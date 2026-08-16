#!/usr/bin/env python
"""Verify a generated customer delivery against its sources. Exits non-zero on violation.

    python scripts/test_customer_delivery.py deliveries/acme-timber/20260812

`test_shared_baseline.py` does this for a processed layer; this does it for a delivery. The
two check different things and neither substitutes for the other -- that one asserts a layer
is shaped right, this one asserts the delivery faithfully carries what the layer said.

The central check is PASSTHROUGH. Every number in `values.csv` is recomputed here from the
source NetCDF with a Gaussian weighting implemented independently of
`spatial_extract.extract_by_point`, and must match bit for bit. That is deliberate: calling
the same function the delivery called would prove only that it is deterministic. An
independent reimplementation is what actually catches a percentile inverted twice, a slope
scaled by 10, a wrong cell, or a broken longitude wrap.

Checks:
  1. Referential integrity across the star schema.
  2. Passthrough of every metric, recomputed independently.
  3. Manifest SHA-256 still matches every source file.
  4. Percentile orientation: layers declaring higher_is_better must anti-correlate value
     with percentile, and every percentile must sit in [1, 100].
  5. slopes_agree matches the active-cell rule.
  6. Baseline-decade slopes are NaN.
  7. data_status agrees with whether `value` is finite.
  8. Taxonomy -> registry linkage: every `covered_by` entry names a real layer.
"""

import hashlib
import html
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

#: DELIBERATELY A SECOND, INDEPENDENT COPY of the scenario -> forcing-tier mapping.
#:
#: This was originally imported from `viz_common` "so the two cannot drift apart" -- which
#: was exactly backwards. Importing the writer's classification means a WRONG mapping (say
#: ssp585 typed as `medium`) is used to both produce and check the Climate Score, and the
#: check passes. A test must restate what it expects; that is what makes it a test. If this
#: list and viz_common.SCENARIO_TIER disagree, one of them is wrong and the delivery fails,
#: which is the correct outcome.
TIER = {
    "rcp26": "low", "ssp126": "low",
    "rcp45": "medium", "rcp60": "medium", "ssp245": "medium", "ssp370": "medium",
    "rcp85": "high", "ssp585": "high",
}

SIGMA = 0.25
SEARCH_RADIUS = 0.5
#: float32 round-trips through CSV as repr(float(x)); allow only representation noise.
TOL = 1e-9

failures: list[str] = []
checks_run = 0


def check(condition: bool, message: str) -> None:
    global checks_run
    checks_run += 1
    if not condition:
        failures.append(message)


def independent_point_extract(ds, var, lat, lon):
    """Gaussian distance-weighted point value, reimplemented from the spec.

    Deliberately NOT importing spatial_extract. Mirrors the documented method: weights
    exp(-0.5 (d/sigma)^2) over cell centres within SEARCH_RADIUS, longitude separation
    wrapped at the antimeridian, NaN cells dropped and remaining weights renormalized.
    """
    lon = ((lon + 180.0) % 360.0) - 180.0
    lats, lons = ds.lat.values, ds.lon.values

    dlat = np.abs(lats - lat)
    dlon = np.abs(lons - lon)
    dlon = np.minimum(dlon, 360.0 - dlon)

    lat_sel, lon_sel = dlat <= SEARCH_RADIUS, dlon <= SEARCH_RADIUS
    if not lat_sel.any() or not lon_sel.any():
        return {}

    dlon_grid, dlat_grid = np.meshgrid(dlon[lon_sel], dlat[lat_sel])
    dist = np.sqrt(dlat_grid**2 + dlon_grid**2)
    weights = np.exp(-0.5 * (dist / SIGMA) ** 2)
    weights = weights / weights.sum()

    block = ds[var].isel(
        lat=np.where(lat_sel)[0], lon=np.where(lon_sel)[0]
    ).values  # (decade, lat, lon)

    out = {}
    for i, decade in enumerate(ds.decade.values):
        panel = block[i]
        valid = ~np.isnan(panel)
        if valid.any():
            w = weights[valid] / weights[valid].sum()
            out[int(decade)] = float(np.sum(panel[valid] * w))
        else:
            out[int(decade)] = np.nan
    return out


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


#: Vocabulary that belongs in our repository and not in a customer's filing. A report is a
#: document somebody else reads under scrutiny; a stray module name or the word UNVERIFIED
#: reads as either noise or an alarm, and neither is what it means internally.
INTERNAL_VOCABULARY = (
    "UNVERIFIED", "this pipeline", "this repository", "codebase", "backlog",
    "(nan)", "process_", "viz_common", "RISK_BANDS", "delivery.py", "TODO",
    "config/", "spatial_extract", "OUTPUT-SPEC", "GUARDRAILS",
)


def independent_vulnerable_counts(values, assets, hazard_layers, threshold):
    """Recompute (tier, decade) -> vulnerable asset count, longhand.

    Deliberately does NOT import report_common.vulnerability_frame. The rule, restated:
    hazard layers only; mean percentile across the native scenarios inside a tier, per
    hazard; MAX across hazards; vulnerable when that maximum reaches the threshold; an asset
    with no finite hazard percentile is not assessed and is not counted either way.
    """
    v = values[values["layer_id"].isin(hazard_layers)].copy()
    v["tier"] = v["scenario"].map(TIER)
    out = {}
    for (tier, decade), grp in v.groupby(["tier", "decade"]):
        n_vuln = n_assessed = 0
        for _asset, ag in grp.groupby("asset_id"):
            worst = None
            for _layer, lg in ag.groupby("layer_id"):
                pcts = [p for p in lg["percentile"] if pd.notna(p)]
                if not pcts:
                    continue
                mean_p = sum(pcts) / len(pcts)
                worst = mean_p if worst is None else max(worst, mean_p)
            if worst is None:
                continue
            n_assessed += 1
            if worst >= threshold:
                n_vuln += 1
        out[(tier, int(decade))] = (n_vuln, n_assessed)
    return out


def check_reports(delivery: Path, manifest: dict, values, assets, layers) -> None:
    """Stage 3 and 4 artifacts, when present."""
    global checks_run
    import yaml

    stages = manifest.get("stages", {})
    cav_path = delivery / "caveats.json"
    compliance = delivery / "report_compliance.html"
    bespoke = delivery / "report_bespoke.html"

    # Stage record and artifact must agree, in BOTH directions.
    #
    # The second direction is the one that bites: re-running the extract rewrites the
    # manifest and resets downstream stages, while the previous run's reports stay in the
    # folder looking finished. A file whose stage is not `built` is a document that does not
    # describe the data next to it.
    for stage, path in (("caveats", cav_path),
                        ("compliance_report", compliance),
                        ("bespoke_report", bespoke)):
        status = stages.get(stage, {}).get("status")
        if status == "built":
            check(path.exists(),
                  f"manifest records stage {stage!r} as built but {path.name} is missing")
        elif path.exists():
            failures.append(
                f"{path.name} exists but stage {stage!r} is {status!r} -- it was built "
                f"against an earlier extract and no longer describes this delivery. "
                f"Rebuild it.")
            checks_run += 1

    if not cav_path.exists():
        if compliance.exists() or bespoke.exists():
            failures.append(
                "a report exists but caveats.json does not. Stage 4 runs BEFORE the "
                "reports and is their input; a report built without it cannot have "
                "carried the required disclosures.")
            checks_run += 1
        return

    payload = json.loads(cav_path.read_text())
    caveats = payload.get("caveats", [])
    check(bool(caveats), "caveats.json contains no caveats")
    ids = [c["id"] for c in caveats]
    check(len(ids) == len(set(ids)), "caveats.json has duplicate ids -- citations are by id")
    check(all(c.get("severity") in {"must_disclose", "should_note", "informational"}
              for c in caveats),
          "a caveat carries an unknown severity")
    check(all(c.get("text", "").strip() for c in caveats), "a caveat has empty text")
    check(stages.get("caveats", {}).get("status") == "built",
          "caveats.json exists but the caveats stage is not recorded as built")
    must = [c["id"] for c in caveats if c["severity"] == "must_disclose"]
    check(bool(must), "no caveat is marked must_disclose -- every delivery has at least "
                      "the coverage and resolution disclosures")
    print(f"  caveats.json: {len(caveats)} caveats, {len(must)} must-disclose")

    cfg_path = delivery / "report_config.yaml"
    cfg = yaml.safe_load(cfg_path.read_text()) if cfg_path.exists() else {}
    check(bool(cfg), "report_config.yaml is missing or empty")

    taxonomy = yaml.safe_load((PROJECT_ROOT / "config" / "hazard_taxonomy.yaml").read_text())
    non_hazard = set(taxonomy.get("non_hazard_layers") or {})
    hazard_layers = [l for l in layers.index if l not in non_hazard]  # layers is indexed by layer_id
    check(len(hazard_layers) < len(layers) or not non_hazard,
          "no layer was excluded as a non-hazard, but the taxonomy declares some")

    for path in (compliance, bespoke):
        if not path.exists():
            continue
        html_text = path.read_text()

        missing = [cid for cid in must if cid not in html_text]
        check(not missing,
              f"{path.name} omits must-disclose caveat(s): {', '.join(missing)}")

        stamp = re.search(r'class="stamp">build ([0-9a-f]{8})<', html_text)
        check(stamp is not None, f"{path.name} carries no build stamp")

        visible = re.sub(r"<[^>]+>", " ", html_text)
        leaked = [t for t in INTERNAL_VOCABULARY if t in visible]
        check(not leaked,
              f"{path.name} leaks internal vocabulary into a customer document: "
              f"{', '.join(leaked)}")

        check("<h1>" in html_text and "</html>" in html_text,
              f"{path.name} is not a complete HTML document")
        check("<!--" not in html_text,
              f"{path.name} contains an HTML comment -- narrative guidance must be stripped")
        print(f"  {path.name} checked ({path.stat().st_size // 1024} KB)")

    # While the vulnerability method is deferred, NO report may publish a determination.
    #
    # This is the check that keeps a considered "we have not decided" from quietly becoming
    # a number again. A percentile threshold produces a vulnerable-asset count very easily,
    # and the count looks authoritative; the whole point of deferring is that nothing
    # connects a global exposure rank to susceptibility to harm.
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    from utils.report_common import TBD_SECTIONS  # noqa: E402

    if "vulnerability_metric" in TBD_SECTIONS:
        for path in (compliance, bespoke):
            if not path.exists():
                continue
            # Test each TEXT NODE, not the flattened page. Flattening replaces tags with
            # spaces, so adjacent table cells run together and "Sections 2, 4, 7" followed
            # by a cell reading "Assets vulnerable" matched a count-of-vulnerable-assets
            # pattern that nothing had published. Per-node testing removes that class of
            # false positive without losing sensitivity: a real determination -- a caption,
            # a callout sentence, a determination cell -- is always within one node.
            nodes = [
                " ".join(html.unescape(n).split())
                for n in re.split(r"<[^>]+>", path.read_text())
            ]
            nodes = [n for n in nodes if n]

            # A literal phrase list was the wrong shape here and missed a must-disclose
            # caveat telling the customer "this report therefore discloses counts and
            # percentages of assets only" -- a direct contradiction of the deferral, in both
            # reports, past a passing verifier.
            patterns = [
                (r"\d+\s*(?:of\s*\d+\s*)?(?:assets?\s+)?vulnerable", "a count of vulnerable assets"),
                (r"vulnerable[^.]{0,40}?\b\d+(?:\.\d+)?\s*%", "a percentage vulnerable"),
                (r"^\s*Not vulnerable\s*$", "a per-asset not-vulnerable determination"),
                (r"assets? vulnerable to[^.]{0,60}?threshold", "a vulnerability table caption"),
                (r"sensitivity of the vulnerable count", "a threshold-sensitivity table"),
                (r"(?:disclos|report|present)\w*\s+(?:only\s+)?(?:the\s+)?"
                 r"(?:count|percentage|amount|share)s?\b[^.]{0,80}?\basset",
                 "a claim that counts/percentages of assets are disclosed"),
            ]
            published = sorted({
                why for pat, why in patterns
                for n in nodes if re.search(pat, n, re.I)
            })
            check(not published,
                  f"{path.name} publishes or claims a vulnerability determination "
                  f"({'; '.join(published)}) while the method is recorded as deferred in "
                  f"TBD_SECTIONS. Either agree the method and remove the TBD entry, or do "
                  f"not print the figure.")
            check("method to be agreed" in visible.lower()
                  or "not yet determined" in visible.lower(),
                  f"{path.name} does not state that the vulnerability method is deferred")
        print("  vulnerability metric correctly reported as deferred, not published")
    else:
        # Method agreed: the printed counts must match an independent recomputation.
        threshold = float((cfg.get("vulnerability") or {}).get("threshold", 80))
        expected = independent_vulnerable_counts(values, assets, hazard_layers, threshold)
        html_text = compliance.read_text() if compliance.exists() else ""
        m = re.search(r"threshold[^<]*?(\d+)th percentile\.(.*?)</table>", html_text, re.S)
        check(m is not None, "compliance report has no vulnerability table to check")
        if m:
            check(int(m.group(1)) == int(threshold),
                  f"compliance report states threshold {m.group(1)} but report_config.yaml "
                  f"says {threshold:.0f}")
            horizons = {
                str(v.get("label", k)): int(v["decade"])
                for k, v in (cfg.get("horizons") or {}).items()
                if isinstance(v, dict) and "decade" in v
            }
            tiers = {"Low forcing": "low", "Medium forcing": "medium", "High forcing": "high"}
            compared = 0
            for row in re.findall(r"<tr>(.*?)</tr>", m.group(2), re.S):
                cells = [re.sub(r"<[^>]+>", "", c).strip()
                         for c in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)]
                if len(cells) < 4 or cells[0] not in horizons or cells[1] not in tiers:
                    continue
                decade, tier = horizons[cells[0]], tiers[cells[1]]
                exp_vuln, exp_assessed = expected.get((tier, decade), (0, 0))
                check(int(cells[2]) == exp_assessed,
                      f"compliance report says {cells[2]} assets assessed at "
                      f"{tier}/{decade}s; independently computed {exp_assessed}")
                check(int(cells[3]) == exp_vuln,
                      f"compliance report says {cells[3]} assets vulnerable at "
                      f"{tier}/{decade}s; independently computed {exp_vuln}")
                compared += 1
            check(compared > 0, "no vulnerability rows could be parsed and compared")
            print(f"  vulnerability metric independently recomputed for {compared} "
                  f"horizon x tier cells")

    # The bespoke report's narrative must still be complete and every citation resolve.
    if bespoke.exists():
        narrative = delivery / "narrative.md"
        check(narrative.exists(),
              "report_bespoke.html exists but narrative.md does not -- the report cannot "
              "be rebuilt or audited")
        if narrative.exists():
            sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
            from utils.report_common import check_narrative, load_delivery  # noqa: E402
            problems = check_narrative(load_delivery(delivery), narrative.read_text())
            check(not problems,
                  "narrative.md no longer validates: " + "; ".join(problems[:5]))
        html_text = bespoke.read_text()
        anchors = {int(x) for x in re.findall(r'id="ref-(\d+)"', html_text)}
        used = {int(x) for x in re.findall(r'href="#ref-(\d+)"', html_text)}
        check(used <= anchors,
              f"bespoke report cites reference(s) with no entry: {sorted(used - anchors)}")
        check(anchors == set(range(1, len(anchors) + 1)) if anchors else True,
              "bespoke report reference numbering is not contiguous from 1")


def check_taxonomy_registry_link() -> None:
    """Every `covered_by` entry names a layer the registry actually has.

    This is a property of the repo config rather than of THIS delivery, and it belongs in the
    verifier anyway because its failure mode is invisible everywhere else and ends up in a
    customer document. `report_common.coverage_summary()` intersects a family's `covered_by`
    with the layers a delivery actually carries, so an entry naming a layer the registry does
    not have can never match anything: the family reports as NOT ASSESSED in every report
    while `hazard_taxonomy.yaml` asserts it is covered. Both files look correct read on their
    own, which is why this needs a machine.

    It is not hypothetical. `permafrost-3b` shipped 2026-08-14 -- processed, documented, and
    named in `families.permafrost-thaw.covered_by` -- and was not added to the registry until
    2026-08-16. For two days the taxonomy claimed a hazard the registry had never heard of,
    and every report built in that window would have filed permafrost thaw under "hazards not
    assessed" while the taxonomy said otherwise.

    Read from the YAML rather than through `load_registry()` / `load_hazard_taxonomy()`, for
    the reason spelled out at TIER above: asking the writer's own loader what the writer wrote
    proves that the loader is deterministic, not that the two files agree.

    NOT CHECKED HERE, DELIBERATELY: the reverse direction. A registered layer that no family
    lists is legitimate -- `csoil` is deliberately uncovered pending a hazard-vs-asset-condition
    decision, and only one rung of a threshold ladder is ever the `preferred` one. Asserting
    the reverse would fail on layers that are correct. It is worth auditing by hand when a
    layer ships; it is not a contract.
    """
    registry = yaml.safe_load((PROJECT_ROOT / "config" / "layer_registry.yaml").read_text())
    taxonomy = yaml.safe_load((PROJECT_ROOT / "config" / "hazard_taxonomy.yaml").read_text())
    known = set(registry.get("layers") or {})
    check(bool(known), "config/layer_registry.yaml declares no layers")

    for family, spec in (taxonomy.get("families") or {}).items():
        for layer_id in spec.get("covered_by") or []:
            check(layer_id in known,
                  f"hazard_taxonomy family {family!r} is covered_by {layer_id!r}, which is "
                  f"not in config/layer_registry.yaml -- the family will report as NOT "
                  f"ASSESSED in every report while the taxonomy claims it is covered")

    for layer_id in taxonomy.get("non_hazard_layers") or {}:
        check(layer_id in known,
              f"hazard_taxonomy non_hazard_layers names {layer_id!r}, which is not in "
              f"config/layer_registry.yaml -- it cannot be excluded from a hazard count "
              f"it can never appear in")


def main(delivery: Path) -> int:
    global checks_run
    for name in ("locations.csv", "assets.csv", "layers.csv", "values.csv", "manifest.json"):
        if not (delivery / name).exists():
            print(f"FAIL: {delivery / name} is missing")
            return 1

    locations = pd.read_csv(delivery / "locations.csv")
    assets = pd.read_csv(delivery / "assets.csv")
    layers = pd.read_csv(delivery / "layers.csv").set_index("layer_id")
    values = pd.read_csv(delivery / "values.csv")
    manifest = json.loads((delivery / "manifest.json").read_text())

    marker = delivery / "DELIVERY-INCOMPLETE.md"
    if marker.exists():
        first = marker.read_text().strip().splitlines()
        failures.append(
            "delivery is marked incomplete and must not be shipped: "
            + (first[2] if len(first) > 2 else marker.name))

    print(f"Delivery: {delivery}")
    print(f"  {len(locations)} locations, {len(assets)} assets, "
          f"{len(layers)} layers, {len(values)} value rows")

    # 0. config linkage -- checked before the delivery, because it is upstream of it -------
    check_taxonomy_registry_link()

    # 1. referential integrity ----------------------------------------------------------
    check(assets["location_id"].isin(locations["location_id"]).all(),
          "assets.csv references a location_id absent from locations.csv")
    check(values["asset_id"].isin(assets["asset_id"]).all(),
          "values.csv references an asset_id absent from assets.csv")
    check(values["layer_id"].isin(layers.index).all(),
          "values.csv references a layer_id absent from layers.csv")
    check(locations["location_id"].is_unique, "locations.csv has duplicate location_id")
    check(assets["asset_id"].is_unique, "assets.csv has duplicate asset_id")
    check(not values.duplicated(["asset_id", "layer_id", "scenario", "decade"]).any(),
          "values.csv has duplicate (asset_id, layer_id, scenario, decade) keys")

    # 1b. COMPLETENESS -- the delivery must contain every row it should ------------------
    #
    # Every other check validates the rows that are PRESENT. Delete every rcp85 row from
    # values.csv (and its climate_score rows) and everything still recomputes correctly:
    # a whole scenario can vanish and the verifier passes. So reconstruct the expected key
    # set from the inputs and the source files, and compare.
    expected_keys = set()
    asset_layer_map = {
        r["asset_id"]: [x for x in str(r["layer_ids"]).split(";") if x]
        for _, r in assets.iterrows()
    }
    for asset_id, layer_ids in asset_layer_map.items():
        for layer_id in layer_ids:
            if layer_id not in layers.index:
                failures.append(f"{asset_id} claims layer {layer_id!r}, absent from layers.csv")
                continue
            scenarios = [x for x in str(layers.at[layer_id, "scenarios"]).split(";") if x]
            decades = [int(x) for x in str(layers.at[layer_id, "decades"]).split(";") if x]
            for sc in scenarios:
                for dec in decades:
                    expected_keys.add((asset_id, layer_id, sc, dec))
    actual_keys = {(r["asset_id"], r["layer_id"], r["scenario"], int(r["decade"]))
                   for _, r in values.iterrows()}
    missing_keys = expected_keys - actual_keys
    extra_keys = actual_keys - expected_keys
    check(not missing_keys,
          f"values.csv is MISSING {len(missing_keys)} expected row(s), e.g. "
          f"{sorted(missing_keys)[:3]}")
    check(not extra_keys,
          f"values.csv has {len(extra_keys)} row(s) that should not exist, e.g. "
          f"{sorted(extra_keys)[:3]}")
    check(int(manifest.get("counts", {}).get("value_rows", -1)) == len(values),
          f"manifest says {manifest.get('counts', {}).get('value_rows')} value rows, "
          f"values.csv has {len(values)}")
    check(int(manifest.get("counts", {}).get("locations", -1)) == len(locations),
          "manifest location count disagrees with locations.csv")
    check(int(manifest.get("counts", {}).get("assets", -1)) == len(assets),
          "manifest asset count disagrees with assets.csv")

    # 3. source hashes -------------------------------------------------------------------
    for src in manifest["source_files"]:
        p = PROJECT_ROOT / src["path"]
        if not p.exists():
            failures.append(f"source file vanished: {src['path']}")
            continue
        checks_run += 1
        actual = sha256(p)
        if actual != src.get("sha256"):
            failures.append(
                f"source changed since delivery: {src['path']}\n"
                f"    manifest {src.get('sha256')}\n    actual   {actual}"
            )

    # 2. passthrough ---------------------------------------------------------------------
    coords = (
        assets.merge(locations, on="location_id")
        .set_index("asset_id")[["lat", "lon"]]
        .to_dict("index")
    )
    metric_to_var = {
        "value": "median", "lower_ci": "lower_ci", "upper_ci": "upper_ci",
        "percentile": "percentile", "ols_slope": "ols_slope", "sen_slope": "sen_slope",
        "n_members": "n_members", "n_models": "n_models",
    }
    counts = {"n_members", "n_models"}

    n_compared = 0
    # Bind to the EXACT file the manifest says was read, not to whatever a glob returns
    # first. A folder holding two files that match `*_ssp585_processed.nc` would otherwise
    # be validated against an arbitrary one, nondeterministically by directory order.
    source_by_key = {(src["layer_id"], src["scenario"]): PROJECT_ROOT / src["path"]
                     for src in manifest["source_files"]}
    for (layer_id, scenario), grp in values.groupby(["layer_id", "scenario"]):
        src_path = source_by_key.get((layer_id, scenario))
        if src_path is None:
            failures.append(
                f"{layer_id}/{scenario} has rows in values.csv but no manifest source entry")
            continue
        if not src_path.exists():
            failures.append(f"manifest source missing on disk: {src_path}")
            continue
        with xr.open_dataset(src_path) as ds:
            for asset_id, arows in grp.groupby("asset_id"):
                lat, lon = coords[asset_id]["lat"], coords[asset_id]["lon"]
                for metric, var in metric_to_var.items():
                    expected = independent_point_extract(ds, var, lat, lon)
                    for _, r in arows.iterrows():
                        exp = expected.get(int(r["decade"]), np.nan)
                        got = r[metric]
                        if metric in counts and np.isfinite(exp):
                            exp = int(round(exp))
                        got_nan = pd.isna(got)
                        exp_nan = pd.isna(exp)
                        checks_run += 1
                        n_compared += 1
                        if got_nan != exp_nan:
                            failures.append(
                                f"{asset_id} {layer_id} {scenario} {r['decade']} {metric}: "
                                f"delivered {got!r}, independent recompute {exp!r}")
                        elif not got_nan and abs(float(got) - float(exp)) > TOL:
                            failures.append(
                                f"{asset_id} {layer_id} {scenario} {r['decade']} {metric}: "
                                f"delivered {got!r}, independent recompute {exp!r}")

    # 4. percentile orientation and range ------------------------------------------------
    pct = values["percentile"].dropna()
    check(pct.between(1, 100).all(),
          f"percentile outside [1, 100]: min {pct.min()}, max {pct.max()}")
    for layer_id, grp in values.groupby("layer_id"):
        direction = str(layers.at[layer_id, "percentile_direction"])
        sub = grp[["value", "percentile"]].dropna()
        if len(sub) < 3 or sub["value"].nunique() < 2:
            continue
        r = sub["value"].corr(sub["percentile"])
        if direction == "higher_is_better":
            check(r < 0,
                  f"{layer_id} declares higher_is_better but value/percentile correlate "
                  f"{r:+.3f} (expected negative -- the inversion may have been applied "
                  f"twice, or not at all)")
        elif direction == "higher_is_worse":
            check(r > 0,
                  f"{layer_id} declares higher_is_worse but value/percentile correlate "
                  f"{r:+.3f} (expected positive)")

    # 5. slopes_agree --------------------------------------------------------------------
    for _, r in values.iterrows():
        ols, sen, agree = r["ols_slope"], r["sen_slope"], r["slopes_agree"]
        if pd.isna(ols) or pd.isna(sen) or (ols == 0 and sen == 0):
            expected = None
        elif ols == 0 or sen == 0:
            expected = False
        else:
            expected = (ols > 0) == (sen > 0)
        actual = None if pd.isna(agree) else bool(agree)
        checks_run += 1
        if actual != expected:
            failures.append(
                f"{r['asset_id']} {r['layer_id']} {r['scenario']} {r['decade']} "
                f"slopes_agree={actual} but ols={ols}, sen={sen} implies {expected}")

    # 6. baseline slopes are NaN ---------------------------------------------------------
    for layer_id, grp in values.groupby("layer_id"):
        baseline = int(layers.at[layer_id, "baseline_decade"])
        base_rows = grp[grp["decade"] <= baseline]
        check(base_rows["ols_slope"].isna().all() and base_rows["sen_slope"].isna().all(),
              f"{layer_id}: slopes are finite at or before the baseline decade {baseline}; "
              f"the expanding window has no span there")

    # 7. data_status ---------------------------------------------------------------------
    check((values.loc[values["data_status"] == "OK", "value"].notna()).all(),
          "a row marked OK has no value")
    check((values.loc[values["data_status"] != "OK", "value"].isna()).all(),
          "a row with a value is not marked OK")

    # 8. Climate Score recomputed independently ------------------------------------------
    score_path = delivery / "climate_score.csv"
    if score_path.exists():
        scores = pd.read_csv(score_path)
        asset_layers = {
            r["asset_id"]: {x for x in str(r["layer_ids"]).split(";") if x}
            for _, r in assets.iterrows()
        }
        vv = values[values["percentile"].notna()].copy()
        vv["tier"] = vv["scenario"].map(TIER)
        check(vv["tier"].notna().all(),
              "a scenario in values.csv has no forcing tier mapping")

        # Two-stage mean, written from the spec rather than imported: average scenario
        # codes WITHIN a hazard, then average the hazards. A flat one-stage mean would
        # double-weight any hazard contributing two codes to one tier.
        stage1 = (vv.groupby(["asset_id", "layer_id", "tier", "decade"])["percentile"]
                  .mean().reset_index())
        stage2 = (stage1.groupby(["asset_id", "tier", "decade"])["percentile"]
                  .agg(["mean", "count"]).reset_index())
        expect = {(r["asset_id"], r["tier"], int(r["decade"])): (r["mean"], int(r["count"]))
                  for _, r in stage2.iterrows()}

        check(len(scores) == len(expect),
              f"climate_score.csv has {len(scores)} rows, independent recompute has "
              f"{len(expect)}")
        for _, r in scores.iterrows():
            k = (r["asset_id"], r["scenario_tier"], int(r["decade"]))
            checks_run += 1
            if k not in expect:
                failures.append(f"climate_score row {k} has no counterpart in values.csv")
                continue
            exp_mean, exp_n = expect[k]
            if abs(float(r["climate_score"]) - exp_mean) > 0.005:
                failures.append(
                    f"climate_score {k}: stored {r['climate_score']}, recomputed "
                    f"{exp_mean:.4f}")
            if int(r["n_hazards"]) != exp_n:
                failures.append(
                    f"climate_score {k}: n_hazards {r['n_hazards']}, recomputed {exp_n}")
            # An asset must never be scored on a hazard outside its catalog set.
            listed = {x for x in str(r["hazards"]).split(";") if x}
            if not listed <= asset_layers.get(r["asset_id"], set()):
                failures.append(
                    f"climate_score {k} names hazard(s) not in the asset's set: "
                    f"{listed - asset_layers.get(r['asset_id'], set())}")
        s = scores["climate_score"]
        check(s.between(1, 100).all(),
              f"climate_score outside [1,100]: min {s.min()}, max {s.max()}")
        check((scores["n_hazards"] >= 1).all(), "a climate_score row rests on 0 hazards")

        # BASELINE INVARIANT. The shared 2020s panel is bit-identical across scenarios, so
        # a given asset's baseline Climate Score MUST be the same in every tier where it
        # has the same hazard set. A difference means composition changed, not risk --
        # which is precisely the artifact that made a portfolio's high-tier 2020s mean read
        # 39.9 against 42.1 for low and medium before the balanced panel landed.
        base_dec = int(layers["baseline_decade"].dropna().astype(int).max())
        base = scores[scores["decade"] == base_dec]
        for (asset_id, n_haz), grp in base.groupby(["asset_id", "n_hazards"]):
            if len(grp) < 2:
                continue
            checks_run += 1
            spread = grp["climate_score"].max() - grp["climate_score"].min()
            if spread > 0.02:
                failures.append(
                    f"{asset_id}: baseline ({base_dec}s) climate_score differs across "
                    f"tiers on the same {n_haz}-hazard set "
                    f"({grp['climate_score'].tolist()}) -- the shared 2020s panel is "
                    f"bit-identical across scenarios, so this is a composition artifact")
        # BALANCED-PANEL INVARIANT, the aggregation form of the check above.
        #
        # The dashboard averages Climate Scores up to a location and to the portfolio. Any
        # such rollup must first restrict to assets present in every tier being compared,
        # or an asset dropping out of one tier reads as a risk difference. This has bitten
        # twice: once at portfolio level (high 39.9 vs low 42.1) and once per location
        # (Shasta high 62.3 vs low 51.7), both at the BASELINE decade where the tiers are
        # required to be identical. The guard is that equality, evaluated on the same
        # balanced panels the dashboard builds.
        nexp = {r["asset_id"]: len([x for x in str(r["layer_ids"]).split(";") if x])
                for _, r in assets.iterrows()}
        complete = {(r["asset_id"], r["scenario_tier"], int(r["decade"]))
                    for _, r in scores.iterrows()
                    if int(r["n_hazards"]) >= nexp.get(r["asset_id"], 0)}
        all_tiers = sorted(scores["scenario_tier"].unique())
        all_decs = sorted(scores["decade"].unique())
        loc_of = dict(zip(assets["asset_id"], assets["location_id"]))

        groups = {"__portfolio": list(assets["asset_id"])}
        for aid, loc in loc_of.items():
            groups.setdefault(loc, []).append(aid)

        for gname, members in groups.items():
            decs = [d for d in all_decs
                    if any(all((m, t, d) in complete for t in all_tiers) for m in members)]
            if not decs:
                continue
            panel = [m for m in members
                     if all((m, t, d) in complete for t in all_tiers for d in decs)]
            if not panel or base_dec not in decs:
                continue
            sub = scores[(scores["asset_id"].isin(panel)) & (scores["decade"] == base_dec)]
            by_tier = sub.groupby("scenario_tier")["climate_score"].mean()
            checks_run += 1
            if len(by_tier) > 1 and (by_tier.max() - by_tier.min()) > 0.02:
                failures.append(
                    f"balanced-panel baseline ({base_dec}s) differs across tiers for "
                    f"{gname}: {by_tier.round(2).to_dict()} on a panel of {len(panel)} "
                    f"asset(s) -- a rollup is mixing different asset sets across tiers")

        print(f"  climate_score.csv: {len(scores)} rows recomputed independently")
    else:
        print("  no climate_score.csv")

    # 9. dashboard, when one has been generated ------------------------------------------
    dash = delivery / "dashboard.html"
    if dash.exists():
        html = dash.read_text()
        payload = re.search(r"const DATA = (\{.*?\});\n", html, re.S)
        check(payload is not None, "dashboard.html has no embedded DATA payload")
        if payload:
            try:
                embedded = json.loads(payload.group(1))
            except json.JSONDecodeError as exc:
                failures.append(f"dashboard payload is not valid JSON: {exc}")
                embedded = None
            if embedded:
                # The dashboard dedupes values to (location, layer, scenario, decade)
                # because extraction is coordinate-driven -- two assets at one site read
                # identically. Assert that, rather than trusting it.
                key = ["location_id", "layer_id", "scenario", "decade"]
                merged = values.merge(
                    assets[["asset_id", "location_id"]], on="asset_id")
                collisions = merged.groupby(key)["value"].nunique(dropna=False)
                check((collisions <= 1).all(),
                      "two assets at the same site disagree on a value for the same "
                      "layer/scenario/decade -- the dashboard's location-level dedup "
                      "would silently drop one")
                check(len(embedded["values"]) == len(merged.drop_duplicates(key)),
                      f"dashboard payload has {len(embedded['values'])} rows, deduped "
                      f"delivery has {len(merged.drop_duplicates(key))}")
                pcts = [r["p"] for r in embedded["values"] if r["p"] is not None]
                check(not pcts or (min(pcts) >= 1 and max(pcts) <= 100),
                      "dashboard payload carries an out-of-range percentile")
        style = re.search(r"<style>(.*?)</style>", html, re.S)
        check(style is not None and "{{" not in style.group(1)
              and "}}" not in style.group(1),
              "dashboard CSS contains leaked f-string escapes")
        ids_used = set(re.findall(r'\$\("([a-z0-9\-]+)"\)', html))
        ids_defined = set(re.findall(r'id="([a-z0-9\-]+)"', html))
        check(not (ids_used - ids_defined),
              f"dashboard JS references missing element ids: {ids_used - ids_defined}")
        print(f"  dashboard.html checked ({dash.stat().st_size // 1024} KB)")
    else:
        # Not optional. A delivery is CSVs *and* dashboard -- shipping bare CSVs leaves the
        # customer, and us, with no way to look at what was extracted.
        failures.append(
            "dashboard.html is missing. Every delivery must ship one; rebuild with "
            "`python scripts/generate_delivery_dashboard.py <delivery>`")
        checks_run += 1

    # 10. delivery stages -----------------------------------------------------------------
    stages = manifest.get("stages", {})
    check(bool(stages), "manifest records no delivery stages")
    for req in ("inputs", "extract", "dashboard"):
        st = stages.get(req, {}).get("status")
        check(st in {"built", "confirmed"},
              f"stage {req!r} is {st!r}, expected built/confirmed -- delivery incomplete")
    pending = [k for k, v in stages.items() if v.get("status") == "not_started"]
    if pending:
        print(f"  stages not started (expected, not built yet): {', '.join(sorted(pending))}")

    # 11. Stage 3 and 4 artifacts ----------------------------------------------------------
    #
    # Same philosophy as everywhere else in this file: RESTATE what is expected rather than
    # importing it. The vulnerable-asset count is recomputed here from values.csv with the
    # rule written out longhand, and compared against the number the compliance report
    # actually printed. That is the end-to-end proof for the headline disclosure metric --
    # everything else about a report can be right while the table says something else.
    check_reports(delivery, manifest, values, assets, layers)

    print(f"  {n_compared} metric values independently recomputed")
    print(f"  {checks_run} checks run")
    if failures:
        print(f"\nFAIL -- {len(failures)} violation(s):")
        for f in failures[:40]:
            print(f"  - {f}")
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
        return 1
    print("\nPASS -- delivery faithfully carries its sources")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(Path(sys.argv[1])))
