#!/usr/bin/env python3
"""Stage 4 -- derive this delivery's caveat set. Runs BEFORE the reports, not after.

    python scripts/generate_delivery_caveats.py deliveries/<customer>/<date>

WHY THIS IS A STAGE AND NOT A PARAGRAPH SOMEBODY REMEMBERS TO WRITE
-------------------------------------------------------------------
Everything below is already recorded somewhere -- in the manifest, in a layer's
`interpretation_caveat`, in the registry, in the report config, in
config/hazard_taxonomy.yaml. Scattered across six files it is inert: nobody assembles it
under deadline, and a caveat nobody assembles is a caveat the customer never reads.

So it is assembled mechanically, once, into `caveats.json`, and BOTH reports are then
required to carry every `must_disclose` entry -- `report_common.assert_caveats_present()`
refuses to render otherwise and the verifier re-checks it afterwards. That ordering is the
point of moving this stage ahead of the reports: the caveat set is an INPUT to a report, not
a summary of one. If each report derived its own list, the two documents describing one
delivery would eventually disagree about what is wrong with it.

IDs ARE SEMANTIC AND STABLE
---------------------------
`HAZARD-COVERAGE`, not `CAV-001`. A narrative cites a caveat by id
(`[data:caveats.json#HAZARD-COVERAGE]`), so a numbering scheme that shifts when a caveat is
added would silently repoint every citation in every previously written narrative.

WHAT THIS DOES NOT DO
---------------------
It does not judge materiality. A caveat's severity says how certainly it must be disclosed,
never how much it matters to this customer's decision -- that judgement needs the business
context and belongs in the bespoke report's narrative, with a citation.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.delivery import DeliveryError, record_stage  # noqa: E402
from utils.report_common import (  # noqa: E402
    CAVEATS_JSON,
    CAVEATS_MD,
    NOT_ASSESSED,
    SEVERITY_INFO,
    SEVERITY_MUST,
    SEVERITY_SHOULD,
    Delivery,
    coverage_summary,
    hazard_layer_ids,
    load_delivery,
    load_hazard_taxonomy,
    vulnerability_frame,
)
from utils.viz_common import SCENARIO_TIER, TIER_LABELS, TIER_ORDER  # noqa: E402


def _cav(
    cid: str,
    severity: str,
    scope: str,
    title: str,
    text: str,
    evidence: str = "",
    affects: Optional[List[str]] = None,
) -> dict:
    return {
        "id": cid,
        "severity": severity,
        "scope": scope,
        "title": title,
        "text": " ".join(text.split()),
        "evidence": " ".join(evidence.split()),
        "affects": affects or [],
    }


# ---------------------------------------------------------------------------------------
# The derivations
# ---------------------------------------------------------------------------------------


def hazard_coverage(delivery: Delivery, taxonomy: dict) -> List[dict]:
    cov = coverage_summary(delivery, taxonomy)
    covered = [c["name"] for c in cov["covered"]]
    uncovered = cov["uncovered"]
    named = [u["name"] for u in uncovered if u.get("materiality_note")]
    out = [
        _cav(
            "HAZARD-COVERAGE",
            SEVERITY_MUST,
            "delivery",
            "This assessment covers three hazard families, not all of them",
            f"""
            The hazards assessed for this portfolio are: {', '.join(covered)}. Every other
            physical hazard was NOT assessed, and its absence from this report is not a
            finding that it is immaterial -- it was not examined. The omissions most likely
            to matter are {', '.join(named)}. A reader comparing sites within this report is
            comparing them on the hazards listed, and on nothing else. The full list of
            hazards outside this assessment is set out in its own section; note that the
            grouping of hazards into families is ours rather than any standard's, so the
            COUNT of families is not a meaningful coverage ratio and none is quoted.
            """,
            evidence=(
                "config/hazard_taxonomy.yaml, scoped to the layers this delivery actually "
                f"extracted: {', '.join(hazard_layer_ids(delivery, taxonomy))}"
            ),
            affects=[u["id"] for u in uncovered],
        )
    ]
    if not cov["cdp_labels_verified"]:
        out.append(
            _cav(
                # Neutral identifier. Caveat ids are printed in the customer document for
                # traceability, so an id is customer-facing text and "UNVERIFIED" reads as
                # an alarm rather than as the precise bookkeeping term it is internally.
                "CDP-HAZARD-LABELS",
                SEVERITY_SHOULD,
                "delivery",
                "The CDP hazard labels in the mapping appendix are a mapping, not a transcription",
                """
                The hazard families in this assessment were mapped onto CDP's question 3.1.1
                hazard vocabulary by name. The dropdown in the current questionnaire was not
                read option by option, so a label here may not match a label there exactly.
                Before submitting, open the questionnaire and reconcile the two lists.
                """,
                evidence=cov["cdp_labels_note"],
            )
        )
    return out


def qa_and_catalog(delivery: Delivery) -> List[dict]:
    out: List[dict] = []
    qa = delivery.manifest.get("qa_review", {})
    unreviewed = sorted(k for k, v in qa.items() if v in (None, "", "NOT CONFIRMED"))
    if unreviewed:
        out.append(
            _cav(
                "QA-UNREVIEWED",
                SEVERITY_MUST,
                "delivery",
                "No layer in this delivery has completed human QA review",
                f"""
                {len(unreviewed)} of {len(qa)} layers carry no recorded QA review:
                {', '.join(unreviewed)}. The automated checks confirm that a dataset is
                SHAPED correctly -- right dimensions, right statistics, an internally
                consistent baseline. They cannot confirm that the underlying model output is
                about what its name says. That distinction is not hypothetical: a crop layer
                has passed every automated check in this programme and was withdrawn on
                inspection, because the model had not simulated the crop anywhere in its main
                growing regions. Treat the values in this report as provisional until a
                reviewer has examined the layers and their maps.
                """,
                evidence="manifest.json qa_review; config/layer_registry.yaml qa_reviewed_on",
                affects=unreviewed,
            )
        )
    catalog_sourced = delivery.assets[delivery.assets["layer_source"] == "catalog"]
    if len(catalog_sourced):
        types = sorted(set(catalog_sourced["catalog_entry"]))
        out.append(
            _cav(
                "CATALOG-UNCONFIRMED",
                SEVERITY_MUST,
                "delivery",
                "The asset-to-hazard mapping is a draft that has not been signed off",
                f"""
                Which hazards were assessed for each asset came from our asset catalog
                ({', '.join(types)}). Every seeded entry carries confirmed_on: null, meaning
                no one has approved it. That mapping is a claim about which hazards have a
                transmission channel to a given asset -- a warehouse is deliberately not
                assessed for drought, for instance -- and the customer is better placed than
                we are to say whether it is right for their operations. Confirm it.
                """,
                evidence="config/asset_catalog.yaml -- every entry confirmed_on: null",
                affects=types,
            )
        )
    return out


def method_caveats(delivery: Delivery) -> List[dict]:
    """Caveats about how a number was made. These are true of every delivery."""
    out = [
        _cav(
            "SPATIAL-FOOTPRINT",
            SEVERITY_MUST,
            "delivery",
            "Each value describes roughly a 1° × 1° area, not the site itself",
            """
            The underlying model grid is 0.5° (about 55 km at the equator) and each site's
            value is a Gaussian-weighted blend of the four surrounding grid-cell centres --
            a footprint of about 1° × 1°, or 111 km north-south. This is a screening
            instrument for comparing sites and horizons, not a site-specific engineering
            assessment, and it cannot resolve local topography, drainage, defensible space
            or flood defences. Coordinate precision matters accordingly: shifting the sites
            in this portfolio by a quarter of a degree changed their 2090s burnt-area values
            by between 44% and 569%. A wrong coordinate is not a rounding error.
            """,
            evidence=(
                "manifest.json extraction block: Gaussian weighting, sigma 0.25°, search "
                "radius 0.5°, unmodelled neighbouring cells excluded and the remaining "
                "weights rescaled. Both figures above are reproducible on demand -- 20,000 "
                "random sites all resolve to a four-cell blend, and the coordinate shifts "
                "were measured on this portfolio's own sites."
            ),
        ),
        _cav(
            "PERCENTILE-BASIS",
            SEVERITY_MUST,
            "delivery",
            "Percentile ranks a site against global land, and is a level rather than a change",
            """
            Every percentile in this report ranks that site's decadal value against the
            distribution of all global land cells in the shared 2020s baseline. A percentile
            of 80 means "worse than 80% of global land was in the 2020s"; it does not mean
            the site has worsened by 80%, and it does not mean the site has worsened at all.
            Level and trend are separate axes and are reported separately: a site can sit at
            the 95th percentile and be almost flat, while another climbs steeply and remains
            less exposed in absolute terms.
            """,
            evidence="OUTPUT-SPEC.md — percentile-of-score against the shared 2020s global land baseline",
        ),
    ]

    rounds = {
        "2b/RCP" if s.startswith("rcp") else "3b/SSP"
        for s in delivery.values["scenario"].unique()
    }
    if len(rounds) > 1:
        by_round: Dict[str, List[str]] = {}
        for lid in delivery.layer_ids:
            scen = sorted(delivery.values[delivery.values["layer_id"] == lid]["scenario"].unique())
            key = "ISIMIP2b (RCP)" if scen and scen[0].startswith("rcp") else "ISIMIP3b (SSP)"
            by_round.setdefault(key, []).append(f"{lid} ({', '.join(scen)})")
        detail = "; ".join(f"{k}: {', '.join(v)}" for k, v in sorted(by_round.items()))
        out.append(
            _cav(
                "SCENARIO-TIERS",
                SEVERITY_MUST,
                "delivery",
                "Hazards come from two model generations, grouped into forcing tiers that are only approximately comparable",
                f"""
                This portfolio draws on both ISIMIP rounds, and no scenario code spans them:
                {detail}. To compare hazards at all, scenarios are grouped into low, medium
                and high forcing tiers -- rcp26 with ssp126, rcp60 with ssp370, rcp85 with
                ssp585. Those pairings are approximate. They are built on different CMIP
                generations with different climate sensitivities and different socioeconomic
                assumptions, so a low-tier comparison across two hazards is a comparison of
                broadly similar futures, not of one scenario.
                """,
                evidence="values.csv scenario column; viz_common.SCENARIO_TIER",
                affects=sorted(rounds),
            )
        )
    return out


def layer_caveats(delivery: Delivery) -> List[dict]:
    """Per-layer caveats the processed files themselves declare."""
    out: List[dict] = []
    for _, row in delivery.layers.iterrows():
        lid = row["layer_id"]
        notes = [
            str(row.get(k) or "").strip()
            for k in ("interpretation_caveat", "delivery_note")
        ]
        note = " ".join(n for n in notes if n and n.lower() != "nan")
        if note:
            out.append(
                _cav(
                    f"LAYER-{lid.upper()}",
                    SEVERITY_SHOULD,
                    f"layer:{lid}",
                    f"How to read {lid} ({row.get('hazard', '')})",
                    note,
                    evidence=f"layers.csv, read from {row.get('source_folder', '')}",
                    affects=[lid],
                )
            )

        # A relative-baseline layer scores "unusual for this place", not "bad". Read as an
        # absolute risk it inverts the conclusion for exactly the sites a reader cares most
        # about, and it does so while every number is correct -- so it is MUST-DISCLOSE,
        # not should-note. Promoted 2026-08-13 after `cropfailure-3b` was found to rank
        # Iowa at the 99.3rd percentile of cropland against the Sahel's 69.4; the two
        # drought layers were already in this class and their own notes already said the
        # caveat "must be stated in any customer narrative" while the machinery filed it
        # as optional.
        if str(row.get("relative_baseline") or "").strip().lower() == "yes":
            note = str(row.get("relative_baseline_note") or "").strip()
            out.append(
                _cav(
                    f"RELATIVE-BASELINE-{lid.upper()}",
                    SEVERITY_MUST,
                    f"layer:{lid}",
                    f"{lid} scores departure from a local historical baseline, "
                    f"not absolute {str(row.get('hazard', '')).lower()} risk",
                    note,
                    evidence=f"layers.csv relative_baseline_note, "
                             f"read from {row.get('source_folder', '')}",
                    affects=[lid],
                )
            )

        models = str(row.get("impact_models") or "")
        n_models = len([m for m in models.split(",") if m.strip()])
        if n_models == 1:
            out.append(
                _cav(
                    f"SINGLE-MODEL-{lid.upper()}",
                    SEVERITY_SHOULD,
                    f"layer:{lid}",
                    f"{lid} has one impact model, so its confidence interval understates uncertainty",
                    f"""
                    {lid} is produced by a single impact model ({models}). Its confidence
                    interval therefore reflects spread between climate forcings only and
                    carries no structural model uncertainty at all. The interval is narrower
                    than the true uncertainty, and by an unknown amount -- there is no second
                    model to measure the gap against.
                    """,
                    evidence=f"layers.csv impact_models = {models}",
                    affects=[lid],
                )
            )

        smoothing = str(row.get("spatial_smoothing") or "")
        if smoothing and not smoothing.lower().startswith("none"):
            out.append(
                _cav(
                    f"SMOOTHING-{lid.upper()}",
                    SEVERITY_SHOULD,
                    f"layer:{lid}",
                    f"{lid} was spatially smoothed before extraction, so its footprint is wider than 1°",
                    f"""
                    This layer was smoothed at processing time as well as at extraction, so
                    the two stages compound and its effective footprint is materially wider
                    than the ~1° blend that applies to the other layers. Values for this
                    hazard should not be described as site-specific.
                    """,
                    evidence=f"layers.csv spatial_smoothing = {smoothing}",
                    affects=[lid],
                )
            )
    return out


def site_caveats(delivery: Delivery) -> List[dict]:
    out: List[dict] = []
    status_counts = delivery.manifest.get("counts", {}).get("data_status", {})

    vf = vulnerability_frame(delivery)
    unassessed_assets = sorted(
        set(vf.loc[vf["status"] == NOT_ASSESSED, "asset_id"])
    ) if not vf.empty else []
    if unassessed_assets:
        labels = [delivery.asset_label(a) for a in unassessed_assets]
        out.append(
            _cav(
                "SITES-NOT-ASSESSED",
                SEVERITY_MUST,
                "delivery",
                "Some assets returned no hazard values and are unassessed, not unexposed",
                f"""
                {len(unassessed_assets)} of {len(delivery.assets)} assets have no finite
                hazard value in at least one scenario-decade combination:
                {'; '.join(labels)}. A blank is not a zero. These assets sit outside the
                modelled land domain or outside a layer's own mask, and they must be
                excluded from any statement about the portfolio rather than counted as
                low-risk. Check first whether the coordinate is right -- a reversed
                latitude/longitude pair produces exactly this signature.
                """,
                evidence=f"values.csv data_status: {json.dumps(status_counts)}",
                affects=unassessed_assets,
            )
        )

    derived = delivery.locations[
        ~delivery.locations["coord_source"].astype(str).str.strip().isin(["supplied", ""])
    ]
    if len(derived):
        rows = [f"{r['name']} ({r['coord_source']})" for _, r in derived.iterrows()]
        out.append(
            _cav(
                "COORDS-DERIVED",
                SEVERITY_SHOULD,
                "delivery",
                "Some coordinates were derived rather than supplied by the customer",
                f"""
                {len(derived)} location(s) carry a coordinate we produced rather than
                received: {'; '.join(rows)}. Given the ~1° extraction footprint, a city
                centroid standing in for a large estate is a real approximation and not a
                rounding difference. Confirm these against the customer's own records.
                """,
                evidence="locations.csv coord_source",
                affects=list(derived["location_id"]),
            )
        )
    return out


def score_caveats(delivery: Delivery, taxonomy: dict) -> List[dict]:
    out: List[dict] = []
    if delivery.climate_score.empty:
        return out
    s = delivery.climate_score

    non_hazard = [l for l in delivery.layer_ids if l in (taxonomy.get("non_hazard_layers") or {})]
    if non_hazard:
        out.append(
            _cav(
                "SCORE-MIXED-AXES",
                SEVERITY_MUST,
                "delivery",
                "The Climate Score averages a hazard-exposure axis with an asset-condition axis",
                f"""
                The Climate Score is the unweighted mean of the risk percentiles of an
                asset's layers. For this portfolio that average includes
                {', '.join(non_hazard)}, which measures the asset's own productivity rather
                than a hazard acting on it. Its percentile is inverted at processing time so
                that low productivity reads as high risk, which puts it on the same 1-100
                axis -- but a productivity response and a cyclone exposure are different
                kinds of quantity, and averaging them is a presentational choice, not a
                measurement. Hazard counts elsewhere in this report exclude it; the Climate
                Score does not.
                """,
                evidence="config/hazard_taxonomy.yaml non_hazard_layers; manifest.json climate_score",
                affects=non_hazard,
            )
        )

    incomplete = int((s["n_hazards"] < s["n_hazards_expected"]).sum())
    if incomplete:
        out.append(
            _cav(
                "SCORE-INCOMPLETE-ROWS",
                SEVERITY_SHOULD,
                "delivery",
                "Some Climate Score rows average fewer layers than the asset carries",
                f"""
                {incomplete} of {len(s)} Climate Score rows were computed over fewer layers
                than the asset is mapped to, because a layer publishes no scenario in that
                forcing tier or no value in that decade. Those rows are not comparable with
                complete ones: an average over two layers and an average over four are
                different quantities with the same name. The n_hazards and
                n_hazards_expected columns in climate_score.csv distinguish them.
                """,
                evidence=f"manifest.json climate_score.incomplete_rows = {incomplete}",
            )
        )

    # Which assets fall out of a cross-tier comparison, and because of which layer.
    tiers_by_layer = {
        lid: {SCENARIO_TIER.get(sc) for sc in
              delivery.values[delivery.values["layer_id"] == lid]["scenario"].unique()}
        for lid in delivery.layer_ids
    }
    missing_tier = {
        lid: sorted(set(TIER_ORDER) - t) for lid, t in tiers_by_layer.items()
        if set(TIER_ORDER) - t
    }
    if missing_tier:
        dropped = sorted(
            {
                r["asset_id"]
                for _, r in delivery.assets.iterrows()
                for lid in str(r["layer_ids"]).split(";")
                if lid in missing_tier
            }
        )
        detail = "; ".join(
            f"{lid} publishes no {', '.join(TIER_LABELS[t].lower() for t in ts)} scenario"
            for lid, ts in sorted(missing_tier.items())
        )
        out.append(
            _cav(
                "PANEL-UNBALANCED",
                SEVERITY_SHOULD,
                "delivery",
                "Comparing forcing tiers drops assets whose hazards lack a scenario in every tier",
                f"""
                {detail}. Any chart or statement comparing one forcing tier against another
                is therefore restricted to the assets that are complete in all three tiers,
                which excludes {len(dropped)} of {len(delivery.assets)} assets
                ({', '.join(dropped)}). Without that restriction the high tier would be
                averaged over a different set of assets than the low tier, and the two would
                differ for reasons that have nothing to do with climate. Single-tier
                statements use the full portfolio.
                """,
                evidence="values.csv scenario coverage per layer",
                affects=dropped,
            )
        )
    return out


def config_caveats(delivery: Delivery) -> List[dict]:
    out: List[dict] = []
    cfg = delivery.config

    if str(cfg["horizons"].get("source", "")).lower() == "default":
        labels = ", ".join(
            f"{k} = {cfg['horizons'][k]['decade']}s" for k in ("short", "medium", "long")
            if k in cfg["horizons"]
        )
        out.append(
            _cav(
                "HORIZONS-DEFAULT",
                SEVERITY_MUST,
                "delivery",
                "The short, medium and long term horizons are ours, not the customer's",
                f"""
                This report uses {labels}. IFRS S2 requires the reporting entity to define
                these itself and to explain how they connect to its own planning and capital
                allocation horizons. These are placeholders chosen to span the available
                projection period. Replace them with the customer's real horizons -- for a
                timberland holding on a 25-year rotation, or an asset with a 10-year hold
                period, the correct buckets are materially different and would change which
                results lead.
                """,
                evidence="report_config.yaml horizons.source = default",
            )
        )

    out.append(
        _cav(
            "VULNERABILITY-METRIC-DEFERRED",
            SEVERITY_MUST,
            "delivery",
            "The headline disclosure metric — assets vulnerable to physical risk — is not reported",
            """
            IFRS S2 paragraph 29(c) and ESRS E1-9 both require the amount and percentage of
            assets vulnerable to physical climate risk. This assessment does NOT report it.
            What is measured here is exposure — where a site sits in the global distribution
            of a modelled hazard — and vulnerability is a statement about susceptibility to
            harm. Converting one into the other is a judgement about each asset type rather
            than a calculation, and it has not been made. Reporting a provisional figure
            would put a number nobody chose deliberately into a disclosure, so the section
            states the position instead. Exposure level, portfolio rank and direction of
            change are reported and are measured rather than inferred.
            """,
            evidence="report_config.yaml vulnerability block; see the report's section 6",
        )
    )

    values = pd.to_numeric(delivery.assets.get("asset_value"), errors="coerce") \
        if "asset_value" in delivery.assets.columns else pd.Series(dtype=float)
    n_valued = int(values.notna().sum())
    if n_valued < len(delivery.assets):
        out.append(
            _cav(
                "ASSET-VALUES-ABSENT",
                SEVERITY_MUST,
                "delivery",
                "Asset carrying amounts were not supplied, so the monetary disclosure is incomplete",
                f"""
                {len(delivery.assets) - n_valued} of {len(delivery.assets)} assets carry no
                value. IFRS S2 paragraph 29(c) and ESRS E1-9 both require a monetary AMOUNT,
                and no climate model can supply it. This is the SECOND obstacle to that
                metric rather than the only one: the method for determining which assets are
                vulnerable has also not been agreed, so this report reports no count,
                percentage or amount -- see VULNERABILITY-METRIC-DEFERRED. Supplying
                Asset_Value, Currency, Valuation_Date and Value_Basis with the site list
                removes this obstacle but not the other. The valuation basis matters as much
                as the figure, since book, insured, market and replacement values give four
                different answers to the same question.
                """,
                evidence=f"assets.csv asset_value populated for {n_valued} of {len(delivery.assets)}",
            )
        )
    return out


def build_caveats(delivery: Delivery) -> List[dict]:
    taxonomy = load_hazard_taxonomy()
    caveats: List[dict] = []
    caveats += hazard_coverage(delivery, taxonomy)
    caveats += qa_and_catalog(delivery)
    caveats += method_caveats(delivery)
    caveats += config_caveats(delivery)
    caveats += site_caveats(delivery)
    caveats += score_caveats(delivery, taxonomy)
    caveats += layer_caveats(delivery)

    order = {SEVERITY_MUST: 0, SEVERITY_SHOULD: 1, SEVERITY_INFO: 2}
    caveats.sort(key=lambda c: (order[c["severity"]], c["id"]))

    ids = [c["id"] for c in caveats]
    if len(ids) != len(set(ids)):
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        raise DeliveryError(
            f"Duplicate caveat id(s): {', '.join(dupes)}. Ids are cited from narratives and "
            f"must be unique within a delivery."
        )
    return caveats


# ---------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------

SEVERITY_TEXT = {
    SEVERITY_MUST: (
        "Must be disclosed in every report built from this delivery. The report build "
        "refuses to render without these and the verifier re-checks them afterwards."
    ),
    SEVERITY_SHOULD: (
        "Should be noted where relevant. These change how a specific number is read rather "
        "than what the delivery as a whole can claim."
    ),
    SEVERITY_INFO: "Context for whoever maintains the delivery.",
}


def write_markdown(delivery: Delivery, caveats: List[dict]) -> str:
    lines = [
        f"# Caveats and anomalies — {delivery.customer}",
        "",
        f"Delivery `{delivery.path.name}` · generated "
        f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "",
        "Derived mechanically from `manifest.json`, `layers.csv`, `values.csv`, "
        "`climate_score.csv`, `report_config.yaml` and `config/hazard_taxonomy.yaml`. "
        "This file is an **input** to the reports, not a summary of them: every "
        "`must_disclose` entry below is required to appear in each report, by ID.",
        "",
    ]
    for severity in (SEVERITY_MUST, SEVERITY_SHOULD, SEVERITY_INFO):
        group = [c for c in caveats if c["severity"] == severity]
        if not group:
            continue
        lines += [f"## {severity} ({len(group)})", "", SEVERITY_TEXT[severity], ""]
        for c in group:
            lines += [f"### `{c['id']}` — {c['title']}", "", c["text"], ""]
            if c["evidence"]:
                lines += [f"> **Evidence:** {c['evidence']}", ""]
            if c["affects"]:
                lines += [f"> **Affects:** {', '.join(c['affects'])}", ""]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("delivery", type=Path, help="deliveries/<customer>/<date>")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    try:
        delivery = load_delivery(args.delivery)
        caveats = build_caveats(delivery)
    except DeliveryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    counts = {
        s: sum(1 for c in caveats if c["severity"] == s)
        for s in (SEVERITY_MUST, SEVERITY_SHOULD, SEVERITY_INFO)
    }
    payload = {
        "customer": delivery.customer,
        "delivery": delivery.path.name,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
        "severity_meaning": SEVERITY_TEXT,
        "caveats": caveats,
    }
    (delivery.path / CAVEATS_JSON).write_text(json.dumps(payload, indent=2) + "\n")
    (delivery.path / CAVEATS_MD).write_text(write_markdown(delivery, caveats))

    record_stage(
        delivery.path,
        "caveats",
        "built",
        f"{len(caveats)} caveats ({counts[SEVERITY_MUST]} must-disclose)",
    )

    if not args.quiet:
        print(f"  caveats: {len(caveats)} "
              f"({counts[SEVERITY_MUST]} must-disclose, {counts[SEVERITY_SHOULD]} should-note)")
        for c in caveats:
            if c["severity"] == SEVERITY_MUST:
                print(f"    MUST  {c['id']:24} {c['title']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
