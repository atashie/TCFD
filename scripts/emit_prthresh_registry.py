"""Emit the `layer_registry.yaml` block for the precipitation layers, READ FROM THE FILES.

PRINTS the YAML for review; it does not write the registry. That file carries human
decisions (status, delivery_note, qa_reviewed_on) alongside measured ones, and a script that
rewrote it in place would be one bad run away from clobbering the other twenty-three layers.

Everything quoted is read from the processed files' own global attributes, so the registry
cannot claim a `recommended_slope` or a `relative_baseline` the layer does not carry.

    .venv/bin/python3 scripts/emit_prthresh_registry.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import xarray as xr
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from process_prthresh import METRICS, PROCESSED  # noqa: E402

#: THREE headline metrics, one per question, because this family asks three different ones
#: and the user chose to carry both an absolute and a relative framing (2026-08-17).
#:   Rx1day  intensity  -- the best single pluvial proxy, and needs no threshold at all
#:   R20mm   absolute   -- sits at the MEDIAN cell's own wet-day p95 (19.0 mm measured)
#:   R95pD   relative   -- "extreme for here", which is what local drainage was built to
#: The other seven are built, contract-passing and `status: alternate`: discoverable, not
#: chosen by default. Routing any of them is a separate decision needing user sign-off.
PREFERRED = {"Rx1day", "R20mm", "R95pD"}

LAYER_IDS = {m: f"pluvial-{m.lower()}" for m in METRICS}


def delivery_note(metric, a, relative, n_scen):
    parts = [
        f"Constructed from ISIMIP3b bias-adjusted daily `pr`, 14 GCMs x {n_scen} SSPs, with "
        f"NO impact model in the chain -- the confidence band carries GCM spread and "
        f"interannual variability only.",
        "A RAINFALL HAZARD, NEVER A FLOOD OUTCOME. Pluvial flooding is rainfall against "
        "DRAINAGE CAPACITY and drainage is not in this data, so a well-drained urban cell "
        "and an undrained rural one with the same rainfall read identically. Distinct from "
        "the `flood-3b-*` CaMa-Flood layers, which are RIVER flooding routed through a "
        "hydrological model; a site can be exposed to either without the other.",
        "THIS ENSEMBLE IS NOT THE HEAT LADDER'S: CESM2-WACCM and IITM-ESM publish `pr` but "
        "no `tasmax`, so a site's rainfall and temperature numbers come from different model "
        "sets. Say so if a narrative combines them.",
    ]
    if relative:
        parts.append(
            "RELATIVE METRIC -- see the must-disclose RELATIVE-BASELINE caveat. The threshold "
            "differs in every cell (3.3-165.2 mm/day, median 19.0), the 2020s panel is a "
            "control by construction, and the ranking INVERTS against the absolute metrics: "
            "Singapore 14.9 d/yr against Cherrapunji's 6.8, while on R20mm Cherrapunji reads "
            "100.8 against Singapore's 23.3. Deliver it beside R20mm or Rx1day, never alone.")
    else:
        parts.append(
            "ABSOLUTE metric: the threshold means the same thing everywhere, which is what "
            "makes it comparable across a portfolio and what makes it near-silent across the "
            "arid belt. Pair with R95pD where 'extreme for this place' is the question.")
    if a.get("sparsity_caveat", "").startswith("SPARSE"):
        parts.append("SPARSE: most sites read exactly 0 and tie at percentile tier 1. Use a "
                     "lower depth or Rx1day for a measure that discriminates globally.")
    if metric in ("R95pD", "R99pD"):
        parts.append("MASKED in 599 land cells (0.91%) where the pooled 2020s baseline holds "
                     "fewer than 840 wet days and no percentile can be defined. Those cells "
                     "are NaN, not 0 -- the absolute metrics still cover them.")
    if metric == "Rx5day":
        parts.append("Rx5day windows CROSS chunk boundaries and are assigned to the year they "
                     "END in; verified against Rx1day (a 5-day window contains its own "
                     "wettest day) with 0 violations in every scenario.")
    return " ".join(parts)


def main():
    block, missing = {}, []
    for metric, (stem, hazard, measure, units, central, relative) in METRICS.items():
        folder = f"{stem}_{metric}_annual"
        files = sorted((PROCESSED / folder).glob(f"{metric}_*_processed.nc"))
        if not files:
            missing.append(metric)
            continue
        slopes, rels = set(), set()
        for f in files:
            with xr.open_dataset(f) as d:
                slopes.add(d.attrs["recommended_slope"])
                rels.add("relative_baseline_note" in d.attrs)
        if len(slopes) != 1:
            raise SystemExit(f"{metric}: scenarios disagree on recommended_slope: {slopes}")
        if len(rels) != 1:
            raise SystemExit(f"{metric}: scenarios disagree on relative_baseline: {rels}")
        ds = xr.open_dataset(files[0])
        a = ds.attrs
        entry = {
            "folder": folder,
            "file_prefix": metric,
            "hazard": hazard,
            "hazard_measure": measure,
            "status": "preferred" if metric in PREFERRED else "alternate",
            "qa_reviewed_on": None,
            "recommended_slope": a["recommended_slope"],
            "recommended_slope_rationale": a["recommended_slope_rationale"],
            "relative_baseline": bool(relative),
        }
        if relative:
            entry["relative_baseline_note"] = a["relative_baseline_note"]
        entry["delivery_note"] = delivery_note(metric, a, relative, len(files))
        block[LAYER_IDS[metric]] = entry
        ds.close()

    if missing:
        print(f"# NOT BUILT, omitted: {', '.join(missing)}\n", file=sys.stderr)
    print(yaml.dump(block, sort_keys=False, width=96, default_flow_style=False))


if __name__ == "__main__":
    sys.exit(main())
