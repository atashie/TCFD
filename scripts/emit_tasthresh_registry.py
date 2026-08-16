"""Emit the `layer_registry.yaml` block for the threshold-count ladder, READ FROM THE FILES.

PRINTS the YAML for review; it does not write the registry. The registry is a curated file
holding human decisions (status, delivery_note) alongside measured ones, and a script that
rewrote it in place would be one bad run away from clobbering the other fifteen layers.

Everything quoted here is read from the processed files' own global attributes, so the
registry cannot claim a `recommended_slope` the layer does not carry. Run after a build and
diff the output against what is already in the registry:

    .venv/bin/python3 scripts/emit_tasthresh_registry.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import xarray as xr
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from process_tasthresh import PROCESSED, RUNGS  # noqa: E402

#: Which rungs are cleared for delivery. hd35 and FD are the headline rungs recorded in
#: config/isimip_search_catalog.yaml search_results.chronic_heat.chosen_2026_08_14; the
#: other seven are built and contract-passing but not selected, exactly as `heatwave-2b` is.
#: `status: alternate` keeps them discoverable rather than invisible. No rung reaches a
#: delivery until an asset type in config/asset_catalog.yaml routes to it, which is a
#: separate decision needing user sign-off.
PREFERRED = {"hd35", "FD"}

LAYER_IDS = {
    "hd30": "heatdays-hd30", "hd35": "heatdays-hd35",
    "hd40": "heatdays-hd40", "hd45": "heatdays-hd45",
    "TR20": "tropicalnights-tr20", "TR25": "tropicalnights-tr25",
    "ID": "icedays-id", "FD": "frostdays-fd", "FDm10": "frostdays-fdm10",
}


def main():
    block, missing = {}, []
    for rung, (stem, hazard, measure, src, cmp_, thr) in RUNGS.items():
        folder = f"{stem}_{rung}_annual"
        files = sorted((PROCESSED / folder).glob(f"{rung}_*_processed.nc"))
        if not files:
            missing.append(rung)
            continue
        ds = xr.open_dataset(files[0])
        a = ds.attrs
        # Sanity: every scenario of a rung must agree about how to read itself. This is the
        # defect that forced the second rebuild, so it is asserted rather than assumed.
        slopes = set()
        for f in files:
            with xr.open_dataset(f) as d2:
                slopes.add(d2.attrs["recommended_slope"])
        if len(slopes) != 1:
            raise SystemExit(f"{rung}: scenarios disagree on recommended_slope: {slopes}")

        block[LAYER_IDS[rung]] = {
            "folder": folder,
            "file_prefix": rung,
            "hazard": hazard,
            "hazard_measure": measure,
            "status": "preferred" if rung in PREFERRED else "alternate",
            "qa_reviewed_on": None,
            "recommended_slope": a["recommended_slope"],
            "recommended_slope_rationale": a["recommended_slope_rationale"],
            # An ABSOLUTE threshold: 35 C is 35 C everywhere. This is the deliberate
            # complement to `heatwave-3b`, which IS relative_baseline: true.
            "relative_baseline": False,
            "delivery_note": delivery_note(rung, a, len(files)),
        }
        ds.close()

    if missing:
        print(f"# NOT YET BUILT, omitted: {', '.join(missing)}\n", file=sys.stderr)
    print(yaml.dump(block, sort_keys=False, width=96, default_flow_style=False))


def delivery_note(rung, a, n_scen):
    parts = [
        f"Absolute dry-bulb threshold count from ISIMIP3b bias-adjusted daily "
        f"{'tasmax' if rung.startswith(('hd', 'ID')) else 'tasmin'}, 12 GCMs x {n_scen} SSPs, "
        f"no impact model in the chain -- the confidence band carries GCM spread and "
        f"interannual variability only.",
        "THE ABSOLUTE COMPLEMENT TO `heatwave-3b`, NOT A SUBSTITUTE FOR IT. That layer "
        "scores departure from each cell's own preindustrial normal and ranks Chicago above "
        "Delhi; this one counts crossings of a fixed physical threshold. Measured over the "
        "65,797 shared land cells, the 2020s rank correlation is +0.554 and the two agree "
        "on only 47.2% of their worst-decile cells, so screening on one gives a materially "
        "different site list than screening on the other. Deliver them together and say "
        "which question each answers.",
    ]
    if a.get("sparsity_caveat", "").startswith("SPARSE"):
        parts.append("SPARSE: most sites read exactly 0 and tie at percentile tier 1. Use a "
                     "lower rung of the ladder for a measure that discriminates globally.")
    if a.get("saturation_caveat", "").startswith("MATERIAL"):
        parts.append("SATURATES at the calendar ceiling in the late century, where both "
                     "slope estimators go to zero and AGREE -- a flat trend there means no "
                     "headroom left, not stability.")
    if rung.startswith(("FD", "ID")):
        parts.append("DIRECTION: higher_is_worse on the frost/ice-day COUNT (user decision "
                     "2026-08-14), so this layer reports a hazard FALLING almost everywhere. "
                     "A declining percentile here is the physical result, not a data defect. "
                     "Note the opposite reading is also legitimate for some assets -- loss of "
                     "reliable frost harms crops needing vernalisation and winter roads "
                     "needing frozen ground -- and this layer does NOT represent that.")
    return " ".join(parts)


if __name__ == "__main__":
    sys.exit(main())
