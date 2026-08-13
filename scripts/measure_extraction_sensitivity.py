#!/usr/bin/env python3
"""Reproduce the two extraction-footprint measurements that reports quote as fact.

    python scripts/measure_extraction_sensitivity.py

WHY THIS EXISTS
---------------
Two numbers appear in customer-facing reports as measured facts:

    "moving a site 0.25 deg changed a burnt-area value by 166%"
    "measured across 20,000 random sites, 100% resolve to a 4-cell blend"

Both were genuinely measured, and neither had a retained receipt: no script, no seed, no
result artifact. A number in a filing whose evidence exists only in someone's memory of a
session is indistinguishable, to a reviewer, from a number that was invented. An external
review flagged exactly that, and it was a fair hit.

So this script reproduces both on demand, deterministically. It reads the processed layers
directly and depends on no delivery. Run it whenever the extraction parameters change; if a
figure moves, update the caveat text and this docstring together.

The 4-cell result is a property of the SEARCH GEOMETRY (a 0.5 deg radius on a 0.5 deg grid,
applied per axis), not of any particular layer -- the 3x3 stencil requires a site to sit
exactly on a cell centre, and any drift collapses it to 2x2. The coordinate-sensitivity
figure IS layer- and site-dependent, which is the point of it: it says a coordinate error is
not a rounding error, and the size varies.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.utils.delivery import (  # noqa: E402
    EXTRACT_SEARCH_RADIUS,
    EXTRACT_SIGMA,
    DeliveryError,
    discover_scenarios,
    load_registry,
    scenario_path,
)
from scripts.utils.spatial_extract import extract_by_point  # noqa: E402

#: Fixed so the run is reproducible. Any seed gives the same answer for the blend count --
#: the result is geometric, not statistical -- but a stated seed is what makes it checkable.
SEED = 20260813
N_SITES = 20_000

#: Sites the sensitivity figure is quoted for. Kept here so the quoted number and the site
#: that produced it cannot drift apart.
SENSITIVITY_SITES = [
    ("Shasta Timberland", 40.8962, -121.7530),
    ("Mobile Distribution Hub", 30.6954, -88.0399),
    ("Bavaria Forest Block", 48.1370, 11.5760),
]
OFFSET_DEG = 0.25


def measure_blend(rng) -> dict:
    """How many grid cells does a random site's Gaussian window actually draw on?"""
    lats = np.arange(-89.75, 90.0, 0.5)
    lons = np.arange(-179.75, 180.0, 0.5)
    counts = np.zeros(N_SITES, dtype=int)
    site_lats = rng.uniform(-89.0, 89.0, N_SITES)
    site_lons = rng.uniform(-179.0, 179.0, N_SITES)
    for i, (la, lo) in enumerate(zip(site_lats, site_lons)):
        lat_d = np.abs(lats - la)
        lon_d = np.abs(lons - lo)
        lon_d = np.minimum(lon_d, 360.0 - lon_d)
        counts[i] = int((lat_d <= EXTRACT_SEARCH_RADIUS).sum()) * int(
            (lon_d <= EXTRACT_SEARCH_RADIUS).sum()
        )
    uniq, n = np.unique(counts, return_counts=True)
    return {int(u): int(c) for u, c in zip(uniq, n)}


def measure_sensitivity(registry, layer_id: str) -> list:
    spec = registry.get(layer_id)
    scenario = discover_scenarios(registry, spec)[-1]
    out = []
    with xr.open_dataset(scenario_path(registry, spec, scenario)) as ds:
        decade = int(ds.decade.values[-1])

        def at(la, lo):
            # extract_by_point returns {variable: {decade: value}} and needs the decade
            # dimension present, so it takes the whole dataset rather than a panel.
            got = extract_by_point(ds, la, lo, variables=["median"],
                                   sigma=EXTRACT_SIGMA,
                                   search_radius=EXTRACT_SEARCH_RADIUS)
            return got.get("median", {}).get(decade)

        for name, lat, lon in SENSITIVITY_SITES:
            b = at(lat, lon)
            if b is None or not np.isfinite(b) or b == 0:
                out.append((name, scenario, decade, b, None, None))
                continue
            worst = None
            for dlat, dlon in ((OFFSET_DEG, 0), (-OFFSET_DEG, 0), (0, OFFSET_DEG), (0, -OFFSET_DEG)):
                v = at(lat + dlat, lon + dlon)
                if v is None or not np.isfinite(v):
                    continue
                pct = abs(v - b) / abs(b) * 100.0
                if worst is None or pct > worst[0]:
                    worst = (pct, v, (dlat, dlon))
            out.append((name, scenario, decade, b, worst[0] if worst else None,
                        worst[2] if worst else None))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--layer", default="wildfire", help="layer to measure sensitivity on")
    args = ap.parse_args()

    print(f"Extraction parameters: sigma={EXTRACT_SIGMA} deg, "
          f"search_radius={EXTRACT_SEARCH_RADIUS} deg\n")

    rng = np.random.default_rng(SEED)
    blend = measure_blend(rng)
    total = sum(blend.values())
    print(f"1. CELL-BLEND COUNT  ({N_SITES:,} random sites, seed {SEED})")
    for cells, n in sorted(blend.items()):
        print(f"     {cells:2d} cells : {n:6,d}  ({100.0 * n / total:.2f}%)")
    print()

    try:
        registry = load_registry()
        rows = measure_sensitivity(registry, args.layer)
    except (DeliveryError, FileNotFoundError, OSError) as exc:
        print(f"2. COORDINATE SENSITIVITY -- skipped: {exc}")
        return 0

    print(f"2. COORDINATE SENSITIVITY  (layer {args.layer}, offset {OFFSET_DEG} deg)")
    for name, scenario, decade, base, pct, direction in rows:
        if pct is None:
            print(f"     {name:26} {scenario} {decade}s  base={base}  (no finite neighbour)")
        else:
            print(f"     {name:26} {scenario} {decade}s  base={base:.4g}  "
                  f"worst shift {pct:6.1f}%  at {direction}")
    print("\n   'Worst shift' is the largest change over the four cardinal offsets. This is a\n"
          "   per-site, per-layer figure: quote it as an illustration that a coordinate error\n"
          "   is not a rounding error, never as a universal sensitivity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
