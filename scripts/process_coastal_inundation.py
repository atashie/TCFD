"""Stage 3 of the coastal layer: OUTPUT-SPEC panels from hypsometry x sea-level delta.

Slides each 0.5 degree cell's coastal-land elevation histogram (Stage 2) against that
cell's per-member annual sea-level delta (Stage 1) to produce the contract panels.

WHAT THE LAYER REPORTS, AND WHY IT IS A CHANGE RATHER THAN A COUNT
------------------------------------------------------------------
`median` is the AREA FRACTION of a cell's ocean-connected coastal land lying at or below
projected sea level -- the same dimensionless [0,1] exposure share as `floodedarea` and
`driedarea`, so it slots into the existing contract without a new value class.

It is deliberately not an absolute inundation count. GEBCO is a digital SURFACE model and
in coastal zones SRTM-class surface models run +2.49 m (Australia) to +3.67 m (US) high,
while the sea-level signal here is +0.21 m (rcp26) to +0.34 m (rcp60) by the 2090s. The
DEM's systematic error is an order of magnitude larger than the climate signal, so an
absolute threshold would publish the DEM's bias. Evaluating the SAME histogram at the
baseline and at each future decade puts the bias on both sides: it still decides where on
the hypsometric curve the calculation sits, but it no longer decides the answer.
Correcting that bias is known to TRIPLE exposure estimates (Kulp & Strauss 2019), which is
the size of the effect being sidestepped, not eliminated.

THE PERCENTILE IS A DECLARED DEVIATION FROM OUTPUT-SPEC
-------------------------------------------------------
The contract defines `percentile` as percentile-of-score against the shared 2020s baseline
distribution -- a RELATIVE rank. This layer instead uses an ABSOLUTE calibration (user
decision 2026-08-14):

    elevation above projected sea level <= 0 m   -> 100
    elevation above projected sea level >= 10 m  ->   1
    in between                                   -> 2..99, filled by the empirical
                                                    distribution of coastal elevations

10 m is the Low Elevation Coastal Zone bound of McGranahan, Balk & Anderson (2007) -- "the
contiguous area along the coast that is less than 10 metres above sea level", defined as
hydrologically connected to the sea, which is also why Stage 2 computes connectivity.

Two consequences worth stating plainly. First, this percentile means something absolute: 0
m reads 100 regardless of what any other cell does, so unlike `drought-*` and
`cropfailure-3b` this layer is NOT scoring departure from a local baseline and carries no
relative_baseline caveat. Second, it is therefore NOT on the same axis as the other
layers' percentiles, and the Climate Score averages percentiles across layers -- that has
to be stated wherever the score is reported.

The per-cell percentile is computed per sub-cell and then area-weighted within the 0.5
degree cell, rather than by calibrating the cell's mean elevation. A cell that is half at
2 m and half at 40 m is not the same risk as a cell uniformly at 21 m, and averaging the
elevation first would say it was.

Run:   .venv/bin/python3 scripts/process_coastal_inundation.py --self-test
       .venv/bin/python3 scripts/process_coastal_inundation.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

HYPSO = Path("data/interim/coastal/hypsometry_halfdeg.nc")
DELTA = Path("data/interim/slr_delta")
OUT = Path("data/processed/coastal-inundation")

LECZ_M = 10.0          # McGranahan et al. 2007
P_AT_SEA_LEVEL = 100
P_AT_LECZ = 1
P_HI, P_LO = 99, 2     # the band filled by the elevation distribution


def calibration_curve(ref_hist: np.ndarray, edges: np.ndarray):
    """Map elevation-above-projected-sea-level -> percentile, from a fixed reference.

    Built ONCE from the 2020s baseline distribution and reused for every decade and
    scenario. Rebuilding it per decade would renormalise the axis each time, so a cell
    whose elevation never changed would drift in percentile as the rest of the world moved
    -- which is the opposite of what an absolute calibration is for.
    """
    centres = 0.5 * (edges[:-1] + edges[1:])
    band = (centres > 0) & (centres < LECZ_M)
    w = ref_hist[band].astype("float64")
    if w.sum() == 0:
        raise ValueError("reference distribution is empty inside (0, LECZ)")
    cdf = np.cumsum(w) / w.sum()
    # p falls from P_HI just above sea level to P_LO just below the LECZ bound.
    p_band = P_HI - (P_HI - P_LO) * cdf

    def to_percentile(e: np.ndarray) -> np.ndarray:
        out = np.empty_like(e, dtype="float64")
        out[:] = np.interp(e, centres[band], p_band,
                           left=P_HI, right=P_LO)
        out[e <= 0] = P_AT_SEA_LEVEL
        out[e >= LECZ_M] = P_AT_LECZ
        return out

    return to_percentile, centres


def cell_panels(hist, n_below, n_above, centres, water_levels, to_percentile):
    """Exposure fraction and area-weighted percentile for one cell at many water levels.

    hist          (n_bins,)  counts by elevation bin
    water_levels  (k,)       sea-level delta in m, one per (member, year)
    returns       (k,) exposed fraction, (k,) area-weighted percentile
    """
    total = hist.sum() + n_below + n_above
    if total == 0:
        return (np.full(water_levels.shape, np.nan),
                np.full(water_levels.shape, np.nan))

    # Cumulative counts, so P(elev <= w) is a lookup rather than a per-level scan.
    cum = np.concatenate([[0.0], np.cumsum(hist)])
    idx = np.clip(np.searchsorted(centres, water_levels, side="right"), 0, len(hist))
    exposed = (cum[idx] + n_below) / total

    # Area-weighted mean percentile: sum over bins of count * p(centre - w).
    e_rel = centres[None, :] - water_levels[:, None]        # (k, n_bins)
    p = to_percentile(e_rel.ravel()).reshape(e_rel.shape)
    num = (p * hist[None, :]).sum(axis=1)
    # Land below the histogram floor is far under water -> 100; land above the ceiling is
    # far above the LECZ bound -> 1. Neither is dropped from the denominator.
    num += n_below * P_AT_SEA_LEVEL + n_above * P_AT_LECZ
    return exposed, num / total


def self_test():
    """Verify the math on histograms whose answers can be worked out by hand."""
    edges = np.arange(-30.0, 30.0 + 1e-9, 0.05)
    centres = 0.5 * (edges[:-1] + edges[1:])
    n = len(centres)

    # Reference: elevations uniform on (0, 10] -> the calibration is LINEAR in elevation.
    ref = np.zeros(n)
    ref[(centres > 0) & (centres < LECZ_M)] = 1.0
    to_p, _ = calibration_curve(ref, edges)

    print("calibration on a uniform reference (should be linear 99 -> 2):")
    for e in (-1.0, 0.0, 0.01, 2.5, 5.0, 7.5, 9.99, 10.0, 25.0):
        print(f"  e={e:6.2f} m -> p={float(to_p(np.array([e]))[0]):7.3f}")
    mid = float(to_p(np.array([5.0]))[0])
    assert abs(mid - 50.5) < 1.5, f"midpoint should be ~50.5, got {mid}"
    assert float(to_p(np.array([-0.5]))[0]) == 100
    assert float(to_p(np.array([12.0]))[0]) == 1

    # A cell with land uniform on (0, 20]: half of it inside the LECZ.
    hist = np.zeros(n)
    hist[(centres > 0) & (centres < 20.0)] = 1.0
    total = hist.sum()

    w = np.array([0.0, 0.5, 1.0])
    exp, pct = cell_panels(hist, 0, 0, centres, w, to_p)
    print("\ncell uniform on (0,20], exposure at rising sea level:")
    for k, ww in enumerate(w):
        print(f"  w={ww:.2f} m -> exposed={exp[k]:.4f} (expect {ww/20:.4f})  "
              f"percentile={pct[k]:.2f}")
    assert abs(exp[0] - 0.0) < 1e-6 and abs(exp[1] - 0.025) < 2e-3, exp
    assert exp[2] > exp[1] > exp[0], "exposure must rise with sea level"
    assert pct[2] > pct[1] > pct[0], "percentile must rise with sea level"

    # Sanity: a cell entirely above the LECZ never moves.
    high = np.zeros(n)
    high[centres > 20.0] = 1.0
    e2, p2 = cell_panels(high, 0, 0, centres, w, to_p)
    print(f"\ncell entirely above 20 m: exposed={e2}  percentile={p2}")
    assert np.allclose(e2, 0) and np.allclose(p2, 1), "high ground must stay at 1"

    # Sanity: a cell entirely below present sea level is fully exposed and pinned at 100.
    low = np.zeros(n)
    low[centres < -1.0] = 1.0
    e3, p3 = cell_panels(low, 0, 0, centres, w, to_p)
    print(f"cell entirely below -1 m: exposed={e3}  percentile={p3}")
    assert np.allclose(e3, 1) and np.allclose(p3, 100), "drowned land must pin at 100"

    # The denominator must include the tails, or exposure can exceed 1.
    e4, p4 = cell_panels(hist, n_below=int(total), n_above=int(total), centres=centres,
                         water_levels=np.array([0.0]), to_percentile=to_p)
    assert abs(e4[0] - 1 / 3) < 1e-6, f"tails must count in the denominator, got {e4[0]}"
    print(f"\ntail accounting: exposed={e4[0]:.4f} (expect 0.3333)  percentile={p4[0]:.2f}")
    print("\nself-test PASSED")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    if not HYPSO.exists():
        raise SystemExit(f"{HYPSO} not built yet -- run build_coastal_hypsometry.py first")
    raise SystemExit("full run not wired yet; Stage 2 output is required first")


if __name__ == "__main__":
    main()
