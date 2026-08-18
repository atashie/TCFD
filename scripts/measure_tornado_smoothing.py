#!/usr/bin/env python3
"""Measure whether the CONUS tornado layers should be spatially smoothed, and at what
decay length. Reports a table; adopts nothing.

WHY THIS IS A MEASUREMENT AND NOT A SETTING
    CLAUDE.md: thin ensembles get spatial smoothing, but "the decay length is a
    per-layer measurement, not a constant, and so is whether to smooth at all". The
    tie-breaker there is a SPLIT-HALF TEST -- halve the ensemble, and if roughness
    barely moves while the halves still correlate, the roughness is real structure
    rather than sampling noise.

    This layer has no ensemble to halve, so the record is halved instead: ODD years
    against EVEN years. That split is deliberate. Splitting early-vs-late would
    confound the test with the reporting inhomogeneity documented in the layer's own
    `reporting_bias_caveat` (reports rose 2.8x since the 1950s), and the test would
    then measure the observing system rather than the sampling noise. Odd/even is
    balanced across the whole record, which removes the gross early-vs-late imbalance.
    It does NOT make the test immune to reporting bias: both halves carry the SAME
    spatial population, radar-coverage and rating-practice biases, and an abrupt
    observing-system change lands in both. It controls the secular trend, nothing more.

THE CRITERION -- AND THE ONE THAT DOES NOT WORK
    The obvious criterion is "smooth both halves, correlate them, take the sigma that
    maximises agreement". IT IS DEGENERATE, and measured to be so here before it was
    discarded: correlation rose MONOTONICALLY with sigma on all four rungs, out to the
    largest length tested (rho 0.81 -> 0.99 on `all`, 0.38 -> 0.98 on `f3plus`). Of
    course it did -- blurring any two fields toward a constant makes them agree. That
    criterion always returns the largest sigma offered, so it selects nothing. It is
    retained below only as a diagnostic, labelled, because the UNSMOOTHED value of it
    is informative even though its trend is not.

    The criterion actually used is HELD-OUT PREDICTIVE LIKELIHOOD, which is
    cross-validation and does have a genuine optimum:

        smooth the ODD half to get a rate field, then score how well that rate
        predicts the EVEN half's ACTUAL counts under a Poisson model, and vice versa.

    Too little smoothing predicts badly because the estimate is noise; too much
    predicts badly because real structure has been flattened. The optimum sits
    between, and its location is the decay length the data supports.

    A floor is applied to the predicted rate. Without one, an unsmoothed estimate that
    says "zero" in a cell where the held-out half recorded a tornado scores minus
    infinity, and sigma = 0 could never be compared with anything. The floor is
    reported so its effect is auditable rather than hidden.

    Note what this test can and cannot settle. It measures SAMPLING noise. It cannot
    see a bias shared by both halves -- if reporting density inflates a metro cell,
    it inflates it in odd and even years alike, and this test will call that
    agreement "real structure". It is a resolution question, not a validity one.

USAGE
    python3 scripts/measure_tornado_smoothing.py
    python3 scripts/measure_tornado_smoothing.py --rungs all f2plus
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.special import gammaln
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from process_tornado_spc import (  # noqa: E402
    LAT_MAX, LAT_MIN, LON_MAX, LON_MIN, RUNGS, SPC_CSV,
    build_conus_mask, build_grid, cell_area_km2, count_grid,
)

# Candidate decay lengths. 0 is the null (no smoothing) and must be in the list --
# a sweep that omits it can only ever recommend smoothing.
SIGMAS_KM = [0, 10, 15, 25, 40, 60, 90, 130, 180, 250]


def nan_aware_gaussian(field: np.ndarray, valid: np.ndarray,
                       sigma_y: float, sigma_x: float) -> np.ndarray:
    """Normalised convolution: smooth only with the weight that actually landed on
    valid cells, so the coast does not bleed zeros inland."""
    # Validity must include finiteness of the FIELD, not just the mask. A NaN sitting
    # inside the mask would otherwise be handed to gaussian_filter and contaminate its
    # whole neighbourhood, and the function's name would be a lie. (Current rate arrays
    # are finite in-mask, so this is defensive rather than a live fix.)
    valid = valid & np.isfinite(field)
    if sigma_y <= 0 and sigma_x <= 0:
        out = field.copy()
        out[~valid] = np.nan
        return out
    filled = np.where(valid, field, 0.0)
    num = gaussian_filter(filled, sigma=(sigma_y, sigma_x), mode="constant", cval=0.0)
    den = gaussian_filter(valid.astype(float), sigma=(sigma_y, sigma_x),
                          mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 1e-9, num / den, np.nan)
    out[~valid] = np.nan
    return out


def roughness(field: np.ndarray, valid: np.ndarray) -> float:
    """Mean |adjacent difference| relative to the mean level. Scale-free, so it is
    comparable across rungs whose absolute rates differ by an order of magnitude."""
    d = []
    for ax in (0, 1):
        a = np.diff(field, axis=ax)
        v = np.diff(valid.astype(int), axis=ax) == 0
        m = valid.take(np.arange(field.shape[ax] - 1), axis=ax) & v
        d.append(np.abs(a[m & np.isfinite(a)]))
    diffs = np.concatenate(d) if d else np.array([0.0])
    level = np.nanmean(field[valid])
    return float(np.mean(diffs) / level) if level else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolution", type=float, default=0.25)
    ap.add_argument("--rungs", nargs="*", default=list(RUNGS))
    args = ap.parse_args()

    res = args.resolution
    df = pd.read_csv(SPC_CSV, low_memory=False)
    df = df[(df.slat > LAT_MIN) & (df.slat < LAT_MAX)
            & (df.slon > LON_MIN) & (df.slon < LON_MAX)]

    lat_edges, lon_edges, lat_c, lon_c = build_grid(res)
    nlat, nlon = len(lat_c), len(lon_c)
    area = cell_area_km2(lat_edges, lon_edges)
    mask = build_conus_mask(lat_c, lon_c)

    km_per_cell_y = res * 111.0
    km_per_cell_x = res * 111.0 * float(np.cos(np.deg2rad(np.mean(lat_c))))
    print(f"Split-half smoothing test -- {res} deg, {int(mask.sum()):,} CONUS cells")
    print(f"cell ~{km_per_cell_y:.1f} km N-S, ~{km_per_cell_x:.1f} km E-W at mean latitude")
    print("Halves: ODD years vs EVEN years -- controls the secular reporting trend;\n       shared spatial bias is present in BOTH halves and is invisible here\n")

    odd = df[df.yr % 2 == 1]
    even = df[df.yr % 2 == 0]

    for rung in args.rungs:
        threshold = RUNGS[rung][0]
        sub_o = odd if threshold is None else odd[odd.mag >= threshold]
        sub_e = even if threshold is None else even[even.mag >= threshold]

        # EXPOSURE IS CALENDAR TIME, NOT OBSERVED TIME. `sub.yr.nunique()` counts only
        # years containing a qualifying report, so a rung with a nationally event-free
        # year would silently shrink its own denominator and inflate the rate. Verified
        # 2026-08-18 that all four rungs currently have 38/38 in both halves, i.e. this
        # was inert -- but it is wrong in principle and would bite a thinner subset.
        yrs = range(int(df.yr.min()), int(df.yr.max()) + 1)
        n_yr_o = sum(1 for y in yrs if y % 2 == 1)
        n_yr_e = sum(1 for y in yrs if y % 2 == 0)
        ca, _ = count_grid(sub_o, res, nlat, nlon, "track")
        cb, _ = count_grid(sub_e, res, nlat, nlon, "track")
        ra = ca / ((area / 1e4) * n_yr_o)
        rb = cb / ((area / 1e4) * n_yr_e)

        exp_o = (area / 1e4) * n_yr_o
        exp_e = (area / 1e4) * n_yr_e
        floor = float(np.nanmean(ra[mask])) / 1000.0

        print(f"[{rung}] odd n={len(sub_o):,} ({n_yr_o} yr)  even n={len(sub_e):,} ({n_yr_e} yr)"
              f"   rate floor {floor:.5f}")
        print(f"    {'sigma_km':>9} {'cells':>7} {'NB logL/cell':>14} "
              f"{'plugin logL/cell':>17} {'rho(smoothed)':>14} {'roughness':>10}")
        best = (-np.inf, None)
        best_plugin = (-np.inf, None)
        for s_km in SIGMAS_KM:
            sy = s_km / km_per_cell_y
            sx = s_km / km_per_cell_x
            fa = nan_aware_gaussian(ra, mask, sy, sx)
            fb = nan_aware_gaussian(rb, mask, sy, sx)

            # Cross-validated Poisson predictive log-likelihood, both directions.
            ll = 0.0
            n = 0
            for pred, cnt, expo in ((fa, cb, exp_e), (fb, ca, exp_o)):
                ok = mask & np.isfinite(pred)
                lam = np.clip(pred[ok], floor, None) * expo[ok]
                k = cnt[ok]
                ll += float(np.sum(k * np.log(lam) - lam))
                n += int(ok.sum())
            ll_per_cell = ll / n

            # FLOOR-FREE ALTERNATIVE: smooth the training COUNTS and the training
            # EXPOSURE, giving a Gamma(alpha, beta) posterior per cell, and integrate it
            # out to a negative-binomial predictive for the held-out count. No pseudocount
            # is needed, because the posterior never assigns zero probability to an event.
            # This is the criterion to trust; the plug-in column is kept for comparison.
            nb = 0.0
            nb_n = 0
            for ctr, cte, T_tr, T_te in ((ca, cb, exp_o, exp_e), (cb, ca, exp_e, exp_o)):
                a_s = nan_aware_gaussian(ctr.astype(float), mask, sy, sx)
                b_s = nan_aware_gaussian(T_tr, mask, sy, sx)
                okk = mask & np.isfinite(a_s) & np.isfinite(b_s)
                alpha = a_s[okk] + 0.5
                beta = b_s[okk]
                kk = cte[okk]
                pr = beta / (beta + T_te[okk])
                nb += float(np.sum(gammaln(kk + alpha) - gammaln(alpha) - gammaln(kk + 1)
                                   + alpha * np.log(pr) + kk * np.log1p(-pr)))
                nb_n += int(okk.sum())
            nb_per_cell = nb / nb_n

            ok = mask & np.isfinite(fa) & np.isfinite(fb)
            rho = float(spearmanr(fa[ok], fb[ok]).statistic)
            rgh = roughness(fa, mask)
            if nb_per_cell > best[0]:
                best = (nb_per_cell, s_km)
            if ll_per_cell > best_plugin[0]:
                best_plugin = (ll_per_cell, s_km)
            print(f"    {s_km:>9} {sy:>7.2f} {nb_per_cell:>14.4f} {ll_per_cell:>17.4f} "
                  f"{rho:>14.4f} {rgh:>10.3f}")
        print(f"    -> NB predictive (floor-free, PRIMARY) peaks at sigma = {best[1]} km")
        print(f"       plug-in Poisson peaks at sigma = {best_plugin[1]} km -- the plug-in "
              f"OVER-SMOOTHS because it ignores estimation uncertainty")
        print(f"       [rho column is MONOTONE BY CONSTRUCTION -- read its sigma=0 value, "
              f"not its trend]\n")

    print("READ THIS BEFORE ADOPTING ANYTHING:")
    print("  * A peak at sigma = 0 means the field is already resolved and smoothing")
    print("    would only cost resolution.")
    print("  * A peak at sigma > 0 gives the decay length the DATA supports, for THIS")
    print("    rung. It is not transferable to another rung, another window, or another")
    print("    layer -- each is its own measurement.")
    print("  * This test sees sampling noise only. A bias present in both halves (metro")
    print("    reporting density) reads as 'real structure' here. It is not validity.")
    print("  * Nothing was written. No layer was modified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
