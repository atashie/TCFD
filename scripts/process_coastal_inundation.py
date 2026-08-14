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

import sys

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))

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


def cell_panels(hist, n_below, n_safe, hist_expo, n_below_expo, n_land,
                centres, water_levels, to_percentile):
    """Exposure fraction and area-weighted percentile for one cell at many water levels.

    TWO histograms, because the two questions need different connectivity. `hist_expo`
    holds land the sea can actually reach (connected at 2 m) and answers "what floods".
    `hist` holds the LECZ pool (connected at 10 m) and answers "how close to sea level is
    this place", which has to keep grading coastal plain at 5 m instead of scoring it 1.
    Both are shares of `n_land`, the cell's full coastal land, so neither is a fraction of
    a subset the method chose for itself.
    """
    if n_land == 0:
        nan = np.full(water_levels.shape, np.nan)
        return nan, nan.copy()

    # Cumulative counts, so P(elev <= w) is a lookup rather than a per-level scan.
    cum = np.concatenate([[0.0], np.cumsum(hist_expo)])
    idx = np.clip(np.searchsorted(centres, water_levels, side="right"), 0, len(hist_expo))
    exposed = (cum[idx] + n_below_expo) / n_land

    # Area-weighted mean percentile: sum over bins of count * p(centre - w).
    e_rel = centres[None, :] - water_levels[:, None]        # (k, n_bins)
    p = to_percentile(e_rel.ravel()).reshape(e_rel.shape)
    num = (p * hist[None, :]).sum(axis=1)
    # Pool land below the histogram floor is far under water -> 100. Everything the sea
    # cannot reach -- above the LECZ bound, or below it inside a closed basin -> 1.
    num += n_below * P_AT_SEA_LEVEL + n_safe * P_AT_LECZ
    return exposed, num / n_land


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
    exp, pct = cell_panels(hist, 0, 0, hist, 0, total, centres, w, to_p)
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
    e2, p2 = cell_panels(high, 0, 0, high, 0, high.sum(), centres, w, to_p)
    print(f"\ncell entirely above 20 m: exposed={e2}  percentile={p2}")
    assert np.allclose(e2, 0) and np.allclose(p2, 1), "high ground must stay at 1"

    # Sanity: a cell entirely below present sea level is fully exposed and pinned at 100.
    low = np.zeros(n)
    low[centres < -1.0] = 1.0
    e3, p3 = cell_panels(low, 0, 0, low, 0, low.sum(), centres, w, to_p)
    print(f"cell entirely below -1 m: exposed={e3}  percentile={p3}")
    assert np.allclose(e3, 1) and np.allclose(p3, 100), "drowned land must pin at 100"

    # The denominator must include the tails, or exposure can exceed 1.
    e4, p4 = cell_panels(hist, int(total), int(total), hist, int(total), int(3 * total),
                         centres, np.array([0.0]), to_p)
    assert abs(e4[0] - 1 / 3) < 1e-6, f"tails must count in the denominator, got {e4[0]}"
    print(f"\ntail accounting: exposed={e4[0]:.4f} (expect 0.3333)  percentile={p4[0]:.2f}")

    # Closed-basin case: land is in the LECZ pool but the sea cannot reach it, so it must
    # be graded on elevation yet never exposed. This is the Salton Sink.
    sink = np.zeros(n); sink[centres < -1.0] = 1.0
    e5, p5 = cell_panels(sink, 0, 0, np.zeros(n), 0, int(sink.sum()), centres, w, to_p)
    print(f"closed basin below sea level: exposed={e5}  percentile={p5}")
    assert np.allclose(e5, 0), "unreachable land must never be exposed"
    assert np.allclose(p5, 100), "but it is still below sea level, so it grades 100"
    print("\nself-test PASSED")


SCENARIOS = ["rcp26", "rcp60"]
NLAT, NLON = 360, 720
DECADES = list(range(2010, 2100, 10))
BASELINE_DECADE = 2020
VARIABLE = "coastalinundation"


def build_annual(scen: str):
    """(exposed, pct) each (member, year, cell), plus the coords needed to write out."""
    hyp = xr.open_dataset(HYPSO)
    cell_ix = hyp["cell_flat_index"].values
    hist = hyp["hist"].values.astype("float64")
    hist_expo = hyp["hist_expo"].values.astype("float64")
    n_below = hyp["n_below"].values.astype("float64")
    n_below_expo = hyp["n_below_expo"].values.astype("float64")
    n_safe = hyp["n_above"].values.astype("float64")
    n_land = hyp["n_land"].values.astype("float64")
    bin_left = hyp["bin_left"].values
    bw = float(bin_left[1] - bin_left[0])
    centres = bin_left + bw / 2
    edges = np.append(bin_left, bin_left[-1] + bw)

    # The calibration reference is the POOLED baseline distribution over every coastal
    # cell -- one global axis, fixed once, so a percentile means the same thing in every
    # decade and every scenario.
    to_p, _ = calibration_curve(hist.sum(axis=0), edges)

    d = xr.open_dataset(DELTA / f"slr_delta_{scen}.nc")
    land_ix = d["land_flat_index"].values
    years = d["year"].values
    n_mem = d.sizes["member"]

    # Join hypsometry cells to sea-level cells on the shared flat index.
    pos = -np.ones(NLAT * NLON, dtype="int64")
    pos[land_ix] = np.arange(land_ix.size)
    sel = pos[cell_ix]
    keep = sel >= 0
    print(f"[{scen}] hypsometry cells {cell_ix.size:,}; matched to a sea level "
          f"{int(keep.sum()):,}")

    delta = d["slr_delta"].values[:, :, sel[keep]]        # (member, year, cell)
    hist, hist_expo = hist[keep], hist_expo[keep]
    n_below, n_below_expo = n_below[keep], n_below_expo[keep]
    n_safe, n_land = n_safe[keep], n_land[keep]
    out_ix = cell_ix[keep]
    n_cells = out_ix.size

    exposed = np.full((n_mem, len(years), n_cells), np.nan, "float32")
    pct = np.full((n_mem, len(years), n_cells), np.nan, "float32")
    for c in range(n_cells):
        w = delta[:, :, c].ravel()
        ok = np.isfinite(w)
        if not ok.any():
            continue
        e, p = cell_panels(hist[c], n_below[c], n_safe[c], hist_expo[c],
                           n_below_expo[c], n_land[c], centres, w[ok], to_p)
        eo = np.full(w.shape, np.nan); eo[ok] = e
        po = np.full(w.shape, np.nan); po[ok] = p
        exposed[:, :, c] = eo.reshape(n_mem, len(years))
        pct[:, :, c] = po.reshape(n_mem, len(years))
        if c % 2000 == 0:
            print(f"  cell {c:,}/{n_cells:,}", end="\r")
    print()
    hyp.close(); d.close()
    return exposed, pct, years, out_ix, [str(m) for m in
                                         xr.open_dataset(DELTA / f"slr_delta_{scen}.nc")
                                         ["member"].values]


def measure():
    """Zero-inflation of the exposure field -- decides the decadal-statistic branch.

    CLAUDE.md requires the third branch (pooled_mean_zero_inflated) to be adopted only on
    a measurement, never to improve contrast, so this reports what the median branch would
    publish before anything is written.
    """
    for scen in SCENARIOS:
        exposed, pct, years, out_ix, mem = build_annual(scen)
        fin = np.isfinite(exposed)
        z = (exposed == 0) & fin
        print(f"\n[{scen}] exposure field")
        print(f"  finite cell-member-years : {int(fin.sum()):,}")
        print(f"  exactly zero             : {100 * z.sum() / fin.sum():.2f}%")
        base = (years >= 2020) & (years < 2030)
        panel_med = np.nanmedian(exposed[:, base], axis=(0, 1))
        panel_mean = np.nanmean(exposed[:, base], axis=(0, 1))
        print(f"  2020s panel, MEDIAN branch: {int((panel_med > 0).sum()):,} exposed cells "
              f"of {out_ix.size:,}")
        print(f"  2020s panel, MEAN   branch: {int((panel_mean > 0).sum()):,} exposed cells")
        f9 = years >= 2090
        print(f"  2090s mean exposure       : {np.nanmean(exposed[:, f9]):.5f}")
        print(f"  2090s mean percentile     : {np.nanmean(pct[:, f9]):.2f}")


def run():
    """Write the OUTPUT-SPEC panels, one file per scenario.

    DECADAL STATISTIC: median + IQR, the contract default (user decision 2026-08-14). The
    third branch was measured and declined: 55.3% of the exposure field is exactly zero and
    the median branch publishes 3,493 exposed cells against the mean branch's 10,620, a 67%
    erasure. Real, but well short of the two layers that took the deviation -- `let` at
    97.8% zero / 93% erased and `cropfailure-3b` at 96.1% / 96.6% -- and above `burntarea`,
    which took the median branch at 29.2%. Recorded so the decision is auditable, not
    because the count is a precedent.
    """
    from utils.decadal_stats import pooled_decadal_stat, expanding_slopes

    panels = {s: build_annual(s) for s in SCENARIOS}
    ref = panels[SCENARIOS[0]]
    years, out_ix, members = ref[2], ref[3], ref[4]
    n_cells = out_ix.size

    # SHARED 2020s BASELINE. OUTPUT-SPEC requires the baseline panel to be bit-identical
    # across scenarios, so each member's 2020s window is averaged ACROSS scenarios before
    # pooling. Pooling each scenario's own 2020s would leave the panels merely similar and
    # the contract test checks identity.
    base_sel = (years >= BASELINE_DECADE) & (years < BASELINE_DECADE + 10)
    shared_exp = np.nanmean(np.stack([panels[s][0][:, base_sel] for s in SCENARIOS]), axis=0)
    shared_pct = np.nanmean(np.stack([panels[s][1][:, base_sel] for s in SCENARIOS]), axis=0)
    base_years = years[base_sel]

    # Internally everything is ASCENDING in latitude, because the sea-level grid, the GEBCO
    # grid and the flat cell indices all are. The product convention for a PROCESSED file is
    # DESCENDING, so the flip happens once, here, at write time -- not in the indexing, where
    # it would have to be right in four places instead of one.
    lat = np.arange(-90 + 0.25, 90, 0.5)
    lon = np.arange(-180 + 0.25, 180, 0.5)
    jj, ii = np.unravel_index(out_ix, (NLAT, NLON))
    OUT.mkdir(parents=True, exist_ok=True)

    for scen in SCENARIOS:
        exposed, pct, _, _, _ = panels[scen]
        grids = {k: np.full((len(DECADES), NLAT, NLON), np.nan, "float32")
                 for k in ("median", "lower_ci", "upper_ci", "percentile",
                           "ols_slope", "sen_slope")}
        # float32 with NaN off-mask, matching the shipped layers. Writing 0 off-mask makes
        # the whole ocean read "zero members", and the contract test checks the WHOLE
        # array, not just the masked part.
        n_mem = np.full((len(DECADES), NLAT, NLON), np.nan, "float32")
        n_mod = np.full((len(DECADES), NLAT, NLON), np.nan, "float32")

        for di, dec in enumerate(DECADES):
            if dec == BASELINE_DECADE:
                a_exp, a_pct, yrs = shared_exp, shared_pct, base_years
            else:
                a_exp, a_pct, yrs = exposed, pct, years
            stat, lo, hi = pooled_decadal_stat(a_exp, yrs, dec, boolean=False)
            pstat, _, _ = pooled_decadal_stat(a_pct, yrs, dec, boolean=False)
            sl = expanding_slopes(exposed, years, dec, BASELINE_DECADE)

            sel = (yrs >= dec) & (yrs < dec + 10)
            cnt = np.isfinite(a_exp[:, sel]).any(axis=1).sum(axis=0).astype("int16")
            # A layer whose mask varies by decade must mask its slopes to that decade's
            # median mask, or a cell absent from the decade gets a finite trend against a
            # NaN median and the contract check rejects it -- correctly.
            # decadal_stats returns PER YEAR; the shipped layers declare per DECADE, so
            # multiply by 10 exactly once and say so in slope_units. Fitting against a
            # decade index AND multiplying would inflate every trend 10x.
            ols = np.where(np.isfinite(stat), sl.ols_slope * 10.0, np.nan)
            sen = np.where(np.isfinite(stat), sl.sen_slope * 10.0, np.nan)

            for name, vals in (("median", stat), ("lower_ci", lo), ("upper_ci", hi),
                               ("percentile", pstat), ("ols_slope", ols), ("sen_slope", sen)):
                grids[name][di, jj, ii] = vals
            cntf = np.where(np.isfinite(stat) & (cnt > 0), cnt.astype("float32"), np.nan)
            n_mem[di, jj, ii] = cntf
            # ONE sea-level model (Mengel et al. 2016) forced by four GCMs, so the CI is
            # inter-GCM spread only and carries no structural model uncertainty.
            n_mod[di, jj, ii] = np.where(np.isfinite(cntf), 1.0, np.nan)
            print(f"  [{scen}] {dec}s done", end="\r")
        print()

        fin = np.isfinite(exposed)
        b_i = DECADES.index(BASELINE_DECADE)
        stats = {
            "zero_frac": float(((exposed == 0) & fin).sum() / fin.sum()),
            "median_exposed": int(np.nansum(
                np.nanmedian(shared_exp, axis=(0, 1)) > 0)),
            "mean_exposed": int(np.nansum(np.nanmean(shared_exp, axis=(0, 1)) > 0)),
            "pct_at_floor": float(
                np.nansum(grids["percentile"][b_i] <= 1.0)
                / max(np.isfinite(grids["percentile"][b_i]).sum(), 1)),
        }
        stats["erasure"] = 1 - stats["median_exposed"] / max(stats["mean_exposed"], 1)
        flip = slice(None, None, -1)
        ds = xr.Dataset(
            {k: (("decade", "lat", "lon"), v[:, flip, :]) for k, v in grids.items()}
            | {"n_members": (("decade", "lat", "lon"), n_mem[:, flip, :]),
               "n_models": (("decade", "lat", "lon"), n_mod[:, flip, :])},
            coords={"decade": DECADES, "lat": lat[flip], "lon": lon},
            attrs=_attrs(scen, members, stats),
        )
        f = OUT / f"{VARIABLE}_{scen}_processed.nc"
        ds.to_netcdf(f, encoding={k: {"zlib": True, "complevel": 4} for k in grids})
        print(f"[{scen}] wrote {f} ({f.stat().st_size / 1e6:.0f} MB)")
        ds.close()


def write_members(scen: str = "rcp60"):
    """Per-member diagnostic panel for the dashboard's Members tab.

    Not part of the OUTPUT-SPEC contract -- a diagnostic. Every statistic in the QA table
    is invariant under spatial rearrangement, so the table cannot see a spatial defect; the
    per-member panel is what caught a ~4x5 degree member elsewhere in this repo after it had
    passed 37 algebraic checks twice. With MIROC5 masked per cell here, this is also the
    only view that shows WHERE each member was allowed to contribute.
    """
    exposed, _, years, out_ix, members = build_annual(scen)
    f = years >= 2090
    lat = np.arange(-90 + 0.25, 90, 0.5)
    lon = np.arange(-180 + 0.25, 180, 0.5)
    jj, ii = np.unravel_index(out_ix, (NLAT, NLON))
    grid = np.full((len(members), NLAT, NLON), np.nan, "float32")
    with np.errstate(invalid="ignore"):
        for k in range(len(members)):
            grid[k, jj, ii] = np.nanmean(exposed[k][f], axis=0)
    ds = xr.Dataset(
        {"value": (("member", "lat", "lon"), grid[:, ::-1, :])},
        coords={"member": members, "lat": lat[::-1], "lon": lon},
        attrs={"member_field": f"2090s mean exposed area fraction, {scen}",
               "note": "NaN where the member was dropped by the consensus mask or had no "
                       "sea-level cell within max_borrow_km"},
    )
    p = OUT / f"{VARIABLE}_members.nc"
    ds.to_netcdf(p, encoding={"value": {"zlib": True, "complevel": 4}})
    print(f"wrote {p} ({p.stat().st_size / 1e6:.1f} MB)")
    for k, m in enumerate(members):
        v = grid[k][np.isfinite(grid[k])]
        print(f"  {m:14s} cells {v.size:6,}  mean {v.mean():.5f}  p95 {np.percentile(v,95):.5f}"
              f"  max {v.max():.5f}  zero% {100*(v==0).mean():.1f}")
    ds.close()


def _attrs(scen, members, stats):
    """Global attributes.

    THE DELIVERY ATTRIBUTE ALLOWLIST IS CLOSED (scripts/utils/delivery.py). A caveat
    written under any other name is silently dropped on the way to layers.csv and never
    reaches a report, so the two that matter are placed on allowlisted attributes rather
    than on descriptive ones of my own naming:

      sparsity_caveat      -> promoted to MUST-DISCLOSE. This layer is pinned at its floor
                              over most of the domain and the floor is partly an artefact
                              of the DEM, not of the coast.
      interpretation_caveat -> should-note: defences, the absence of surge, and the
                              scenario ceiling.
    """
    return {
        "title": "Coastal inundation exposure from sea-level rise and terrain elevation",
        "variable": VARIABLE,
        "scenario": scen,
        "units": "1",
        "long_name": "area fraction of the cell's coastal land at or below projected sea level",
        "spatial_resolution_degrees": 0.5,
        "baseline_decade": BASELINE_DECADE,
        "baseline_source": "shared_across_all_scenarios",
        "window_years": 10,
        "slope_units": "1 decade-1",
        "decadal_statistic": "pooled_median",
        "decadal_statistic_rationale":
            "Median + IQR, the contract default (user decision 2026-08-14). Measured: "
            f"{stats['zero_frac']:.1%} of the exposure field is exactly zero and the median "
            f"branch publishes {stats['median_exposed']:,} exposed cells in the 2020s against "
            f"the mean branch's {stats['mean_exposed']:,} -- a {stats['erasure']:.0%} erasure. "
            "Short of `let` (97.8% zero / 93% erased) and `cropfailure-3b` (96.1% / 96.6%), "
            "above `burntarea` (29.2%, median branch).",
        "field_nature": "continuous, zero-inflated area fraction in [0, 1]",
        "percentile_direction": "higher_is_worse",
        "percentile_zero_fraction": float(stats["pct_at_floor"]),
        "percentile_method":
            "DEVIATION FROM OUTPUT-SPEC: an ABSOLUTE calibration, not percentile-of-score "
            "against the 2020s baseline. Elevation above projected sea level <=0 m -> 100, "
            ">=10 m -> 1, between -> 2..99 by the empirical distribution of coastal "
            "elevations, fixed once from the 2020s and reused. Area-weighted over sub-cells. "
            "Because it is absolute, this layer's percentile is NOT on the same axis as "
            "other layers' and any Climate Score averaging them must say so.",
        "lecz_bound_m": LECZ_M,
        "lecz_source": "McGranahan, Balk & Anderson (2007), Environment and Urbanization",
        "relative_baseline": "false -- 0 m reads 100 regardless of other cells",
        "spatial_smoothing": "none -- the ensemble is thin (4 members) but the field is not "
                             "a sparse track field like `let`; coastal structure is real "
                             "terrain at 15 arcsec, and a kernel would bleed exposure across "
                             "the coastline into cells with no low-lying land.",
        "normalization": "none -- all four members are the same sea-level model under "
                         "different GCM forcing, in metres, on one grid. Model democracy.",
        "percentile_baseline": "NOT a 2020s score distribution. The calibration reference is "
                               "the pooled empirical distribution of coastal-land elevations "
                               "inside the LECZ, built once from the 2020s panel and reused "
                               "for every decade and scenario so a percentile means the same "
                               "thing over time. See percentile_method.",
        "mask_rule": "Land = GEBCO_2026 TID == 0 (the publisher's own verdict), NOT an "
                     "elevation threshold -- thresholding flood-fills below-sea-level land "
                     "behind sub-pixel dikes into the ocean and deletes the Netherlands. "
                     "Cells are published only where a sea-level member solved within "
                     "250 km; 419 coastal cells and the whole Persian Gulf interior fall "
                     "outside that and are NaN rather than filled.",
        "value_note": "Area fraction of the cell's coastal land at or below projected mean "
                      "sea level. Denominator is ALL land in the cell, including land above "
                      "the LECZ and land the sea cannot reach, so the value is a share of "
                      "real coastal land and not of a subset the method chose.",
        "ci_definition": "25th/75th percentile of the pooled (year x member) decade sample, "
                         "the contract's continuous branch. With ONE sea-level model the "
                         "spread is inter-GCM plus interannual only and carries NO "
                         "structural model uncertainty.",
        "slope_definition": "OLS and Theil-Sen over the expanding (year x member) stack from "
                            "the start of the 2020s baseline through the end of the target "
                            "decade, per OUTPUT-SPEC, then multiplied by 10 to express per "
                            "decade. READ ols_slope: sen_slope is exactly zero on 100% of "
                            "active cells.",
        "source_dataset": "ISIMIP2b InputData/sealevelrise `total` p50 (Mengel et al. PNAS "
                          "2016; Bamber & Riva 2010) x GEBCO_2026 15 arcsec terrain "
                          "(land base SRTM15+ v2.8, land mask TID==0)",
        "impact_models": "mengel2016 regional sea-level model (ONE model)",
        "gcms": ", ".join(members),
        "members_by_scenario": "; ".join(f"{sc}: {', '.join(members)}" for sc in SCENARIOS),
        "ensemble_uniform_across_scenarios": "true",
        "n_models_note": "ONE sea-level model forced by four CMIP5 GCMs -- the CI is "
                         "inter-GCM spread only and carries no structural model uncertainty.",
        "member_qc": "MIROC5 is dropped per cell where it deviates more than 0.5 m from the "
                     "median of all members in any decade of either scenario. Unmasked it "
                     "gives a FALLING Mediterranean (-0.753 m) and Black Sea (-1.053 m) by "
                     "the 2090s. See scripts/check_slr_member_outliers.py",
        "resolution_caveat":
            "SCREENING RESOLUTION -- THIS IS A REGIONAL RESULT, NOT A SITE RESULT. Every "
            "value is published on a 0.5 degree grid, about 55 km across at the equator, so "
            "the number returned for a site is the statistic for the whole cell containing "
            "it. Terrain is read at 15 arcsec (~460 m) and the full elevation distribution "
            "inside each cell is used -- it is not averaged away -- but sea level itself is "
            "one value per 0.5 degree cell, and the published result is one value per cell. "
            "Coastal inundation turns on metres of elevation across hundreds of metres of "
            "ground: a quayside and a hillside two kilometres apart sit in the same cell and "
            "receive the SAME number. Hydrological connectivity is decided at ~1.85 km, so "
            "any barrier, inlet or channel narrower than that is not resolved. Treat this "
            "How much it matters is measurable: a site on a cell boundary falls into one of "
            "two adjacent cells, and at 52.5N 4.75E -- the Dutch coast -- the two give 0.334 "
            "and 0.367 exposed fraction (percentile 73 and 79), while at 33.0S 71.5W the "
            "percentile moves from 3 to 8. Same site, same data, different 55 km box. Treat "
            "this layer as a FIRST-PASS SCREEN that ranks which coastlines and which sites "
            "deserve investigation. It cannot support a site-level conclusion, a design "
            "elevation, an asset-level financial estimate or any statement about an "
            "individual building or berth, and a site-specific study at metre-scale "
            "elevation is required before any of those.",
        "spatial_support_km": 55.0,
        "sparsity_caveat":
            f"This layer reads EXACTLY ZERO over most of its domain -- {stats['zero_frac']:.0%} "
            f"of the exposure field and {stats['pct_at_floor']:.0%} of cells sit at the "
            "percentile floor -- and the zero is only partly a statement about the coast. "
            "Elevation comes from GEBCO, a digital SURFACE model that includes vegetation "
            "canopy and buildings; SRTM-class coastal bias runs +2.49 m (Australia) to "
            "+3.67 m (US) against a sea-level signal of +0.21 to +0.34 m by the 2090s, so "
            "ABSOLUTE exposure is systematically UNDERSTATED and correcting the bias is "
            "measured to TRIPLE exposure estimates (Kulp & Strauss 2019). The worked example "
            "is South Florida, which reads 0.00000 at every decade although 99.9% of its "
            "land lies below 10 m. A zero here means 'no land below projected sea level in "
            "this DEM', NOT 'no coastal risk', and a flat slope at the floor is censored "
            "rather than measured. Read the percentile, which still ranks South Florida at "
            "48, before concluding anything from a zero.",
        "interpretation_caveat":
            "Permanent-inundation potential from MEAN sea level, not coastal flood "
            "frequency: no tide, no storm surge and no wave setup are included, and most "
            "coastal damage occurs at extreme water levels well above mean sea level. "
            "COASTAL DEFENCES ARE NOT REPRESENTED -- dikes, levees and seawalls are narrower "
            "than a 460 m grid cell and absent from the DEM, so the most heavily defended "
            "coasts read exactly like undefended ones; the Netherlands is the standing "
            "example at 33% of its coastal land below projected sea level by the 2090s. "
            "Scenario coverage is rcp26 and rcp60 only -- no high-forcing member exists for "
            "this product, so the high end of sea-level rise is not representable. Local "
            "land subsidence is not included and exceeds the climate signal in several "
            "deltas. Values are borrowed from the nearest cell where the ocean model solved, "
            "a median of 74-98 km away.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--members-only", action="store_true",
                    help="write only the per-member diagnostic panel")
    ap.add_argument("--measure", action="store_true",
                    help="report zero-inflation and both branches, write nothing")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    if not HYPSO.exists():
        raise SystemExit(f"{HYPSO} not built yet -- run build_coastal_hypsometry.py first")
    if a.members_only:
        write_members()
        return
    if a.measure:
        measure()
        return
    run()


if __name__ == "__main__":
    main()
