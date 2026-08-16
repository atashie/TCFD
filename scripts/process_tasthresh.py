"""Stage 2: turn the ISIMIP3b threshold-count ladder into TCFD/CDP layers.

Reads the stage-1 interim counts (data/interim/tasthresh/{GCM}_{scenario}_counts.nc,
written by download_reduce_tasthresh_isimip3b.py) and emits one OUTPUT-SPEC layer per
rung. Nine rungs, one ingest: the daily data was read once and tested against every
threshold, so the rungs cost nothing extra and no rung needs a second 1.34 TB.

    hd30 hd35 hd40 hd45   days with tasmax above 30/35/40/45 C      heat-stress-chronic
    TR20 TR25             days with tasmin above 20/25 C            heat-stress-chronic
    ID                    days with tasmax below 0 C  (ETCCDI ID)   cold-frost
    FD FDm10              days with tasmin below 0/-10 C (ETCCDI FD) cold-frost

WHAT THIS LAYER IS, AND WHY IT EXISTS ALONGSIDE `heatwave-3b`. Every threshold here is
ABSOLUTE: 35 C is 35 C in Singapore and in Stockholm. `heatwave-3b` scores departure from
each cell's OWN preindustrial distribution. The two answer different questions and neither
substitutes for the other -- this one says how often a fixed physical threshold is crossed,
that one says how unusual the year was for that place.

That is not a stylistic distinction, it is a different site list. Measured 2026-08-16 over
the 65,797 shared land cells, 2020s panel, ssp585:

    Spearman(hd35, heatwave-3b)  +0.554        top-decile agreement   47.2%
    Spearman(hd30, heatwave-3b)  +0.631        (3,105 of 6,580 cells)
    Spearman(FD,   heatwave-3b)  -0.662

    site          hd35 2020s   hd35 2090s   heatwave-3b percentile (2020s)
    Delhi              130.6        228.6                    53.0
    Chicago             10.0         71.8                    75.4
    Phoenix AZ         144.0        205.5                    60.8
    Singapore            4.8        268.8                    90.0
    Kuwait City        188.9        230.7                    97.0
    Stockholm            0.0          0.9                    22.1

Delhi is hot and used to it; Chicago is not hot and is not used to it. Screening on one
layer hands you a different half of the worst decile than screening on the other. Singapore
is the sharpest case in the other direction: 4.8 hot days now and 268.8 by the 2090s, a
56x rise, because its temperature distribution is narrow and sits just below 35 C -- a
threshold count reads near zero right up until the distribution crosses, so a low value
here is NOT evidence of a safe margin. Ship the two layers together, and say which question
each answers, or a reader will use one to answer the other's.

---------------------------------------------------------------------------
THE MEASUREMENTS THIS PROCESSOR RESTS ON (all 12 members, ssp585, 2026-08-16)
---------------------------------------------------------------------------

1. THE LAND MASK IS `landseamask_no-ant.nc`, NOT THE FULL MASK. This is measured, not a
   convention borrowed from a sibling layer. The full ISIMIP3b mask has 92,889 cells and
   the no-Antarctica mask 65,797; the 27,092-cell difference is Antarctica, 29% of the
   full mask. Antarctica is not a neutral addition to a temperature-threshold layer, it
   is a saturated block that breaks two things at once:

       rung / decade      ceiling share (>=364 d)   full mask   no-ant   Antarctica alone
       FD    2020s                                    30.25%     2.01%        98.85%
       ID    2020s                                    28.19%     1.08%        94.04%
       FDm10 2020s                                    21.25%     0.00%        72.83%

   Carrying Antarctica censors ~30% of the FD pool at the calendar ceiling, which is the
   regime OUTPUT-SPEC warns about: at a bound BOTH slope estimators go to ~0 and AGREE, so
   the dual-slope disagreement rule stops warning. It also fills the top of the frost
   percentile with a continent nobody sites an asset on, pushing genuinely severe
   mid-latitude frost down the ranking. Dropped -- and the drop is what makes the cold
   rungs readable, not a cosmetic choice.

   NOTE this is the ISIMIP3b W5E5 mask, 65,797 cells, against `heatwave-3b`'s 67,420. The
   1,623-cell gap is that layer's own member-coverage union, not a different landmask
   decision; percentiles remain comparable to within that.

2. THE WHOLE LADDER TAKES `pooled_mean_zero_inflated`, AND THE ALTERNATIVE WAS MEASURED.
   The median branch does not merely lose precision on the hot rungs, it erases them --
   cells whose pooled median is exactly 0 while their mean is strictly positive, i.e.
   places that DO cross the threshold, published as never crossing it:

       rung   decade   cells erased   % of land   their true mean (d/yr)
       hd45   2090s        29,058       44.16%           0.87
       hd40   2050s        25,597       38.90%           0.70
       hd35   2020s        18,428       28.01%           0.52
       hd30   2020s        10,418       15.83%           0.60
       TR25   2090s        17,772       27.01%           0.87
       FD     2020s         5,350        8.13%           0.30

   Going the other way costs almost nothing where the median works: on FD the 2020s pooled
   median is 106.0 days against a mean of 104.3, a 1.6% difference. So the mean is taken
   for EVERY rung -- an expected annual count, which is also the physically additive
   quantity -- rather than splitting the ladder across two branches and having hd35 report
   an expected count while FD reports a typical year. A customer reads these rungs side by
   side; two statistics under one family would be a trap. DECLARED, per OUTPUT-SPEC, in
   `decadal_statistic` + `decadal_statistic_rationale`. Not taken to improve contrast.

3. THE CENSORED RUNGS ARE THE HOT ONES, NOT THE COLD ONES. On the no-ant mask at ssp585
   2090s, hd30 sits at the 364-day ceiling for 11.80% of pooled observations and TR20 for
   11.36%, while FD is 0.63% and FDm10 exactly 0.00%. A `saturation_caveat` is therefore
   emitted for the rungs that earn it, per panel, from the measurement -- not asserted for
   the family. (An earlier reading of this had it backwards; it was measured on the full
   mask, where Antarctica supplied the censoring.)

4. THE CALENDAR RISK DID NOT MATERIALISE, AND THAT IS A MEASUREMENT TOO. A 360_day member
   gets 360 chances a year to cross a threshold and a Gregorian one 365.25 -- a ~1.5% bias
   that is CONSTANT PER MEMBER, so it never averages out. All 32 members measured
   `proleptic_gregorian`, all Kelvin, all 86 years: ISIMIP's bias adjustment puts every GCM
   on a real calendar first. No correction is applied because none is needed. Stage 1's
   `*_days_per_year_HETEROGENEOUS` attribute is an artifact of dividing each chunk's day
   count by its own unique-year count (a 6-year chunk gives 365.33, a 10-year one 365.30)
   and is NOT evidence of a mixed calendar -- read `*_calendar`, which is uniform.

5. THE 12 GCMs ARE NOT 12 INDEPENDENT MODELS, AND ARE NOT POOLED BY FAMILY EITHER. Two
   lineage pairs exist: CNRM-CM6-1/CNRM-ESM2-1 share ARPEGE-Climat, NEMO and ISBA-CTRIP,
   and KACE-1-0-G runs the HadGEM3 atmosphere it shares with UKESM1-0-LL. Family pooling
   was tested by correlating each member's residual from the ensemble mean, and it does
   not survive: on FD the CNRM pair ranks 1 of 66 (+0.764, z=+2.14), but on hd35 it ranks
   14 of 66 (+0.228) behind KACE x UKESM (+0.546). The duplication is real, rung-dependent,
   and does not partition cleanly, so `n_models` counts GCMs and the non-independence is
   DECLARED in the attributes instead of being silently resolved one way or the other.
   Read the CI as spanning 12 correlated members, not 12 independent draws.

REFERENCE SITES (GUARDRAILS 12) are printed for every rung before writing. A count layer
fails visibly -- Phoenix must have hd40 and Yakutsk must have FDm10 -- so the check is
worth more here than on a field whose plausible range is wide.

Run:
    .venv/bin/python3 scripts/process_tasthresh.py --plan
    .venv/bin/python3 scripts/process_tasthresh.py --run [--rung hd35] [--jobs 8]
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.decadal_stats import (  # noqa: E402
    MIN_SLOPE_OBS,
    expanding_slopes,
    pooled_decadal_stat,
)

INTERIM = Path("data/interim/tasthresh")
PROCESSED = Path("data/processed")
MASK_PATH = Path("data/masks/ISIMIP3b_landseamask_no-ant.nc")

DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2015, 2099
SLOPE_PER_DECADE = 10.0            # expanding_slopes returns per YEAR
TWO_TIER_ZERO_THRESHOLD = 0.02
CEILING_DAYS = 364.0               # "every day of the year", allowing for a short year
SATURATION_REPORT_THRESHOLD = 0.05  # declare a saturation_caveat above this share
SPARSITY_REPORT_THRESHOLD = 0.50    # declare a sparsity_caveat above this baseline-zero share
SLOPE_MEM_BUDGET_BYTES = 400 * 1024 * 1024

#: rung -> (folder stem, hazard, hazard_measure, source variable, comparison, threshold C)
RUNGS = {
    "hd30": ("heatdays-isimip3b", "Chronic heat",
             "Days per year with daily maximum temperature above 30 C",
             "tasmax", "gt", 30.0),
    "hd35": ("heatdays-isimip3b", "Chronic heat",
             "Days per year with daily maximum temperature above 35 C",
             "tasmax", "gt", 35.0),
    "hd40": ("heatdays-isimip3b", "Chronic heat",
             "Days per year with daily maximum temperature above 40 C",
             "tasmax", "gt", 40.0),
    "hd45": ("heatdays-isimip3b", "Chronic heat",
             "Days per year with daily maximum temperature above 45 C",
             "tasmax", "gt", 45.0),
    "TR20": ("tropicalnights-isimip3b", "Chronic heat",
             "Days per year with daily minimum temperature above 20 C (tropical nights)",
             "tasmin", "gt", 20.0),
    "TR25": ("tropicalnights-isimip3b", "Chronic heat",
             "Days per year with daily minimum temperature above 25 C",
             "tasmin", "gt", 25.0),
    "ID": ("icedays-isimip3b", "Cold and frost",
           "Days per year with daily maximum temperature below 0 C (ice days)",
           "tasmax", "lt", 0.0),
    "FD": ("frostdays-isimip3b", "Cold and frost",
           "Days per year with daily minimum temperature below 0 C (frost days)",
           "tasmin", "lt", 0.0),
    "FDm10": ("frostdays-isimip3b", "Cold and frost",
              "Days per year with daily minimum temperature below -10 C (hard frost)",
              "tasmin", "lt", -10.0),
}

#: (name, lat, lon, climate descriptor, {rung: what a working layer MUST show}).
#: A count layer fails visibly, so the expectation is written PER RUNG and printed only
#: for the rung being processed. A single expectation string shared across nine rungs is
#: not falsifiable -- it renders "hd40 well above zero" beside an FD value, where it is
#: neither true nor false, and a reader skims past it. Where a rung is not diagnostic at a
#: site, no claim is made rather than a vague one.
REFERENCE_SITES = [
    ("Phoenix AZ",  33.45, -112.07, "hot desert",
     {"hd35": ">100", "hd40": ">30", "FD": "<15", "TR25": ">30"}),
    ("Kuwait City", 29.38,   47.98, "hottest inhabited",
     {"hd35": ">150", "hd40": ">60", "FD": "0", "TR25": ">60"}),
    ("Singapore",    1.35,  103.82, "equatorial maritime",
     {"hd35": "<30 now, RISING STEEPLY later", "TR20": ">350", "FD": "0", "hd40": "0"}),
    ("Delhi",       28.61,   77.21, "monsoon subtropics",
     {"hd35": ">100", "hd40": ">20", "TR25": ">60", "FD": "0"}),
    ("Yakutsk",     62.03,  129.73, "extreme continental",
     {"FD": ">200", "FDm10": ">150", "ID": ">150", "hd35": "~0"}),
    ("Winnipeg",    49.90,  -97.14, "cold continental",
     {"FD": ">160", "FDm10": ">80", "ID": ">60", "hd35": "<10"}),
    ("Chicago",     41.88,  -87.63, "humid temperate",
     {"FD": ">80", "hd35": "small but NON-ZERO", "ID": ">20"}),
    ("Stockholm",   59.33,   18.07, "maritime boreal",
     {"FD": ">80", "hd30": "~0", "hd35": "0"}),
    ("Nairobi",     -1.29,   36.82, "tropical highland (1,795 m)",
     {"hd30": "LOW despite the tropics", "hd35": "0", "FD": "0"}),
    ("Manaus",      -3.12,  -60.02, "humid tropics",
     {"TR20": ">300", "hd35": ">30", "FD": "0"}),
]

_CUBE = None


def log(msg=""):
    print(msg, flush=True)


def load_mask():
    """Land cells from the no-Antarctica W5E5 mask.

    Read through `.filled(np.nan)`. A masked array's `.data` still holds 1e20 in the fill
    cells, so `np.asarray(mask) > 0.5` marks the WHOLE GLOBE as land -- which is exactly
    what happened on the first read of this file (259,200 "land" cells).
    """
    with Dataset(MASK_PATH) as ds:
        raw = ds.variables["mask"][:]
        arr = np.asarray(raw.filled(np.nan) if np.ma.isMaskedArray(raw) else raw,
                         dtype="f8").squeeze()
        mlat = np.asarray(ds.variables["lat"][:])
    return (np.isfinite(arr) & (arr > 0.5)), mlat


def discover_members():
    """{scenario: {gcm: path}} from what is on disk (GUARDRAILS 3 -- never hardcoded)."""
    out: dict[str, dict[str, Path]] = {}
    for p in sorted(INTERIM.glob("*_counts.nc")):
        if p.name.startswith("SMOKE_"):
            continue
        with Dataset(p) as ds:
            gcm, scen = ds.gcm, ds.scenario
            if getattr(ds, "partial_smoke_test", None):
                log(f"  SKIPPING {p.name}: marked partial_smoke_test")
                continue
        out.setdefault(scen, {})[gcm] = p
    return out


def make_pct_fn(baseline_flat):
    """Two-tier percentile-of-score against the 2020s baseline. Higher count = higher risk.

    Every rung of this ladder is materially zero-inflated somewhere (a threshold nobody
    crosses is a structural zero, not a missing value), so the tier is chosen from the
    measured baseline rather than assumed. Zeros take tier 1; positives rank against the
    NON-ZERO baseline into [2, 100]. No inversion: on both the hot and the cold rungs more
    days is worse -- the cold direction is the user's decision of 2026-08-14 ("more frost
    = worse"), and it means a frost-free tropical cell scores 1, correctly, because it has
    no frost risk at all.
    """
    frac_zero = float(np.mean(baseline_flat == 0.0))
    two_tier = frac_zero >= TWO_TIER_ZERO_THRESHOLD
    if two_tier:
        nz = np.sort(baseline_flat[baseline_flat > 0])
        n_nz = nz.size

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            res = np.ones(vals.shape, np.float32)
            pos = vals > 0
            if n_nz:
                frac = np.searchsorted(nz, vals[pos], side="left") / n_nz
                res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
            out[fin] = res
            return out.reshape(arr.shape)
    else:
        srt = np.sort(baseline_flat)
        n = srt.size

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            if n:
                frac = np.searchsorted(srt, flat[fin], side="right") / n
                out[fin] = np.clip(100.0 * frac, 1.0, 100.0).astype(np.float32)
            return out.reshape(arr.shape)

    return pct, ("two_tier" if two_tier else "single_tier"), frac_zero


def scatter(vec, land_idx, shape):
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = vec
    return out


def slope_chunk_cells(n_members, n_years, jobs=1):
    """Cells per block, sized so the memory budget is a TOTAL across workers."""
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    per_cell = 4 * pairs * 4
    budget = SLOPE_MEM_BUDGET_BYTES // max(jobs, 1)
    return max(4, min(512, int(budget // max(per_cell, 1))))


def _slope_block(task):
    s, e, years, decade, chunk = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, BASELINE_DECADE,
                           window_years=WINDOW_YEARS, chunk_cells=chunk)
    # Plain arrays, not the SlopeResult: its __getattr__ = dict.__getitem__ turns pickle's
    # __getstate__ probe into a KeyError and kills the pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def constant_cells(cube, years, decade):
    """Cells whose every finite observation in the expanding window is the same value.

    Both estimators are EXACTLY 0 there, by algebra rather than approximation: every
    pairwise difference y_j - y_i is 0, so the Theil-Sen median of the pairwise slopes is
    0, and the OLS numerator sum((x-xbar)(y-ybar)) is 0 because y is constant. So these
    cells can be filled directly instead of fitted.

    This is not an optimisation of convenience -- Theil-Sen is quadratic in stacked points
    (960 points = 455,040 pairs per cell) and this ladder is dominated by cells that never
    cross their threshold: hd45 is 92.2% exact-zero on land, hd40 75.4%. Fitting a
    half-million pairwise slopes to discover they are all zero is the bulk of the work on
    the sparse rungs.

    A cell with too few observations is NOT included: those must reach expanding_slopes so
    it can apply MIN_SLOPE_OBS and return NaN, which is a different answer from 0.
    """
    mask = (years >= BASELINE_DECADE) & (years <= decade + WINDOW_YEARS - 1)
    win = cube[:, mask, :]
    flat = win.reshape(-1, win.shape[2])
    fin = np.isfinite(flat)
    n_obs = fin.sum(axis=0)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        lo = np.nanmin(flat, axis=0)
        hi = np.nanmax(flat, axis=0)
    return (n_obs >= MIN_SLOPE_OBS) & np.isfinite(lo) & (hi == lo)


def compute_slopes(cube, years, decade, chunk, jobs, n_land):
    if decade == BASELINE_DECADE:
        return expanding_slopes(cube, years, decade, BASELINE_DECADE,
                                window_years=WINDOW_YEARS, chunk_cells=chunk)

    flat_cells = constant_cells(cube, years, decade)
    live = np.flatnonzero(~flat_cells)
    ols = np.full(n_land, np.nan, np.float32)
    sen = np.full(n_land, np.nan, np.float32)
    ols[flat_cells] = 0.0
    sen[flat_cells] = 0.0
    if live.size == 0:
        return {"ols_slope": ols, "sen_slope": sen}

    sub = np.ascontiguousarray(cube[:, :, live])
    if jobs <= 1:
        r = expanding_slopes(sub, years, decade, BASELINE_DECADE,
                             window_years=WINDOW_YEARS, chunk_cells=chunk)
        ols[live], sen[live] = r["ols_slope"], r["sen_slope"]
        return {"ols_slope": ols, "sen_slope": sen}

    global _CUBE
    keep = _CUBE
    _CUBE = sub                       # workers inherit the SUBSET through the fork
    try:
        edges = np.linspace(0, live.size, max(jobs * 8, 1) + 1).astype(int)
        tasks = [(int(a), int(b), years, decade, chunk)
                 for a, b in zip(edges[:-1], edges[1:]) if b > a]
        with mp.get_context("fork").Pool(jobs) as pool:
            for s, e, o, sn in pool.imap_unordered(_slope_block, tasks):
                ols[live[s:e]], sen[live[s:e]] = o, sn
    finally:
        _CUBE = keep
    return {"ols_slope": ols, "sen_slope": sen}


def site_indices(lats, lons, rung):
    """Grid indices plus THIS rung's expectation for each reference site."""
    out = []
    for name, la, lo, climate, expects in REFERENCE_SITES:
        lo360 = lo % 360 if lons.max() > 180 else lo
        exp = expects.get(rung)
        out.append((name, int(np.argmin(np.abs(lats - la))),
                    int(np.argmin(np.abs(lons - lo360))), climate, exp))
    return out


def process_rung(rung, members, mask, lats, lons, jobs, dry=False):
    folder_stem, hazard, measure, src_var, cmp_, thr = RUNGS[rung]
    out_dir = PROCESSED / f"{folder_stem}_{rung}_annual"
    scenarios = sorted(members)
    LAT, LON = mask.shape
    land_idx = np.flatnonzero(mask.ravel())
    n_land = land_idx.size

    log("=" * 74)
    log(f"{rung}: {measure}")
    log(f"  -> {out_dir}")
    log("=" * 74)

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    annual, members_by_scen = {}, {}
    for s in scenarios:
        gcms = sorted(members[s])
        members_by_scen[s] = gcms
        arr = np.full((len(gcms), n_years, n_land), np.nan, np.float32)
        for i, g in enumerate(gcms):
            with Dataset(members[s][g]) as ds:
                yr = np.asarray(ds.variables["year"][:], dtype=int)
                vals = np.asarray(ds.variables[rung][:], dtype="float32")
            for k, y in enumerate(yr):
                if int(y) in y_index:
                    arr[i, y_index[int(y)]] = vals.reshape(len(yr), -1)[k, land_idx]
            del vals
        annual[s] = arr
        log(f"  loaded {s}: {len(gcms)} members -> {arr.shape}")

    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("  WARNING: ensemble composition DIFFERS across scenarios. The shared 2020s "
            "baseline is only valid for a uniform ensemble; declaring members_by_scenario.")
        for s in scenarios:
            log(f"    {s}: {len(members_by_scen[s])} -- {','.join(members_by_scen[s])}")

    # ---- Shared 2020s baseline, pooled over every (year, member, scenario) ------
    base_pool = np.concatenate([annual[s] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, window_years=WINDOW_YEARS, central="mean")

    # What the median branch WOULD have published -- measured, not asserted (OUTPUT-SPEC
    # requires the deviation to be declared against the number it replaced).
    m_med, _, _ = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, window_years=WINDOW_YEARS, central="median")
    erased = int(np.count_nonzero((m_med == 0) & (b_med > 0)))
    erased_mean = float(np.mean(b_med[(m_med == 0) & (b_med > 0)])) if erased else 0.0

    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    ceil_base = float(np.mean(base_pool[:, bwin, :] >= CEILING_DAYS))
    del base_pool, m_med

    b_lo = np.clip(b_lo, 0.0, None)
    b_hi = np.clip(b_hi, 0.0, 366.0)
    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)

    log(f"\n  2020s shared baseline: land n={baseline_flat.size:,}  "
        f"mean={np.nanmean(b_med):.2f} d/yr  exact-zero={frac_zero:.2%}  "
        f"percentile={pct_mode}")
    log(f"  median branch would erase {erased:,} cells ({erased / n_land:.2%} of land) "
        f"whose true mean is {erased_mean:.2f} d/yr")
    log(f"  at the {CEILING_DAYS:.0f}-day ceiling: {ceil_base:.2%} of the baseline pool")

    # ---- Sub-grid variability, measured on THIS rung ---------------------------
    # The resolution caveat quotes a number rather than asserting coarseness. Adjacent
    # 0.5 deg cells are ~55 km apart; if they routinely differ by tens of days, terrain
    # inside a single cell plausibly does too, so the published cell mean matches neither
    # the valley floor nor the ridge.
    grid2d = scatter(b_med, land_idx, (LAT, LON))
    adj = np.concatenate([np.abs(np.diff(grid2d, axis=0)).ravel(),
                          np.abs(np.diff(grid2d, axis=1)).ravel()])
    adj = adj[np.isfinite(adj)]
    adj_p99 = float(np.percentile(adj, 99)) if adj.size else 0.0
    adj_max = float(adj.max()) if adj.size else 0.0
    log(f"  adjacent-cell |difference| on the baseline: p99 {adj_p99:.1f}, "
        f"max {adj_max:.1f} days per year")

    # ---- Reference sites, BEFORE writing anything (GUARDRAILS 12) --------------
    grid_lookup = np.full(LAT * LON, -1, np.int64)
    grid_lookup[land_idx] = np.arange(n_land)
    log(f"\n  Reference sites, 2020s baseline, {rung} (days per year). An expectation is "
        f"shown only where this rung is diagnostic:")
    for name, i, j, climate, expect in site_indices(lats, lons, rung):
        k = grid_lookup[i * LON + j]
        val = "OCEAN/masked" if k < 0 else f"{b_med[k]:7.1f}"
        flag = f"expect {expect}" if expect else "-"
        log(f"    {name:<14} {val:>12}   {climate:<28} {flag}")

    if dry:
        log("\n  --plan: stopping before the slope stage.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    sat_by_panel, sen_zero_by_panel = {}, {}
    # Every scenario is computed BEFORE anything is written, because `recommended_slope`
    # is a property of the RUNG and must not differ between that rung's own files. It is
    # decided from the sen==0 share, which measured 51-53% on hd35 -- straddling the 0.5
    # cut, so a per-scenario decision would have had ssp126 recommend ols_slope and ssp585
    # sen_slope for one layer, while the registry carries exactly one. Holding three
    # scenarios of output costs ~200 MB.
    outputs = {}

    for s in scenarios:
        global _CUBE
        cube = annual[s]
        _CUBE = cube
        mem = members_by_scen[s]
        chunk = slope_chunk_cells(len(mem), n_years, jobs)
        out = {k: np.full((len(DECADES), LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}
        log(f"\n  --- {s} ---")

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, window_years=WINDOW_YEARS, central="mean")
                lo = np.clip(lo, 0.0, None)
                hi = np.clip(hi, 0.0, 366.0)
                pc = pct(med)

            sl = compute_slopes(cube, years, d, chunk, jobs, n_land)

            # Slopes must not survive where the median does not: the expanding window
            # reaches back before the decade, so a cell could otherwise carry a trend
            # across a decade it has no observation in.
            gone = ~np.isfinite(med)
            leak = int(np.count_nonzero(
                gone & (np.isfinite(sl["ols_slope"]) | np.isfinite(sl["sen_slope"]))))
            if leak:
                sl = {"ols_slope": np.where(gone, np.nan, sl["ols_slope"]),
                      "sen_slope": np.where(gone, np.nan, sl["sen_slope"])}
                log(f"    masked {leak} slope cell(s) outside the {d}s median window")

            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
            n_mem[n_mem == 0] = np.nan
            # One GCM is one model here: this layer has no impact model, and the two
            # lineage pairs do not partition cleanly (see the module docstring).
            n_mod = n_mem.copy()

            for key, vec in (("median", med), ("lower_ci", lo), ("upper_ci", hi),
                             ("percentile", pc),
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            sat = float(np.mean(cube[:, win, :] >= CEILING_DAYS))
            sat_by_panel[(s, d)] = sat
            if d != BASELINE_DECADE:
                # ACTIVE cells only -- either slope non-zero -- which is the repo's
                # convention for this diagnostic (DATASET-ATTRIBUTES "Which slope to read")
                # and the only denominator that answers the question being asked. Measuring
                # over all FINITE cells instead counts every cell that never crosses the
                # threshold, where both slopes are correctly and exactly 0. Those are not
                # Sen failures, and including them inflated the share enough to flip the
                # recommendation on four rungs: ID and FDm10 read 0.542 over finite cells
                # against 0.182 over active ones.
                o, sn = sl["ols_slope"], sl["sen_slope"]
                act = np.isfinite(o) & np.isfinite(sn) & ((o != 0.0) | (sn != 0.0))
                sen_zero_by_panel[(s, d)] = (
                    float(np.mean(sn[act] == 0.0)) if act.any() else np.nan)
            log(f"    {d}s: mean={np.nanmean(out['median'][i]):6.2f} d/yr  "
                f"ceiling={sat:6.2%}  "
                + ("slopes=NaN (baseline)" if d == BASELINE_DECADE else
                   f"ols={np.nanmean(out['ols_slope'][i]):+.3f}  "
                   f"sen={np.nanmean(out['sen_slope'][i]):+.3f} d/dec  "
                   f"sen==0 on {sen_zero_by_panel[(s, d)]:.1%}"))

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), "baseline ols must be NaN"
                assert np.all(np.isnan(out["sen_slope"][i])), "baseline sen must be NaN"
                continue
            fin = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~fin
                assert not extra.any(), (
                    f"{k} finite where median is NaN at {d}s ({int(extra.sum())} cells)")

        outputs[s] = out
        del cube
        _CUBE = None

    # ---- Rung-level slope recommendation, then write ---------------------------
    # THE SLOPE IS DECIDED BY WHICH ESTIMATOR CAN FAIL HERE, NOT BY A THRESHOLD.
    # OUTPUT-SPEC's two estimators fail in opposite regimes, so the choice is settled by
    # asking which regime this layer is actually in -- and for this family the answer is
    # measured and the same for all nine rungs:
    #   * ols_slope's failure mode REQUIRES UNEVEN MEMBER COVERAGE (it absorbs between-member
    #     level offsets as trend). Measured: n_members is exactly 12 in EVERY cell of every
    #     rung and scenario. The bias has no mechanism here. Contrast `csoil`, where CLASSIC
    #     contributes 2 GCMs against 5 for the others and ols is genuinely biased.
    #   * sen_slope's failure mode -- collapse to exactly 0 where most year-pairs are tied --
    #     IS present on every rung, from 18.2% of active cells (ID, FDm10) to 85.4% (hd45),
    #     because a threshold count is zero wherever the threshold is not crossed.
    # So ols_slope is the read across the family: one estimator, for the same reason the
    # ladder takes one statistic. Both slopes are in every file and their DISAGREEMENT is
    # still the per-cell robustness signal; this only names the headline.
    finite_sen0 = [v for v in sen_zero_by_panel.values() if np.isfinite(v)]
    worst_sen0 = max(finite_sen0) if finite_sen0 else float("nan")
    best_sen0 = min(finite_sen0) if finite_sen0 else float("nan")
    slope_read = "ols_slope"
    log(f"\n  {rung}: sen_slope == 0 on up to {worst_sen0:.1%} of active cells across all "
        f"scenarios -> recommended_slope = {slope_read}")
    for s in scenarios:
        write_scenario(outputs[s], rung, s, out_dir, lats, lons, members_by_scen[s],
                       members_by_scen, scenarios, uniform, pct_mode, frac_zero, erased,
                       erased_mean, n_land, sat_by_panel, ceil_base, slope_read,
                       worst_sen0, adj_p99, adj_max)
    del outputs

    write_members(annual, members_by_scen, scenarios, years, bwin, rung,
                  out_dir, lats, lons, land_idx, (LAT, LON))


def saturation_text(rung, sat_by_panel, scenarios):
    worst = max(sat_by_panel.items(), key=lambda kv: kv[1])
    if worst[1] < SATURATION_REPORT_THRESHOLD:
        return (f"none -- the highest share of pooled observations at the "
                f"{CEILING_DAYS:.0f}-day calendar ceiling is {worst[1]:.2%} "
                f"({worst[0][0]} {worst[0][1]}s). Both slope estimators are readable.")
    return (
        f"MATERIAL CENSORING. {worst[1]:.1%} of pooled (year x member) observations sit at "
        f"the {CEILING_DAYS:.0f}-day calendar ceiling at {worst[0][0]} {worst[0][1]}s -- the "
        f"cell is above the threshold on essentially every day of the year and CANNOT go "
        f"higher. In that regime BOTH slope estimators are censored toward zero and AGREE, "
        f"so the dual-slope disagreement rule gives no warning: agreement there is ambiguous "
        f"between 'no trend' and 'maximally exposed, permanently'. Read `median` alongside "
        f"any slope, and treat a flat slope on a saturated cell as saturation, not stability. "
        f"Declared, never re-estimated.")


def sparsity_text(rung, frac_zero, thr, cmp_):
    """A rung nobody crosses is a sparse layer, and sparsity is a delivery-visible fact.

    `heatwave-2b` carries this attribute at 65.7% of land at zero; the top rungs here are
    far sparser than that, and a customer whose whole portfolio reads 0 needs to be told
    that is the layer's normal state rather than an extraction failure.
    """
    if frac_zero < SPARSITY_REPORT_THRESHOLD:
        return (f"none -- {frac_zero:.1%} of land is exactly zero on the 2020s baseline, "
                f"below the {SPARSITY_REPORT_THRESHOLD:.0%} at which sparsity is reported.")
    return (
        f"SPARSE LAYER. {frac_zero:.1%} of land never reaches {thr:g} C on the 2020s "
        f"baseline, so a majority of sites read exactly 0 and rank at percentile tier 1. "
        f"That is the measure working as defined, NOT a missing value and NOT an extraction "
        f"failure -- most of the world does not cross this threshold today. Two "
        f"consequences: a portfolio can be entirely zero and still gain exposure later "
        f"(this rung grows fastest where it starts at zero), and the percentile is "
        f"uninformative between zero sites, which are tied rather than ranked. For a "
        f"screening measure that discriminates across the whole land surface, use a lower "
        f"rung of the same ladder.")


def resolution_text(rung, adj_p99, adj_max):
    """0.5 deg support against a hazard that turns on elevation.

    Measured rather than asserted: the caveat quotes the observed spread between ADJACENT
    cells on this rung's own baseline. If neighbouring 55 km cells differ by tens of days,
    then sub-cell terrain inside one cell plausibly does too, and a site value cannot be
    read as that site's own count.
    """
    return (
        f"SCREENING ONLY -- 0.5 deg (~55 km) SUPPORT AGAINST A HAZARD THAT TURNS ON "
        f"ELEVATION. A site's value is its cell's value, and this rung's own baseline "
        f"shows ADJACENT cells differing by up to {adj_p99:.0f} days per year at the 99th "
        f"percentile (max {adj_max:.0f}). Temperature falls roughly 6.5 C per km of "
        f"elevation, so a cell spanning a valley floor and a ridge spans a threshold count "
        f"that varies over its own footprint by a comparable amount -- and the single "
        f"published number is the cell mean, matching neither. Urban heat island is not "
        f"represented either, so a city centre runs hotter than its cell. This layer RANKS "
        f"which sites deserve investigation; it cannot support a design temperature, an "
        f"asset-level exceedance count, or a site-specific engineering threshold.")


def write_scenario(out, rung, s, out_dir, lats, lons, mem, members_by_scen, scenarios,
                   uniform, pct_mode, frac_zero, erased, erased_mean, n_land,
                   sat_by_panel, ceil_base, slope_read, worst_sen0, adj_p99, adj_max):
    folder_stem, hazard, measure, src_var, cmp_, thr = RUNGS[rung]
    sat_txt = saturation_text(rung, sat_by_panel, scenarios)
    spa_txt = sparsity_text(rung, frac_zero, thr, cmp_)
    sat_applies = sat_txt.startswith("MATERIAL")
    spa_applies = spa_txt.startswith("SPARSE")

    ds_out = xr.Dataset(
        {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
        coords={"decade": DECADES, "lat": lats, "lon": lons},
        attrs={
            "variable": rung,
            "scenario": s,
            "long_name": measure,
            "units": "days per year",
            "spatial_resolution_degrees": round(float(abs(lats[1] - lats[0])), 6),
            "output_spec": "OUTPUT-SPEC.md",
            "decadal_statistic": "pooled_mean_zero_inflated",
            "field_nature": ("count of days per year, integer-valued on [0, 366]; "
                             "CONTINUOUS for statistical purposes, but zero-inflated -- "
                             "a threshold nobody crosses is a STRUCTURAL zero, not a "
                             "missing value"),
            "value_note": (
                f"median = MEAN over the pooled (year x member) sample inside the decade "
                f"window across {len(mem)} ISIMIP3b bias-adjusted GCMs. The variable is "
                f"named `median` by the OUTPUT-SPEC contract; the estimator here is the "
                f"MEAN, per decadal_statistic. Source: daily {src_var} tested "
                f"{'>' if cmp_ == 'gt' else '<'} {thr:g} C, counted per calendar year."),
            "ci_definition": (
                "lower_ci/upper_ci = mean -/+ 1 SD of the same pooled (year x member) "
                "sample, clipped to [0, 366]. The band therefore carries BOTH interannual "
                "variability and inter-model spread and is not a pure model-spread band. "
                "It is also NOT symmetric in meaning near the bounds: a cell whose mean is "
                "near 0 or near 366 has its band clipped by the calendar, not by agreement."),
            "slope_definition": (
                "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both are "
                "fitted over an EXPANDING window from the start of the 2020s baseline "
                "through the end of the target decade, stacking every (year, member) "
                "observation as an independent point. The baseline panel is NaN. The two "
                "fail in OPPOSITE regimes -- sen collapses to exactly 0 on zero-inflated "
                "fields, ols absorbs member level offsets as trend -- so their disagreement "
                "flags a cell whose trend is not robust. See saturation_caveat for where "
                "that rule stops working."),
            "slope_units": "days per year, per decade",
            "recommended_slope": slope_read,
            "recommended_slope_rationale": (
                f"ols_slope, decided by which estimator CAN fail on this layer rather than "
                f"by a threshold. (1) ols_slope's documented failure mode requires UNEVEN "
                f"member coverage, and coverage here is exactly {len(mem)} members in EVERY "
                f"cell of every rung and scenario -- so the bias has no mechanism. (2) "
                f"sen_slope collapses to exactly 0 where most year-pairs are tied, which on "
                f"a threshold count happens wherever the threshold is rarely crossed: "
                f"measured {best_sen0:.1%}-{worst_sen0:.1%} of ACTIVE cells (either slope "
                f"non-zero) across all {len(scenarios)} scenarios and every post-baseline "
                f"panel. One estimator across the ladder, for the same reason it takes one "
                f"statistic. BOTH slopes are in this file and their DISAGREEMENT remains the "
                f"per-cell signal that a trend is not robust -- this names the headline, not "
                f"the only number worth reading. NOTE the denominator: shares quoted over "
                f"ALL FINITE cells run far higher (ID reads 54.2% there against 18.2% here) "
                f"because cells that never cross the threshold have both slopes correctly "
                f"at 0 and are not Sen failures."),
            "percentile_baseline": (
                f"{pct_mode}: each cell's count is ranked against the 2020s ensemble land "
                f"distribution. Zeros take tier 1; positive counts rank against the "
                f"NON-ZERO baseline into [2, 100]."),
            "percentile_zero_fraction": round(frac_zero, 5),
            "percentile_direction": "higher_is_worse",
            "baseline_decade": BASELINE_DECADE,
            "baseline_source": "shared_across_all_scenarios",
            "members_by_scenario": ";".join(
                f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
            "ensemble_uniform_across_scenarios": str(uniform),
            "window_years": WINDOW_YEARS,
            "n_members": len(mem),
            "gcms": ",".join(mem),
            "impact_models": ("none -- this layer is derived directly from bias-adjusted "
                              "climate forcing, with no impact model in the chain"),
            "normalization": (
                "none. Every member reports the same physical quantity in the same unit "
                "(days per year) on the same grid and calendar, so members are "
                "equal-weighted raw (model democracy). Rescaling would destroy a unit that "
                "a reader can check against a weather station."),
            "spatial_smoothing": (
                f"none. {len(mem)} members x {WINDOW_YEARS} years = "
                f"{len(mem) * WINDOW_YEARS} draws per cell-decade, and the input is "
                f"bias-adjusted W5E5 forcing that is already spatially coherent; the thin-"
                f"ensemble smoothing case does not arise here."),
            "land_mask": (
                f"ISIMIP3b W5E5 landseamask_no-ant.nc, {n_land:,} cells. ANTARCTICA IS "
                f"EXCLUDED, on measurement: it is 27,092 cells (29% of the full 92,889-cell "
                f"mask) and carries 98.85% of FD's calendar-ceiling censoring, which would "
                f"censor ~30% of the cold-rung pool and fill the top of the frost percentile "
                f"with a continent holding no assets. NOTE the counts field is FINITE OVER "
                f"THE WHOLE GLOBE INCLUDING OCEAN -- `isfinite` is not a mask on this layer, "
                f"the landseamask is."),
            "decadal_statistic_rationale": (
                f"pooled_mean_zero_inflated, a DECLARED deviation from pooled_median, taken "
                f"on measurement. The median branch would publish exactly 0 for {erased:,} "
                f"land cells ({erased / n_land:.2%}) whose true expected count is "
                f"{erased_mean:.2f} days per year -- places that DO cross the threshold, "
                f"reported as never crossing it. Going the other way is nearly free where "
                f"the median works (on FD the 2020s pooled median is 106.0 days against a "
                f"mean of 104.3, 1.6%). The mean is applied to EVERY rung of the ladder so "
                f"that a customer reading hd35 beside FD is reading one statistic, not two. "
                f"Baseline pool at the calendar ceiling: {ceil_base:.2%}."),
            # saturation_caveat / sparsity_caveat are set ONLY WHEN THEY APPLY, and the
            # negative measurement goes to a *_measured attribute instead. Both names are on
            # delivery.LAYER_ATTRS_EXPORTED, and generate_delivery_caveats promotes them to
            # MUST-DISCLOSE on the test `if note and note.lower() != "nan"` -- i.e. on the
            # attribute being NON-EMPTY, not on it saying anything. Writing "none -- 0.69% at
            # the ceiling" would therefore publish a must-disclose caveat whose body says
            # "none", in both reports, for seven of the nine rungs. That is the same defect
            # load_registry() already guards against for a blank relative_baseline_note: a
            # heading promising a caveat with nothing under it. `heatwave-3b` carries the
            # attribute because it saturates; `cropfailure-3b` omits all three.
            # resolution_caveat is unconditional and that is correct -- every rung is 0.5 deg
            # support against a hazard that turns on elevation.
            "resolution_caveat": resolution_text(rung, adj_p99, adj_max),
            **({"saturation_caveat": sat_txt} if sat_applies
               else {"saturation_measured": sat_txt}),
            **({"sparsity_caveat": spa_txt} if spa_applies
               else {"sparsity_measured": spa_txt}),
            "source_dataset": (
                f"ISIMIP3b bias-adjusted (W5E5v2) daily {src_var}, "
                f"InputData/climate/atmosphere/bias-adjusted/global/daily (group-I GCMs) "
                f"and SecondaryInputData/... (extended GCMs), constructed into an annual "
                f"threshold-exceedance count by scripts/download_reduce_tasthresh_isimip3b.py. "
                f"Each daily chunk was verified against the publisher sidecar sha512 and "
                f"deleted after reduction; see data/interim/tasthresh/download_provenance.csv."),
            "interpretation_caveat": (
                "AN ABSOLUTE DRY-BULB THRESHOLD COUNT -- NOT A HEAT-STRESS INDEX AND NOT A "
                "RELATIVE ANOMALY. (1) It counts days crossing a FIXED air-temperature "
                "threshold, with NO humidity term. Human and livestock heat strain depends "
                "on humidity: 35 C at 80% relative humidity is far more dangerous than 40 C "
                "in dry desert air, and this layer cannot tell them apart. For a "
                "humidity-aware measure use `heatwave-3b` alongside it. (2) It is the "
                "ABSOLUTE counterpart to `heatwave-3b`, which scores departure from each "
                "cell's OWN preindustrial distribution. THE TWO ARE NOT INTERCHANGEABLE AND "
                "THE DIFFERENCE IS MEASURED, not argued: over the 65,797 shared land cells "
                "the 2020s rank correlation is only +0.554 (hd35), and the two layers agree "
                "on just 47.2% of their worst-decile cells -- screen on one and you get a "
                "materially different site list than screening on the other. Site examples, "
                "2020s: Delhi has 130.6 hot days a year and sits at the 53rd percentile of "
                "`heatwave-3b`, while Chicago has 10.0 and sits at the 75th. Neither layer "
                "is wrong; they answer different questions, and reading either as a "
                "substitute for the other inverts the answer. "
                "(3) The threshold is a screening value, not a damage function: an asset's "
                "actual failure temperature depends on its design, and 35 C is defensible "
                "as a shared reference point rather than as any specific asset's limit. "
                "(4) Daily tasmax/tasmin are cell-mean values at ~55 km; urban heat-island "
                "and local terrain effects are NOT represented, and a city centre will run "
                "hotter than the cell that contains it. (5) THE 12 GCMs ARE NOT 12 "
                "INDEPENDENT MODELS, so the confidence band is narrower than it looks: "
                "CNRM-CM6-1 and CNRM-ESM2-1 share ARPEGE-Climat, NEMO and ISBA-CTRIP, and "
                "KACE-1-0-G runs the HadGEM3 atmosphere it shares with UKESM1-0-LL. Family "
                "pooling was tested by correlating each member's residual from the ensemble "
                "mean and REJECTED because the duplication is rung-dependent and does not "
                "partition cleanly: on FD the CNRM pair ranks 1 of 66 pairs (+0.764), on "
                "hd35 it ranks 14 (+0.228), behind KACE x UKESM (+0.546). n_models counts "
                "GCMs; read the CI as spanning correlated members, not independent draws."),
            "created_by": "scripts/process_tasthresh.py",
        },
    )
    enc = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
               "_FillValue": np.float32(np.nan)} for k in out}
    path = out_dir / f"{rung}_{s}_processed.nc"
    ds_out.to_netcdf(path, encoding=enc)
    log(f"    wrote {path.name} ({path.stat().st_size / 1e6:.1f} MB)")


def write_members(annual, members_by_scen, scenarios, years, bwin, rung, out_dir,
                  lats, lons, land_idx, shape):
    """Per-member 2020s diagnostic -- the Members tab, and the contact sheet to LOOK at."""
    names = sorted({m for s in scenarios for m in members_by_scen[s]})
    grid = np.full((len(names), *shape), np.nan, np.float32)
    for mi, name in enumerate(names):
        stack = [annual[s][members_by_scen[s].index(name)][bwin]
                 for s in scenarios if name in members_by_scen[s]]
        if not stack:
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            flat = np.nanmean(np.concatenate(stack, axis=0), axis=0)
        grid[mi] = scatter(flat, land_idx, shape)
    ds = xr.Dataset(
        {"value": (["member", "lat", "lon"], grid)},
        coords={"member": names, "lat": lats, "lon": lons},
        attrs={"variable": rung, "units": "days per year",
               "member_field": (f"mean annual count over the {BASELINE_DECADE}s baseline "
                                "decade, pooled across all scenarios"),
               "note": "Diagnostic only -- not part of the OUTPUT-SPEC contract."},
    )
    p = out_dir / f"{rung}_members.nc"
    ds.to_netcdf(p, encoding={"value": {"dtype": "float32", "zlib": True, "complevel": 4,
                                        "_FillValue": np.float32(np.nan)}})
    log(f"  wrote per-member diagnostic {p.name} ({len(names)} members)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--rung", action="append", default=None,
                    help="process only these rungs (default: all nine)")
    ap.add_argument("--jobs", type=int, default=max(1, (mp.cpu_count() or 2) - 2))
    a = ap.parse_args()

    members = discover_members()
    if not members:
        log(f"ERROR: no interim counts found in {INTERIM}")
        return 1
    mask, mlat = load_mask()
    first = next(iter(next(iter(members.values())).values()))
    with Dataset(first) as ds:
        lats = np.asarray(ds.variables["lat"][:])
        lons = np.asarray(ds.variables["lon"][:])
    if np.allclose(lats, mlat[::-1]):
        mask = mask[::-1]
        log("NOTE: mask latitude axis flipped to match the counts grid")
    elif not np.allclose(lats, mlat):
        log("ERROR: mask and counts are on different grids")
        return 1

    rungs = a.rung or list(RUNGS)
    unknown = [r for r in rungs if r not in RUNGS]
    if unknown:
        log(f"ERROR: unknown rung(s) {unknown}; known: {list(RUNGS)}")
        return 1

    log(f"interim       {INTERIM}")
    log(f"land mask     {MASK_PATH.name} -- {int(mask.sum()):,} cells (Antarctica EXCLUDED)")
    for s in sorted(members):
        log(f"scenario      {s}: {len(members[s])} members "
            f"({','.join(sorted(members[s]))})")
    log(f"rungs         {rungs}")
    log(f"jobs          {a.jobs}")
    complete = {s for s in members if len(members[s]) == 12}
    if len(complete) != len(members):
        log("\nWARNING: not every scenario has all 12 members. The shared 2020s baseline "
            "assumes a uniform ensemble; a subset run will declare members_by_scenario, "
            "but the layer should be rebuilt once stage 1 finishes.")

    if not (a.run or a.plan):
        log("\nNothing to do. Pass --plan or --run.")
        return 0

    for rung in rungs:
        process_rung(rung, members, mask, lats, lons, a.jobs, dry=not a.run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
