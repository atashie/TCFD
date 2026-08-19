"""Stage 2: turn the ISIMIP3b precipitation statistics into TCFD/CDP layers.

Reads the stage-1 interim files (data/interim/prthresh/{GCM}_{scenario}_pr.nc, written by
download_reduce_prthresh_isimip3b.py) and emits one OUTPUT-SPEC layer per metric. Ten
metrics, one ingest: ~764 GB of daily `pr` was read once, so a metric we do not publish
today needs no second pass.

TEN METRICS IN THREE FAMILIES, and the families answer different questions.
  ABSOLUTE frequency, days/yr -- "how often does a damaging depth fall here"
    R10mm R20mm   ETCCDI heavy / very-heavy precipitation days
    R50mm         World Bank CCKP ships r50mm
    R100mm        OUR EXTENSION, not a standard index
  RELATIVE frequency, days/yr -- "how often is a rainy day extreme FOR HERE"
    R95pD R99pD   days above this cell's OWN 2020s wet-day p95/p99. `relative_baseline`.
                  NOT ETCCDI `R95p`, which is a TOTAL in mm on a 1961-1990 base.
  INTENSITY / TOTAL, mm
    Rx1day        annual maximum 1-day precipitation -- the best single pluvial proxy, and
                  the only metric here needing no threshold choice at all
    Rx5day        annual maximum 5-day total: the saturation-driven and compound case
    prcptot       annual total precipitation from wet days
  wetdays ships too: it is the denominator the relative rungs are defined on, so a reader
  who wants to check that R95pD is ~5% of wet days can.

---------------------------------------------------------------------------
THE MEASUREMENTS THIS PROCESSOR RESTS ON (all 42 members, 2026-08-18)
---------------------------------------------------------------------------

1. TWO STATISTIC BRANCHES, SPLIT BY UNIT FAMILY, EACH ON MEASUREMENT. The pooled median
   erases cells that genuinely experience the event -- publishing exactly 0 where the
   expected count is positive:

       metric   zeros(2020s)  cells erased by median   their true mean
       R50mm         73.0%      35,228  (53.5% land)      0.18 d/yr
       R99pD         50.2%      28,086  (42.7%)           0.47
       R100mm        94.3%      26,417  (40.2%)           0.14
       R20mm         26.0%      15,511  (23.6%)           0.34
       R10mm          6.4%       3,063  ( 4.7%)           0.36
       Rx1day         0.02%          0  ( 0.00%)          --

   So the seven DAY COUNTS take `pooled_mean_zero_inflated`: the median would erase half the
   land on R50mm, and going the other way costs ~1% where the median works (R10mm median
   13.00 against mean 13.15; wetdays 86.00 against 86.32). The three mm-valued metrics
   (Rx1day, Rx5day, prcptot) are continuous and essentially zero-free -- the median erases
   NOTHING there -- so they keep OUTPUT-SPEC's default `pooled_median`.

   This is a split WITHIN one ingest, which the heat ladder deliberately avoided. It is
   defensible here and was not there because these are different QUANTITIES in different
   UNITS: a reader of Rx1day (mm) already knows it is a different object from R20mm (d/yr),
   whereas hd35 and FD were both day counts and a split would have been a trap.

2. PRECIPITATION IS NOT A SECOND TEMPERATURE LADDER, which is why both families ship. Each
   cell's 2020s wet-day p95 spans 3.3-165.2 mm/day, median 19.0 (all 42 members). A fixed
   20 mm threshold therefore sits at roughly the MEDIAN cell's own p95 -- which is why it is
   the absolute headline -- while being unreachable in the driest cells and an ordinary wet
   day in the wettest. Daily temperature varies about two-fold globally and a fixed threshold
   means the same thing everywhere; precipitation does not.

3. THE PERCENTILE IS OVER WET DAYS, NOT CALENDAR DAYS, and this is load-bearing. Verified
   over 8 members: R95pD/wetdays = 0.0509, R99pD/wetdays = 0.0093, as constructed. On an
   ALL-DAYS basis the metric breaks twice: below 18.25 wet days a year the p95 of all days is
   exactly 0 mm and "exceedance" collapses into "did it rain at all" (Cairo has 4.7 wet days
   a year), and above that it slides what "extreme" means with local rain frequency -- p95 of
   all days is ~the 45th percentile of WET days in Phoenix against the 92nd in Singapore.

4. THE RELATIVE RUNGS ARE A `relative_baseline` PAIR AND THEIR 2020s PANEL IS A CONTROL BY
   CONSTRUCTION. If the threshold IS the cell's own 2020s p95, then in the 2020s every cell
   sits at ~5% of its wet days by definition, so the baseline panel is essentially
   0.05 x wet-day count -- a map of how often it rains, not of hazard. The signal is the
   CHANGE. Accepted deliberately (user decision 2026-08-17: "that is by design and not
   problematic"), and declared: `relative_baseline: true` promotes it to a MUST-DISCLOSE
   caveat in every report. Worked example from the data: Singapore reads R95pD 14.9 against
   Cherrapunji's 6.8, while on the ABSOLUTE rung Cherrapunji reads R20mm 100.8 against
   Singapore's 23.3. The two families invert, for the same reason Chicago outranks Delhi on
   `heatwave-3b`.

5. 599 OF 65,797 CELLS ARE MASKED FOR THE RELATIVE RUNGS ONLY (0.91%), where the pooled
   baseline holds fewer than 840 wet days (2 per model-year across 14 GCMs x 3 scenarios x
   10 years; user decision 2026-08-17, deliberately liberal). The median cell carries 36,255
   pooled wet days -- a 43x margin -- so the rule bites only in the hyper-arid interior. The
   ABSOLUTE rungs and the intensities are NOT masked and stay global, so an arid site still
   receives R20mm and Rx1day.

6. THE 14 GCMs ARE TWO MORE THAN THE HEAT LADDER HAS, and that is a declared difference, not
   an oversight. CESM2-WACCM and IITM-ESM publish daily `pr` but no `tasmax`, so a site's
   rainfall and heat numbers come from different model sets (user decision 2026-08-17:
   prefer 14 over 12). All 14 are bias-adjusted to W5E5, so their 2020s climatologies are
   constrained to the same observations; the divergence is in the trend.

7. IITM-ESM/ssp370 ENDS AT 2098, verified against the server (its final published chunk is
   `2091_2098`), so that member contributes 9 years to the 2090s pool instead of 10. Every
   member covers 2020-2029 in full, so the SHARED BASELINE IS UNIFORM and only the final
   decade is affected; a pooled statistic absorbs it.

8. `isfinite` IS NOT A MASK. The interim fields are finite over the whole globe including
   ocean, exactly like the temperature counts, so the landseamask does all the work. The
   relative rungs additionally carry -1 as the masked-cell fill, converted to NaN on read --
   a -1 treated as a count would be a silent negative day count.

Run:
    .venv/bin/python3 scripts/process_prthresh.py --plan --metric R20mm
    .venv/bin/python3 scripts/process_prthresh.py --run --jobs 8
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

INTERIM = Path("data/interim/prthresh")
PROCESSED = Path("data/processed")
MASK_PATH = Path("data/masks/ISIMIP3b_landseamask_no-ant.nc")
BASELINE_FILE = INTERIM / "baseline_percentiles.nc"

DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2015, 2099
SLOPE_PER_DECADE = 10.0
TWO_TIER_ZERO_THRESHOLD = 0.02
CEILING_DAYS = 364.0
SATURATION_REPORT_THRESHOLD = 0.05
SPARSITY_REPORT_THRESHOLD = 0.50
SLOPE_MEM_BUDGET_BYTES = 400 * 1024 * 1024

#: metric -> (folder stem, hazard label, measure, units, central branch, relative_baseline)
METRICS = {
    "R10mm":   ("pluvial-isimip3b", "Heavy precipitation",
                "Days per year with precipitation at or above 10 mm",
                "days per year", "mean", False),
    "R20mm":   ("pluvial-isimip3b", "Heavy precipitation",
                "Days per year with precipitation at or above 20 mm",
                "days per year", "mean", False),
    "R50mm":   ("pluvial-isimip3b", "Heavy precipitation",
                "Days per year with precipitation at or above 50 mm",
                "days per year", "mean", False),
    "R100mm":  ("pluvial-isimip3b", "Heavy precipitation",
                "Days per year with precipitation at or above 100 mm",
                "days per year", "mean", False),
    "R95pD":   ("pluvial-isimip3b", "Heavy precipitation (local extreme)",
                "Days per year above this location's own 2020s wet-day 95th percentile",
                "days per year", "mean", True),
    "R99pD":   ("pluvial-isimip3b", "Heavy precipitation (local extreme)",
                "Days per year above this location's own 2020s wet-day 99th percentile",
                "days per year", "mean", True),
    "wetdays": ("pluvial-isimip3b", "Wet days",
                "Days per year with precipitation at or above 1 mm",
                "days per year", "mean", False),
    "Rx1day":  ("pluvial-isimip3b", "Extreme daily rainfall",
                "Annual maximum 1-day precipitation", "mm", "median", False),
    "Rx5day":  ("pluvial-isimip3b", "Extreme multi-day rainfall",
                "Annual maximum 5-day precipitation total", "mm", "median", False),
    "prcptot": ("pluvial-isimip3b", "Annual precipitation",
                "Annual total precipitation from wet days", "mm", "median", False),
}

#: (name, lat, lon, climate descriptor, {metric: expectation}). Printed per metric, only
#: where that metric is diagnostic -- one string shared across ten metrics is not
#: falsifiable, which is the mistake the heat ladder's first version made.
REFERENCE_SITES = [
    ("Cherrapunji", 25.30,  91.70, "wettest inhabited place",
     {"prcptot": ">5000", "R20mm": ">60", "Rx1day": ">150", "R100mm": ">2"}),
    ("Mumbai",      19.08,  72.88, "monsoon coast",
     {"prcptot": ">2000", "R20mm": ">25", "Rx1day": ">100"}),
    ("Singapore",    1.35, 103.82, "equatorial maritime",
     {"prcptot": ">1800", "wetdays": ">180", "R95pD": "HIGH -- narrow distribution"}),
    ("Manaus",      -3.12, -60.02, "humid tropics",
     {"prcptot": ">1700", "wetdays": ">140"}),
    ("London",      51.51,  -0.13, "temperate maritime",
     {"prcptot": "500-900", "wetdays": ">100", "R50mm": "~0"}),
    ("Seattle",     47.61,-122.33, "wet temperate",
     {"prcptot": ">800", "wetdays": ">120"}),
    ("Phoenix AZ",  33.45,-112.07, "hot desert",
     {"prcptot": "<500", "R20mm": "<6", "R50mm": "~0"}),
    ("Cairo",       30.04,  31.24, "hyper-arid",
     {"prcptot": "<80", "R20mm": "~0", "R95pD": "~0 or MASKED"}),
]

_CUBE = None


def log(msg=""):
    print(msg, flush=True)


def load_mask():
    """Land cells. Read through .filled(np.nan): a masked array's .data holds _FillValue,
    and `np.asarray(mask) > 0.5` marked all 259,200 cells as land the first time."""
    with Dataset(MASK_PATH) as ds:
        raw = ds.variables["mask"][:]
        arr = np.asarray(raw.filled(np.nan) if np.ma.isMaskedArray(raw) else raw,
                         dtype="f8").squeeze()
        mlat = np.asarray(ds.variables["lat"][:])
    return (np.isfinite(arr) & (arr > 0.5)), mlat


def discover_members():
    out: dict[str, dict[str, Path]] = {}
    for p in sorted(INTERIM.glob("*_pr.nc")):
        if p.name.startswith("SMOKE_"):
            continue
        with Dataset(p) as ds:
            out.setdefault(ds.scenario, {})[ds.gcm] = p
    return out


def make_pct_fn(baseline_flat):
    """Two-tier percentile-of-score against the 2020s baseline. More rain = more risk.

    Every metric here is materially zero-inflated somewhere -- a depth nobody reaches is a
    structural zero, not a missing value -- so the tier is chosen from the measured baseline
    rather than assumed. No inversion: on every metric in this family a higher value is a
    larger rainfall hazard.
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
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    per_cell = 4 * pairs * 4
    budget = SLOPE_MEM_BUDGET_BYTES // max(jobs, 1)
    return max(4, min(512, int(budget // max(per_cell, 1))))


def _slope_block(task):
    s, e, years, decade, chunk = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, BASELINE_DECADE,
                           window_years=WINDOW_YEARS, chunk_cells=chunk)
    return s, e, res["ols_slope"], res["sen_slope"]


def constant_cells(cube, years, decade):
    """Cells whose every finite observation in the expanding window is identical.

    Both slopes are EXACTLY 0 there by algebra -- every pairwise difference is 0, and the OLS
    numerator sum((x-xbar)(y-ybar)) is 0 -- so they are filled directly instead of fitted.
    Theil-Sen is quadratic in stacked points and this family is dominated by cells that never
    reach their threshold (R100mm is 94.3% zeros). Verified bit-identical to the fitted path
    on the temperature ladder. Cells with too few observations are EXCLUDED so
    expanding_slopes can apply MIN_SLOPE_OBS and return NaN, a different answer from 0.
    """
    mask = (years >= BASELINE_DECADE) & (years <= decade + WINDOW_YEARS - 1)
    flat = cube[:, mask, :].reshape(-1, cube.shape[2])
    fin = np.isfinite(flat)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        lo = np.nanmin(flat, axis=0)
        hi = np.nanmax(flat, axis=0)
    return (fin.sum(axis=0) >= MIN_SLOPE_OBS) & np.isfinite(lo) & (hi == lo)


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
    _CUBE = sub
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


def site_indices(lats, lons, metric):
    out = []
    for name, la, lo, climate, expects in REFERENCE_SITES:
        lo360 = lo % 360 if lons.max() > 180 else lo
        out.append((name, int(np.argmin(np.abs(lats - la))),
                    int(np.argmin(np.abs(lons - lo360))), climate, expects.get(metric)))
    return out


def saturation_text(metric, sat_by_panel, is_count):
    """Only meaningful for the day counts; an mm-valued metric has no calendar ceiling."""
    if not is_count:
        return ("not applicable -- this metric is a depth in mm, not a day count, so it has "
                "no calendar ceiling to saturate against.")
    worst = max(sat_by_panel.items(), key=lambda kv: kv[1])
    if worst[1] < SATURATION_REPORT_THRESHOLD:
        return (f"none -- the highest share of pooled observations at the "
                f"{CEILING_DAYS:.0f}-day calendar ceiling is {worst[1]:.2%} "
                f"({worst[0][0]} {worst[0][1]}s). Both slope estimators are readable.")
    return (
        f"MATERIAL CENSORING. {worst[1]:.1%} of pooled (year x member) observations sit at the "
        f"{CEILING_DAYS:.0f}-day calendar ceiling at {worst[0][0]} {worst[0][1]}s -- the cell "
        f"is above the threshold on essentially every day of the year and CANNOT go higher. "
        f"There BOTH slope estimators are censored toward zero and AGREE, so the dual-slope "
        f"disagreement rule gives no warning: agreement is ambiguous between 'no trend' and "
        f"'maximally exposed, permanently'. Read `median` alongside any slope.")


def sparsity_text(metric, frac_zero):
    if frac_zero < SPARSITY_REPORT_THRESHOLD:
        return (f"none -- {frac_zero:.1%} of land is exactly zero on the 2020s baseline, "
                f"below the {SPARSITY_REPORT_THRESHOLD:.0%} at which sparsity is reported.")
    return (
        f"SPARSE LAYER. {frac_zero:.1%} of land never reaches this depth on the 2020s "
        f"baseline, so a majority of sites read exactly 0 and tie at percentile tier 1. That "
        f"is the measure working as defined, NOT a missing value and NOT an extraction "
        f"failure -- most of the world does not see rainfall this heavy today. Two "
        f"consequences: a portfolio can read entirely zero and still gain exposure later "
        f"(this metric grows fastest where it starts at zero), and the percentile cannot "
        f"discriminate between zero sites, which are tied rather than ranked. For a screening "
        f"measure that discriminates across the whole land surface, use a lower depth "
        f"(R10mm/R20mm) or the intensity metrics (Rx1day), which are non-zero everywhere.")


def resolution_text(metric, adj_p99, adj_max, units):
    return (
        f"SCREENING ONLY -- 0.5 deg (~55 km) SUPPORT AGAINST A HAZARD THAT TURNS ON TERRAIN "
        f"AND ON DRAINAGE. A site's value is its cell's value, and this metric's own baseline "
        f"shows ADJACENT cells differing by up to {adj_p99:.1f} {units} at the 99th percentile "
        f"(max {adj_max:.1f}). Rainfall is more spatially variable than temperature, not less: "
        f"orographic gradients concentrate rain on windward slopes over a few kilometres, and "
        f"convective cells that drive urban surface flooding are often smaller than a single "
        f"grid box, so a cell mean systematically UNDERSTATES the peak intensity anywhere "
        f"inside it. This layer RANKS which sites deserve investigation; it cannot support a "
        f"drainage design, a depth-duration-frequency curve, or a site-specific return period.")


def process_metric(metric, members, mask, lats, lons, jobs, dry=False):
    stem, hazard, measure, units, central, relative = METRICS[metric]
    out_dir = PROCESSED / f"{stem}_{metric}_annual"
    scenarios = sorted(members)
    LAT, LON = mask.shape
    land_idx = np.flatnonzero(mask.ravel())
    n_land = land_idx.size
    is_count = units == "days per year"

    log("=" * 78)
    log(f"{metric}: {measure}")
    log(f"  -> {out_dir}   [{units}, decadal statistic = "
        f"{'pooled_mean_zero_inflated' if central == 'mean' else 'pooled_median'}"
        f"{', RELATIVE BASELINE' if relative else ''}]")
    log("=" * 78)

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    # THE RELATIVE RUNGS MUST BE MASKED HERE, AND STAGE 1 CANNOT DO IT.
    # Stage 1 counts `hit = (value >= threshold) & usable`, so a cell with no definable
    # threshold accumulates a count of 0 -- a perfectly legitimate-looking integer. The -1
    # fill in the interim file therefore only ever marks OCEAN, never the 599 cells the
    # 840-wet-day rule excluded. Left alone those cells publish "0 days above the local
    # 95th percentile", which is indistinguishable from a genuine zero, ranks at percentile
    # tier 1 (LOWEST risk), and is wrong: the truth is "not computable here". Caught
    # 2026-08-18 after R95pD's first build reported 0 of 599 masked cells as NaN.
    usable = None
    if relative:
        with Dataset(BASELINE_FILE) as bl:
            usable = np.asarray(bl.variables["usable"][:], "i1").astype(bool).ravel()[
                land_idx]
        log(f"  relative rung: masking {int((~usable).sum()):,} land cells with no definable "
            f"threshold (pooled baseline wet days < 840)")

    annual, members_by_scen = {}, {}
    for s in scenarios:
        gcms = sorted(members[s])
        members_by_scen[s] = gcms
        arr = np.full((len(gcms), n_years, n_land), np.nan, np.float32)
        for i, g in enumerate(gcms):
            with Dataset(members[s][g]) as ds:
                yr = np.asarray(ds.variables["year"][:], dtype=int)
                vals = np.asarray(ds.variables[metric][:], dtype="float32")
            # -1 marks ocean in the interim file. Left as -1 it would be a silent NEGATIVE
            # day count that survives every algebraic check (lower_ci <= median holds at -1).
            vals = np.where(vals < 0, np.nan, vals)
            flatv = vals.reshape(len(yr), -1)[:, land_idx]
            if usable is not None:
                flatv[:, ~usable] = np.nan
            for k, y in enumerate(yr):
                if int(y) in y_index:
                    arr[i, y_index[int(y)]] = flatv[k]
            del vals, flatv
        annual[s] = arr
        log(f"  loaded {s}: {len(gcms)} members -> {arr.shape}")

    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("  WARNING: ensemble composition DIFFERS across scenarios; declaring "
            "members_by_scenario so QA groups by composition.")

    base_pool = np.concatenate([annual[s] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, window_years=WINDOW_YEARS, central=central)
    # What the OTHER branch would have published -- OUTPUT-SPEC requires the deviation to be
    # declared against the number it replaced, not merely asserted.
    other = "median" if central == "mean" else "mean"
    o_med, _, _ = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, window_years=WINDOW_YEARS, central=other)
    erased = int(np.count_nonzero((o_med == 0) & (b_med > 0))) if central == "mean" else 0
    erased_mean = (float(np.mean(b_med[(o_med == 0) & (b_med > 0)]))
                   if erased else 0.0)
    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    ceil_base = (float(np.nanmean(base_pool[:, bwin, :] >= CEILING_DAYS))
                 if is_count else 0.0)
    del base_pool, o_med

    b_lo = np.clip(b_lo, 0.0, None)
    b_hi = np.clip(b_hi, 0.0, 366.0 if is_count else None)
    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)

    log(f"\n  2020s shared baseline: land n={baseline_flat.size:,}  "
        f"mean={np.nanmean(b_med):.2f} {units}  exact-zero={frac_zero:.2%}  "
        f"percentile={pct_mode}")
    if central == "mean":
        log(f"  the median branch would erase {erased:,} cells "
            f"({erased / n_land:.2%} of land) whose true mean is {erased_mean:.3f} {units}")
    if is_count:
        log(f"  at the {CEILING_DAYS:.0f}-day ceiling: {ceil_base:.2%} of the baseline pool")

    grid2d = scatter(b_med, land_idx, (LAT, LON))
    adj = np.concatenate([np.abs(np.diff(grid2d, axis=0)).ravel(),
                          np.abs(np.diff(grid2d, axis=1)).ravel()])
    adj = adj[np.isfinite(adj)]
    adj_p99 = float(np.percentile(adj, 99)) if adj.size else 0.0
    adj_max = float(adj.max()) if adj.size else 0.0
    log(f"  adjacent-cell |difference| on the baseline: p99 {adj_p99:.2f}, "
        f"max {adj_max:.2f} {units}")

    # A masked relative cell must be NaN, never 0 -- assert it rather than trust it, because
    # the failure mode is a legitimate-looking integer that no contract check rejects.
    if relative:
        n_bad = int(np.count_nonzero(np.isfinite(b_med[~usable])))
        assert n_bad == 0, (
            f"{metric}: {n_bad} of {int((~usable).sum())} cells with no definable threshold "
            f"carry a finite baseline value -- they would publish as a real count")
        log(f"  verified: all {int((~usable).sum()):,} unmaskable cells are NaN, not 0")

    grid_lookup = np.full(LAT * LON, -1, np.int64)
    grid_lookup[land_idx] = np.arange(n_land)
    log(f"\n  Reference sites, 2020s baseline, {metric} ({units}). An expectation is shown "
        f"only where this metric is diagnostic:")
    for name, i, j, climate, expect in site_indices(lats, lons, metric):
        k = grid_lookup[i * LON + j]
        val = "OCEAN/masked" if k < 0 else (
            "MASKED (relative)" if not np.isfinite(b_med[k]) else f"{b_med[k]:9.1f}")
        flag = f"expect {expect}" if expect else "-"
        log(f"    {name:<13} {val:>18}   {climate:<24} {flag}")

    if dry:
        log("\n  --plan: stopping before the slope stage.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    sat_by_panel, sen_zero_by_panel, outputs = {}, {}, {}

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
                    cube, years, d, window_years=WINDOW_YEARS, central=central)
                lo = np.clip(lo, 0.0, None)
                hi = np.clip(hi, 0.0, 366.0 if is_count else None)
                pc = pct(med)

            sl = compute_slopes(cube, years, d, chunk, jobs, n_land)
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
            n_mod = n_mem.copy()   # one GCM is one model; no impact model in the chain

            for key, vec in (("median", med), ("lower_ci", lo), ("upper_ci", hi),
                             ("percentile", pc),
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            if is_count:
                sat_by_panel[(s, d)] = float(np.nanmean(cube[:, win, :] >= CEILING_DAYS))
            if d != BASELINE_DECADE:
                o, sn = sl["ols_slope"], sl["sen_slope"]
                act = np.isfinite(o) & np.isfinite(sn) & ((o != 0.0) | (sn != 0.0))
                sen_zero_by_panel[(s, d)] = (
                    float(np.mean(sn[act] == 0.0)) if act.any() else np.nan)
            log(f"    {d}s: mean={np.nanmean(out['median'][i]):8.2f} {units}  "
                + ("slopes=NaN (baseline)" if d == BASELINE_DECADE else
                   f"ols={np.nanmean(out['ols_slope'][i]):+.3f}  "
                   f"sen={np.nanmean(out['sen_slope'][i]):+.3f} /dec  "
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

    # Rung-level slope recommendation, decided from the FAILURE MODES (see write_scenario).
    finite = [v for v in sen_zero_by_panel.values() if np.isfinite(v)]
    worst_sen0 = max(finite) if finite else float("nan")
    best_sen0 = min(finite) if finite else float("nan")
    log(f"\n  {metric}: sen_slope == 0 on {best_sen0:.1%}-{worst_sen0:.1%} of active cells "
        f"-> recommended_slope = ols_slope (uniform member coverage; see rationale)")
    for s in scenarios:
        write_scenario(outputs[s], metric, s, out_dir, lats, lons, members_by_scen[s],
                       members_by_scen, scenarios, uniform, pct_mode, frac_zero, erased,
                       erased_mean, n_land, sat_by_panel, ceil_base, best_sen0, worst_sen0,
                       adj_p99, adj_max, central, other, is_count, relative)
    del outputs
    write_members(annual, members_by_scen, scenarios, years, bwin, metric, out_dir,
                  lats, lons, land_idx, (LAT, LON), units)


def write_scenario(out, metric, s, out_dir, lats, lons, mem, members_by_scen, scenarios,
                   uniform, pct_mode, frac_zero, erased, erased_mean, n_land, sat_by_panel,
                   ceil_base, best_sen0, worst_sen0, adj_p99, adj_max, central, other,
                   is_count, relative):
    stem, hazard, measure, units, _c, _r = METRICS[metric]
    stat_name = ("pooled_mean_zero_inflated" if central == "mean" else "pooled_median")

    if central == "mean":
        rationale = (
            f"pooled_mean_zero_inflated, a DECLARED deviation from pooled_median taken on "
            f"measurement. The median branch would publish exactly 0 for {erased:,} land "
            f"cells ({erased / n_land:.2%}) whose true expected value is {erased_mean:.3f} "
            f"{units} -- places that DO experience the event, reported as never experiencing "
            f"it. Going the other way is nearly free where the median works (on R10mm the "
            f"2020s pooled median is 13.00 days against a mean of 13.15, and on wetdays "
            f"86.00 against 86.32). Applied to all SEVEN day-count metrics so a reader "
            f"comparing R20mm with R50mm is reading one statistic. The three mm-valued "
            f"metrics (Rx1day, Rx5day, prcptot) keep pooled_median: they are continuous and "
            f"essentially zero-free, so the median erases nothing there and no deviation is "
            f"warranted. Baseline exact-zero share {frac_zero:.2%}.")
    else:
        rationale = (
            f"pooled_median, OUTPUT-SPEC's default, retained ON MEASUREMENT rather than by "
            f"default. This metric is a continuous depth in mm with a baseline exact-zero "
            f"share of {frac_zero:.2%}, and the median erases NOTHING -- zero cells have a "
            f"median of 0 alongside a positive mean. The seven day-count metrics from the "
            f"same ingest DO take pooled_mean_zero_inflated, because there the median erases "
            f"up to 53.5% of land (R50mm). The split follows the measured nature of each "
            f"field, and the two families are already in different units.")

    attrs = {
        "variable": metric,
        "scenario": s,
        "long_name": measure,
        "units": units,
        "spatial_resolution_degrees": round(float(abs(lats[1] - lats[0])), 6),
        "output_spec": "OUTPUT-SPEC.md",
        "decadal_statistic": stat_name,
        "field_nature": (
            "count of days per year, integer-valued on [0, 366]; CONTINUOUS for statistical "
            "purposes but zero-inflated -- a depth nobody reaches is a STRUCTURAL zero, not "
            "a missing value" if is_count else
            "continuous non-negative depth in mm; effectively zero-free over land"),
        "value_note": (
            f"median = {'MEAN' if central == 'mean' else 'MEDIAN'} over the pooled "
            f"(year x member) sample inside the decade window across {len(mem)} ISIMIP3b "
            f"bias-adjusted GCMs. The variable is named `median` by the OUTPUT-SPEC "
            f"contract; the estimator is given by decadal_statistic. Source: daily `pr` "
            f"(kg m-2 s-1, converted x86400 to mm/day -- scale MEASURED from values, not "
            f"from the declared unit)."),
        "ci_definition": (
            f"lower_ci/upper_ci = "
            + ("mean -/+ 1 SD" if central == "mean" else "25th/75th percentile (IQR)")
            + " of the same pooled (year x member) sample, clipped at 0"
            + (" and 366" if is_count else "") + ". The band carries BOTH interannual "
            "variability and inter-model spread; it is not a pure model-spread band."),
        "slope_definition": (
            "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both are fitted "
            "over an EXPANDING window from the start of the 2020s baseline through the end "
            "of the target decade, stacking every (year, member) observation as an "
            "independent point. The baseline panel is NaN. The two fail in OPPOSITE regimes, "
            "so their disagreement flags a cell whose trend is not robust."),
        "slope_units": f"{units} per decade",
        # THE SLOPE FOLLOWS THE SAME BOUNDARY AS THE STATISTIC, and both follow measurement.
        # ols's failure mode needs uneven coverage, which does not exist here (14 members in
        # every cell), so ols is SAFE everywhere -- but safe is not best. Sen collapses on
        # the zero-inflated day counts (41.8%-100% of active cells at exactly zero) and does
        # NOT collapse on the mm-valued metrics (0.0%-1.5%), where it is also structurally
        # the better estimator: Rx1day/Rx5day are ANNUAL MAXIMA, a heavy-tailed series in
        # which one freak year drags a least-squares fit and a median of pairwise slopes
        # shrugs. Counts -> mean + ols; mm -> median + sen.
        "recommended_slope": "ols_slope" if is_count else "sen_slope",
        "recommended_slope_rationale": (
            (f"ols_slope, decided by which estimator CAN fail on this layer rather than by a "
             f"threshold. (1) ols_slope's documented failure mode requires UNEVEN member "
             f"coverage, and coverage here is {len(mem)} members in every cell -- so the bias "
             f"has no mechanism. (2) sen_slope collapses to exactly 0 where most year-pairs "
             f"are tied, measured on {best_sen0:.1%}-{worst_sen0:.1%} of ACTIVE cells (either "
             f"slope non-zero) across every post-baseline panel of all {len(scenarios)} "
             f"scenarios. The three mm-valued metrics from this same ingest take sen_slope "
             f"instead, because there it does not collapse.")
            if is_count else
            (f"sen_slope, on measurement. Theil-Sen has NOT collapsed on this metric: it is "
             f"exactly zero on only {best_sen0:.1%}-{worst_sen0:.1%} of ACTIVE cells across "
             f"every post-baseline panel of all {len(scenarios)} scenarios -- cleaner than "
             f"either other sen_slope layer in the product (conifer-npp 2.1%, csoil 6.8%), so "
             f"the robust estimator is the read. It is also structurally right here: an "
             f"annual maximum is a heavy-tailed extreme-value series, where one freak year "
             f"drags a least-squares fit and a median of pairwise slopes does not. CONTRAST "
             f"the seven day-count metrics from the same ingest, which take ols_slope because "
             f"sen collapses on 41.8%-100% of their active cells.")
        ) + (" BOTH slopes are in this file and their DISAGREEMENT remains the per-cell "
             "signal that a trend is not robust."),
        "percentile_baseline": (
            f"{pct_mode}: each cell's value is ranked against the 2020s ensemble land "
            f"distribution. Zeros take tier 1; positive values rank against the NON-ZERO "
            f"baseline into [2, 100]."),
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
        "impact_models": ("none -- derived directly from bias-adjusted climate forcing, with "
                          "no impact model in the chain, so the confidence band carries GCM "
                          "spread and interannual variability only"),
        "normalization": (
            "none. Every member reports the same physical quantity in the same unit on the "
            "same grid and calendar, so members are equal-weighted raw (model democracy)."),
        "spatial_smoothing": (
            f"none. {len(mem)} members x {WINDOW_YEARS} years = {len(mem) * WINDOW_YEARS} "
            f"draws per cell-decade, and the input is bias-adjusted W5E5 forcing that is "
            f"already spatially coherent."),
        "land_mask": (
            f"ISIMIP3b W5E5 landseamask_no-ant.nc, {n_land:,} cells, Antarctica EXCLUDED. "
            f"NOTE the interim field is FINITE OVER THE WHOLE GLOBE INCLUDING OCEAN -- "
            f"`isfinite` is not a mask on this layer, the landseamask is."),
        "decadal_statistic_rationale": rationale,
        "source_dataset": (
            "ISIMIP3b bias-adjusted (W5E5v2) daily `pr`, "
            "InputData/climate/atmosphere/bias-adjusted/global/daily (5 group-I GCMs) and "
            "SecondaryInputData/... (9 extended GCMs), constructed into annual statistics by "
            "scripts/pr_baseline_percentiles.py (stage 0, per-cell wet-day percentiles) + "
            "scripts/download_reduce_prthresh_isimip3b.py (stage 1). ~764 GB streamed, each "
            "chunk sha512-verified against the publisher sidecar and deleted after "
            "reduction; see data/interim/prthresh/download_provenance.csv."),
        "resolution_caveat": resolution_text(metric, adj_p99, adj_max, units),
        "interpretation_caveat": (
            "A RAINFALL HAZARD, NOT A FLOOD OUTCOME. (1) Pluvial (surface-water) flooding is "
            "rainfall intensity against DRAINAGE CAPACITY, and drainage is not in this data "
            "at all. A well-drained urban cell and an undrained rural one with the same "
            "rainfall read identically. This layer ranks rainfall hazard; it never states "
            "that a site floods. (2) It is also NOT river flooding -- for channel overflow "
            "use the `flood-3b-*` CaMa-Flood layers, which route water through a hydrological "
            "model. The two are different hazards and a site can be exposed to either "
            "without the other. (3) There is NO IMPACT MODEL in the chain, so the confidence "
            "band carries GCM spread and interannual variability only, with no impact-model "
            "structural uncertainty. (4) THE 14 GCMs ARE NOT 14 INDEPENDENT MODELS: "
            "CNRM-CM6-1 and CNRM-ESM2-1 share ARPEGE-Climat, NEMO and ISBA-CTRIP, and "
            "KACE-1-0-G runs the HadGEM3 atmosphere it shares with UKESM1-0-LL. n_models "
            "counts GCMs; read the CI as spanning correlated members. (5) THIS ENSEMBLE IS "
            "NOT THE HEAT LADDER'S: CESM2-WACCM and IITM-ESM publish `pr` but no `tasmax`, "
            "so a site's rainfall and temperature numbers come from different model sets. "
            "(6) Sub-daily intensity is absent -- a 30-minute cloudburst that overwhelms a "
            "storm drain can sit inside an unremarkable daily total."),
        "created_by": "scripts/process_prthresh.py",
    }

    sat = saturation_text(metric, sat_by_panel, is_count)
    spa = sparsity_text(metric, frac_zero)
    if sat.startswith("MATERIAL"):
        attrs["saturation_caveat"] = sat
    else:
        attrs["saturation_measured"] = sat
    if spa.startswith("SPARSE"):
        attrs["sparsity_caveat"] = spa
    else:
        attrs["sparsity_measured"] = spa

    if relative:
        attrs["relative_baseline_note"] = (
            "THIS METRIC SCORES DEPARTURE FROM THE SITE'S OWN HISTORY, NOT ABSOLUTE RAINFALL. "
            "It counts days exceeding the 95th (or 99th) percentile of that cell's OWN 2020s "
            "wet-day rainfall, so the threshold is different in every location: 3.3 mm/day in "
            "the driest cells and 165.2 in the wettest, median 19.0. A high score means 'rain "
            "here is moving furthest from what this place is used to', NEVER 'it rains hard "
            "here'. THE RANKING INVERTS AGAINST THE ABSOLUTE METRICS, measured on this "
            "ensemble: Singapore reads 14.9 days a year above its own p95 against "
            "Cherrapunji's 6.8, while on the absolute rung Cherrapunji reads 100.8 days above "
            "20 mm against Singapore's 23.3. Both are correct and they answer different "
            "questions. BY CONSTRUCTION the 2020s panel is a control, not a hazard map: every "
            "cell sits at ~5% (or ~1%) of its own wet days in the baseline decade by "
            "definition, so the baseline is essentially a map of how often it rains and the "
            "SIGNAL IS THE CHANGE from it. Pair this with the absolute metrics (R20mm, "
            "Rx1day) before drawing any conclusion, and never quote it as a rainfall depth.")

    ds_out = xr.Dataset(
        {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
        coords={"decade": DECADES, "lat": lats, "lon": lons},
        attrs=attrs,
    )
    enc = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
               "_FillValue": np.float32(np.nan)} for k in out}
    path = out_dir / f"{metric}_{s}_processed.nc"
    ds_out.to_netcdf(path, encoding=enc)
    log(f"    wrote {path.name} ({path.stat().st_size / 1e6:.1f} MB)")


def write_members(annual, members_by_scen, scenarios, years, bwin, metric, out_dir,
                  lats, lons, land_idx, shape, units):
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
        attrs={"variable": metric, "units": units,
               "member_field": (f"mean over the {BASELINE_DECADE}s baseline decade, pooled "
                                "across all scenarios"),
               "note": "Diagnostic only -- not part of the OUTPUT-SPEC contract."},
    )
    p = out_dir / f"{metric}_members.nc"
    ds.to_netcdf(p, encoding={"value": {"dtype": "float32", "zlib": True, "complevel": 4,
                                        "_FillValue": np.float32(np.nan)}})
    log(f"  wrote per-member diagnostic {p.name} ({len(names)} members)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--metric", action="append", default=None,
                    help="process only these metrics (default: all ten)")
    ap.add_argument("--jobs", type=int, default=max(1, (mp.cpu_count() or 2) - 2))
    a = ap.parse_args()

    members = discover_members()
    if not members:
        log(f"ERROR: no interim files found in {INTERIM}")
        return 1
    if not BASELINE_FILE.exists():
        log(f"ERROR: {BASELINE_FILE} missing -- the relative rungs were built against it and "
            "its provenance belongs with the layers. Run stage 0.")
        return 1
    mask, mlat = load_mask()
    first = next(iter(next(iter(members.values())).values()))
    with Dataset(first) as ds:
        lats = np.asarray(ds.variables["lat"][:])
        lons = np.asarray(ds.variables["lon"][:])
    if np.allclose(lats, mlat[::-1]):
        mask = mask[::-1]
        log("NOTE: mask latitude axis flipped to match the interim grid")
    elif not np.allclose(lats, mlat):
        log("ERROR: mask and interim files are on different grids")
        return 1

    metrics = a.metric or list(METRICS)
    unknown = [m for m in metrics if m not in METRICS]
    if unknown:
        log(f"ERROR: unknown metric(s) {unknown}; known: {list(METRICS)}")
        return 1

    log(f"interim       {INTERIM}")
    log(f"land mask     {MASK_PATH.name} -- {int(mask.sum()):,} cells (Antarctica EXCLUDED)")
    for s in sorted(members):
        log(f"scenario      {s}: {len(members[s])} members")
    log(f"metrics       {metrics}")
    log(f"jobs          {a.jobs}")
    if not (a.run or a.plan):
        log("\nNothing to do. Pass --plan or --run.")
        return 0
    for m in metrics:
        process_metric(m, members, mask, lats, lons, a.jobs, dry=not a.run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
