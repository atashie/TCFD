"""Process `leh` (ISIMIP2b/RCP heatwave exposure) into the TCFD output contract.

`leh` (Lange2020, ISIMIP2b DerivedOutputData) is the heatwave member of the Lange 2020
extreme-event exposure family. Each cell carries a per-year flag for whether it saw a
heatwave.

WHY THIS LAYER EXISTS ALONGSIDE `heatwave-3b`: IT IS A DIFFERENT INDEX, NOT AN OLDER
VERSION. This one is the only ISIMIP heatwave product with an ABSOLUTE, health-anchored
criterion.

    ISIMIP2b  Lange2020 / `HWMId-humidex` -> HWMId (relative) AND Humidex >= 45 (absolute)
    ISIMIP3b  Zantout2025 / `HWMID-NONE`  -> HWMId ONLY; the "NONE" is the humidity term

Lange et al. 2020 (Earth's Future, doi:10.1029/2020EF001616) combine a RELATIVE indicator
with an ABSOLUTE one so that the events counted "would also adversely affect human health".
The reference implementation states the combination explicitly: land is exposed "if BOTH a
relative indicator based on temperature (Russo et al. 2015, 2017) AND an absolute indicator
based on temperature and relative humidity (Masterton & Richardson, 1979) exceed their
respective threshold values". The Humidex threshold is 45, inside Environment Canada's
40-45 "great discomfort; avoid exertion" band.

READ THAT AS AN INTERSECTION, NOT AS AN ABSOLUTE INDEX. The Humidex gate makes the counted
events health-relevant, but the relative HWMId gate still has to open first, so `leh` ADDS
a health filter on top of a per-cell-distribution framing rather than replacing it. Both
framings therefore apply to this layer, and both are declared below.

MEASURED data nature (GUARDRAILS 9 -- from the values), all 8 members,
`scripts/check_leh_nature.py` 2026-08-14:

  * STRICTLY BINARY {0,1}: exactly 2 unique values per member, range [0,1]. NOT inherited --
    the sibling `let` in this same publication is a CONTINUOUS fraction in [0, 0.952).
    -> pooled MEAN = heatwave frequency.
  * Annual exact-zero share 97.93%-99.29% per member. Extremely sparse; see below.
  * 4 members per scenario, composition IDENTICAL across rcp26/rcp60, so the shared 2020s
    baseline is valid. soc/sens uniformly nosoc/co2. No rcp85 exists for this product.

THE PRODUCT HAS NO LAND MASK AND THIS PROCESSOR MUST SUPPLY ONE. Every one of the 259,200
grid cells is finite in every member -- the publisher ZERO-FILLS ocean rather than masking
it, exactly like its 3b cousin `cropfailure`. So `isfinite` carries no footprint here, and
a processor that trusted it would publish the entire ocean as a confident zero and rank it
in the percentile distribution. The ISIMIP2b land-sea mask is applied explicitly:
67,420 land cells.

  NOTE THE MASK VARIABLE NAME DIFFERS BY ROUND: ISIMIP2b publishes `LSM` (with a singleton
  time axis), ISIMIP3b publishes `mask`. Assuming `mask` here silently disabled the land
  denominator during the value check and quoted every sparsity share against the whole
  globe (89.1% "silent") instead of against land (59.7%).

TIME DECODING: the axis is `years since 1661-1-1` on a **360_day** calendar, and xarray
CANNOT decode that combination -- `decode_times=True` raises
"unable to decode time units 'years since 1661-1-1 00:00:00' with calendar '360_day'".
This is the THIRD time convention across three publications of the same concept
(Heinicke2026 `days since 1601-01-01`; Zantout2025 `days since 2015-01-01`, mid-year
stamps, alternating 365/365.5 steps). The convention is per-PUBLICATION, always.
Here the values are integer years since an epoch, so the epoch is parsed FROM THE UNITS
STRING and the resulting span is asserted against the filename -- never hardcoded, and
never divided by 365.

THE DEFINING PROPERTY OF THIS LAYER: THE ABSOLUTE GATE IS SILENT OVER MOST OF THE
EXTRATROPICS, AND THAT IS THE INDEX WORKING, NOT A DEFECT.

Measured over 2006-2099 across both RCPs and all 4 GCMs, on the 67,420-cell land mask:
40,265 cells (59.7%) are NEVER non-zero in any member or any year. By latitude belt, the
share of land that never registers a single exposed year:

    60-90N (boreal)        99%          tropics              18%
    35-60N (temperate)     81%          23.5-35N             30%
    35-60S                 68%          23.5-35S             19%

Reference sites, mean 2006-2099 -> mean over the 2090s:
    Jacobabad 0.553->0.713   Kuwait 0.548->0.575   Lagos 0.363->0.500
    Delhi 0.134->0.225       Singapore 0.134->0.250  Phoenix 0.112->0.213
    Cordoba 0.036->0.038     Sydney 0.019->0.025
    Chicago 0.0027->0.0125   Sao Paulo 0.0027->0.0125
    PARIS, FRANKFURT AND YAKUTSK ARE EXACTLY ZERO IN EVERY MEMBER AND EVERY YEAR.

This is the opposite failure mode to `heatwave-3b`, which SATURATES (45.9% of ssp585 2090s
cells pinned at 1.0). Between them: the absolute index is silent where the relative one
saturates, and vice versa. Neither is wrong; they measure different things.

A ZERO HERE MEANS "NEVER CROSSES HUMIDEX 45", NOT "LOW HEAT RISK", and the distinction is
the whole reason `interpretation_caveat` and `sparsity_caveat` are written the way they
are. A delivery that renders this layer as a low percentile across Europe without saying
which of the two it means would be reporting a threshold artifact as a finding. Contrast
the withdrawn sugarcane layers, where the model simply did not run: here the model ran
everywhere and the threshold is genuinely not crossed.

TWO DECISIONS SPECIFIC TO THIS LAYER. Neither inherited; both measured here.

1. NO MINIMUM-MODEL MASK, AND THE RULE IS VOID RATHER THAN RELAXED -- one index model
   (hwmid-humidex), so `n_models` is 1 at every published cell and the `led` layer's
   ">= 2 impact models" rule cannot be expressed. Nor is there a GCM-count fallback: after
   the land mask is applied the four GCM members cover identical cells. `n_members` and
   `n_models` are emitted per cell so the degeneracy is visible in the file.

   CAVEAT recorded rather than masked: with one index model the CI carries interannual and
   inter-GCM spread ONLY. There is NO inter-model structural uncertainty in it -- and on
   this layer that understates uncertainty unusually badly, because the Humidex threshold
   is a single hard cut with no second formulation to disagree with it.

2. SPATIAL SMOOTHING: see `spatial_smoothing` in the output attributes -- decided from the
   per-member contact sheet and the measured spatial structure, NOT inherited from `let`
   (4 members, L=2.5) merely because the member count matches.

Percentile: ranked against the shared 2020s land distribution, `higher_is_worse`, two-tier
because the baseline is heavily zero-inflated. Tier 1 therefore holds the large silent
majority -- correct, and the reason `percentile_tier1_share` is published.

Expect `sen_slope` to be exactly 0 nearly everywhere: on a binary field every pairwise
difference is -1, 0 or +1 and same-valued pairs dominate. READ `ols_slope`.

Usage:
    python scripts/process_leh_isimip2b.py [--scenarios rcp26] [--skip-slopes]
                                           [--limit-cells N] [--members-only] [--jobs N]
                                           [--smooth-l L]
"""

import argparse
import glob
import multiprocessing as mp
import os
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.decadal_stats import (  # noqa: E402
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

VAR = "leh"
OUT_VAR = "leh"
LAYER_ID = "heatwave-isimip2b_leh_annual"
#: ISIMIP2b runs from 2006, so a full 2010s panel exists and the baseline sits at index 1 --
#: unlike the ISIMIP3b layers, whose baseline IS index 0. `test_shared_baseline.py` locates
#: it from the `baseline_decade` attribute rather than assuming index 0.
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2010, 2099

CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False
SLOPE_PER_DECADE = 10.0
MIN_MODELS = 1
STAT_NAME = "pooled_mean_boolean"
MODEL_FAMILY = {}

SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2
_CUBE = None


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) from a Lange2020 filename.

    lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099.nc4

    ISIMIP2b Lange2020 carries a LEADING publication token AND the `.nc4` extension. Read
    from the END so the parse is stable across the round's grammar differences.
    """
    p = os.path.basename(fpath).split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def decode_years(ds, fpath):
    """Integer years, decoded from the units epoch and CHECKED against the filename.

    xarray cannot decode `years since ... ` on a `360_day` calendar at all (it raises), so
    this is done by hand -- but the epoch is PARSED from the units string rather than
    hardcoded, and the resulting span is asserted against the filename, which is part of the
    ISIMIP grammar. If either the epoch or the span moves, this raises instead of silently
    mis-binning a decade.
    """
    units = str(ds["time"].attrs.get("units", ""))
    m = re.match(r"\s*years?\s+since\s+(-?\d+)", units)
    if not m:
        raise ValueError(f"{os.path.basename(fpath)}: expected a 'years since <epoch>' time "
                         f"axis for this publication; got {units!r}")
    epoch = int(m.group(1))
    t = np.asarray(ds["time"].values, dtype="float64")
    if not np.allclose(t, np.round(t)):
        raise ValueError(f"{os.path.basename(fpath)}: 'years since' axis carries "
                         "non-integer offsets; the convention changed")
    yrs = (epoch + np.round(t)).astype(int)

    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    y0, y1 = int(p[-2]), int(p[-1])
    if yrs.min() != y0 or yrs.max() != y1 or yrs.size != (y1 - y0 + 1):
        raise ValueError(f"{os.path.basename(fpath)}: decoded {yrs.min()}-{yrs.max()} "
                         f"({yrs.size} records) but the filename declares {y0}-{y1}")
    if np.unique(yrs).size != yrs.size:
        raise ValueError(f"{os.path.basename(fpath)}: duplicate years after decoding")
    return yrs


def load_land_mask(root, shape):
    """The ISIMIP2b land-sea mask, as a 2-D boolean.

    Required, not optional: this product zero-fills ocean, so without an explicit mask the
    entire ocean would be published as a confident zero and would enter the percentile
    reference distribution. The mask VARIABLE NAME is round-specific -- 2b is `LSM`, 3b is
    `mask` -- so it is resolved rather than assumed.
    """
    sys.path.insert(0, str(root / "scripts"))
    from utils.land_mask import get_isimip_landmask  # noqa: E402
    with xr.open_dataset(get_isimip_landmask("2b")) as lm:
        name = "LSM" if "LSM" in lm.variables else "mask"
        arr = np.squeeze(lm[name].values)
    land = np.nan_to_num(arr) > 0
    if land.shape != shape:
        raise ValueError(f"land mask {land.shape} does not match the data grid {shape}")
    return land


def load_member(fpath, land):
    """Load one member as (years, (year, lat, lon)) with ocean set to NaN by the LAND MASK.

    The file itself is finite everywhere; the mask applied here is what makes the output a
    land product.
    """
    ds = xr.open_dataset(fpath, decode_times=False)
    da = ds[VAR]
    yrs = decode_years(ds, fpath)
    vals = da.values.astype("float32")
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan
    vals[:, ~land] = np.nan

    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    return yrs[keep], vals[keep]


def smooth_annual(vals, decay_l, land):
    """5x5 exponential-decay spatial smoothing applied per member per YEAR.

    Applied to ANNUAL maps, not decadal means, so the CI and both slopes see the smoothed
    series -- smoothing the decadal mean alone would leave the CI describing an unsmoothed
    field. Weights are normalised over the FINITE neighbours only, so land cells beside the
    coast are not diluted toward zero by ocean.
    """
    off = np.arange(-2, 3)
    dy, dx = np.meshgrid(off, off, indexing="ij")
    w = np.exp(-np.sqrt(dy**2 + dx**2) / decay_l)
    out = np.full_like(vals, np.nan)
    finite = np.isfinite(vals)
    num = np.zeros_like(vals)
    den = np.zeros_like(vals)
    filled = np.where(finite, vals, 0.0)
    for k in range(5):
        for j in range(5):
            sy, sx = int(dy[k, j]), int(dx[k, j])
            num += w[k, j] * np.roll(np.roll(filled, sy, axis=1), sx, axis=2)
            den += w[k, j] * np.roll(np.roll(finite.astype(np.float32), sy, axis=1),
                                     sx, axis=2)
    good = den > 0
    out[good] = num[good] / den[good]
    out[:, ~land] = np.nan
    return out.astype(np.float32)


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline land distribution.

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros), which on this
    layer it emphatically is: a cell that never crosses the Humidex threshold across the
    pooled baseline sample is the LOWEST risk this index can express, so it maps to 1 and
    positives rank against the NON-ZERO baseline into [2,100]. That tier-1 block is large by
    construction here -- see `percentile_tier1_share` -- and it means "below the health
    threshold", not "no heat".
    """
    frac_zero = float(np.mean(baseline_flat == 0.0))
    two_tier = frac_zero >= TWO_TIER_ZERO_THRESHOLD

    def _invert(sc):
        return np.clip(101.0 - sc, 1.0, 100.0) if higher_is_better else sc

    if two_tier:
        nz_sort = np.sort(baseline_flat[baseline_flat > 0])
        n_nz = len(nz_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            res = np.ones(vals.shape, np.float32)
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="right") / n_nz
                res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
            out[fin] = _invert(res)
            return out.reshape(arr.shape)
    else:
        all_sort = np.sort(baseline_flat)
        n = len(all_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            if n > 0:
                frac = np.searchsorted(all_sort, vals, side="right") / n
                out[fin] = _invert(np.clip(100.0 * frac, 1.0, 100.0)).astype(np.float32)
            return out.reshape(arr.shape)

    return pct, ("two_tier" if two_tier else "single_tier"), frac_zero


def scatter(flat_land, land_idx, shape):
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


def slope_chunk_cells(n_members, n_years, max_pairs):
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    if max_pairs is not None:
        pairs = min(pairs, max_pairs)
    per_cell = 4 * pairs * 4
    return max(4, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def _slope_block(task):
    s, e, years, decade, baseline, window, chunk, max_pairs = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk, max_pairs=max_pairs)
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, max_pairs, jobs, n_land):
    if decade == BASELINE_DECADE or jobs <= 1:
        return expanding_slopes(cube, years, decade, BASELINE_DECADE,
                                window_years=WINDOW_YEARS, chunk_cells=chunk,
                                max_pairs=max_pairs)
    n_blocks = max(jobs * 8, 1)
    edges = np.linspace(0, n_land, n_blocks + 1).astype(int)
    tasks = [(int(a), int(b), years, decade, BASELINE_DECADE, WINDOW_YEARS, chunk, max_pairs)
             for a, b in zip(edges[:-1], edges[1:]) if b > a]
    ols = np.full(n_land, np.nan, np.float32)
    sen = np.full(n_land, np.nan, np.float32)
    ctx = mp.get_context("fork")
    with ctx.Pool(jobs) as pool:
        for s, e, o, sn in pool.imap_unordered(_slope_block, tasks):
            ols[s:e] = o
            sen[s:e] = sn
    return {"ols_slope": ols, "sen_slope": sen}


def main():
    global _CUBE
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", nargs="*", default=None)
    ap.add_argument("--limit-cells", type=int, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--members-only", action="store_true")
    ap.add_argument("--skip-slopes", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--min-models", type=int, default=MIN_MODELS)
    ap.add_argument("--smooth-l", type=float, default=None,
                    help="5x5 exponential-decay length. Omitted = NO smoothing, which is "
                         "this layer's measured default; do not copy `let`'s L=2.5 just "
                         "because the member count matches.")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    out_dir = root / "data" / "processed" / LAYER_ID
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc4")))
    if not files:
        log(f"ERROR: no {VAR} files in {raw_dir}")
        return 2

    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    write_scenarios = args.scenarios or scenarios
    bad = [s for s in write_scenarios if s not in scenarios]
    if bad:
        log(f"ERROR: unknown scenario(s) {bad}; found {scenarios}")
        return 2

    log("=" * 74)
    log(f"leh (ISIMIP2b/RCP heatwave exposure) -> TCFD contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")
    log("INDEX NOTE: HWMId AND Humidex>=45 -- an INTERSECTION. The absolute gate is silent "
        "over most of the extratropics under rcp26/rcp60; that is the index, not a defect.")

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    LAT, LON = len(lats), len(lons)

    land = load_land_mask(root, (LAT, LON))
    log(f"\nLand mask (ISIMIP2b `LSM`, applied EXPLICITLY -- the product zero-fills the "
        f"globe): {int(land.sum()):,} of {LAT * LON:,} cells")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    members_by_scen = {s: [] for s in scenarios}
    for f in files:
        members_by_scen[meta[f]["scenario"]].append(meta[f]["member"])
    for s in scenarios:
        members_by_scen[s] = sorted(members_by_scen[s])

    models = sorted({meta[f]["model"] for f in files})
    gcms = sorted({meta[f]["gcm"] for f in files})
    socs = sorted({meta[f]["soc"] for f in files})
    senss = sorted({meta[f]["sens"] for f in files})

    keep2d = land
    land_idx = np.flatnonzero(keep2d.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        land_idx = land_idx[np.linspace(0, land_idx.size - 1,
                                        args.limit_cells).astype(int)]
    n_land = land_idx.size

    t0 = time.time()
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        yrs, cube = load_member(f, land)
        if args.smooth_l:
            cube = smooth_annual(cube, args.smooth_l, land)
        flat = cube.reshape(cube.shape[0], -1)
        for k, y in enumerate(yrs):
            yi = y_index.get(int(y))
            if yi is not None:
                annual[s][slot[s][m], yi] = flat[k, land_idx]
        del cube, flat
        log(f"  loaded {info['model']:<14} {info['gcm']:<15} {s}  [{time.time() - t0:.0f}s]")
    log(f"Packed {len(files)} members over {n_land:,} land cells "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident)")

    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> is_boolean_field={boolean}")
    if boolean and args.smooth_l:
        log("  NOTE: smoothing was applied, so the packed cube is no longer strictly "
            "binary; nature was classified on the SMOOTHED values above.")
    fin0 = annual[scenarios[0]][np.isfinite(annual[scenarios[0]])]
    log(f"  annual cell-year exact-zero fraction (land only): {float((fin0 == 0).mean()):.2%}")
    del fin0
    if not boolean and not args.smooth_l:
        log("  ERROR: `leh` measured BINARY at ingest 2026-08-14 and no smoothing was "
            "applied, so a continuous read means the inputs changed -- re-run "
            "check_leh_nature.py before processing.")
        return 3
    stat_name = STAT_NAME if boolean else "pooled_mean_smoothed_boolean"
    log(f"  decadal statistic: {stat_name} -- pooled MEAN = heatwave frequency")

    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios.")
    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    base_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years[bwin], BASELINE_DECADE, boolean=True, window_years=WINDOW_YEARS)
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, CI_CAP)
    b_hi = np.clip(b_hi, CI_FLOOR, CI_CAP)

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    tier1 = float(np.mean(b_pct[np.isfinite(b_pct)] <= 1.0))
    log(f"\nShared {BASELINE_DECADE}s baseline: land n={baseline_flat.size:,}, "
        f"exact-zero {frac_zero:.2%}, percentile mode={pct_mode}, tier-1 share {tier1:.2%}, "
        f"mean frequency {np.nanmean(b_med):.4f}, max {np.nanmax(b_med):.4f}")

    mem_names = members_by_scen[scenarios[0]]
    mem_grid = np.full((len(mem_names), LAT, LON), np.nan, np.float32)
    for mi, mname in enumerate(mem_names):
        stack = [annual[s][members_by_scen[s].index(mname)][bwin]
                 for s in scenarios if mname in members_by_scen[s]]
        if not stack:
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            flat = np.nanmean(np.concatenate(stack, axis=0), axis=0)
        mem_grid[mi] = scatter(flat, land_idx, (LAT, LON))
    mem_ds = xr.Dataset(
        {"value": (["member", "lat", "lon"], mem_grid)},
        coords={"member": mem_names, "lat": lats, "lon": lons},
        attrs={
            "variable": OUT_VAR, "units": "1",
            "member_field": (f"heatwave frequency over the {BASELINE_DECADE}s baseline "
                             "decade (mean of the binary annual flag), pooled across "
                             "scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. All four "
                     "members are the SAME index model (hwmid-humidex) driven by four "
                     "different CMIP5 GCMs, so differences are GCM spread, not model "
                     "structural spread. The field is SPARSE by construction (absolute "
                     "Humidex>=45 gate): large silent areas are expected, and what to look "
                     "for is block structure, seams, and whether the ACTIVE regions agree "
                     "between GCMs."),
        },
    )
    mem_path = out_dir / f"{OUT_VAR}_members.nc"
    mem_ds.to_netcdf(mem_path, encoding={"value": {"dtype": "float32", "zlib": True,
                                                   "complevel": 4,
                                                   "_FillValue": np.float32(np.nan)}})
    log(f"  wrote per-member diagnostic {mem_path.name} ({len(mem_names)} members)")

    # ---- How many GCMs agree a cell is active AT ALL -------------------------- #
    # Not a level-agreement check: on a sparse threshold index the first question is
    # whether the members even agree the threshold is ever crossed. Measured on the
    # baseline decade, pooled across scenarios, per member.
    base_active = np.zeros((len(mem_names), n_land), bool)
    for mi, mname in enumerate(mem_names):
        for sc in scenarios:
            if mname in members_by_scen[sc]:
                blk = annual[sc][members_by_scen[sc].index(mname)][bwin]
                base_active[mi] |= (np.nan_to_num(blk) > 0).any(axis=0)
    n_agree = base_active.sum(axis=0)
    agree_tiers = {k: int((n_agree == k).sum()) for k in range(len(mem_names) + 1)}
    n_mem_tot = len(mem_names)
    active_any = int((n_agree > 0).sum())
    solo_share = (agree_tiers.get(1, 0) / active_any) if active_any else float("nan")
    agree_txt = "; ".join(f"{k}/{n_mem_tot}: {v:,} ({100.0 * v / n_land:.1f}% of land)"
                          for k, v in agree_tiers.items() if v)
    log(f"\nGCM agreement on ACTIVITY over the {BASELINE_DECADE}s baseline: {agree_txt}")
    log(f"  of the {active_any:,} cells active in any member, {solo_share:.1%} rest on a "
        f"SINGLE GCM")
    del mem_grid
    if args.members_only:
        log("\n--members-only: done.")
        return 0

    for s_drop in [s for s in scenarios if s not in write_scenarios]:
        del annual[s_drop]

    chunk = slope_chunk_cells(max(len(m) for m in members_by_scen.values()),
                              n_years, args.max_pairs)
    jobs = max(1, args.jobs)
    log(f"Theil-Sen chunk_cells={chunk}, jobs={jobs}")

    for s in write_scenarios:
        log(f"\n{'=' * 74}\nAssembling scenario {s}\n{'=' * 74}")
        mem = members_by_scen[s]
        cube = annual[s]
        _CUBE = cube
        fams = sorted({family_of(m.rsplit("_", 1)[0]) for m in mem})
        fam_idx = {f: k for k, f in enumerate(fams)}

        nd = len(DECADES)
        out = {k: np.full((nd, LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}
        silent_by_decade = {}

        for i, d in enumerate(DECADES):
            td = time.time()
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, boolean=True, window_years=WINDOW_YEARS)
                lo = np.clip(lo, CI_FLOOR, CI_CAP)
                hi = np.clip(hi, CI_FLOOR, CI_CAP)
                pc = pct(med)

            mf = med[np.isfinite(med)]
            silent_by_decade[d] = float((mf == 0.0).mean()) if mf.size else float("nan")

            if args.skip_slopes:
                sl = dict(ols_slope=np.full(n_land, np.nan, np.float32),
                          sen_slope=np.full(n_land, np.nan, np.float32))
            else:
                sl = compute_slopes(cube, years, d, chunk, args.max_pairs, jobs, n_land)

            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
            fam_present = np.zeros((len(fams), n_land), bool)
            for mi, m in enumerate(mem):
                fam_present[fam_idx[family_of(m.rsplit("_", 1)[0])]] |= present[mi]
            n_mod = fam_present.sum(axis=0).astype(np.float32)
            n_mem[n_mem == 0] = np.nan
            n_mod[np.isnan(n_mem)] = np.nan

            for key, vec in (("median", med), ("lower_ci", lo), ("upper_ci", hi),
                             ("percentile", pc),
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                warnings.filterwarnings("ignore", message="All-NaN slice encountered")
                if d <= BASELINE_DECADE:
                    slope_txt = "slopes=NaN (baseline or before)"
                else:
                    slope_txt = (f"ols={np.nanmean(out['ols_slope'][i]):+.3e} "
                                 f"sen={np.nanmean(out['sen_slope'][i]):+.3e} /dec")
                tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
                log(f"  {d}s: {tag:<15} mean={np.nanmean(out['median'][i]):.5f}  "
                    f"silent={silent_by_decade[d]:6.2%}  {slope_txt}  "
                    f"[{time.time() - td:.0f}s]")

        for i, d in enumerate(DECADES):
            if d <= BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), \
                    f"ols must be NaN at {d}s (no elapsed period from the baseline)"
                assert np.all(np.isnan(out["sen_slope"][i])), f"sen must be NaN at {d}s"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), (
                    f"{k} finite where median is NaN at {d}s ({extra.sum()} cells)"
                    " -- ocean leak")
            assert np.all(out["lower_ci"][i][med_finite]
                          <= out["median"][i][med_finite] + 1e-5), f"lower_ci > median {d}s"
            assert np.all(out["upper_ci"][i][med_finite]
                          >= out["median"][i][med_finite] - 1e-5), f"upper_ci < median {d}s"
        # The whole point of the explicit land mask: nothing may be published off it.
        for i in range(len(DECADES)):
            assert not (np.isfinite(out["median"][i]) & ~land).any(), \
                "median finite off the ISIMIP2b land mask -- the ocean zero-fill leaked"

        fin_sen = out["sen_slope"][-1][np.isfinite(out["sen_slope"][-1])]
        sen_zero = float((fin_sen == 0).mean()) if fin_sen.size else float("nan")
        silent_txt = "; ".join(f"{d}s {v:.1%}" for d, v in silent_by_decade.items())

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Land area fraction exposed to heatwave (HWMId AND Humidex>=45)",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "field_nature": "boolean",
                "value_note": (
                    "median = MEAN over the pooled (year x member) sample inside the decade "
                    f"window, across {len(mem)} members (1 index model x 4 CMIP5 GCMs), i.e. "
                    "the HEATWAVE FREQUENCY: the fraction of model-years the cell was "
                    "flagged as exposed. The raw ISIMIP field is a BINARY {0,1} annual flag "
                    "(measured: exactly 2 unique values in all 8 members). NOT inherited -- "
                    "the sibling `let` in this same publication is a continuous fraction."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clamped to [{CI_FLOOR}, {CI_CAP}]. WITH ONE "
                    "INDEX MODEL THIS CI CARRIES INTERANNUAL AND INTER-GCM SPREAD ONLY. On "
                    "this layer that understates uncertainty unusually badly: the Humidex "
                    ">= 45 criterion is a single hard threshold with no second formulation "
                    "in the ensemble to disagree with it, so the sensitivity of the whole "
                    "field to that one number is invisible in the CI. Do not compare CI "
                    "widths across layers."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope, both over "
                    "an EXPANDING window from the start of the 2020s baseline through the "
                    "end of the target decade, stacking every (year, member) observation. "
                    "The 2010s and 2020s panels are NaN. This field is BINARY and extremely "
                    "sparse, so same-valued pairs dominate every window: sen_slope is "
                    f"exactly 0 on {sen_zero:.1%} of finite cells in the final panel. READ "
                    "ols_slope. Note a slope of 0 over a permanently silent cell is a TRUE "
                    "zero here, not a censored one -- the opposite situation to "
                    "`heatwave-3b`, where zeros come from saturation at the ceiling."),
                "slope_units": "1 decade-1",
                "sparsity_caveat": (
                    "THE ABSOLUTE HUMIDEX GATE IS SILENT OVER MOST OF THE EXTRATROPICS, AND "
                    "A ZERO HERE MEANS 'NEVER CROSSES HUMIDEX 45', NOT 'LOW HEAT RISK'. "
                    "Measured over 2006-2099 across both RCPs and all 4 GCMs on the "
                    "67,420-cell land mask: 40,265 cells (59.7%) are NEVER non-zero in any "
                    "member or year. Share of land that never registers a single exposed "
                    "year, by belt: 60-90N 99%, 35-60N 81%, 35-60S 68%, 23.5-35N 30%, "
                    "23.5-35S 19%, tropics 18%. Paris, Frankfurt and Yakutsk are EXACTLY "
                    "ZERO in every member and every year. This is the absolute threshold "
                    "behaving as designed -- the model ran everywhere, the threshold is "
                    "genuinely not crossed -- and NOT the withdrawn-sugarcane failure mode, "
                    "where the model did not run. But a delivery that renders a low "
                    "percentile across Europe without saying which of the two it means is "
                    "reporting a threshold artifact as a finding. "
                    f"Per-decade share of land at exactly zero on this scenario: "
                    f"{silent_txt}."),
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal frequency ranked against the shared "
                    f"{BASELINE_DECADE}s ensemble land distribution. Not inverted -- more "
                    "heatwave exposure is worse. The tier-1 block is large BY CONSTRUCTION "
                    "on this layer and means 'below the health threshold this index "
                    "measures', not 'no heat'."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_tier1_share": round(tier1, 5),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_treatment": f"UNIFORM {','.join(socs)} -- no compromise needed",
                "co2_treatment": f"UNIFORM {','.join(senss)}",
                "normalization": (
                    "none -- every member reports the same dimensionless binary exposure "
                    "flag from the same index code; only the driving GCM differs."),
                "spatial_smoothing": (
                    "none -- decided from the per-member contact sheet, NOT inherited from "
                    "`let` (also 4 members, L=2.5). `let` needed a kernel because its raw "
                    "field is one-cell-wide storm tracks; this field's active regions are "
                    "broad contiguous zones (Sahel, Gulf, South Asia, maritime SE Asia) and "
                    "its sparsity is a THRESHOLD boundary, not sampling noise. Smoothing "
                    "would bleed non-zero frequency across that boundary into cells the "
                    "index says never cross Humidex 45, which is precisely the claim this "
                    "layer must not blur."
                    if not args.smooth_l else
                    f"5x5 exponential decay, L={args.smooth_l}, applied per member per YEAR "
                    "(so the CI and both slopes see the smoothed series)."),
                "land_mask_rule": (
                    "THE SOURCE PRODUCT HAS NO MASK -- it is finite on all 259,200 grid "
                    "cells and zero-fills ocean, so `isfinite` carries no footprint. The "
                    "ISIMIP2b land-sea mask (`LSM`, not the 3b `mask`) is applied "
                    f"EXPLICITLY here: {int(land.sum()):,} land cells published, everything "
                    "else NaN. Without this the entire ocean would publish as a confident "
                    "zero and would enter the percentile reference distribution."),
                "minimum_models": args.min_models,
                "mask_rule": (
                    "THE MINIMUM-MODEL RULE IS VOID, not relaxed -- one index model "
                    "(hwmid-humidex), so n_models is 1 at every cell and the `led` layer's "
                    "'>= 2 impact models' rule cannot be expressed. No GCM-count fallback "
                    "either: after the land mask is applied the four GCM members cover "
                    "identical cells."),
                "interpretation_caveat": (
                    "THIS IS AN INTERSECTION OF A RELATIVE AND AN ABSOLUTE CRITERION, AND "
                    "BOTH APPLY. A cell counts as exposed only when the relative HWMId "
                    "exceeds its threshold AND the absolute Humidex reaches 45 -- inside "
                    "Environment Canada's 40-45 'great discomfort; avoid exertion' band -- "
                    "so counted events are health-relevant by construction (Lange et al. "
                    "2020, doi:10.1029/2020EF001616). But the health filter sits ON TOP OF a "
                    "per-cell-distribution framing, it does not replace it: the relative "
                    "gate still has to open first. So this layer is neither a pure absolute "
                    "heat-stress index nor a pure departure-from-normal index, and it should "
                    "not be described as either. Its ISIMIP3b successor `HWMID-NONE` DROPS "
                    "the Humidex criterion entirely and is relative-only."),
                "ensemble_agreement": (
                    "THE ACTIVE AREA IS SUBSTANTIALLY SINGLE-GCM, AND THAT IS A SEPARATE "
                    "CAVEAT FROM THE SPARSITY. On a sparse threshold index the first "
                    "question is not whether members agree on a level but whether they "
                    "agree the threshold is ever crossed at all. Measured on the "
                    f"{BASELINE_DECADE}s baseline, cells active in exactly k of "
                    f"{n_mem_tot} GCMs: {agree_txt}. Of the {active_any:,} cells active in "
                    f"ANY member, {solo_share:.1%} rest on a SINGLE GCM. Per-member global "
                    "land mean frequency spans 0.0204 (gfdl-esm2m) to 0.0440 "
                    "(ipsl-cm5a-lr), a 2.16x spread with CV 28.5%. "
                    "FOR CONTRAST the ISIMIP3b `heatwave` layer, which drops the absolute "
                    "criterion, has 74.8% of land active in all 5 of its GCMs and only 1.2% "
                    "in two or fewer, at 1.58x spread and CV 17.5%. The absolute gate does "
                    "not merely shrink the footprint -- it makes the surviving footprint "
                    "much more GCM-dependent, because a hard threshold amplifies small "
                    "differences in the driving climate. Filter on n_members before "
                    "reporting a site-level value from a thinly-agreed cell."),
                "ensemble_note": (
                    "1 index model (hwmid-humidex) x 4 CMIP5 GCMs = 4 members per scenario, "
                    "COMPLETE (1 x 4 x 2 = 8 files enumerated 2026-08-13, no missing "
                    "combination). Only rcp26 and rcp60 exist for this product -- there is "
                    "NO high-forcing member for this hazard, even though the ISIMIP2b "
                    "climate forcing itself publishes rcp85. That matters here: the Humidex "
                    "gate is reported to open rarely in temperate regions below 2 degC of "
                    "warming, and rcp26/rcp60 are exactly that range."),
                "source_dataset": (
                    "ISIMIP2b DerivedOutputData/Lange2020/HWMId-humidex. SEPARATE LAYER from "
                    "the ISIMIP3b `heatwave` layer and NOT superseded by it: different GCM "
                    "generation, different scenarios, and a DIFFERENT INDEX. All 8 files "
                    "sha512-verified against the publisher sidecars at ingest 2026-08-14."),
                "description": (
                    "ISIMIP2b/RCP heatwave exposure frequency processed to the TCFD output "
                    "contract (OUTPUT-SPEC.md) with a shared 2020s baseline; "
                    f"{len(mem)}-member 1-index-model x 4-GCM CMIP5 ensemble, "
                    f"{stat_name} decadal statistic, no smoothing, explicit ISIMIP2b land "
                    f"mask, {pct_mode} percentile, higher_is_worse. Health-anchored and "
                    "SPARSE -- read sparsity_caveat before using it."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)")
        log(f"    land at exactly zero by decade: {silent_txt}")
        log(f"    final panel: sen_slope exactly 0 on {sen_zero:.1%} of finite cells")
        _CUBE = None

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
