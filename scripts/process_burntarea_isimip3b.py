"""Process ISIMIP3b burntarea-total (wildfire burnt-area) to the TCFD output contract.

burntarea-total = the percent of each grid cell burned, reported by ISIMIP3b fire /
vegetation models in PERCENT. This is the direct biophysical fire signal, NOT the Lange2020
`lew` / Zantout2025 `wildfire` *exposure* framing.

Ensemble: ALL ISIMIP3b burntarea-total, 5 impact models x their CMIP6 GCMs =
**22 members per scenario** x {ssp126, ssp370, ssp585} (user decision 2026-08-08: maximize
members across BOTH GCMs and impact models, and prefer the SSP round over RCP).

    model             GCMs  soc                   cadence   ann_mean %  zero%   land cells
    classic             2   2015soc               monthly    3.57-3.72  44-45%  ~66,660
    elm-eca             5   2015soc-from-histsoc  monthly    2.85-3.45  37-39%   60,723
    lpjml5-7-10-fire    5   2015soc               monthly    1.89-2.03  20%      67,420
    visit               5   2015soc               monthly    3.18-3.45  39-41%   58,714
    mc2-usfs-r87g5c1    5   nat                   ANNUAL     3.26-3.62  43-44%   58,919

MONTHLY -> ANNUAL = **SUM** (user instruction 2026-08-08: "annual totals"). Burnt area
ACCUMULATES; the mean would under-scale fire 12x. `mc2-usfs` already publishes annual and is
read directly, then pooled with the summed members.

The sum is NaN-preserving: `np.nansum` maps an all-NaN ocean cell to 0.0, which then reads
as "land that never burns" and silently triples the apparent land mask (measured: 259,200
"land" cells instead of ~67,000). A year contributes only if all 12 of its months are
finite.

Value-check, 2026-08-08 (GUARDRAILS §9 -- measured, not inherited):

* All 5 models declare `units="%"`, `long_name="Burnt Area Fraction"`, a `days since`
  time axis on a `365_day` calendar, and `_FillValue=1e20` for ocean. Unusually uniform for
  this variable -- the ISIMIP2b generation mislabelled lpj-guess's long_name and mixed
  `days since` with `years since`.
* Annual-total land means span 1.89-3.72 % (a ~2x spread) -- the SAME unit at COMPARABLE
  magnitude, so **no normalization**: equal-weight model democracy in raw %, and the
  inter-member spread becomes the CI.
* Within-model cross-GCM spread is 1.0-1.2x, i.e. **no mis-scaling**. This is why `classic`
  is taken from its `2015soc` run: its `2015soc-from-histsoc` run is documented mis-scaled
  ACROSS GCMs within the one model (gfdl fraction-scaled, ukesm percent, identical
  metadata). Verified clean here at 1.0x.
* Strongly **zero-inflated** over land (20-45% exact zeros) -> two-tier percentile, and
  `ols_slope` is the slope to read (`sen_slope` collapses to 0 where most year-pairs are
  0 -> 0). Not boolean: the values are continuous, so the median/IQR branch applies.
* Annual totals EXCEED 100% for `classic` (max ~151) and `elm-eca` (max ~575). A cell that
  reburns within a year legitimately accumulates >100% burnt area, so values are NOT
  clipped at 100. `elm-eca` reaching exactly 100.00 in single months repeatedly is flagged
  in `outlier_note` rather than silently corrected.

Direction: higher = worse (more fire), so the percentile is NOT inverted.

Output: `burntarea_{scenario}_processed.nc`, dims (decade, lat, lon), decades 2020s-2090s
(ISIMIP3b starts 2015, so there is no full 2010s decade), carrying
{median, lower_ci, upper_ci, percentile, ols_slope, sen_slope, n_members, n_models}
per OUTPUT-SPEC.md.

Usage:
    python scripts/process_burntarea_isimip3b.py                # full run
    python scripts/process_burntarea_isimip3b.py --limit-cells 200 --scenarios ssp370
    python scripts/process_burntarea_isimip3b.py --max-pairs 200000   # fast, approximate
"""

import argparse
import glob
import os
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.utils.decadal_stats import (  # noqa: E402
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

VAR = "burntarea-total"
OUT_VAR = "burntarea"
LAYER_ID = "wildfire-isimip3b_burntarea-total_annual"
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
#: ISIMIP3b runs start 2015, so 2015-2019 can never form a full 2010s decade and the
#: expanding window starts at the 2020s baseline -- those years are simply unused.
MIN_YEAR, MAX_YEAR = 2020, 2099
CI_FLOOR = 0.0          # burnt area is non-negative
CI_CAP = None           # NO upper cap: cumulative annual burnt area may exceed 100%
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False  # more fire = worse -> percentile NOT inverted

SLOPE_PER_DECADE = 10.0

#: All five models are structurally distinct codes, so each is its own family. Explicit so
#: a future orchidee/orchidee-dgvm-style pair cannot silently get two votes.
MODEL_FAMILY = {}

#: Peak slope memory is ~4 x n_pairs x chunk_cells x itemsize, and this ensemble has ~1.5M
#: pairs per cell at the widest panel -- the module default of 512 would want ~12.7 GB.
#: Size the chunk from the actual pair count instead.
SLOPE_MEM_BUDGET_BYTES = 700 * 1024**2


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, cadence, member) from an ISIMIP3b filename.

    e.g. visit_gfdl-esm4_w5e5_ssp370_2015soc_default_burntarea-total_global_monthly_2015_2100.nc

    Parsed from the END for the variable/cadence fields, because the variable token itself
    contains hyphens; model/gcm/scenario are safe from the front (no underscores in names).
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], soc=p[4],
                cadence=p[-3], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Integer calendar years from a CF time axis (xarray cannot decode 365_day)."""
    t = ds.time
    units = t.attrs.get("units", "")
    cal = t.attrs.get("calendar", "365_day")
    m = re.search(r"(years|days|months)\s+since\s+(\d+)", units)
    if not m:
        raise ValueError(f"cannot parse time units {units!r}")
    step, base = m.group(1), int(m.group(2))
    vals = t.values.astype("float64")
    if step == "years":
        yrs = base + vals
    elif step == "months":
        yrs = base + vals / 12.0
    else:
        dpy = 360.0 if "360" in cal else 365.0 if "365" in cal else 365.25
        yrs = base + vals / dpy
    return np.floor(yrs + 1e-6).astype(int)


def annualize(fpath, cache_dir):
    """Load one member as annual totals -> (years, (year, lat, lon) float32).

    Monthly members are SUMMED over each calendar year (burnt area accumulates). A year is
    emitted only where all 12 months are finite, so an all-NaN ocean cell stays NaN instead
    of becoming a finite 0.0.

    Reading + summing 6 GB of monthly data is the slow part of this script, so the result is
    cached as .npz; delete the cache dir to force a re-read.
    """
    base = os.path.basename(fpath).replace(".nc", "")
    cache = cache_dir / f"{base}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["years"], z["cube"]

    info = parse_name(fpath)
    ds = xr.open_dataset(fpath, decode_times=False, decode_cf=False)
    da = ds[VAR]
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    yrs = years_from_time(ds)
    vals = da.values.astype("float32")
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan

    uy = np.array(sorted(set(int(y) for y in yrs)))
    keep = (uy >= MIN_YEAR) & (uy <= MAX_YEAR)
    uy = uy[keep]
    out = np.full((uy.size, vals.shape[1], vals.shape[2]), np.nan, np.float32)

    for k, y in enumerate(uy):
        sel = np.where(yrs == y)[0]
        blk = vals[sel]
        if info["cadence"] == "monthly":
            if blk.shape[0] != 12:
                # Partial year: leave NaN rather than emit an under-counted total.
                continue
            full = np.isfinite(blk).all(axis=0)
            tot = np.nansum(blk, axis=0)
            out[k] = np.where(full, tot, np.nan)
        else:
            out[k] = blk[0] if blk.shape[0] == 1 else np.nanmean(blk, axis=0)

    del vals
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, years=uy, cube=out)
    return uy, out


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline land distribution.

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros), which it is
    for burnt area: a cell that never burns is the LOWEST fire risk, so zeros -> 1 and
    positives are ranked against the NON-ZERO baseline into [2, 100]. Ranking zeros inside
    the full distribution would instead hand every unburnt cell the ~20-45th percentile.
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
            res = np.ones(vals.shape, np.float32)   # never-burns -> raw 1
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
    """Put a land-cell vector back on the (lat, lon) grid; ocean stays NaN."""
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


def slope_chunk_cells(n_members, n_years, max_pairs):
    """Chunk width that keeps Theil-Sen peak memory inside the budget."""
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    if max_pairs is not None:
        pairs = min(pairs, max_pairs)
    per_cell = 4 * pairs * 4          # 4 float32 operands of n_pairs each
    return max(8, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", nargs="*", default=None,
                    help="subset of scenarios (default: all found)")
    ap.add_argument("--limit-cells", type=int, default=None,
                    help="process only the first N land cells (benchmarking)")
    ap.add_argument("--max-pairs", type=int, default=None,
                    help="cap Theil-Sen pairs; default exact. Iteration only.")
    ap.add_argument("--members-only", action="store_true",
                    help="write only the per-member diagnostic and exit (fast, uses cache)")
    ap.add_argument("--skip-slopes", action="store_true",
                    help="write NaN slopes (fast structural check)")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    cache_dir = raw_dir / "_annual_cache"
    out_dir = root / "data" / "processed" / LAYER_ID
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_*_2015_2100.nc")))
    if not files:
        log(f"ERROR: no {VAR} member files in {raw_dir}")
        return 1
    meta = {f: parse_name(f) for f in files}
    # The shared 2020s baseline is pooled over EVERY scenario, so all scenarios are always
    # loaded even when only a subset is written. Selecting with --scenarios must not change
    # the baseline -- otherwise parallel per-scenario runs would each emit a different
    # "shared" baseline and the panels would not be bit-identical across files.
    scenarios = sorted({m["scenario"] for m in meta.values()})
    write_scenarios = ([s for s in scenarios if s in args.scenarios]
                       if args.scenarios else list(scenarios))
    if not write_scenarios:
        log(f"ERROR: --scenarios {args.scenarios} matched none of {scenarios}")
        return 1

    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    soc_by_model = {m["model"]: m["soc"] for m in meta.values()}
    cad_by_model = {m["model"]: m["cadence"] for m in meta.values()}

    log("=" * 74)
    log("Processing ISIMIP3b burntarea-total (wildfire) -> TCFD output contract")
    log("=" * 74)
    log(f"Files: {len(files)} | scenarios: {scenarios} | writing: {write_scenarios}")
    log(f"Impact models ({len(models)}): {models}")
    log(f"GCMs ({len(gcms)}): {gcms}")
    log(f"soc by model: {soc_by_model}")
    log(f"cadence by model: {cad_by_model}")
    log("Monthly -> annual = SUM (burnt area accumulates); mc2-usfs is already annual.")
    log("No normalization (raw %, model democracy); no spatial smoothing; higher = worse.")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    # ---- Pass A: land mask only, one member in memory at a time -------------
    # Holding all 66 members as full 360x720 grids costs 83 MB each = 5.5 GB, which on a
    # 16 GB machine drives the process into swap (measured: 19.5M swapouts, a 2070s panel
    # that should take ~11 min took 9.2 HOURS). Scan for the mask first, then pack straight
    # into the land-cell array and never materialise a full-grid stack.
    log(f"\nPass A: land mask (cache: {cache_dir})...")
    t0 = time.time()
    LAT = LON = None
    finite_any = None
    members_by_scen = {s: [] for s in scenarios}
    for f in files:
        info = meta[f]
        uy, cube = annualize(f, cache_dir)
        if LAT is None:
            LAT, LON = cube.shape[1], cube.shape[2]
            finite_any = np.zeros((LAT, LON), bool)
        finite_any |= np.isfinite(cube).any(axis=0)
        members_by_scen[info["scenario"]].append(info["member"])
        del cube
    for s in scenarios:
        members_by_scen[s] = sorted(members_by_scen[s])
    land_idx = np.flatnonzero(finite_any.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        # Evenly spaced across the land mask, not the first N -- a head slice is all
        # polar cells and benchmarks/QA on it are unrepresentative.
        land_idx = land_idx[np.linspace(0, land_idx.size - 1,
                                        args.limit_cells).astype(int)]
    n_land = land_idx.size
    log(f"  land cells (union over members): {n_land:,} of {LAT * LON:,} "
        f"[{time.time() - t0:.0f}s]")

    with xr.open_dataset(files[0], decode_times=False, decode_cf=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values

    # ---- Pass B: pack (member, year, land_cell) directly --------------------
    log("Pass B: packing land-cell cubes...")
    t0 = time.time()
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        uy, cube = annualize(f, cache_dir)
        flat = cube.reshape(cube.shape[0], -1)
        for k, y in enumerate(uy):
            yi = y_index.get(int(y))
            if yi is not None:
                annual[s][slot[s][m], yi] = flat[k, land_idx]
        del cube, flat
    log(f"  packed {len(files)} members in {time.time() - t0:.0f}s "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident)")

    # ---- Field nature: measured, never assumed (GUARDRAILS §9) --------------
    boolean = is_boolean_field(annual[scenarios[0]])
    stat_name = "pooled_mean_boolean" if boolean else "pooled_median"
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> decadal statistic = {stat_name}")
    if boolean:
        log("  WARNING: burnt-area % is continuous; a boolean classification means the "
            "input is not what this processor expects. Check the members.")

    for s in scenarios:
        cov = np.isfinite(annual[s]).any(axis=1)
        log(f"  {s}: {len(members_by_scen[s])} members, "
            f"per-member land cells {cov.sum(axis=1).min():,}-{cov.sum(axis=1).max():,}")

    # ---- Shared 2020s baseline ---------------------------------------------
    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios; shared baseline is "
            "only valid for a uniform ensemble. Declaring members_by_scenario.")
    # Slice to the 10-year baseline window BEFORE concatenating: pooling the full 80-year
    # cubes would copy 1.44 GB when only 180 MB is ever read.
    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    base_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years[bwin], BASELINE_DECADE, boolean=boolean,
        window_years=WINDOW_YEARS)
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, CI_CAP)
    b_hi = np.clip(b_hi, CI_FLOOR, CI_CAP)

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"\nShared 2020s baseline: land n={baseline_flat.size:,}, "
        f"exact-zero fraction={frac_zero:.2%}, percentile mode={pct_mode}, "
        f"global-mean burnt={np.nanmean(b_med):.4f} %, "
        f"max={np.nanmax(b_med):.2f} %")

    # ---- Per-member diagnostic for the dashboard "Members" tab --------------
    # (member, lat, lon) mean annual value over the shared baseline decade, pooled across
    # scenarios. Plotting every LSM x GCM on ONE scale is how a member that runs
    # systematically hot/cold -- or whose spatial distribution is unlike its siblings --
    # becomes visible; the flat member pool gives each one an equal vote, so an outlier
    # matters. Cheap to build from the annualised cache, so it is emitted on every run.
    mem_names = members_by_scen[scenarios[0]]
    mem_grid = np.full((len(mem_names), LAT, LON), np.nan, np.float32)
    for mi, mname in enumerate(mem_names):
        stack = []
        for s in scenarios:
            if mname in members_by_scen[s]:
                stack.append(annual[s][members_by_scen[s].index(mname)][bwin])
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
            "variable": OUT_VAR,
            "units": "%",
            "member_field": (f"mean annual burnt area over the {BASELINE_DECADE}s "
                             "baseline decade, pooled across all scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Consumed "
                     "by scripts/generate_maps.py to render the Members tab."),
        },
    )
    mem_path = out_dir / f"{OUT_VAR}_members.nc"
    mem_ds.to_netcdf(mem_path, encoding={"value": {"dtype": "float32", "zlib": True,
                                                   "complevel": 4,
                                                   "_FillValue": np.float32(np.nan)}})
    log(f"  wrote per-member diagnostic {mem_path.name} ({len(mem_names)} members)")
    del mem_grid
    if args.members_only:
        log("\n--members-only: done.")
        return 0

    # The baseline is now fixed, so scenarios we are not writing can be released.
    for s_drop in [s for s in scenarios if s not in write_scenarios]:
        del annual[s_drop]

    chunk = slope_chunk_cells(max(len(m) for m in members_by_scen.values()),
                              n_years, args.max_pairs)
    log(f"Theil-Sen chunk_cells={chunk} (memory-bounded to "
        f"~{SLOPE_MEM_BUDGET_BYTES / 1024**2:.0f} MB/panel)")

    # ---- Per-scenario assembly --------------------------------------------
    for s in write_scenarios:
        log(f"\n{'=' * 74}\nAssembling scenario {s}\n{'=' * 74}")
        mem = members_by_scen[s]
        cube = annual[s]
        fams = sorted({family_of(m.split("_")[0]) for m in mem})
        fam_idx = {f: k for k, f in enumerate(fams)}

        nd = len(DECADES)
        out = {k: np.full((nd, LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}

        for i, d in enumerate(DECADES):
            td = time.time()
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, boolean=boolean, window_years=WINDOW_YEARS)
                lo = np.clip(lo, CI_FLOOR, CI_CAP)
                hi = np.clip(hi, CI_FLOOR, CI_CAP)
                pc = pct(med)

            if args.skip_slopes:
                sl = dict(ols_slope=np.full(n_land, np.nan, np.float32),
                          sen_slope=np.full(n_land, np.nan, np.float32))
            else:
                sl = expanding_slopes(cube, years, d, BASELINE_DECADE,
                                      window_years=WINDOW_YEARS,
                                      chunk_cells=chunk, max_pairs=args.max_pairs)

            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
            fam_present = np.zeros((len(fams), n_land), bool)
            for mi, m in enumerate(mem):
                fam_present[fam_idx[family_of(m.split("_")[0])]] |= present[mi]
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
                slope_txt = (f"ols={np.nanmean(out['ols_slope'][i]):+.4f} "
                             f"sen={np.nanmean(out['sen_slope'][i]):+.4f} %/dec"
                             if d != BASELINE_DECADE else "slopes=NaN (baseline)")
            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
            log(f"  {d}s: {tag:<15} mean={np.nanmean(out['median'][i]):.4f}%  "
                f"{slope_txt}  [{time.time() - td:.0f}s]")

        # ---- GUARDRAIL: slope and median masks must agree -------------------
        # A bare np.zeros() baseline panel would make the whole OCEAN a finite zero and QA
        # would not catch it (it only checks that FINITE baseline slopes are zero).
        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), "baseline ols must be NaN"
                assert np.all(np.isnan(out["sen_slope"][i])), "baseline sen must be NaN"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), (
                    f"{k} finite where median is NaN at {d}s ({extra.sum()} cells)"
                    " -- ocean leak")
            assert np.all(out["lower_ci"][i][med_finite] <= out["median"][i][med_finite]
                          + 1e-5), f"lower_ci > median at {d}s"
            assert np.all(out["upper_ci"][i][med_finite] >= out["median"][i][med_finite]
                          - 1e-5), f"upper_ci < median at {d}s"

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Annual burnt area (percent of grid cell burned per year)",
                "units": "%",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "field_nature": "boolean_01" if boolean else "continuous",
                "value_note": (
                    "median = MEDIAN over the pooled (year x member) sample inside the "
                    f"decade window, across {len(mem)} ISIMIP3b members "
                    "(classic, elm-eca, lpjml5-7-10-fire, visit, mc2-usfs x their CMIP6 "
                    "GCMs) in raw percent burnt area."),
                "annual_aggregation": (
                    "Monthly members (classic, elm-eca, lpjml5-7-10-fire, visit) were "
                    "SUMMED over each calendar year -- burnt area ACCUMULATES, so a mean "
                    "would under-scale fire 12x. A year contributes only where all 12 "
                    "months are finite, so ocean stays NaN instead of summing to 0. "
                    "mc2-usfs publishes annual values and is read directly."),
                "ci_definition": (
                    "lower_ci/upper_ci = 25th/75th percentile (IQR) of the same pooled "
                    "(year x member) sample, floored at 0 and NOT capped at 100. The IQR "
                    "carries BOTH interannual variability and inter-model spread; it is "
                    "not a pure model-spread band."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "fitted over an EXPANDING window from the start of the 2020s baseline "
                    "through the end of the target decade, stacking every (year, member) "
                    "observation as an independent point. The baseline panel is NaN (no "
                    "elapsed period). The estimators fail in OPPOSITE regimes -- sen "
                    "collapses to exactly 0 on zero-inflated fields, ols absorbs member "
                    "level offsets as trend when coverage is uneven -- so disagreement "
                    "between them means a cell's trend is not robust. This field is "
                    f"{frac_zero:.1%} exact zeros in the baseline, i.e. ZERO-INFLATED, so "
                    "ols_slope is the slope to read here."),
                "slope_units": "% decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal median ranked against the shared "
                    "2020s ensemble land distribution. Zeros (cells that never burn) map "
                    "to 1; positives rank against the NON-ZERO baseline into [2,100]. Not "
                    "inverted -- more fire is worse."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "decade_note": (
                    "ISIMIP3b runs start 2015, so there is no full 2010s decade; the "
                    "layer begins at the 2020s baseline (2020s-2090s)."),
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_by_model": ";".join(f"{k}:{v}" for k, v in sorted(soc_by_model.items())),
                "soc_treatment": (
                    "MIXED, deliberately (user decision 2026-08-08: maximize ensemble "
                    "depth). No single soc token spans the ensemble -- elm-eca publishes "
                    "ONLY 2015soc-from-histsoc, visit ONLY 2015soc, mc2-usfs ONLY nat -- "
                    "so a uniform filter would drop 5-15 of 22 members. classic is taken "
                    "from 2015soc specifically to AVOID its 2015soc-from-histsoc run, "
                    "which is documented mis-scaled across GCMs within the one model; "
                    "measured cross-GCM ratio here is 1.0x (clean)."),
                "co2_treatment": (
                    "UNIFORM transient CO2: every member is the 'default' sens run."),
                "normalization": (
                    "none -- all 5 models report the SAME unit (% burnt area) with "
                    "COMPARABLE annual-total land means (1.89-3.72 %), so they are "
                    "equal-weighted in raw % (model democracy). Note classic contributes "
                    "only 2 GCMs vs 5 each for the others, so it is underweighted in the "
                    "flat member pool."),
                "outlier_note": (
                    "Cumulative annual burnt area may legitimately exceed 100% where a "
                    "cell reburns within the year, so values are NOT clipped: measured "
                    "annual maxima are classic ~151%, elm-eca ~575%, lpjml/visit ~100%, "
                    "mc2-usfs 100%. elm-eca reports exactly 100.00% in single months "
                    "repeatedly, which is flagged rather than corrected. Per-member land "
                    "masks differ (58,714-67,420 cells), so n_members/n_models are "
                    "emitted per cell and a level step between multi- and single-model "
                    "cells is expected at mask edges."),
                "spatial_smoothing": f"none ({len(mem)}-member ensemble is thick)",
                "source_dataset": (
                    "ISIMIP3b OutputData/fire (classic, elm-eca, lpjml5-7-10-fire, visit) "
                    "+ OutputData/biomes (mc2-usfs, the only annual publisher). The fire "
                    "and biomes sectors publish byte-identical burntarea files, so each "
                    "model is ingested from exactly one sector."),
                "description": (
                    "Wildfire burnt-area percent processed to the TCFD output contract "
                    f"(OUTPUT-SPEC.md) with a shared 2020s baseline; {len(models)}-model "
                    f"x CMIP6-GCM ensemble ({len(mem)} members) in raw %, monthly summed "
                    "to annual totals, no normalization, no spatial smoothing, "
                    "higher_is_worse."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)")

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
