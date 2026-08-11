"""Process `let` (tropical-cyclone exposure) into the TCFD output contract.

`let` (Lange 2020, ISIMIP2b DerivedOutputData, impact model ke-tg-meanfield) is the
tropical-cyclone member of the extreme-event exposure family: each cell holds the
FRACTION OF ITS LAND AREA exposed to tropical-cyclone winds that year.

Rebuilt 2026-08-11 onto [OUTPUT-SPEC.md](../OUTPUT-SPEC.md). The 2026-07-24 version was a
"family B" processor -- it emitted the retired single `trend`, used mean +/- 1 SD as the
CI, and shipped no `ols_slope` / `sen_slope` / `n_members` / `n_models` and no
`{var}_members.nc`. It is fully superseded by this file.

MEASURED data nature (GUARDRAILS 9 -- from the values, never the name):

  * continuous in [0, 0.952]; NEVER exactly 1; `is_boolean_field` -> False
  * 97.84% of finite cell-years are EXACTLY ZERO
  * land mask identical across all 8 files (68,249 cells), and time-invariant
  * 4 members/scenario (1 impact model x 4 CMIP5 GCMs), uniform across rcp26/rcp60,
    so the shared 2020s baseline is valid
  * soc/sens uniformly nosoc/co2 -- no harmonization compromise

TWO DECISIONS SPECIFIC TO THIS LAYER (user, 2026-08-11). Neither transfers elsewhere.

1. THE DECADAL STATISTIC IS THE POOLED **MEAN**, not the contract's default median.

   OUTPUT-SPEC sends a zero-inflated *continuous* field down the median/IQR branch. At
   97.84% annual zero-inflation that branch does not merely lose precision, it erases the
   layer. Measured on the 2020s panel, after smoothing:

       branch                     exposed cells    exact-zero land
       pooled median + IQR             2,684            96.07%
       pooled mean +/- 1 SD           15,122            77.84%

   93% of exposed land reads exactly 0 under the median, and the two-tier percentile then
   assigns tier 1 to 96% of land. `let` is the continuous analogue of `led`: both estimate
   an EXPECTED ANNUAL EXPOSED FRACTION, and the spec already grants `led` the mean because
   a median of Bernoulli draws is meaningless. So this layer passes `central="mean"` and
   DECLARES it as `decadal_statistic: pooled_mean_zero_inflated`. The CI is correspondingly
   mean -/+ 1 SD of the same pooled (year x member) sample.

   For contrast: `burntarea` is 29.2% zeros and the median branch handles it fine. This is
   a regime, not a slippery slope -- measure the median's exact-zero share before reaching
   for this branch on any other layer.

2. AGGRESSIVE SPATIAL SMOOTHING, applied per member per YEAR on the raw data.

   Raw `let` is sparse: individual storm tracks, one cell wide, rather than a smooth
   projected impact field. The 4-member ensemble is far too thin to average that out. The
   previous kernel (5x5, L=0.7, described in its own code as "sharp") kept 32.1% of the
   weight on the centre cell and removed only 63% of the track structure.

   This layer uses the SAME 5x5 / 111 km footprint with a much flatter profile, L=2.5,
   which puts 8.1% on the centre. Measured on the 2020s panel, roughness = mean
   |cell - its 4-neighbour mean| over exposed land, normalized by the mean:

       kernel                        centre wt   roughness   exposed cells
       none (raw)                      100.0%      0.389        10,794
       5x5 L=0.7 (previous, "sharp")    32.1%      0.142        15,122
       5x5 L=2.5 (THIS LAYER)            8.1%      0.044        15,122
       7x7 L=2.0                         6.8%      0.035        16,664

   L=2.5 was chosen over widening the window because it costs NO extra spatial reach: it
   adds no newly-exposed cell relative to the previous kernel, whereas every radius
   increase asserts that land further from any track is exposed (7x7 adds 1,542 cells,
   9x9 adds 3,940). The radius is physically anchored -- 2 cells = 111 km at the equator,
   about the radius of hurricane-force winds. Do not widen past 3 cells (167 km, the outer
   damaging-wind envelope); beyond that the kernel is smearing, not borrowing strength.

   Smoothing is normalized over non-NaN LAND neighbours only, longitude-wrapped at the
   antimeridian, latitude-clamped at the grid edge, and cos(lat)-scaled in longitude so the
   footprint stays ~isotropic in km. Ocean is never filled. Mass is conserved to +3.1%
   (global mean 0.004322 -> 0.004456); the residual is coastal redistribution.

   NOTE it is applied to each member's ANNUAL map before any pooling, as specified. Because
   smoothing is linear and the land mask is time-invariant, this yields a `median` panel
   identical to smoothing the decadal mean -- but the CI and `sen_slope` genuinely differ,
   which is why it is done here and not later.

Percentile: two-tier, as the baseline is heavily zero-inflated -- zeros -> 1, positives
ranked against the NON-ZERO 2020s baseline into [2,100]. `higher_is_worse` (more exposure
= more risk), so it is NOT inverted.

Expect `sen_slope` to collapse toward exactly 0 over much of the domain -- that is the
documented zero-inflation failure mode, and it is why both slopes are emitted. READ
`ols_slope` ON THIS LAYER.

Usage:
    python scripts/process_let_cyclone.py [--scenarios rcp26] [--skip-slopes]
                                          [--limit-cells N] [--members-only]
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.decadal_stats import (  # noqa: E402
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

VAR = "let"
OUT_VAR = "let"
LAYER_ID = "cyclone_let_annual"
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
#: Source files span 2006-2099; 2006-2009 cannot form a full 2010s decade and are unused.
MIN_YEAR, MAX_YEAR = 2010, 2099

#: `let` is a fraction of a cell's land area, so it is genuinely bounded by [0, 1] --
#: unlike burnt area, which may exceed 100% when a cell reburns. Measured max is 0.952.
CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False          # more exposure = worse -> percentile NOT inverted
SLOPE_PER_DECADE = 10.0

#: Decadal central value: "mean" (this layer) or "median" (the contract default).
#: See decision 1 in the module docstring. Changing this WITHOUT re-measuring the
#: median's exact-zero share is a contract violation, not a tuning knob.
CENTRAL = "mean"
STAT_NAME = "pooled_mean_zero_inflated"

#: Smoothing kernel -- see decision 2. L is the exponential length scale in grid cells.
SMOOTH_L = 2.5
SMOOTH_RADIUS = 2                 # 5x5 window = 111 km at the equator

#: Single impact model, so every member shares one family and n_models is 1 everywhere.
MODEL_FAMILY = {}

SLOPE_MEM_BUDGET_BYTES = 700 * 1024**2


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) from a Lange 2020 filename.

    e.g. lange2020_ke-tg-meanfield_gfdl-esm2m_ewembi_rcp26_nosoc_co2_let_global_annual_2006_2099.nc4

    DerivedOutputData filenames carry a LEADING PUBLICATION TOKEN (`lange2020_`) that
    shifts every forward index by one relative to OutputData names. The variable and
    cadence are therefore read from the END, which is index-stable either way.
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[1], gcm=p[2], scenario=p[4], soc=p[5], sens=p[6],
                variable=p[-5], cadence=p[-3], member=f"{p[1]}_{p[2]}")


def years_from_time(ds):
    """Integer calendar years from the CF axis ('years since 1661-1-1', 360_day)."""
    t = ds.time
    units = t.attrs.get("units", "")
    cal = t.attrs.get("calendar", "360_day")
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


# --------------------------------------------------------------------------- #
# Spatial smoothing
# --------------------------------------------------------------------------- #

def make_kernel(lats, L=SMOOTH_L, radius=SMOOTH_RADIUS):
    """Per-row weight vectors for every (di, dj) offset in the window.

    Distance is cos(lat)-scaled in longitude -- d = sqrt(di^2 + (dj*cos(lat))^2) -- so a
    longitude step counts for less as cells narrow toward the poles and the footprint
    stays roughly isotropic in km. Returns [(di, dj, w[nlat]), ...].
    """
    coslat = np.cos(np.deg2rad(lats)).astype(np.float64)
    offs = []
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            d = np.sqrt(di * di + (dj * coslat) ** 2)
            offs.append((di, dj, np.exp(-d / L)))
    return offs


def _shift(a, di, dj):
    """a[..., i+di, j+dj] with longitude wrap; rows off the lat edge are zeroed."""
    nlat, nlon = a.shape[-2:]
    i_src = np.arange(nlat) + di
    lat_ok = (i_src >= 0) & (i_src < nlat)
    i_src = np.clip(i_src, 0, nlat - 1)
    j_src = (np.arange(nlon) + dj) % nlon
    out = a[..., i_src, :][..., :, j_src]
    if not lat_ok.all():
        out = out.copy()
        out[..., ~lat_ok, :] = 0.0
    return out


def smooth_cube(cube, kernel, valid2d):
    """Smooth every annual map in a (year, lat, lon) cube. Ocean stays NaN.

    `valid2d` is the member's time-invariant land mask, so the weight normalizer is
    computed ONCE rather than per year -- half the work, and it also guarantees every
    year is normalized identically (a per-year normalizer would let a member's own NaNs
    change the effective kernel from year to year).
    """
    f0 = np.where(np.isfinite(cube), cube, 0.0).astype(np.float32)
    ws = np.zeros_like(f0)
    wt = np.zeros(valid2d.shape, np.float64)
    vf = valid2d.astype(np.float64)
    for di, dj, w in kernel:
        wc = w[:, None]
        ws += (wc * _shift(f0, di, dj)).astype(np.float32)
        wt += wc * _shift(vf, di, dj)
    out = np.full(cube.shape, np.nan, np.float32)
    good = wt > 0
    out[:, good] = ws[:, good] / wt[good].astype(np.float32)
    out[:, ~valid2d] = np.nan          # never fill ocean
    return out


def load_member(fpath, kernel=None):
    """Load one member as (years, (year, lat, lon)), smoothed if a kernel is given."""
    ds = xr.open_dataset(fpath, decode_times=False, decode_cf=False)
    da = ds[VAR]
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    yrs = years_from_time(ds)
    vals = da.values.astype("float32")
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan

    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    vals, yrs = vals[keep], yrs[keep]

    valid = np.isfinite(vals)
    mask2d = valid.any(axis=0)
    # A time-varying land mask would invalidate the single-normalizer smoothing above.
    if not np.array_equal(valid, np.broadcast_to(mask2d, valid.shape)):
        raise ValueError(f"{os.path.basename(fpath)}: land mask varies by year")

    if kernel is not None:
        vals = smooth_cube(vals, kernel, mask2d)
    return yrs, vals, mask2d


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline land distribution.

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros), which it
    emphatically is here: a cell no cyclone reaches is the LOWEST exposure risk, so zeros
    -> 1 and positives are ranked against the NON-ZERO baseline into [2,100]. Ranking the
    zeros inside the full distribution would hand every unexposed inland cell a middling
    percentile and crush the real coastal gradient.
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
            res = np.ones(vals.shape, np.float32)      # never exposed -> raw 1
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
    per_cell = 4 * pairs * 4
    return max(8, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", nargs="*", default=None,
                    help="subset of scenarios to WRITE (the baseline always pools all)")
    ap.add_argument("--limit-cells", type=int, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--members-only", action="store_true")
    ap.add_argument("--skip-slopes", action="store_true")
    ap.add_argument("--no-smoothing", action="store_true",
                    help="diagnostic only -- disables the kernel this layer depends on")
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
    log(f"let (tropical-cyclone exposure) -> TCFD output contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")

    with xr.open_dataset(files[0], decode_times=False, decode_cf=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    LAT, LON = len(lats), len(lons)

    kernel = None if args.no_smoothing else make_kernel(lats)
    if kernel is None:
        log("WARNING: --no-smoothing: the sparse track structure will NOT be resolved.")
    else:
        r = int(np.argmin(np.abs(lats)))
        tot = sum(float(w[r]) for _, _, w in kernel)
        ctr = sum(float(w[r]) for di, dj, w in kernel if di == 0 and dj == 0)
        log(f"Smoothing: {2*SMOOTH_RADIUS+1}x{2*SMOOTH_RADIUS+1} w=exp(-d/{SMOOTH_L}), "
            f"cos(lat)-scaled, centre weight {100*ctr/tot:.1f}% at the equator")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    # ---- Load + smooth every member ---------------------------------------- #
    t0 = time.time()
    members_by_scen = {s: [] for s in scenarios}
    for f in files:
        members_by_scen[meta[f]["scenario"]].append(meta[f]["member"])
    for s in scenarios:
        members_by_scen[s] = sorted(members_by_scen[s])

    raw_gm, sm_gm = [], []
    cubes, finite_any = {}, np.zeros((LAT, LON), bool)
    for f in files:
        info = meta[f]
        yrs, cube, mask2d = load_member(f, kernel)
        if kernel is not None:
            _, raw_cube, _ = load_member(f, None)
            raw_gm.append(np.nanmean(raw_cube))
            sm_gm.append(np.nanmean(cube))
            del raw_cube
        finite_any |= mask2d
        cubes[f] = (yrs, cube)
        log(f"  loaded {info['gcm']:<13} {info['scenario']}  "
            f"{mask2d.sum():,} land cells  [{time.time() - t0:.0f}s]")

    if raw_gm:
        log(f"\nSmoothing mass check: raw global-mean {np.mean(raw_gm):.6f} -> "
            f"smoothed {np.mean(sm_gm):.6f} "
            f"({100*(np.mean(sm_gm)/np.mean(raw_gm) - 1):+.1f}%; "
            "smoothing redistributes, it must not inflate)")

    land_idx = np.flatnonzero(finite_any.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        land_idx = land_idx[np.linspace(0, land_idx.size - 1,
                                        args.limit_cells).astype(int)]
    n_land = land_idx.size
    log(f"Land cells (union over members): {n_land:,} of {LAT*LON:,}")

    # ---- Pack (member, year, land_cell) ------------------------------------ #
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        yrs, cube = cubes[f]
        flat = cube.reshape(cube.shape[0], -1)
        for k, y in enumerate(yrs):
            yi = y_index.get(int(y))
            if yi is not None:
                annual[s][slot[s][m], yi] = flat[k, land_idx]
        del cube, flat
    del cubes
    log(f"Packed {len(files)} members "
        f"({sum(a.nbytes for a in annual.values())/1024**2:.0f} MB resident)")

    # ---- Field nature: measured, never assumed (GUARDRAILS 9) -------------- #
    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> is_boolean_field={boolean}")
    if boolean:
        log("  ERROR: `let` must be continuous; a boolean read means the input is wrong.")
        return 3
    zf = float(np.mean(annual[scenarios[0]][np.isfinite(annual[scenarios[0]])] == 0))
    log(f"  annual cell-year exact-zero fraction: {zf:.2%}")
    log(f"  decadal statistic: {STAT_NAME} (central='{CENTRAL}') -- see the module "
        "docstring, decision 1")

    # ---- Shared 2020s baseline --------------------------------------------- #
    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios; the shared "
            "baseline is only valid for a uniform ensemble.")
    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    base_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years[bwin], BASELINE_DECADE, window_years=WINDOW_YEARS,
        central=CENTRAL)
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, CI_CAP)
    b_hi = np.clip(b_hi, CI_FLOOR, CI_CAP)

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"\nShared 2020s baseline: land n={baseline_flat.size:,}, "
        f"exact-zero {frac_zero:.2%}, percentile mode={pct_mode}, "
        f"exposed cells {int((baseline_flat > 0).sum()):,}, "
        f"global-mean exposed fraction {np.nanmean(b_med):.6f}, "
        f"max {np.nanmax(b_med):.4f}")

    # ---- Per-member diagnostic for the dashboard "Members" tab ------------- #
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
            "variable": OUT_VAR,
            "units": "1",
            "member_field": (f"mean annual exposed land-area fraction over the "
                             f"{BASELINE_DECADE}s baseline decade, pooled across "
                             "scenarios, AFTER spatial smoothing"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Consumed "
                     "by scripts/generate_maps.py to render the Members tab. All four "
                     "members share one impact model (ke-tg-meanfield), so differences "
                     "here are pure inter-GCM storm-field spread."),
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

    for s_drop in [s for s in scenarios if s not in write_scenarios]:
        del annual[s_drop]

    chunk = slope_chunk_cells(max(len(m) for m in members_by_scen.values()),
                              n_years, args.max_pairs)
    log(f"Theil-Sen chunk_cells={chunk}")

    models = sorted({meta[f]["model"] for f in files})
    gcms = sorted({meta[f]["gcm"] for f in files})
    socs = sorted({meta[f]["soc"] for f in files})
    senss = sorted({meta[f]["sens"] for f in files})

    # ---- Per-scenario assembly --------------------------------------------- #
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
                    cube, years, d, window_years=WINDOW_YEARS, central=CENTRAL)
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
                if d <= BASELINE_DECADE:
                    slope_txt = "slopes=NaN (at/before baseline)"
                else:
                    slope_txt = (f"ols={np.nanmean(out['ols_slope'][i]):+.3e} "
                                 f"sen={np.nanmean(out['sen_slope'][i]):+.3e} /dec")
                tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
                log(f"  {d}s: {tag:<15} mean={np.nanmean(out['median'][i]):.6f}  "
                    f"{slope_txt}  [{time.time() - td:.0f}s]")

        # ---- GUARDRAIL: slope and median masks must agree ------------------- #
        for i, d in enumerate(DECADES):
            if d <= BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), \
                    f"ols must be NaN at {d}s (no elapsed period from the baseline)"
                assert np.all(np.isnan(out["sen_slope"][i])), \
                    f"sen must be NaN at {d}s"
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

        sen_zero = float(np.mean(out["sen_slope"][-1][np.isfinite(out["sen_slope"][-1])] == 0))

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Land area fraction exposed to tropical cyclones",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": STAT_NAME,
                "field_nature": "continuous_zero_inflated",
                "value_note": (
                    "median = MEAN over the pooled (year x member) sample inside the "
                    f"decade window, across {len(mem)} members (ke-tg-meanfield x 4 CMIP5 "
                    "GCMs), of the annual fraction of each cell's LAND AREA exposed to "
                    "tropical cyclones. Raw ISIMIP values are fractional in [0,1) -- NOT "
                    "a binary flag (measured max 0.952, never exactly 1)."),
                "decadal_statistic_rationale": (
                    "DECLARED DEVIATION from the OUTPUT-SPEC default. The spec sends a "
                    "zero-inflated continuous field to the median/IQR branch, but this "
                    "field is 97.84% exact zeros at annual resolution, where the pooled "
                    "median is 0 on 96.07% of land and leaves only 2,684 exposed cells "
                    "against 15,122 under the mean -- 93% of exposed land erased. `let` is "
                    "the continuous analogue of `led`, which the spec already grants the "
                    "mean for the same reason: both estimate an expected annual exposed "
                    "fraction, and a median of near-Bernoulli draws is meaningless. "
                    "User decision 2026-08-11. Contrast burntarea at 29.2% zeros, which "
                    "takes the median branch without difficulty."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clamped to [{CI_FLOOR}, {CI_CAP}] since the "
                    "value is an area fraction. The band carries BOTH interannual "
                    "variability and inter-GCM spread; with a single impact model there is "
                    "no inter-model component at all."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "fitted over an EXPANDING window from the start of the 2020s baseline "
                    "through the end of the target decade, stacking every (year, member) "
                    "observation as an independent point. The 2010s and 2020s panels are "
                    "NaN (no elapsed period from the baseline). The estimators fail in "
                    "OPPOSITE regimes, so disagreement means a cell's trend is not robust. "
                    f"This field is {frac_zero:.1%} exact zeros in the baseline, i.e. "
                    "STRONGLY ZERO-INFLATED, and sen_slope is exactly 0 on "
                    f"{sen_zero:.1%} of finite cells in the final panel -- READ ols_slope "
                    "ON THIS LAYER."),
                "slope_units": "1 decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal value ranked against the shared "
                    "2020s ensemble land distribution. Zeros (never reached by a cyclone) "
                    "map to 1; positives rank against the NON-ZERO baseline into [2,100]. "
                    "Not inverted -- more exposure is worse."),
                "percentile_zero_fraction": round(frac_zero, 5),
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
                    "none -- a single impact model reporting a dimensionless area "
                    "fraction; there are no cross-model scales to reconcile."),
                "spatial_smoothing": (
                    "none" if kernel is None else
                    f"{2*SMOOTH_RADIUS+1}x{2*SMOOTH_RADIUS+1} window, w=exp(-d/{SMOOTH_L}) "
                    "grid-cells, cos(lat)-scaled longitude distance, normalized over "
                    "non-NaN land neighbours, longitude-wrapped, applied to EACH member's "
                    "ANNUAL map before any pooling. Ocean is never filled."),
                "smoothing_rationale": (
                    "Raw `let` is sparse storm TRACKS one cell wide, and a 4-member "
                    "ensemble cannot average them into a projected impact field. L=2.5 "
                    "puts 8.1% of the weight on the centre cell (the previous L=0.7 kernel "
                    "kept 32.1% and left the tracks visible), cutting the roughness index "
                    "from 0.389 raw / 0.142 previous to 0.044 while adding NO newly "
                    "exposed cell. The 2-cell radius is 111 km at the equator, about the "
                    "radius of hurricane-force winds; do not widen past 3 cells. Measured "
                    "mass change +3.1% (coastal redistribution). User decision 2026-08-11."),
                "ensemble_note": (
                    "SINGLE impact model (ke-tg-meanfield), so n_models is 1 everywhere "
                    "and the CI is inter-GCM spread only -- there is no structural "
                    "impact-model uncertainty represented in this layer."),
                "source_dataset": (
                    "ISIMIP2b DerivedOutputData/Lange2020/KE-TG-meanfield. No ISIMIP3b/SSP "
                    "re-issue of TC exposure exists (the family was re-issued as "
                    "Heinicke2026 driedarea/floodedarea and Zantout2025 heatwave/wildfire/"
                    "cropfailure; ke-tg-meanfield appears in no 3a or 3b directory). "
                    "ISIMIP3b's newer TC product is per-storm wind footprints under "
                    "InputData/climate/tropical_cyclones/MIT -- a different quantity, "
                    "noted for a future delivery."),
                "description": (
                    "Tropical-cyclone exposed land-area fraction processed to the TCFD "
                    "output contract (OUTPUT-SPEC.md) with a shared 2020s baseline; "
                    f"{len(mem)}-member single-model x CMIP5-GCM ensemble, pooled-mean "
                    "decadal statistic (declared deviation for extreme zero-inflation), "
                    "aggressive 5x5 spatial smoothing, two-tier percentile, "
                    "higher_is_worse."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)  "
            f"sen_slope exactly 0 on {sen_zero:.1%} of finite cells (final panel)")

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
