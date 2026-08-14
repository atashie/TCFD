#!/usr/bin/env python3
"""Process ISIMIP3b TipESM2025 `fldfrcmax` (CaMa-Flood inundation) onto the TCFD contract.

    python scripts/process_fldfrcmax_isimip3b.py --prepare-mask       # once, first
    python scripts/process_fldfrcmax_isimip3b.py --protection flopros # second, writes the reference
    python scripts/process_fldfrcmax_isimip3b.py --protection 40yr
    python scripts/process_fldfrcmax_isimip3b.py --protection none

THE ORDER IS NOT A CONVENIENCE. `--prepare-mask` establishes the permanent-water mask and
the land grid; the `flopros` run then writes the percentile reference every other variant
ranks against. Running `none` or `40yr` first fails loudly rather than silently inventing a
reference.

WHAT THIS LAYER IS
------------------
Annual MAXIMUM flooded fraction of a 15 arcmin (0.25 deg, ~28 km) cell, from CaMa-Flood
v4.0.0 routing seven ISIMIP3b GHMs. It is an ABSOLUTE inundation field, which is the whole
reason it was chosen: the 0.5 deg `floodedarea` exposure flag scores departure from a
preindustrial reference and therefore reads 0.000 across the Amazon floodplain in all 45 of
its members. Measured here on the full 258-member ingest, the Amazon main stem box reads
0.11496 unprotected -- the third-wettest of eight reference boxes.

Because the value is absolute, this layer does NOT carry the `relative_baseline` caveat that
`drought-2b`, `drought-3b` and `cropfailure-3b` require. It carries a different one; see
PROTECTION below.

THREE VARIANTS, THREE LAYERS, ONE PERCENTILE REFERENCE (user decisions 2026-08-14)
---------------------------------------------------------------------------------
`none`, `40yr` and `flopros` ship as three separate registry entries -- OUTPUT-SPEC has no
protection dimension -- with `flopros` the default for customer deliveries.

  * EVERY RAW VALUE IS THE VARIANT'S OWN: median, lower_ci, upper_ci, ols_slope and
    sen_slope are computed from that variant's members and its own within-decade values.
  * ONLY THE PERCENTILE REFERENCE IS SHARED: all three rank against the flopros 2020s
    baseline distribution, so a percentile means the same thing in all three.

Measured cost of the shared reference (8 members, ssp126): `none` ranked against flopros
2020s has a median percentile of 89.7 with 11.1% of cells at >=99, against a median of 50.0
and 1.0% on its own baseline. So the unprotected variant loses top-end resolution -- that is
the accepted price of cross-variant comparability, and it must be stated wherever the
`none` percentile is shown.

What is NOT shared is the 2020s VALUES. Substituting flopros 2020s values into the `none`
layer would put 0.0055 in its 2020s panel and 0.039 in its 2030s -- a 7x discontinuity with
no physical process behind it.

PROTECTION IS A MODELLING CHOICE, AND ZERO DOES NOT MEAN SAFE
------------------------------------------------------------
`flopros` applies the empirical FLOPROS flood-protection standards. Measured: the Rhine
(NL/DE) box falls from 0.05909 unprotected to 0.00122, and the Lower Mississippi from
0.14505 to 0.00358. A defended site reading ~0 means PROTECTED TO STANDARD, not NO HAZARD,
and protection is held constant into the future, so it embeds a no-further-adaptation
assumption. That is a must-disclose caveat for any flopros delivery.

`40yr` and `flopros` are NOT ordered versions of each other -- the ranking reverses by
region (Amazon 0.01077 vs 0.01607; Ganges 0.06853 vs 0.05664), because a uniform 40-year
standard exceeds the empirical standard in some countries and falls short in others.

MEASURED PROPERTIES (all 258 ingested members, 2026-08-14)
----------------------------------------------------------
  * CONTINUOUS [0,1], not binary.
  * Zero share DIFFERS BY VARIANT: none 17.05%, flopros 83.11%, 40yr 87.27%. The decadal
    statistic branch is therefore a per-variant measurement, not a layer-wide choice, and
    this script refuses to deviate from the median branch without an explicit --central.
  * MASK IS UNIFORM: finite over exactly 24.34% of the grid in every member, min == max.
    Unlike `led` / `driedarea` / `cropfailure` there is NO heterogeneous coverage, so there
    is no minimum-model rule to resolve. The script asserts this rather than assuming it.
  * PERMANENT WATER: cells fully flooded in every year of every member -- 511 unprotected
    (median; 494-531 across members), ~230 protected, ~373 of them in the Caspian. They are
    inland water bodies, not flood risk, and they would occupy the top of the percentile
    ranking. Masked out (user decision 2026-08-14), identified from the `none` variant
    because it is the most permissive and the physical set does not depend on defences.

TWO DEFECTS IN THE SOURCE, BOTH SILENT
--------------------------------------
  1. LATITUDE ORIENTATION DIFFERS BETWEEN PROTECTION VARIANTS. `none` is descending in all
     96 files, `40yr`/`flopros` ascending in all 162 -- same publication, same member. The
     variable's dim order is ('time','lat','lon') everywhere, so the arrays are NOT
     transposed; it is purely the coordinate order. A `.sel(lat=slice(hi, lo))` returns zero
     cells on whichever half does not match, and positional indexing builds an upside-down
     layer that passes every contract check. Every read here goes through
     `load_member()`, which sorts to DESCENDING lat to match every other shipped layer.
  2. Time is `days since 1661-01-01` (proleptic_gregorian) -- a fourth convention in this
     repository. Decoded with cftime through xarray; never by days/365 arithmetic.
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing as mp
import os
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.decadal_stats import (  # noqa: E402
    expanding_slopes, is_boolean_field, pooled_decadal_stat,
)

VAR = "fldfrcmax"
RAW_LAYER_ID = "flood-isimip3b_fldfrcmax_annual"
OUT_LAYER_TMPL = "flood-isimip3b_fldfrcmax-{prot}_annual"
PROTECTIONS = ("none", "40yr", "flopros")
REFERENCE_PROTECTION = "flopros"

DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

CI_FLOOR, CI_CAP = 0.0, 1.0          # a flooded FRACTION is bounded
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False             # more inundation = worse
SLOPE_PER_DECADE = 10.0
SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

#: Above this share of exact zeros in the baseline-decade median panel, the median branch
#: is collapsing and the script REFUSES to publish it without an explicit --central.
MEDIAN_COLLAPSE_THRESHOLD = 0.50

#: Sanity boxes. Caspian must be wetter than Sahara in COORDINATE space or a file's data
#: disagrees with its own axis. (lat_south, lat_north, lon_west, lon_east)
CASPIAN = (37.0, 46.0, 48.0, 54.0)
SAHARA = (23.0, 26.0, 8.0, 12.0)

_CUBE = None                          # set before forking the slope pool


def log(msg):
    print(msg, flush=True)


def artefact_dir(root: Path) -> Path:
    d = root / "data" / "processed" / "_fldfrcmax_shared"
    d.mkdir(parents=True, exist_ok=True)
    return d


def parse_name(fpath: str) -> dict:
    """Fields from the END. TipESM2025 has NO bias-adjustment token, so its field count
    differs from Zimmer2023's -- never carry an offset between publications."""
    # cama-flood_h08_gfdl-esm4_ssp585_2015soc_default_fldfrcmax_flopros_15arcmin_global_annual_2015_2100.nc
    #            [-12]   [-11]     [-10]   [-9]     [-8]     [-7]      [-6]     [-5]   [-4]  [-3]  [-2] [-1]
    # The `variable != VAR` assertion below is what caught these offsets being one out on
    # the first write -- keep it.
    p = os.path.basename(fpath).split("_")
    info = dict(model=p[-12], gcm=p[-11], scenario=p[-10], soc=p[-9], sens=p[-8],
                variable=p[-7], protection=p[-6], resolution=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def load_member(fpath: str):
    """(years, (year, lat, lon) DESCENDING lat, land_mask2d), NaN off-mask.

    Sorting is not cosmetic here: this publication publishes `none` with descending lat and
    the protected variants with ascending lat, so anything positional is wrong for half the
    files. Descending is chosen to match every other shipped layer.
    """
    ds = xr.open_dataset(fpath)                      # cftime decode
    da = ds[VAR].sortby("lat", ascending=False).sortby("lon")
    yrs = ds.time.dt.year.values.astype(int)
    vals = da.values.astype("float32")
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    lats, lons = da["lat"].values, da["lon"].values
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan

    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    vals, yrs = vals[keep], yrs[keep]

    valid = np.isfinite(vals)
    mask2d = valid.any(axis=0)
    return yrs, vals, mask2d, lats, lons


def box_mean(grid2d, lats, lons, b):
    la_lo, la_hi, lo_w, lo_e = b
    i = (lats >= la_lo) & (lats <= la_hi)
    j = (lons >= lo_w) & (lons <= lo_e)
    sub = grid2d[np.ix_(i, j)]
    return float(np.nanmean(sub)) if np.isfinite(sub).any() else np.nan


def assert_orientation(grid2d, lats, lons, name):
    """A file whose data disagrees with its own coordinates would build an upside-down
    layer that satisfies every line of OUTPUT-SPEC. Two box means catch it."""
    casp = box_mean(grid2d, lats, lons, CASPIAN)
    sah = box_mean(grid2d, lats, lons, SAHARA)
    if not (np.isnan(casp) and np.isnan(sah)) and not casp > sah:
        raise ValueError(
            f"{name}: data disagrees with its declared latitude axis "
            f"(Caspian {casp:.4f} <= Sahara {sah:.4f}) -- the array is flipped.")


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the SHARED flopros 2020s land distribution.

    Two-tier when the reference is materially zero-inflated: a cell that never floods
    across the pooled baseline sample is the lowest risk and maps to 1; positives rank
    against the NON-ZERO reference into [2,100]. The mode follows the measured zero share.
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
    # Plain arrays, not the SlopeResult dict subclass -- its __getattr__ turns pickle's
    # __getstate__ probe into a KeyError and kills the pool.
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


def files_for(raw_dir: Path, protection: str):
    pat = str(raw_dir / f"*_{VAR}_{protection}_*_global_annual_*.nc")
    return sorted(glob.glob(pat))


# --------------------------------------------------------------------------------------
# Stage 1 -- the shared mask
# --------------------------------------------------------------------------------------
def prepare_mask(root: Path) -> int:
    """Permanent-water mask and the published land grid, both from the `none` variant.

    `none` is used because permanent water is a physical property of the cell and does not
    depend on assumed defences -- and because it is the most permissive variant, so it
    identifies the largest such set (511 cells median, against ~230 under protection).
    """
    raw_dir = root / "data" / "raw" / RAW_LAYER_ID
    files = files_for(raw_dir, "none")
    if not files:
        log(f"ERROR: no `none` members in {raw_dir}")
        return 2

    log("=" * 78)
    log(f"prepare-mask: permanent water + land grid from {len(files)} `none` members")
    log("=" * 78)

    always = None
    land = None
    coverage = []
    lats = lons = None
    t0 = time.time()
    for n, f in enumerate(files, 1):
        yrs, vals, mask2d, la, lo = load_member(f)
        if lats is None:
            lats, lons = la, lo
        assert_orientation(np.nanmean(vals, axis=0), lats, lons, os.path.basename(f))
        a = np.all(np.nan_to_num(vals, nan=0.0) >= 1.0, axis=0) & mask2d
        always = a if always is None else (always & a)
        land = mask2d if land is None else (land & mask2d)
        coverage.append(int(mask2d.sum()))
        del vals
        if n % 16 == 0 or n == len(files):
            log(f"  {n}/{len(files)} members  [{time.time() - t0:.0f}s]")

    lo_cov, hi_cov = min(coverage), max(coverage)
    log(f"\nCoverage per member: min {lo_cov:,} max {hi_cov:,} cells")
    if lo_cov != hi_cov:
        log("  WARNING: coverage is NOT uniform across members. The measured ingest said it "
            "was (24.34% in every member); a minimum-model rule may now be needed.")
    else:
        log("  UNIFORM -- no minimum-model rule required for this layer.")

    keep = land & ~always
    log(f"Permanent water (>=1.0 in EVERY year of EVERY `none` member): {int(always.sum()):,} cells")
    log(f"  of which inside the Caspian box: "
        f"{int(always[np.ix_((lats >= CASPIAN[0]) & (lats <= CASPIAN[1]), (lons >= CASPIAN[2]) & (lons <= CASPIAN[3]))].sum()):,}")
    log(f"Published land cells: {int(keep.sum()):,} (was {int(land.sum()):,} before masking)")

    out = artefact_dir(root) / "fldfrcmax_mask.npz"
    np.savez_compressed(
        out, land=keep, permanent_water=always, lats=lats, lons=lons,
        n_members=len(files), coverage_min=lo_cov, coverage_max=hi_cov,
    )
    log(f"\nwrote {out}")
    return 0


# --------------------------------------------------------------------------------------
# Stage 2 -- process one protection variant
# --------------------------------------------------------------------------------------
def main() -> int:
    global _CUBE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prepare-mask", action="store_true")
    ap.add_argument("--protection", choices=PROTECTIONS)
    ap.add_argument("--scenarios", nargs="*", default=None)
    ap.add_argument("--central", choices=("median", "mean"), default=None,
                    help="decadal statistic branch; omit to auto-measure and be told")
    ap.add_argument("--limit-cells", type=int, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--skip-slopes", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    if args.prepare_mask:
        return prepare_mask(root)
    if not args.protection:
        log("ERROR: pass --prepare-mask or --protection {none,40yr,flopros}")
        return 2

    prot = args.protection
    raw_dir = root / "data" / "raw" / RAW_LAYER_ID
    out_dir = root / "data" / "processed" / OUT_LAYER_TMPL.format(prot=prot)
    out_dir.mkdir(parents=True, exist_ok=True)
    shared = artefact_dir(root)
    out_var = f"{VAR}-{prot}"

    mask_path = shared / "fldfrcmax_mask.npz"
    if not mask_path.exists():
        log(f"ERROR: {mask_path} missing. Run --prepare-mask first: the permanent-water "
            "mask must be identical across all three variants.")
        return 2
    mz = np.load(mask_path)
    keep2d, lats, lons = mz["land"], mz["lats"], mz["lons"]
    LAT, LON = len(lats), len(lons)
    land_idx = np.flatnonzero(keep2d.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        land_idx = land_idx[np.linspace(0, land_idx.size - 1, args.limit_cells).astype(int)]
    n_land = land_idx.size

    ref_path = shared / "fldfrcmax_percentile_reference.npz"
    is_reference_run = prot == REFERENCE_PROTECTION
    if not is_reference_run and not ref_path.exists():
        log(f"ERROR: {ref_path} missing. Process --protection {REFERENCE_PROTECTION} first: "
            "all three variants rank against its 2020s baseline (user decision 2026-08-14).")
        return 2

    files = files_for(raw_dir, prot)
    if not files:
        log(f"ERROR: no `{prot}` members in {raw_dir}")
        return 2
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    write_scenarios = args.scenarios or scenarios

    log("=" * 78)
    log(f"fldfrcmax [{prot}] -> TCFD contract  ({out_var})")
    log("=" * 78)
    log(f"{len(files)} files | scenarios {scenarios} | {n_land:,} published cells "
        f"| grid {LAT}x{LON} @ {abs(lats[1]-lats[0]):.2f} deg")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    y_index = {int(y): i for i, y in enumerate(years)}
    members_by_scen = {s: sorted({meta[f]["member"] for f in files
                                  if meta[f]["scenario"] == s}) for s in scenarios}
    models = sorted({meta[f]["model"] for f in files})
    if len({tuple(v) for v in members_by_scen.values()}) != 1:
        log("ERROR: ensemble composition differs across scenarios; a shared 2020s baseline "
            "would not be valid.")
        return 3
    log(f"members/scenario {len(members_by_scen[scenarios[0]])} | impact models {len(models)}: "
        f"{', '.join(models)}")

    # ---- Pass A: the 2020s window only, pooled across scenarios ---------------------- #
    # Full cubes for three scenarios would be ~8 GB on a 17 GB machine; the baseline needs
    # only ten years, so it is built first and the scenarios are then streamed one at a time.
    t0 = time.time()
    bwin_years = np.arange(BASELINE_DECADE, BASELINE_DECADE + WINDOW_YEARS)
    base_stack = []
    for f in files:
        yrs, vals, _, la, lo = load_member(f)
        sel = np.isin(yrs, bwin_years)
        flat = vals[sel].reshape(sel.sum(), -1)[:, land_idx]
        base_stack.append(flat)
        del vals
    base_pool = np.stack(base_stack, axis=0)          # (member*scenario, year, cell)
    del base_stack
    log(f"\nBaseline pool {base_pool.shape} ({base_pool.nbytes/1e9:.2f} GB) "
        f"[{time.time() - t0:.0f}s]")

    # ---- Field nature and the statistic branch: measured, then declared -------------- #
    boolean = is_boolean_field(base_pool)
    fin = base_pool[np.isfinite(base_pool)]
    zero_share = float((fin == 0).mean())
    log(f"Field nature: {'BOOLEAN' if boolean else 'CONTINUOUS'} | annual exact-zero "
        f"share {zero_share:.2%}")
    if boolean:
        log("  ERROR: measured BOOLEAN. fldfrcmax was continuous at ingest -- re-run "
            "scripts/check_fldfrcmax_nature.py before processing.")
        return 3
    del fin

    med_test, _, _ = pooled_decadal_stat(base_pool, bwin_years, BASELINE_DECADE,
                                         window_years=WINDOW_YEARS, central="median")
    med_zero = float((med_test[np.isfinite(med_test)] == 0).mean())
    mean_test, _, _ = pooled_decadal_stat(base_pool, bwin_years, BASELINE_DECADE,
                                          window_years=WINDOW_YEARS, central="mean")
    mean_zero = float((mean_test[np.isfinite(mean_test)] == 0).mean())
    log(f"\nDecadal-statistic measurement on the {BASELINE_DECADE}s panel:")
    log(f"  median branch -> {med_zero:.2%} of published cells exactly zero")
    log(f"  mean   branch -> {mean_zero:.2%} of published cells exactly zero")
    central = args.central
    if central is None:
        if med_zero >= MEDIAN_COLLAPSE_THRESHOLD:
            log(f"\n  The median branch erases {med_zero:.1%} of published cells, above the "
                f"{MEDIAN_COLLAPSE_THRESHOLD:.0%} collapse threshold. That is a DECLARED "
                "deviation, not a default: re-run with --central mean to take the "
                "zero-inflated branch, or --central median to publish the median anyway.")
            return 4
        central = "median"
    stat_name = ("pooled_median" if central == "median" else "pooled_mean_zero_inflated")
    log(f"  branch: {stat_name} (central={central})")
    del med_test, mean_test

    b_med, b_lo, b_hi = pooled_decadal_stat(base_pool, bwin_years, BASELINE_DECADE,
                                            window_years=WINDOW_YEARS, central=central)
    b_lo = np.clip(b_lo, CI_FLOOR, CI_CAP)
    b_hi = np.clip(b_hi, CI_FLOOR, CI_CAP)
    del base_pool

    # ---- The percentile reference: written by flopros, read by everyone -------------- #
    # A cell-limited run must never write the reference: it would be built on a subsampled
    # cell set and would silently mis-rank every later variant. Such runs write nothing.
    dry_run = bool(args.limit_cells)
    if is_reference_run and dry_run:
        log("\n--limit-cells set: timing/validation run only. NOT writing the percentile "
            "reference and NOT writing any NetCDF; ranking against this run's own baseline.")
    if is_reference_run:
        if not dry_run:
            np.savez_compressed(ref_path, baseline=b_med, land_idx=land_idx,
                                central=central, stat_name=stat_name,
                                n_members=len(members_by_scen[scenarios[0]]),
                                scenarios=np.array(scenarios))
            log(f"\nwrote percentile reference {ref_path.name} "
                f"({np.isfinite(b_med).sum():,} cells, central={central})")
        ref_flat = b_med[np.isfinite(b_med)]
        ref_note = (f"this layer's OWN {BASELINE_DECADE}s baseline, which is also the shared "
                    "reference for the other protection variants")
    else:
        rz = np.load(ref_path, allow_pickle=False)
        if not np.array_equal(rz["land_idx"], land_idx):
            log("ERROR: the reference was built on a different cell set. Re-run "
                "--prepare-mask and reprocess flopros.")
            return 3
        ref = rz["baseline"]
        ref_flat = ref[np.isfinite(ref)]
        ref_note = (f"the {REFERENCE_PROTECTION} {BASELINE_DECADE}s baseline, NOT this "
                    "variant's own (user decision 2026-08-14)")
        log(f"\nreference: {REFERENCE_PROTECTION} {BASELINE_DECADE}s "
            f"({ref_flat.size:,} cells, central={str(rz['central'])})")

    pct, pct_mode, ref_zero = make_pct_fn(ref_flat)
    b_pct = pct(b_med)
    log(f"  percentile mode={pct_mode}, reference exact-zero {ref_zero:.2%}, "
        f"own {BASELINE_DECADE}s mean {np.nanmean(b_med):.5f}")

    # ---- Pass B: one scenario at a time ---------------------------------------------- #
    chunk = slope_chunk_cells(len(members_by_scen[scenarios[0]]), years.size, args.max_pairs)
    log(f"Theil-Sen chunk_cells={chunk}, jobs={args.jobs}")

    for s in write_scenarios:
        log(f"\n{'=' * 78}\nScenario {s}\n{'=' * 78}")
        mem = members_by_scen[s]
        cube = np.full((len(mem), years.size, n_land), np.nan, np.float32)
        slot = {m: i for i, m in enumerate(mem)}
        ts = time.time()
        for f in [f for f in files if meta[f]["scenario"] == s]:
            yrs, vals, _, _, _ = load_member(f)
            flat = vals.reshape(vals.shape[0], -1)
            for k, y in enumerate(yrs):
                yi = y_index.get(int(y))
                if yi is not None:
                    cube[slot[meta[f]["member"]], yi] = flat[k, land_idx]
            del vals, flat
        log(f"  packed {len(mem)} members ({cube.nbytes/1e9:.2f} GB) [{time.time()-ts:.0f}s]")
        _CUBE = cube

        out = {k: np.full((len(DECADES), LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}

        for i, d in enumerate(DECADES):
            td = time.time()
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(cube, years, d, window_years=WINDOW_YEARS,
                                                  central=central)
                lo = np.clip(lo, CI_FLOOR, CI_CAP)
                hi = np.clip(hi, CI_FLOOR, CI_CAP)
                pc = pct(med)

            if args.skip_slopes:
                sl = dict(ols_slope=np.full(n_land, np.nan, np.float32),
                          sen_slope=np.full(n_land, np.nan, np.float32))
            else:
                sl = compute_slopes(cube, years, d, chunk, args.max_pairs, args.jobs, n_land)

            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
            mod_present = np.zeros((len(models), n_land), bool)
            for mi, m in enumerate(mem):
                mod_present[models.index(m.rsplit("_", 1)[0])] |= present[mi]
            n_mod = mod_present.sum(axis=0).astype(np.float32)
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
                txt = ("slopes=NaN (baseline)" if d <= BASELINE_DECADE else
                       f"ols={np.nanmean(out['ols_slope'][i]):+.3e} "
                       f"sen={np.nanmean(out['sen_slope'][i]):+.3e} /dec")
                log(f"  {d}s: mean={np.nanmean(out['median'][i]):.5f}  {txt}  "
                    f"[{time.time()-td:.0f}s]")

        # ---- contract guardrails ---------------------------------------------------- #
        for i, d in enumerate(DECADES):
            if d <= BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), f"ols must be NaN at {d}s"
                assert np.all(np.isnan(out["sen_slope"][i])), f"sen must be NaN at {d}s"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), f"{k} finite where median is NaN at {d}s -- leak"
            assert np.all(out["lower_ci"][i][med_finite]
                          <= out["median"][i][med_finite] + 1e-5), f"lower_ci > median {d}s"
            assert np.all(out["upper_ci"][i][med_finite]
                          >= out["median"][i][med_finite] - 1e-5), f"upper_ci < median {d}s"
        assert_orientation(out["median"][-1], lats, lons, f"{out_var} {s} output")

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": out_var,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Annual maximum flooded fraction of the cell",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                # Both required by test_shared_baseline.py; omitting them failed the
                # contract test on the first flopros build.
                "slope_units": "1 decade-1",
                "baseline_source": "shared_across_all_scenarios",
                "grid": f"{LAT}x{LON} at {abs(lats[1]-lats[0]):.2f} deg (15 arcmin, ~28 km) "
                        "-- NOT the 0.5 deg of the other shipped layers",
                "spatial_resolution_degrees": round(float(abs(lats[1] - lats[0])), 6),
                "decadal_statistic": stat_name,
                "decadal_statistic_rationale": (
                    f"Measured on the {BASELINE_DECADE}s panel of THIS variant: the median "
                    f"branch leaves {med_zero:.2%} of published cells exactly zero, the mean "
                    f"branch {mean_zero:.2%}. Measured across all three variants 2026-08-14 "
                    "(median / mean exact-zero share of the panel): none 13.31% / 0.08%, "
                    "flopros 92.68% / 2.69%, 40yr 97.25% / 0.12%. "
                    "ALL THREE VARIANTS TAKE THE MEAN BRANCH, and for `none` that is a "
                    "deviation its own distribution does NOT require -- its median is "
                    "healthy at 13.31%, comfortably inside the range where `burntarea` "
                    "ships the median (29.2% annual zeros). The reason is the SHARED "
                    "PERCENTILE REFERENCE: all three variants rank against the flopros "
                    f"{BASELINE_DECADE}s panel, so a mixed choice would rank `none` MEDIANS "
                    "against a distribution of flopros MEANS -- a second incomparability on "
                    "top of the protection difference the shared reference exists to "
                    "neutralise. Uniform estimator, declared (user decision 2026-08-14). "
                    "This is NOT the branch being taken to improve contrast."),
                "field_nature": "continuous",
                "protection_level": prot,
                "protection_note": (
                    "PROTECTION IS A MODELLING CHOICE. `flopros` applies empirical flood-"
                    "protection standards held constant into the future, so it embeds a "
                    "no-further-adaptation assumption; a defended site reading ~0 means "
                    "PROTECTED TO STANDARD, not NO HAZARD (measured: Rhine NL/DE 0.05909 "
                    "unprotected -> 0.00122 flopros). `40yr` assumes a uniform 40-year "
                    "standard everywhere. The two are NOT ordered versions of each other -- "
                    "the ranking reverses by region."),
                "percentile_reference": ref_note,
                "percentile_mode": pct_mode,
                "percentile_direction": "higher_is_worse",
                "percentile_note": (
                    "ALL THREE protection variants rank against the flopros 2020s baseline "
                    "so a percentile means the same thing across them (user decision "
                    "2026-08-14). Raw values -- median, CIs and both slopes -- are this "
                    "variant's own. Measured cost: the `none` variant's percentile compresses "
                    "at the top (median 89.7, 11.1% of cells at >=99, against 50.0 and 1.0% "
                    "on its own baseline)."),
                "mask_rule": (
                    "Uniform member coverage (24.34% of grid in every member), so NO "
                    "minimum-model rule. Permanent water is EXCLUDED: cells at >=1.0 in "
                    "every year of every `none` member (inland water bodies, mostly the "
                    "Caspian) would otherwise occupy the top of the percentile ranking."),
                "relative_baseline_note": (
                    "This layer is ABSOLUTE inundation, not a departure from a historical "
                    "reference, so it does NOT carry the relative-baseline caveat that "
                    "drought-2b/3b and cropfailure-3b require."),
                "source": ("ISIMIP3b DerivedOutputData/TipESM2025 (CaMa-Flood v4.0.0), "
                           "MIT/{scenario}/{GHM}/{gcm}/"),
                "ensemble": (f"{len(mem)} members, {len(models)} impact models: "
                             f"{', '.join(models)}"),
                "smoothing": "none -- deep ensemble (27-32 members)",
                "slope_note": (
                    "Both slopes are published per OUTPUT-SPEC. Theil-Sen pairs are "
                    f"SUBSAMPLED at max_pairs={args.max_pairs} with a fixed seed (rng_seed=0, "
                    "deterministic). Measured 2026-08-14 on a 4,000-cell run of this "
                    "variant: `ols_slope` is IDENTICAL with and without the cap -- the cap "
                    "touches only the pairwise estimator -- while `sen_slope` differs at the "
                    "1e-7 level, which is immaterial on a field where Sen has already "
                    "collapsed (all decades read |sen| < 3e-6 against |ols| ~ 1e-4). Full "
                    "pairs cost 204s per 4,000 cells against 7s capped, i.e. ~32 h versus "
                    "~1 h for the three variants. Read `ols_slope` on this layer."),
                "processed_on": time.strftime("%Y-%m-%d"),
            },
        )
        enc = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                   "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{out_var}_{s}_processed.nc"
        if dry_run:
            log(f"  --limit-cells: NOT writing {path.name} (subsampled cells)")
        else:
            ds_out.to_netcdf(path, encoding=enc)
            log(f"  wrote {path.name}")
        del cube, out
        _CUBE = None

    log(f"\nDone: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
