"""Process `driedarea` (ISIMIP3b/SSP drought exposure) into the TCFD output contract.

`driedarea` (Heinicke2026, ISIMIP3b DerivedOutputData) is the SSP re-issue of the Lange 2020
drought-exposure concept, renamed by hazard word rather than `le*` code. Each cell carries a
per-year flag for whether it was in drought.

This is a SEPARATE LAYER from `led` (ISIMIP2b, rcp26/rcp60), not a replacement: different
impact models, a different GCM generation, and different scenarios. Both ship.

MEASURED data nature (GUARDRAILS 9 -- from the values, never the name), all 45 members:

  * STRICTLY BINARY {0,1}: exactly 2 unique values per member, exact-0 plus exact-1 = 100%.
    `is_boolean_field` -> True, so the decadal statistic is the pooled MEAN = drought
    FREQUENCY. Verified independently of `led`; the two agree, but that was not assumed.
  * 75.55%-95.79% of finite cell-years are exact zeros, depending on member.
  * `long_name` is the generic "Exposed Area Share" -- it does not even say drought, which
    is another reason the nature had to come from values.
  * units "1"; grid 360x720; 86 records, one per year, 2015-2100.
  * THE MASK IS CLEAN. This mattered: the sibling `floodedarea` -- same publication, same
    models -- is non-NaN over 94.7% of the globe INCLUDING OCEAN. `driedarea` is 20.4-22.2%,
    i.e. genuine land only. A per-product mask assumption would have been wrong; it is
    per-VARIABLE.
  * Land masks are TIME-INVARIANT per member but differ BY GCM WITHIN a model (h08 spans
    52,820-53,735 cells across its five GCMs), not merely between models.
  * 15 members per scenario, composition IDENTICAL across ssp126/ssp370/ssp585 -- the
    enumerated matrix is complete (3 x 5 x 3 = 45, no gaps), so the shared 2020s baseline
    is valid. soc/sens uniformly 2015soc/default -- no harmonization compromise.

TIME DECODING: the axis is `days since 1601-01-01` on a `proleptic_gregorian` calendar --
NOT the `years since` axis the ISIMIP2b Lange2020 files use. It is decoded with cftime via
xarray (`decode_times=True`) and the year read off `.dt.year`; dividing days by 365 would
drift ~0.4 years over the 400-year offset and silently mis-bin decades.

TWO DECISIONS SPECIFIC TO THIS LAYER. Neither is inherited from `led`; both were measured
here first.

1. NO MINIMUM-MODEL MASK -- the full union of 63,455 cells is published.

   The `led` layer publishes only where >= 2 of its 8 models have data, because its
   single-model cells read 1.63x the all-model level: one high model (orchidee, 0.0660)
   undiluted by the low ones (h08 0.0085, watergap2 0.0123) across a 7.8x inter-model
   spread. That mechanism is ABSENT here, measured on the 2020s panel:

       coverage tier    cells      2020s mean freq
       1 model          5,396          0.0745
       2 models         8,629          0.0599
       3 models        49,430          0.0722

   The step is 1.03x, inter-model spread is only 2.69x (h08 0.0582, jules-w2 0.1118,
   watergap2-2e 0.0416), and the single-model cells are split across all three models
   (1,130 / 2,968 / 1,298) rather than dominated by one outlier. There is no composition
   artifact to remove, so a mask would cost 5,396 cells (8.5% of the union) and buy nothing.
   `--min-models` remains available, and `n_models` is emitted per cell either way.

   CAVEAT recorded rather than masked: a 1-model cell still carries NO inter-model
   uncertainty in its CI, whatever its level. Filter on `n_models` if that matters.

2. NO SPATIAL SMOOTHING. 15 members x 10 years = up to 150 Bernoulli draws per cell-decade.
   This is a broad, spatially coherent field -- not the one-cell-wide storm tracks that
   forced an aggressive kernel onto the 4-member `let` layer -- so smoothing would blur real
   aridity and permafrost gradients for no variance benefit.

Percentile: ranked against the shared 2020s land distribution, `higher_is_worse`. Tier mode
is chosen from the MEASURED zero-fraction of that baseline, not assumed.

Expect `sen_slope` to collapse to exactly 0 across the domain: on a binary field every
pairwise difference is -1, 0 or +1 and same-valued pairs dominate. READ `ols_slope`.

Usage:
    python scripts/process_driedarea_isimip3b.py [--scenarios ssp126] [--skip-slopes]
                                                 [--limit-cells N] [--members-only]
                                                 [--jobs N] [--min-models K]
"""

import argparse
import glob
import multiprocessing as mp
import os
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

VAR = "driedarea"
OUT_VAR = "driedarea"
LAYER_ID = "drought-isimip3b_driedarea_annual"
#: ISIMIP3b starts in 2015, so the baseline decade IS decade index 0 -- unlike the ISIMIP2b
#: `led`/`let` layers, which carry a full 2010s panel first. 2015-2019 cannot form a decade.
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

#: The decadal value is a FREQUENCY (fraction of model-years in drought), bounded by [0, 1].
CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False          # more drought = worse -> percentile NOT inverted
SLOPE_PER_DECADE = 10.0

#: See decision 1 -- measured, not inherited from `led`. 1 == publish the full union.
MIN_MODELS = 1

#: Measured boolean -> the contract's `pooled_mean_boolean` branch (mean +/- 1 SD).
STAT_NAME = "pooled_mean_boolean"

#: Three distinct models, no two configurations of one code, so family == model.
MODEL_FAMILY = {}

SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

#: Set in the parent before forking the slope pool; workers inherit it copy-on-write.
_CUBE = None


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) from a Heinicke2026 filename.

    e.g. h08_gfdl-esm4_w5e5_ssp126_2015soc_default_driedarea_global_annual_2015_2100.nc

    NOTE these filenames carry NO leading publication token, unlike ISIMIP2b's
    `lange2020_...`. DerivedOutputData grammar is per-PUBLICATION, so the fields are read
    from the END, which is index-stable under either convention.
    """
    p = os.path.basename(fpath).split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def load_member(fpath):
    """Load one member as (years, (year, lat, lon), land_mask2d), NaN over ocean.

    Time is decoded with cftime through xarray rather than by days-per-year arithmetic:
    the axis is `days since 1601-01-01` (proleptic_gregorian), and dividing by 365 drifts
    by ~0.4 yr over that offset, which is enough to push a January stamp into the wrong
    decade bin.
    """
    ds = xr.open_dataset(fpath)                    # decode_times=True
    da = ds[VAR]
    yrs = ds.time.dt.year.values.astype(int)
    vals = da.values.astype("float32")
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan

    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    vals, yrs = vals[keep], yrs[keep]

    valid = np.isfinite(vals)
    mask2d = valid.any(axis=0)
    if not np.array_equal(valid, np.broadcast_to(mask2d, valid.shape)):
        raise ValueError(f"{os.path.basename(fpath)}: land mask varies by year")
    return yrs, vals, mask2d


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline land distribution.

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros): a cell never
    in drought across 150 pooled model-years is the LOWEST risk, so it maps to 1 and
    positives rank against the NON-ZERO baseline into [2,100]. Mode is chosen from the
    measured share, not assumed.
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
            res = np.ones(vals.shape, np.float32)      # never in drought -> raw 1
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
    """Chunk width that keeps Theil-Sen peak memory inside the per-worker budget."""
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    if max_pairs is not None:
        pairs = min(pairs, max_pairs)
    per_cell = 4 * pairs * 4
    return max(4, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def _slope_block(task):
    """Worker: expanding slopes for one contiguous block of land cells."""
    s, e, years, decade, baseline, window, chunk, max_pairs = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk, max_pairs=max_pairs)
    # Return plain arrays, NOT the SlopeResult: it is a dict subclass whose
    # `__getattr__ = dict.__getitem__` turns pickle's probe for `__getstate__` into a
    # KeyError instead of an AttributeError, which kills the whole pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, max_pairs, jobs, n_land):
    """expanding_slopes over all land cells, fanned across `jobs` forked workers."""
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
    ap.add_argument("--scenarios", nargs="*", default=None,
                    help="subset of scenarios to WRITE (the baseline always pools all)")
    ap.add_argument("--limit-cells", type=int, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--members-only", action="store_true")
    ap.add_argument("--skip-slopes", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--min-models", type=int, default=MIN_MODELS,
                    help="publish a cell only where this many impact models have data; "
                         "1 = full union (this layer's measured default, see docstring)")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    out_dir = root / "data" / "processed" / LAYER_ID
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc")))
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
    log(f"driedarea (ISIMIP3b/SSP drought exposure) -> TCFD contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")

    with xr.open_dataset(files[0]) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    LAT, LON = len(lats), len(lons)

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

    # ---- Pass 1: per-MODEL land masks -> publication mask -------------------- #
    t0 = time.time()
    model_mask = {m: np.zeros((LAT, LON), bool) for m in models}
    for f in files:
        with xr.open_dataset(f) as ds:
            da = ds[VAR]
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
            v = da.values.astype("float32")
        if fill is not None:
            v = np.where(np.isclose(v, np.float32(fill), rtol=1e-6), np.nan, v)
        model_mask[meta[f]["model"]] |= np.isfinite(v).any(axis=0)
        del v

    nmod_static = np.sum([model_mask[m] for m in models], axis=0).astype(np.int16)
    union = nmod_static > 0
    keep2d = nmod_static >= args.min_models
    log(f"\nLand mask (per-model coverage, {time.time() - t0:.0f}s):")
    for k in range(1, len(models) + 1):
        c = int((nmod_static == k).sum())
        if c:
            log(f"  exactly {k} model(s): {c:>7,}")
    log(f"  union {int(union.sum()):,} -> publishing >= {args.min_models} models: "
        f"{int(keep2d.sum()):,} cells; dropped {int((union & ~keep2d).sum()):,}")

    land_idx = np.flatnonzero(keep2d.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        land_idx = land_idx[np.linspace(0, land_idx.size - 1,
                                        args.limit_cells).astype(int)]
    n_land = land_idx.size

    # ---- Pass 2: pack (member, year, land_cell) per scenario ----------------- #
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        yrs, cube, mask2d = load_member(f)
        flat = cube.reshape(cube.shape[0], -1)
        for k, y in enumerate(yrs):
            yi = y_index.get(int(y))
            if yi is not None:
                annual[s][slot[s][m], yi] = flat[k, land_idx]
        del cube, flat
        log(f"  loaded {info['model']:<14} {info['gcm']:<15} {s}  "
            f"{int(mask2d.sum()):,} land cells  [{time.time() - t0:.0f}s]")
    log(f"Packed {len(files)} members over {n_land:,} cells "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident)")

    # ---- Field nature: measured, never assumed (GUARDRAILS 9) ---------------- #
    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> is_boolean_field={boolean}")
    if not boolean:
        log("  ERROR: `driedarea` measured BINARY at ingest 2026-08-11. A continuous read "
            "means the inputs changed -- re-run the value check before processing.")
        return 3
    fin0 = annual[scenarios[0]][np.isfinite(annual[scenarios[0]])]
    log(f"  annual cell-year exact-zero fraction: {float((fin0 == 0).mean()):.2%}")
    log(f"  decadal statistic: {STAT_NAME} -- pooled MEAN = drought frequency")
    del fin0

    # ---- Shared 2020s baseline ----------------------------------------------- #
    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios; the shared "
            "baseline is only valid for a uniform ensemble.")
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
    log(f"\nShared {BASELINE_DECADE}s baseline: land n={baseline_flat.size:,}, "
        f"exact-zero {frac_zero:.2%}, percentile mode={pct_mode}, "
        f"mean frequency {np.nanmean(b_med):.4f}, max {np.nanmax(b_med):.4f}")

    # ---- Per-member diagnostic for the dashboard "Members" tab --------------- #
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
            "member_field": (f"drought frequency over the {BASELINE_DECADE}s baseline "
                             "decade (mean of the binary annual flag), pooled across "
                             "scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Inter-model "
                     "spread is 2.69x (watergap2-2e 0.0416, h08 0.0582, jules-w2 0.1118), "
                     "so members are expected to differ in level; look for spatial-pattern "
                     "outliers and mask differences, not level agreement."),
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
    jobs = max(1, args.jobs)
    log(f"Theil-Sen chunk_cells={chunk}, jobs={jobs}")

    # ---- Per-scenario assembly ----------------------------------------------- #
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
                    slope_txt = "slopes=NaN (baseline)"
                else:
                    slope_txt = (f"ols={np.nanmean(out['ols_slope'][i]):+.3e} "
                                 f"sen={np.nanmean(out['sen_slope'][i]):+.3e} /dec")
                tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
                log(f"  {d}s: {tag:<15} mean={np.nanmean(out['median'][i]):.4f}  "
                    f"{slope_txt}  [{time.time() - td:.0f}s]")

        # ---- GUARDRAIL: slope and median masks must agree -------------------- #
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
            assert np.nanmin(out["n_models"][i]) >= args.min_models, \
                f"a published cell at {d}s has fewer than {args.min_models} models"

        fin_sen = out["sen_slope"][-1][np.isfinite(out["sen_slope"][-1])]
        sen_zero = float((fin_sen == 0).mean()) if fin_sen.size else float("nan")

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Land area fraction exposed to drought",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": STAT_NAME,
                "field_nature": "boolean",
                "value_note": (
                    "median = MEAN over the pooled (year x member) sample inside the decade "
                    f"window, across {len(mem)} members (3 GHM x 5 CMIP6 GCMs), i.e. the "
                    "DROUGHT FREQUENCY: the fraction of model-years the cell was flagged as "
                    "in drought. The raw ISIMIP field is a BINARY {0,1} annual flag "
                    "(measured: exactly 2 unique values in all 45 members). Its published "
                    "long_name is the generic 'Exposed Area Share', which names neither the "
                    "hazard nor the per-cell nature -- both were taken from the values."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clamped to [{CI_FLOOR}, {CI_CAP}] since the "
                    "value is a frequency. Carries interannual, inter-GCM AND inter-model "
                    "spread; inter-model is 2.69x here (watergap2-2e 0.0416, h08 0.0582, "
                    "jules-w2 0.1118 global-mean frequency)."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "fitted over an EXPANDING window from the start of the 2020s baseline "
                    "through the end of the target decade, stacking every (year, member) "
                    "observation as an independent point. The 2020s panel is NaN (no "
                    "elapsed period). The estimators fail in OPPOSITE regimes, so "
                    "disagreement means a cell's trend is not robust. This field is BINARY, "
                    "so every pairwise difference is -1, 0 or +1 and same-valued pairs "
                    f"dominate: sen_slope is exactly 0 on {sen_zero:.1%} of finite cells in "
                    "the final panel. READ ols_slope ON THIS LAYER."),
                "slope_units": "1 decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal frequency ranked against the shared "
                    f"{BASELINE_DECADE}s ensemble land distribution. Not inverted -- more "
                    "drought is worse."),
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
                    "none -- every member reports the same dimensionless binary exposure "
                    "flag. The inter-model level spread is structural uncertainty, not a "
                    "unit mismatch, and belongs in the CI."),
                "spatial_smoothing": (
                    "none -- 15 members x 10 years gives up to 150 Bernoulli draws per "
                    "cell-decade, and this is a broad spatially coherent field rather than "
                    "the one-cell-wide storm tracks that forced a kernel onto `let`. "
                    "Smoothing would blur real aridity and permafrost gradients."),
                "minimum_models": args.min_models,
                "mask_rule": (
                    f"Publishing where >= {args.min_models} impact model(s) have data: "
                    f"{int(keep2d.sum()):,} of {int(union.sum()):,} union land cells. "
                    "MEASURED on the 2020s panel, the single-model cells show NO level "
                    "step -- 1 model 0.0745 (n=5,396), 2 models 0.0599 (n=8,629), 3 models "
                    "0.0722 (n=49,430), a ratio of 1.03x -- because inter-model spread is "
                    "only 2.69x and the solo cells are split across all three models "
                    "(1,130/2,968/1,298), not dominated by one outlier. So no minimum-model "
                    "cut is applied. This DIFFERS from the ISIMIP2b `led` layer, which "
                    "masks at >= 2 models because its solo cells read 1.63x across a 7.8x "
                    "inter-model spread; the rule was re-measured here, not inherited. "
                    "Caveat: a 1-model cell still carries no inter-model uncertainty in its "
                    "CI regardless of level -- filter on n_models if that matters."),
                "interpretation_caveat": (
                    "THIS IS DEPARTURE FROM A FIXED REFERENCE, NOT ARIDITY. The exposure "
                    "flag is a soil-moisture threshold relative to a fixed historical "
                    "reference distribution, so a permanently arid cell with a stable "
                    "regime scores LOW while a wet cell whose regime has shifted scores "
                    "HIGH. Do not read this layer as a dryness map."),
                "coverage_geography": (
                    "THIN COVERAGE CONCENTRATES IN THE ARID BELT. The 1- and 2-model cells "
                    "sit at a median |latitude| of 34 and 33 degrees, against 46 for the "
                    "3-model cells: hydrological models disagree about desert cells (no "
                    "runoff, no soil column), so the subtropical dry zone is where they "
                    "drop out. Mean n_models by region on the 2020s panel -- Sahara 1.57, "
                    "Arabian 1.87, Gobi/central Asia 2.10, Australian interior 2.28, "
                    "against Amazon 3.00, Europe/Mediterranean 2.77, global 2.69. "
                    "Reassuringly the LARGEST projected changes land on the best-covered "
                    "cells (ssp585 2090s minus 2020s: Amazon +0.436 at n_models 3.00, "
                    "Europe/Med +0.323 at 2.77) while the thin arid cells barely move "
                    "(Sahara +0.034, Australian interior +0.024, Arabian -0.017). So the "
                    "headline signal does not rest on thin coverage -- but a site-level "
                    "query in a desert should read n_members/n_models before the value."),
                "single_model_dominance": (
                    "THE HIGH-LATITUDE SIGNAL IS LARGELY ONE MODEL. In the 60-80N band the "
                    "2020s frequency reads jules-w2 0.2397, watergap2-2e 0.0453, h08 0.0272 "
                    "-- an 8.8x spread against a 2.69x global one -- so the ensemble's "
                    "0.1028 there is jules-w2 carrying the band. Above 60N jules-w2 is "
                    "2.14x the ensemble mean, versus 1.52x globally. This is visible as a "
                    "bright Arctic band in that model's panels on the Members tab. It is "
                    "structural model disagreement, not a defect, and the CI widens "
                    "accordingly -- but do not report a boreal drought trend from this "
                    "layer without saying it is one model of three."),
                "ensemble_note": (
                    "3 impact models x 5 CMIP6 GCMs = 15 members per scenario, COMPLETE "
                    "(3 x 5 x 3 = 45 files enumerated 2026-08-11, no missing combination). "
                    "Land masks differ BY GCM WITHIN a model, not only between models "
                    "(h08 spans 52,820-53,735 cells across its five GCMs), which is why "
                    "n_members is emitted per cell alongside n_models."),
                "source_dataset": (
                    "ISIMIP3b DerivedOutputData/Heinicke2026 -- the SSP re-issue of the "
                    "Lange 2020 drought-exposure concept, named by hazard word rather than "
                    "le* code. SEPARATE LAYER from the ISIMIP2b `led` drought layer "
                    "(8 models x 4 CMIP5 GCMs, rcp26/rcp60), not a replacement: different "
                    "models, GCM generation and scenarios. Sibling `floodedarea` in this "
                    "same publication has a broken land mask (non-NaN over 94.7% of the "
                    "globe including ocean); `driedarea` was checked and is clean."),
                "description": (
                    "ISIMIP3b/SSP drought exposure frequency processed to the TCFD output "
                    "contract (OUTPUT-SPEC.md) with a shared 2020s baseline; "
                    f"{len(mem)}-member 3-model x 5-GCM CMIP6 ensemble, pooled-mean boolean "
                    f"decadal statistic, no smoothing, full-union mask, {pct_mode} "
                    "percentile, higher_is_worse."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)  "
            f"sen_slope exactly 0 on {sen_zero:.1%} of finite cells (final panel)")
        _CUBE = None

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
