"""Process `heatwave` (ISIMIP3b/SSP heatwave exposure) into the TCFD output contract.

`heatwave` (Zantout2025, ISIMIP3b DerivedOutputData) is the SSP re-issue of the Lange 2020
heatwave-exposure concept (`leh`), renamed by hazard word rather than `le*` code. Each cell
carries a per-year flag for whether it saw a heatwave.

IT IS NOT A NEWER VERSION OF `leh`. THE INDEX ITSELF CHANGED, and that is the first thing to
know about this layer:

    ISIMIP2b  Lange2020 / `HWMId-humidex` -> HWMId (relative) AND Humidex >= 45 (absolute)
    ISIMIP3b  Zantout2025 / `HWMID-NONE`  -> HWMId ONLY; the "NONE" is the humidity term

Lange et al. 2020 (Earth's Future, doi:10.1029/2020EF001616) pairs a RELATIVE indicator --
the Heat Wave Magnitude Index daily, which depends only on daily maximum temperature -- with
an ABSOLUTE one, the Humidex, so that the events counted "would also adversely affect human
health". Zantout et al. 2025 (Nat. Commun., doi:10.1038/s41467-025-65600-7, Methods) drops
the absolute criterion: a cell is exposed "if the HWMId of that year exceeds the 97.5th
percentile of the picontrol distribution of that grid cell".

So THIS LAYER MEASURES UNPRECEDENTED HEAT RELATIVE TO EACH CELL'S OWN PREINDUSTRIAL CONTROL,
NOT DANGEROUS-TO-HUMANS HEAT. It carries no humidity term at all, so it cannot evidence a
wet-bulb, heat-stress or workforce-safety claim, and it is expected to flag MORE area than
`leh` in arid and high-latitude cells where the Humidex gate would never have opened.

MEASURED data nature (GUARDRAILS 9 -- from the values, never the name), all 15 members,
`scripts/check_heatwave_nature.py` 2026-08-13:

  * STRICTLY BINARY {0,1}: exactly 2 unique values per member, range [0,1] throughout.
    `is_boolean_field` -> True, so the decadal statistic is the pooled MEAN = heatwave
    FREQUENCY. The published `long_name` is "heatwave area share" and `units` is "1", which
    reads CONTINUOUS and is a MISNOMER -- the field is an occurrence flag, so its decadal
    mean is a frequency, not an area. Same trap as `floodedarea` in the sibling publication.
  * Annual exact-zero share 19.6%-69.1% (ensemble median 45.3%) -- nowhere near `let`'s
    97.84%, so the `pooled_mean_zero_inflated` third branch is NOT in scope here.
  * Mask: 67,420 finite cells, IDENTICAL in all 15 members and time-invariant. 26.0% of the
    grid, i.e. a genuine land footprint -- NOT the `floodedarea` failure mode (94.7% of the
    globe including open ocean). 1,623 cells fall off the ISIMIP3b land-sea mask (945
    coastal fringe, 678 isolated, half of those tropical islands the 0.5deg mask does not
    resolve); the shipped `driedarea` carries 1,236 off the same mask and `led` 46, so this
    is a mask-CONVENTION difference, not a defect. The 27,092 ISIMIP3b land cells the product
    never covers are 100% Antarctica.
  * 5 members per scenario, composition IDENTICAL across ssp126/ssp370/ssp585 (complete
    1 x 5 x 3 = 15 matrix, no gaps), so the shared 2020s baseline is valid. soc/sens
    uniformly 2015soc/default.

TIME DECODING: the axis is `days since 2015-01-01` on `proleptic_gregorian`, stamped
mid-year with ALTERNATING 365 / 365.5 day steps -- the publisher is tracking leap years, so
cftime decoding through xarray is correct and yields exactly 86 unique years 2015-2100 with
no duplicates and no gaps (verified 2026-08-13). Note this differs from the Heinicke2026
`driedarea` axis (`days since 1601-01-01`) and from the ISIMIP2b Lange2020 `years since`
axis: the time convention is per-publication. Dividing days by 365 would drift.

THE DEFINING PROPERTY OF THIS LAYER: THE FLAG SATURATES.

Because exposure is "above THIS cell's own preindustrial 97.5th percentile", warming pushes
cells permanently over the threshold and the flag pins at 1. Measured on the ingested files,
pooled (year x member) decadal exposure frequency and the share of finite cells at exactly
1.0:

    scenario   2020s    2050s    2090s     cells at 1.0 in the 2090s
    ssp126     0.316    0.473    0.462            1.2%
    ssp370     0.310    0.607    0.846           32.6%
    ssp585     0.316    0.650    0.903           45.9%

Two consequences, BOTH OF WHICH INVERT THE USUAL READING, and neither of which any other
shipped layer has:

  1. THE TREND IS CENSORED AT THE CEILING. A cell exposed in every year of the decade in
     every member has slope ~0 -- and that means "permanently and maximally exposed", not
     "no trend". Both estimators fail the same way here, so the dual-slope disagreement rule
     does NOT rescue it: agreement between ols_slope and sen_slope at ~0 is exactly what a
     saturated cell produces.
  2. THE PERCENTILE LOSES DISCRIMINATION AT THE TOP. Cells at 1.0 all rank 100, so the
     ssp585 2090s panel carries one tie block holding nearly half the published cells.

NO NEW STATISTIC IS INVENTED HERE TO PAPER OVER THAT. The contract's fields are emitted as
specified and the censoring is DECLARED instead, in `saturation_caveat`, with the measured
per-panel shares in `saturation_by_decade`. A saturated cell-decade is identifiable from the
published fields alone -- `median == 1.0` with a zero-width CI (`lower_ci == upper_ci == 1`,
because the pooled sample has no variance) -- so a downstream consumer can mask or flag it
without a private convention. How the tie block should be RANKED, and whether a
time-to-saturation measure should be published alongside, are open decisions and are
deliberately left open rather than answered here.

TWO DECISIONS SPECIFIC TO THIS LAYER. Neither is inherited; both were measured here first.

1. NO MINIMUM-MODEL MASK, AND THE RULE IS VOID RATHER THAN RELAXED. There is exactly ONE
   impact model (hwmid-none), so `n_models` is 1 at every published cell and the `led`
   layer's ">= 2 impact models" rule cannot be expressed at all -- it is not that the
   threshold was lowered, it is that the quantity it thresholds does not vary. Nor is there
   a GCM-count fallback to apply: the finite mask is IDENTICAL across all five GCM members
   (67,420 cells each), so every published cell sits in 5 of 5 members and no coverage tier
   exists to cut on. `n_members` and `n_models` are still emitted per cell so the degeneracy
   is visible in the file rather than only in this docstring.

   CAVEAT recorded rather than masked: with one impact model the CI carries interannual and
   inter-GCM spread ONLY. It contains NO inter-model structural uncertainty, so it is
   narrower than a multi-model layer's CI for reasons that have nothing to do with
   confidence. Do not compare CI widths across layers.

2. NO SPATIAL SMOOTHING. 5 members x 10 years = up to 50 Bernoulli draws per cell-decade,
   and this is a broad, spatially coherent field driven by a large-scale temperature signal
   -- not the one-cell-wide storm tracks that forced an aggressive kernel onto the 4-member
   `let` layer. Confirmed by eye on the per-member contact sheet before this was decided:
   all five members show smooth continental-scale structure with no block structure, seams,
   banding, hemisphere flip or mask defect, and the same geography member to member. GCM
   level spread is 1.58x (gfdl-esm4 0.2583 -> ukesm1-0-ll 0.4073 global land mean).

WHAT LOOKING AT THE CONTACT SHEET ALSO SHOWED, and what a value-check table could never
have: the field is BRIGHT ACROSS THE TROPICS AND DARK AT HIGH LATITUDES in the 2020s
baseline, in every member. Measured, |lat| <= 23.5 reads 0.561 against 0.149 for |lat| > 50
-- 3.77x -- before any warming signal has accumulated. This is the relative threshold being
crossed sooner where interannual variance is LOW, and the tropics have the lowest
temperature variance on the planet. It is the same relative-baseline artifact the drought
layers carry, running the other way, and it is recorded in `latitudinal_artifact` because a
reader will otherwise take it for a statement about heat severity.

Percentile: ranked against the shared 2020s land distribution, `higher_is_worse`. Tier mode
is chosen from the MEASURED zero-fraction of that baseline, not assumed.

Unlike `driedarea`, do NOT expect `sen_slope` to be zero everywhere for zero-inflation
reasons -- this field is far less zero-inflated. Expect it to be zero for the OPPOSITE
reason in the high-forcing panels: saturation. Read both, and read `saturation_by_decade`
before reading either.

GUARDRAILS 12 reference sites (2015-2100 mean -> 2090s mean, all 15 members, none NaN,
none identically zero): Kuwait 0.957->0.982, Singapore 0.889->0.952, Sao Paulo 0.826->0.897,
Cordoba 0.688->0.818, Chicago 0.688->0.836, Jacobabad 0.646->0.855, Lagos 0.650->0.739,
Phoenix 0.609->0.800, Delhi 0.548->0.776, Paris 0.406->0.636, Frankfurt 0.345->0.576,
Yakutsk 0.271->0.509, Sydney 0.226->0.461. Chicago reading ABOVE Delhi is the relative
baseline working as designed, not a defect.

Usage:
    python scripts/process_heatwave_isimip3b.py [--scenarios ssp126] [--skip-slopes]
                                                [--limit-cells N] [--members-only]
                                                [--jobs N]
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

VAR = "heatwave"
OUT_VAR = "heatwave"
LAYER_ID = "heatwave-isimip3b_heatwave_annual"
#: ISIMIP3b starts in 2015, so the baseline decade IS decade index 0 -- unlike the ISIMIP2b
#: `led`/`let` layers, which carry a full 2010s panel first. 2015-2019 cannot form a decade.
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

#: The decadal value is a FREQUENCY (fraction of model-years with a heatwave), in [0, 1].
CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False          # more heatwave = worse -> percentile NOT inverted
SLOPE_PER_DECADE = 10.0

#: See decision 1. One impact model exists, so this rule is VOID here, not relaxed.
MIN_MODELS = 1

#: Measured boolean -> the contract's `pooled_mean_boolean` branch (mean +/- 1 SD).
STAT_NAME = "pooled_mean_boolean"

#: A single impact model; family == model.
MODEL_FAMILY = {}

#: A cell-decade is SATURATED when the pooled (year x member) sample is degenerate at 1 --
#: every observation flagged. Reported per panel; see the saturation discussion above.
SATURATION_EPS = 1e-6

SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

#: Set in the parent before forking the slope pool; workers inherit it copy-on-write.
_CUBE = None


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) from a Zantout2025 filename.

    e.g. zantout2025_hwmid-none_gfdl-esm4_w5e5_ssp126_2015soc_default_heatwave_global_
         annual_2015_2100.nc

    Zantout2025 filenames DO carry a leading publication token, where the 3b sibling
    Heinicke2026 does not, so a forward index reads the wrong column and would report the
    bias-adjustment (`w5e5`) as the scenario. Read from the END, which is index-stable
    under either convention.
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

    Time is decoded with cftime through xarray rather than by days-per-year arithmetic. The
    axis here is `days since 2015-01-01` (proleptic_gregorian) stamped mid-year with
    alternating 365/365.5 steps -- the publisher IS tracking leap years, so decoding yields
    exactly one record per year 2015-2100. Dividing by 365 would drift into the wrong bin.
    """
    ds = xr.open_dataset(fpath)                    # decode_times=True
    da = ds[VAR]
    yrs = ds.time.dt.year.values.astype(int)
    vals = da.values.astype("float32")
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    ds.close()

    if len(set(yrs.tolist())) != yrs.size:
        raise ValueError(f"{os.path.basename(fpath)}: decoded time axis has duplicate "
                         "years -- the calendar convention changed, re-check before binning")

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

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros). On this layer
    the baseline is NOT zero-inflated -- it is the TOP of the range that degenerates in the
    later high-forcing panels, which a two-tier rule does not address and this function does
    not pretend to. Mode is chosen from the measured share, not assumed.
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
            res = np.ones(vals.shape, np.float32)      # never a heatwave -> raw 1
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


def saturated_fraction(med_flat):
    """Share of finite cells whose decadal frequency is pinned at 1.

    The whole point of publishing this number: at 1.0 the pooled sample has no variance, so
    both slopes go to ~0 and the percentile ties at 100. That reads as "stable and extreme"
    when it means "maximally exposed every year of the decade".
    """
    fin = med_flat[np.isfinite(med_flat)]
    if fin.size == 0:
        return float("nan")
    return float((fin >= 1.0 - SATURATION_EPS).mean())


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
                         "VOID on this layer -- there is one impact model (see docstring)")
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
    log(f"heatwave (ISIMIP3b/SSP heatwave exposure) -> TCFD contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")
    log("INDEX NOTE: HWMID-NONE = HWMId with NO humidity term. This is unprecedented-heat "
        "exposure relative to each cell's own picontrol, not a heat-stress index.")

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

    # ---- Pass 1: per-MEMBER land masks -> publication mask ------------------- #
    # Per-MEMBER, not per-model: with one impact model a per-model survey would collapse
    # to a single row and hide whether the five GCMs actually agree on coverage. They do
    # (measured identical), and that is worth asserting rather than assuming.
    t0 = time.time()
    model_mask = {m: np.zeros((LAT, LON), bool) for m in models}
    member_cells = {}
    for f in files:
        with xr.open_dataset(f) as ds:
            da = ds[VAR]
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
            v = da.values.astype("float32")
        if fill is not None:
            v = np.where(np.isclose(v, np.float32(fill), rtol=1e-6), np.nan, v)
        m2d = np.isfinite(v).any(axis=0)
        model_mask[meta[f]["model"]] |= m2d
        member_cells.setdefault(meta[f]["member"], set()).add(int(m2d.sum()))
        del v

    nmod_static = np.sum([model_mask[m] for m in models], axis=0).astype(np.int16)
    union = nmod_static > 0
    keep2d = nmod_static >= args.min_models
    cell_counts = sorted({c for cs in member_cells.values() for c in cs})
    log(f"\nLand mask ({time.time() - t0:.0f}s): {len(models)} impact model(s), "
        f"{len(member_cells)} members")
    log(f"  finite cells per member: {cell_counts}"
        + ("  <-- IDENTICAL across all members" if len(cell_counts) == 1 else
           "  <-- members DIFFER; n_members carries it"))
    log(f"  union {int(union.sum()):,} -> publishing {int(keep2d.sum()):,} cells "
        f"(the >=N-models rule is void: N_models is 1 everywhere)")

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
        log(f"  loaded {info['model']:<12} {info['gcm']:<15} {s}  "
            f"{int(mask2d.sum()):,} land cells  [{time.time() - t0:.0f}s]")
    log(f"Packed {len(files)} members over {n_land:,} cells "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident)")

    # ---- Field nature: measured, never assumed (GUARDRAILS 9) ---------------- #
    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> is_boolean_field={boolean}")
    if not boolean:
        log("  ERROR: `heatwave` measured BINARY at ingest 2026-08-13. A continuous read "
            "means the inputs changed -- re-run check_heatwave_nature.py before processing.")
        return 3
    fin0 = annual[scenarios[0]][np.isfinite(annual[scenarios[0]])]
    log(f"  annual cell-year exact-zero fraction: {float((fin0 == 0).mean()):.2%}  "
        f"exact-one: {float((fin0 == 1).mean()):.2%}")
    log(f"  decadal statistic: {STAT_NAME} -- pooled MEAN = heatwave frequency")
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
        f"mean frequency {np.nanmean(b_med):.4f}, max {np.nanmax(b_med):.4f}, "
        f"saturated {saturated_fraction(b_med):.2%}")

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
            "member_field": (f"heatwave frequency over the {BASELINE_DECADE}s baseline "
                             "decade (mean of the binary annual flag), pooled across "
                             "scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. All five "
                     "members are the SAME impact model (hwmid-none) driven by five "
                     "different CMIP6 GCMs, so differences here are GCM spread, not model "
                     "structural spread. Look for spatial-pattern outliers, block structure "
                     "and mask differences."),
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
        sat_by_decade = {}

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

            sat_by_decade[d] = saturated_fraction(med)

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
                    f"sat={sat_by_decade[d]:6.2%}  {slope_txt}  [{time.time() - td:.0f}s]")

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
        fin_pct = out["percentile"][-1][np.isfinite(out["percentile"][-1])]
        tie_top = float((fin_pct >= 99.5).mean()) if fin_pct.size else float("nan")
        sat_txt = "; ".join(f"{d}s {v:.1%}" for d, v in sat_by_decade.items())

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Land area fraction exposed to heatwave",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": STAT_NAME,
                "field_nature": "boolean",
                "value_note": (
                    "median = MEAN over the pooled (year x member) sample inside the decade "
                    f"window, across {len(mem)} members (1 index model x 5 CMIP6 GCMs), i.e. "
                    "the HEATWAVE FREQUENCY: the fraction of model-years the cell was "
                    "flagged as exposed. The raw ISIMIP field is a BINARY {0,1} annual flag "
                    "(measured: exactly 2 unique values in all 15 members). Its published "
                    "long_name 'heatwave area share' and units '1' read CONTINUOUS and are a "
                    "MISNOMER -- the nature was taken from the values, not the metadata."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clamped to [{CI_FLOOR}, {CI_CAP}] since the "
                    "value is a frequency. WITH ONE IMPACT MODEL THIS CI CARRIES "
                    "INTERANNUAL AND INTER-GCM SPREAD ONLY -- there is no inter-model "
                    "structural uncertainty in it. It is therefore narrower than a "
                    "multi-model layer's CI for a reason unrelated to confidence; do not "
                    "compare CI widths across layers. A cell at frequency 1.0 has a "
                    "ZERO-WIDTH CI because the pooled sample is degenerate, not because the "
                    "estimate is certain -- see saturation_caveat."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "fitted over an EXPANDING window from the start of the 2020s baseline "
                    "through the end of the target decade, stacking every (year, member) "
                    "observation as an independent point. The 2020s panel is NaN (no "
                    "elapsed period). The estimators normally fail in OPPOSITE regimes so "
                    "disagreement flags a fragile trend -- BUT ON THIS LAYER THAT RULE HAS "
                    "AN EXCEPTION: a saturated cell drives BOTH to ~0 and they AGREE at "
                    "zero. Agreement near zero here is therefore ambiguous between 'no "
                    "trend' and 'pinned at the ceiling'; disambiguate with median == 1.0. "
                    f"sen_slope is exactly 0 on {sen_zero:.1%} of finite cells in the final "
                    "panel."),
                "slope_units": "1 decade-1",
                "saturation_caveat": (
                    "THE FLAG SATURATES, AND THIS IS THE LAYER'S DEFINING LIMITATION. "
                    "Exposure is defined as the annual HWMId exceeding the 97.5th percentile "
                    "of THIS cell's own preindustrial control distribution, so under warming "
                    "cells cross the threshold permanently and the frequency pins at 1.0. "
                    "At 1.0 the pooled sample has no variance: the CI collapses to zero "
                    "width, BOTH slopes go to ~0, and the percentile ties at 100. A "
                    "saturated cell-decade therefore reads as 'extreme and no longer "
                    "changing' when it means 'exposed in every year of the decade in every "
                    "member, with no headroom left to measure'. IDENTIFY SATURATED CELLS AS "
                    "median == 1.0 (equivalently a zero-width CI at 1.0) and treat their "
                    "slopes and percentile ranks as censored, not as measurements. How the "
                    "tied top block should be ranked, and whether a time-to-saturation "
                    "measure should be published alongside, are OPEN DECISIONS deliberately "
                    "not answered in this build. "
                    "WORKED EXAMPLE, measured on this ssp585 file -- the censoring INVERTS "
                    "the trend ranking between regions. Amazon (10S-0, 70-55W) goes "
                    "median 0.601 -> 1.000 and 0% -> 100% saturated, and its ols_slope FALLS "
                    "from +0.160 dec-1 in the 2030s to +0.046 in the 2090s. Siberia "
                    "(60-70N, 80-120E) never saturates (0% throughout, median 0.119 -> "
                    "0.750) and its ols_slope RISES from +0.069 to +0.098. So on the 2090s "
                    "panel Siberia out-trends the Amazon 2.1x -- which reads as 'the Amazon "
                    "has stabilised, Siberia is the hotspot' and means 'the Amazon is "
                    "exposed every year in every member and has no headroom left'. Sahara "
                    "(18-28N, 0-25E) shows the same shape, +0.085 -> +0.032 as saturation "
                    "goes 5% -> 100%. NEVER RANK REGIONS BY SLOPE ON THIS LAYER WITHOUT "
                    "READING median AND THE SATURATION SHARE ALONGSIDE."),
                "saturation_by_decade": sat_txt,
                "percentile_tie_at_top_final_panel": round(tie_top, 5),
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal frequency ranked against the shared "
                    f"{BASELINE_DECADE}s ensemble land distribution. Not inverted -- more "
                    "heatwave exposure is worse. The two-tier zero rule does NOT apply here "
                    f"(baseline exact-zero share is {frac_zero:.2%}); the degeneracy on this "
                    "layer is at the TOP of the range, which no tier rule addresses -- see "
                    "saturation_caveat and percentile_tie_at_top_final_panel."),
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
                    "flag from the same index code; only the driving GCM differs."),
                "spatial_smoothing": (
                    "none -- 5 members x 10 years gives up to 50 Bernoulli draws per "
                    "cell-decade, and this is a broad, spatially coherent field driven by a "
                    "large-scale temperature signal, not the one-cell-wide storm tracks that "
                    "forced a kernel onto `let`. Confirmed on the per-member contact sheet "
                    "before the decision was made."),
                "minimum_models": args.min_models,
                "mask_rule": (
                    f"Publishing the full footprint: {int(keep2d.sum()):,} cells. THE "
                    "MINIMUM-MODEL RULE IS VOID ON THIS LAYER, not relaxed -- there is one "
                    "impact model (hwmid-none), so n_models is 1 at every cell and the `led` "
                    "layer's '>= 2 impact models' rule cannot be expressed. There is no "
                    "GCM-count fallback either: the finite mask is IDENTICAL across all five "
                    "GCM members (67,420 cells each, time-invariant), so every published "
                    "cell sits in 5 of 5 members and no coverage tier exists to cut on. "
                    "1,623 of those cells fall off the ISIMIP3b land-sea mask (945 coastal "
                    "fringe, 678 isolated, largely tropical islands the 0.5deg mask does not "
                    "resolve); the shipped `driedarea` carries 1,236 off the same mask and "
                    "`led` 46, so this is a mask-CONVENTION difference, not an ocean leak "
                    "(the field covers 26.0% of the grid; `floodedarea`'s real leak was "
                    "94.7%). Antarctica is absent from the product entirely."),
                "interpretation_caveat": (
                    "THIS IS UNPRECEDENTED HEAT RELATIVE TO A PER-CELL PREINDUSTRIAL "
                    "CONTROL, NOT DANGEROUS-TO-HUMANS HEAT, AND NOT ABSOLUTE HEAT. A cell is "
                    "exposed when its annual HWMId exceeds the 97.5th percentile of its OWN "
                    "picontrol distribution, so a permanently hot site whose regime is "
                    "stable scores LOWER than a temperate site warming past its own "
                    "historical spread. Measured on this ensemble: Chicago 0.688 reads ABOVE "
                    "Delhi 0.548. THE INDEX CARRIES NO HUMIDITY TERM -- the ISIMIP2b "
                    "predecessor `leh` additionally required Humidex >= 45 so that counted "
                    "events 'would also adversely affect human health', and Zantout et al. "
                    "2025 drops that criterion. This layer therefore CANNOT evidence a "
                    "wet-bulb, heat-stress, workforce-safety or equipment-derating claim, "
                    "and it is expected to flag MORE area than `leh` in arid and "
                    "high-latitude cells where the Humidex gate would never have opened."),
                "latitudinal_artifact": (
                    "THE TROPICS READ ~3.8x THE HIGH LATITUDES IN THE BASELINE ITSELF, AND "
                    "THAT IS THE INDEX, NOT THE CLIMATE. Measured on the shared 2020s panel: "
                    "|lat| <= 23.5 mean frequency 0.561 against 0.149 for |lat| > 50, a "
                    "3.77x ratio, with the band profile running 60-90N 0.156, 30-60N 0.228, "
                    "10-30N 0.559, 10S-10N 0.557, 10-30S 0.442, 30-60S 0.153. A relative "
                    "threshold is crossed sooner where INTERANNUAL VARIANCE IS LOW, and the "
                    "tropics have the lowest temperature variance on the planet, so a small "
                    "absolute warming clears the 97.5th picontrol percentile there while a "
                    "much larger one does not in Siberia or Canada. Visible directly on the "
                    "per-member contact sheet as a bright tropical band in all five GCMs. "
                    "Do NOT present this gradient as tropical heat being ~4x more severe "
                    "than boreal heat; it is a statement about departure from local normal."),
                "ensemble_note": (
                    "1 index model (hwmid-none) x 5 CMIP6 GCMs = 5 members per scenario, "
                    "COMPLETE (1 x 5 x 3 = 15 files enumerated 2026-08-13, no missing "
                    "combination). All five members share one index implementation, so the "
                    "spread here is GCM spread only and understates total uncertainty: "
                    "there is no second heatwave-index formulation to disagree with it."),
                "source_dataset": (
                    "ISIMIP3b DerivedOutputData/Zantout2025 (doi:10.48364/ISIMIP.920810) -- "
                    "the SSP re-issue of the Lange 2020 heatwave-exposure concept, named by "
                    "hazard word rather than le* code. SEPARATE LAYER from the ISIMIP2b "
                    "`leh` layer, and NOT a drop-in newer version: different GCM generation, "
                    "different scenarios, and a DIFFERENT INDEX (see interpretation_caveat). "
                    "All 15 files sha512-verified against the publisher sidecars at ingest."),
                "description": (
                    "ISIMIP3b/SSP heatwave exposure frequency processed to the TCFD output "
                    "contract (OUTPUT-SPEC.md) with a shared 2020s baseline; "
                    f"{len(mem)}-member 1-index-model x 5-GCM CMIP6 ensemble, pooled-mean "
                    f"boolean decadal statistic, no smoothing, full-footprint mask, "
                    f"{pct_mode} percentile, higher_is_worse. Saturates at high forcing -- "
                    "read saturation_caveat before using the slopes or the percentile."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)")
        log(f"    saturated cells by decade: {sat_txt}")
        log(f"    final panel: sen_slope exactly 0 on {sen_zero:.1%} of finite cells, "
            f"percentile >= 99.5 on {tie_top:.1%}")
        _CUBE = None

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
