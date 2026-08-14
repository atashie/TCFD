"""Process `cropfailure` (ISIMIP3b/SSP crop-failure exposure) into the TCFD output contract.

`cropfailure` (Zantout2025, ISIMIP3b DerivedOutputData) is the SSP re-issue of the Lange
2020 crop-failure-exposure concept (`lec`), renamed by hazard word rather than `le*` code.
Each cell carries a per-year share of the cell's cropland exposed to an unprecedented
crop-failure event.

This is a SEPARATE LAYER from `lec` (ISIMIP2b, rcp26/rcp60, GEPIC + PEPIC), not a
replacement: 8 crop models vs 2, CMIP6 vs CMIP5, SSPs vs RCPs.

Ensemble: 8 crop models x 5 CMIP6 GCMs x 3 SSPs = 120 files, 700 MiB, enumerated
2026-08-13 with zero gaps. 40 members per scenario, composition IDENTICAL across
ssp126/ssp370/ssp585, so the shared 2020s baseline is valid. soc/sens uniformly
2015soc/default -- no harmonization compromise, no model dropped.

MEASURED data nature (GUARDRAILS 9 -- from the values, never the name), all 120 members:

  * CONTINUOUS fraction in [0, 1], NOT binary: 82,962-130,168 distinct values per member.
    This differs from its 2b relatives' behaviour and was not inherited from any of them.
  * `long_name` "cropfailure area share", `units` "1".
  * 96.14% of finite cell-years are exact zeros -- BUT that figure is an artifact of
    zero-fill (see the mask section). Within the cropland footprint the annual exact-zero
    share is 46.37%-80.63%, median 60.89%.

THE PUBLISHER ZERO-FILLS THE ENTIRE GLOBE. This is the single most important fact about
the raw files and it is invisible to a finite-mask survey:

  * Every member is non-NaN over ~100% of the 259,200-cell grid (only 236 land cells are
    NaN anywhere). Ocean, Antarctica, Greenland and the Sahara all read exact 0, not NaN.
  * So `np.isfinite(...)` carries NO footprint here, and using it as the land mask would
    put 87% ocean zeros into the percentile baseline population, report n_members=40 over
    open ocean, and hand every ocean cell the lowest-risk percentile tier.
  * The footprint therefore comes from where the field is ever NON-ZERO -- the cropland
    each model actually simulates. Union across all 8 models: 39,890 cells.

  This is NOT the `floodedarea` ocean leak, and the distinction was measured rather than
  assumed. `floodedarea` carries real values over open ocean; here, 39,406 of the 39,890
  footprint cells (98.8%) fall inside the official ISIMIP3b land-sea mask. Of the 484 that
  do not, 63% are directly adjacent to land (coastline disagreement at 0.5 degrees) and the
  remainder are small islands the generic mask does not resolve as land -- Lofoten
  (68.25N 13.25E), Shetland (60.25N -1.75E). Their values are marginal (median 0.0016).
  They come from exactly two models, `ldndc` (476) and `lpjml` (453); the other six have
  none. The model-derived footprint is published as-is: the crop models are the authority
  on where they simulate crops, and dropping real island cropland to satisfy a coarser
  generic mask would be the worse error. `n_models` is emitted per cell either way.

  EACH MEMBER IS MASKED TO ITS OWN CROPLAND FOOTPRINT before pooling, which is the second
  consequence of the zero-fill and is easy to miss because nothing errors without it. A
  model that grows nothing in a cell is not reporting "0% of the cropland here failed"; it
  has no opinion about that cell. Left unmasked, its structural zero is pooled as an
  observation, and `n_members` / `n_models` -- which are CONTRACT variables that downstream
  delivery filters on -- would read 40 and 8 in every published cell while this layer's own
  `mask_rule` tells the reader to filter on them. Measured: only 38.1% of footprint cells
  have all 40 members simulating crops (1-10 members: 563 cells, 11-30: 5,471, 31-40:
  33,856). Effect on the value is 1.03x globally (0.00535 -> 0.00550) but 4.4x on the
  1-10-member cells and 1.5x on 11-30 -- negligible for headline numbers, decisive for a
  site-level query in marginal cropland, which is exactly where such a query is fragile.

FIVE DECISIONS SPECIFIC TO THIS LAYER. None is inherited; each was measured here.

1. DECADAL STATISTIC = `pooled_mean_zero_inflated` (mean +/- 1 SD on a CONTINUOUS field).

   This is the OUTPUT-SPEC third branch and a DECLARED deviation, taken on measurement and
   confirmed by the user 2026-08-13. Measured on the 2020s panel over the footprint, using
   the SAME pooling the published statistic uses (scenarios concatenated along the member
   axis, as every shipped layer does -- not averaged first):

       branch          exposed cells   exact-zero share
       pooled_median          1,351          96.61%
       pooled_mean           39,872           0.05%

   The median branch erases 96.6% of exposed cropland -- Iowa, the Pampas, Beauce, the
   Sahel and Java all read exactly 0 at baseline under it, so their decadal change would be
   measured from a floor of zero. `cropfailure` estimates an EXPECTED ANNUAL EXPOSED SHARE,
   and a median of near-Bernoulli draws does not estimate that. Same justification as
   `let`, re-measured; `burntarea` at 29.2% did not qualify and took the median branch.

   The pooling convention matters to this number and is worth stating once: averaging each
   member across scenarios BEFORE taking the median gives a much milder 80.08% / 7,948
   exposed. Both say the same thing, but the figure recorded in
   `decadal_statistic_rationale` is recomputed at run time under the convention actually
   published, so the counterfactual stays comparable to the statistic it justifies rather
   than to a differently-pooled cousin. Do not hand-edit it to match this docstring; if the
   two disagree, the attribute is right and this text is stale.

2. NO MINIMUM-MODEL MASK -- the full 39,890-cell union is published (`--min-models 1`).

   The coverage tiers show a steep monotone level gradient on the 2020s panel:

       1 model     46 cells   0.00001        5 models    820 cells   0.00018
       2 models   469 cells   0.00007        6 models  1,453 cells   0.00030
       3 models   136 cells   0.00016        7 models  6,201 cells   0.00166
       4 models   591 cells   0.00007        8 models 30,174 cells   0.00670

   That is the OPPOSITE direction from `led`'s artifact, and the difference decides the
   rule. `led` masks at >= 2 models because its solo cells read 1.63x the all-model level:
   one high model undiluted, i.e. a composition artifact to remove. Here low-coverage cells
   are LOW, not high -- a cell only one of eight models grows crops in is marginal cropland
   with near-zero failure exposure. There is no inflation to correct, and a mask would
   delete real (if marginal) cropland. Re-measured, not inherited from `led` OR `driedarea`.

   CAVEAT recorded rather than masked: a 1-model cell still carries NO inter-model
   uncertainty in its CI. Filter on `n_models` if that matters.

3. NO SPATIAL SMOOTHING, established by a SPLIT-HALF TEST rather than by draw count.

   Roughness (mean |cell - 4-neighbour mean| / mean) is 0.347 on the 2020s panel, close
   enough to `let`'s raw 0.389 that the draw count alone (40 members x 10 yr = 400, vs
   `let`'s 40) would not have settled it. Splitting the ensemble into two disjoint
   20-member halves:

       roughness   full 0.347   half A 0.351   half B 0.359
       half A vs half B spatial correlation: Pearson 0.977, Spearman 0.991 (n=39,872)

   Halving the ensemble barely moves roughness and the halves agree almost exactly, so the
   roughness is REAL SPATIAL STRUCTURE, not sampling noise -- unlike `let`, where a thin
   ensemble produced one-cell-wide storm tracks. It is concentrated at the cropland
   boundary: roughness is 0.790 on footprint EDGE cells against 0.304 in the interior.
   Smoothing would smear cropland exposure into neighbouring non-cropland.

4. PERCENTILE: single_tier, `higher_is_worse`. The two-tier mode exists for baselines that
   are materially zero-inflated (>2% exact zeros); under the mean branch this baseline is
   0.05% exact zeros, so tiering is unnecessary and would misrepresent a continuous field.
   Ranked against the shared 2020s footprint distribution -- note the reference population
   is CROPLAND, not all land, so a percentile here is a rank among cropland.

5. NO NORMALIZATION. Inter-model spread on the 2020s footprint mean is 6.69x
   (promet 0.00162, lpjml 0.00220, ldndc 0.00226, pepic 0.00399, epic-iiasa 0.00566,
   crover 0.00679, isam 0.00942, cygma1p74 0.01081). All 8 report the same dimensionless
   share, so this is structural model uncertainty and belongs in the CI, not a unit
   mismatch to rescale away.

TIME DECODING: the `time` coordinate carries `long_name`/`standard_name`/`axis` and NO
`units` attribute at all, so xarray cannot decode it (`.dt` raises AttributeError) and
there is no epoch to read. The values are a bare contiguous integer sequence, 165..250 for
a file declaring 2015_2100 -- `years since 1850`, undeclared. Rather than hardcode 1850,
years are taken from the filename span (which IS part of the ISIMIP grammar) and the axis
is asserted to match it in length and to step by exactly 1. Contrast `driedarea`, whose
axis is a properly declared `days since 1601-01-01`; the handling is per-publication.

GUARDRAILS 12 reference-site check (`scripts/check_cropfailure_nature.py`), mean over
2015-2100 across all 120 members -- non-zero and non-NaN at every site:

    Iowa, US Corn Belt      0.1303      Ukraine steppe wheat    0.0256
    Mato Grosso, BR soy     0.0958      Beauce, FR wheat        0.0257
    Pampas, AR              0.0355      WA wheatbelt, AU        0.0727
    Punjab, IN wheat/rice   0.1490      Guinea savanna, NG      0.0275
    North China Plain       0.1042      Ethiopian highlands     0.0426
    Sahel millet, ML        0.0117      Java rice, ID           0.0088

HOW THIS LAYER MUST BE READ, and the one thing a reader will get wrong. The value is
DEPARTURE FROM A FIXED HISTORICAL REFERENCE, not absolute agricultural risk, and that
difference REVERSES the ranking an untrained reader expects. A reliably productive region
with a narrow historical range crosses its own unprecedentedness threshold more readily
than a region that has always been marginal and variable. Measured here: Iowa sits at the
99.3rd percentile of world cropland, the Sahel at 69.4.

That is the measure working as defined, not a defect -- the intensity of change relative to
local norms genuinely is greater in Iowa, and a farming system built around a stable
climate has further to fall. It is not a claim that Iowan agriculture is more fragile than
Sahelian agriculture in absolute terms, and the layer says nothing about food security,
yield levels, or capacity to absorb a bad year.

Because a reader can be wrong about this while every number in front of them is correct,
the caveat is MACHINE-ENFORCED rather than left to prose: `cropfailure-3b` declares
`relative_baseline: true` in config/layer_registry.yaml, `generate_delivery_caveats.py`
promotes that to a MUST-DISCLOSE caveat, and both reports refuse to render without it. The
customer-facing wording lives in the registry's `relative_baseline_note` so it can be
revised without reprocessing 700 MiB. NOTE: files built before 2026-08-13 carry an earlier,
substantively identical `interpretation_caveat` that does not quote the Iowa/Sahel figures;
the registry note is the authoritative customer-facing text either way.

Expect `sen_slope` to collapse toward exactly 0 across much of the domain: the field is
zero-inflated, so most year-pairs are 0->0. READ `ols_slope` ON THIS LAYER.

Usage:
    python scripts/process_cropfailure_isimip3b.py [--scenarios ssp126] [--skip-slopes]
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

VAR = "cropfailure"
OUT_VAR = "cropfailure"
LAYER_ID = "cropfailure-isimip3b_cropfailure_annual"
#: ISIMIP3b starts in 2015, so the baseline decade IS decade index 0 -- 2015-2019 cannot
#: form a decade. Same as `driedarea`; unlike the ISIMIP2b layers, which carry a 2010s panel.
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

#: The decadal value is an exposed-area SHARE, bounded by [0, 1].
CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False          # more crop failure = worse -> percentile NOT inverted
SLOPE_PER_DECADE = 10.0

#: See decision 2 -- measured here, NOT inherited. 1 == publish the full footprint.
MIN_MODELS = 1

#: See decision 1 -- OUTPUT-SPEC third branch, declared deviation, user-confirmed.
STAT_NAME = "pooled_mean_zero_inflated"
CENTRAL = "mean"

#: Eight distinct crop models, no two configurations of one code, so family == model.
MODEL_FAMILY = {}

SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

#: Set in the parent before forking the slope pool; workers inherit it copy-on-write.
_CUBE = None


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) from a Zantout2025 filename.

    e.g. zantout2025_crover_gfdl-esm4_w5e5_ssp126_2015soc_default_cropfailure_global_annual_2015_2100.nc

    NOTE these filenames DO carry a leading publication token, where the 3b sibling
    Heinicke2026 (`driedarea`) does not. DerivedOutputData grammar is per-PUBLICATION, so
    the fields are read from the END, which is index-stable under either convention.
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
    """Years per record, WITHOUT trusting the time axis to declare itself.

    Zantout2025's `time` carries no `units` attribute, so there is no epoch to decode and
    xarray leaves it as bare floats (165..250). The span comes from the filename and the
    axis is checked against it; a member whose axis disagrees raises rather than silently
    binning into the wrong decade.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    y0, y1 = int(p[-2]), int(p[-1])
    t = np.asarray(ds["time"].values, dtype="float64")
    n = y1 - y0 + 1
    if t.size != n:
        raise ValueError(f"{os.path.basename(fpath)}: {t.size} records but the filename "
                         f"declares {y0}-{y1} ({n} years)")
    steps = np.unique(np.diff(t))
    if steps.size != 1 or not np.isclose(steps[0], 1.0):
        raise ValueError(f"{os.path.basename(fpath)}: time axis is not a contiguous "
                         f"annual sequence (steps {steps})")
    return np.arange(y0, y1 + 1, dtype=int)


def read_values(fpath):
    """Raw (year, lat, lon) values with the declared fill replaced by NaN."""
    with xr.open_dataset(fpath, decode_times=False) as ds:
        da = ds[VAR]
        yrs = decode_years(ds, fpath)
        vals = da.values.astype("float32")
        fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    if fill is not None and np.isfinite(fill):
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan
    return yrs, vals


def load_member(fpath):
    """Load one member as (years, (year, lat, lon)) restricted to the analysis window.

    Deliberately does NOT derive or assert a land mask. This product zero-fills the whole
    globe, so `isfinite` is ~everything and carries no information; the footprint is built
    separately from where the field is ever non-zero. It also does not assert a
    time-invariant mask: 90 of 120 members have a handful of NaN cells that move between
    years (236 land cells total), which is noise in the fill, not a moving footprint.
    """
    yrs, vals = read_values(fpath)
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    return yrs[keep], vals[keep]


def member_footprint(fpaths):
    """Cells where this member ever records non-zero exposure, over its FULL record.

    A member that never fails in a cell across 86 years is not saying "0% of the cropland
    here failed" -- it is saying it has no cropland here. Because the publisher zero-fills,
    those two are identical within any single cell-year and can only be separated across
    the whole record. Even very low-risk cropland resolves: the Pampas, Beauce and Sahel
    reference sites all carry small non-zeros.

    Taken as a UNION across the member's three scenario files. The cropland mask is static
    (uniform 2015soc), so it should not vary by scenario, and a union is the inclusive
    reading if it ever does.
    """
    acc = None
    for f in fpaths:
        _, v = read_values(f)
        nz = (np.nan_to_num(v) > 0).any(axis=0)
        acc = nz if acc is None else (acc | nz)
    return acc


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline CROPLAND distribution.

    Two-tier when the baseline is materially zero-inflated (>2% exact zeros). Under the
    mean branch this baseline measures 0.05% exact zeros, so the single-tier path is taken
    -- but the mode is CHOSEN FROM THE MEASURED SHARE at run time, never hardcoded.
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
    """Put a footprint-cell vector back on the (lat, lon) grid; everything else NaN."""
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
    """Worker: expanding slopes for one contiguous block of footprint cells."""
    s, e, years, decade, baseline, window, chunk, max_pairs = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk, max_pairs=max_pairs)
    # Return plain arrays, NOT the SlopeResult: it is a dict subclass whose
    # `__getattr__ = dict.__getitem__` turns pickle's probe for `__getstate__` into a
    # KeyError instead of an AttributeError, which kills the whole pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, max_pairs, jobs, n_land):
    """expanding_slopes over all footprint cells, fanned across `jobs` forked workers."""
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
                    help="publish a cell only where this many crop models simulate crops; "
                         "1 = full footprint (this layer's measured default, see docstring)")
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
    log(f"cropfailure (ISIMIP3b/SSP crop-failure exposure) -> TCFD contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")

    with xr.open_dataset(files[0], decode_times=False) as ds0:
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

    # ---- Pass 1: per-MODEL CROPLAND footprint -------------------------------- #
    # NOT the finite mask. This product zero-fills the globe, so `isfinite` is ~100%
    # everywhere and says nothing; what carries the footprint is where a model ever
    # records non-zero exposure, i.e. the cropland it actually simulates.
    t0 = time.time()
    all_members = sorted({meta[f]["member"] for f in files})
    files_by_member = {m: [f for f in files if meta[f]["member"] == m]
                       for m in all_members}
    mem_mask = {m: member_footprint(files_by_member[m]) for m in all_members}

    model_mask = {m: np.zeros((LAT, LON), bool) for m in models}
    for mem, mk in mem_mask.items():
        model_mask[mem.rsplit("_", 1)[0]] |= mk

    nmod_static = np.sum([model_mask[m] for m in models], axis=0).astype(np.int16)
    union = nmod_static > 0
    keep2d = nmod_static >= args.min_models
    log(f"\nCropland footprint (per-model non-zero coverage, {time.time() - t0:.0f}s):")
    for k in range(1, len(models) + 1):
        c = int((nmod_static == k).sum())
        if c:
            log(f"  exactly {k} model(s): {c:>7,}")
    log(f"  union {int(union.sum()):,} -> publishing >= {args.min_models} model(s): "
        f"{int(keep2d.sum()):,} cells; dropped {int((union & ~keep2d).sum()):,}")

    try:
        from utils.land_mask import get_isimip_landmask
        land = np.nan_to_num(
            xr.open_dataset(get_isimip_landmask("3b"))["mask"].values) > 0
        off = int((union & ~land).sum())
        log(f"  vs ISIMIP3b land-sea mask: {int((union & land).sum()):,} on land, "
            f"{off:,} off ({100.0 * off / max(int(union.sum()), 1):.1f}%) -- coastline "
            f"disagreement and unresolved small islands, NOT an ocean leak; see docstring")
    except Exception as e:  # noqa: BLE001
        log(f"  (land-mask cross-check unavailable: {type(e).__name__}: {e})")

    land_idx = np.flatnonzero(keep2d.ravel())
    if args.limit_cells and args.limit_cells < land_idx.size:
        land_idx = land_idx[np.linspace(0, land_idx.size - 1,
                                        args.limit_cells).astype(int)]
    n_land = land_idx.size

    # ---- Pass 2: pack (member, year, footprint_cell) per scenario ------------ #
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    # Each member is masked to ITS OWN cropland footprint. Without this, the global
    # zero-fill makes every member finite in every cell, so a model that grows nothing in
    # a cell contributes a structural 0 that (a) drags the pooled mean down and (b) makes
    # n_members read 40 and n_models 8 everywhere -- flatly contradicting this layer's own
    # mask_rule, which tells a reader to filter on n_models. Measured effect: global mean
    # 0.00535 -> 0.00550 (1.03x), but 4.4x on cells with 1-10 simulating members and 1.5x
    # on 11-30. Only 38.1% of footprint cells have all 40 members simulating crops.
    mem_keep = {m: mem_mask[m].ravel()[land_idx] for m in all_members}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        yrs, cube = load_member(f)
        flat = cube.reshape(cube.shape[0], -1)
        drop = ~mem_keep[m]
        for k, y in enumerate(yrs):
            yi = y_index.get(int(y))
            if yi is not None:
                row = flat[k, land_idx]
                row[drop] = np.nan
                annual[s][slot[s][m], yi] = row
        del cube, flat
    nmem_static = np.sum([mem_keep[m] for m in all_members], axis=0)
    log(f"Packed {len(files)} members over {n_land:,} footprint cells "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident) "
        f"[{time.time() - t0:.0f}s]")
    log("  per-member footprints applied (structural zeros -> NaN); simulating-member "
        "count per cell:")
    for lo, hi in ((1, 10), (11, 20), (21, 30), (31, 39), (40, 40)):
        c = int(((nmem_static >= lo) & (nmem_static <= hi)).sum())
        if c:
            log(f"    {lo:>2}-{hi:<2} members: {c:>7,}")

    # ---- Field nature: measured, never assumed (GUARDRAILS 9) ---------------- #
    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> is_boolean_field={boolean}")
    if boolean:
        log("  ERROR: `cropfailure` measured CONTINUOUS at ingest 2026-08-13 "
            "(82,962-130,168 distinct values per member). A boolean read means the inputs "
            "changed -- re-run scripts/check_cropfailure_nature.py before processing.")
        return 3
    fin0 = annual[scenarios[0]][np.isfinite(annual[scenarios[0]])]
    log(f"  annual cell-year exact-zero fraction WITHIN the footprint: "
        f"{float((fin0 == 0).mean()):.2%}")
    log(f"  decadal statistic: {STAT_NAME} (central={CENTRAL}) -- see docstring decision 1")
    del fin0

    # ---- Shared 2020s baseline ----------------------------------------------- #
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

    # The declared deviation must be justified by NUMBERS, recomputed at run time rather
    # than quoted from the docstring: what would the median branch have published here?
    med_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    alt_med, _, _ = pooled_decadal_stat(
        med_pool, years[bwin], BASELINE_DECADE, window_years=WINDOW_YEARS,
        central="median")
    del med_pool
    alt_zero = float(np.mean(alt_med[np.isfinite(alt_med)] == 0.0))
    this_zero = float(np.mean(b_med[np.isfinite(b_med)] == 0.0))
    alt_exposed = int(np.sum(alt_med > 0))
    this_exposed = int(np.sum(b_med > 0))
    del alt_med
    rationale = (
        f"OUTPUT-SPEC third branch, DECLARED deviation (user decision 2026-08-13). "
        f"Measured on this run's {BASELINE_DECADE}s panel over the cropland footprint: "
        f"the pooled-MEDIAN branch leaves {alt_zero:.2%} of cells at exactly zero "
        f"({alt_exposed:,} exposed), the pooled-MEAN branch {this_zero:.2%} "
        f"({this_exposed:,} exposed) -- the median erases "
        f"{100.0 * (1 - alt_exposed / max(this_exposed, 1)):.1f}% of exposed cropland. "
        f"`cropfailure` estimates an EXPECTED ANNUAL EXPOSED SHARE and a median of "
        f"near-Bernoulli draws does not estimate that. Same justification as `let`, "
        f"re-measured here; `burntarea` at 29.2% annual zeros did NOT qualify.")

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"\nShared {BASELINE_DECADE}s baseline: footprint n={baseline_flat.size:,}, "
        f"exact-zero {frac_zero:.2%}, percentile mode={pct_mode}, "
        f"mean share {np.nanmean(b_med):.5f}, max {np.nanmax(b_med):.4f}")
    log(f"  median-branch counterfactual: {alt_zero:.2%} exact-zero "
        f"({alt_exposed:,} exposed) vs mean-branch {this_exposed:,}")

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
            "member_field": (f"mean exposed cropland share over the {BASELINE_DECADE}s "
                             "baseline decade, pooled across scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Inter-model "
                     "spread is 6.69x (promet 0.00162 to cygma1p74 0.01081 on the 2020s "
                     "footprint mean), so members are expected to differ in level; look "
                     "for spatial-pattern outliers and footprint differences, not level "
                     "agreement."),
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
                    cube, years, d, window_years=WINDOW_YEARS, central=CENTRAL)
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
                log(f"  {d}s: {tag:<15} mean={np.nanmean(out['median'][i]):.5f}  "
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
        fin_ols = out["ols_slope"][-1][np.isfinite(out["ols_slope"][-1])]
        sen_zero = float((fin_sen == 0).mean()) if fin_sen.size else float("nan")
        active = (fin_sen != 0) | (fin_ols != 0)
        sen_zero_active = (float((fin_sen[active] == 0).mean())
                           if active.any() else float("nan"))

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Cropland area share exposed to unprecedented crop failure",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": STAT_NAME,
                "decadal_statistic_rationale": rationale,
                "field_nature": "continuous",
                "value_note": (
                    "median = MEAN over the pooled (year x member) sample inside the decade "
                    f"window, across {len(mem)} members (8 crop models x 5 CMIP6 GCMs), "
                    "i.e. the EXPECTED ANNUAL SHARE of the cell's cropland exposed to an "
                    "unprecedented crop-failure event. The raw ISIMIP field is a CONTINUOUS "
                    "fraction in [0,1] (measured: 82,962-130,168 distinct values per "
                    "member), NOT the binary flag some members of this hazard family use. "
                    "Published long_name is 'cropfailure area share', units '1'."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clamped to [{CI_FLOOR}, {CI_CAP}] since the "
                    "value is a share. Carries interannual, inter-GCM AND inter-model "
                    "spread; inter-model is 6.69x here (promet 0.00162, lpjml 0.00220, "
                    "ldndc 0.00226, pepic 0.00399, epic-iiasa 0.00566, crover 0.00679, "
                    "isam 0.00942, cygma1p74 0.01081 on the 2020s footprint mean)."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "fitted over an EXPANDING window from the start of the 2020s baseline "
                    "through the end of the target decade, stacking every (year, member) "
                    "observation as an independent point. The 2020s panel is NaN (no "
                    "elapsed period). The estimators fail in OPPOSITE regimes, so "
                    "disagreement means a cell's trend is not robust. This field is "
                    f"zero-inflated, so sen_slope is exactly 0 on {sen_zero:.1%} of finite "
                    f"cells in the final panel ({sen_zero_active:.1%} of ACTIVE cells, "
                    "i.e. those where either slope is non-zero -- the active figure is the "
                    "honest one). READ ols_slope ON THIS LAYER."),
                "slope_units": "1 decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal exposed share ranked against the "
                    f"shared {BASELINE_DECADE}s ensemble CROPLAND distribution. The "
                    "reference population is the cropland footprint, NOT all land, so a "
                    "percentile here is a rank among cropland. Not inverted -- more crop "
                    "failure is worse."),
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
                    "none -- all 8 models report the same dimensionless exposed-area share. "
                    "The 6.69x inter-model level spread is structural uncertainty, not a "
                    "unit mismatch, and belongs in the CI."),
                "spatial_smoothing": (
                    "none. Roughness on the 2020s panel is 0.347, close to `let`'s raw "
                    "0.389, so the draw count alone (40 members x 10 yr = 400) did not "
                    "settle it. A SPLIT-HALF TEST did: two disjoint 20-member halves give "
                    "roughness 0.351 and 0.359 and correlate at Pearson 0.977 / Spearman "
                    "0.991 (n=39,872), so the roughness is REAL SPATIAL STRUCTURE, not "
                    "sampling noise. It concentrates at the cropland boundary -- 0.790 on "
                    "footprint edge cells vs 0.304 in the interior -- so smoothing would "
                    "smear cropland exposure into neighbouring non-cropland."),
                "minimum_models": args.min_models,
                "mask_rule": (
                    "THE PUBLISHER ZERO-FILLS THE WHOLE GLOBE: every member is non-NaN over "
                    "~100% of the 259,200-cell grid, with ocean, Antarctica, Greenland and "
                    "the Sahara all reading exact 0 rather than NaN. `isfinite` therefore "
                    "carries NO footprint and was NOT used as the mask; using it would put "
                    "87% ocean zeros into the percentile baseline and report n_members=40 "
                    "over open ocean. The footprint is where the field is ever NON-ZERO -- "
                    f"the cropland each model actually simulates -- union {int(union.sum()):,} "
                    f"cells, publishing >= {args.min_models} model(s) = {int(keep2d.sum()):,}. "
                    "EACH MEMBER IS ALSO MASKED TO ITS OWN FOOTPRINT before pooling, so a "
                    "model that grows nothing in a cell contributes no observation there "
                    "rather than a structural zero. Without that, n_members/n_models would "
                    "read 40/8 in EVERY cell and this caveat would be unusable; measured, "
                    "only 38.1% of footprint cells have all 40 members simulating crops "
                    "(1-10 members: 563 cells, 11-30: 5,471, 31-40: 33,856). Effect on the "
                    "value is 1.03x globally (0.00535 -> 0.00550) but 4.4x on the "
                    "1-10-member cells and 1.5x on 11-30. "
                    "MEASURED coverage tiers on the 2020s panel show a monotone level "
                    "gradient (1 model 0.00001 at n=46, rising to 8 models 0.00670 at "
                    "n=30,174), which is the OPPOSITE direction from `led`'s artifact: low-"
                    "coverage cells here are marginal cropland with near-zero exposure, not "
                    "one high model undiluted. There is no composition artifact to remove, "
                    "so no minimum-model cut is applied. The rule was re-measured, NOT "
                    "inherited from `led` (>= 2) or `driedarea` (full union). Caveat: a "
                    "1-model cell carries no inter-model uncertainty in its CI -- filter on "
                    "n_models if that matters."),
                "footprint_vs_land_mask": (
                    "39,406 of 39,890 footprint cells (98.8%) fall inside the official "
                    "ISIMIP3b land-sea mask. The 484 that do not are NOT an ocean leak: 63% "
                    "are directly adjacent to land (coastline disagreement at 0.5 degrees) "
                    "and the rest are small islands the generic mask does not resolve -- "
                    "Lofoten (68.25N 13.25E), Shetland (60.25N -1.75E). They come from two "
                    "models only, ldndc (476) and lpjml (453), with marginal values (median "
                    "0.0016). Contrast `floodedarea`, which carries REAL values over open "
                    "ocean. The model-derived footprint is published as-is: the crop models "
                    "are the authority on where they simulate crops."),
                "interpretation_caveat": (
                    "THIS IS DEPARTURE FROM A FIXED HISTORICAL REFERENCE, NOT ABSOLUTE "
                    "AGRICULTURAL RISK, AND THE DIFFERENCE REVERSES THE RANKING A READER "
                    "EXPECTS. The exposure flag is an unprecedentedness threshold against a "
                    "fixed reference distribution, so a reliably productive region with a "
                    "narrow historical range crosses it more readily than a region that has "
                    "always been marginal and variable. MEASURED on this layer: Iowa sits "
                    "at the 99.3rd percentile of world cropland and the Sahel at 69.4. That "
                    "is the measure working as defined -- the intensity of change relative "
                    "to local norms IS greater in Iowa, and a farming system built around a "
                    "stable climate has further to fall -- and it is NOT a claim that Iowan "
                    "agriculture is more fragile than Sahelian agriculture in absolute "
                    "terms. The layer says nothing about food security, yield levels, or "
                    "capacity to absorb a bad year. Read a high score as 'this place is "
                    "moving furthest from what it is built for' and pair it with local "
                    "baseline productivity and adaptive capacity before drawing any "
                    "conclusion about consequences. Do not read it as a yield map."),
                "crop_composition_undocumented": (
                    "THE INDEX CARRIES NO CROP TOKEN, and the sidecar `specifiers` block "
                    "names only the variable `cropfailure`, so WHICH crops are aggregated "
                    "into it, and how they are weighted, is not readable from the archive. This is material for site-level use: a "
                    "crop-aggregated index over crops a site does not grow is not that "
                    "site's risk. Resolve against the Zantout2025 publication before "
                    "quoting this layer for a specific crop or a specific farm."),
                "ensemble_note": (
                    "8 crop models x 5 CMIP6 GCMs = 40 members per scenario, COMPLETE "
                    "(8 x 5 x 3 = 120 files enumerated 2026-08-13, no missing combination). "
                    "Composition is IDENTICAL across ssp126/ssp370/ssp585, so the shared "
                    "2020s baseline is valid. This is the deepest ensemble of any shipped "
                    "layer in this product."),
                "source_dataset": (
                    "ISIMIP3b DerivedOutputData/Zantout2025 -- the SSP re-issue of the Lange "
                    "2020 crop-failure-exposure concept, named by hazard word rather than "
                    "le* code. SEPARATE LAYER from the ISIMIP2b `lec` product (GEPIC + "
                    "PEPIC x 4 CMIP5 GCMs, rcp26/rcp60), not a replacement. Every input "
                    "file is verified against the publisher's sha512, published at "
                    "{stem}.json -- NOT {stem}.nc.json, which 404s and reads as 'no "
                    "sidecars'. The ingest made exactly that error and ran on "
                    "Content-Length alone until it was caught and corrected 2026-08-14."),
                "reference_sites": (
                    "GUARDRAILS 12, mean over 2015-2100 across all 120 members, all "
                    "non-zero and non-NaN: Iowa 0.1303, Punjab 0.1490, North China Plain "
                    "0.1042, Mato Grosso 0.0958, WA wheatbelt 0.0727, Ethiopian highlands "
                    "0.0426, Pampas 0.0355, Guinea savanna 0.0275, Beauce 0.0257, Ukraine "
                    "0.0256, Sahel Mali 0.0117, Java 0.0088."),
                "description": (
                    "ISIMIP3b/SSP crop-failure exposure processed to the TCFD output "
                    "contract (OUTPUT-SPEC.md) with a shared 2020s baseline; "
                    f"{len(mem)}-member 8-model x 5-GCM CMIP6 ensemble, "
                    f"{STAT_NAME} decadal statistic, no smoothing, full cropland-footprint "
                    f"mask, {pct_mode} percentile, higher_is_worse."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / 1024**2:.1f} MB)  "
            f"sen_slope exactly 0 on {sen_zero:.1%} of finite / {sen_zero_active:.1%} of "
            f"active cells (final panel)")
        _CUBE = None

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
