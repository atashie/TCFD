#!/usr/bin/env python3
"""Process ISIMIP2b TEMPERATE NEEDLELEAF EVERGREEN cveg / npp into the TCFD format.

Timber-growth productivity for maritime/temperate conifers (the Sitka-spruce use case).
Two tracks share this script; pick one on the command line:

    python scripts/process_timber_tempnle.py cveg     -> timber_cveg-tempnle_annual
    python scripts/process_timber_tempnle.py npp      -> timber_npp-tempnle_annual

WHY ISIMIP2b AND NOT 3b. There is NO PFT-resolved `cwood` anywhere in ISIMIP2b (all 11
biomes models and the permafrost sector were grepped, 2026-07-28), and no climate-zone-
resolved conifer PFT carries `cwood` in ANY round. ISIMIP3b offers cwood only for
all-climate conifer classes (CLASSIC `evgndltr` 2 GCMs + JULES `ndlevg` 5 GCMs), which pool
boreal, temperate and subtropical conifers under one parameter set -- measured, 34.6% and
29.3% of their global PFT wood sits in the boreal 55-70N band. The user chose PFT
specificity over the wood-only pool, accepting RCP scenarios. Cost of cveg for cwood:
wood is p50 77% (CLASSIC) to 90% (JULES) of conifer cveg, so cveg is a slightly damped,
slightly noisier wood surrogate -- not a different signal.

ENSEMBLE (2026-07-28; see config/isimip_search_catalog.yaml -> temperate_needleleaf_evergreen_2b)

    model           PFT code                              cveg  npp   GCMs  scenarios
    orchidee        tendev                                 y     y     4     rcp26/60/85
    orchidee-dgvm   tendev                                 y     y     2     rcp26/60
    clm45           needleleaf-evergreen-tree-temperate    y     y     ragged (see below)
    lpjml           temperate-needleleaved-evergreen-tree  -     y     4     rcp26/60/85
    ---- EXCLUDED ----
    caraib          ndevteclt                              y     y     4     rcp26/60

Pooled members per scenario: cveg 8/8/6, npp 12/12/10 at rcp26/60/85. rcp85 is thinner by
construction -- orchidee-dgvm and caraib do not publish it.

THREE THINGS THAT MAKE THIS LAYER UNLIKE THE OTHERS, all measured, none assumable:

1. TWO VALUE CONVENTIONS. A PFT field is reported either as a PER-TILE DENSITY (on the
   PFT's own tile area, so sum_i(frac_i * value_i) recovers the cell total) or PER-GRIDCELL
   (already cover-scaled). Verified: orchidee / orchidee-dgvm / clm45 are per-tile;
   lpjml and caraib are per-gridcell. Test that settles it: a per-tile field EXCEEDS the
   all-PFT total (clm45 does so in 94.4% of cells, ratio p50 2.28); a per-gridcell field
   never does and collapses toward 0 where cover -> 0.
   This layer harmonizes onto PER-TILE (stand-level "what would conifers carry here"),
   which is the right framing for a plantation question and the only one that can retain
   clm45 -- clm45 publishes NO cover fraction at all, so it can never be converted the
   other way. lpjml is divided by its own cover with a COVER_FLOOR: its cover is healthy
   in the target regions (p50 27.3% UK/Ireland, 56.6% PNW-BC-SEAK -> only 4x / 2x
   amplification) but globally 61.8% of its non-zero cells sit under 1% cover, where
   dividing would amplify ~528x. Floored cells drop out of n_members rather than carry a
   fabricated value.

2. THE MODELS BARELY SHARE A MASK, and level disagreement is between multi-model and
   single-model CELLS, not between models. On a common mask the models agree well --
   2020s medians span 2.35x (cveg) and 1.83x (npp), right at the csoil precedent -> NO
   NORMALIZATION. But the pooled mean where every model reports is 6.29 vs 0.52 where only
   one does (cveg; 12x) and 2.37e-08 vs 4.61e-09 (npp; 5x), because each model's periphery
   is marginal habitat. Pooling the union would print that step into the maps as hard mask
   edges AND distort the percentile (a periphery cell would rank low for having fewer
   contributing models, not for low growth potential).
   -> MIN_FAMILIES = 2: a cell is emitted only where at least TWO DISTINCT MODEL FAMILIES
   report. Coverage: cveg 17,217 cells (6.6% of land; UK 275, PNW 390), npp 19,485 (7.5%;
   UK 297, PNW 398). n_members / n_models ship per cell so a consumer can still filter to
   the all-model subset. User decision 2026-07-28.

3. ORCHIDEE AND ORCHIDEE-DGVM ARE ONE MODEL, two configurations (with/without dynamic
   vegetation; the non-DGVM file's own global attrs read `model:
   ORCHIDEE-MICT(w/tDGVM)` with a `nodgvm` processing path -- contradictory metadata).
   Counting them separately would give ORCHIDEE 6 of 8 cveg votes at rcp26 and would let
   the pair alone satisfy a naive ">=2 models" rule. So members are aggregated by FAMILY:
   the ensemble value is the mean of family means, and the CI is the spread ACROSS FAMILY
   MEANS. That is a deliberate departure from led/let/burntarea/csoil, which use a flat
   inter-member SD -- documented in ci_definition. User decision 2026-07-28.

Also verified per member across all staged files: units are `kg m-2` (cveg) / `kg m-2 s-1`
(npp) everywhere; every member is annual, 94 years 2006-2099, on a `years since 1661`
365_day axis. Metadata is unreliable and is NOT trusted: orchidee's cveg long_name reads
"crop biomass yield" (an ncrename artifact visible in the file's own history), its npp
reads "positive" yet holds negatives, and caraib carries no long_name at all.
Effective resolution measured on STRICTLY POSITIVE 2x2 blocks (these fields are
zero-inflated, so all-zero blocks otherwise read as false evidence of a coarse grid):
orchidee 3.4%, orchidee-dgvm 2.8%, clm45 0.0%, lpjml 0.3% -> all native 0.5 deg.
caraib 49% (and 35% of 4x4) -> effectively ~1-2 deg, one of its three disqualifiers.

npp is converted from the published kg m-2 s-1 to kg m-2 yr-1 using the members' own
365_day calendar, so the value class is a readable annual carbon flux. npp legitimately
goes NEGATIVE (measured min -3.7e-08 kg m-2 s-1 = net carbon loss), so the npp CI is
NOT floored; the cveg CI is floored at 0 because a carbon stock cannot be negative.

ISIMIP2b runs 2006-2099, so a full 2010s decade exists -- but the layer still begins at
the 2020s baseline (decades 2020s-2090s) to stay aligned with the other TCFD layers.

Output: {track}_{scenario}_processed.nc with {median, percentile, trend, lower_ci,
upper_ci, n_members, n_models, pft_cover} on (decade, lat, lon).
"""

import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402
from utils.trend_significance import (  # noqa: E402
    METHOD as SIGNIFICANCE_METHOD, TREND_METHOD, AnnualEnsembleMean, mk_expanding,
    significance_definition, theilsen_decadal, trend_definition_decadal,
)
from utils.contact_sheet import render_contact_sheet  # noqa: E402

# --- Track configuration ----------------------------------------------------
TRACKS = {
    "cveg": dict(
        layer_id="timber_cveg-tempnle_annual",
        long_name="Vegetation carbon stock, temperate needleleaf evergreen conifers",
        units="kg m-2",
        unit_scale=1.0,
        ci_floor=0.0,          # a carbon STOCK cannot be negative
        statistic="decadal_mean_conifer_vegetation_carbon_kg_per_m2",
    ),
    "npp": dict(
        layer_id="timber_npp-tempnle_annual",
        long_name="Net primary productivity, temperate needleleaf evergreen conifers",
        units="kg m-2 yr-1",
        # Published as kg m-2 s-1 on a 365_day calendar (verified for every member), so the
        # exact factor for these files is 365 * 86400 -- not 365.2425, which would import a
        # calendar the models do not use.
        unit_scale=365.0 * 86400.0,
        ci_floor=None,         # npp is legitimately negative (net carbon loss)
        statistic="decadal_mean_conifer_npp_kg_per_m2_per_year",
    ),
}

# Model families. orchidee + orchidee-dgvm are ONE model in two configurations.
FAMILY = {"orchidee": "orchidee", "orchidee-dgvm": "orchidee",
          "clm45": "clm45", "lpjml": "lpjml", "caraib": "caraib"}

EXCLUDED_MODELS = {
    "caraib": ("effectively ~1-2 deg (49% of non-zero 2x2 blocks and 35% of 4x4 blocks "
               "exactly constant, vs 0.0-3.4% for the retained models); per-gridcell "
               "convention; and no rcp85 for the PFT-resolved fields"),
}
# Models whose PFT values are PER-GRIDCELL and must be divided by cover to reach per-tile.
PER_GRIDCELL = {"lpjml", "caraib"}
COVER_FLOOR = 0.01          # below 1% cover, do not divide -- drop the cell for that member

DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2006, 2099
MIN_FAMILIES = 2            # emit a cell only where >=2 distinct model families report
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = True     # more standing carbon / more growth = better; risk = decline


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract member identity from an ISIMIP2b filename.

    e.g. orchidee_gfdl-esm2m_ewembi_rcp60_2005soc_co2_cveg-tendev_global_annual_2006_2099.nc4
         [0]=model [1]=gcm [2]=bias [3]=scenario [4]=soc [5]=sens [-5]=variable
         [-4]=region [-3]=cadence [-2]=start [-1]=end
    Parsed from the END: model and GCM names contain hyphens but never underscores, while
    the PFT-qualified variable itself contains hyphens too
    (`cveg-needleleaf-evergreen-tree-temperate`), so a fixed forward index would break.
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], soc=p[4], sens=p[5],
                var=p[-5], cadence=p[-3], member=f"{p[0]}_{p[1]}",
                family=FAMILY.get(p[0], p[0]))


def years_from_time(ds):
    """Parse integer calendar years from the CF time axis.

    Every ISIMIP2b member here uses `years since {base}` on a 365_day calendar (verified
    across all staged files; note lpjml writes `years since 1661-01-01` and the others
    `years since 1661-1-1`, so the base year is regex-extracted rather than string-matched).
    `years since` is not valid udunits for cftime, so it is handled arithmetically with a
    FLOOR -- never a round, which would misassign boundary steps. Months and days are
    supported too so the function does not silently mis-handle a future member at another
    cadence.
    """
    t = ds.time
    units = t.attrs.get("units", "")
    cal = t.attrs.get("calendar") or "standard"
    m = re.search(r"(years|months|days|hours|minutes|seconds)\s+since\s+(\d+)", units)
    if not m:
        raise ValueError(f"cannot parse time units {units!r}")
    step, base = m.group(1), int(m.group(2))
    if step in ("years", "months"):
        vals = t.values.astype("float64")
        return np.floor(base + (vals if step == "years" else vals / 12.0)).astype(int)
    import cftime  # only needed for date-based axes; 2b members never reach here
    dates = cftime.num2date(t.values, units, calendar=cal)
    return np.array([d.year for d in np.atleast_1d(dates)], dtype=int)


def load_raw(fpath):
    """Load one file as annual (year, lat, lon) over MIN..MAX_YEAR, unscaled."""
    ds = xr.open_dataset(fpath, decode_times=False)
    var = [k for k in ds.data_vars if ds[k].ndim >= 3][0]
    yrs = years_from_time(ds)
    da = ds[var].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    if da.year.size != np.unique(da.year.values).size:      # sub-annual -> annual mean
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            da = da.groupby("year").mean("year", skipna=True)
    return da.sortby("year")


def cover_fraction(da):
    """Normalize a `pft-*` field to a 0-1 fraction.

    The declared unit CANNOT be trusted: classic and orchidee label these '%' but store
    0-1, while jules, lpjml and caraib store true percent. Decided from the VALUES -- if
    anything exceeds 1.5 the field must be percent (a fraction cannot).
    """
    mx = float(np.nanmax(da.values)) if da.size else 0.0
    return da / 100.0 if mx > 1.5 else da


def decade_mean_map(da, decade):
    """Mean over a decade window -> one map per member per decade."""
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return da.isel(year=sel).mean("year").values.astype(np.float32)


def anchored_trend(stack, i, b_idx):
    """Baseline-anchored trend: (value[i] - value[baseline]) / elapsed decades.

    Each decade's panel is the mean rate of change FROM the 2020s baseline TO that decade,
    so it is exactly the change map divided by elapsed decades -- spatially coherent by
    construction and consistent with the change map (GUARDRAILS S10). Chosen over a
    within-decade Theil-Sen slope because both tracks are noisy year-to-year at the PFT
    level (npp especially), matching the burntarea and csoil layers. The baseline decade
    has no elapsed change -> 0.

    The baseline zero is written ONLY where the baseline median is finite; off-mask cells
    stay NaN. GUARDRAILS S10 requires the trend's finite mask to match the median's, and a
    bare np.zeros() here is a known latent defect in process_burntarea_fire.py and
    process_csoil_soilcarbon.py (WORKFLOW-ISSUES 2026-07-28, finding 6): it makes the whole
    OCEAN a finite zero in the baseline trend, which QA does not catch because it only
    checks that finite baseline trends equal zero -- never that the two masks agree.
    """
    span = i - b_idx
    if span == 0:
        return np.where(np.isfinite(stack[b_idx]), 0.0, np.nan).astype(np.float32)
    return ((stack[i] - stack[b_idx]) / span).astype(np.float32)


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the 2020s baseline, returned as a RISK score.

    Higher standing carbon / higher productivity is BETTER, so the raw rank is INVERTED
    via (101 - raw): a low or declining value gets a HIGH risk percentile.

    Tier is decided from the baseline's exact-zero fraction, not assumed. These PFT fields
    are heavily zero-inflated at source (per-member exact zeros: orchidee 50.0%, lpjml
    ~72%, orchidee-dgvm ~60%, clm45 ~6%), so the two-tier branch is the likely one: a zero
    means the PFT carries no biomass there, which is a real floor rather than a small
    positive value, and ranking zeros inside the continuous distribution would compress
    everything above them.

    For npp, values can be NEGATIVE (net carbon loss). In the two-tier branch they fall in
    the non-positive tier alongside zeros, i.e. the worst productivity, which after
    inversion is the highest risk -- the intended reading.
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
            res = np.ones(vals.shape, np.float32)     # zero / negative -> raw 1
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="left") / n_nz
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


def family_pool(maps_by_member, families_by_member, min_families=MIN_FAMILIES,
                ci_floor=None):
    """Pool member maps by MODEL FAMILY, equal weight per family.

    Returns (value, lower, upper, n_members, n_families). The central value is the mean of
    family means -- not a flat member mean -- so orchidee's two configurations and its 4
    GCMs cannot outvote clm45's 1-2 members (decision 3 in the module docstring). The CI is
    +/- 1 standard deviation ACROSS FAMILY MEANS, i.e. between-model disagreement, which
    is the dominant uncertainty here (2.35x / 1.83x level spread between models).

    Cells with fewer than `min_families` contributing families are set to NaN: pooling them
    would mix multi-model consensus with single-model marginal habitat that sits ~12x lower
    and would print mask edges into the maps.
    """
    fams = sorted(set(families_by_member))
    fam_means, fam_counts = [], []
    for fam in fams:
        idx = [i for i, f in enumerate(families_by_member) if f == fam]
        sub = np.stack([maps_by_member[i] for i in idx], 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fam_means.append(np.nanmean(sub, axis=0))
            fam_counts.append(np.sum(np.isfinite(sub), axis=0).astype(np.float32))
    fam_stack = np.stack(fam_means, 0)                       # (family, lat, lon)
    n_members = np.sum(np.stack(fam_counts, 0), axis=0)
    n_families = np.sum(np.isfinite(fam_stack), axis=0).astype(np.float32)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        value = np.nanmean(fam_stack, axis=0)
        sd = np.nanstd(fam_stack, axis=0)

    thin = n_families < min_families
    value[thin] = np.nan
    sd[thin] = np.nan
    n_members[thin] = 0
    n_families[thin] = 0

    lower = value - sd
    upper = value + sd
    if ci_floor is not None:
        lower = np.clip(lower, ci_floor, None)
        upper = np.clip(upper, ci_floor, None)
    return (value.astype(np.float32), lower.astype(np.float32), upper.astype(np.float32),
            n_members, n_families)


def main():
    track = (sys.argv[1] if len(sys.argv) > 1 else "").strip().lower()
    if track not in TRACKS:
        log(f"usage: {os.path.basename(__file__)} {{{'|'.join(TRACKS)}}}")
        return 2
    cfg = TRACKS[track]
    layer_id = cfg["layer_id"]

    value_files = [str(p) for p in storage.stage_raw(
        layer_id, f"*_{track}-*_global_annual_{MIN_YEAR}_{MAX_YEAR}.nc4")]
    cover_files = [str(p) for p in storage.stage_raw(
        layer_id, f"*_pft-*_global_annual_{MIN_YEAR}_{MAX_YEAR}.nc4")]
    if not value_files:
        log(f"ERROR: no {track} member files in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(layer_id)}")
        log("Ingest members with storage.ingest_raw() first.")
        return 1

    kept = [f for f in value_files if parse_name(f)["model"] not in EXCLUDED_MODELS]
    if len(kept) != len(value_files):
        for mdl, why in sorted(EXCLUDED_MODELS.items()):
            n = sum(1 for f in value_files if parse_name(f)["model"] == mdl)
            if n:
                log(f"EXCLUDING {n} {mdl} files -- {why}")
        log(f"  {len(kept)} of {len(value_files)} staged value files retained")
    value_files = kept
    if not value_files:
        log("ERROR: every staged file was excluded; nothing to process.")
        return 1

    meta = {f: parse_name(f) for f in value_files}
    cover_idx = {}
    for f in cover_files:
        i = parse_name(f)
        cover_idx[(i["model"], i["gcm"], i["scenario"])] = f

    scenarios = sorted({m["scenario"] for m in meta.values()})      # dynamic (GUARDRAILS S3)
    models = sorted({m["model"] for m in meta.values()})
    families = sorted({m["family"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    pfts = sorted({m["var"].split("-", 1)[1] for m in meta.values()})

    log("=" * 74)
    log(f"Processing {track} for TEMPERATE NEEDLELEAF EVERGREEN -> TCFD format")
    log("=" * 74)
    log(f"Members: {len(value_files)} | scenarios: {scenarios}")
    log(f"Models: {models}  -> families: {families}")
    log(f"GCMs: {gcms}")
    log(f"PFT codes pooled: {pfts}")
    log(f"Units: {cfg['units']} (scale x{cfg['unit_scale']:g} from source)")
    log(f"Convention: harmonized to PER-TILE density; per-gridcell models converted by "
        f"cover (floor {COVER_FLOOR:.0%}): {sorted(PER_GRIDCELL & set(models))}")
    log(f"Mask rule: >= {MIN_FAMILIES} distinct model families per cell.")
    log("No normalization (models agree to 2.35x/1.83x on a common mask); no smoothing.")
    log("Direction: higher_is_better (risk = DECLINE) -> percentile inverted.")

    da0 = load_raw(value_files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0

    # ---- Pass 1: per-member decadal maps, harmonized to per-tile density -------
    dec = {s: {} for s in scenarios}
    fam_of = {}
    cover_dec = {s: {} for s in scenarios}
    # The p-value is tested on the ensemble-mean ANNUAL series (GUARDRAILS S15), so the
    # harmonized annual values are accumulated as members stream past. NOTE this is a
    # FLAT member mean, whereas `median` pools family-mean-of-family-means -- a
    # deliberate simplification (user decision 2026-07-30), recorded in
    # significance_pooling. `trend` uses only the decadal medians (S10).
    annual_acc = {s: AnnualEnsembleMean(MIN_YEAR, MAX_YEAR, (LAT, LON))
                  for s in scenarios}
    for f in sorted(value_files):
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_raw(f)
        if da.year.size == 0:
            log(f"  WARNING: {os.path.basename(f)} has no years in "
                f"{MIN_YEAR}-{MAX_YEAR}; skipping (check the time-axis parse)")
            continue

        cov_da = None
        key = (info["model"], info["gcm"], info["scenario"])
        if info["model"] in PER_GRIDCELL:
            if key not in cover_idx:
                log(f"  ERROR: {info['model']} is per-gridcell but no pft- cover file was "
                    f"staged for {key}; refusing to pool it uncorrected. Skipping.")
                continue
            cov_da = cover_fraction(load_raw(cover_idx[key]))
        elif key in cover_idx:
            cov_da = cover_fraction(load_raw(cover_idx[key]))    # diagnostic only

        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        cmaps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        n_floored = 0
        for i, d in enumerate(DECADES):
            raw = decade_mean_map(da, d)
            if raw is None:
                continue
            if cov_da is not None:
                cmap = decade_mean_map(cov_da, d)
                cmaps[i] = cmap
                if info["model"] in PER_GRIDCELL:
                    ok = np.isfinite(cmap) & (cmap >= COVER_FLOOR)
                    n_floored += int(np.sum(np.isfinite(raw) & (raw > 0) & ~ok))
                    raw = np.where(ok, raw / np.maximum(cmap, COVER_FLOOR), np.nan)
            maps[i] = raw * cfg["unit_scale"]

        # Same harmonization as above but per YEAR: per-gridcell models are divided by
        # their own annual cover with the 1% floor before the unit scale (S13).
        ann = da.values.astype(np.float32)
        if cov_da is not None and info["model"] in PER_GRIDCELL:
            cov_y = cov_da.reindex(year=da.year.values).values
            ok_y = np.isfinite(cov_y) & (cov_y >= COVER_FLOOR)
            ann = np.where(ok_y, ann / np.maximum(cov_y, COVER_FLOOR), np.nan)
        annual_acc[s].add(da.year.values, ann * cfg["unit_scale"])

        dec[s][member] = maps
        cover_dec[s][member] = cmaps
        fam_of[member] = info["family"]
        extra = f"  [/cover, {n_floored:,} cells floored out]" if info["model"] in PER_GRIDCELL else ""
        log(f"  loaded {info['model']:<15} {info['gcm']:<13} {s}{extra}")
        del da, cov_da

    if not any(dec[s] for s in scenarios):
        log("ERROR: no members loaded successfully.")
        return 1

    # ---- Shared 2020s baseline (per member, averaged across its scenarios) -----
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec[s]})
    shared_member_maps = []
    for member in all_members:
        per_scen = [dec[s][member][b_idx] for s in scenarios if member in dec[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_member_maps.append(np.nanmean(np.stack(per_scen, 0), axis=0))

    # ---- Per-member contact sheet (GUARDRAILS S11) -----------------------------
    # Rendered while each member's own field still exists separately -- once pooled, a bad
    # member is diluted and invisible. This is the only view that can show a spatial defect.
    sheet_path = None
    stage = storage.staging_dir(layer_id)
    out_dir = stage / "data"
    try:
        sheet_path = render_contact_sheet(
            {m: shared_member_maps[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", layer_id, BASELINE_DECADE, units=cfg["units"],
            note=(f"Ensemble as processed: {len(all_members)} members, "
                  f"{len(families)} model families, harmonized to per-tile density. "
                  f"Excluded and NOT shown: {', '.join(sorted(EXCLUDED_MODELS))}. "
                  f"Per-member panels are PRE-mask: the published layer keeps only cells "
                  f"with >= {MIN_FAMILIES} families."))
        log(f"\n  contact sheet: {sheet_path}  <-- REVIEW THIS before trusting the layer")
    except Exception as e:                                       # noqa: BLE001
        log(f"\n  WARNING: contact sheet failed ({type(e).__name__}: {e}); the per-member "
            f"visual check of GUARDRAILS S11 has NOT been produced")

    # The percentile REFERENCE distribution is fitted once, over every member, so a
    # percentile means the same thing in every scenario.
    shared_value, shared_lo, shared_hi, shared_nmem, shared_nfam = family_pool(
        shared_member_maps, [fam_of[m] for m in all_members], ci_floor=cfg["ci_floor"])

    baseline_flat = shared_value[np.isfinite(shared_value)].ravel()
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)

    # Each scenario's 2020s PANEL is pooled over exactly that scenario's members, using
    # those members' scenario-independent baseline maps.
    #
    # WHY, and this is a correctness fix rather than a preference: ensemble composition
    # VARIES BY SCENARIO here (orchidee-dgvm publishes no rcp85; clm45's GCMs differ per
    # RCP -- rcp26 hadgem2+miroc5, rcp60 ipsl+gfdl, rcp85 gfdl+hadgem2). Pooling an
    # all-member baseline against per-scenario decades therefore differences two different
    # ensembles and manufactures trend. Measured in the first dry run: rcp85 read
    # -0.72 kg m-2/dec at the 2030s and then "recovered", purely because dgvm was in the
    # baseline but absent from rcp85 and sits HIGHER than orchidee on the retained cells.
    # Anchoring each scenario to its own composition removes it. The cost is that the 2020s
    # panel is no longer bit-identical across scenarios (as it is for csoil, whose ensemble
    # was uniform); that is the right trade, since generate_maps derives the change map by
    # differencing the PUBLISHED medians and so must see composition-matched decades.
    idx_of = {m: i for i, m in enumerate(all_members)}
    base_by_scen = {}
    for s in scenarios:
        ms = sorted(dec[s])
        if not ms:
            continue
        base_by_scen[s] = family_pool([shared_member_maps[idx_of[m]] for m in ms],
                                      [fam_of[m] for m in ms], ci_floor=cfg["ci_floor"])

    land = np.isfinite(shared_value)
    log(f"\nShared 2020s baseline: {len(all_members)} members -> {len(families)} families")
    log(f"  cells retained (>= {MIN_FAMILIES} families): {int(land.sum()):,} "
        f"({100*land.sum()/land.size:.2f}% of grid); "
        f"all-family {100*np.mean(shared_nfam[land] == len(families)):.1f}%")
    log(f"  exact-zero fraction={frac_zero:.2%} -> percentile mode={pct_mode}; "
        f"global-mean {track}={np.nanmean(shared_value):.6g} {cfg['units']}")

    # ---- Per-scenario assembly ------------------------------------------------
    for s in scenarios:
        members = sorted(dec[s])
        if not members:
            log(f"\nWARNING: scenario {s} has no members; skipping.")
            continue
        log(f"\n{'='*74}\nAssembling scenario {s} -- {len(members)} members, "
            f"{len(set(fam_of[m] for m in members))} families\n{'='*74}")

        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)
        percentile = np.full_like(median, np.nan)
        nmem = np.zeros_like(median)
        nfam = np.zeros_like(median)
        cover = np.full_like(median, np.nan)

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                # Each member's baseline is its cross-scenario 2020s mean (so it is
                # scenario-independent), but the POOL is this scenario's members, matching
                # the composition of every other decade in this file. See base_by_scen.
                bv, blo, bhi, bnm, bnf = base_by_scen[s]
                median[i], lower[i], upper[i] = bv, blo, bhi
                percentile[i] = pct(bv)
                nmem[i], nfam[i] = bnm, bnf
            else:
                v, lo, hi, nm, nf = family_pool(
                    [dec[s][m][i] for m in members], [fam_of[m] for m in members],
                    ci_floor=cfg["ci_floor"])
                median[i], lower[i], upper[i] = v, lo, hi
                percentile[i] = pct(v)
                nmem[i], nfam[i] = nm, nf
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cover[i] = np.nanmean(
                    np.stack([cover_dec[s][m][i] for m in members], 0), axis=0)

        # Theil-Sen on the DECADAL medians (S10); p-value on the ANNUAL series (S15).
        trend = theilsen_decadal(median, DECADES, window_years=WINDOW_YEARS,
                                 baseline_decade=BASELINE_DECADE).astype(np.float32)
        years_a, mean_annual = annual_acc[s].result()
        tpval, ttau, tnobs = mk_expanding(years_a, mean_annual, DECADES,
                                          window_years=WINDOW_YEARS,
                                          baseline_decade=BASELINE_DECADE)
        for arr in (trend, tpval, ttau):
            arr[~np.isfinite(median)] = np.nan
        tnobs[~np.isfinite(median)] = 0
        for i, d in enumerate(DECADES):
            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(members)} members"
            log(f"  {d}s: {tag:<15} global-mean={np.nanmean(median[i]):.6g}  "
                f"trend={np.nanmean(trend[i]):+.4g} {cfg['units']}/dec  "
                f"cells={int(np.isfinite(median[i]).sum()):,}")

        nmem[nmem == 0] = np.nan        # off-mask: NaN, not a misleading 0
        nfam[~np.isfinite(nmem)] = np.nan
        cover[~np.isfinite(median)] = np.nan

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "trend_pvalue": (["decade", "lat", "lon"], tpval.astype(np.float32)),
                "trend_tau": (["decade", "lat", "lon"], ttau.astype(np.float32)),
                "trend_n_obs": (["decade", "lat", "lon"], tnobs.astype(np.float32)),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
                "n_members": (["decade", "lat", "lon"], nmem),
                "n_models": (["decade", "lat", "lon"], nfam),
                "pft_cover": (["decade", "lat", "lon"], cover),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": f"{track}-tempnle",
                "scenario": s,
                "long_name": cfg["long_name"],
                "units": cfg["units"],
                "statistic": cfg["statistic"],
                "pft_concept": ("temperate needleleaf evergreen conifers. Pooled across each "
                                f"model's own code for that concept: {', '.join(pfts)}. These "
                                "are NOT the same PFT -- each is a distinct model's class with "
                                "its own parameter set -- but they are the same climate-zone "
                                "concept, which is the axis this layer selects on."),
                "value_note": (f"median = mean of MODEL-FAMILY means (not a flat member mean) "
                               f"of the decadal-mean {track} on each model's own PFT tile area. "
                               f"Families: {', '.join(families)}. Only cells with >= "
                               f"{MIN_FAMILIES} contributing families are emitted."),
                "why_isimip2b": ("ISIMIP2b publishes NO PFT-resolved cwood (all 11 biomes "
                                 "models and the permafrost sector checked 2026-07-28), and no "
                                 "climate-zone-resolved conifer PFT carries cwood in ANY round. "
                                 "ISIMIP3b has cwood only for all-climate conifer classes "
                                 "(CLASSIC evgndltr, JULES ndlevg), which pool boreal with "
                                 "temperate conifers -- 34.6% and 29.3% of their global PFT wood "
                                 "sits in the boreal 55-70N band. This layer trades the wood-only "
                                 "pool for climate-zone specificity (user decision 2026-07-28). "
                                 "cveg vs cwood: wood is p50 77-90% of conifer cveg, so cveg is "
                                 "a damped wood surrogate, not a different signal."),
                "value_convention": ("PER-TILE DENSITY (value on the PFT's own tile area, "
                                     "independent of how much of the cell it covers) -- the "
                                     "stand-level 'what would conifers carry here' framing, "
                                     "appropriate to a plantation question. Verified per model: "
                                     "orchidee / orchidee-dgvm / clm45 publish per-tile "
                                     "(a per-tile field EXCEEDS the all-PFT total; clm45 does so "
                                     "in 94.4% of cells, ratio p50 2.28); lpjml and caraib "
                                     "publish PER-GRIDCELL (already cover-scaled) and are "
                                     f"divided by their own cover fraction, floored at "
                                     f"{COVER_FLOOR:.0%}. Per-tile was chosen over per-gridcell "
                                     "because clm45 publishes NO cover fraction at all and could "
                                     "not otherwise be retained. lpjml's cover is healthy where "
                                     "it matters (p50 27.3% UK/Ireland, 56.6% PNW-BC-SEAK -> 4x "
                                     "and 2x amplification) but globally 61.8% of its non-zero "
                                     "cells are under 1% cover, where dividing would amplify "
                                     "~528x; those cells are dropped from n_members instead."),
                "mask_rule": (f">= {MIN_FAMILIES} distinct model families per cell. The models "
                              "barely share a PFT mask, and the pooled level differs far more "
                              "between multi-model and single-model CELLS than between models: "
                              "6.29 vs 0.52 kg m-2 for cveg (12x) and 2.37e-08 vs 4.61e-09 "
                              "kg m-2 s-1 for npp (5x), because each model's periphery is "
                              "marginal habitat. Pooling the union would print that step into "
                              "the maps as hard mask edges and would distort the percentile "
                              "(a periphery cell ranking low for having fewer contributing "
                              "models rather than for low growth potential). Coverage retained: "
                              "cveg 17,217 cells (6.6% of land; UK/Ireland 275, PNW-BC-SEAK 390) "
                              "and npp 19,485 (7.5%; 297 and 398). n_members / n_models are "
                              "emitted so an all-model subset stays recoverable. User decision "
                              "2026-07-28."),
                "model_family_note": ("orchidee and orchidee-dgvm are ONE model in two "
                                      "configurations (with/without dynamic vegetation); the "
                                      "non-DGVM files' own global attrs read "
                                      "'ORCHIDEE-MICT(w/tDGVM)' with a 'nodgvm' processing path, "
                                      "which is contradictory metadata. They are therefore "
                                      "aggregated into ONE family: counting them separately "
                                      "would give orchidee 6 of 8 cveg votes at rcp26 and would "
                                      "let the pair alone satisfy the >=2-model rule. "
                                      "n_models counts FAMILIES, not raw model names."),
                "normalization": ("none. The apparent cross-model spread collapses once the "
                                 "value convention is harmonized and models are compared on a "
                                 "COMMON mask: 2020s medians span 2.35x (cveg) and 1.83x (npp), "
                                 "at or below the csoil precedent of 1.8x that established "
                                 "model democracy in raw units. Measured on own-masks the same "
                                 "numbers read 15.5x and 4.6x, and BEFORE convention "
                                 "harmonization 10.5x and 177x -- those larger figures are "
                                 "artifacts of comparing per-tile against per-gridcell values on "
                                 "different cell populations, not model disagreement."),
                "spatial_smoothing": (f"none. {len(all_members)} members across {len(families)} "
                                      "families is thick enough, and the >=2-family mask already "
                                      "removes the single-model periphery that smoothing would "
                                      "otherwise bleed across."),
                "trend_definition": trend_definition_decadal(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "trend_method": TREND_METHOD,
                "significance_method": SIGNIFICANCE_METHOD,
                "significance_definition": significance_definition(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "significance_pooling": (
                    "FLAT mean across members within each year -- note `median` "
                    "pools family-mean-of-family-means with a >=2-family mask, so "
                    "the significance series is the same members equally weighted, "
                    "not the identical estimator (user decision 2026-07-30). "
                    "Coverage is still taken from median, so no cell gains a "
                    "p-value that median masks out."),
                "trend_units": f"{cfg['units']} decade-1",
                "ci_definition": ("lower/upper_ci = ensemble value +/- 1 standard deviation "
                                  "ACROSS MODEL-FAMILY MEANS, i.e. between-model disagreement, "
                                  "which dominates here (2.35x / 1.83x level spread). This is a "
                                  "DELIBERATE DEPARTURE from led / let / burntarea / csoil, "
                                  "which use a flat inter-member SD: with orchidee supplying 6 of "
                                  "8 cveg members at rcp26, a flat member SD would be dominated "
                                  "by within-family GCM spread and would understate model "
                                  "disagreement. Consequence: on a 2-family cell the SD is over "
                                  "2 values, so the CI is coarse -- use n_models to qualify it."
                                  + ("  Floored at 0 (a carbon stock cannot be negative)."
                                     if cfg["ci_floor"] is not None else
                                     "  NOT floored: npp is legitimately negative (measured min "
                                     "-3.7e-08 kg m-2 s-1 = net carbon loss), so flooring would "
                                     "push lower_ci above the value in loss cells.")),
                "coverage_note": ("n_members = contributing model x GCM simulations; n_models = "
                                  "contributing model FAMILIES (orchidee + orchidee-dgvm count "
                                  "once). NaN off-mask. Cells below the family threshold are "
                                  "masked, not reported with a low count."),
                "pft_cover_note": ("pft_cover = ensemble-mean fractional cover (0-1) of this PFT, "
                                   "a CONFIDENCE/CONTEXT field, not part of the 8 value classes. "
                                   "It answers 'do the models actually place this conifer here', "
                                   "which the per-tile median deliberately does not encode. "
                                   "INCOMPLETE BY CONSTRUCTION: clm45 publishes no cover, so "
                                   "cover averages only the models that do. Low cover with a high "
                                   "median is not a contradiction -- it is the plantation case "
                                   "(models simulate natural vegetation, so UK/Ireland conifer "
                                   "cover is low while the per-tile density is well defined)."),
                "percentile_baseline": (f"{pct_mode}: raw score ranks each cell against the 2020s "
                                       "ensemble land distribution, then INVERTED to a risk "
                                       "percentile (101 - raw) because more standing carbon / "
                                       "higher productivity is better -> low or declining value = "
                                       "high risk percentile."
                                       + ("  Zero and NEGATIVE values share the lowest tier, so "
                                          "net-carbon-loss cells map to the highest risk."
                                          if track == "npp" else "")),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_better",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": ("per_member_shared_across_scenarios, pooled over THIS "
                                    "scenario's members. Each member's 2020s map is the mean "
                                    "of its own 2020s across the scenarios it appears in, so a "
                                    "member's baseline is scenario-independent; but the POOL is "
                                    "this scenario's composition, not all members. Necessary "
                                    "because composition VARIES BY SCENARIO -- orchidee-dgvm "
                                    "publishes no rcp85, and clm45's GCM pair differs per RCP "
                                    "(rcp26 hadgem2+miroc5, rcp60 ipsl+gfdl, rcp85 "
                                    "gfdl+hadgem2). Differencing an all-member baseline against "
                                    "per-scenario decades manufactures trend: it read "
                                    "-0.72 kg m-2/dec at rcp85 2030s with a spurious later "
                                    "'recovery', purely because dgvm was in the baseline, is "
                                    "absent from rcp85, and sits higher than orchidee on the "
                                    "retained cells. Consequence to be aware of: the 2020s panel "
                                    "is NOT bit-identical across scenarios (it is for csoil, "
                                    "whose ensemble was uniform). The percentile REFERENCE "
                                    "distribution is still fitted once over all members, so a "
                                    "percentile means the same thing in every scenario."),
                "decade_note": ("ISIMIP2b runs 2006-2099, so a full 2010s decade DOES exist, but "
                                "the layer still begins at the 2020s baseline (2020s-2090s) to "
                                "stay aligned with the other TCFD layers (user decision "
                                "2026-07-28). The 2090s window is 2090-2099."),
                "window_years": WINDOW_YEARS,
                # Recorded as a per-scenario BREAKDOWN that is byte-identical in every
                # scenario file, rather than as a scalar count that differs between them.
                # The manifest guard in utils/layer_publish.py rightly rejects attributes
                # that disagree across scenarios (it is how a non-uniformly PROCESSED
                # ensemble gets caught); this ensemble is non-uniform in the SOURCE DATA
                # (orchidee-dgvm publishes no rcp85), which is a fact to record, not a
                # reason to loosen the guard.
                # Member IDENTITY, not just a count: rcp26 and rcp60 both have 8 members
                # but they are not the SAME 8 (clm45 contributes hadgem2-es+miroc5 at rcp26
                # and ipsl-cm5a-lr+gfdl-esm2m at rcp60), so a count-based signature would
                # claim their baseline panels must match bit-for-bit when they legitimately
                # cannot. Consumed by the composition grouping in generate_qa_report.py.
                "members_by_scenario": "; ".join(
                    f"{sc}=[{','.join(sorted(dec[sc]))}]"
                    for sc in scenarios if dec[sc]),
                "member_counts_by_scenario": ", ".join(
                    f"{sc}:{len(dec[sc])}m/{len(set(fam_of[m] for m in dec[sc]))}f"
                    for sc in scenarios if dec[sc]),
                "n_members_total": len(all_members),
                "n_model_families": len(families),
                "impact_models": ",".join(sorted({meta[f]["model"] for f in value_files})),
                "model_families": ",".join(families),
                "gcms": ",".join(gcms),
                "scenario_depth": ("rcp85 is thinner by construction: orchidee-dgvm and the "
                                   "excluded caraib do not publish the PFT-resolved fields at "
                                   "rcp85. Pooled members per scenario -- cveg 8/8/6 and "
                                   "npp 12/12/10 at rcp26/60/85."),
                "soc_sens": ("land use held fixed for every member (2005soc). CO2: transient "
                             "'co2' everywhere EXCEPT clm45 at rcp60/gfdl-esm2m, which publishes "
                             "only the fixed-2005-CO2 run and is accepted as a fallback rather "
                             "than losing a GCM. Mixed CO2 is documented precedent (csoil)."),
                "model_notes": ("Verified per member across all staged files 2026-07-28: units "
                                "'kg m-2' (cveg) / 'kg m-2 s-1' (npp) everywhere; every member "
                                "annual, 94 years 2006-2099, 'years since 1661' on a 365_day "
                                "calendar (lpjml writes '1661-01-01', others '1661-1-1', so the "
                                "base year is regex-extracted). METADATA IS NOT TRUSTED: "
                                "orchidee's cveg long_name reads 'crop biomass yield' (an "
                                "ncrename artifact visible in the file's own history attribute), "
                                "its npp long_name claims 'positive' yet the field holds "
                                "negatives, and caraib carries no long_name at all. Source "
                                "exact-zero fractions per member: orchidee 50.0%, lpjml ~72%, "
                                "orchidee-dgvm ~60%, clm45 ~6%. Land masks differ sharply -- "
                                "clm45 writes only ~9.7k cells (its PFT tile extent) vs "
                                "orchidee ~37k -- which is why the mask rule exists."),
                "effective_resolution": ("Measured from the VALUES on STRICTLY POSITIVE 2x2 "
                                         "blocks, not from the declared dims: orchidee 3.4%, "
                                         "orchidee-dgvm 2.8%, clm45 0.0%, lpjml 0.3% of blocks "
                                         "constant -> all native 0.5 deg. The test MUST exclude "
                                         "all-zero blocks: these fields are 50-72% zeros at "
                                         "source, and counting zero blocks inflates the same "
                                         "models to 49-66% and reads as a false coarse-grid "
                                         "signal. caraib measures 49% (and 35% of 4x4 blocks) "
                                         "even after that correction -> genuinely ~1-2 deg."),
                "excluded_sources": ("caraib is EXCLUDED from the pooled ensemble though its "
                                     "ndevteclt files are ingested (user decision 2026-07-28), "
                                     "for three independent reasons: effectively ~1-2 deg "
                                     "resolution (GUARDRAILS S11); the per-gridcell value "
                                     "convention; and no rcp85. Its raw files stay in S3 staging "
                                     "as a sensitivity panel. caraib's class codes carry NO "
                                     "long_name, so they were identified by cover-weighted "
                                     "biogeography: ndevteclt peaks at 40-55N (48.2%) and is the "
                                     "genuine temperate class, while ndevtecdt -- despite 'te' in "
                                     "the code -- peaks at 55-70N (65.0%) and is BOREAL. "
                                     "Also excluded: the 'ndsu*' codes, which are summergreen "
                                     "(deciduous), not evergreen; every non-2005soc land-use "
                                     "variant; and lpjml/orchidee's PERMAFROST-sector copies of "
                                     "the same PFT fields, which would double-weight them."),
                "source_dataset": (f"ISIMIP2b OutputData/biomes ({track}-<temperate needleleaf "
                                   f"evergreen code>, annual 2006-2099, EWEMBI-bias-adjusted "
                                   f"CMIP5)"),
                "description": (f"{cfg['long_name']} processed to the TCFD 8-value-class format "
                                f"with a shared 2020s baseline. {len(families)}-family x "
                                f"CMIP5-GCM ensemble harmonized to per-tile density, no "
                                f"normalization, no spatial smoothing, >= {MIN_FAMILIES}-family "
                                f"mask, higher_is_better (risk = declining growth potential)."),
            },
        )
        path = out_dir / f"{track}-tempnle_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        layer_id,
        stage,
        created_by=f"scripts/process_timber_tempnle.py {track}",
        # A same-day, same-commit rerun publishes alongside its predecessor as -b, -c, ...
        # rather than overwriting it, so a superseded version stays auditable.
        on_exists="bump",
        notes=(f"ISIMIP2b temperate needleleaf evergreen {track}; "
               f"{len(families)} model families ({', '.join(families)}) x CMIP5 GCMs, "
               f"harmonized to per-tile density, >= {MIN_FAMILIES}-family mask. "
               f"caraib EXCLUDED (~1-2 deg, per-gridcell, no rcp85). "
               f"See WORKFLOW-ISSUES.md 2026-07-28; GUARDRAILS S8-S11."),
    )

    log("\nGenerating QA/QC report and maps...")
    finalize_layer(layer_id, version=version, extra_maps=[sheet_path])
    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
