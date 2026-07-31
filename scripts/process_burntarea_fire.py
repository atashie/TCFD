"""Process burntarea (wildfire burnt-area fraction) into the TCFD 8-value-class format.

burntarea = the fraction of each grid cell burned per reporting interval, reported by
ISIMIP vegetation/fire models in PERCENT [0, 100]. This is the direct biophysical fire
signal (NOT the Lange 2020 "exposure to unprecedented wildfire" family member lew, and
not the `ffire` carbon-emissions flux).

SOURCE: ISIMIP3b `biomes`, variable `burntarea-total`, ssp126/ssp370/ssp585.
Moved from ISIMIP2b/RCP to ISIMIP3b/SSP on 2026-07-28: newer round and scenario family
are preferred wherever the newer data is viable. The 2b ensemble was thicker
(4 models x 4 GCMs) but built on CMIP5 forcing; see WORKFLOW-ISSUES.md 2026-07-28.

Ensemble: 12 members per scenario from 3 models x their own GCM sets --

    model             cadence  GCMs  soc/sens          effective grid
    mc2-usfs-r87g5c1  annual     5   nat/default       0.5 deg (clean)
    visit             monthly    5   2015soc/default   0.5 deg (clean)
    classic           monthly    2   2015soc/default   1.0 deg (2x2 blocks)

CLASSIC's soc variant is pinned to `2015soc`, NOT the `2015soc-from-histsoc` that the
csoil layer uses for the same model. That is deliberate and must not be "corrected" to
match: for burntarea the `2015soc-from-histsoc` set is MIXED-SCALE. Its gfdl-esm4 files
are written ~100x low (a 0-1 fraction: monthly mean 0.0032, max 0.32) while its
ukesm1-0-ll files are percent (monthly mean 0.31, max 35.4) -- both declaring
units='%' and long_name='Burnt Area Fraction'. The `2015soc` set is uniformly percent
across both GCMs and all three SSPs (3.58-3.74 %/yr). Verified 2026-07-28; a first
build on the mixed set was caught only because one GCM's max was 1.18% against its
sibling's 114%. **The soc/sens variant must be value-checked per variable -- a variant
that is sound for one variable can be broken for another in the same model.**

CLASSIC is retained despite running effectively at 1.0 deg -- 100.0% of its 2x2 cell
blocks are exactly constant (offset (0,1)), so it was replicated onto the 0.5 deg grid
from a 1 deg native grid. A 1 deg member still carries real information (user decision
2026-07-28). elm-eca, which also publishes burntarea-total monthly for 5 GCMs, is
EXCLUDED: it is effectively ~4x5 deg, far too coarse to be useful here. See
EXCLUDED_MODELS and GUARDRAILS S11.

MONTHLY MEMBERS ARE ANNUALIZED BY **SUM**, NOT MEAN. burntarea accumulates over its
reporting interval: verified against CLASSIC, which publishes the same run at both daily
and monthly cadence -- each published monthly value equals the SUM of that month's daily
values (agreement to 1e-6, correlation 1.00000000 in all 12 months). So the annual burnt
fraction is the sum of the 12 monthly fractions. Taking the mean instead -- as the csoil
layer correctly does for a *stock* -- would under-scale this layer by 12x. See
GUARDRAILS S9 and WORKFLOW-ISSUES.md 2026-07-28.

Value-checked 2026-07-28 (area-weighted global burnt area, first 20 future years):

    model      units  annual %/yr (area-wtd)  Mkm2/yr   long_name
    mc2-usfs     %          4.24               5.63     "Burnt Area Fraction"
    visit        %          3.45               4.56     "Burnt Area Fraction"
    classic      %          4.88               7.26     "Burnt Area Fraction"

For context GFED4 observes ~348 Mha/yr = ~3.5 Mkm2/yr burned globally, so all three
members are the right order of magnitude -- unlike the 2b models lpj-guess (1.09) and
lpjml (0.94), which run ~3.5x low, or 2b clm45/orchidee (~0.00), which declare '%' while
sitting on a 0-1 FRACTION scale. Trust the values, never the declared unit.

All three retained members share the SAME unit (% burnt area) AND agree to within 1.6x
on the global mean, so -- unlike the water-index TWS -- NO normalization is applied: the
models are equal-weighted in raw %, and the inter-member spread becomes the CI ("model
democracy"). Decadal statistic = MEAN (expected annual burnt-area %). Higher = worse.
12-member ensemble => NO spatial smoothing (contrast let). Shared 2020s baseline across
all three scenarios; ISIMIP3b starts in 2015, so the layer begins at the 2020s baseline
with no full 2010s decade (same as the csoil layer).

The three models do NOT share a land mask (mc2-usfs 58,919 cells; visit 58,714;
classic 66,644), so n_members / n_models are emitted per cell and decade. Thin-coverage
cells are RETAINED, not masked: a decade pools 10 years per member, so even a 2-member
cell rests on 20 annual samples (user decision 2026-07-27).

KNOWN OPEN ISSUE -- `visit` HIGH-LATITUDE BIAS (2026-07-28, undecided). The 5 `visit`
members have an INVERTED zonal profile: 2090s ssp585 land-mean burnt %/yr is 4.76 at
45-60N, 3.36 at 60-70N, 5.95 at 70-75N and 25.94 above 75N, against mc2-usfs 0.65-1.45
and classic 0.00-1.15. visit thus burns more above 75N than in the tropics (7.3%), which
is not physical; its worst cells are visit-only Arctic islands saturated near 100%/yr.
The global area-weighted total is unaffected (visit 4.56 Mkm2/yr, closest of the three to
GFED4) because polar cells carry almost no area -- exactly why no aggregate check catches
it, and why the layer declares `zonal_expectation` so the QA report warns. **All 12
members are retained and nothing is masked**, pending review of these results. Treat
values poleward of ~70N as unreliable meanwhile. See `known_issues` in the output attrs.

Output: burntarea_{scenario}_processed.nc with variables
{median, percentile, trend, lower_ci, upper_ci, n_members, n_models} on
(decade, lat, lon), units %.
"""

import os
import re
import sys
import warnings
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.contact_sheet import render_contact_sheet  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402
from utils.trend_significance import (  # noqa: E402
    METHOD as SIGNIFICANCE_METHOD, TREND_METHOD, AnnualEnsembleMean, mk_expanding,
    significance_definition, theilsen_decadal, trend_definition_decadal,
)

VAR = "burntarea-total"
LAYER_ID = "wildfire_burntarea_annual"
RAW_PATTERN = "*_burntarea-total_global_*_2015_2100.nc"
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099
# Reference ceiling for reporting only -- NOT a clamp. Annual burnt area is a
# CUMULATIVE fraction, so a cell that reburns within a year legitimately exceeds
# 100%: summing monthly members yields a small real tail (visit 0.15% of annual
# values, all <=120%; classic 0.05%, a few to ~182%; decadal means peak at ~107%).
# Clamping upper_ci to 100 would therefore drive it BELOW the median on exactly
# those cells and break CI ordering, so the CI is floored at 0 and left unbounded
# above.
VREF = 100.0
TWO_TIER_ZERO_THRESHOLD = 0.02   # use two-tier percentile if >2% of baseline is exact 0
MONTHS_PER_YEAR = 12

# elm-eca declares the 0.5 deg ISIMIP grid but is effectively ~4x5 deg: seams recur
# every 10 columns / 8 rows. It slipped a full tabular value-check and 37 algebraic
# checks on the csoil layer because distribution statistics are invariant under
# spatial rearrangement -- only rendering the member revealed it (GUARDRAILS S11).
# Detection traps that made it hard: an origin-aligned 2x2 test found only 3.3%; an
# exact-tie test scored it clean because its blocks are SMOOTH inside rather than
# constant; a modulo test capped at k=2 inverted at k=10; and a variogram was
# confounded by its much larger sill.
EXCLUDED_MODELS = {"elm-eca": "effectively ~4x5 deg, not 0.5 deg -- far too coarse"}


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (model, gcm, scenario, member) from a standard ISIMIP filename.

    e.g. visit_gfdl-esm4_w5e5_ssp585_2015soc_default_burntarea-total_global_monthly_2015_2100.nc
         [0]=model  [1]=gcm  [2]=bias  [3]=scenario  ...
    (model/gcm names contain hyphens but no underscores, so the split is safe.)
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Parse integer calendar years from a CF time axis (xarray cannot decode these
    non-standard-calendar files). All ISIMIP3b members here use 'days since 1601' on a
    365_day calendar.

    Decoding is delegated to cftime, which applies each file's own calendar exactly,
    instead of dividing by an assumed days-per-year. That matters for the two MONTHLY
    members: under days/365 arithmetic December lands at ~year+0.96, and naive ROUNDING
    pushes it into the following year -- silently misassigning one month in twelve and
    corrupting both the first and last decade. cftime returns the true calendar year.

    'years since' / 'months since' are not valid udunits for cftime (ISIMIP2b uses them,
    ISIMIP3b does not), so they are handled arithmetically with a FLOOR -- never a
    round -- for the same reason.
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
        yrs = base + (vals if step == "years" else vals / 12.0)
        return np.floor(yrs).astype(int)
    dates = cftime.num2date(t.values, units, calendar=cal)
    return np.array([d.year for d in np.atleast_1d(dates)], dtype=int)


def load_member(fpath):
    """Load one member as annual (year, lat, lon), restricted to MIN..MAX_YEAR.

    visit and classic publish burntarea-total only MONTHLY, so those files are collapsed
    to one value per year by SUMMING the 12 months -- burnt area accumulates over its
    reporting interval (verified against CLASSIC's daily output; see module docstring).
    A year that does not carry all 12 months is DROPPED rather than summed, because a
    partial sum is a silent under-count that no downstream check would catch.

    ``skipna=False`` is deliberate: with skipna=True an all-NaN ocean cell sums to 0.0,
    which is finite, so the ocean would be admitted as zero-burn land and drag every
    spatial statistic down.
    """
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    if da.year.size != np.unique(da.year.values).size:      # monthly -> annual
        uy, counts = np.unique(da.year.values, return_counts=True)
        partial = uy[counts != MONTHS_PER_YEAR]
        if partial.size:
            log(f"    NOTE: dropping {partial.size} incomplete year(s) "
                f"{partial.min()}-{partial.max()} (fewer than {MONTHS_PER_YEAR} months)")
            da = da.isel(year=np.where(np.isin(da.year.values, uy[counts == MONTHS_PER_YEAR]))[0])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            da = da.groupby("year").sum("year", skipna=False)
    return da.sortby("year")


def decade_mean_map(da, decade):
    """Mean over a decade window = per-member annual burnt-area-% map (lat, lon)."""
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return da.isel(year=sel).mean("year").values.astype(np.float32)


def anchored_trend(med_stack, i, b_idx):
    """Baseline-anchored trend: the mean rate of change of the decadal median
    FROM the baseline decade TO decade i, in units of value per DECADE:

        trend[i] = (median[i] - median[baseline]) / (i - baseline)

    So the 2090s panel shows the full baseline->2090s trend, the 2050s panel the
    baseline->2050s trend, etc. Anchored at the 2020s origin, so it is exactly the
    (decade - baseline) change map divided by the elapsed decades -> spatially
    coherent by construction and consistent with the change map (a within-decade
    slope of the annual series is NOT used: fire is far too noisy year-to-year,
    which yields a spotty sign-flipping field). The baseline decade itself has no
    elapsed change -> 0 (identical across scenarios, preserving bit-identity)."""
    span = i - b_idx
    if span == 0:
        return np.zeros(med_stack.shape[1:], np.float32)
    return ((med_stack[i] - med_stack[b_idx]) / span).astype(np.float32)


def make_pct_fn(baseline_flat):
    """Percentile-of-score against the 2020s baseline (higher = worse).

    Chooses the tier empirically from the baseline's exact-zero fraction:

    * SINGLE-tier: rank each value against the FULL 2020s land distribution -> [1, 100].

    * TWO-tier (zero-inflated, like let): value 0 -> 1; value > 0 -> ranked against the
      NON-ZERO baseline cells -> [2, 100]. Used if the baseline is materially
      zero-inflated (>2% exact zeros), which would otherwise crush the gradient.
      burntarea is strongly zero-inflated (all three models report large true-zero
      areas), so this is the expected branch here.
    """
    frac_zero = float(np.mean(baseline_flat == 0.0))
    two_tier = frac_zero >= TWO_TIER_ZERO_THRESHOLD

    if two_tier:
        nz_sort = np.sort(baseline_flat[baseline_flat > 0])
        n_nz = len(nz_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            res = np.ones(vals.shape, np.float32)  # zero-burn cells -> 1
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="left") / n_nz
                res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
            out[fin] = res
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
                out[fin] = np.clip(100.0 * frac, 1.0, 100.0).astype(np.float32)
            return out.reshape(arr.shape)

    return pct, ("two_tier" if two_tier else "single_tier"), frac_zero


def main():
    files = [str(p) for p in storage.stage_raw(LAYER_ID, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no burntarea member files found in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
        log("Upload members with storage.ingest_raw() first.")
        return

    # Drop excluded models (see EXCLUDED_MODELS) and say so loudly -- a silently
    # dropped member is indistinguishable from one that was never there.
    kept = [f for f in files if parse_name(f)["model"] not in EXCLUDED_MODELS]
    dropped = len(files) - len(kept)
    if dropped:
        for mdl, why in sorted(EXCLUDED_MODELS.items()):
            n = sum(1 for f in files if parse_name(f)["model"] == mdl)
            if n:
                log(f"EXCLUDING {n} {mdl} files -- {why}")
        log(f"  {len(kept)} of {len(files)} staged files retained")
    files = kept

    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log("=" * 70)
    log("Processing burntarea (wildfire, ISIMIP3b) -> TCFD 8-value-class format")
    log("=" * 70)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models} | GCMs: {gcms}")
    log("Monthly members annualized by SUM (burnt area accumulates); annual members as-is.")
    log("No normalization (raw %, equal-weight model democracy); no spatial smoothing.")

    da0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0

    # ---- Pass 1: per-member decadal-mean maps ---------------------------------
    dec = {s: {} for s in scenarios}            # dec[scen][member] = (n_dec, lat, lon)
    obs_max = 0.0
    # trend_pvalue is tested on the ensemble-mean ANNUAL series (GUARDRAILS S15), which
    # is gone once each member is reduced to decadal maps. `trend` needs only the
    # decadal medians (S10). Sorted order matches backfill_trend_significance.py.
    annual_acc = {s: AnnualEnsembleMean(MIN_YEAR, MAX_YEAR, (LAT, LON))
                  for s in scenarios}
    for f in sorted(files):
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
        if da.year.size:
            annual_acc[s].add(da.year.values, da.values)
        if da.year.size == 0:
            log(f"  WARNING: {os.path.basename(f)} has no years in "
                f"{MIN_YEAR}-{MAX_YEAR}; skipping (check time-axis parse)")
            continue
        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        for i, d in enumerate(DECADES):
            raw = decade_mean_map(da, d)
            if raw is not None:
                maps[i] = raw
        dec[s][member] = maps
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mx = float(np.nanmax(maps))
        obs_max = max(obs_max, mx)
        log(f"  loaded {info['model']:<18} {info['gcm']:<14} {s}  "
            f"years={da.year.size:<3d} max={mx:.2f}%")
        del da

    log(f"\nObserved max decadal-mean annual burnt area across all members: {obs_max:.3f}%")
    if obs_max > VREF:
        log(f"  NOTE: exceeds {VREF:.0f}% -- expected. Annual burnt area is cumulative, so a "
            f"cell that reburns within a year exceeds 100%. Values are NOT clamped.")

    # ---- Shared 2020s baseline (ensemble mean across all members/scenarios) ---
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec[s]})
    shared_2020 = []
    for member in all_members:
        per_scen = [dec[s][member][b_idx] for s in scenarios if member in dec[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_2020.append(np.nanmean(np.stack(per_scen, 0), axis=0))
    shared_2020 = np.stack(shared_2020, 0).astype(np.float32)  # (member, lat, lon)

    # ---- Per-member contact sheet (GUARDRAILS S11) -----------------------------
    # Rendered here because this is the only point where every member's own field
    # still exists separately; once pooled, an individual bad member is diluted and
    # invisible. CLASSIC's 1 deg blocks should be plainly visible in its two panels.
    sheet_path = None
    try:
        sheet_path = render_contact_sheet(
            {m: shared_2020[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", LAYER_ID, BASELINE_DECADE, units="%",
            note=(f"Ensemble as processed: {len(all_members)} members. "
                  f"classic runs at an effective 1.0 deg (2x2 constant blocks) and is "
                  f"retained deliberately. "
                  + (f"Excluded and NOT shown: "
                     f"{', '.join(sorted(EXCLUDED_MODELS))}." if EXCLUDED_MODELS else "")))
        log(f"  contact sheet: {sheet_path}  <-- REVIEW THIS before trusting the layer")
    except Exception as e:                                      # noqa: BLE001
        log(f"  WARNING: contact sheet failed ({type(e).__name__}: {e}); "
            f"the per-member visual check of GUARDRAILS S11 has NOT been produced")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared_median = np.nanmean(shared_2020, axis=0)
        shared_sd = np.nanstd(shared_2020, axis=0)
    shared_lo = np.clip(shared_median - shared_sd, 0, None)
    shared_hi = shared_median + shared_sd

    # Per-cell coverage: the 3 models do not share a land mask, so record how many
    # members and how many DISTINCT models back each cell. Thin cells are retained,
    # not masked -- a decade pools 10 years per member, so a 2-member cell still
    # rests on 20 annual samples. Emitting the counts makes the CI auditable.
    shared_nmem = np.sum(np.isfinite(shared_2020), axis=0).astype(np.float32)
    member_models = np.array([m.split("_")[0] for m in all_members])
    shared_nmodel = np.zeros_like(shared_nmem)
    for mdl in np.unique(member_models):
        has = np.any(np.isfinite(shared_2020[member_models == mdl]), axis=0)
        shared_nmodel += has.astype(np.float32)
    land = shared_nmem > 0
    n_land = int(land.sum())

    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"land cells n={len(baseline_flat):,}, exact-zero fraction={frac_zero:.2%}, "
        f"percentile mode={pct_mode}, global-mean burnt%={np.nanmean(shared_median):.4f}")
    log(f"  coverage: {n_land:,} land cells; "
        f"all-member {100*np.mean(shared_nmem[land]==len(all_members)):.2f}%, "
        f">=2 models {100*np.mean(shared_nmodel[land] >= 2):.2f}%, "
        f"single-model {100*np.mean(shared_nmodel[land] == 1):.2f}%, "
        f"single-member {int(np.sum(shared_nmem[land] == 1))} cells")

    # ---- Phase A: per-scenario median / CI / percentile (no trend yet) --------
    (med_by_scen, lo_by_scen, hi_by_scen, pct_by_scen,
     nmem_by_scen, nmodel_by_scen, members_by_scen) = {}, {}, {}, {}, {}, {}, {}
    for s in scenarios:
        members = sorted(dec[s])
        members_by_scen[s] = members
        mdl_of = np.array([m.split("_")[0] for m in members])
        stack = np.stack([dec[s][m] for m in members], 0)  # (member, dec, lat, lon)
        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)
        percentile = np.full_like(median, np.nan)
        nmem = np.zeros_like(median)
        nmodel = np.zeros_like(median)
        for i, d in enumerate(DECADES):
            layer = stack[:, i, :, :]  # (member, lat, lon)
            nmem[i] = np.sum(np.isfinite(layer), axis=0)
            for mdl in np.unique(mdl_of):
                nmodel[i] += np.any(np.isfinite(layer[mdl_of == mdl]), axis=0)
            if d == BASELINE_DECADE:
                # Baseline is the SHARED cross-scenario map, so its coverage is the
                # shared one too -- not this scenario's slice.
                median[i], lower[i], upper[i], percentile[i] = (
                    shared_median, shared_lo, shared_hi, shared_pct)
                nmem[i], nmodel[i] = shared_nmem, shared_nmodel
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(layer, axis=0)
                sd = np.nanstd(layer, axis=0)
            lower[i] = np.clip(median[i] - sd, 0, None)
            upper[i] = median[i] + sd
            percentile[i] = pct(median[i])
        nmem[nmem == 0] = np.nan       # off-land: NaN, not a misleading 0
        nmodel[~np.isfinite(nmem)] = np.nan
        med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s] = (
            median, lower, upper, percentile)
        nmem_by_scen[s], nmodel_by_scen[s] = nmem, nmodel

    # ---- Phase B: baseline-anchored trend + write -----------------------------
    for s in scenarios:
        log(f"\n{'='*70}\nAssembling scenario {s}\n{'='*70}")
        members = members_by_scen[s]
        median, lower, upper, percentile = (
            med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s])
        nmem, nmodel = nmem_by_scen[s], nmodel_by_scen[s]
        # Theil-Sen on the DECADAL medians (S10) -- NOT the annual series: burnt area
        # is zero over most of the land in most years, and a median of pairwise annual
        # slopes is exactly 0 once over half the year-pairs are 0-to-0.
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
            log(f"  {d}s: {tag:<15}  global-mean burnt%={np.nanmean(median[i]):.4f}  "
                f"trend={np.nanmean(trend[i]):+.4f}%/dec")

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
                "n_models": (["decade", "lat", "lon"], nmodel),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": "Annual burnt area (percent of grid cell)",
                "units": "%",
                "statistic": "decadal_mean_annual_burnt_area_percent",
                "value_note": ("median = equal-weighted ensemble mean over 3 fire models "
                               "(mc2-usfs, visit, classic) x their GCMs of the annual burnt-"
                               "area percent; raw ISIMIP3b values in [0,100] %."),
                "annualization": ("monthly members (visit, classic) are annualized by SUMMING "
                                  "the 12 months -- burnt area ACCUMULATES over its reporting "
                                  "interval. Verified against classic, which publishes the same "
                                  "run daily and monthly: each monthly value equals the sum of "
                                  "that month's daily values (1e-6 agreement, r=1.00000000). "
                                  "A mean would under-scale the layer 12x. Years without all "
                                  "12 months are dropped, not partially summed. mc2-usfs is "
                                  "published annually and is used as-is."),
                "normalization": ("none -- all 3 models report the SAME unit (% burnt area) AND "
                                  "agree within 1.6x on the global mean, so they are equal-"
                                  "weighted in raw %; inter-member spread is retained as the CI "
                                  "(model-democracy decision). NOTE the unit string alone is not "
                                  "sufficient grounds: 2b clm45/orchidee also declare '%' while "
                                  "sitting ~1000x lower on a 0-1 fraction scale."),
                "model_notes": ("classic runs at an EFFECTIVE 1.0 deg -- 100.0% of its 2x2 cell "
                                "blocks are exactly constant -- and is retained deliberately "
                                "(a 1 deg member still carries real information). mc2-usfs is "
                                "strongly zero-inflated (~44% exact zeros) and publishes only "
                                "the 'nat' natural-vegetation run, so its land use differs from "
                                "visit/classic (2015soc). elm-eca is EXCLUDED as effectively "
                                "~4x5 deg. Land masks differ across models -- see n_members / "
                                "n_models."),
                "excluded_models": ("; ".join(f"{k}: {v}" for k, v in sorted(EXCLUDED_MODELS.items()))
                                    or "none"),
                "spatial_smoothing": "none (12-member ensemble is thick enough)",
                "trend_definition": trend_definition_decadal(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "trend_method": TREND_METHOD,
                "significance_method": SIGNIFICANCE_METHOD,
                "significance_definition": significance_definition(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "significance_pooling": (
                    "flat mean across members within each year; the p-value is "
                    "tested on that ANNUAL series while `trend` is fitted on the "
                    "DECADAL medians (GUARDRAILS S10)"),
                "trend_units": "% decade-1",
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-member standard "
                                  "deviation (across the 12 model x GCM members), floored at 0 "
                                  "and UNBOUNDED above -- annual burnt area is cumulative, so a "
                                  "reburning cell legitimately exceeds 100% and clamping the "
                                  "upper CI to 100 would push it below the median there. "
                                  "Spread reflects fire-model + GCM uncertainty. "
                                  "The spread is taken across MEMBERS after each member's decade "
                                  "is reduced over its 10 years. Use n_members / n_models to "
                                  "qualify the CI on thin cells."),
                "coverage_note": ("n_members / n_models give, per cell and decade, how many of "
                                  "the 12 members and how many of the 3 distinct models supplied "
                                  "a value. The models do not share a land mask. Thin-coverage "
                                  "cells are RETAINED, not masked: a decade pools 10 years per "
                                  "member, so even a 2-member cell rests on 20 annual samples."),
                "percentile_baseline": (f"{pct_mode}: ranked against the 2020s ensemble-mean "
                                        "land distribution (single-tier -> [1,100]; two-tier "
                                        "reserves 1 for exact-zero cells and maps positives "
                                        "to [2,100])."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "zonal_expectation": ("low_latitude_dominated -- observed burned area is "
                                      "overwhelmingly tropical/subtropical savanna and is "
                                      "near zero poleward of ~70 deg. Declaring this makes "
                                      "generate_qa_report.py warn if a polar latitude band "
                                      "exceeds the tropical band."),
                "known_issues": (
                    "OPEN (2026-07-28, no decision yet -- do not treat as resolved): the 5 "
                    "`visit` members carry a systematic HIGH-LATITUDE BIAS with an INVERTED "
                    "zonal profile. 2090s ssp585 land-mean burnt %/yr by band -- "
                    "visit: 45-60N 4.76, 60-70N 3.36, 70-75N 5.95, >75N 25.94; "
                    "mc2-usfs: 0.65 / 1.35 / 1.45 / 0.67; classic: 1.15 / 0.10 / 0.02 / 0.00. "
                    "visit therefore reports MORE burning above 75N (25.9%) than in the "
                    "tropics (7.3% at 23S-0), which is not physical: the worst cells are "
                    "visit-only Arctic islands (Franz Josef Land, Severnaya Zemlya, Novaya "
                    "Zemlya) saturated near 100%/yr, e.g. 81.25N/56.75E reads 0.000% in the "
                    "2020s and 100.04% in the 2090s on a single model. Poleward of 75N the "
                    "pooled median has 320 cells above 20%/yr. The GLOBAL area-weighted "
                    "figure is unaffected (visit 4.56 Mkm2/yr, the closest of the three to "
                    "GFED4's ~3.5) because polar cells carry negligible area -- which is "
                    "precisely why an aggregate check cannot see this. All 12 members are "
                    "RETAINED and nothing is masked, pending review of these results; "
                    "affected cells are auditable via n_members / n_models. Treat "
                    "burnt-area values poleward of ~70N as unreliable for now."),
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "baseline_note": ("ISIMIP3b begins in 2015, so the layer starts at the 2020s "
                                  "baseline; there is no full 2010s decade."),
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "n_models": len(models),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_sens": ("mc2-usfs nat/default; visit 2015soc/default; classic "
                             "2015soc/default (transient CO2 throughout). classic is pinned to "
                             "2015soc -- NOT the 2015soc-from-histsoc used by the csoil layer -- "
                             "because for burntarea that variant is MIXED-SCALE: its gfdl-esm4 "
                             "files are written ~100x low (0-1 fraction) while its ukesm1-0-ll "
                             "files are percent, both declaring units='%'. Value-check the "
                             "soc/sens variant per variable, never inherit it."),
                "source_dataset": "ISIMIP3b OutputData/biomes (burntarea-total)",
                "observational_context": ("GFED4 observes ~348 Mha/yr = ~3.5 Mkm2/yr burned "
                                          "globally; the 3 retained members give 4.6-7.3 "
                                          "Mkm2/yr, i.e. the right order of magnitude."),
                "description": ("Wildfire burnt-area processed to TCFD 8-value-class format "
                                "with shared 2020s baseline; 3-model x GCM ensemble (12 members "
                                "per scenario) in raw %, no normalization, no spatial smoothing. "
                                "ISIMIP3b/SSP supersedes the earlier ISIMIP2b/RCP build."),
            },
        )
        path = out_dir / f"burntarea_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID,
        stage,
        created_by="scripts/process_burntarea_fire.py",
        notes=("ISIMIP3b burntarea-total, ssp126/370/585; 12 members/scenario from "
               "mc2-usfs (annual, 5 GCMs) + visit (monthly, 5) + classic (monthly, 2); "
               "monthly annualized by SUM; elm-eca excluded (~4x5 deg). "
               "See WORKFLOW-ISSUES.md 2026-07-28; GUARDRAILS S9-S11."),
    )

    # Every ingest-and-process run leaves reviewable HTML evidence behind.
    # Never raises: the data is already published and gated, and these
    # artifacts are regenerable (see scripts/utils/finalize.py).
    log("\nGenerating QA/QC report and maps...")
    finalize_layer(LAYER_ID, version=version,
                   extra_maps=[sheet_path] if sheet_path else None)

    log("\nDone.")


if __name__ == "__main__":
    main()
