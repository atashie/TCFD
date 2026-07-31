"""Process csoil-total (soil organic carbon stock) into the TCFD 8-value-class format.

csoil-total = the total soil organic carbon pool of each land grid cell, reported
by ISIMIP3b biomes (vegetation) models in kg C m-2. This is the direct subsurface
carbon-STORAGE signal (distinct from the vegetation pools cveg/croot/cvegbg and
from the net-sink FLUX nbp; and distinct from the Lange 2020 exposure family).

Ensemble: 4 of the 5 ISIMIP3b biomes models that publish csoil-total x their CMIP6
GCMs x {ssp126, ssp370, ssp585} = 17 members per scenario. Every member was
value-checked 2026-07-27/28; all report the SAME unit (kg C m-2) with COMPARABLE
central magnitudes (2020s medians span 1.8x):

    model             GCMs  cadence   CO2 run    median   p99.9    max   %zero  eff.res
    classic            2    annual    default*     5.74    53.9   70.2   4.42   1.0 deg
    mc2-usfs-r87g5c1   5    annual    default*     7.60    26.1   37.2   0.00   0.5 deg
    visit              5    MONTHLY   default*     8.84   152.3  180.4   0.00   0.5 deg
    jules-es-vn6p3     5    annual    2015co2**   10.34    51.2   66.0   0.00   0.5 deg
    ---- EXCLUDED ----
    elm-eca            5    MONTHLY   default*     7.67   135.6  162.4   2.82   ~4x5 deg

  * default = transient CO2 (SSP-consistent, includes CO2 fertilization).
 ** 2015co2 (FIXED CO2) is the ONLY run JULES publishes for csoil-total, so the
    ensemble MIXES CO2 treatments. Contrary to this layer's earlier note, fixing
    CO2 does NOT mute JULES's trend -- it gives JULES the STRONGEST and most
    scenario-sensitive response of the five, because without fertilization nothing
    offsets warming-driven decomposition. Land use is fixed for every member.

Global-mean 2090s-vs-2020s change per model (own land mask), ssp126/370/585:

    classic          +0.65 / +0.77 / +0.77 %      transient CO2
    mc2-usfs         +0.40 / -0.11 / -0.05 %      transient CO2
    elm-eca          +2.01 / +2.43 / +2.60 %      transient CO2
    visit            +6.70 / +7.31 / +7.14 %      transient CO2
    jules-es-vn6p3   +1.57 / -2.92 / -4.37 %      FIXED 2015 CO2

So the SIGN of the global-mean change is contested: of the 4 RETAINED models 2
accumulate, 1 is flat, 1 loses strongly. Ensemble membership -- not processing --
drives the headline number: the original 3-model / 12-member layer gave
+1.06 / -1.43 / -2.17 %, and that result is reproducible from these same files to
within 0.05pp. Going to 4 models / 17 members cuts JULES's weight from 42% to 29%
and raises the ensemble mean (ssp585 measured at +0.68% vs -2.17%), so ssp370 and
ssp585 no longer read as net loss. Scenario ordering stays monotonic. Read the sign
of the global mean as contested across models, not settled.

elm-eca is EXCLUDED (user decision 2026-07-28): it declares 0.5 deg but is
effectively ~4 deg lat x 5 deg lon, and combined with the fattest tail of the five it
rendered as large bright rectangles that dominated the maps. Dropping it cost 5 land
cells (0.01%). See EXCLUDED_MODELS for the measurement and the detection traps.

visit publishes csoil-total ONLY monthly and is annualized by the mean of each year's
12 months. Verified immaterial: within-year CV 0.108%, |Dec - annual mean|/mean 0.103%
-- ~3 orders of magnitude below the inter-model spread. See load_member().

All 4 retained models share the SAME unit and comparable centres, so -- like
burntarea and unlike the water-index TWS -- NO normalization is applied: members are
equal-weighted in raw kg C m-2 and the inter-member spread becomes the CI ("model
democracy"). Caveats: CLASSIC contributes only 2 GCMs vs 5 each for the others, and
is the one member coarser than 0.5 deg (natively 1 deg, 2/17 = 12% of the weight);
and visit carries a much fatter peat tail than the other three (p99.9 152 vs 26-54),
concentrated in the subtropics (median latitude of its top 0.1% is 31.9N) -- deep
tropical peat, not boreal permafrost -- so the mean in peat cells stays
model-dependent.

Decadal statistic = MEAN over the decade window (soil carbon is a smooth, slowly
varying STOCK). Direction: higher = BETTER (more stored carbon), so the risk is
LOSS -- the percentile is INVERTED (low stock / decline -> high risk percentile)
and change/trend maps read red = soil-carbon decline. 17-member ensemble is thick
enough => NO spatial smoothing. Shared 2020s baseline across ssp126/370/585.

Coverage: the 4 models do NOT share a land mask (58,714 to 67,647 cells of
259,200), so most land carries all 17 members and a thin remainder as few as
one model. Those cells are retained, not masked; n_members / n_models are emitted
per cell so the CI can be qualified. See the coverage block in main().

Time axis: all members use "days since 1601" but different calendars (365_day for
classic/mc2-usfs/visit; proleptic_gregorian for jules) -> years are
decoded with cftime rather than days/365 arithmetic, which would misassign
December of each monthly member. ISIMIP3b csoil starts in 2015, so there is no
full 2010s decade; the layer begins at the 2020s baseline (decades 2020s-2090s).

Output: csoil_{scenario}_processed.nc with variables {median, percentile, trend,
lower_ci, upper_ci, n_members, n_models} on (decade, lat, lon), units kg m-2.
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
from utils.finalize import finalize_layer  # noqa: E402
from utils.trend_significance import (  # noqa: E402
    METHOD as SIGNIFICANCE_METHOD, TREND_METHOD, AnnualEnsembleMean, mk_expanding,
    significance_definition, theilsen_decadal, trend_definition_decadal,
)
from utils.contact_sheet import render_contact_sheet  # noqa: E402

VAR = "csoil-total"
LAYER_ID = "soilcarbon_csoil_annual"
# Matches BOTH cadences: 3 models publish csoil-total annually, 2 (elm-eca, visit)
# only monthly. Monthly members are annualized on load -- see load_member().
RAW_PATTERN = "*_csoil-total_global_*_2015_2100.nc"

# --- Excluded members -------------------------------------------------------
# elm-eca is EXCLUDED although it is ingested and matches RAW_PATTERN. It declares the
# standard 0.5 deg / 360x720 ISIMIP grid but its csoil-total field is effectively
# ~4 deg lat x 5 deg lon: gradient seams recur every 10 columns (5.0 deg, 62 of 75
# detected) and every 8 rows (4.0 deg), with only fine within-block variation on top.
# Measured 2026-07-28 after a user spotted large boxes in the published maps.
#
# It matters because elm-eca ALSO has by far the fattest tail (max ~160 vs 35-70
# kg C m-2 for the others), so its coarse boxes are biased high and render as bright
# rectangles that dominated the ensemble visually. At 5 of 22 members it alone inflated
# the 2020s global mean from 9.32 to 10.99 kg C m-2 (+18%).
#
# Dropping it costs 5 land cells (0.01%) because the other four models cover almost the
# same mask, and it takes the 5 deg seam ratio from 0.871 to 0.957. Kept in raw staging
# for provenance rather than deleted. User decision 2026-07-28.
#
# Note the trap: neither the declared dimensions, an exact-tie test (its blocks are
# smooth inside, only 2.9% exact ties), nor an origin-aligned modulo-gradient test
# capped at 2 deg detects this. See GUARDRAILS S9 and the seam-spacing check in
# scripts/generate_qa_report.py.
EXCLUDED_MODELS = {"elm-eca": "effectively ~4x5 deg, not 0.5 deg; fattest tail"}
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2015, 2099   # ISIMIP3b csoil covers 2015-2100 (2100 unused)
CI_FLOOR = 0.0                    # soil carbon is non-negative; no upper cap
TWO_TIER_ZERO_THRESHOLD = 0.02    # use two-tier percentile if >2% of baseline is exact 0
HIGHER_IS_BETTER = True           # more stored carbon = better; risk = loss -> invert pct


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (model, gcm, scenario, member) from a standard ISIMIP3b filename.

    e.g. classic_gfdl-esm4_w5e5_ssp585_2015soc-from-histsoc_default_csoil-total_global_annual_2015_2100.nc
         [0]=model [1]=gcm [2]=bias [3]=scenario [4]=soc [5]=sens [6]=variable
         [7]=region [8]=cadence [9]=start [10]=end
    (model/gcm names contain hyphens but no underscores, so the split is safe.)
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], soc=p[4], sens=p[5],
                cadence=p[8], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Parse integer calendar years from a CF time axis (xarray cannot decode
    these non-standard-calendar files). All 5 members use 'days since 1601' but
    disagree on the calendar: 365_day for classic/mc2-usfs/elm-eca/visit,
    proleptic_gregorian for jules-es-vn6p3.

    Decoding is delegated to cftime, which applies each file's own calendar
    exactly, instead of dividing by an assumed days-per-year. That matters for the
    two MONTHLY members: under days/365 arithmetic December lands at ~year+0.96,
    and naive ROUNDING pushes it into the following year -- silently misassigning
    one month in twelve and corrupting both the first and last decade. cftime
    returns the true calendar year.

    'years since' / 'months since' are not valid udunits for cftime (ISIMIP2b uses
    them, ISIMIP3b does not), so they are handled arithmetically with a FLOOR --
    never a round -- for the same reason.
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

    Members arrive at two cadences (elm-eca and visit publish csoil-total only
    monthly), so a monthly file is collapsed to one value per year by taking the
    ANNUAL MEAN of its 12 months. Mean -- not a December snapshot -- because:
      * it is the same statistic already used to reduce years to a decade, so the
        monthly members are reduced exactly like the annual ones rather than being
        handed a different estimator, and
      * soil carbon is a large, slowly-turning stock whose within-year swing is
        tiny next to the inter-model spread (value-checked; see module docstring),
        so the mean is a low-variance estimate of the year's stock and carries no
        seasonal phase bias, whereas a December snapshot would import each model's
        northern-winter state.
    Units are unchanged by the mean (kg C m-2 stays kg C m-2).
    """
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    if da.year.size != np.unique(da.year.values).size:      # monthly -> annual
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            da = da.groupby("year").mean("year", skipna=True)
    return da.sortby("year")


def decade_mean_map(da, decade):
    """Mean over a decade window = per-member soil-carbon stock map (lat, lon)."""
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
    coherent by construction and consistent with the change map (GUARDRAILS S10).
    The baseline decade itself has no elapsed change -> 0 (identical across
    scenarios, preserving bit-identity)."""
    span = i - b_idx
    if span == 0:
        return np.zeros(med_stack.shape[1:], np.float32)
    return ((med_stack[i] - med_stack[b_idx]) / span).astype(np.float32)


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the 2020s baseline, returned as a RISK score.

    The raw score ranks each value against the 2020s land distribution (higher raw
    value -> higher raw score). Because higher soil carbon is BETTER, the score is
    INVERTED to a risk percentile via (101 - raw) -> a cell with a LOW stock (or
    that declines below the baseline distribution) gets a HIGH risk percentile.

    Tier is chosen empirically from the baseline's exact-zero fraction:

    * SINGLE-tier (expected here): rank against the FULL 2020s land distribution
      -> [1, 100]. Correct when zeros are negligible. Only CLASSIC (~4.4% desert/
      ice zeros, 2 of 17 members) emits exact zeros at all -- the other three
      retained models emit none -- and the baseline is the ensemble MEAN, so a cell
      is 0 there only when every contributing member is 0. The realized fraction is
      far below CLASSIC's own rate. The mode is decided from the data, not assumed.

    * TWO-tier (zero-inflated): value 0 -> lowest stock; value > 0 -> ranked
      against the NON-ZERO baseline. Used only if the baseline is materially
      zero-inflated (>2% exact zeros). After inversion, exact-zero (barren)
      cells map to the highest risk.
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
            res = np.ones(vals.shape, np.float32)  # zero-stock cells -> raw 1
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


def main():
    files = [str(p) for p in storage.stage_raw(LAYER_ID, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no csoil-total member files found in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
        log("Upload members with storage.ingest_raw() first.")
        return

    # Drop excluded models (see EXCLUDED_MODELS) and say so loudly -- a silently
    # thinned ensemble is exactly the kind of thing that should never be implicit.
    kept = [f for f in files if parse_name(f)["model"] not in EXCLUDED_MODELS]
    dropped = len(files) - len(kept)
    if dropped:
        for mdl, why in sorted(EXCLUDED_MODELS.items()):
            n = sum(1 for f in files if parse_name(f)["model"] == mdl)
            log(f"EXCLUDING {n} {mdl} files -- {why}")
        log(f"  {len(kept)} of {len(files)} staged files retained")
    files = kept
    if not files:
        log("ERROR: every staged file was excluded; nothing to process.")
        return

    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log("=" * 66)
    log("Processing csoil-total (soil organic carbon) -> TCFD 8-value-class format")
    log("=" * 66)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models} | GCMs: {gcms}")
    log("No normalization (raw kg C m-2, equal-weight model democracy); no spatial smoothing.")
    log("Direction: higher_is_better (risk = soil-carbon LOSS) -> percentile inverted.")

    da0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0

    # ---- Pass 1: per-member decadal-mean maps ---------------------------------
    dec = {s: {} for s in scenarios}            # dec[scen][member] = (n_dec, lat, lon)
    # trend_pvalue is tested on the ensemble-mean ANNUAL series (GUARDRAILS S15), which
    # is gone once each member is reduced to decadal maps -- so accumulate it here.
    # `trend` itself needs only the decadal medians (S10). Sorted order matches
    # backfill_trend_significance.py so both paths agree bit-for-bit.
    annual_acc = {s: AnnualEnsembleMean(MIN_YEAR, MAX_YEAR, (LAT, LON))
                  for s in scenarios}
    for f in sorted(files):
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
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
        log(f"  loaded {info['model']:<18} {info['gcm']:<13} {s}")
        del da

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
    # invisible. Costs nothing extra to read -- the arrays are already in memory.
    sheet_path = None
    try:
        sheet_path = render_contact_sheet(
            {m: shared_2020[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", LAYER_ID, BASELINE_DECADE, units="kg C m-2",
            note=(f"Ensemble as processed: {len(all_members)} members. "
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
    shared_lo = np.clip(shared_median - shared_sd, CI_FLOOR, None)
    shared_hi = np.clip(shared_median + shared_sd, CI_FLOOR, None)

    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"land cells n={len(baseline_flat):,}, exact-zero fraction={frac_zero:.2%}, "
        f"percentile mode={pct_mode}, global-mean csoil={np.nanmean(shared_median):.4f} kg m-2")

    # ---- Coverage: the 4 models do NOT share a land mask ----------------------
    # Land-cell counts differ by model (visit 58,714 .. jules 67,647 of 259,200),
    # and each model's GCMs share its mask exactly, so the ensemble is taken over a
    # VARYING number of members per cell: most land carries all 17 members, while
    # the thin remainder (arid / ice-margin cells holding ~2% of global soil
    # carbon) is carried by as few as one model. Those cells are KEPT, not masked --
    # a decade supplies 10 annual samples per member, so the underlying spread
    # stays estimable even on a 2-member cell.
    #
    # The caveat is specific to THIS layer's CI estimator rather than to the data:
    # decade_mean_map() averages each member's 10 years into one decadal mean
    # BEFORE the spread is taken, so lower/upper_ci is a spread over member means
    # (n = members), not over year x member samples (n = 10 x members). On a
    # classic-only cell that is 2 samples, and on the 2 single-member cells the SD
    # collapses to 0. Kept as-is for cross-layer consistency (led / let / burntarea
    # all define the CI as mean +/- 1 inter-member SD); the per-cell member and
    # distinct-model counts are emitted so this is auditable rather than implicit.
    shared_nmem = np.sum(np.isfinite(shared_2020), axis=0).astype(np.float32)
    member_models = np.array([m.split("_")[0] for m in all_members])
    shared_nmodel = np.zeros_like(shared_nmem)
    for mdl in np.unique(member_models):
        has = np.any(np.isfinite(shared_2020[member_models == mdl]), axis=0)
        shared_nmodel += has.astype(np.float32)
    land = shared_nmem > 0
    n_land = int(land.sum())
    log(f"  coverage: {n_land:,} land cells; all-member {100*np.mean(shared_nmem[land]==len(all_members)):.2f}%, "
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
            lower[i] = np.clip(median[i] - sd, CI_FLOOR, None)
            upper[i] = np.clip(median[i] + sd, CI_FLOOR, None)
            percentile[i] = pct(median[i])
        nmem[nmem == 0] = np.nan       # off-land: NaN, not a misleading 0
        nmodel[~np.isfinite(nmem)] = np.nan
        med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s] = (
            median, lower, upper, percentile)
        nmem_by_scen[s], nmodel_by_scen[s] = nmem, nmodel

    # ---- Phase B: baseline-anchored trend + write -----------------------------
    for s in scenarios:
        log(f"\n{'='*66}\nAssembling scenario {s}\n{'='*66}")
        members = members_by_scen[s]
        median, lower, upper, percentile = (
            med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s])
        nmem, nmodel = nmem_by_scen[s], nmodel_by_scen[s]
        # Theil-Sen on the DECADAL medians (S10), p-value on the ANNUAL series (S15).
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
            log(f"  {d}s: {tag:<15}  global-mean csoil={np.nanmean(median[i]):.4f}  "
                f"trend={np.nanmean(trend[i]):+.4f} kg m-2/dec")

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
                "long_name": "Soil organic carbon stock (total soil carbon pool)",
                "units": "kg m-2",
                "statistic": "decadal_mean_soil_carbon_kg_per_m2",
                "value_note": ("median = equal-weighted ensemble mean over 4 ISIMIP3b biomes "
                               "models (classic, jules-es-vn6p3, mc2-usfs, visit) x their "
                               "CMIP6 GCMs of the total soil-carbon stock; raw ISIMIP3b "
                               "values in kg C m-2. A 5th model (elm-eca) publishes "
                               "csoil-total and was deliberately EXCLUDED -- see "
                               "excluded_sources."),
                "cadence_note": ("MIXED source cadence, harmonized on load: classic, "
                                 "jules-es-vn6p3 and mc2-usfs publish csoil-total ANNUALLY; "
                                 "visit publishes it only MONTHLY and is reduced to annual by "
                                 "taking the mean of each year's 12 months. Verified 2026-07-27 "
                                 "that the choice is immaterial: visit's within-year coefficient "
                                 "of variation is 0.108% and |December - annual mean| / mean is "
                                 "0.103% -- three orders of magnitude below the ~1.8x "
                                 "inter-model spread. The annual mean is unit-preserving "
                                 "(kg C m-2 -> kg C m-2) and matches the estimator already used "
                                 "to reduce years to decades."),
                "normalization": ("none -- all 4 retained models report the SAME unit (kg C "
                                  "m-2) with COMPARABLE magnitudes (2020s medians 5.74 classic / "
                                  "7.60 mc2-usfs / 8.84 visit / 10.34 jules = 1.8x spread), so "
                                  "they are equal-weighted in raw kg C m-2 and the inter-member "
                                  "spread is retained as the CI (model-democracy decision). Two "
                                  "caveats: CLASSIC contributes only 2 GCMs vs 5 each for the "
                                  "others, so it is underweighted in the flat 17-member mean; "
                                  "and disagreement is concentrated in the UPPER TAIL rather "
                                  "than the centre -- visit carries a much fatter peat tail "
                                  "(99.9th pct ~152, max ~180 kg C m-2) than classic / jules / "
                                  "mc2-usfs (~54 / ~51 / ~26, max ~70 / ~66 / ~37), "
                                  "concentrated in the subtropics (median latitude of its top "
                                  "0.1% is 31.9N) -- deep tropical peat, not boreal permafrost. "
                                  "Because members are pooled with a MEAN, that tail makes the "
                                  "pooled value more tail-sensitive than the 1.8x median "
                                  "agreement suggests, and the ensemble mean in peat cells stays "
                                  "materially model-dependent with a correspondingly wide CI. A "
                                  "cross-member MEDIAN would be robust to this; the mean is "
                                  "retained for consistency with the other TCFD layers."),
                "co2_treatment": ("MIXED: classic, mc2-usfs and visit use transient CO2 "
                                  "(default run, SSP-consistent, includes CO2 fertilization); "
                                  "jules-es-vn6p3 publishes ONLY the fixed-2015-CO2 run for "
                                  "csoil-total. CORRECTION (measured 2026-07-27, superseding the "
                                  "earlier note that JULES's trend is 'muted'): holding CO2 at "
                                  "2015 does NOT mute the trend, it AMPLIFIES the loss. JULES has "
                                  "the strongest and most scenario-sensitive response of the four "
                                  "(global-mean 2090s vs 2020s: +1.57% ssp126, -2.92% ssp370, "
                                  "-4.37% ssp585), because without fertilization there is no "
                                  "offsetting rise in litter input and warming-driven "
                                  "decomposition dominates. The four transient-CO2 models run "
                                  "flat-to-positive (classic +0.7%, mc2-usfs ~0%, visit +7.1% at "
                                  "ssp585). JULES is therefore the ensemble's main source of "
                                  "soil-carbon LOSS, and cutting its weight from 5/12 = 42% to "
                                  "5/17 = 29% by expanding the ensemble made the ensemble mean "
                                  "MORE POSITIVE, not less damped -- see ensemble_change."),
                "ensemble_change": ("This 17-member layer supersedes a 12-member version "
                                    "(classic + jules + mc2-usfs only). Recomputing that older "
                                    "3-model ensemble from the same raw files reproduces its "
                                    "published global-mean change to within 0.05pp (+1.06% / "
                                    "-1.43% / -2.17% vs the recorded +1.1% / -1.4% / -2.2%), so "
                                    "the difference is due to ensemble membership, not to a "
                                    "processing change. Adding visit (5 members) and dropping "
                                    "elm-eca raises the level relative to the 12-member run "
                                    "(ssp585 measured at +0.68% vs -2.17%), so ssp370/ssp585 no "
                                    "longer read as net loss. The scenario ORDERING is unchanged "
                                    "and monotonic (less carbon retained at higher forcing). Read "
                                    "the sign of the global mean as contested across models "
                                    "rather than settled: of the 4 retained models 2 gain, 1 is "
                                    "flat, 1 (jules) loses strongly."),
                "model_notes": ("Verified per-member 2026-07-27 across all 66 files: units are "
                                "'kg m-2' and long_name 'Carbon Mass in Soil Pool' for all 66 staged "
                                "no mislabelling (contrast burntarea's lpj-guess). Time axis is "
                                "'days since 1601' for all, but the calendar differs (365_day for "
                                "classic/mc2-usfs/elm-eca/visit; proleptic_gregorian for jules) "
                                "-> years are decoded with cftime, not days/365 arithmetic, which "
                                "would misassign December of each monthly member. Every member "
                                "spans 2015-2100 with exactly 1 or 12 steps per year. Land masks "
                                "DIFFER by model (visit 58,714 / mc2-usfs 58,919 / classic ~66,678 "
                                "/ jules 67,647 cells of 259,200) and each "
                                "model's GCMs share its mask exactly -- see n_members / n_models. "
                                "Exact zeros: classic ~4.4%, the other three 0%; no "
                                "negatives anywhere. mc2-usfs is a natural-vegetation biome model "
                                "(nat soc, no land use) and the most compressed."),
                "effective_resolution": ("Measured from the VALUES, not the declared dims -- a "
                                        "natively coarse model replicated onto the ISIMIP grid "
                                        "still reports 360x720 at 0.5 deg. Retained members: "
                                        "jules-es-vn6p3, mc2-usfs and visit are genuine 0.5 deg; "
                                        "classic is NATIVELY 1.0 deg x 1.0 deg, replicated 2x2 "
                                        "with a one-cell longitude offset (100% of 2x2 blocks "
                                        "constant at offset lat=0/lon=1; 99.35% of longitude runs "
                                        "exactly 2 cells), and carries 2/17 = 12% of the weight, "
                                        "leaving a mild 1-deg signature. No GCM contributes block "
                                        "structure. The excluded elm-eca was effectively "
                                        "~4 deg lat x 5 deg lon -- see excluded_sources. Detection "
                                        "requires care: exact-tie tests miss blocks that are "
                                        "smooth inside, and an origin-aligned modulo test misses "
                                        "offset blocks; use seam SPACING, as the check in "
                                        "scripts/generate_qa_report.py now does."),
                "excluded_sources": ("elm-eca is EXCLUDED from the ensemble entirely (user "
                                     "decision 2026-07-28) despite publishing csoil-total for 5 "
                                     "GCMs. It declares 0.5 deg but its field is effectively "
                                     "~4 deg lat x 5 deg lon: gradient seams recur every 10 "
                                     "columns (5.0 deg; 62 of 75 detected) and every 8 rows "
                                     "(4.0 deg), with only fine within-block variation on top. "
                                     "Combined with the fattest tail of the five (max ~160 vs "
                                     "35-70 kg C m-2) its coarse boxes were biased high and "
                                     "rendered as large bright rectangles that dominated the "
                                     "maps; at 5 of 22 members it alone inflated the 2020s global "
                                     "mean from 9.32 to 10.99 kg C m-2 (+18%). Dropping it cost 5 "
                                     "land cells (0.01%) and took the 5-deg seam ratio from 0.871 "
                                     "to 0.957. Its raw files remain in S3 staging for provenance. "
                                     "Separately, elm-eca also publishes a csoil-total copy under "
                                     "the ISIMIP3b "
                                     "PERMAFROST sector that is BYTE-IDENTICAL to its biomes copy "
                                     "(same sha512, 237,500,430 bytes; verified 2026-07-27). Only "
                                     "the biomes copy is ingested -- harvesting both would "
                                     "silently double-weight elm-eca. Also excluded: the "
                                     "nondep / nondep2015co2 / ssp585ndep nitrogen-deposition "
                                     "sensitivity runs, and duplicate soc variants, so each model "
                                     "contributes exactly one run per GCM x scenario."),
                "spatial_smoothing": "none (17-member ensemble is thick enough)",
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
                "trend_units": "kg m-2 decade-1",
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-member standard "
                                  "deviation (across the 17 model x GCM members), floored at 0. "
                                  "Spread reflects model + GCM (+ mixed-CO2) uncertainty. NOTE "
                                  "the spread is taken over member DECADAL MEANS -- each member's "
                                  "10 annual values are averaged first -- so n = members, not "
                                  "10 x members. Kept this way for consistency with the other "
                                  "TCFD layers (led / let / burntarea). Consequence on the thin "
                                  "part of the land mask: where only classic reports, the SD is "
                                  "over its 2 GCMs (one model), and on the 2 single-member cells "
                                  "it collapses to 0. Use n_members / n_models to qualify the CI."),
                "coverage_note": ("n_members / n_models give, per cell and decade, how many of "
                                  "the 17 member simulations and how many of the 4 distinct "
                                  "impact models actually reported there (NaN off-land). The "
                                  "models do not share a land mask, so most land carries all "
                                  "17 members while the remainder -- arid / ice-margin cells "
                                  "holding ~2% of global soil carbon -- is carried by as few as "
                                  "one model. Those cells are RETAINED, not masked: a decade "
                                  "supplies 10 annual samples per member, so the underlying "
                                  "spread remains estimable even on a 2-member cell (user "
                                  "decision 2026-07-27). The counts are emitted so any consumer "
                                  "needing a uniform-confidence subset can filter explicitly."),
                "percentile_baseline": (f"{pct_mode}: raw score ranks each cell against the 2020s "
                                        "ensemble-mean land distribution, then INVERTED to a risk "
                                        "percentile (101 - raw) because higher soil carbon is "
                                        "better -> low stock / decline = high risk percentile."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_better",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "decade_note": ("ISIMIP3b csoil covers 2015-2100, so there is no full 2010s "
                                "decade; the layer begins at the 2020s baseline (2020s-2090s)."),
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_sens": ("land use held fixed for every member: 2015soc-from-histsoc "
                             "(classic, jules), 2015soc (visit), nat (mc2-usfs, a "
                             "natural-vegetation model with no land-use forcing). CO2: transient "
                             "'default' (classic, mc2-usfs, visit), fixed '2015co2' "
                             "(jules -- its only published csoil-total run)."),
                "source_dataset": ("ISIMIP3b OutputData/biomes (csoil-total; annual for classic/"
                                   "jules-es-vn6p3/mc2-usfs, monthly annualized for visit)"),
                "description": ("Soil organic carbon stock processed to TCFD 8-value-class format "
                                "with shared 2020s baseline; 4-model x CMIP6-GCM ensemble "
                                "(17 members per scenario) in raw kg C m-2, no normalization, no "
                                "spatial smoothing, higher_is_better (risk = loss)."),
            },
        )
        path = out_dir / f"csoil_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID,
        stage,
        created_by="scripts/process_csoil_soilcarbon.py",
        notes=("17-member ensemble: 4 of the 5 ISIMIP3b models publishing csoil-total "
               "(classic 2 GCMs, jules-es-vn6p3 5, mc2-usfs 5, visit 5); visit "
               "annualized from monthly. elm-eca EXCLUDED -- declares 0.5 deg but is "
               "effectively ~4x5 deg. See WORKFLOW-ISSUES.md 2026-07-28; "
               "GUARDRAILS §8-§9, §11."),
    )

    # Every ingest-and-process run leaves reviewable HTML evidence behind.
    # Never raises: the data is already published and gated, and these
    # artifacts are regenerable (see scripts/utils/finalize.py).
    log("\nGenerating QA/QC report and maps...")
    finalize_layer(LAYER_ID, version=version, extra_maps=[sheet_path])

    log("\nDone.")


if __name__ == "__main__":
    main()
