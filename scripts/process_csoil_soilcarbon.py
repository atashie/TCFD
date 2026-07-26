"""Process csoil-total (soil organic carbon stock) into the TCFD 6-value-class format.

csoil-total = the total soil organic carbon pool of each land grid cell, reported
by ISIMIP3b biomes (vegetation) models in kg C m-2. This is the direct subsurface
carbon-STORAGE signal (distinct from the vegetation pools cveg/croot/cvegbg and
from the net-sink FLUX nbp; and distinct from the Lange 2020 exposure family).

Ensemble: 3 ISIMIP3b models x their CMIP6 GCMs x {ssp126, ssp370, ssp585} =
12 members per scenario. The three models were value-checked (2026-07-25), all in
the SAME unit (kg C m-2) with COMPARABLE magnitudes (medians within ~2x):

    model             GCMs  CO2 run     median  mean   max   %zero  long_name
    classic            2    default*    5.76    6.67   70.3  4.48   "Carbon Mass in Soil Pool"
    jules-es-vn6p3     5    2015co2**   10.28   10.69  67.9  0.00   "Carbon Mass in Soil Pool"
    mc2-usfs-r87g5c1   5    default*    7.67    8.42   37.8  0.00   "Carbon Mass in Soil Pool"

   * default        = transient CO2 (SSP-consistent, includes CO2 fertilization).
  ** 2015co2 (FIXED CO2) is the ONLY run JULES-ES-VN6P3 publishes for csoil-total,
     so the ensemble MIXES CO2 treatments -- JULES's soil-carbon TREND is muted
     (no fertilization signal) relative to the two transient models. Land use is
     held fixed for all (2015soc-from-histsoc for classic/jules; nat for mc2-usfs,
     a natural-vegetation biome model with no land-use forcing). User decision
     2026-07-25: retain all 12 members (accept the mixed CO2) rather than drop
     JULES and thin the ensemble to 2 structurally-diverse models.

All three share the SAME unit and comparable magnitudes, so -- like burntarea and
unlike the water-index TWS -- NO normalization is applied: members are equal-
weighted in raw kg C m-2 and the inter-member spread becomes the CI ("model
democracy"; note CLASSIC contributes only 2 GCMs vs 5 each for the others).

Decadal statistic = MEAN over the decade window (soil carbon is a smooth, slowly
varying STOCK). Direction: higher = BETTER (more stored carbon), so the risk is
LOSS -- the percentile is INVERTED (low stock / decline -> high risk percentile)
and change/trend maps read red = soil-carbon decline. Thick 12-member ensemble
=> NO spatial smoothing. Shared 2020s baseline across ssp126/370/585.

Time axis: all members use "days since 1601" but different calendars (365_day for
classic/mc2-usfs; proleptic_gregorian for jules) -> calendar-aware year parse.
ISIMIP3b csoil starts in 2015, so there is no full 2010s decade; the layer begins
at the 2020s baseline (decades 2020s-2090s).

Output: csoil_{scenario}_processed.nc with variables
{median, percentile, trend, lower_ci, upper_ci} on (decade, lat, lon), units kg m-2.
"""

import os
import re
import glob
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

VAR = "csoil-total"
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
         [0]=model  [1]=gcm  [2]=bias  [3]=scenario  [4]=soc  [5]=sens  [6]=variable ...
    (model/gcm names contain hyphens but no underscores, so the split is safe.)
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Parse integer calendar years from a CF time axis (xarray cannot decode
    these non-standard-calendar files). All members use 'days since 1601' but
    disagree on the calendar:
      classic/mc2-usfs -> 365_day             (divide by 365)
      jules-es-vn6p3   -> proleptic_gregorian (divide by 365.25)
    so we handle years / days / months and use the calendar's days-per-year."""
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
    else:  # days
        dpy = 360.0 if "360" in cal else 365.0 if "365" in cal else 365.25
        yrs = base + vals / dpy
    return np.round(yrs).astype(int)


def load_member(fpath):
    """Load one member file as (year, lat, lon), restricted to MIN..MAX_YEAR."""
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    return da


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
      -> [1, 100]. Correct when zeros are negligible (pooled ensemble is ~0.7%
      zeros: only CLASSIC has ~4.5% desert/ice zeros, over 2 of 12 members).

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
    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / "soilcarbon_csoil_annual"
    out_dir = root / "data" / "processed" / "soilcarbon_csoil_annual"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / "*_csoil-total_global_annual_2015_2100.nc")))
    if not files:
        log(f"ERROR: no csoil-total member files found in {raw_dir}")
        return
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log("=" * 66)
    log("Processing csoil-total (soil organic carbon) -> TCFD 6-value-class format")
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
    for f in files:
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
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

    # ---- Phase A: per-scenario median / CI / percentile (no trend yet) --------
    med_by_scen, lo_by_scen, hi_by_scen, pct_by_scen, members_by_scen = {}, {}, {}, {}, {}
    for s in scenarios:
        members = sorted(dec[s])
        members_by_scen[s] = members
        stack = np.stack([dec[s][m] for m in members], 0)  # (member, dec, lat, lon)
        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)
        percentile = np.full_like(median, np.nan)
        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                median[i], lower[i], upper[i], percentile[i] = (
                    shared_median, shared_lo, shared_hi, shared_pct)
                continue
            layer = stack[:, i, :, :]  # (member, lat, lon)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(layer, axis=0)
                sd = np.nanstd(layer, axis=0)
            lower[i] = np.clip(median[i] - sd, CI_FLOOR, None)
            upper[i] = np.clip(median[i] + sd, CI_FLOOR, None)
            percentile[i] = pct(median[i])
        med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s] = (
            median, lower, upper, percentile)

    # ---- Phase B: baseline-anchored trend + write -----------------------------
    for s in scenarios:
        log(f"\n{'='*66}\nAssembling scenario {s}\n{'='*66}")
        members = members_by_scen[s]
        median, lower, upper, percentile = (
            med_by_scen[s], lo_by_scen[s], hi_by_scen[s], pct_by_scen[s])
        trend = np.full_like(median, np.nan)
        for i, d in enumerate(DECADES):
            trend[i] = anchored_trend(median, i, b_idx)  # baseline->decade rate
            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(members)} members"
            log(f"  {d}s: {tag:<15}  global-mean csoil={np.nanmean(median[i]):.4f}  "
                f"trend(base->{d}s)={np.nanmean(trend[i]):+.4f} kg m-2/dec")

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": "Soil organic carbon stock (total soil carbon pool)",
                "units": "kg m-2",
                "statistic": "decadal_mean_soil_carbon_kg_per_m2",
                "value_note": ("median = equal-weighted ensemble mean over 3 ISIMIP3b biomes "
                               "models (classic, jules-es-vn6p3, mc2-usfs) x their CMIP6 GCMs "
                               "of the total soil-carbon stock; raw ISIMIP3b values in kg C m-2."),
                "normalization": ("none -- all 3 models report the SAME unit (kg C m-2) with "
                                  "COMPARABLE magnitudes (2020s medians ~5.8/7.7/10.3), so they "
                                  "are equal-weighted in raw kg C m-2 and the inter-member spread "
                                  "is retained as the CI (model-democracy decision). Note CLASSIC "
                                  "contributes only 2 GCMs vs 5 each for jules/mc2-usfs, so it is "
                                  "underweighted in the flat 12-member mean."),
                "co2_treatment": ("MIXED: classic & mc2-usfs use transient CO2 (default run, "
                                  "SSP-consistent, includes CO2 fertilization); jules-es-vn6p3 "
                                  "publishes ONLY the fixed-2015-CO2 run for csoil-total, so its "
                                  "soil-carbon trend is muted (no fertilization). Retained to "
                                  "keep 12 members (user decision 2026-07-25)."),
                "model_notes": ("long_name is consistent & correct across members ('Carbon Mass "
                                "in Soil Pool'). Time axis is 'days since 1601' for all, but the "
                                "calendar differs (365_day for classic/mc2-usfs, "
                                "proleptic_gregorian for jules) -> calendar-aware year parse. "
                                "classic has ~4.5% exact-zero (desert/ice) cells; mc2-usfs is a "
                                "coarse natural-vegetation biome model (nat soc, no land use) and "
                                "the most compressed (max ~38 vs ~68-70 kg m-2)."),
                "spatial_smoothing": "none (12-member ensemble is thick enough)",
                "trend_definition": ("baseline-anchored rate: trend[decade] = (median[decade] "
                                     "- median[2020s]) / (elapsed decades), units kg C m-2 per "
                                     "decade. Each decade's panel is the trend FROM the 2020s "
                                     "baseline TO that decade (2090s = the full baseline->2090s "
                                     "trend). Built on decadal means and anchored at the "
                                     "baseline, so it is exactly the (decade-2020s) change map / "
                                     "elapsed decades -- spatially coherent by construction and "
                                     "consistent with the change map (GUARDRAILS S10). The 2020s "
                                     "baseline has no elapsed change -> trend 0 (identical across "
                                     "scenarios). A NEGATIVE trend = soil-carbon LOSS = the "
                                     "adverse direction."),
                "trend_units": "kg m-2 decade-1",
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-member standard "
                                  "deviation (across the 12 model x GCM members), floored at 0. "
                                  "Spread reflects model + GCM (+ mixed-CO2) uncertainty."),
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
                "soc_sens": ("land use fixed: 2015soc-from-histsoc (classic/jules) / nat "
                             "(mc2-usfs); CO2: transient 'default' (classic/mc2-usfs), fixed "
                             "'2015co2' (jules)."),
                "source_dataset": "ISIMIP3b OutputData/biomes (csoil-total, annual)",
                "description": ("Soil organic carbon stock processed to TCFD 6-value-class format "
                                "with shared 2020s baseline; 3-model x CMIP6-GCM ensemble in raw "
                                "kg C m-2, no normalization, no spatial smoothing, higher_is_"
                                "better (risk = loss)."),
            },
        )
        path = out_dir / f"csoil_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  saved {path}")

    log("\nDone.")


if __name__ == "__main__":
    main()
