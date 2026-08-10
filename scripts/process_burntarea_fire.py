"""SUPERSEDED 2026-08-10 -- do not use for new work.

    The shipped wildfire hazard layer is now `process_burntarea_isimip3b.py`
    (ISIMIP3b `burntarea-total`, ssp126/370/585, 22 members/scenario across 5 impact
    models x 5 CMIP6 GCMs). It supersedes this ISIMIP2b/RCP layer on every axis the
    selection criteria weigh: newer experiment generation (SSP over RCP), a deeper
    ensemble in BOTH dimensions (22 vs 12 members, 5 vs 3 impact models), and a
    high-forcing scenario (ssp585).

    This module is kept for provenance -- it documents the ISIMIP2b generation's
    per-member metadata defects (see below), which are real findings about that data
    and not reproduced in 3b. Do not extend it; do not read its framing decisions as
    precedent for a new layer (GUARDRAILS 9: measure, never inherit).

Process burntarea (wildfire burnt-area fraction) into the TCFD 6-value-class format.

burntarea = the fraction of each grid cell burned per year, reported by ISIMIP2b
vegetation/fire models in PERCENT [0, 100]. This is the direct biophysical fire
signal (NOT the Lange 2020 "exposure to unprecedented wildfire" family member lew).

Ensemble: 3 annual ISIMIP2b models x 4 CMIP5 GCMs x {rcp26, rcp60, rcp85} =
12 members per scenario. The three fire models were value-checked (2026-07-24):

    model      units  min   max   mean   %zero   long_name
    lpj-guess    %    0.1   97    0.71   0.00    "Fire Return Interval (burntarea)" *
    lpjml        %    0     80    0.57   1.05    "burnt area fraction"
    mc2-usfs     %    0    100    3.84   44.8    "Burnt Area Fraction"

  * lpj-guess's long_name is a MISLABEL: the values are unambiguously burnt-area
    percent (same units, same shape as the others; a 0.1-"year" return interval
    would mean burning 10x/yr). lpj-guess also floors reporting at 0.1% (never a
    true zero). mc2-usfs (a coarse biome model) is strongly zero-inflated and
    runs ~5-7x hotter in the mean.

All three share the SAME unit (% burnt area), so -- unlike the water-index TWS --
NO normalization is applied: the models are equal-weighted in raw %, and the
inter-member spread becomes the CI (user decision 2026-07-24, "model democracy").

Decadal statistic = MEAN (expected annual burnt-area %). Higher = worse (more
fire). Thick 12-member ensemble => NO spatial smoothing (contrast let). Shared
2020s baseline across rcp26/60/85. Percentile tier (single vs two-tier zero-
inflated) is chosen empirically from the pooled baseline's exact-zero fraction.

Output: burntarea_{scenario}_processed.nc with variables
{median, percentile, trend, lower_ci, upper_ci} on (decade, lat, lon), units %.
"""

import os
import re
import glob
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

VAR = "burntarea"
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2010, 2099
VMAX = 100.0                     # burnt area is a percent in [0, 100]
TWO_TIER_ZERO_THRESHOLD = 0.02   # use two-tier percentile if >2% of baseline is exact 0


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (model, gcm, scenario, member) from a standard ISIMIP filename.

    e.g. lpj-guess_gfdl-esm2m_ewembi_rcp26_2005soc_co2_burntarea_global_annual_2006_2099.nc4
         [0]=model  [1]=gcm  [2]=bias  [3]=scenario  ...
    (model/gcm names contain hyphens but no underscores, so the split is safe.)
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Parse integer calendar years from a CF time axis (xarray cannot decode
    these 365_day-calendar files). Members disagree on the unit:
      lpj-guess/lpjml -> 'years since 1661'   (offset is already in years)
      mc2-usfs        -> 'days since 1661'     (365_day: divide by 365)
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

    * SINGLE-tier (default here): rank each value against the FULL 2020s land
      distribution -> [1, 100]. Correct when zeros are negligible (lpj-guess
      floors at 0.1%, so the pooled ensemble mean is rarely exactly 0 on land).

    * TWO-tier (zero-inflated, like let): value 0 -> 1; value > 0 -> ranked
      against the NON-ZERO baseline cells -> [2, 100]. Used only if the baseline
      is materially zero-inflated (>2% exact zeros), which would otherwise crush
      the gradient.
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
    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / "wildfire_burntarea_annual"
    out_dir = root / "data" / "processed" / "wildfire_burntarea_annual"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / "*_2006_2099.nc4")))
    if not files:
        log(f"ERROR: no burntarea member files found in {raw_dir}")
        return
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log("=" * 66)
    log("Processing burntarea (wildfire) -> TCFD 6-value-class format")
    log("=" * 66)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models} | GCMs: {gcms}")
    log("No normalization (raw %, equal-weight model democracy); no spatial smoothing.")

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
        log(f"  loaded {info['model']:<16} {info['gcm']:<13} {s}")
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
    shared_lo = np.clip(shared_median - shared_sd, 0, VMAX)
    shared_hi = np.clip(shared_median + shared_sd, 0, VMAX)

    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"land cells n={len(baseline_flat):,}, exact-zero fraction={frac_zero:.2%}, "
        f"percentile mode={pct_mode}, global-mean burnt%={np.nanmean(shared_median):.4f}")

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
            lower[i] = np.clip(median[i] - sd, 0, VMAX)
            upper[i] = np.clip(median[i] + sd, 0, VMAX)
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
            log(f"  {d}s: {tag:<15}  global-mean burnt%={np.nanmean(median[i]):.4f}  "
                f"trend(base->{d}s)={np.nanmean(trend[i]):+.4f}%/dec")

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
                "long_name": "Annual burnt area (percent of grid cell)",
                "units": "%",
                "statistic": "decadal_mean_annual_burnt_area_percent",
                "value_note": ("median = equal-weighted ensemble mean over 3 fire models "
                               "(lpj-guess, lpjml, mc2-usfs) x 4 GCMs of the annual burnt-"
                               "area percent; raw ISIMIP2b values in [0,100] %."),
                "normalization": ("none -- all 3 models report the SAME unit (% burnt area) "
                                  "so they are equal-weighted in raw %; inter-member spread "
                                  "is retained as the CI (model-democracy decision)."),
                "model_notes": ("lpj-guess long_name is mislabeled 'Fire Return Interval' but "
                                "the data are burnt-area % (verified); lpj-guess floors "
                                "reporting at 0.1%; mc2-usfs is a coarse biome model, "
                                "strongly zero-inflated and ~5-7x hotter in the mean."),
                "spatial_smoothing": "none (12-member ensemble is thick enough)",
                "trend_definition": ("baseline-anchored rate: trend[decade] = (median[decade] "
                                     "- median[2020s]) / (elapsed decades), units % per decade. "
                                     "Each decade's panel is the trend FROM the 2020s baseline "
                                     "TO that decade (so 2090s = the full baseline->2090s "
                                     "trend). Built on decadal means (each a 12-member x 10-yr "
                                     "average) and anchored at the baseline, so it is exactly "
                                     "the (decade-2020s) change map / elapsed decades -- "
                                     "spatially coherent by construction. A within-decade slope "
                                     "of the annual series is NOT used (fire is too noisy year-"
                                     "to-year -> spotty sign-flipping field). The 2020s baseline "
                                     "has no elapsed change -> trend 0 (identical across "
                                     "scenarios)."),
                "trend_units": "% decade-1",
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-member standard "
                                  "deviation (across the 12 model x GCM members), clamped "
                                  "to [0,100]. Spread reflects fire-model + GCM uncertainty."),
                "percentile_baseline": (f"{pct_mode}: ranked against the 2020s ensemble-mean "
                                        "land distribution (single-tier -> [1,100]; two-tier "
                                        "reserves 1 for exact-zero cells and maps positives "
                                        "to [2,100])."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_sens": "2005soc/nosoc, co2 (transient CO2, fixed/no land-use)",
                "source_dataset": "ISIMIP2b OutputData/biomes (burntarea, annual)",
                "description": ("Wildfire burnt-area processed to TCFD 6-value-class format "
                                "with shared 2020s baseline; 3-model x 4-GCM ensemble in "
                                "raw %, no normalization, no spatial smoothing."),
            },
        )
        path = out_dir / f"burntarea_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  saved {path}")

    log("\nDone.")


if __name__ == "__main__":
    main()
