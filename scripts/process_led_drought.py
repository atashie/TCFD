"""Process led (drought exposure) into the TCFD 6-value-class annualized format.

led (Lange 2020, ISIMIP2b DerivedOutputData) is a BINARY per-cell annual flag:
1 = grid cell exposed to drought that year (monthly soil moisture below the
2.5th percentile of the preindustrial baseline), 0 = not. The "fraction of land
area exposed" emerges only on aggregation.

Consequences handled here (differ from the continuous-variable reference,
process_qg.py):
  * Decadal statistic is the MEAN over years/members (= drought frequency),
    NOT the median (median of a binary series collapses to 0/1).
  * Confidence bounds are the 25th/75th percentile of per-member decadal
    frequencies (spread across the 8-model x 4-GCM ensemble).
  * Time is "years since 1661" on a 360_day calendar and cannot be decoded by
    xarray; years are parsed manually.
  * Ensemble member = (impact_model, GCM), parsed from filename fields [1],[2]
    (lange2020_{model}_{gcm}_ewembi_{scenario}_..._led_global_annual_YYYY_YYYY).

Implements the standard TCFD pipeline otherwise: shared 2020s baseline across
scenarios, percentile-of-score against the 2020s distribution (higher = worse
for drought), and a within-decade trend slope. Output files:
  led_{scenario}_processed.nc  with variables
  {median, percentile, trend, lower_ci, upper_ci} on (decade, lat, lon).
"""

import os
import re
import sys
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402

VAR = "led"
LAYER_ID = "drought_led_annual"
RAW_PATTERN = "*_2006_2099.nc4"
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10  # 62 members >> 10, so standard 10-year decade windows
MIN_YEAR, MAX_YEAR = 2010, 2099


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (impact_model, gcm, scenario, member) from a lange2020 filename."""
    p = os.path.basename(fpath).split("_")
    return dict(model=p[1], gcm=p[2], scenario=p[4], member=f"{p[1]}_{p[2]}")


def years_from_time(ds):
    """Parse integer calendar years from a 'years since YYYY-...' time axis."""
    units = ds.time.attrs.get("units", "")
    m = re.search(r"years since (\d+)", units)
    base = int(m.group(1)) if m else 0
    return (base + ds.time.values).astype(int)


def load_member(fpath):
    """Load one member file as (year, lat, lon), restricted to MIN..MAX_YEAR."""
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    return da


def decade_freq_map(da, decade):
    """Mean over a decade window = per-member drought frequency map (lat, lon)."""
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None
    return da.isel(year=sel).mean("year").values.astype(np.float32)


def slope_over_years(annual_stack, years):
    """OLS slope per cell across `years` for a (nyear, lat, lon) stack."""
    if len(years) < 3:
        return np.full(annual_stack.shape[1:], np.nan, np.float32)
    t = years - years.mean()
    var_t = float(np.mean(t ** 2))
    if var_t <= 0:
        return np.full(annual_stack.shape[1:], np.nan, np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yc = annual_stack - np.nanmean(annual_stack, axis=0, keepdims=True)
        cov = np.nanmean(t[:, None, None] * yc, axis=0)
    return (cov / var_t).astype(np.float32)


def make_pct_fn(baseline_flat):
    """percentile-of-score against the 2020s baseline (higher value = worse)."""
    bsort = np.sort(baseline_flat)
    n = len(bsort)

    def pct(arr):
        flat = arr.ravel()
        out = np.full(flat.shape, np.nan, np.float32)
        fin = np.isfinite(flat)
        ranks = np.searchsorted(bsort, flat[fin], side="right")
        out[fin] = np.clip(ranks / n * 100.0, 1, 100)
        return out.reshape(arr.shape)

    return pct


def main():
    files = [str(p) for p in storage.stage_raw(LAYER_ID, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no led member files found in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
        log("Upload members with storage.ingest_raw() first.")
        return
    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    log("=" * 64)
    log("Processing led (drought exposure) -> TCFD 6-value-class format")
    log("=" * 64)
    log(f"Members: {len(files)} | scenarios discovered: {scenarios}")

    # grid geometry from the first file
    da0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0

    # ---- Pass 1: load every member, compute per-member decadal freq maps -----
    # dec_freq[scenario][member] = (n_decades, lat, lon)
    dec_freq = {s: {} for s in scenarios}
    # ensemble annual accumulators for trend: sum / count over members per year
    years_all = np.arange(MIN_YEAR, MAX_YEAR + 1)
    yidx = {y: i for i, y in enumerate(years_all)}
    ens_sum = {s: np.zeros((len(years_all), LAT, LON), np.float64) for s in scenarios}
    ens_cnt = {s: np.zeros((len(years_all), LAT, LON), np.float64) for s in scenarios}

    for f in files:
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
        # per-member decadal frequency maps
        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        for i, d in enumerate(DECADES):
            fm = decade_freq_map(da, d)
            if fm is not None:
                maps[i] = fm
        dec_freq[s][member] = maps
        # accumulate ensemble annual (for trend)
        v = da.values.astype(np.float64)  # (year, lat, lon), NaN over ocean
        fin = np.isfinite(v)
        rows = np.array([yidx[y] for y in da.year.values])
        ens_sum[s][rows] += np.where(fin, v, 0.0)
        ens_cnt[s][rows] += fin
        log(f"  loaded {info['model']:<11} {info['gcm']:<13} {s}")
        del da, v, fin

    ens_annual = {}
    for s in scenarios:
        with np.errstate(invalid="ignore", divide="ignore"):
            ea = ens_sum[s] / ens_cnt[s]
        ea[ens_cnt[s] == 0] = np.nan
        ens_annual[s] = ea.astype(np.float32)

    # ---- Shared 2020s baseline across scenarios --------------------------------
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec_freq[s]})
    shared_2020 = []
    for member in all_members:
        per_scen = [dec_freq[s][member][b_idx] for s in scenarios if member in dec_freq[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_2020.append(np.nanmean(np.stack(per_scen, 0), axis=0))
    shared_2020 = np.stack(shared_2020, 0).astype(np.float32)  # (member, lat, lon)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared_median = np.nanmean(shared_2020, axis=0)
        shared_sd = np.nanstd(shared_2020, axis=0)
    shared_lo = np.clip(shared_median - shared_sd, 0, 1)
    shared_hi = np.clip(shared_median + shared_sd, 0, 1)

    # Percentile baseline = 2020s ensemble-mean spatial distribution (the same
    # statistic we score), so the 2020s decade centres near the 50th percentile.
    # Scoring against the per-member (binary-inflated) distribution would push
    # every cell into the top quantiles for this rare-event field.
    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"baseline sample n={len(baseline_flat):,}, "
        f"global-mean freq={np.nanmean(shared_median):.4f}")

    # shared 2020s trend: pool both scenarios' 2020s window annual frequency
    win = np.arange(BASELINE_DECADE, BASELINE_DECADE + WINDOW_YEARS)
    wrows = np.array([yidx[y] for y in win])
    pooled_sum = sum(ens_sum[s][wrows] for s in scenarios)
    pooled_cnt = sum(ens_cnt[s][wrows] for s in scenarios)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_annual = pooled_sum / pooled_cnt
    pooled_annual[pooled_cnt == 0] = np.nan
    shared_trend = slope_over_years(pooled_annual.astype(np.float32), win)

    # ---- Per-scenario assembly + write ----------------------------------------
    for s in scenarios:
        log(f"\n{'='*64}\nAssembling scenario {s}\n{'='*64}")
        members = sorted(dec_freq[s])
        stack = np.stack([dec_freq[s][m] for m in members], 0)  # (member, dec, lat, lon)

        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        percentile = np.full_like(median, np.nan)
        trend = np.full_like(median, np.nan)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                median[i], percentile[i] = shared_median, shared_pct
                lower[i], upper[i], trend[i] = shared_lo, shared_hi, shared_trend
                log(f"  {d}s: shared baseline (identical across scenarios)")
                continue
            layer = stack[:, i, :, :]  # (member, lat, lon)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(layer, axis=0)
                sd = np.nanstd(layer, axis=0)
            lower[i] = np.clip(median[i] - sd, 0, 1)
            upper[i] = np.clip(median[i] + sd, 0, 1)
            percentile[i] = pct(median[i])
            win_d = np.arange(d, d + WINDOW_YEARS)
            rows = np.array([yidx[y] for y in win_d if y in yidx])
            trend[i] = slope_over_years(ens_annual[s][rows], win_d[: len(rows)])
            log(f"  {d}s: {layer.shape[0]} members  "
                f"global-mean freq={np.nanmean(median[i]):.4f}")

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
                "long_name": "Land area fraction exposed to drought",
                "units": "1",
                "statistic": "decadal_mean_exposure_frequency",
                "value_note": ("median = ensemble-mean fraction of model-years in "
                               "drought (built from binary per-cell annual flags)"),
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-model "
                                  "standard deviation, clamped to [0,1]"),
                "percentile_baseline": ("percentile-of-score of the decadal ensemble-"
                                        "mean vs the 2020s ensemble-mean spatial "
                                        "distribution (2020s centres near 50)"),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "impact_models": ",".join(sorted({meta[f]["model"] for f in files})),
                "gcms": ",".join(sorted({meta[f]["gcm"] for f in files})),
                "source_dataset": "lange2020 ISIMIP2b DerivedOutputData",
                "description": ("Drought exposure processed to TCFD 6-value-class "
                                "format with shared 2020s baseline."),
            },
        )
        path = out_dir / f"led_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID,
        stage,
        created_by="scripts/process_led_drought.py",
        notes="See WORKFLOW-ISSUES.md 2026-07-24; GUARDRAILS §8-§10.",
    )

    # Every ingest-and-process run leaves reviewable HTML evidence behind.
    # Never raises: the data is already published and gated, and these
    # artifacts are regenerable (see scripts/utils/finalize.py).
    log("\nGenerating QA/QC report and maps...")
    finalize_layer(LAYER_ID, version=version)

    log("\nDone.")


if __name__ == "__main__":
    main()
