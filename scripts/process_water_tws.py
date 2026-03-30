"""Process TWS (Total Water Storage) monthly ISIMIP3b data into water index format.

Produces a NetCDF file with dimensions (lat, lon, scenario, value_type, decade) where:
  - value_types 0-11: Per-month ensemble means (Jan-Dec) within each decade
  - value_type 12: Annual mean (= mean of monthly means)
  - value_types 13-19: Annual quantiles Q05, Q15, Q25, Q50, Q75, Q85, Q95

Architecture: Latitude-chunked processing to keep memory under 2 GB.
Data source: SSP projection files only (no historical model runs).

Usage:
    python scripts/process_water_tws.py
    python scripts/process_water_tws.py --data-dir /path/to/raw/data
    python scripts/process_water_tws.py --output /path/to/output.nc
"""

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import xarray as xr

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_water_tws import (
    VARIABLE, VARIABLE_LONG_NAME, VARIABLE_UNITS, AGGREGATION_METHOD,
    IMPACT_MODELS, GCMS, SCENARIOS, DECADES,
    MIN_YEAR, MAX_YEAR, N_VALUE_TYPES, N_LAT, N_LON,
    LAT_CHUNK_SIZE, VALUE_TYPE_NAMES,
    NORMALIZE_MODELS, NORM_TARGET_MEAN, NORM_TARGET_SD, NORM_REF_YEARS,
    get_decade_years,
)

# Quantile levels for vt13-19
ANNUAL_QUANTILES = [5, 15, 25, 50, 75, 85, 95]


def log(msg: str):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def _filter_standard_run(nc_files: List[Path]) -> List[Path]:
    """Filter NetCDF files to keep only the standard social forcing variant.

    Preference order:
      1. Files with '2015soc_default' (standard run)
      2. Files with '2015soc-from-histsoc_default' (transition variant)
    Excludes: '1850soc' (pre-industrial) and '2015co2' (fixed CO2 sensitivity).
    """
    standard = [f for f in nc_files if "_2015soc_default_" in f.name]
    if standard:
        return standard

    histsoc = [f for f in nc_files if "_2015soc-from-histsoc_default_" in f.name]
    if histsoc:
        return histsoc

    fallback = [f for f in nc_files
                if "2015soc" in f.name and "1850soc" not in f.name and "2015co2" not in f.name]
    return fallback if fallback else nc_files


def discover_files(data_dir: Path) -> Dict[str, Dict[str, Dict[str, List[Path]]]]:
    """Discover downloaded NetCDF files organized as model/gcm_scenario/*.nc

    Filters to keep only the standard social forcing variant per (model, gcm, scenario).
    Returns: {model: {gcm: {scenario: [paths]}}}
    """
    inventory: Dict[str, Dict[str, Dict[str, List[Path]]]] = {}

    for model_dir in sorted(data_dir.iterdir()) if data_dir.exists() else []:
        if not model_dir.is_dir():
            continue
        model = model_dir.name
        inventory[model] = {}

        for sub in sorted(model_dir.iterdir()):
            if not sub.is_dir():
                continue
            parts = sub.name.split("_", 1)
            if len(parts) != 2:
                continue
            gcm, scenario = parts
            if gcm not in inventory[model]:
                inventory[model][gcm] = {}
            nc_files = sorted(sub.glob("*.nc"))
            nc_files = _filter_standard_run(nc_files)
            if nc_files:
                inventory[model][gcm][scenario] = nc_files

    if not inventory:
        nc_files = sorted(data_dir.glob("**/*.nc"))
        for fpath in nc_files:
            name = fpath.stem.lower()
            model = gcm = scenario = None

            for m in IMPACT_MODELS:
                if m in name:
                    model = m
                    break
            for g in GCMS:
                if g.replace("-", "") in name.replace("-", ""):
                    gcm = g
                    break
            for s in SCENARIOS:
                if s in name:
                    scenario = s
                    break

            if model and gcm and scenario:
                inventory.setdefault(model, {}).setdefault(gcm, {}).setdefault(scenario, []).append(fpath)

        for model in inventory:
            for gcm in inventory[model]:
                for scenario in inventory[model][gcm]:
                    inventory[model][gcm][scenario] = _filter_standard_run(inventory[model][gcm][scenario])

    return inventory


def count_ensemble_members(inventory: Dict) -> int:
    """Count unique (model, gcm) combinations."""
    members = set()
    for model in inventory:
        for gcm in inventory[model]:
            members.add((model, gcm))
    return len(members)


# ---------------------------------------------------------------------------
# Model normalization
# ---------------------------------------------------------------------------

def compute_normalization_stats(
    inventory: Dict,
    cache_path: Optional[Path] = None,
) -> Dict[str, Dict[str, float]]:
    """Compute per-model median and IQR of annual-mean TWS from reference period.

    Streams one year at a time per file and spatially subsamples (every 4th cell)
    to keep memory manageable.

    Returns: {model: {"median": float, "iqr": float, "p25": float, "p75": float, "n": int}}
    """

    # Check cache first
    if cache_path and cache_path.exists():
        log(f"Loading normalization stats from cache: {cache_path}")
        with open(cache_path) as f:
            stats = json.load(f)
        for model, s in stats.items():
            log(f"  {model}: median={s['median']:.1f}, IQR={s['iqr']:.1f}")
        return stats

    ref_start, ref_end = NORM_REF_YEARS
    stride = 4  # spatial subsampling for efficiency

    stats = {}

    for model in sorted(inventory.keys()):
        log(f"  Computing stats for {model}...")
        all_annual_values = []

        for gcm in sorted(inventory[model].keys()):
            for scenario in sorted(inventory[model][gcm].keys()):
                files = inventory[model][gcm][scenario]
                for fpath in files:
                    try:
                        ds = xr.open_dataset(fpath)
                    except Exception:
                        continue

                    var_name = VARIABLE if VARIABLE in ds.data_vars else list(ds.data_vars)[0]
                    da = ds[var_name]

                    try:
                        time_years = da.time.dt.year
                    except Exception:
                        ds.close()
                        continue

                    # Process one year at a time
                    for yr in range(ref_start, ref_end + 1):
                        yr_mask = time_years == yr
                        if yr_mask.values.sum() == 0:
                            continue
                        monthly = da.isel(time=yr_mask).values  # (months, lat, lon)
                        # Subsample spatially
                        monthly_sub = monthly[:, ::stride, ::stride]
                        # Annual mean
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            ann_mean = np.nanmean(monthly_sub, axis=0)  # (lat_sub, lon_sub)
                        valid = ann_mean[~np.isnan(ann_mean)]
                        if len(valid) > 0:
                            all_annual_values.append(valid)

                    ds.close()

        if all_annual_values:
            pooled = np.concatenate(all_annual_values)
            p25 = float(np.percentile(pooled, 25))
            p75 = float(np.percentile(pooled, 75))
            median = float(np.median(pooled))
            iqr = p75 - p25

            stats[model] = {
                "median": median,
                "iqr": iqr,
                "p25": p25,
                "p75": p75,
                "n": int(len(pooled)),
            }
            log(f"    {model}: median={median:.1f}, IQR={iqr:.1f} "
                f"(P25={p25:.1f}, P75={p75:.1f}, n={len(pooled):,})")
        else:
            log(f"    WARNING: no data for {model}")

    # Save cache
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(stats, f, indent=2)
        log(f"  Normalization stats cached: {cache_path}")

    return stats


def apply_normalization(
    data: np.ndarray,
    model_median: float,
    model_iqr: float,
) -> np.ndarray:
    """Apply robust z-score normalization: target_mean + (value - median) / IQR * target_sd."""
    if model_iqr <= 0:
        return data
    return NORM_TARGET_MEAN + (data - model_median) / model_iqr * NORM_TARGET_SD


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_monthly_data(
    files: List[Path], year_start: int, year_end: int,
) -> Optional[xr.DataArray]:
    """Load and concatenate monthly data from multiple files, filtering to year range."""
    arrays = []
    for fpath in files:
        try:
            ds = xr.open_dataset(fpath, use_cftime=True)
        except Exception:
            try:
                ds = xr.open_dataset(fpath)
            except Exception as e:
                log(f"    WARNING: Cannot open {fpath.name}: {e}")
                continue

        var_name = VARIABLE if VARIABLE in ds.data_vars else list(ds.data_vars)[0]
        da = ds[var_name]

        try:
            years = da.time.dt.year
            mask = (years >= year_start) & (years <= year_end)
            da = da.isel(time=mask)
        except Exception:
            pass

        if len(da.time) > 0:
            arrays.append(da)
        ds.close()

    if not arrays:
        return None

    combined = xr.concat(arrays, dim="time")
    combined = combined.sortby("time")
    return combined


# ---------------------------------------------------------------------------
# Core processing: per lat-chunk
# ---------------------------------------------------------------------------

def process_chunk(
    inventory: Dict[str, Dict[str, Dict[str, List[Path]]]],
    lat_start: int,
    lat_end: int,
    lats: np.ndarray,
    lons: np.ndarray,
    norm_stats: Optional[Dict[str, Dict[str, float]]] = None,
) -> np.ndarray:
    """Process a latitude chunk for all scenarios, decades, and value_types.

    If norm_stats is provided, applies per-model robust z-score normalization
    before ensemble averaging.

    Returns array of shape (lat_chunk, n_lon, n_scenarios, n_value_types, n_decades).
    """
    n_lat_chunk = lat_end - lat_start
    n_lon = len(lons)
    n_scenarios = len(SCENARIOS)
    n_decades = len(DECADES)

    result = np.full(
        (n_lat_chunk, n_lon, n_scenarios, N_VALUE_TYPES, n_decades),
        np.nan, dtype=np.float32,
    )

    for s_idx, scenario in enumerate(SCENARIOS):
        log(f"    Scenario {scenario}...")

        # Collect all ensemble members' monthly data for this scenario
        ensemble_monthly = []

        for model in inventory:
            for gcm in inventory[model]:
                if scenario not in inventory[model][gcm]:
                    continue
                files = inventory[model][gcm][scenario]
                da = load_monthly_data(files, MIN_YEAR, MAX_YEAR)
                if da is None:
                    continue

                chunk_lats = lats[lat_start:lat_end]
                try:
                    da_chunk = da.sel(lat=chunk_lats, method="nearest")
                except Exception:
                    da_chunk = da.isel(lat=slice(lat_start, lat_end))

                # Apply per-model normalization if stats available
                if norm_stats and model in norm_stats:
                    ms = norm_stats[model]
                    da_chunk = da_chunk.copy(data=apply_normalization(
                        da_chunk.values, ms["median"], ms["iqr"]))

                ensemble_monthly.append(da_chunk)

        if not ensemble_monthly:
            log(f"      No data for {scenario}")
            continue

        n_members = len(ensemble_monthly)
        log(f"      {n_members} ensemble members loaded")

        for d_idx, decade in enumerate(DECADES):
            year_start, year_end = get_decade_years(decade)

            # --- Value types 0-11: Monthly ensemble means ---
            for month in range(12):
                month_values = []
                for da in ensemble_monthly:
                    try:
                        years_mask = (da.time.dt.year >= year_start) & (da.time.dt.year <= year_end)
                        month_mask = da.time.dt.month == (month + 1)
                        subset = da.isel(time=years_mask & month_mask)
                        if len(subset.time) > 0:
                            month_values.append(subset.values)
                    except Exception:
                        continue

                if month_values:
                    stacked = np.concatenate(month_values, axis=0)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        result[:, :, s_idx, month, d_idx] = np.nanmean(stacked, axis=0)

            # --- Value type 12: Annual mean (= mean of vt0-11) ---
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result[:, :, s_idx, 12, d_idx] = np.nanmean(
                    result[:, :, s_idx, 0:12, d_idx], axis=2)

            # --- Value types 13-19: Annual quantiles Q05..Q95 ---
            # Pool raw annual values across all ensemble members and years
            member_annual_values = []
            for da in ensemble_monthly:
                try:
                    years_mask = (da.time.dt.year >= year_start) & (da.time.dt.year <= year_end)
                    subset = da.isel(time=years_mask)
                    if len(subset.time) == 0:
                        continue

                    if AGGREGATION_METHOD == "mean":
                        annual = subset.groupby("time.year").mean(dim="time")
                    else:
                        annual = subset.groupby("time.year").sum(dim="time")

                    # annual has shape (n_years_this_member, lat_chunk, lon)
                    member_annual_values.append(annual.values)
                except Exception:
                    continue

            if member_annual_values:
                # Pool: (n_members * n_years, lat_chunk, lon)
                pooled = np.concatenate(member_annual_values, axis=0)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    for q_idx, q_level in enumerate(ANNUAL_QUANTILES):
                        result[:, :, s_idx, 13 + q_idx, d_idx] = np.nanpercentile(
                            pooled, q_level, axis=0)

    return result


# ---------------------------------------------------------------------------
# Output writer
# ---------------------------------------------------------------------------

def write_output(
    result: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    output_path: Path,
    norm_stats: Optional[Dict] = None,
):
    """Write the processed result to NetCDF matching the target format."""
    attrs = {
        "Conventions": "CF-1.8",
        "title": f"Water Index Underlying Data - {VARIABLE_LONG_NAME}",
        "variable": VARIABLE,
        "variable_long_name": VARIABLE_LONG_NAME,
        "units": "normalized (robust z-score)" if norm_stats else VARIABLE_UNITS,
        "units_original": VARIABLE_UNITS,
        "source": "ISIMIP3b",
        "scenarios": ", ".join(SCENARIOS),
        "data_source": "SSP projections only (no historical model runs)",
        "aggregation_method": AGGREGATION_METHOD,
        "impact_models": ", ".join(IMPACT_MODELS),
        "gcms": ", ".join(GCMS),
        "value_type_names": str(VALUE_TYPE_NAMES),
    }

    if norm_stats:
        attrs["normalization_method"] = "robust_zscore"
        attrs["normalization_formula"] = (
            f"target_mean + (value - model_median) / model_IQR * target_sd"
        )
        attrs["normalization_target_mean"] = float(NORM_TARGET_MEAN)
        attrs["normalization_target_sd"] = float(NORM_TARGET_SD)
        attrs["normalization_ref_period"] = f"{NORM_REF_YEARS[0]}-{NORM_REF_YEARS[1]}"
        for model, s in norm_stats.items():
            attrs[f"normalization_{model}_median"] = s["median"]
            attrs[f"normalization_{model}_iqr"] = s["iqr"]

    ds = xr.Dataset(
        {
            VARIABLE: (
                ["lat", "lon", "scenario", "value_type", "decade"],
                result,
            ),
        },
        coords={
            "lat": ("lat", lats, {"units": "degrees_north", "long_name": "latitude"}),
            "lon": ("lon", lons, {"units": "degrees_east", "long_name": "longitude"}),
            "scenario": ("scenario", SCENARIOS),
            "value_type": ("value_type", list(range(N_VALUE_TYPES))),
            "decade": ("decade", DECADES),
        },
        attrs=attrs,
    )

    encoding = {
        VARIABLE: {
            "dtype": "float32",
            "zlib": True,
            "complevel": 4,
            "_FillValue": np.float32(np.nan),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output_path, encoding=encoding)
    log(f"Output saved: {output_path}")

    size_mb = output_path.stat().st_size / (1024 * 1024)
    log(f"File size: {size_mb:.1f} MB")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Process TWS monthly ISIMIP3b data into water index format"
    )
    parser.add_argument(
        "--data-dir", type=Path, default=PROJECT_ROOT / "data" / "raw" / "water_tws",
        help="Directory with downloaded raw NetCDF files",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path(r"C:\Cai_data\WaterIndex\waterIndexUnderlyingData_tws_ssp.nc"),
        help="Output NetCDF file path",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=LAT_CHUNK_SIZE,
        help=f"Latitude rows per chunk (default: {LAT_CHUNK_SIZE})",
    )
    parser.add_argument(
        "--norm-cache", type=Path,
        default=PROJECT_ROOT / "data" / "normalization_stats_tws.json",
        help="Path to cache normalization statistics JSON",
    )
    parser.add_argument(
        "--no-normalize", action="store_true",
        help="Disable per-model normalization",
    )
    args = parser.parse_args()

    log("=" * 70)
    log(f"Processing {VARIABLE_LONG_NAME} ({VARIABLE}) - Monthly to Water Index")
    log("=" * 70)

    log(f"\nDiscovering files in {args.data_dir}...")
    inventory = discover_files(args.data_dir)

    if not inventory:
        log("ERROR: No data files found!")
        log(f"Expected directory structure: {args.data_dir}/{{model}}/{{gcm}}_{{scenario}}/*.nc")
        sys.exit(1)

    n_members = count_ensemble_members(inventory)
    log(f"\nInventory:")
    for model in sorted(inventory):
        for gcm in sorted(inventory[model]):
            scenarios = sorted(inventory[model][gcm].keys())
            n_files = sum(len(inventory[model][gcm][s]) for s in scenarios)
            log(f"  {model} / {gcm}: {', '.join(scenarios)} ({n_files} files)")
    log(f"\nTotal ensemble members: {n_members}")

    # --- Model normalization pre-pass ---
    do_normalize = NORMALIZE_MODELS and not args.no_normalize
    norm_stats = None

    if do_normalize:
        log(f"\n--- Model Normalization ---")
        log(f"Reference period: {NORM_REF_YEARS[0]}-{NORM_REF_YEARS[1]}")
        log(f"Target: mean={NORM_TARGET_MEAN}, SD={NORM_TARGET_SD}")
        norm_stats = compute_normalization_stats(inventory, cache_path=args.norm_cache)
        if not norm_stats:
            log("WARNING: No normalization stats computed, proceeding without normalization")
            norm_stats = None
    else:
        log("\nNormalization disabled.")

    lats = np.arange(89.75, -90, -0.5, dtype=np.float32)[:N_LAT]
    lons = np.arange(-179.75, 180, 0.5, dtype=np.float32)[:N_LON]

    full_result = np.full(
        (N_LAT, N_LON, len(SCENARIOS), N_VALUE_TYPES, len(DECADES)),
        np.nan, dtype=np.float32,
    )

    chunk_size = args.chunk_size
    n_chunks = (N_LAT + chunk_size - 1) // chunk_size

    log(f"\nProcessing {N_LAT} latitude rows in {n_chunks} chunks of {chunk_size}...")

    for chunk_idx in range(n_chunks):
        lat_start = chunk_idx * chunk_size
        lat_end = min(lat_start + chunk_size, N_LAT)

        log(f"\n{'='*50}")
        log(f"Chunk {chunk_idx + 1}/{n_chunks}: lats [{lat_start}:{lat_end}] "
            f"({lats[lat_start]:.2f} to {lats[lat_end - 1]:.2f})")
        log(f"{'='*50}")

        chunk_result = process_chunk(
            inventory, lat_start, lat_end, lats, lons, norm_stats,
        )

        full_result[lat_start:lat_end, :, :, :, :] = chunk_result

    log(f"\nWriting output to {args.output}...")
    write_output(full_result, lats, lons, args.output, norm_stats)

    log("\nSummary statistics:")
    for s_idx, scenario in enumerate(SCENARIOS):
        data = full_result[:, :, s_idx, 12, :]
        valid = ~np.isnan(data)
        log(f"  {scenario}: {valid.sum():,} valid cells, "
            f"range [{np.nanmin(data):.4f}, {np.nanmax(data):.4f}]")

    log("\n" + "=" * 70)
    log("Processing complete!")
    log("=" * 70)


if __name__ == "__main__":
    main()
