#!/usr/bin/env python3
"""Process fish b30cm (large fish biomass >30cm) data into decadal percentile metrics.

This script processes ISIMIP2b marine fisheries b30cm data following process-metrics
skill requirements:
- Temporal scope: 2010s-2090s (data starts 2006)
- Shared 2020s baseline across all scenarios
- Adaptive windowing based on number of models
- Percentile calculation against 2020s reference distribution

The b30cm variable represents biomass density of consumers with asymptotic length >30cm,
capturing large fish species like tuna and salmon.

Implements shared 2020s baseline: all scenarios share identical 2020s values
computed as the average across ALL available models and scenarios.
"""

import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr

# Configuration following process-metrics skill
DEFAULT_DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
MIN_DATA_POINTS = 100
BASELINE_DECADE = 2020
MIN_YEAR = 2006  # b30cm data starts at 2006
MAX_YEAR = 2099

# RCP scenarios (ISIMIP2b)
TARGET_SCENARIOS = ["rcp26", "rcp45", "rcp60", "rcp85"]


def log(msg: str):
    """Print with immediate flush."""
    print(msg, flush=True)


def parse_filename(fpath: Path) -> Dict[str, str]:
    """Parse ISIMIP2b fish filename to extract metadata.

    Example: apecosm_ipsl-cm5a-lr_nobc_rcp26_wo-diaz_no-fishing_no-oa_b30cm_global_monthly_2006_2100.nc4
    """
    parts = fpath.stem.split("_")
    return {
        'model': parts[0],           # apecosm, ecoocean, macroecological, zoomss
        'gcm': parts[1],             # gfdl-esm2m, ipsl-cm5a-lr, cesm1-bgc
        'scenario': parts[3],        # rcp26, rcp45, rcp60, rcp85
        'diazotroph': parts[4],      # w-diaz, wo-diaz
        'variable': parts[7],        # b30cm
        'timestep': parts[9],        # monthly, annual
    }


def calculate_window_size(num_models: int) -> int:
    """Calculate adaptive window size based on number of models."""
    window_years = max(10, int(np.ceil(MIN_DATA_POINTS / num_models)))
    return window_years


def get_decade_window(decade_start: int, window_years: int) -> Tuple[int, int]:
    """Get temporal window for a decade, centered on decade midpoint."""
    if window_years == 10:
        return decade_start, decade_start + 9

    decade_center = decade_start + 4.5
    half_window = window_years / 2
    window_start = int(np.floor(decade_center - half_window + 0.5))
    window_end = window_start + window_years - 1

    # Ensure we don't go before data start
    window_start = max(window_start, MIN_YEAR)

    return window_start, window_end


def load_and_aggregate_file(fpath: Path) -> Tuple[xr.DataArray, Dict]:
    """Load a NetCDF file and aggregate monthly to annual.

    Returns:
        Tuple of (annual DataArray, metadata dict)
    """
    log(f"    Loading {fpath.name}...")

    meta = parse_filename(fpath)

    # Try normal loading first, fall back to decode_times=False for non-standard calendars
    try:
        ds = xr.open_dataset(fpath)
        time_is_decoded = True
    except ValueError as e:
        if "decode" in str(e).lower() or "calendar" in str(e).lower():
            log(f"      Using decode_times=False for non-standard calendar")
            ds = xr.open_dataset(fpath, decode_times=False)
            time_is_decoded = False
        else:
            raise

    # Get b30cm variable
    var_name = 'b30cm' if 'b30cm' in ds.data_vars else list(ds.data_vars)[0]
    da = ds[var_name]

    # Handle time coordinate based on whether it was decoded
    if time_is_decoded:
        years = da.time.dt.year
    else:
        # For annual data with "years since" units, the time values are year offsets
        # Parse start year from file name (e.g., _2006_2100.nc4)
        parts = fpath.stem.split("_")
        start_year = int(parts[-2])
        end_year = int(parts[-1])

        # For annual data, time values are typically 0, 1, 2, ... or actual years
        time_vals = da.time.values
        if np.max(time_vals) < 500:  # Probably offsets from reference year
            # Check the time units attribute
            time_units = ds['time'].attrs.get('units', '')
            if 'years since' in time_units:
                # Extract reference year
                ref_year = int(time_units.split('years since')[1].split('-')[0].strip())
                years = xr.DataArray(ref_year + time_vals.astype(int), dims=['time'])
            else:
                # Assume years from start_year
                years = xr.DataArray(np.arange(start_year, start_year + len(time_vals)), dims=['time'])
        else:
            # Values are actual years
            years = xr.DataArray(time_vals.astype(int), dims=['time'])

    # Filter to relevant years
    mask = (years >= MIN_YEAR) & (years <= MAX_YEAR)
    da_filtered = da.isel(time=mask)
    years_filtered = years.isel(time=mask)

    if len(da_filtered.time) == 0:
        ds.close()
        return None, meta

    # Aggregate monthly to annual mean (or keep annual as-is)
    if meta['timestep'] == 'monthly' and time_is_decoded:
        annual = da_filtered.groupby("time.year").mean(dim="time")
    else:
        # Annual data or non-decoded time - assign year coordinate
        annual = da_filtered.rename({'time': 'year'})
        annual['year'] = years_filtered.values

    ds.close()
    return annual, meta


def collect_shared_2020s_baseline(
    all_files: List[Path],
    window_years: int
) -> Tuple[xr.DataArray, np.ndarray, np.ndarray, np.ndarray]:
    """Collect 2020s data from ALL files and create shared baseline.

    Args:
        all_files: All b30cm NetCDF files
        window_years: Adaptive window size

    Returns:
        Tuple of:
        - Combined 2020s data averaged across all models/scenarios
        - Flattened baseline values for percentile calculation
        - lat coordinates
        - lon coordinates
    """
    baseline_start, baseline_end = get_decade_window(BASELINE_DECADE, window_years)
    log(f"\nCollecting shared 2020s baseline from ALL files (window: {baseline_start}-{baseline_end})...")

    all_2020s_data = []  # List of (model, gcm, scenario, annual_data)
    lat = None
    lon = None

    for fpath in all_files:
        annual, meta = load_and_aggregate_file(fpath)
        if annual is None:
            continue

        # Store coordinates
        if lat is None:
            lat = annual.lat.values
            lon = annual.lon.values

        # Extract 2020s window
        years = annual.year.values
        mask = (years >= baseline_start) & (years <= baseline_end)
        years_in_window = years[mask]

        if len(years_in_window) > 0:
            data_2020s = annual.sel(year=years_in_window)
            all_2020s_data.append({
                'model': meta['model'],
                'gcm': meta['gcm'],
                'scenario': meta['scenario'],
                'data': data_2020s,
            })
            log(f"      {meta['model']}/{meta['gcm']}/{meta['scenario']}: {len(years_in_window)} years")

    if not all_2020s_data:
        log("  ERROR: No 2020s data found!")
        return None, np.array([]), lat, lon

    log(f"\n  Collected 2020s data from {len(all_2020s_data)} model-gcm-scenario combinations")

    # Average across all to create shared baseline
    # First, compute median per combination
    medians = []
    for item in all_2020s_data:
        median_2d = item['data'].median(dim='year').values
        medians.append(median_2d)

    # Stack and average
    stacked = np.stack(medians, axis=0)
    shared_baseline_2d = np.nanmean(stacked, axis=0)

    log(f"  Shared 2020s baseline computed from {len(medians)} combinations")
    log(f"  Value range: {np.nanmin(shared_baseline_2d):.6f} - {np.nanmax(shared_baseline_2d):.6f} g C/m²")

    # Create flattened baseline for percentile calculation (subsample for efficiency)
    baseline_sample = shared_baseline_2d[::5, ::5]  # Subsample spatially
    baseline_flat = baseline_sample.flatten()
    baseline_flat = baseline_flat[~np.isnan(baseline_flat)]
    log(f"  Baseline sample size for percentiles: {len(baseline_flat):,}")

    # Create DataArray for shared baseline
    shared_baseline_da = xr.DataArray(
        shared_baseline_2d,
        dims=['lat', 'lon'],
        coords={'lat': lat, 'lon': lon},
    )

    return shared_baseline_da, baseline_flat, lat, lon


def process_scenario(
    scenario_files: List[Path],
    scenario: str,
    output_dir: Path,
    window_years: int,
    shared_baseline_2d: xr.DataArray,
    baseline_flat: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
):
    """Process a single scenario with shared 2020s baseline.

    Args:
        scenario_files: Files for this scenario
        scenario: Scenario name (rcp26, rcp45, rcp60, rcp85)
        output_dir: Output directory
        window_years: Adaptive window size
        shared_baseline_2d: Shared 2020s baseline (lat, lon)
        baseline_flat: Flattened baseline for percentile calculation
        lat, lon: Coordinate arrays
    """
    log(f"\n{'='*60}")
    log(f"Processing scenario: {scenario} (with SHARED 2020s baseline)")
    log(f"{'='*60}")

    # Load all files for this scenario
    all_models_data = []

    for fpath in scenario_files:
        annual, meta = load_and_aggregate_file(fpath)
        if annual is not None:
            all_models_data.append({
                'model': meta['model'],
                'gcm': meta['gcm'],
                'diazotroph': meta['diazotroph'],
                'data': annual,
            })
            log(f"      {meta['model']}/{meta['gcm']}: years {annual.year.values[0]}-{annual.year.values[-1]}")

    if not all_models_data:
        log(f"  No data loaded for {scenario}")
        return

    log(f"\n  Loaded {len(all_models_data)} model-gcm combinations")

    # Initialize output arrays
    n_decades = len(DEFAULT_DECADES)
    n_lat, n_lon = len(lat), len(lon)

    result = {
        "median": np.full((n_decades, n_lat, n_lon), np.nan, dtype=np.float32),
        "percentile": np.full((n_decades, n_lat, n_lon), np.nan, dtype=np.float32),
        "lower_ci": np.full((n_decades, n_lat, n_lon), np.nan, dtype=np.float32),
        "upper_ci": np.full((n_decades, n_lat, n_lon), np.nan, dtype=np.float32),
        "model_spread": np.full((n_decades, n_lat, n_lon), np.nan, dtype=np.float32),
    }

    # Process each decade
    for d_idx, decade_start in enumerate(DEFAULT_DECADES):
        window_start, window_end = get_decade_window(decade_start, window_years)

        if decade_start == BASELINE_DECADE:
            # Use SHARED baseline for 2020s
            log(f"  Decade {decade_start}s: Using SHARED baseline (identical across scenarios)")

            decade_median = shared_baseline_2d.values
            result["median"][d_idx] = decade_median

            # Percentile (should be ~50 since it's the baseline)
            # For "higher_is_better" variables: invert so high value → low percentile (safe)
            if len(baseline_flat) > 100:
                flat_median = decade_median.flatten()
                percentiles = np.zeros_like(flat_median)
                for i, val in enumerate(flat_median):
                    if np.isnan(val):
                        percentiles[i] = np.nan
                    else:
                        pct = np.sum(baseline_flat <= val) / len(baseline_flat) * 100
                        pct = 100 - pct  # Invert: high biomass → low percentile (safe)
                        percentiles[i] = max(1, min(100, pct))
                result["percentile"][d_idx] = percentiles.reshape(n_lat, n_lon)

            # For CI and spread, use actual model data if available
            decade_model_data = []
            for item in all_models_data:
                years = item['data'].year.values
                mask = (years >= window_start) & (years <= window_end)
                if np.any(mask):
                    window_data = item['data'].sel(year=years[mask])
                    decade_model_data.append(window_data.median(dim='year').values)

            if decade_model_data:
                stacked = np.stack(decade_model_data, axis=0)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result["lower_ci"][d_idx] = np.nanpercentile(stacked, 25, axis=0)
                    result["upper_ci"][d_idx] = np.nanpercentile(stacked, 75, axis=0)
                    result["model_spread"][d_idx] = np.nanstd(stacked, axis=0)

        else:
            # For 2010s and 2030s-2090s: use scenario-specific data
            log(f"  Decade {decade_start}s: Using {scenario}-specific data (window: {window_start}-{window_end})")

            decade_model_data = []

            for item in all_models_data:
                years = item['data'].year.values
                mask = (years >= window_start) & (years <= window_end)
                years_in_window = years[mask]

                if len(years_in_window) > 0:
                    window_data = item['data'].sel(year=years_in_window)
                    model_median = window_data.median(dim='year').values
                    decade_model_data.append(model_median)

            if not decade_model_data:
                log(f"    No data available for {decade_start}s")
                continue

            log(f"    {len(decade_model_data)} models with data")

            # Stack across models
            stacked = np.stack(decade_model_data, axis=0)

            # Ensemble median
            decade_median = np.nanmedian(stacked, axis=0)
            result["median"][d_idx] = decade_median

            # Percentile against 2020s baseline
            # For "higher_is_better" variables: invert so high value → low percentile (safe)
            if len(baseline_flat) > 100:
                flat_median = decade_median.flatten()
                percentiles = np.zeros_like(flat_median)
                for i, val in enumerate(flat_median):
                    if np.isnan(val):
                        percentiles[i] = np.nan
                    else:
                        pct = np.sum(baseline_flat <= val) / len(baseline_flat) * 100
                        pct = 100 - pct  # Invert: high biomass → low percentile (safe)
                        percentiles[i] = max(1, min(100, pct))
                result["percentile"][d_idx] = percentiles.reshape(n_lat, n_lon)

            # Confidence intervals (25th-75th percentile across models)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result["lower_ci"][d_idx] = np.nanpercentile(stacked, 25, axis=0)
                result["upper_ci"][d_idx] = np.nanpercentile(stacked, 75, axis=0)
                result["model_spread"][d_idx] = np.nanstd(stacked, axis=0)

    # Create output dataset
    log("\nCreating output dataset...")

    ds_out = xr.Dataset(
        {
            "median": (["decade", "lat", "lon"], result["median"], {
                "long_name": "Ensemble median large fish biomass (>30cm)",
                "units": "g C m-2",
            }),
            "percentile": (["decade", "lat", "lon"], result["percentile"], {
                "long_name": "Percentile rank against 2020s baseline",
                "units": "percent",
            }),
            "lower_ci": (["decade", "lat", "lon"], result["lower_ci"], {
                "long_name": "25th percentile across models",
                "units": "g C m-2",
            }),
            "upper_ci": (["decade", "lat", "lon"], result["upper_ci"], {
                "long_name": "75th percentile across models",
                "units": "g C m-2",
            }),
            "model_spread": (["decade", "lat", "lon"], result["model_spread"], {
                "long_name": "Standard deviation across models",
                "units": "g C m-2",
            }),
        },
        coords={
            "decade": DEFAULT_DECADES,
            "lat": lat,
            "lon": lon,
        },
        attrs={
            "title": "Large Fish Biomass (>30cm) - Decadal Projections",
            "variable": "b30cm",
            "scenario": scenario,
            "long_name": "Biomass density of consumers with asymptotic length >30cm",
            "units": "g C m-2",
            "window_years": window_years,
            "baseline_decade": str(BASELINE_DECADE),
            "baseline_source": "shared_across_all_scenarios",
            "source": "ISIMIP2b FishMIP marine-fishery_global",
            "n_models": len(all_models_data),
            "models": ", ".join(sorted(set(item['model'] for item in all_models_data))),
            "gcms": ", ".join(sorted(set(item['gcm'] for item in all_models_data))),
            "description": "Processed with shared 2020s baseline. "
                          "2020s values are identical across scenarios. "
                          "b30cm captures large fish species like tuna and salmon.",
            "percentile_direction": "higher_is_better",  # High biomass = low risk
        }
    )

    # Save
    output_path = output_dir / f"b30cm_{scenario}_processed.nc"
    ds_out.to_netcdf(output_path)
    log(f"  Saved: {output_path}")
    log(f"  File size: {output_path.stat().st_size / 1e6:.1f} MB")


def main():
    """Main processing function."""
    # Get paths
    project_root = Path(__file__).parent.parent
    input_dir = project_root / "data" / "raw" / "fish-b30cm"
    output_dir = project_root / "data" / "processed" / "fish-b30cm_b30cm_monthly"

    # Allow command line override
    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    output_dir.mkdir(parents=True, exist_ok=True)

    log("=" * 60)
    log("Processing Fish b30cm (Large Fish Biomass >30cm) Data")
    log("Following process-metrics skill requirements")
    log("=" * 60)
    log(f"\nInput directory: {input_dir}")
    log(f"Output directory: {output_dir}")

    # Find all files
    all_files = sorted(list(input_dir.glob("*.nc4")) + list(input_dir.glob("*.nc")))
    log(f"\nFound {len(all_files)} NetCDF files")

    if not all_files:
        log("ERROR: No NetCDF files found!")
        return

    # Group by scenario
    scenario_files = {scenario: [] for scenario in TARGET_SCENARIOS}
    for f in all_files:
        meta = parse_filename(f)
        if meta['scenario'] in TARGET_SCENARIOS:
            scenario_files[meta['scenario']].append(f)

    for scenario in TARGET_SCENARIOS:
        log(f"  {scenario.upper()}: {len(scenario_files[scenario])} files")

    # Calculate adaptive window based on unique model-gcm combinations
    unique_combos = set()
    for f in all_files:
        meta = parse_filename(f)
        unique_combos.add((meta['model'], meta['gcm']))

    num_models = len(unique_combos)
    window_years = calculate_window_size(num_models)

    log(f"\nUnique model-GCM combinations: {num_models}")
    log(f"Adaptive window size: {window_years} years")

    # Collect shared 2020s baseline from ALL files
    log("\n" + "=" * 60)
    log("COLLECTING SHARED 2020s BASELINE FROM ALL SCENARIOS")
    log("=" * 60)

    shared_baseline_2d, baseline_flat, lat, lon = collect_shared_2020s_baseline(
        all_files, window_years
    )

    if shared_baseline_2d is None:
        log("ERROR: Could not create shared baseline!")
        return

    # Process each scenario with shared baseline
    log("\n" + "=" * 60)
    log("PROCESSING SCENARIOS WITH SHARED 2020s BASELINE")
    log("=" * 60)

    for scenario in TARGET_SCENARIOS:
        if scenario_files[scenario]:
            process_scenario(
                scenario_files[scenario], scenario, output_dir, window_years,
                shared_baseline_2d, baseline_flat, lat, lon
            )

    log("\n" + "=" * 60)
    log("PROCESSING COMPLETE!")
    log("=" * 60)
    log(f"\nOutput files saved to: {output_dir}")
    log("\nVerify shared baseline with:")
    log(f"  python scripts/test_shared_baseline.py {output_dir}")


if __name__ == "__main__":
    main()
