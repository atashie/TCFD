#!/usr/bin/env python3
"""Process ISIMIP2b temperature-related mortality data with shared 2020s baseline.

This script processes TRM-Tsukuba heat mortality data (an-tot-heat variable),
implementing shared 2020s baseline: all scenarios share identical 2020s values
computed as the average across ALL available scenarios.

Variable: an-tot-heat (Annual total heat-related mortality)
Units: deaths per year (need to verify from file metadata)
"""

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import xarray as xr
from scipy import stats

# Add scripts directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.land_mask import (
    get_isimip_landmask,
    load_landmask,
    apply_land_mask,
    verify_grid_compatibility,
    get_mask_stats,
)


def extract_years(ds: xr.Dataset) -> np.ndarray:
    """Extract integer years from dataset time coordinate."""
    time_values = ds.time.values
    if hasattr(time_values[0], 'year'):
        return np.array([t.year for t in time_values])
    return ds.time.dt.year.values


def collect_2020s_baseline(gcm_data: Dict, window_start: int = 2020, window_end: int = 2030) -> Dict[str, np.ndarray]:
    """Collect 2020s data from ALL scenarios and average per GCM.

    This ensures all climate projections share identical 2020s baseline values,
    computed as the average across all scenarios (projection + historical/picontrol).

    Args:
        gcm_data: Dict mapping GCM name to scenario data
        window_start: Start year for baseline window
        window_end: End year for baseline window (exclusive)

    Returns:
        Dict mapping GCM name to shared 2020s baseline array (lat, lon)
    """
    shared_baselines = {}

    for gcm, scenarios_data in gcm_data.items():
        all_2020s = []

        print(f"\n  Collecting 2020s baseline for {gcm}...")

        for scenario, sdata in scenarios_data.items():
            data = sdata['data']
            years = sdata['years']

            # Extract 2020s window
            mask = (years >= window_start) & (years < window_end)

            if np.any(mask):
                decade_data = data[mask]
                # Compute median for this scenario's 2020s
                scenario_median = np.nanmedian(decade_data, axis=0)
                all_2020s.append(scenario_median)
                valid_vals = scenario_median[~np.isnan(scenario_median)]
                if len(valid_vals) > 0:
                    print(f"    {scenario}: {np.sum(mask)} years (median range: "
                          f"{np.nanmin(valid_vals):.2e} - {np.nanmax(valid_vals):.2e})")

        if all_2020s:
            # Average across all scenarios to create shared baseline
            stacked = np.stack(all_2020s, axis=0)
            shared_baseline = np.nanmean(stacked, axis=0)
            shared_baselines[gcm] = shared_baseline
            valid_vals = shared_baseline[~np.isnan(shared_baseline)]
            if len(valid_vals) > 0:
                print(f"    SHARED 2020s baseline from {len(all_2020s)} scenarios: "
                      f"range {np.nanmin(valid_vals):.2e} - {np.nanmax(valid_vals):.2e}")
        else:
            print(f"    WARNING: No 2020s data for {gcm}")

    return shared_baselines


def calculate_trend(data: np.ndarray, decades: list) -> np.ndarray:
    """Calculate linear trend (slope) across decades for each grid cell.

    Args:
        data: Array of shape (n_decades, lat, lon)
        decades: List of decade start years

    Returns:
        Array of shape (lat, lon) with trend values (units per decade)
    """
    n_decades, n_lat, n_lon = data.shape
    trend = np.full((n_lat, n_lon), np.nan)

    x = np.arange(n_decades)  # 0, 1, 2, ... for decades

    for i in range(n_lat):
        for j in range(n_lon):
            y = data[:, i, j]
            valid = ~np.isnan(y)
            if np.sum(valid) >= 3:  # Need at least 3 points for meaningful trend
                slope, _, _, _, _ = stats.linregress(x[valid], y[valid])
                trend[i, j] = slope

    return trend


def process_health_mortality_files(input_dir: Path, output_dir: Path):
    """Process health mortality files and create decadal statistics with shared 2020s baseline.

    Args:
        input_dir: Directory containing raw NetCDF files
        output_dir: Directory for processed output
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download/load ISIMIP land-sea mask
    print("Loading ISIMIP2b land-sea mask...")
    mask_path = get_isimip_landmask(version="2b")
    land_mask = load_landmask(mask_path)
    mask_stats = get_mask_stats(land_mask)
    print(f"  Land cells: {mask_stats['land_cells']:,} ({mask_stats['land_fraction']:.1%})")
    print(f"  Ocean cells: {mask_stats['ocean_cells']:,} ({mask_stats['ocean_fraction']:.1%})")

    # Find all NetCDF files
    nc_files = list(input_dir.glob("*.nc4")) + list(input_dir.glob("*.nc"))
    print(f"Found {len(nc_files)} NetCDF files")

    if len(nc_files) == 0:
        print("ERROR: No NetCDF files found!")
        return

    # Group files by GCM
    gcm_data = {}
    units = None
    long_name = None

    for f in nc_files:
        # Extract GCM and scenario from filename
        # Pattern: trm-tsukuba_{gcm}_ewembi_{scenario}_2005soc_an-tot-heat-all_global_annual_{start}_{end}.nc4
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm2m
        scenario = parts[3]  # e.g., historical, picontrol, rcp26, rcp60

        if gcm not in gcm_data:
            gcm_data[gcm] = {}

        # Open with decode_times=False to handle non-standard time units ('years since...')
        ds = xr.open_dataset(f, decode_times=False)
        var_name = list(ds.data_vars)[0]

        # Handle non-standard time units manually
        # Time is in 'years since YYYY-M-D', so we add the reference year to the values
        time_units = ds.time.attrs.get('units', '')
        if 'years since' in time_units:
            # Extract reference year from units like 'years since 1661-1-1 00:00:00'
            ref_year = int(time_units.split('since')[1].strip().split('-')[0])
            years = (ds.time.values + ref_year).astype(int)
        else:
            years = extract_years(ds)

        # Capture units and long_name from first file
        if units is None and 'units' in ds[var_name].attrs:
            units = ds[var_name].attrs['units']
        if long_name is None and 'long_name' in ds[var_name].attrs:
            long_name = ds[var_name].attrs['long_name']

        # Verify grid compatibility on first file
        if len(gcm_data) == 0:
            verify_grid_compatibility(ds.lat.values, ds.lon.values, land_mask)
            print("  Grid compatibility verified with land mask")

        # Apply land mask to data (set ocean cells to NaN)
        raw_data = ds[var_name].values
        masked_data = apply_land_mask(raw_data, land_mask)

        # Append to existing data for this GCM/scenario or create new
        if scenario in gcm_data[gcm]:
            # Concatenate with existing data
            existing = gcm_data[gcm][scenario]
            gcm_data[gcm][scenario] = {
                'data': np.concatenate([existing['data'], masked_data], axis=0),
                'years': np.concatenate([existing['years'], years]),
                'lat': ds.lat.values,
                'lon': ds.lon.values,
            }
        else:
            gcm_data[gcm][scenario] = {
                'data': masked_data,  # (time, lat, lon) with ocean masked
                'years': years,
                'lat': ds.lat.values,
                'lon': ds.lon.values,
            }
        print(f"  Loaded {f.name}: {gcm}/{scenario}, years {years[0]}-{years[-1]}")
        ds.close()

    print(f"\nVariable units: {units}")
    print(f"Variable long_name: {long_name}")

    # Discover available scenarios dynamically
    all_scenarios = set()
    for gcm_scenarios in gcm_data.values():
        all_scenarios.update(gcm_scenarios.keys())
    print(f"\nDiscovered scenarios: {sorted(all_scenarios)}")

    # Identify projection scenarios (exclude historical and picontrol)
    projection_scenarios = sorted([s for s in all_scenarios if s not in ['historical', 'picontrol']])
    print(f"Projection scenarios for output: {projection_scenarios}")

    # Step 1: Calculate shared 2020s baseline across ALL scenarios
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_baselines = collect_2020s_baseline(gcm_data)

    # Step 2: Process each scenario using shared baseline for 2020s
    decades = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]

    print("\n" + "=" * 60)
    print("PROCESSING SCENARIOS WITH SHARED BASELINE")
    print("=" * 60)

    # Collect per-GCM results for ensemble averaging
    gcm_results = {scenario: {} for scenario in projection_scenarios}

    for gcm, scenarios_data in gcm_data.items():
        print(f"\nProcessing GCM: {gcm}")

        # Get dimensions from first scenario
        first_scenario = list(scenarios_data.values())[0]
        lat = first_scenario['lat']
        lon = first_scenario['lon']

        # Get shared baseline for this GCM
        shared_2020s = shared_baselines.get(gcm)
        if shared_2020s is None:
            print(f"  WARNING: No shared baseline for {gcm}, skipping...")
            continue

        # Process each projection scenario separately
        for scenario in projection_scenarios:
            if scenario not in scenarios_data:
                print(f"  Skipping {scenario} (not available)")
                continue

            print(f"\n  Processing {scenario}...")

            # Initialize output arrays for this scenario
            n_decades = len(decades)
            median_data = np.full((n_decades, len(lat), len(lon)), np.nan)

            data = scenarios_data[scenario]['data']
            years = scenarios_data[scenario]['years']

            for d_idx, decade_start in enumerate(decades):
                decade_end = decade_start + 10
                mask = (years >= decade_start) & (years < decade_end)

                if decade_start == 2020:
                    # Use SHARED baseline for 2020s (identical across all scenarios)
                    median_data[d_idx] = shared_2020s
                    print(f"    {decade_start}s: SHARED BASELINE (from all scenarios)")
                elif np.any(mask):
                    # Use scenario-specific data for 2030s-2090s
                    decade_data = data[mask]
                    median_data[d_idx] = np.nanmedian(decade_data, axis=0)
                    print(f"    {decade_start}s: {np.sum(mask)} years (scenario-specific)")

            # Calculate trend
            trend_data = calculate_trend(median_data, decades)

            # Store for ensemble averaging
            gcm_results[scenario][gcm] = {
                'median': median_data,
                'trend': trend_data,
                'lat': lat,
                'lon': lon,
            }

            # Create per-GCM output dataset
            output_ds = xr.Dataset(
                data_vars={
                    'mortality_median': (['decade', 'lat', 'lon'], median_data),
                    'mortality_trend': (['lat', 'lon'], trend_data),
                },
                coords={
                    'decade': decades,
                    'lat': lat,
                    'lon': lon,
                },
                attrs={
                    'description': 'Decadal median heat-related mortality',
                    'units': units or 'unknown',
                    'long_name': long_name or 'Annual total heat-related mortality',
                    'variable': 'an-tot-heat',
                    'gcm': gcm,
                    'scenario': scenario,
                    'source': 'ISIMIP2b TRM-Tsukuba model',
                    'baseline_source': 'shared_across_all_scenarios',
                    'baseline_decade': '2020s',
                    'baseline_note': '2020s values are identical across all scenarios, '
                                     'computed as average of all available scenarios',
                    'land_masked': 'true',
                    'mask_source': 'ISIMIP2b_landseamask_generic.nc4',
                }
            )

            # Output per-GCM file
            output_file = output_dir / f"an-tot-heat_{scenario}_{gcm}_decadal.nc"
            output_ds.to_netcdf(output_file)
            print(f"    Saved: {output_file}")
            output_ds.close()

    # Step 3: Create ensemble files (average across GCMs)
    print("\n" + "=" * 60)
    print("CREATING ENSEMBLE FILES (MULTI-GCM AVERAGE)")
    print("=" * 60)

    for scenario in projection_scenarios:
        gcm_medians = []
        gcm_trends = []

        for gcm, results in gcm_results[scenario].items():
            gcm_medians.append(results['median'])
            gcm_trends.append(results['trend'])
            lat = results['lat']
            lon = results['lon']

        if len(gcm_medians) == 0:
            print(f"  No GCM data for {scenario}, skipping ensemble...")
            continue

        # Stack and compute ensemble statistics
        stacked_medians = np.stack(gcm_medians, axis=0)  # (n_gcm, n_decades, lat, lon)
        stacked_trends = np.stack(gcm_trends, axis=0)  # (n_gcm, lat, lon)

        ensemble_median = np.nanmean(stacked_medians, axis=0)
        ensemble_trend = np.nanmean(stacked_trends, axis=0)

        # GCM spread (lower and upper bounds)
        lower_ci = np.nanmin(stacked_medians, axis=0)
        upper_ci = np.nanmax(stacked_medians, axis=0)

        # Calculate percentile rank against 2020s distribution
        # Get 2020s values (first decade)
        baseline_2020s = ensemble_median[0]  # (lat, lon)

        # Calculate percentile for each decade
        percentile_data = np.full_like(ensemble_median, np.nan)
        baseline_flat = baseline_2020s.flatten()
        baseline_valid = baseline_flat[~np.isnan(baseline_flat)]

        if len(baseline_valid) > 0:
            for d_idx in range(len(decades)):
                decade_values = ensemble_median[d_idx]
                for i in range(len(lat)):
                    for j in range(len(lon)):
                        val = decade_values[i, j]
                        if not np.isnan(val):
                            percentile_data[d_idx, i, j] = stats.percentileofscore(baseline_valid, val)

        # Create ensemble output dataset
        ensemble_ds = xr.Dataset(
            data_vars={
                'median': (['decade', 'lat', 'lon'], ensemble_median),
                'percentile': (['decade', 'lat', 'lon'], percentile_data),
                'trend': (['lat', 'lon'], ensemble_trend),
                'lower_ci': (['decade', 'lat', 'lon'], lower_ci),
                'upper_ci': (['decade', 'lat', 'lon'], upper_ci),
            },
            coords={
                'decade': decades,
                'lat': lat,
                'lon': lon,
            },
            attrs={
                'description': 'Decadal heat-related mortality ensemble statistics',
                'units': units or 'unknown',
                'long_name': long_name or 'Annual total heat-related mortality',
                'variable': 'an-tot-heat',
                'scenario': scenario,
                'source': 'ISIMIP2b TRM-Tsukuba model',
                'ensemble_members': list(gcm_results[scenario].keys()),
                'n_gcms': len(gcm_results[scenario]),
                'baseline_source': 'shared_across_all_scenarios',
                'baseline_decade': '2020s',
                'land_masked': 'true',
                'mask_source': 'ISIMIP2b_landseamask_generic.nc4',
            }
        )

        output_file = output_dir / f"an-tot-heat_{scenario}_processed.nc"
        ensemble_ds.to_netcdf(output_file)
        print(f"  Saved ensemble: {output_file} ({len(gcm_results[scenario])} GCMs)")
        ensemble_ds.close()

    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE!")
    print("=" * 60)
    print("\nVerify shared baseline with:")
    print(f"  python scripts/test_shared_baseline.py {output_dir}")
    print("\nGenerate QA report with:")
    print(f"  python scripts/generate_maps.py an-tot-heat {output_dir} reports")


if __name__ == "__main__":
    input_dir = Path("data/raw/health-mortality_an-tot-heat_annual")
    output_dir = Path("data/processed/health-mortality_an-tot-heat_annual")

    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    process_health_mortality_files(input_dir, output_dir)
