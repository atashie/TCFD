#!/usr/bin/env python3
"""Process timber cwood data with cftime handling and shared 2020s baseline.

This script processes ISIMIP3b cwood data for deciduous broadleaf trees,
handling the cftime.DatetimeNoLeap calendar used by CLASSIC model.

Implements shared 2020s baseline: all scenarios share identical 2020s values
computed as the average across ALL available scenarios.
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr


def extract_years(ds: xr.Dataset) -> np.ndarray:
    """Extract integer years from dataset time coordinate."""
    time_values = ds.time.values
    if hasattr(time_values[0], 'year'):
        return np.array([t.year for t in time_values])
    return ds.time.dt.year.values


def collect_2020s_baseline(gcm_data: Dict) -> Dict[str, np.ndarray]:
    """Collect 2020s data from ALL scenarios and average per GCM.

    This ensures all climate projections share identical 2020s baseline values,
    computed as the average across all scenarios (projection + historical/picontrol).

    Args:
        gcm_data: Dict mapping GCM name to scenario data

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

            # Extract 2020s window (2020-2029)
            mask = (years >= 2020) & (years < 2030)

            if np.any(mask):
                decade_data = data[mask]
                # Compute median for this scenario's 2020s
                scenario_median = np.nanmedian(decade_data, axis=0)
                all_2020s.append(scenario_median)
                print(f"    {scenario}: {np.sum(mask)} years (median range: "
                      f"{np.nanmin(scenario_median):.2f} - {np.nanmax(scenario_median):.2f})")

        if all_2020s:
            # Average across all scenarios to create shared baseline
            stacked = np.stack(all_2020s, axis=0)
            shared_baseline = np.nanmean(stacked, axis=0)
            shared_baselines[gcm] = shared_baseline
            print(f"    SHARED 2020s baseline from {len(all_2020s)} scenarios: "
                  f"range {np.nanmin(shared_baseline):.2f} - {np.nanmax(shared_baseline):.2f}")
        else:
            print(f"    WARNING: No 2020s data for {gcm}")

    return shared_baselines


def process_cwood_files(input_dir: Path, output_dir: Path):
    """Process cwood files and create decadal statistics with shared 2020s baseline.

    Implements shared baseline: all scenarios share identical 2020s values
    computed as the average across ALL available scenarios.

    Args:
        input_dir: Directory containing raw NetCDF files
        output_dir: Directory for processed output
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all cwood files
    nc_files = list(input_dir.glob("*.nc"))
    print(f"Found {len(nc_files)} NetCDF files")

    # Group files by GCM
    gcm_data = {}
    for f in nc_files:
        # Extract scenario from filename
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm4
        scenario = parts[3]  # e.g., ssp126, ssp370, ssp585, historical

        if gcm not in gcm_data:
            gcm_data[gcm] = {}

        ds = xr.open_dataset(f)
        var_name = list(ds.data_vars)[0]
        years = extract_years(ds)

        gcm_data[gcm][scenario] = {
            'data': ds[var_name].values,  # (time, lat, lon)
            'years': years,
            'lat': ds.lat.values,
            'lon': ds.lon.values,
        }
        print(f"  Loaded {f.name}: {scenario}, years {years[0]}-{years[-1]}")
        ds.close()

    # Step 1: Calculate shared 2020s baseline across ALL scenarios
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_baselines = collect_2020s_baseline(gcm_data)

    # Step 2: Process each scenario using shared baseline for 2020s
    decades = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
    # Only process projection scenarios (not historical/picontrol)
    projection_scenarios = ['ssp126', 'ssp370', 'ssp585']

    print("\n" + "=" * 60)
    print("PROCESSING SCENARIOS WITH SHARED BASELINE")
    print("=" * 60)

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

            # Initialize output array for this scenario
            n_decades = len(decades)
            output_data = np.full((n_decades, len(lat), len(lon)), np.nan)

            data = scenarios_data[scenario]['data']
            years = scenarios_data[scenario]['years']

            for d_idx, decade_start in enumerate(decades):
                decade_end = decade_start + 10
                mask = (years >= decade_start) & (years < decade_end)

                if decade_start == 2020:
                    # Use SHARED baseline for 2020s (identical across all scenarios)
                    output_data[d_idx] = shared_2020s
                    print(f"    {decade_start}s: SHARED BASELINE (from all scenarios)")
                elif np.any(mask):
                    # Use scenario-specific data for 2030s-2090s
                    decade_data = data[mask]
                    output_data[d_idx] = np.nanmedian(decade_data, axis=0)
                    print(f"    {decade_start}s: {np.sum(mask)} years (scenario-specific)")

            # Create output dataset for this scenario
            output_ds = xr.Dataset(
                data_vars={
                    'cwood_median': (['decade', 'lat', 'lon'], output_data)
                },
                coords={
                    'decade': decades,
                    'lat': lat,
                    'lon': lon,
                },
                attrs={
                    'description': 'Decadal median wood carbon for deciduous broadleaf trees',
                    'units': 'kg m-2',
                    'variable': 'cwood-dcdcldbdltr',
                    'gcm': gcm,
                    'scenario': scenario,
                    'source': 'ISIMIP3b CLASSIC model',
                    'pft': 'dcdcldbdltr (deciduous cold-climate broadleaf trees)',
                    'baseline_source': 'shared_across_all_scenarios',
                    'baseline_decade': '2020s',
                    'baseline_note': '2020s values are identical across all scenarios, '
                                     'computed as average of all available scenarios',
                }
            )

            # Output per-scenario file
            output_file = output_dir / f"cwood_dcdcldbdltr_{scenario}_{gcm}_decadal.nc"
            output_ds.to_netcdf(output_file)
            print(f"    Saved: {output_file}")
            output_ds.close()

    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE!")
    print("=" * 60)
    print("\nVerify shared baseline with:")
    print(f"  python scripts/test_shared_baseline.py {output_dir}")


if __name__ == "__main__":
    input_dir = Path("data/raw/timber-cwood-future")
    output_dir = Path("data/processed/timber-cwood_cwood_annual")

    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    process_cwood_files(input_dir, output_dir)
