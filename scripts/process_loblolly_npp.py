#!/usr/bin/env python3
"""Process CLM45 npp data for temperate needleleaf trees (loblolly pine proxy).

This script processes ISIMIP2b CLM45 npp data for needleleaf-evergreen-tree-temperate PFT.
Handles incomplete GCM/scenario coverage and implements shared 2020s baseline.

Data characteristics:
- Model: CLM45 (ISIMIP2b)
- Variable: npp (net primary productivity)
- PFT: needleleaf-evergreen-tree-temperate
- Period: 2006-2099
- Scenarios: RCP (rcp26, rcp60, rcp85) - incomplete coverage
"""

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import xarray as xr
from scipy import stats


def extract_years(ds: xr.Dataset) -> np.ndarray:
    """Extract integer years from dataset time coordinate.

    Handles non-standard calendars (365_day, noleap) used by CLM45.
    """
    time_values = ds.time.values

    # If time is already decoded as cftime objects
    if hasattr(time_values[0], 'year'):
        return np.array([t.year for t in time_values])

    # If time values are numeric (years since reference)
    # Check if we can use dt accessor
    try:
        return ds.time.dt.year.values
    except AttributeError:
        pass

    # If time is numeric, check the units attribute
    if 'units' in ds.time.attrs:
        units = ds.time.attrs['units']
        if 'years since' in units:
            # Parse reference year and add time values
            ref_year = int(units.split('since')[1].strip().split('-')[0])
            return (ref_year + time_values).astype(int)

    # Fallback: assume time values are years directly
    return time_values.astype(int)


def calculate_trend(data: np.ndarray, decades: list) -> tuple:
    """Calculate linear trend across decades using Theil-Sen estimator.

    Args:
        data: Array of shape (n_decades, lat, lon)
        decades: List of decade start years

    Returns:
        Tuple of (slope, p_value) arrays of shape (lat, lon)
    """
    n_decades, n_lat, n_lon = data.shape
    slope = np.full((n_lat, n_lon), np.nan)
    p_value = np.full((n_lat, n_lon), np.nan)

    x = np.array(decades)

    for i in range(n_lat):
        for j in range(n_lon):
            y = data[:, i, j]
            valid = ~np.isnan(y)
            if np.sum(valid) >= 3:  # Need at least 3 points for trend
                try:
                    result = stats.theilslopes(y[valid], x[valid])
                    slope[i, j] = result.slope
                    # Approximate p-value using Mann-Kendall
                    _, p = stats.kendalltau(x[valid], y[valid])
                    p_value[i, j] = p
                except:
                    pass

    return slope, p_value


def collect_2020s_baseline(file_data: Dict) -> np.ndarray:
    """Collect 2020s data from ALL files and compute shared baseline.

    Args:
        file_data: Dict mapping filename to data dict

    Returns:
        Shared 2020s baseline array (lat, lon)
    """
    all_2020s = []

    print("\n  Collecting 2020s baseline from all files...")

    for fname, fdata in file_data.items():
        data = fdata['data']
        years = fdata['years']

        # Extract 2020s window (2020-2029)
        mask = (years >= 2020) & (years < 2030)

        if np.any(mask):
            decade_data = data[mask]
            # Compute median for this file's 2020s
            file_median = np.nanmedian(decade_data, axis=0)
            all_2020s.append(file_median)
            print(f"    {fname}: {np.sum(mask)} years")

    if all_2020s:
        # Average across all files to create shared baseline
        stacked = np.stack(all_2020s, axis=0)
        shared_baseline = np.nanmean(stacked, axis=0)
        print(f"  SHARED 2020s baseline from {len(all_2020s)} files")
        print(f"    Range: {np.nanmin(shared_baseline):.6f} - {np.nanmax(shared_baseline):.6f} kg C m-2 s-1")
        return shared_baseline
    else:
        raise ValueError("No 2020s data found in any file!")


def process_npp_files(input_dir: Path, output_dir: Path):
    """Process npp files and create decadal statistics with shared 2020s baseline.

    Args:
        input_dir: Directory containing raw NetCDF files
        output_dir: Directory for processed output
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all npp files (nc4 extension for ISIMIP2b)
    nc_files = list(input_dir.glob("*.nc4")) + list(input_dir.glob("*.nc"))
    print(f"Found {len(nc_files)} NetCDF files")

    if not nc_files:
        print("ERROR: No NetCDF files found!")
        return

    # Load all files and extract metadata
    file_data = {}
    scenarios_available = set()
    gcms_available = set()

    for f in nc_files:
        # Parse CLM45 filename: clm45_{gcm}_ewembi_{scenario}_2005soc_co2_npp-..._global_annual_{start}_{end}.nc4
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm2m
        scenario = parts[3]  # e.g., rcp26, rcp60, rcp85

        ds = xr.open_dataset(f, decode_times=False)
        var_name = list(ds.data_vars)[0]
        years = extract_years(ds)

        file_data[f.name] = {
            'data': ds[var_name].values,  # (time, lat, lon)
            'years': years,
            'lat': ds.lat.values,
            'lon': ds.lon.values,
            'gcm': gcm,
            'scenario': scenario,
        }

        scenarios_available.add(scenario)
        gcms_available.add(gcm)

        print(f"  Loaded {f.name}: {gcm} {scenario}, years {years[0]}-{years[-1]}")
        ds.close()

    print(f"\nScenarios available: {sorted(scenarios_available)}")
    print(f"GCMs available: {sorted(gcms_available)}")

    # Get dimensions from first file
    first_file = list(file_data.values())[0]
    lat = first_file['lat']
    lon = first_file['lon']

    # Step 1: Calculate shared 2020s baseline
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_2020s = collect_2020s_baseline(file_data)

    # Step 2: Process by scenario (ensemble across available GCMs)
    decades = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]

    print("\n" + "=" * 60)
    print("PROCESSING SCENARIOS WITH SHARED BASELINE")
    print("=" * 60)

    for scenario in sorted(scenarios_available):
        print(f"\nProcessing scenario: {scenario}")

        # Find all files for this scenario
        scenario_files = {k: v for k, v in file_data.items() if v['scenario'] == scenario}
        n_gcms = len(scenario_files)
        print(f"  Available GCMs: {n_gcms}")

        if n_gcms == 0:
            continue

        # Initialize arrays for ensemble statistics
        n_decades = len(decades)
        all_gcm_medians = []

        # Process each GCM
        for fname, fdata in scenario_files.items():
            gcm = fdata['gcm']
            data = fdata['data']
            years = fdata['years']

            # Initialize output for this GCM
            gcm_output = np.full((n_decades, len(lat), len(lon)), np.nan)

            for d_idx, decade_start in enumerate(decades):
                decade_end = decade_start + 10
                mask = (years >= decade_start) & (years < decade_end)

                if decade_start == 2020:
                    # Use SHARED baseline for 2020s
                    gcm_output[d_idx] = shared_2020s
                elif np.any(mask):
                    decade_data = data[mask]
                    gcm_output[d_idx] = np.nanmedian(decade_data, axis=0)

            all_gcm_medians.append(gcm_output)

            # Save per-GCM file
            gcm_ds = xr.Dataset(
                data_vars={
                    'npp_median': (['decade', 'lat', 'lon'], gcm_output)
                },
                coords={
                    'decade': decades,
                    'lat': lat,
                    'lon': lon,
                },
                attrs={
                    'description': 'Decadal median net primary productivity for temperate needleleaf trees',
                    'units': 'kg C m-2 s-1',
                    'variable': 'npp-needleleaf-evergreen-tree-temperate',
                    'gcm': gcm,
                    'scenario': scenario,
                    'source': 'ISIMIP2b CLM45 model',
                    'pft': 'needleleaf-evergreen-tree-temperate (temperate conifers)',
                    'baseline_source': 'shared_across_all_scenarios',
                }
            )

            gcm_file = output_dir / f"npp_needleleaf-evergreen-tree-temperate_{scenario}_{gcm}_decadal.nc"
            gcm_ds.to_netcdf(gcm_file)
            print(f"    Saved: {gcm_file.name}")
            gcm_ds.close()

        # Create ensemble file if multiple GCMs
        if n_gcms >= 1:
            stacked = np.stack(all_gcm_medians, axis=0)  # (n_gcms, n_decades, lat, lon)

            # Ensemble statistics
            ensemble_median = np.nanmean(stacked, axis=0)  # Mean across GCMs
            ensemble_lower = np.nanmin(stacked, axis=0)    # GCM spread lower
            ensemble_upper = np.nanmax(stacked, axis=0)    # GCM spread upper

            # Calculate trend
            trend_slope, trend_pvalue = calculate_trend(ensemble_median, decades)

            # Calculate percentile rank against 2020s global distribution
            baseline_2020s = ensemble_median[1]  # decade index 1 = 2020s
            baseline_flat = baseline_2020s.flatten()
            baseline_valid = baseline_flat[~np.isnan(baseline_flat)]

            percentile = np.full_like(ensemble_median, np.nan)
            if len(baseline_valid) > 0:
                for d_idx in range(n_decades):
                    decade_values = ensemble_median[d_idx]
                    for i in range(len(lat)):
                        for j in range(len(lon)):
                            val = decade_values[i, j]
                            if not np.isnan(val):
                                # Rank this cell's value against 2020s global distribution
                                pct = stats.percentileofscore(baseline_valid, val)
                                # Invert: high value → low percentile (safe) for "higher_is_better"
                                percentile[d_idx, i, j] = 100 - pct

            ensemble_ds = xr.Dataset(
                data_vars={
                    'median': (['decade', 'lat', 'lon'], ensemble_median),
                    'percentile': (['decade', 'lat', 'lon'], percentile),
                    'trend': (['lat', 'lon'], trend_slope),
                    'trend_pvalue': (['lat', 'lon'], trend_pvalue),
                    'lower_ci': (['decade', 'lat', 'lon'], ensemble_lower),
                    'upper_ci': (['decade', 'lat', 'lon'], ensemble_upper),
                },
                coords={
                    'decade': decades,
                    'lat': lat,
                    'lon': lon,
                },
                attrs={
                    'description': 'Ensemble decadal net primary productivity for temperate needleleaf trees',
                    'units': 'kg C m-2 s-1',
                    'variable': 'npp-needleleaf-evergreen-tree-temperate',
                    'scenario': scenario,
                    'source': 'ISIMIP2b CLM45 model ensemble',
                    'pft': 'needleleaf-evergreen-tree-temperate (loblolly pine proxy)',
                    'n_gcms': n_gcms,
                    'gcms': ', '.join(sorted([f['gcm'] for f in scenario_files.values()])),
                    'baseline_source': 'shared_across_all_scenarios',
                    'baseline_decade': '2020s',
                    'trend_method': 'Theil-Sen estimator (kg C m-2 s-1 per decade)',
                    'percentile_direction': 'higher_is_better',  # More productivity is good
                    'percentile_note': 'Normalized to 2020s baseline (50 = same as 2020s)',
                }
            )

            ensemble_file = output_dir / f"npp-needleleaf-evergreen-tree-temperate_{scenario}_processed.nc"
            ensemble_ds.to_netcdf(ensemble_file)
            print(f"  Ensemble saved: {ensemble_file.name} ({n_gcms} GCMs)")
            ensemble_ds.close()

    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print("\nFiles created:")
    for f in sorted(output_dir.glob("*.nc")):
        print(f"  {f.name}")


if __name__ == "__main__":
    input_dir = Path("data/raw/loblolly-temperate_npp-needleleaf-evergreen-tree-temperate_annual")
    output_dir = Path("data/processed/loblolly-temperate_npp-needleleaf-evergreen-tree-temperate_annual")

    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    process_npp_files(input_dir, output_dir)
