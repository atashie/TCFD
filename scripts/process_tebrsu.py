#!/usr/bin/env python3
"""Process ORCHIDEE tebrsu (temperate broadleaf summergreen) data.

This script processes ISIMIP2b ORCHIDEE data for the tebrsu PFT
(temperate deciduous broadleaf trees - birch proxy for continental USA).

Data characteristics:
- Model: ORCHIDEE (ISIMIP2b)
- PFT: tebrsu (TEmperate BRoadleaf SUmmergreen)
- Period: 2006-2099
- Scenarios: RCP (rcp26, rcp60, rcp85)
- GCMs: GFDL-ESM2M, HadGEM2-ES, IPSL-CM5A-LR, MIROC5
"""

import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import xarray as xr
from scipy import stats


def load_mask(mask_file: Path) -> np.ndarray:
    """Load spatial mask from NetCDF file.

    Args:
        mask_file: Path to mask NetCDF file with 'mask' variable (True=keep, False=mask)

    Returns:
        Boolean array (lat, lon) where True = keep, False = mask
    """
    ds = xr.open_dataset(mask_file)
    mask = ds['mask'].values
    ds.close()
    return mask


def extract_years(ds: xr.Dataset) -> np.ndarray:
    """Extract integer years from dataset time coordinate.

    ORCHIDEE uses 365_day calendar with 'years since 1661-1-1' units.
    """
    time_values = ds.time.values

    # Check time units attribute
    if 'units' in ds.time.attrs:
        units = ds.time.attrs['units']
        if 'years since' in units:
            # Parse reference year and add time values
            ref_year = int(units.split('since')[1].strip().split('-')[0])
            return (ref_year + time_values).astype(int)

    # If time is already decoded as cftime objects
    if hasattr(time_values[0], 'year'):
        return np.array([t.year for t in time_values])

    # Try dt accessor
    try:
        return ds.time.dt.year.values
    except AttributeError:
        pass

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
        return shared_baseline
    else:
        raise ValueError("No 2020s data found in any file!")


def compute_zero_mask(file_data: Dict) -> np.ndarray:
    """Compute mask excluding cells with 0 value in both 2020s AND 2090s.

    Args:
        file_data: Dict mapping filename to data dict

    Returns:
        Boolean array (lat, lon) where True = keep, False = mask
    """
    all_2020s = []
    all_2090s = []

    print("\n  Computing zero-value mask...")

    for fname, fdata in file_data.items():
        data = fdata['data']
        years = fdata['years']

        # Get 2020s data
        mask_2020s = (years >= 2020) & (years < 2030)
        if np.any(mask_2020s):
            median_2020s = np.nanmedian(data[mask_2020s], axis=0)
            all_2020s.append(median_2020s)

        # Get 2090s data
        mask_2090s = (years >= 2090) & (years < 2100)
        if np.any(mask_2090s):
            median_2090s = np.nanmedian(data[mask_2090s], axis=0)
            all_2090s.append(median_2090s)

    if not all_2020s or not all_2090s:
        print("  WARNING: Could not compute mask - missing 2020s or 2090s data")
        return None

    # Average across all files
    avg_2020s = np.nanmean(np.stack(all_2020s, axis=0), axis=0)
    avg_2090s = np.nanmean(np.stack(all_2090s, axis=0), axis=0)

    # Mask = True where value is NOT zero in BOTH decades
    # i.e., keep cells where either 2020s > 0 OR 2090s > 0
    zero_2020s = (avg_2020s == 0) | np.isnan(avg_2020s)
    zero_2090s = (avg_2090s == 0) | np.isnan(avg_2090s)

    # Mask out cells that are zero in BOTH
    mask = ~(zero_2020s & zero_2090s)

    n_masked = np.sum(~mask)
    n_kept = np.sum(mask)
    print(f"  Cells with value=0 in both 2020s AND 2090s: {n_masked:,}")
    print(f"  Cells kept: {n_kept:,}")

    return mask


def process_tebrsu_files(input_dir: Path, output_dir: Path, variable: str,
                         mask_file: Optional[Path] = None):
    """Process tebrsu files and create decadal statistics with shared 2020s baseline.

    Args:
        input_dir: Directory containing raw NetCDF files
        output_dir: Directory for processed output
        variable: Variable name (cveg or npp)
        mask_file: Optional path to mask NetCDF file (True=keep, False=mask)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load mask if provided
    spatial_mask = None
    if mask_file and mask_file.exists():
        print(f"Loading mask from: {mask_file}")
        spatial_mask = load_mask(mask_file)
        print(f"  Mask shape: {spatial_mask.shape}")
        print(f"  Cells to keep: {np.sum(spatial_mask)}")
        print(f"  Cells to mask: {np.sum(~spatial_mask)}")

    # Find all NetCDF files
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
        # Parse ORCHIDEE filename: orchidee_{gcm}_ewembi_{scenario}_2005soc_co2_{var}-tebrsu_global_annual_{start}_{end}.nc4
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm2m
        scenario = parts[3]  # e.g., rcp26, rcp60, rcp85

        ds = xr.open_dataset(f, decode_times=False)
        var_name = list(ds.data_vars)[0]
        years = extract_years(ds)

        # Get data and apply mask if provided
        data = ds[var_name].values  # (time, lat, lon)
        if spatial_mask is not None:
            # Apply mask: set masked cells to NaN
            data = np.where(spatial_mask, data, np.nan)

        file_data[f.name] = {
            'data': data,
            'years': years,
            'lat': ds.lat.values,
            'lon': ds.lon.values,
            'gcm': gcm,
            'scenario': scenario,
            'units': ds[var_name].attrs.get('units', 'unknown'),
        }

        scenarios_available.add(scenario)
        gcms_available.add(gcm)

        print(f"  Loaded {f.name}: {gcm} {scenario}, years {years[0]}-{years[-1]}")
        ds.close()

    print(f"\nScenarios available: {sorted(scenarios_available)}")
    print(f"GCMs available: {sorted(gcms_available)}")

    # Compute zero-value mask if no external mask provided
    if spatial_mask is None:
        print("\n" + "=" * 60)
        print("COMPUTING ZERO-VALUE MASK")
        print("=" * 60)
        spatial_mask = compute_zero_mask(file_data)

        # Apply mask to file_data
        if spatial_mask is not None:
            for fname in file_data:
                file_data[fname]['data'] = np.where(
                    spatial_mask,
                    file_data[fname]['data'],
                    np.nan
                )

    # Get dimensions and units from first file
    first_file = list(file_data.values())[0]
    lat = first_file['lat']
    lon = first_file['lon']
    units = first_file['units']

    # Step 1: Calculate shared 2020s baseline
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_2020s = collect_2020s_baseline(file_data)
    print(f"    Range: {np.nanmin(shared_2020s):.6e} - {np.nanmax(shared_2020s):.6e} {units}")

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

        # Create ensemble file
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

            # Build attributes dict
            attrs = {
                'description': f'Ensemble decadal {variable} for temperate broadleaf summergreen trees',
                'units': units,
                'variable': f'{variable}-tebrsu',
                'scenario': scenario,
                'source': 'ISIMIP2b ORCHIDEE model ensemble',
                'pft': 'tebrsu (TEmperate BRoadleaf SUmmergreen - birch/deciduous hardwood proxy)',
                'n_gcms': n_gcms,
                'gcms': ', '.join(sorted([f['gcm'] for f in scenario_files.values()])),
                'baseline_source': 'shared_across_all_scenarios',
                'baseline_decade': '2020s',
                'trend_method': f'Theil-Sen estimator ({units} per decade)',
                'percentile_direction': 'higher_is_better',
                'percentile_note': 'Normalized to 2020s baseline (50 = same as 2020s)',
            }

            # Add mask info if mask was applied
            if spatial_mask is not None:
                attrs['mask_applied'] = 'true'
                attrs['mask_criteria'] = 'value=0 in both 2020s AND 2090s across all scenarios'
                attrs['cells_masked'] = int(np.sum(~spatial_mask))
                attrs['cells_kept'] = int(np.sum(spatial_mask))

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
                attrs=attrs
            )

            ensemble_file = output_dir / f"{variable}-tebrsu_{scenario}_processed.nc"
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
    if len(sys.argv) < 3:
        print("Usage: python process_tebrsu.py <variable> <input_dir> [output_dir] [mask_file]")
        print("  variable: cveg or npp")
        print("  input_dir: Directory containing raw NetCDF files")
        print("  output_dir: Directory for processed output (optional)")
        print("  mask_file: Path to mask NetCDF file (optional)")
        sys.exit(1)

    variable = sys.argv[1]
    input_dir = Path(sys.argv[2])

    if len(sys.argv) > 3:
        output_dir = Path(sys.argv[3])
    else:
        output_dir = Path(f"data/processed/tebrsu_{variable}_annual")

    mask_file = None
    if len(sys.argv) > 4:
        mask_file = Path(sys.argv[4])

    process_tebrsu_files(input_dir, output_dir, variable, mask_file)
