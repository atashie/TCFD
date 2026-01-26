#!/usr/bin/env python3
"""Process CLASSIC evgndltr (evergreen needleleaf) NPP and cveg data.

This script processes ISIMIP3b CLASSIC model data for the evgndltr (evergreen
needleleaf trees) PFT - the most generic conifer class combining all needleleaf
evergreen trees (both boreal and temperate).

Data characteristics:
- Model: CLASSIC (ISIMIP3b)
- Variables: npp (net primary productivity), cveg (vegetation carbon)
- PFT: evgndltr (all evergreen needleleaf combined)
- Period: 2015-2100
- Scenarios: SSP (ssp126, ssp370, ssp585)
- NPP: Monthly resolution (aggregated to annual using mean)
- cveg: Annual resolution
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr
from scipy import stats


def calculate_trend(data: np.ndarray, decades: list) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate linear trend across decades using Theil-Sen estimator.

    DEPRECATED: Use calculate_cumulative_trends for per-decade trends.

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
                except Exception:
                    pass

    return slope, p_value


def calculate_cumulative_trends(
    annual_data: np.ndarray,
    years: np.ndarray,
    decades: list,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate cumulative trend from earliest year to end of each decade.

    For each decade, calculates the Theil-Sen slope from the earliest available
    year to the end of that decade. This gives cumulative trends where:
    - 2020s trend = slope from 2010 to 2029
    - 2030s trend = slope from 2010 to 2039
    - etc.

    Args:
        annual_data: Array of shape (n_years, lat, lon) with annual values
        years: Array of years corresponding to annual_data
        decades: List of decade start years [2020, 2030, ...] (NOT 2010s)

    Returns:
        Tuple of (trend, trend_pvalue) arrays of shape (n_decades, lat, lon)
    """
    n_decades = len(decades)
    n_lat, n_lon = annual_data.shape[1], annual_data.shape[2]

    trend = np.full((n_decades, n_lat, n_lon), np.nan)
    trend_pvalue = np.full((n_decades, n_lat, n_lon), np.nan)

    earliest_year = int(years.min())

    # Pre-identify valid cells (have at least some non-NaN values)
    # to avoid processing ocean/masked cells
    valid_cell_mask = np.any(~np.isnan(annual_data), axis=0)
    valid_cells = np.argwhere(valid_cell_mask)
    n_valid = len(valid_cells)
    print(f"    Processing {n_valid:,} valid grid cells...")

    for d_idx, decade_start in enumerate(decades):
        decade_end = decade_start + 9  # e.g., 2020 -> 2029

        # Get all data from earliest year to end of this decade
        mask = (years >= earliest_year) & (years <= decade_end)
        subset_years = years[mask]
        subset_data = annual_data[mask]

        n_years_subset = len(subset_years)
        if n_years_subset < 3:
            continue

        print(f"      Decade {decade_start}s: {earliest_year}-{decade_end} ({n_years_subset} years)...", end=" ", flush=True)
        cells_processed = 0

        # Only process valid cells
        for idx, (i, j) in enumerate(valid_cells):
            y = subset_data[:, i, j]
            valid = ~np.isnan(y)
            if np.sum(valid) >= 3:
                try:
                    result = stats.theilslopes(y[valid], subset_years[valid])
                    trend[d_idx, i, j] = result.slope
                    _, p = stats.kendalltau(subset_years[valid], y[valid])
                    trend_pvalue[d_idx, i, j] = p
                    cells_processed += 1
                except Exception:
                    pass

        print(f"done ({cells_processed:,} cells with trends)")

    return trend, trend_pvalue


def aggregate_monthly_to_annual(ds: xr.Dataset, var_name: str) -> xr.DataArray:
    """Aggregate monthly data to annual using mean.

    Args:
        ds: Dataset with monthly time coordinate
        var_name: Name of variable to aggregate

    Returns:
        DataArray with annual time coordinate
    """
    data = ds[var_name]
    # Group by year and take mean
    annual = data.groupby('time.year').mean(dim='time')
    return annual


def load_spatial_mask(mask_file: Path) -> np.ndarray:
    """Load spatial mask from NetCDF file.

    Args:
        mask_file: Path to mask NetCDF file

    Returns:
        Boolean mask array (lat, lon) where True = valid, False = masked
    """
    if mask_file.exists():
        ds = xr.open_dataset(mask_file)
        mask = ds['mask'].values
        ds.close()
        print(f"  Loaded spatial mask: {mask.sum():,} valid cells")
        return mask
    return None


def collect_2020s_baseline(file_data: Dict, is_monthly: bool = False,
                           spatial_mask: np.ndarray = None) -> np.ndarray:
    """Collect 2020s data from ALL files and compute shared baseline.

    Args:
        file_data: Dict mapping filename to data dict
        is_monthly: Whether data is monthly (needs aggregation)
        spatial_mask: Optional boolean mask (True=valid, False=masked)

    Returns:
        Shared 2020s baseline array (lat, lon)
    """
    all_2020s = []

    print("\n  Collecting 2020s baseline from all files...")

    for fname, fdata in file_data.items():
        data = fdata['data']  # (time, lat, lon) - already aggregated if monthly
        years = fdata['years']

        # Apply spatial mask if provided
        if spatial_mask is not None:
            data = data.copy()
            data[:, ~spatial_mask] = np.nan

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
        valid = shared_baseline[~np.isnan(shared_baseline)]
        print(f"    Range: {np.nanmin(valid):.6e} - {np.nanmax(valid):.6e}")
        return shared_baseline
    else:
        raise ValueError("No 2020s data found in any file!")


def process_variable(
    input_dir: Path,
    output_dir: Path,
    var_name: str,
    file_pattern: str,
    is_monthly: bool = False,
    units: str = "",
    long_name: str = "",
    spatial_mask: np.ndarray = None,
):
    """Process a single variable with shared 2020s baseline.

    Args:
        input_dir: Directory containing raw NetCDF files
        output_dir: Directory for processed output
        var_name: Variable name (e.g., 'npp', 'cveg')
        file_pattern: Glob pattern for files
        is_monthly: Whether data is monthly (needs aggregation)
        units: Variable units for metadata
        long_name: Long name for metadata
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all files for this variable
    nc_files = list(input_dir.glob(file_pattern))
    print(f"\n{'='*60}")
    print(f"Processing {var_name.upper()}")
    print(f"{'='*60}")
    print(f"Found {len(nc_files)} files matching {file_pattern}")

    if not nc_files:
        print("ERROR: No files found!")
        return

    # Load all files and extract metadata
    file_data = {}
    scenarios_available = set()
    gcms_available = set()

    for f in nc_files:
        # Parse CLASSIC filename: classic_{gcm}_w5e5_{scenario}_2015soc_default_{var}_global_{res}_{start}_{end}.nc
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm4
        scenario = parts[3]  # e.g., ssp126, ssp370, ssp585

        ds = xr.open_dataset(f)

        # Find the actual variable name in the dataset
        data_var = [v for v in ds.data_vars if var_name in v.lower()][0]

        if is_monthly:
            # Aggregate monthly to annual using mean
            print(f"  Aggregating {f.name} from monthly to annual (mean)...")
            annual_data = aggregate_monthly_to_annual(ds, data_var)
            data = annual_data.values
            years = annual_data.year.values
        else:
            data = ds[data_var].values
            years = ds.time.dt.year.values

        file_data[f.name] = {
            'data': data,  # (time, lat, lon)
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

    # Apply spatial mask info
    if spatial_mask is not None:
        print(f"\n  Applying spatial mask: {spatial_mask.sum():,} valid cells")

    # Step 1: Calculate shared 2020s baseline
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_2020s = collect_2020s_baseline(file_data, is_monthly, spatial_mask)

    # Step 2: Process by scenario (ensemble across available GCMs)
    decades = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]

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
        all_gcm_annual = []  # For cumulative trend calculation

        # Process each GCM
        for fname, fdata in scenario_files.items():
            gcm = fdata['gcm']
            data = fdata['data'].copy()
            years = fdata['years']

            # Apply spatial mask if provided
            if spatial_mask is not None:
                data[:, ~spatial_mask] = np.nan

            # Store annual data for trend calculation
            all_gcm_annual.append({'data': data, 'years': years})

            # Initialize output for this GCM
            gcm_output = np.full((n_decades, len(lat), len(lon)), np.nan)

            for d_idx, decade_start in enumerate(decades):
                decade_end = decade_start + 10
                time_mask = (years >= decade_start) & (years < decade_end)

                if decade_start == 2020:
                    # Use SHARED baseline for 2020s
                    gcm_output[d_idx] = shared_2020s
                elif np.any(time_mask):
                    decade_data = data[time_mask]
                    gcm_output[d_idx] = np.nanmedian(decade_data, axis=0)

            all_gcm_medians.append(gcm_output)

        # Create ensemble file
        stacked = np.stack(all_gcm_medians, axis=0)  # (n_gcms, n_decades, lat, lon)

        # Ensemble statistics
        ensemble_median = np.nanmean(stacked, axis=0)  # Mean across GCMs
        ensemble_lower = np.nanmin(stacked, axis=0)    # GCM spread lower
        ensemble_upper = np.nanmax(stacked, axis=0)    # GCM spread upper

        # Calculate cumulative per-decade trends from annual data
        # First, create ensemble mean of annual data
        first_gcm = all_gcm_annual[0]
        ensemble_years = first_gcm['years']
        ensemble_annual = np.stack([g['data'] for g in all_gcm_annual], axis=0)
        ensemble_annual = np.nanmean(ensemble_annual, axis=0)  # (n_years, lat, lon)

        print(f"  Calculating cumulative trends from {int(ensemble_years.min())}-{int(ensemble_years.max())}...")
        trend_slope, trend_pvalue = calculate_cumulative_trends(
            ensemble_annual, ensemble_years, decades
        )

        # Calculate percentile rank against 2020s global distribution
        baseline_2020s = ensemble_median[0]  # decade index 0 = 2020s
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
                'trend': (['decade', 'lat', 'lon'], trend_slope),
                'trend_pvalue': (['decade', 'lat', 'lon'], trend_pvalue),
                'lower_ci': (['decade', 'lat', 'lon'], ensemble_lower),
                'upper_ci': (['decade', 'lat', 'lon'], ensemble_upper),
            },
            coords={
                'decade': decades,
                'lat': lat,
                'lon': lon,
            },
            attrs={
                'description': f'Ensemble decadal {long_name} for evergreen needleleaf trees',
                'units': units,
                'variable': f'{var_name}-evgndltr',
                'scenario': scenario,
                'source': 'ISIMIP3b CLASSIC model ensemble',
                'pft': 'evgndltr (all evergreen needleleaf trees - boreal + temperate combined)',
                'n_gcms': n_gcms,
                'gcms': ', '.join(sorted([f['gcm'] for f in scenario_files.values()])),
                'baseline_source': 'shared_across_all_scenarios',
                'baseline_decade': '2020s',
                'trend_method': f'Cumulative Theil-Sen slope from earliest year to end of each decade ({units}/year)',
                'trend_note': 'Per-decade cumulative trends: 2020s=2010-2029, 2030s=2010-2039, etc.',
                'percentile_direction': 'higher_is_better',  # More productivity/carbon is good
                'percentile_note': 'Normalized to 2020s baseline (50 = same as 2020s)',
                'aggregation_method': 'mean' if is_monthly else 'none (already annual)',
                'spatial_mask': 'applied' if spatial_mask is not None else 'none',
                'spatial_mask_note': 'Masked cells with zero values in both 2020s and 2090s' if spatial_mask is not None else '',
            }
        )

        ensemble_file = output_dir / f"{var_name}-evgndltr_{scenario}_processed.nc"
        ensemble_ds.to_netcdf(ensemble_file)
        print(f"  Ensemble saved: {ensemble_file.name} ({n_gcms} GCMs)")
        ensemble_ds.close()

    print(f"\nOutput files for {var_name}:")
    for f in sorted(output_dir.glob(f"{var_name}*.nc")):
        print(f"  {f.name}")


def main(input_dir: Path, output_dir: Path, mask_file: Path = None):
    """Process both NPP and cveg for evergreen needleleaf trees."""

    print("=" * 70)
    print("PROCESSING ISIMIP3b CLASSIC EVERGREEN NEEDLELEAF (evgndltr) DATA")
    print("=" * 70)

    # Load spatial mask if provided
    spatial_mask = None
    if mask_file and mask_file.exists():
        spatial_mask = load_spatial_mask(mask_file)

    # Process NPP (monthly -> annual mean)
    process_variable(
        input_dir=input_dir,
        output_dir=output_dir,
        var_name="npp",
        file_pattern="*npp-evgndltr*.nc",
        is_monthly=True,
        units="kg m-2 s-1",
        long_name="Net Primary Productivity",
        spatial_mask=spatial_mask,
    )

    # Process cveg (already annual)
    process_variable(
        input_dir=input_dir,
        output_dir=output_dir,
        var_name="cveg",
        file_pattern="*cveg-evgndltr*.nc",
        is_monthly=False,
        units="kg m-2",
        long_name="Carbon Mass in Vegetation",
        spatial_mask=spatial_mask,
    )

    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print("\nAll files created:")
    for f in sorted(output_dir.glob("*.nc")):
        print(f"  {f.name}")


if __name__ == "__main__":
    input_dir = Path("data/raw/evgndltr_npp-cveg")
    output_dir = Path("data/processed/evgndltr_npp-cveg")
    mask_file = output_dir / "evgndltr_valid_mask.nc"  # Default mask location

    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])
    if len(sys.argv) > 3:
        mask_file = Path(sys.argv[3])

    main(input_dir, output_dir, mask_file)
