#!/usr/bin/env python3
"""Combine per-GCM oak cwood files into ensemble statistics per scenario.

Creates generate_maps.py compatible output: {variable}_{scenario}_processed.nc
with median, percentile, lower_ci, upper_ci, and trend variables.

Supports latitude masking to exclude anomalous bands before processing.
"""

import numpy as np
import xarray as xr
from pathlib import Path
from scipy import stats
from typing import Optional, Tuple, List


def apply_latitude_mask(data: np.ndarray, lat: np.ndarray,
                        mask_bands: Optional[List[Tuple[float, float]]] = None) -> np.ndarray:
    """Apply latitude mask to data array.

    Args:
        data: Array with shape (..., lat, lon) where lat is second-to-last dimension
        lat: Latitude coordinate array
        mask_bands: List of (min_lat, max_lat) tuples to mask out

    Returns:
        Data array with masked latitudes set to NaN
    """
    if mask_bands is None:
        return data

    masked_data = data.copy()
    for min_lat, max_lat in mask_bands:
        # Find indices to mask
        mask_idx = (lat >= min_lat) & (lat <= max_lat)
        # Apply mask - set to NaN
        if masked_data.ndim == 3:  # (decade, lat, lon)
            masked_data[:, mask_idx, :] = np.nan
        elif masked_data.ndim == 4:  # (gcm, decade, lat, lon)
            masked_data[:, :, mask_idx, :] = np.nan
        elif masked_data.ndim == 2:  # (lat, lon)
            masked_data[mask_idx, :] = np.nan

    return masked_data


def calculate_trend(median: np.ndarray, decades: np.ndarray) -> np.ndarray:
    """Calculate decadal trend (Theil-Sen slope) for each grid cell.

    Uses Theil-Sen estimator which is robust to outliers (up to 29.3% can be outliers).
    This is the recommended method for climate data.

    Args:
        median: Array of shape (n_decades, lat, lon)
        decades: Array of decade values (e.g., [2020, 2030, ...])

    Returns:
        Trend array of shape (lat, lon) with units of [value/decade]
    """
    n_decades, n_lat, n_lon = median.shape
    trend = np.full((n_lat, n_lon), np.nan)

    for i in range(n_lat):
        for j in range(n_lon):
            values = median[:, i, j]
            valid_mask = ~np.isnan(values)

            if np.sum(valid_mask) >= 3:  # Need at least 3 points for regression
                valid_decades = decades[valid_mask]
                valid_values = values[valid_mask]
                # Theil-Sen: median of slopes between all pairs of points
                result = stats.theilslopes(valid_values, valid_decades)
                slope = result.slope  # slope is per year
                # Convert to change per decade (multiply by 10)
                trend[i, j] = slope * 10

    return trend


def combine_gcm_files(input_dir: Path, output_dir: Path,
                      variable: str = "cwood-dcdcldbdltr",
                      mask_bands: Optional[List[Tuple[float, float]]] = None):
    """Combine GCM-specific files into ensemble statistics per scenario.

    Args:
        input_dir: Directory containing per-GCM decadal files
        output_dir: Directory for output files
        variable: Variable name for output files
        mask_bands: List of (min_lat, max_lat) tuples to mask out
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all decadal files (exclude any existing _processed files)
    nc_files = [f for f in input_dir.glob("*.nc") if "_processed" not in f.name]
    print(f"Found {len(nc_files)} NetCDF files")

    if mask_bands:
        print(f"\nLatitude masking enabled:")
        for min_lat, max_lat in mask_bands:
            print(f"  Masking {min_lat}° to {max_lat}°")

    # Group files by scenario
    scenario_files = {}
    for f in nc_files:
        parts = f.stem.split("_")
        # Pattern: cwood_dcdcldbdltr_ssp126_gfdl-esm4_decadal
        scenario = parts[2]  # e.g., ssp126
        if scenario not in scenario_files:
            scenario_files[scenario] = []
        scenario_files[scenario].append(f)

    print(f"Scenarios found: {list(scenario_files.keys())}")

    for scenario, files in scenario_files.items():
        print(f"\nProcessing {scenario} with {len(files)} GCM files...")

        # Load all GCM data
        gcm_data = []
        for f in files:
            ds = xr.open_dataset(f)
            var_name = list(ds.data_vars)[0]  # e.g., cwood_median
            gcm_data.append(ds[var_name].values)
            lat = ds.lat.values
            lon = ds.lon.values
            decades = ds.decade.values
            units = ds.attrs.get('units', 'kg m-2')
            ds.close()

        # Stack GCMs: (n_gcms, n_decades, lat, lon)
        stacked = np.stack(gcm_data, axis=0)
        n_gcms = stacked.shape[0]
        print(f"  Stacked shape: {stacked.shape} ({n_gcms} GCMs)")

        # Apply latitude mask if specified
        if mask_bands:
            stacked = apply_latitude_mask(stacked, lat, mask_bands)
            masked_count = np.sum(np.isnan(stacked[0, 0, :, :])) - np.sum(np.isnan(gcm_data[0][0]))
            print(f"  Applied mask: ~{int(masked_count)} additional NaN cells per timestep")

        # Calculate ensemble statistics
        median = np.nanmedian(stacked, axis=0)  # (decade, lat, lon)
        lower_ci = np.nanpercentile(stacked, 25, axis=0)
        upper_ci = np.nanpercentile(stacked, 75, axis=0)

        # Calculate percentile (vs 2020s baseline)
        # 2020s is first decade
        baseline = median[0]  # (lat, lon)
        percentile = np.full_like(median, np.nan)

        # Flatten for percentile calculation
        baseline_flat = baseline.flatten()
        valid_baseline = baseline_flat[~np.isnan(baseline_flat)]
        print(f"  Baseline (2020s): {len(valid_baseline)} valid cells")

        for d_idx in range(len(decades)):
            decade_data = median[d_idx].flatten()
            pct_values = np.full_like(decade_data, np.nan)

            for i, val in enumerate(decade_data):
                if not np.isnan(val) and len(valid_baseline) > 0:
                    pct = np.sum(valid_baseline <= val) / len(valid_baseline) * 100
                    pct_values[i] = max(1, min(100, pct))

            percentile[d_idx] = pct_values.reshape(median[d_idx].shape)

        # Calculate trend (slope across decades)
        print(f"  Calculating trends...")
        trend = calculate_trend(median, decades)
        valid_trends = np.sum(~np.isnan(trend))
        print(f"  Trend calculated for {valid_trends} cells")

        # Build mask description for metadata
        mask_desc = "None"
        if mask_bands:
            mask_desc = "; ".join([f"{min_lat}° to {max_lat}°" for min_lat, max_lat in mask_bands])

        # Create output dataset
        output_ds = xr.Dataset(
            data_vars={
                'median': (['decade', 'lat', 'lon'], median),
                'percentile': (['decade', 'lat', 'lon'], percentile),
                'lower_ci': (['decade', 'lat', 'lon'], lower_ci),
                'upper_ci': (['decade', 'lat', 'lon'], upper_ci),
                'trend': (['lat', 'lon'], trend),
            },
            coords={
                'decade': decades,
                'lat': lat,
                'lon': lon,
            },
            attrs={
                'description': 'Decadal ensemble statistics for wood carbon (oak proxy)',
                'units': units,
                'trend_units': f'{units} per decade',
                'variable': variable,
                'scenario': scenario,
                'long_name': 'Carbon Mass in Wood - Deciduous Broadleaf Trees',
                'source': 'ISIMIP3b CLASSIC model',
                'n_gcms': n_gcms,
                'baseline_source': 'shared_across_all_scenarios',
                'baseline_decade': '2020s',
                'latitude_mask': mask_desc,
            }
        )

        # Output with expected naming convention
        output_file = output_dir / f"{variable}_{scenario}_processed.nc"
        output_ds.to_netcdf(output_file)
        print(f"  Saved: {output_file}")
        output_ds.close()

    print("\n" + "=" * 60)
    print("ENSEMBLE COMBINATION COMPLETE")
    print("=" * 60)
    print(f"Output files ready for generate_maps.py")


if __name__ == "__main__":
    import sys

    input_dir = Path("data/processed/oak-timber_cwood_annual")
    output_dir = Path("data/processed/oak-timber_cwood_annual")

    # Default latitude mask for oak data: exclude 22.25-23.75°N and 22.25-23.75°S
    # These bands have CLASSIC model artifacts at tropical/temperate boundary
    mask_bands = [
        (22.25, 23.75),    # Northern tropical boundary
        (-23.75, -22.25),  # Southern tropical boundary
    ]

    # Parse command line arguments
    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])
    if len(sys.argv) > 3 and sys.argv[3] == "--no-mask":
        mask_bands = None
        print("Masking disabled via --no-mask flag")

    combine_gcm_files(input_dir, output_dir, mask_bands=mask_bands)
