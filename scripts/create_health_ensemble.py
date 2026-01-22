#!/usr/bin/env python3
"""Create ensemble files from per-GCM health mortality processed files.

This script combines per-GCM files into ensemble statistics without the slow
pixel-by-pixel percentile calculation.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr


def create_ensemble_files(processed_dir: Path):
    """Create ensemble files from per-GCM processed files."""

    # Find all per-GCM files
    gcm_files = list(processed_dir.glob("an-tot-heat_*_*_decadal.nc"))
    print(f"Found {len(gcm_files)} per-GCM files")

    # Group by scenario
    scenarios = {}
    for f in gcm_files:
        parts = f.stem.split("_")
        scenario = parts[1]  # e.g., rcp26, rcp60

        if scenario not in scenarios:
            scenarios[scenario] = []
        scenarios[scenario].append(f)

    print(f"Scenarios: {list(scenarios.keys())}")

    for scenario, files in scenarios.items():
        print(f"\nCreating ensemble for {scenario} ({len(files)} GCMs)...")

        gcm_medians = []
        gcm_trends = []
        gcm_names = []

        for f in files:
            ds = xr.open_dataset(f)
            gcm_medians.append(ds['mortality_median'].values)
            gcm_trends.append(ds['mortality_trend'].values)
            gcm_names.append(ds.attrs.get('gcm', f.stem.split('_')[2]))

            lat = ds.lat.values
            lon = ds.lon.values
            decades = ds.decade.values
            units = ds.attrs.get('units', 'unknown')
            long_name = ds.attrs.get('long_name', 'Heat-related mortality')
            ds.close()

        # Stack and compute ensemble statistics
        stacked_medians = np.stack(gcm_medians, axis=0)  # (n_gcm, n_decades, lat, lon)
        stacked_trends = np.stack(gcm_trends, axis=0)  # (n_gcm, lat, lon)

        ensemble_median = np.nanmean(stacked_medians, axis=0)
        ensemble_trend = np.nanmean(stacked_trends, axis=0)

        # GCM spread (lower and upper bounds)
        lower_ci = np.nanmin(stacked_medians, axis=0)
        upper_ci = np.nanmax(stacked_medians, axis=0)

        # Calculate percentile using vectorized approach
        # Get 2020s values (first decade) as baseline
        baseline_2020s = ensemble_median[0]  # (lat, lon)
        baseline_flat = baseline_2020s.flatten()
        baseline_valid = baseline_flat[~np.isnan(baseline_flat)]
        baseline_sorted = np.sort(baseline_valid)
        n_baseline = len(baseline_sorted)

        print(f"  Baseline 2020s: {n_baseline} valid points")

        # Vectorized percentile calculation
        percentile_data = np.full_like(ensemble_median, np.nan)

        if n_baseline > 0:
            for d_idx in range(len(decades)):
                decade_values = ensemble_median[d_idx].flatten()
                # Use searchsorted for vectorized percentile calculation
                ranks = np.searchsorted(baseline_sorted, decade_values)
                percentiles = 100.0 * ranks / n_baseline
                percentile_data[d_idx] = percentiles.reshape(ensemble_median[d_idx].shape)
                print(f"  Decade {decades[d_idx]}: percentile range {np.nanmin(percentiles):.1f} - {np.nanmax(percentiles):.1f}")

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
                'units': units,
                'long_name': long_name,
                'variable': 'an-tot-heat',
                'scenario': scenario,
                'source': 'ISIMIP2b TRM-Tsukuba model',
                'ensemble_members': gcm_names,
                'n_gcms': len(gcm_names),
                'baseline_source': 'shared_across_all_scenarios',
                'baseline_decade': '2020s',
                'land_masked': 'true',
                'mask_source': 'ISIMIP2b_landseamask_generic.nc4',
            }
        )

        output_file = processed_dir / f"an-tot-heat_{scenario}_processed.nc"
        ensemble_ds.to_netcdf(output_file)
        print(f"  Saved: {output_file}")
        ensemble_ds.close()

    print("\nEnsemble creation complete!")


if __name__ == "__main__":
    processed_dir = Path("data/processed/health-mortality_an-tot-heat_annual")

    if len(sys.argv) > 1:
        processed_dir = Path(sys.argv[1])

    create_ensemble_files(processed_dir)
