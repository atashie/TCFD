#!/usr/bin/env python3
"""Process loblolly pine proxy data (evergreen needleleaf cwood) with shared 2020s baseline.

This script processes ISIMIP3b cwood-evgndltr data for evergreen needleleaf trees,
which serves as a proxy for loblolly pine timber yield.

Implements:
- Shared 2020s baseline (identical across all scenarios)
- Adaptive windowing (17-year windows for 6 datasets)
- Per-pixel trend calculation
- Percentile ranking against 2020s baseline
"""

import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import xarray as xr


# Configuration
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
WINDOW_SIZE = 17  # Adaptive windowing: ceil(100/6) = 17 years


def extract_years(ds: xr.Dataset) -> np.ndarray:
    """Extract integer years from dataset time coordinate."""
    time_values = ds.time.values
    if hasattr(time_values[0], 'year'):
        return np.array([t.year for t in time_values])
    return ds.time.dt.year.values


def get_window_range(decade_start: int, window_size: int) -> Tuple[int, int]:
    """Calculate adaptive window range centered on decade midpoint.

    Args:
        decade_start: Start year of decade (e.g., 2020)
        window_size: Window size in years

    Returns:
        Tuple of (start_year, end_year) for the window
    """
    midpoint = decade_start + 5  # e.g., 2025 for 2020s
    half_window = window_size // 2
    start = midpoint - half_window
    end = start + window_size
    return start, end


def collect_2020s_baseline(gcm_data: Dict) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Collect 2020s data from ALL scenarios and average per GCM.

    Returns:
        Tuple of (shared_baselines dict, reference_distribution array)
    """
    shared_baselines = {}
    all_2020s_values = []

    for gcm, scenarios_data in gcm_data.items():
        gcm_2020s = []

        print(f"\n  Collecting 2020s baseline for {gcm}...")

        for scenario, sdata in scenarios_data.items():
            data = sdata['data']
            years = sdata['years']

            # Get adaptive window for 2020s
            start, end = get_window_range(2020, WINDOW_SIZE)
            mask = (years >= start) & (years < end)

            if np.any(mask):
                decade_data = data[mask]
                # Compute median for this scenario's 2020s
                scenario_median = np.nanmedian(decade_data, axis=0)
                gcm_2020s.append(scenario_median)
                all_2020s_values.append(decade_data)
                print(f"    {scenario}: {np.sum(mask)} years ({start}-{end-1}), "
                      f"median range: {np.nanmin(scenario_median):.2f} - {np.nanmax(scenario_median):.2f}")

        if gcm_2020s:
            # Average across scenarios for shared baseline
            stacked = np.stack(gcm_2020s, axis=0)
            shared_baseline = np.nanmean(stacked, axis=0)
            shared_baselines[gcm] = shared_baseline
            print(f"    SHARED baseline from {len(gcm_2020s)} scenarios")

    # Create reference distribution for percentile calculations
    if all_2020s_values:
        # Flatten all 2020s data for reference distribution
        reference = np.concatenate([v.flatten() for v in all_2020s_values])
        reference = reference[~np.isnan(reference)]
        print(f"\n  Reference distribution: {len(reference)} values, "
              f"range {np.min(reference):.4f} - {np.max(reference):.4f}")
    else:
        reference = np.array([])

    return shared_baselines, reference


def calculate_percentile(value: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Calculate percentile rank of values against reference distribution.

    For "higher_is_better" variables (cwood): inverts so high value → low percentile (safe).
    """
    if len(reference) == 0:
        return np.full_like(value, np.nan)

    # Vectorized percentile calculation
    result = np.full_like(value, np.nan)
    valid = ~np.isnan(value)
    if np.any(valid):
        # For each valid value, count how many reference values are below it
        pct = np.searchsorted(np.sort(reference), value[valid]) / len(reference) * 100
        # Invert: high value → low percentile (safe) for "higher_is_better"
        result[valid] = 100 - pct
    return result


def calculate_trend(data: np.ndarray, years: np.ndarray) -> np.ndarray:
    """Calculate linear trend (slope) per pixel.

    Args:
        data: Array of shape (time, lat, lon)
        years: Array of years

    Returns:
        Trend array of shape (lat, lon) in units per year
    """
    if len(years) < 2:
        return np.full(data.shape[1:], np.nan)

    # Centered time values for regression
    t = years.astype(float) - years.mean()
    var_t = np.var(t)

    if var_t == 0:
        return np.full(data.shape[1:], np.nan)

    # Calculate trend using covariance method
    y_centered = data - np.nanmean(data, axis=0, keepdims=True)
    cov = np.nanmean(t[:, np.newaxis, np.newaxis] * y_centered, axis=0)
    trend = cov / var_t

    return trend


def process_evgndltr_files(input_dir: Path, output_dir: Path):
    """Process evgndltr cwood files with shared 2020s baseline and trends."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all cwood files
    nc_files = list(input_dir.glob("*.nc"))
    print(f"Found {len(nc_files)} NetCDF files")
    print(f"Using adaptive window size: {WINDOW_SIZE} years")

    # Group files by GCM
    gcm_data = {}
    for f in nc_files:
        parts = f.stem.split("_")
        gcm = parts[1]  # e.g., gfdl-esm4
        scenario = parts[3]  # e.g., ssp126, ssp370, ssp585

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

    # Step 1: Calculate shared 2020s baseline
    print("\n" + "=" * 60)
    print("CALCULATING SHARED 2020s BASELINE")
    print("=" * 60)
    shared_baselines, reference_dist = collect_2020s_baseline(gcm_data)

    # Step 2: Process each scenario
    projection_scenarios = ['ssp126', 'ssp370', 'ssp585']

    print("\n" + "=" * 60)
    print("PROCESSING SCENARIOS")
    print("=" * 60)

    for gcm, scenarios_data in gcm_data.items():
        print(f"\nProcessing GCM: {gcm}")

        first_scenario = list(scenarios_data.values())[0]
        lat = first_scenario['lat']
        lon = first_scenario['lon']

        shared_2020s = shared_baselines.get(gcm)
        if shared_2020s is None:
            print(f"  WARNING: No shared baseline for {gcm}, skipping...")
            continue

        for scenario in projection_scenarios:
            if scenario not in scenarios_data:
                print(f"  Skipping {scenario} (not available)")
                continue

            print(f"\n  Processing {scenario}...")

            # Initialize output arrays
            n_decades = len(DECADES)
            median_data = np.full((n_decades, len(lat), len(lon)), np.nan, dtype=np.float32)
            percentile_data = np.full((n_decades, len(lat), len(lon)), np.nan, dtype=np.float32)
            trend_data = np.full((n_decades, len(lat), len(lon)), np.nan, dtype=np.float32)

            data = scenarios_data[scenario]['data']
            years = scenarios_data[scenario]['years']

            for d_idx, decade_start in enumerate(DECADES):
                start, end = get_window_range(decade_start, WINDOW_SIZE)
                mask = (years >= start) & (years < end)

                if decade_start == 2020:
                    # Use SHARED baseline
                    median_data[d_idx] = shared_2020s
                    percentile_data[d_idx] = calculate_percentile(shared_2020s, reference_dist)
                    # Trend from all scenarios' 2020s data
                    if np.any(mask):
                        window_data = data[mask]
                        window_years = years[mask]
                        trend_data[d_idx] = calculate_trend(window_data, window_years)
                    print(f"    {decade_start}s: SHARED BASELINE ({start}-{end-1})")

                elif np.any(mask):
                    # Scenario-specific data for 2030s-2090s
                    window_data = data[mask]
                    window_years = years[mask]

                    decade_median = np.nanmedian(window_data, axis=0)
                    median_data[d_idx] = decade_median
                    percentile_data[d_idx] = calculate_percentile(decade_median, reference_dist)
                    trend_data[d_idx] = calculate_trend(window_data, window_years)

                    print(f"    {decade_start}s: {np.sum(mask)} years ({start}-{end-1}), "
                          f"median range: {np.nanmin(decade_median):.2f} - {np.nanmax(decade_median):.2f}")

            # Create output dataset
            output_ds = xr.Dataset(
                data_vars={
                    'cwood_median': (['decade', 'lat', 'lon'], median_data, {
                        'long_name': 'Decadal median wood carbon',
                        'units': 'kg m-2',
                    }),
                    'cwood_percentile': (['decade', 'lat', 'lon'], percentile_data, {
                        'long_name': 'Percentile rank vs 2020s baseline',
                        'units': 'percentile (0-100)',
                    }),
                    'cwood_trend': (['decade', 'lat', 'lon'], trend_data, {
                        'long_name': 'Linear trend within decade window',
                        'units': 'kg m-2 year-1',
                    }),
                },
                coords={
                    'decade': DECADES,
                    'lat': lat,
                    'lon': lon,
                },
                attrs={
                    'description': 'Decadal statistics for evergreen needleleaf wood carbon (loblolly pine proxy)',
                    'variable': 'cwood-evgndltr',
                    'pft': 'evgndltr (evergreen needleleaf trees)',
                    'gcm': gcm,
                    'scenario': scenario,
                    'source': 'ISIMIP3b CLASSIC model',
                    'baseline_source': 'shared_across_all_scenarios',
                    'baseline_decade': '2020s',
                    'window_size': f'{WINDOW_SIZE} years (adaptive)',
                    'processing_note': 'Loblolly pine timber proxy using evergreen needleleaf PFT',
                    'percentile_direction': 'higher_is_better',  # High wood carbon = low risk
                }
            )

            output_file = output_dir / f"cwood_evgndltr_{scenario}_{gcm}_decadal.nc"
            output_ds.to_netcdf(output_file)
            print(f"    Saved: {output_file}")
            output_ds.close()

    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE!")
    print("=" * 60)
    print(f"\nOutput files in: {output_dir}")
    print("\nVerify shared baseline with:")
    print(f"  python scripts/test_shared_baseline.py {output_dir}")


if __name__ == "__main__":
    input_dir = Path("data/raw/loblolly-pine-proxy_cwood-evgndltr_annual")
    output_dir = Path("data/processed/loblolly-pine-proxy_cwood-evgndltr_annual")

    if len(sys.argv) > 1:
        input_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    process_evgndltr_files(input_dir, output_dir)
