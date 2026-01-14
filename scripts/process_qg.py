"""Process qg (groundwater runoff) data following process-metrics skill requirements.

Implements:
- Temporal scope: 2010s-2090s
- Adaptive windowing: min 100 data points per decade-bin
- Percentile baseline: 2020s reference distribution

Optimized for large global datasets using chunked processing.
"""

import sys
import numpy as np
import xarray as xr
from pathlib import Path
from typing import Dict, List, Tuple
import warnings

# Configuration following process-metrics skill
DEFAULT_DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
MIN_DATA_POINTS = 100
BASELINE_DECADE = 2020
MIN_YEAR = 2010  # Only load data from 2010 onwards
MAX_YEAR = 2099


def log(msg: str):
    """Print with immediate flush."""
    print(msg, flush=True)


def calculate_window_size(num_datasets: int) -> int:
    """Calculate adaptive window size based on number of datasets."""
    window_years = max(10, int(np.ceil(MIN_DATA_POINTS / num_datasets)))
    return window_years


def get_decade_window(decade_start: int, window_years: int) -> Tuple[int, int]:
    """Get temporal window for a decade, centered on decade midpoint."""
    if window_years == 10:
        return decade_start, decade_start + 9

    decade_center = decade_start + 4.5
    half_window = window_years / 2
    window_start = int(np.floor(decade_center - half_window + 0.5))
    window_end = window_start + window_years - 1

    return window_start, window_end


def load_and_aggregate_file(fpath: Path, min_year: int = MIN_YEAR, max_year: int = MAX_YEAR) -> xr.DataArray:
    """Load a NetCDF file, filter to relevant years, and aggregate to annual.

    Args:
        fpath: Path to NetCDF file
        min_year: Minimum year to include
        max_year: Maximum year to include

    Returns:
        DataArray with annual data (year, lat, lon)
    """
    log(f"    Loading {fpath.name}...")

    # Load dataset (no chunking - dask not required)
    ds = xr.open_dataset(fpath)

    # Get the variable
    var_name = "qg" if "qg" in ds.data_vars else list(ds.data_vars)[0]
    da = ds[var_name]

    # Filter to relevant years
    years = da.time.dt.year
    mask = (years >= min_year) & (years <= max_year)
    da_filtered = da.isel(time=mask)

    if len(da_filtered.time) == 0:
        ds.close()
        return None

    # Aggregate monthly to annual mean
    annual = da_filtered.groupby("time.year").mean(dim="time")

    ds.close()
    return annual


def collect_2020s_data_across_scenarios(
    hist_files: List[Path],
    scenario_files: Dict[str, List[Path]],
    window_years: int
) -> Tuple[xr.DataArray, np.ndarray]:
    """Collect and combine 2020s data from all scenarios.

    This ensures all climate projections start with identical 2020s baseline values,
    computed as the average across all scenarios (projection + historical/picontrol).

    Args:
        hist_files: Historical NetCDF files
        scenario_files: Dict mapping scenario name to list of future files
        window_years: Adaptive window size

    Returns:
        Tuple of:
        - Combined 2020s data array (model, year, lat, lon)
        - Flattened baseline values for percentile calculation
    """
    all_2020s_data = []  # List of (gcm, scenario, annual_data)

    # Get 2020s window bounds
    baseline_start, baseline_end = get_decade_window(BASELINE_DECADE, window_years)
    log(f"\nCollecting 2020s baseline from ALL scenarios (window: {baseline_start}-{baseline_end})...")

    # Load 2020s window from each scenario
    for scenario, future_files in scenario_files.items():
        if not future_files:
            continue

        log(f"  Loading 2020s from {scenario}...")
        for fpath in future_files:
            gcm = fpath.stem.split("_")[1]
            # Load only the 2020s window
            annual = load_and_aggregate_file(
                fpath,
                min_year=baseline_start,
                max_year=baseline_end
            )
            if annual is not None and len(annual.year) > 0:
                all_2020s_data.append((gcm, scenario, annual))
                log(f"    {gcm}: {len(annual.year)} years")

    if not all_2020s_data:
        log("  WARNING: No 2020s data collected from any scenario!")
        return None, np.array([])

    # Group by GCM and average across scenarios
    gcm_scenario_data = {}  # gcm -> {scenario -> data}
    for gcm, scenario, data in all_2020s_data:
        if gcm not in gcm_scenario_data:
            gcm_scenario_data[gcm] = {}
        gcm_scenario_data[gcm][scenario] = data

    log(f"\n  GCM contribution to shared baseline:")
    # For each GCM, average its 2020s data across scenarios
    combined_gcms = []
    for gcm, scenario_data in gcm_scenario_data.items():
        scenario_arrays = list(scenario_data.values())

        if len(scenario_arrays) == 1:
            gcm_avg = scenario_arrays[0]
        else:
            # Align years across scenarios and average
            common_years = set(scenario_arrays[0].year.values)
            for arr in scenario_arrays[1:]:
                common_years &= set(arr.year.values)
            common_years = sorted(common_years)

            if common_years:
                aligned = [arr.sel(year=common_years) for arr in scenario_arrays]
                stacked = xr.concat(aligned, dim="scenario")
                gcm_avg = stacked.mean(dim="scenario")
            else:
                gcm_avg = scenario_arrays[0]

        combined_gcms.append(gcm_avg.expand_dims({"model": [gcm]}))
        log(f"    {gcm}: {len(gcm_avg.year)} years from {len(scenario_data)} scenarios")

    # Stack all GCMs
    combined_baseline = xr.concat(combined_gcms, dim="model")
    log(f"\n  Combined 2020s baseline shape: {combined_baseline.shape}")
    log(f"  Models: {list(combined_baseline.model.values)}")

    # Create flattened baseline for percentile calculation
    baseline_sample = combined_baseline.values[:, :, ::10, ::10]  # Subsample spatially
    baseline_flat = baseline_sample.flatten()
    baseline_flat = baseline_flat[~np.isnan(baseline_flat)]
    log(f"  Baseline sample size: {len(baseline_flat):,}")

    return combined_baseline, baseline_flat


def process_scenario(
    hist_files: List[Path],
    future_files: List[Path],
    scenario: str,
    output_dir: Path,
    window_years: int
):
    """Process a single scenario.

    Args:
        hist_files: Historical NetCDF files
        future_files: Future scenario NetCDF files
        scenario: Scenario name
        output_dir: Output directory
        window_years: Adaptive window size
    """
    log(f"\n{'='*60}")
    log(f"Processing scenario: {scenario}")
    log(f"{'='*60}")

    all_data = []  # List of (gcm_name, annual_data)

    # Load historical files (only 2010-2014)
    log("\nLoading historical data (2010-2014)...")
    for fpath in hist_files:
        gcm = fpath.stem.split("_")[1]
        annual = load_and_aggregate_file(fpath, min_year=2010, max_year=2014)
        if annual is not None:
            all_data.append((gcm, "historical", annual))
            log(f"      {gcm}: {len(annual.year)} years ({annual.year.values[0]}-{annual.year.values[-1]})")

    # Load future files (2015-2099)
    log(f"\nLoading {scenario} data (2015-2099)...")
    for fpath in future_files:
        gcm = fpath.stem.split("_")[1]
        annual = load_and_aggregate_file(fpath, min_year=2015, max_year=2099)
        if annual is not None:
            all_data.append((gcm, scenario, annual))
            log(f"      {gcm}: {len(annual.year)} years ({annual.year.values[0]}-{annual.year.values[-1]})")

    if not all_data:
        log("  No data loaded!")
        return

    # Group by GCM and concatenate historical + future
    log("\nCombining historical and future data per GCM...")
    gcm_data = {}
    for gcm, period, data in all_data:
        if gcm not in gcm_data:
            gcm_data[gcm] = []
        gcm_data[gcm].append(data)

    # Concatenate time series for each GCM
    combined_gcms = []
    for gcm, data_list in gcm_data.items():
        if len(data_list) > 1:
            combined = xr.concat(data_list, dim="year")
            combined = combined.sortby("year")
        else:
            combined = data_list[0]
        combined_gcms.append(combined.expand_dims({"model": [gcm]}))
        log(f"    {gcm}: years {combined.year.values[0]}-{combined.year.values[-1]}")

    # Stack all GCMs
    all_models = xr.concat(combined_gcms, dim="model")
    log(f"\nCombined data shape: {all_models.shape}")
    log(f"  Models: {list(all_models.model.values)}")
    log(f"  Years: {all_models.year.values[0]}-{all_models.year.values[-1]}")

    # Get coordinates
    years = all_models.year.values
    lats = all_models.lat.values
    lons = all_models.lon.values

    # Calculate 2020s baseline
    log(f"\nCalculating 2020s baseline (window: {window_years} years)...")
    baseline_start, baseline_end = get_decade_window(BASELINE_DECADE, window_years)
    log(f"  Baseline window: {baseline_start}-{baseline_end}")

    baseline_mask = (years >= baseline_start) & (years <= baseline_end)
    baseline_years = years[baseline_mask]

    if len(baseline_years) == 0:
        log("  WARNING: No data in baseline window!")
        baseline_flat = np.array([])
    else:
        baseline_data = all_models.sel(year=baseline_years)
        # Sample for percentile calculation (use subset to save memory)
        baseline_sample = baseline_data.values[:, :, ::10, ::10]  # Subsample spatially
        baseline_flat = baseline_sample.flatten()
        baseline_flat = baseline_flat[~np.isnan(baseline_flat)]
        log(f"  Baseline sample size: {len(baseline_flat):,}")

    # Initialize output
    n_decades = len(DEFAULT_DECADES)
    result = {
        "median": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "percentile": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "trend": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "lower_ci": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "upper_ci": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
    }

    # Process each decade
    for d_idx, decade_start in enumerate(DEFAULT_DECADES):
        window_start, window_end = get_decade_window(decade_start, window_years)

        # Get years in window
        window_mask = (years >= window_start) & (years <= window_end)
        window_years_arr = years[window_mask]

        if len(window_years_arr) == 0:
            log(f"  Decade {decade_start}s: No data")
            continue

        log(f"  Decade {decade_start}s: {window_years_arr[0]}-{window_years_arr[-1]} ({len(window_years_arr)} years)")

        # Extract data for this window
        window_data = all_models.sel(year=window_years_arr)

        # Calculate statistics
        # Mean across models, then median across time
        model_mean = window_data.mean(dim="model")
        decade_median = model_mean.median(dim="year").values
        result["median"][d_idx] = decade_median

        # Percentile rank against 2020s baseline
        if len(baseline_flat) > 100:
            flat_median = decade_median.flatten()
            percentiles = np.zeros_like(flat_median)
            for i, val in enumerate(flat_median):
                if np.isnan(val):
                    percentiles[i] = np.nan
                else:
                    pct = np.sum(baseline_flat <= val) / len(baseline_flat) * 100
                    percentiles[i] = max(1, min(100, pct))
            result["percentile"][d_idx] = percentiles.reshape(len(lats), len(lons))

        # Confidence intervals
        all_vals = window_data.values.reshape(-1, len(lats), len(lons))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result["lower_ci"][d_idx] = np.nanpercentile(all_vals, 25, axis=0)
            result["upper_ci"][d_idx] = np.nanpercentile(all_vals, 75, axis=0)

        # Trend (simplified)
        if len(window_years_arr) >= 3:
            t = window_years_arr - window_years_arr.mean()
            var_t = np.var(t)
            if var_t > 0:
                y_centered = model_mean.values - model_mean.values.mean(axis=0, keepdims=True)
                cov = np.mean(t[:, np.newaxis, np.newaxis] * y_centered, axis=0)
                result["trend"][d_idx] = cov / var_t

    # Create output dataset
    log("\nCreating output dataset...")
    ds_out = xr.Dataset(
        {
            "median": (["decade", "lat", "lon"], result["median"]),
            "percentile": (["decade", "lat", "lon"], result["percentile"]),
            "trend": (["decade", "lat", "lon"], result["trend"]),
            "lower_ci": (["decade", "lat", "lon"], result["lower_ci"]),
            "upper_ci": (["decade", "lat", "lon"], result["upper_ci"]),
        },
        coords={
            "decade": DEFAULT_DECADES,
            "lat": lats,
            "lon": lons,
        },
        attrs={
            "variable": "qg",
            "scenario": scenario,
            "long_name": "Groundwater runoff",
            "units": "kg m-2 s-1",
            "window_years": window_years,
            "baseline_decade": BASELINE_DECADE,
            "description": "Processed following process-metrics skill requirements",
        }
    )

    # Save
    output_path = output_dir / f"qg_{scenario}_processed.nc"
    ds_out.to_netcdf(output_path)
    log(f"  Saved: {output_path}")


def process_scenario_with_shared_baseline(
    hist_files: List[Path],
    future_files: List[Path],
    scenario: str,
    output_dir: Path,
    window_years: int,
    shared_2020s_data: xr.DataArray,
    shared_baseline_flat: np.ndarray
):
    """Process a single scenario using shared 2020s baseline.

    For the 2020s decade, uses shared_2020s_data (averaged across all scenarios).
    For 2030s-2090s, uses this scenario's data only.

    Args:
        hist_files: Historical NetCDF files
        future_files: Future scenario NetCDF files
        scenario: Scenario name
        output_dir: Output directory
        window_years: Adaptive window size
        shared_2020s_data: Shared 2020s data from all scenarios
        shared_baseline_flat: Flattened baseline for percentile calculation
    """
    log(f"\n{'='*60}")
    log(f"Processing scenario: {scenario} (with SHARED 2020s baseline)")
    log(f"{'='*60}")

    all_data = []  # List of (gcm_name, period, annual_data)

    # Load historical files (only 2010-2014)
    log("\nLoading historical data (2010-2014)...")
    for fpath in hist_files:
        gcm = fpath.stem.split("_")[1]
        annual = load_and_aggregate_file(fpath, min_year=2010, max_year=2014)
        if annual is not None:
            all_data.append((gcm, "historical", annual))
            log(f"      {gcm}: {len(annual.year)} years")

    # Load future files (2015-2099)
    log(f"\nLoading {scenario} data (2015-2099)...")
    for fpath in future_files:
        gcm = fpath.stem.split("_")[1]
        annual = load_and_aggregate_file(fpath, min_year=2015, max_year=2099)
        if annual is not None:
            all_data.append((gcm, scenario, annual))
            log(f"      {gcm}: {len(annual.year)} years")

    if not all_data:
        log(f"  No data loaded for {scenario}")
        return

    # Group by GCM and concatenate historical + future
    log("\nCombining historical and future data per GCM...")
    gcm_data = {}
    for gcm, period, data in all_data:
        if gcm not in gcm_data:
            gcm_data[gcm] = []
        gcm_data[gcm].append(data)

    combined_gcms = []
    for gcm, data_list in gcm_data.items():
        if len(data_list) > 1:
            combined = xr.concat(data_list, dim="year")
            combined = combined.sortby("year")
        else:
            combined = data_list[0]
        combined_gcms.append(combined.expand_dims({"model": [gcm]}))
        log(f"    {gcm}: years {combined.year.values[0]}-{combined.year.values[-1]}")

    all_models = xr.concat(combined_gcms, dim="model")
    years = all_models.year.values
    lats = all_models.lat.values
    lons = all_models.lon.values

    log(f"\nCombined data shape: {all_models.shape}")
    log(f"  Models: {list(all_models.model.values)}")
    log(f"  Years: {years[0]}-{years[-1]}")

    # Initialize result arrays
    n_decades = len(DEFAULT_DECADES)
    result = {
        "median": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "percentile": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "trend": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "lower_ci": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
        "upper_ci": np.full((n_decades, len(lats), len(lons)), np.nan, dtype=np.float32),
    }

    # Process each decade
    for d_idx, decade_start in enumerate(DEFAULT_DECADES):
        window_start, window_end = get_decade_window(decade_start, window_years)

        # KEY CHANGE: For 2020s, use SHARED data
        if decade_start == BASELINE_DECADE and shared_2020s_data is not None:
            log(f"  Decade {decade_start}s: Using SHARED baseline across all scenarios")

            # Use shared 2020s data
            shared_years = shared_2020s_data.year.values
            window_mask = (shared_years >= window_start) & (shared_years <= window_end)
            window_years_arr = shared_years[window_mask]

            if len(window_years_arr) == 0:
                log(f"    No data in shared baseline window!")
                continue

            log(f"    Window: {window_years_arr[0]}-{window_years_arr[-1]} ({len(window_years_arr)} years)")

            window_data = shared_2020s_data.sel(year=window_years_arr)
            model_mean = window_data.mean(dim="model")
            decade_median = model_mean.median(dim="year").values
            result["median"][d_idx] = decade_median

            # Percentile against shared baseline
            if len(shared_baseline_flat) > 100:
                flat_median = decade_median.flatten()
                percentiles = np.zeros_like(flat_median)
                for i, val in enumerate(flat_median):
                    if np.isnan(val):
                        percentiles[i] = np.nan
                    else:
                        pct = np.sum(shared_baseline_flat <= val) / len(shared_baseline_flat) * 100
                        percentiles[i] = max(1, min(100, pct))
                result["percentile"][d_idx] = percentiles.reshape(len(lats), len(lons))

            # Confidence intervals from shared data
            all_vals = window_data.values.reshape(-1, len(lats), len(lons))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result["lower_ci"][d_idx] = np.nanpercentile(all_vals, 25, axis=0)
                result["upper_ci"][d_idx] = np.nanpercentile(all_vals, 75, axis=0)

            # Trend (use shared data years)
            if len(window_years_arr) >= 3:
                t = window_years_arr - window_years_arr.mean()
                var_t = np.var(t)
                if var_t > 0:
                    y_centered = model_mean.values - model_mean.values.mean(axis=0, keepdims=True)
                    cov = np.mean(t[:, np.newaxis, np.newaxis] * y_centered, axis=0)
                    result["trend"][d_idx] = cov / var_t

        else:
            # For 2010s and 2030s-2090s: use THIS scenario's data only
            log(f"  Decade {decade_start}s: Using {scenario}-specific data")

            window_mask = (years >= window_start) & (years <= window_end)
            window_years_arr = years[window_mask]

            if len(window_years_arr) == 0:
                log(f"    No data in window {window_start}-{window_end}")
                continue

            log(f"    Window: {window_years_arr[0]}-{window_years_arr[-1]} ({len(window_years_arr)} years)")

            window_data = all_models.sel(year=window_years_arr)
            model_mean = window_data.mean(dim="model")
            decade_median = model_mean.median(dim="year").values
            result["median"][d_idx] = decade_median

            # Percentile against shared 2020s baseline
            if len(shared_baseline_flat) > 100:
                flat_median = decade_median.flatten()
                percentiles = np.zeros_like(flat_median)
                for i, val in enumerate(flat_median):
                    if np.isnan(val):
                        percentiles[i] = np.nan
                    else:
                        pct = np.sum(shared_baseline_flat <= val) / len(shared_baseline_flat) * 100
                        percentiles[i] = max(1, min(100, pct))
                result["percentile"][d_idx] = percentiles.reshape(len(lats), len(lons))

            # Confidence intervals
            all_vals = window_data.values.reshape(-1, len(lats), len(lons))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result["lower_ci"][d_idx] = np.nanpercentile(all_vals, 25, axis=0)
                result["upper_ci"][d_idx] = np.nanpercentile(all_vals, 75, axis=0)

            # Trend
            if len(window_years_arr) >= 3:
                t = window_years_arr - window_years_arr.mean()
                var_t = np.var(t)
                if var_t > 0:
                    y_centered = model_mean.values - model_mean.values.mean(axis=0, keepdims=True)
                    cov = np.mean(t[:, np.newaxis, np.newaxis] * y_centered, axis=0)
                    result["trend"][d_idx] = cov / var_t

    # Create output dataset
    log("\nCreating output dataset...")
    ds_out = xr.Dataset(
        {
            "median": (["decade", "lat", "lon"], result["median"]),
            "percentile": (["decade", "lat", "lon"], result["percentile"]),
            "trend": (["decade", "lat", "lon"], result["trend"]),
            "lower_ci": (["decade", "lat", "lon"], result["lower_ci"]),
            "upper_ci": (["decade", "lat", "lon"], result["upper_ci"]),
        },
        coords={
            "decade": DEFAULT_DECADES,
            "lat": lats,
            "lon": lons,
        },
        attrs={
            "variable": "qg",
            "scenario": scenario,
            "long_name": "Groundwater runoff",
            "units": "kg m-2 s-1",
            "window_years": window_years,
            "baseline_decade": BASELINE_DECADE,
            "baseline_source": "shared_across_all_scenarios",
            "description": "Processed with shared 2020s baseline across all scenarios",
        }
    )

    # Save
    output_path = output_dir / f"qg_{scenario}_processed.nc"
    ds_out.to_netcdf(output_path)
    log(f"  Saved: {output_path}")


def main():
    """Main processing function."""
    # Get project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "data" / "raw"
    output_dir = project_root / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)

    log("=" * 60)
    log("Processing qg (groundwater runoff) data")
    log("Following process-metrics skill requirements")
    log("=" * 60)

    # Group files by scenario
    files = sorted(data_dir.glob("*.nc"))
    log(f"\nFound {len(files)} NetCDF files")

    hist_files = [f for f in files if "historical" in f.name.lower()]
    ssp126_files = [f for f in files if "ssp126" in f.name.lower()]
    ssp370_files = [f for f in files if "ssp370" in f.name.lower()]
    ssp585_files = [f for f in files if "ssp585" in f.name.lower()]

    log(f"  Historical: {len(hist_files)}")
    log(f"  SSP126: {len(ssp126_files)}")
    log(f"  SSP370: {len(ssp370_files)}")
    log(f"  SSP585: {len(ssp585_files)}")

    # Check if we have data to process
    if not files:
        log("\nERROR: No NetCDF files found in data/raw/")
        log("Please download raw data first using: isimip-pipeline run")
        return

    if not hist_files:
        log("\nERROR: No historical files found. Cannot establish baseline.")
        return

    # Calculate adaptive window
    num_gcms = len(hist_files)
    window_years = calculate_window_size(num_gcms)
    log(f"\nNumber of GCMs: {num_gcms}")
    log(f"Adaptive window size: {window_years} years")

    # Collect shared 2020s baseline from ALL scenarios
    scenario_files = {
        "ssp126": ssp126_files,
        "ssp370": ssp370_files,
        "ssp585": ssp585_files,
    }

    # Add picontrol files if available (for enhanced baseline robustness)
    picontrol_files = [f for f in files if "picontrol" in f.name.lower()]
    if picontrol_files:
        scenario_files["picontrol"] = picontrol_files
        log(f"  Picontrol: {len(picontrol_files)} (will be used for 2020s baseline only)")

    # Collect shared baseline across all scenarios
    shared_2020s_data, shared_baseline_flat = collect_2020s_data_across_scenarios(
        hist_files, scenario_files, window_years
    )

    if shared_2020s_data is None:
        log("\nERROR: Could not collect shared 2020s baseline. Falling back to per-scenario processing.")
        # Fallback to original behavior
        for scenario, future_files in [
            ("ssp126", ssp126_files),
            ("ssp370", ssp370_files),
            ("ssp585", ssp585_files),
        ]:
            if future_files:
                process_scenario(hist_files, future_files, scenario, output_dir, window_years)
    else:
        # Process each scenario with shared baseline
        log("\n" + "=" * 60)
        log("Processing scenarios with SHARED 2020s baseline")
        log("=" * 60)
        for scenario, future_files in [
            ("ssp126", ssp126_files),
            ("ssp370", ssp370_files),
            ("ssp585", ssp585_files),
        ]:
            if future_files:
                process_scenario_with_shared_baseline(
                    hist_files, future_files, scenario, output_dir, window_years,
                    shared_2020s_data, shared_baseline_flat
                )

    log("\n" + "=" * 60)
    log("Processing complete!")
    log("=" * 60)


if __name__ == "__main__":
    main()
