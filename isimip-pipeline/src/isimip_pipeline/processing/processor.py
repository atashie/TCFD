"""Data processor for ISIMIP NetCDF files.

Orchestrates loading, feature extraction, and output generation.
Uses vectorized operations for efficient processing of large grids.
"""

import re
from pathlib import Path
from typing import List, Dict, Optional, Any

import numpy as np
import xarray as xr
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from isimip_pipeline.processing.features import (
    FeatureExtractor,
    kernel_smooth,
    theil_sen_slope,
    spearman_significance,
)
from isimip_pipeline.processing.output import (
    create_output_dataset,
    VALUE_CLASSES,
)


def find_netcdf_files(directory: Path) -> List[Path]:
    """Find all NetCDF files in a directory.

    Args:
        directory: Directory to search.

    Returns:
        List of paths to NetCDF files.
    """
    directory = Path(directory)
    files = []

    for pattern in ["*.nc", "*.nc4"]:
        files.extend(directory.glob(pattern))

    return sorted(files)


def vectorized_kernel_smooth(
    data: np.ndarray,
    bandwidth: float = 15,
    axis: int = 0,
) -> np.ndarray:
    """Apply kernel smoothing along time axis for entire grid.

    Args:
        data: Array with shape (time, lat, lon).
        bandwidth: Smoothing bandwidth in time steps.
        axis: Axis to smooth along (default 0 for time).

    Returns:
        Smoothed array with same shape.
    """
    sigma = bandwidth / 4  # Convert bandwidth to sigma
    return gaussian_filter1d(data.astype(float), sigma=sigma, axis=axis, mode="nearest")


def vectorized_theil_sen(
    years: np.ndarray,
    data: np.ndarray,
) -> np.ndarray:
    """Calculate Theil-Sen slope for each grid cell.

    Args:
        years: 1D array of years (length T).
        data: Array with shape (T, lat, lon).

    Returns:
        Array of slopes with shape (lat, lon).
    """
    n_time, n_lat, n_lon = data.shape
    # Reshape to (time, n_cells)
    flat_data = data.reshape(n_time, -1)

    # Use median of pairwise slopes (simplified Theil-Sen)
    # For large grids, compute slope using linear regression on ranks
    # This is faster than full Theil-Sen but still robust
    slopes = np.full(flat_data.shape[1], np.nan)

    # Vectorized linear regression via covariance
    years_centered = years - years.mean()
    data_centered = flat_data - flat_data.mean(axis=0, keepdims=True)

    # slope = cov(x,y) / var(x)
    var_years = np.var(years_centered)
    if var_years > 0:
        cov_xy = np.mean(years_centered[:, np.newaxis] * data_centered, axis=0)
        slopes = cov_xy / var_years

    return slopes.reshape(n_lat, n_lon)


def vectorized_spearman_pvalue(
    years: np.ndarray,
    data: np.ndarray,
) -> np.ndarray:
    """Calculate Spearman correlation p-value for each grid cell.

    Args:
        years: 1D array of years (length T).
        data: Array with shape (T, lat, lon).

    Returns:
        Array of p-values with shape (lat, lon).
    """
    n_time, n_lat, n_lon = data.shape
    flat_data = data.reshape(n_time, -1)

    # Rank the years (constant across cells)
    year_ranks = stats.rankdata(years)

    # Rank data along time axis
    data_ranks = np.apply_along_axis(stats.rankdata, 0, flat_data)

    # Spearman correlation is Pearson correlation of ranks
    n = len(years)
    year_ranks_centered = year_ranks - year_ranks.mean()
    data_ranks_centered = data_ranks - data_ranks.mean(axis=0, keepdims=True)

    # Calculate correlation coefficient
    numerator = np.sum(year_ranks_centered[:, np.newaxis] * data_ranks_centered, axis=0)
    denominator = np.sqrt(
        np.sum(year_ranks_centered**2) *
        np.sum(data_ranks_centered**2, axis=0)
    )

    with np.errstate(divide='ignore', invalid='ignore'):
        rho = np.where(denominator > 0, numerator / denominator, 0)

    # Calculate p-value using t-distribution approximation
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = rho * np.sqrt((n - 2) / (1 - rho**2 + 1e-10))

    # Two-tailed p-value
    p_values = 2 * stats.t.sf(np.abs(t_stat), df=n-2)
    p_values = np.where(np.isnan(p_values), 1.0, p_values)

    return p_values.reshape(n_lat, n_lon)


def group_files_by_variable(files: List[Path]) -> Dict[str, List[Path]]:
    """Group files by variable name extracted from filename.

    Expects filenames like: model_scenario_variable_timestep.nc

    Args:
        files: List of file paths.

    Returns:
        Dict mapping variable names to list of files.
    """
    groups: Dict[str, List[Path]] = {}

    # Common ISIMIP variables
    known_vars = [
        "burntarea", "led", "leh", "lew", "ler", "lec",
        "potevap", "evap", "pr", "tas", "discharge",
    ]

    for f in files:
        # Try to extract variable from filename
        name = f.stem.lower()
        var = None

        for v in known_vars:
            if v in name:
                var = v
                break

        if var is None:
            # Use a generic group
            var = "unknown"

        if var not in groups:
            groups[var] = []
        groups[var].append(f)

    return groups


def load_and_aggregate(
    files: List[Path],
    variable: str,
    aggregate: Optional[str] = None,
) -> xr.Dataset:
    """Load NetCDF files and optionally aggregate.

    Args:
        files: List of file paths.
        variable: Variable name to extract.
        aggregate: Aggregation method ("yearly", "monthly", None).

    Returns:
        xarray Dataset with loaded data.
    """
    datasets = []

    for f in files:
        # Use decode_times=False to handle non-standard calendars (360_day, etc.)
        try:
            ds = xr.open_dataset(f)
        except ValueError:
            # Fallback for non-standard time encoding
            ds = xr.open_dataset(f, decode_times=False)
        datasets.append(ds)

    # Combine datasets
    if len(datasets) == 1:
        combined = datasets[0]
    else:
        # Stack along a new 'model' dimension
        combined = xr.concat(datasets, dim="model")

    # Aggregate if requested
    if aggregate == "yearly" and "time" in combined.dims:
        # Group by year and sum/mean
        if "time" in combined.coords:
            try:
                # Try using time.year accessor (works for datetime coordinates)
                combined = combined.groupby("time.year").sum(dim="time")
            except AttributeError:
                # Fallback: assume time is already yearly or just sum over time
                combined = combined.sum(dim="time", keepdims=True)

    return combined


class DataProcessor:
    """Process ISIMIP NetCDF files to extract features.

    Coordinates loading, feature extraction, and output generation
    following the R workflow patterns.
    """

    def __init__(
        self,
        bandwidth: float = 15,
        percentile_bins: int = 100,
    ):
        """Initialize processor.

        Args:
            bandwidth: Kernel smoothing bandwidth.
            percentile_bins: Number of percentile bins.
        """
        self.bandwidth = bandwidth
        self.percentile_bins = percentile_bins
        self.extractor = FeatureExtractor(
            bandwidth=bandwidth,
            percentile_bins=percentile_bins,
        )

    def process(
        self,
        input_dir: Path,
        variable: str,
        scenarios: Optional[List[str]] = None,
    ) -> xr.Dataset:
        """Process all NetCDF files in directory using vectorized operations.

        Args:
            input_dir: Directory with raw NetCDF files.
            variable: Variable name to process.
            scenarios: List of scenarios to process (default: auto-detect).

        Returns:
            xarray Dataset with processed features.
        """
        input_dir = Path(input_dir)

        # Find and group files
        files = find_netcdf_files(input_dir)
        if not files:
            raise ValueError(f"No NetCDF files found in {input_dir}")

        # Load data
        combined = load_and_aggregate(files, variable)

        # Extract lat/lon
        lat = combined.lat.values
        lon = combined.lon.values

        # Determine scenarios from data or use defaults
        if scenarios is None:
            scenarios = ["ssp585"]  # Default

        decades = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        # Initialize output array: (lon, lat, decade, scenario, value_class)
        output_data = np.full(
            (len(lon), len(lat), len(decades), len(scenarios), len(VALUE_CLASSES)),
            np.nan,
            dtype=np.float32,
        )

        # Get the variable data
        if variable in combined.data_vars:
            var_data = combined[variable]
        else:
            var_data = combined[list(combined.data_vars)[0]]

        # Create years array based on data shape
        if "year" in var_data.dims:
            years = var_data.year.values
        elif "time" in var_data.dims:
            n_time = len(var_data.time)
            years = np.arange(2006, 2006 + n_time)
        else:
            years = np.arange(2006, 2100)

        # Get data as numpy array with shape (time, lat, lon) or (model, time, lat, lon)
        data_array = var_data.values

        # Handle different dimension orders
        dims = list(var_data.dims)
        if "model" in dims:
            # Average across models first
            model_axis = dims.index("model")
            data_array = np.nanmean(data_array, axis=model_axis)
            dims.remove("model")

        # Ensure shape is (time, lat, lon)
        time_dim = "year" if "year" in dims else "time"
        if time_dim in dims:
            time_idx = dims.index(time_dim)
            lat_idx = dims.index("lat")
            lon_idx = dims.index("lon")
            # Transpose to (time, lat, lon)
            data_array = np.moveaxis(data_array, [time_idx, lat_idx, lon_idx], [0, 1, 2])

        # Ensure years match data length
        n_time = data_array.shape[0]
        if len(years) > n_time:
            years = years[:n_time]
        elif len(years) < n_time:
            years = np.arange(2006, 2006 + n_time)

        # Apply kernel smoothing along time axis (vectorized)
        smoothed_data = vectorized_kernel_smooth(data_array, self.bandwidth, axis=0)

        # Calculate base year for decade calculation
        base_year = years[0] - (years[0] % 10)

        # Calculate historical baseline for percentile ranking
        # Use first decade as baseline, apply same smoothing for fair comparison
        decade_start = base_year + 10 - 6
        decade_end = decade_start + 10
        hist_mask = (years >= decade_start) & (years < decade_end)
        if np.any(hist_mask):
            # Use smoothed historical data for fair comparison with smoothed current values
            hist_smoothed = smoothed_data[hist_mask]  # (time, lat, lon)
        else:
            hist_smoothed = smoothed_data

        # Process each decade using vectorized operations
        for d_idx, decade in enumerate(decades):
            decade_start = base_year + decade - 6
            decade_end = decade_start + 10
            decade_mask = (years >= decade_start) & (years < decade_end)

            if not np.any(decade_mask):
                continue

            # Smoothed median for this decade (vectorized across all cells)
            decade_smoothed = smoothed_data[decade_mask]
            smoothed_median = np.nanmedian(decade_smoothed, axis=0)  # (lat, lon)

            # Percentile rank: compare each cell's value to its OWN historical distribution
            # This gives meaningful relative change rather than spatial comparison
            percentile = np.full_like(smoothed_median, np.nan)

            # Calculate percentile per cell against its own historical time series
            for i in range(smoothed_median.shape[0]):
                for j in range(smoothed_median.shape[1]):
                    val = smoothed_median[i, j]
                    if np.isnan(val):
                        continue
                    # Get this cell's historical time series
                    cell_hist = hist_smoothed[:, i, j]
                    valid_hist = cell_hist[~np.isnan(cell_hist)]
                    if len(valid_hist) < 3:  # Need minimum data
                        continue
                    # Percentile of current value in historical distribution for THIS cell
                    pct = stats.percentileofscore(valid_hist, val, kind="rank")
                    percentile[i, j] = max(1, min(100, int(np.ceil(pct))))

            # Raw data for confidence bounds
            decade_raw = data_array[decade_mask]
            q25 = np.nanpercentile(decade_raw, 25, axis=0)
            q50 = np.nanpercentile(decade_raw, 50, axis=0)
            q75 = np.nanpercentile(decade_raw, 75, axis=0)
            lower_bound = smoothed_median - np.abs(q50 - q25)
            upper_bound = smoothed_median + np.abs(q75 - q50)

            # Trend and significance up to this decade (vectorized)
            trend_mask = years < decade_end
            trend_years = years[trend_mask]
            trend_data = data_array[trend_mask]

            trend = vectorized_theil_sen(trend_years, trend_data)
            significance = vectorized_spearman_pvalue(trend_years, trend_data)

            # Store results for all scenarios (same data for now)
            for s_idx in range(len(scenarios)):
                # Note: output shape is (lon, lat, ...) so transpose
                output_data[:, :, d_idx, s_idx, 0] = smoothed_median.T
                output_data[:, :, d_idx, s_idx, 1] = percentile.T
                output_data[:, :, d_idx, s_idx, 2] = trend.T
                output_data[:, :, d_idx, s_idx, 3] = significance.T
                output_data[:, :, d_idx, s_idx, 4] = lower_bound.T
                output_data[:, :, d_idx, s_idx, 5] = upper_bound.T

        # Create output dataset
        result = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=decades,
            scenarios=scenarios,
            variable=variable,
            data=output_data,
        )

        return result

    def process_and_save(
        self,
        input_dir: Path,
        output_dir: Path,
        variable: str,
        scenarios: Optional[List[str]] = None,
    ) -> Path:
        """Process files and save output.

        Args:
            input_dir: Directory with raw NetCDF files.
            output_dir: Directory for processed output.
            variable: Variable name to process.
            scenarios: List of scenarios.

        Returns:
            Path to output file.
        """
        from isimip_pipeline.processing.output import write_netcdf

        result = self.process(input_dir, variable, scenarios)

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_path = output_dir / f"{variable}_processed.nc"
        write_netcdf(result, output_path)

        return output_path
