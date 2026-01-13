"""Data alignment module for multi-model ISIMIP datasets.

Handles:
- Spatial grid verification and harmonization
- Time coordinate alignment and synchronization
- Calendar conversion (360_day, noleap, standard)
- Optional regridding with xesmf (if available)

This module ensures datasets from different models and sources
have compatible spatial and temporal dimensions before processing.
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import warnings

import numpy as np
import xarray as xr


# Custom exceptions
class SpatialGridMismatchError(Exception):
    """Raised when spatial grids don't match."""
    pass


class TimeCoordinateMismatchError(Exception):
    """Raised when time coordinates can't be aligned."""
    pass


@dataclass
class SpatialGridInfo:
    """Information about a spatial grid."""

    lon: np.ndarray
    lat: np.ndarray

    @property
    def n_lon(self) -> int:
        """Number of longitude points."""
        return len(self.lon)

    @property
    def n_lat(self) -> int:
        """Number of latitude points."""
        return len(self.lat)

    @property
    def lon_min(self) -> float:
        """Minimum longitude."""
        return float(self.lon[0])

    @property
    def lon_max(self) -> float:
        """Maximum longitude."""
        return float(self.lon[-1])

    @property
    def lat_min(self) -> float:
        """Minimum latitude."""
        return float(self.lat[0])

    @property
    def lat_max(self) -> float:
        """Maximum latitude."""
        return float(self.lat[-1])

    @property
    def lon_resolution(self) -> float:
        """Longitude resolution (degrees per grid cell)."""
        if len(self.lon) < 2:
            return 0.0
        return float((self.lon[-1] - self.lon[0]) / (len(self.lon) - 1))

    @property
    def lat_resolution(self) -> float:
        """Latitude resolution (degrees per grid cell)."""
        if len(self.lat) < 2:
            return 0.0
        return float((self.lat[-1] - self.lat[0]) / (len(self.lat) - 1))


def extract_spatial_grid(ds: xr.Dataset) -> SpatialGridInfo:
    """Extract spatial grid information from dataset.

    Args:
        ds: xarray Dataset with lat/lon coordinates.

    Returns:
        SpatialGridInfo with grid parameters.

    Raises:
        ValueError: If lat/lon coordinates not found.
    """
    if "lat" not in ds.coords or "lon" not in ds.coords:
        raise ValueError("Dataset must have 'lat' and 'lon' coordinates")

    lat = ds.coords["lat"].values
    lon = ds.coords["lon"].values

    return SpatialGridInfo(lon, lat)


def verify_spatial_grids(
    ds1: xr.Dataset,
    ds2: xr.Dataset,
    tolerance: float = 0.001,
    check_resolution: bool = True,
) -> None:
    """Verify that two datasets have compatible spatial grids.

    Args:
        ds1: First xarray Dataset.
        ds2: Second xarray Dataset.
        tolerance: Maximum allowed difference in coordinates (degrees).
        check_resolution: Whether to check grid resolution matches.

    Raises:
        SpatialGridMismatchError: If grids are incompatible.
    """
    grid1 = extract_spatial_grid(ds1)
    grid2 = extract_spatial_grid(ds2)

    # Check dimensions match
    if grid1.n_lon != grid2.n_lon:
        raise SpatialGridMismatchError(
            f"Longitude dimension mismatch: {grid1.n_lon} vs {grid2.n_lon}"
        )

    if grid1.n_lat != grid2.n_lat:
        raise SpatialGridMismatchError(
            f"Latitude dimension mismatch: {grid1.n_lat} vs {grid2.n_lat}"
        )

    # Check coordinate bounds
    if abs(grid1.lon_min - grid2.lon_min) > tolerance:
        raise SpatialGridMismatchError(
            f"Longitude minimum mismatch: {grid1.lon_min} vs {grid2.lon_min}"
        )

    if abs(grid1.lon_max - grid2.lon_max) > tolerance:
        raise SpatialGridMismatchError(
            f"Longitude maximum mismatch: {grid1.lon_max} vs {grid2.lon_max}"
        )

    if abs(grid1.lat_min - grid2.lat_min) > tolerance:
        raise SpatialGridMismatchError(
            f"Latitude minimum mismatch: {grid1.lat_min} vs {grid2.lat_min}"
        )

    if abs(grid1.lat_max - grid2.lat_max) > tolerance:
        raise SpatialGridMismatchError(
            f"Latitude maximum mismatch: {grid1.lat_max} vs {grid2.lat_max}"
        )

    # Check resolution if requested
    if check_resolution:
        if abs(grid1.lon_resolution - grid2.lon_resolution) > tolerance:
            raise SpatialGridMismatchError(
                f"Longitude resolution mismatch: {grid1.lon_resolution} vs {grid2.lon_resolution}"
            )

        if abs(grid1.lat_resolution - grid2.lat_resolution) > tolerance:
            raise SpatialGridMismatchError(
                f"Latitude resolution mismatch: {grid1.lat_resolution} vs {grid2.lat_resolution}"
            )


def get_calendar_type(ds: xr.Dataset) -> str:
    """Detect calendar type from dataset.

    Args:
        ds: xarray Dataset.

    Returns:
        Calendar type: "standard", "360_day", "noleap", "gregorian", etc.
    """
    # Check time coordinate attributes
    if "time" in ds.coords:
        time_attrs = ds.coords["time"].attrs
        if "calendar" in time_attrs:
            return str(time_attrs["calendar"])

    # Try to infer from time values
    if "time" in ds.coords:
        try:
            # If we can decode times, it's likely standard
            time_values = ds.coords["time"].values
            if np.issubdtype(time_values.dtype, np.datetime64):
                return "standard"
        except (TypeError, ValueError):
            pass

    # Default to standard calendar
    return "standard"


def convert_calendar(
    ds: xr.Dataset,
    source_calendar: str,
    target_calendar: str,
) -> xr.Dataset:
    """Convert dataset from one calendar to another.

    Args:
        ds: xarray Dataset with time coordinate.
        source_calendar: Source calendar type.
        target_calendar: Target calendar type.

    Returns:
        Dataset with converted time coordinate.

    Raises:
        ValueError: If conversion not supported.
    """
    # If calendars match, return as-is
    if source_calendar == target_calendar:
        return ds

    # Handle common calendar conversions
    if source_calendar == "360_day" and target_calendar == "standard":
        # 360-day calendar: 12 months of 30 days each
        if "time" not in ds.coords:
            return ds

        time_values = ds.coords["time"].values

        # Get units from attributes if available
        time_attrs = ds.coords["time"].attrs
        units = time_attrs.get("units", "days since 2000-01-01")

        # Parse reference date from units
        try:
            ref_date_str = units.split("since ")[-1].strip()
            from datetime import datetime as dt
            ref_date = dt.strptime(ref_date_str.split()[0], "%Y-%m-%d")
        except (IndexError, ValueError):
            from datetime import datetime as dt
            ref_date = dt(2000, 1, 1)

        # Convert 360-day days to standard calendar
        # Each 360-day year = 360 days, so 1 day in 360-day ≈ 1.01389 days standard
        # But we'll use a simpler approach: map day numbers to dates
        try:
            from datetime import timedelta
            new_times = []
            for day_num in time_values:
                # In 360-day calendar, each month = 30 days
                # Convert numpy int64 to Python int for timedelta
                day_int = int(day_num)
                year_fraction = day_int / 360
                years_passed = int(year_fraction)
                days_in_year = int((year_fraction - years_passed) * 365)

                new_date = ref_date + timedelta(days=day_int)
                new_times.append(new_date)

            ds = ds.assign_coords(time=new_times)
        except (TypeError, ValueError) as e:
            warnings.warn(f"Could not convert {source_calendar} to {target_calendar}")
            raise TimeCoordinateMismatchError(
                f"Calendar conversion from {source_calendar} to {target_calendar} failed: {e}"
            )

    elif source_calendar == "noleap" and target_calendar == "standard":
        # Noleap calendar: 365 days per year, no leap years
        if "time" not in ds.coords:
            return ds

        try:
            from datetime import datetime as dt, timedelta
            time_values = ds.coords["time"].values

            # Get reference date
            time_attrs = ds.coords["time"].attrs
            units = time_attrs.get("units", "days since 2000-01-01")
            ref_date_str = units.split("since ")[-1].strip()
            ref_date = dt.strptime(ref_date_str.split()[0], "%Y-%m-%d")

            # Convert noleap days to standard calendar
            new_times = []
            for day_num in time_values:
                new_date = ref_date + timedelta(days=int(day_num))
                new_times.append(new_date)

            ds = ds.assign_coords(time=new_times)
        except (TypeError, ValueError) as e:
            warnings.warn(f"Could not convert {source_calendar} to {target_calendar}")
            raise TimeCoordinateMismatchError(
                f"Calendar conversion from {source_calendar} to {target_calendar} failed: {e}"
            )

    else:
        warnings.warn(f"Calendar conversion {source_calendar}->{target_calendar} not supported")

    return ds


def align_time_coordinates(
    ds1: xr.Dataset,
    ds2: xr.Dataset,
    method: str = "intersection",
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Align time coordinates of two datasets.

    Args:
        ds1: First xarray Dataset.
        ds2: Second xarray Dataset.
        method: "intersection" (default) or "union".
                - intersection: keep only common time periods
                - union: combine all times (pads with NaN)

    Returns:
        Tuple of aligned datasets.

    Raises:
        TimeCoordinateMismatchError: If time coordinates can't be aligned.
    """
    if "time" not in ds1.coords or "time" not in ds2.coords:
        raise TimeCoordinateMismatchError("Both datasets must have 'time' coordinates")

    # Get time values
    time1 = ds1.coords["time"].values
    time2 = ds2.coords["time"].values

    # Check if times are datetime64 or numeric
    are_datetime = np.issubdtype(time1.dtype, np.datetime64)

    if method == "intersection":
        # Find common time period
        if are_datetime:
            common_times = np.intersect1d(time1, time2)
            if len(common_times) == 0:
                raise TimeCoordinateMismatchError(
                    "No overlapping time periods between datasets"
                )
            aligned1 = ds1.sel(time=common_times)
            aligned2 = ds2.sel(time=common_times)
        else:
            # For numeric times, find common range
            min_time = max(time1.min(), time2.min())
            max_time = min(time1.max(), time2.max())

            if min_time > max_time:
                raise TimeCoordinateMismatchError(
                    "No overlapping time periods between datasets"
                )

            aligned1 = ds1.sel(time=slice(min_time, max_time))
            aligned2 = ds2.sel(time=slice(min_time, max_time))

    elif method == "union":
        # Use xarray's built-in reindex to align
        all_times = np.union1d(time1, time2)
        aligned1 = ds1.reindex(time=all_times)
        aligned2 = ds2.reindex(time=all_times)

    else:
        raise ValueError(f"Unknown alignment method: {method}")

    return aligned1, aligned2


def align_datasets(
    datasets: List[xr.Dataset],
    verify_spatial: bool = True,
    time_method: str = "intersection",
) -> List[xr.Dataset]:
    """Align multiple datasets for ensemble processing.

    Args:
        datasets: List of xarray Datasets.
        verify_spatial: Whether to verify spatial grids match.
        time_method: "intersection" or "union" for time alignment.

    Returns:
        List of aligned datasets.

    Raises:
        ValueError: If datasets list is empty.
        SpatialGridMismatchError: If spatial grids incompatible.
        TimeCoordinateMismatchError: If time coordinates can't be aligned.
    """
    if not datasets:
        raise ValueError("Cannot align empty dataset list")

    if len(datasets) == 1:
        return datasets

    # Verify and align spatial grids if requested
    if verify_spatial:
        reference_ds = datasets[0]
        for ds in datasets[1:]:
            verify_spatial_grids(reference_ds, ds)

    # Align time coordinates pairwise
    aligned = [datasets[0]]
    for ds in datasets[1:]:
        # Get previous aligned dataset and current dataset
        prev_aligned = aligned[-1]

        # Check calendars and convert if needed
        cal1 = get_calendar_type(prev_aligned)
        cal2 = get_calendar_type(ds)

        if cal1 != cal2:
            # Convert to standard calendar
            prev_aligned = convert_calendar(prev_aligned, cal1, "standard")
            ds = convert_calendar(ds, cal2, "standard")

        # Align times
        aligned_prev, aligned_curr = align_time_coordinates(
            prev_aligned, ds, method=time_method
        )

        # Update aligned list with new aligned version
        if len(aligned) > 1:
            # Reindex all previous datasets to match new time coordinates
            new_times = aligned_curr.coords["time"].values
            for i in range(len(aligned)):
                aligned[i] = aligned[i].reindex(time=new_times)

        aligned[-1] = aligned_prev
        aligned.append(aligned_curr)

    return aligned


def regrid_to_reference(
    ds: xr.Dataset,
    reference_ds: xr.Dataset,
    method: str = "nearest_s2d",
) -> xr.Dataset:
    """Regrid dataset to reference grid using xesmf (if available).

    Args:
        ds: Dataset to regrid.
        reference_ds: Reference dataset with target grid.
        method: Regridding method ("nearest_s2d", "bilinear", etc).

    Returns:
        Regridded dataset.

    Raises:
        ImportError: If xesmf not available.
        ValueError: If regridding fails.
    """
    try:
        import xesmf as xe
    except ImportError:
        raise ImportError(
            "xesmf required for regridding. "
            "Install with: pip install xesmf"
        )

    # Get lat/lon from reference
    ref_lat = reference_ds.coords["lat"].values
    ref_lon = reference_ds.coords["lon"].values

    # Create regridder
    regridder = xe.Regridder(
        ds,
        reference_ds,
        method=method,
        periodic=True,  # Periodic for lon dimension
    )

    # Apply regridding to all data variables
    regridded = regridder(ds)

    return regridded


def harmonize_datasets(
    datasets: List[xr.Dataset],
    reference_ds: Optional[xr.Dataset] = None,
    verify_spatial: bool = True,
    regrid: bool = False,
    regrid_method: str = "nearest_s2d",
) -> List[xr.Dataset]:
    """Fully harmonize a list of datasets for ensemble processing.

    Includes spatial verification, time alignment, calendar conversion,
    and optional regridding.

    Args:
        datasets: List of xarray Datasets to harmonize.
        reference_ds: Optional reference dataset for regridding.
                      If None, first dataset is used as reference.
        verify_spatial: Whether to verify spatial grids.
        regrid: Whether to regrid all datasets to reference grid.
        regrid_method: Method for regridding if enabled.

    Returns:
        List of harmonized datasets.

    Raises:
        ValueError: If datasets list is empty.
        SpatialGridMismatchError: If spatial grids incompatible.
    """
    if not datasets:
        raise ValueError("Cannot harmonize empty dataset list")

    # Set reference dataset
    if reference_ds is None:
        reference_ds = datasets[0]

    # First: align all datasets
    aligned = align_datasets(datasets, verify_spatial=verify_spatial)

    # Second: regrid if requested
    if regrid:
        aligned = [
            regrid_to_reference(ds, reference_ds, method=regrid_method)
            if ds is not reference_ds
            else ds
            for ds in aligned
        ]

    return aligned
