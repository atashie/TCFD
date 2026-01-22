"""ISIMIP land-sea mask utility module.

Download, cache, and apply official ISIMIP land-sea masks to exclude
ocean or land areas from climate data processing.

Supported masks:
- ISIMIP2b generic (0.5 degree, 720x360)
- ISIMIP3b generic (0.5 degree, 720x360)
"""

from pathlib import Path
from typing import Optional, Union
import urllib.request

import numpy as np
import xarray as xr

# Official ISIMIP land-sea mask URLs
MASK_URLS = {
    "2b": "https://files.isimip.org/ISIMIP2b/InputData/landseamask/ISIMIP2b_landseamask_generic.nc4",
    "3b": "https://files.isimip.org/ISIMIP3b/InputData/geo_conditions/landseamask/landseamask.nc",
}

# Default cache directory
DEFAULT_CACHE_DIR = Path("data/masks")


def get_isimip_landmask(
    version: str = "2b",
    cache_dir: Optional[Path] = None,
    force_download: bool = False,
) -> Path:
    """Download and cache ISIMIP land-sea mask.

    Args:
        version: ISIMIP version ("2b" or "3b")
        cache_dir: Directory to cache masks (default: data/masks/)
        force_download: Re-download even if cached

    Returns:
        Path to the cached mask file
    """
    if version not in MASK_URLS:
        raise ValueError(f"Unknown ISIMIP version: {version}. Use '2b' or '3b'.")

    cache_dir = cache_dir or DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    mask_filename = f"ISIMIP{version}_landseamask.nc4"
    mask_path = cache_dir / mask_filename

    if mask_path.exists() and not force_download:
        print(f"Using cached mask: {mask_path}")
        return mask_path

    url = MASK_URLS[version]
    print(f"Downloading ISIMIP{version} land-sea mask from {url}...")

    urllib.request.urlretrieve(url, mask_path)
    print(f"Saved to: {mask_path}")

    return mask_path


def load_landmask(path: Union[str, Path]) -> xr.DataArray:
    """Load land-sea mask as xarray DataArray.

    Args:
        path: Path to mask NetCDF file

    Returns:
        DataArray with 1 for land, NaN/0 for ocean (2D: lat, lon)
    """
    ds = xr.open_dataset(path)

    # Find the mask variable (different names in different versions)
    mask_var = None
    for var in ["LSM", "mask", "land", "landsea", "lsm"]:
        if var in ds.data_vars:
            mask_var = var
            break

    if mask_var is None:
        # Use the first data variable
        mask_var = list(ds.data_vars)[0]

    mask = ds[mask_var]

    # Squeeze out time dimension if present (ISIMIP2b mask has time=1)
    if "time" in mask.dims:
        mask = mask.squeeze("time", drop=True)

    ds.close()

    return mask


def apply_land_mask(
    data: np.ndarray,
    mask: Union[xr.DataArray, np.ndarray],
    invert: bool = False,
) -> np.ndarray:
    """Apply land-sea mask to data array.

    Sets ocean cells to NaN (or land cells if inverted).

    Args:
        data: Data array of shape (..., lat, lon)
        mask: Mask array of shape (lat, lon), 1=land, NaN/0=ocean
        invert: If True, mask land instead of ocean (for marine data)

    Returns:
        Masked data array with ocean (or land) set to NaN
    """
    if isinstance(mask, xr.DataArray):
        mask_values = mask.values
    else:
        mask_values = mask

    # Ensure mask is 2D (lat, lon)
    if mask_values.ndim > 2:
        # Squeeze extra dimensions
        mask_values = np.squeeze(mask_values)

    # ISIMIP masks use 1 for land, NaN for ocean
    # Create boolean mask: land=True where value is 1 (not NaN)
    land_mask = ~np.isnan(mask_values) & (mask_values > 0)

    if invert:
        # For marine data: keep ocean, mask land
        valid_mask = ~land_mask
    else:
        # For land data: keep land, mask ocean
        valid_mask = land_mask

    # Apply mask - set invalid cells to NaN
    # Use broadcasting for multi-dimensional data (time, lat, lon)
    masked_data = data.astype(np.float64).copy()
    masked_data[..., ~valid_mask] = np.nan

    return masked_data


def verify_grid_compatibility(
    data_lat: np.ndarray,
    data_lon: np.ndarray,
    mask: xr.DataArray,
    tolerance: float = 0.01,
) -> bool:
    """Verify that data and mask grids are compatible.

    Args:
        data_lat: Latitude coordinates from data
        data_lon: Longitude coordinates from data
        mask: Mask DataArray with lat/lon coordinates
        tolerance: Maximum allowed difference in grid points

    Returns:
        True if grids are compatible

    Raises:
        ValueError if grids are incompatible
    """
    mask_lat = mask.lat.values if hasattr(mask, "lat") else mask.coords["lat"].values
    mask_lon = mask.lon.values if hasattr(mask, "lon") else mask.coords["lon"].values

    # Check dimensions match
    if len(data_lat) != len(mask_lat):
        raise ValueError(
            f"Latitude dimension mismatch: data has {len(data_lat)}, "
            f"mask has {len(mask_lat)}"
        )

    if len(data_lon) != len(mask_lon):
        raise ValueError(
            f"Longitude dimension mismatch: data has {len(data_lon)}, "
            f"mask has {len(mask_lon)}"
        )

    # Check coordinate values match
    lat_diff = np.abs(data_lat - mask_lat).max()
    lon_diff = np.abs(data_lon - mask_lon).max()

    if lat_diff > tolerance:
        raise ValueError(
            f"Latitude coordinates differ by up to {lat_diff:.4f} degrees"
        )

    if lon_diff > tolerance:
        raise ValueError(
            f"Longitude coordinates differ by up to {lon_diff:.4f} degrees"
        )

    return True


def get_mask_stats(mask: Union[xr.DataArray, np.ndarray]) -> dict:
    """Get statistics about the mask.

    Args:
        mask: Land-sea mask array (1=land, NaN/0=ocean)

    Returns:
        Dictionary with mask statistics
    """
    if isinstance(mask, xr.DataArray):
        mask_values = mask.values
    else:
        mask_values = mask

    # Squeeze extra dimensions if needed
    if mask_values.ndim > 2:
        mask_values = np.squeeze(mask_values)

    # ISIMIP masks: 1=land, NaN=ocean
    land_cells = np.sum(~np.isnan(mask_values) & (mask_values > 0))
    ocean_cells = np.sum(np.isnan(mask_values) | (mask_values == 0))
    total_cells = mask_values.size

    return {
        "land_cells": int(land_cells),
        "ocean_cells": int(ocean_cells),
        "total_cells": int(total_cells),
        "land_fraction": float(land_cells / total_cells),
        "ocean_fraction": float(ocean_cells / total_cells),
    }
