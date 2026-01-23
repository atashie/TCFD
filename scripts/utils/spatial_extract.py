"""Spatial extraction utilities for ISIMIP processed data.

Extract climate data values by:
- Point (with Gaussian distance-weighted averaging)
- Polygon (with area-weighted averaging)
- Region (Natural Earth boundaries)

All functions work with processed NetCDF files having (decade, lat, lon) dimensions.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import xarray as xr
from shapely.geometry import Point, Polygon, MultiPolygon, box


# Standard ISIMIP grid cell size
CELL_SIZE = 0.5

# Known variable directions for percentile interpretation
# "higher_is_better" means high values are GOOD (low risk) - percentile should be inverted
# "higher_is_worse" means high values are BAD (high risk) - percentile used as-is
KNOWN_HIGHER_IS_BETTER = {
    "cveg", "cwood", "cleaf", "csoil",  # Vegetation carbon
    "npp", "gpp", "ra", "rh",           # Productivity
    "tcb", "b30cm", "b10cm",            # Fish biomass
    "yield", "biom",                     # Crop output
    "soilmoist", "groundwstor",         # Water storage
}

KNOWN_HIGHER_IS_WORSE = {
    "led", "leh", "lew", "ler", "lec",  # Exposure days
    "an-tot-heat", "an-tot-cold",       # Mortality
    "burntarea", "ffire",               # Fire
    "potevap",                          # Evaporative demand
}


def get_percentile_direction(ds: xr.Dataset, variable: str) -> str:
    """Get percentile direction, checking metadata then known lists.

    Args:
        ds: xarray Dataset (may contain percentile_direction attribute)
        variable: Variable name (e.g., "cveg", "cveg-evgndltr")

    Returns:
        "higher_is_worse", "higher_is_better", or "unknown"
    """
    # 1. Check dataset metadata first
    if "percentile_direction" in ds.attrs:
        return ds.attrs["percentile_direction"]

    # 2. Check known variable lists (handle PFT suffixes like "cveg-evgndltr")
    var_base = variable.split("-")[0]
    if var_base in KNOWN_HIGHER_IS_BETTER:
        return "higher_is_better"
    if var_base in KNOWN_HIGHER_IS_WORSE:
        return "higher_is_worse"

    # 3. Unknown - requires user confirmation
    return "unknown"


def apply_percentile_inversion(
    results: Dict[str, Dict[int, float]],
    direction: str,
) -> Dict[str, Dict[int, float]]:
    """Apply percentile inversion if direction is 'higher_is_better'.

    For "higher_is_better" variables, inverts percentile so that:
    - High raw value → Low percentile (good/safe)
    - Low raw value → High percentile (bad/risky)

    Args:
        results: Extraction results dict from extract_by_* functions
        direction: "higher_is_worse", "higher_is_better", or "unknown"

    Returns:
        Results dict with percentile inverted if needed
    """
    if direction != "higher_is_better":
        return results

    if "percentile" not in results:
        return results

    # Invert: 100 - percentile
    inverted = {}
    for decade, value in results["percentile"].items():
        if np.isnan(value):
            inverted[decade] = value
        else:
            inverted[decade] = 100.0 - value

    results["percentile"] = inverted
    return results


def normalize_longitude(lon: float) -> float:
    """Normalize longitude to -180 to 180 range (EPSG:4326).

    Args:
        lon: Longitude value

    Returns:
        Normalized longitude in range [-180, 180]

    Raises:
        ValueError: If longitude cannot be normalized to valid range
    """
    if lon < -180:
        lon += 360
    elif lon > 180:
        lon -= 360

    if lon < -180 or lon > 180:
        raise ValueError(f"Longitude {lon} out of valid range [-180, 180]")

    return lon


def extract_by_point(
    ds: xr.Dataset,
    lat: float,
    lon: float,
    variables: Optional[List[str]] = None,
    search_radius: float = 0.5,
    sigma: float = 0.25,
) -> Dict[str, Dict[int, float]]:
    """Extract values at a point using Gaussian distance-weighted averaging.

    Uses Gaussian weighting from the target point to nearby cell centers.
    This provides smooth interpolation while weighting closer cells more heavily.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        lat: Target latitude (decimal degrees)
        lon: Target longitude (decimal degrees, -180 to 180)
        variables: List of variables to extract (default: all numeric variables)
        search_radius: Search radius in degrees (default: 0.5 = ~4 cells)
        sigma: Gaussian sigma for weighting (default: 0.25 = half cell size)

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Example:
        >>> ds = xr.open_dataset("processed.nc")
        >>> result = extract_by_point(ds, lat=45.5, lon=-122.7)
        >>> print(result["median"][2050])  # Value for 2050 decade
    """
    lon = normalize_longitude(lon)

    lats = ds.lat.values
    lons = ds.lon.values

    # Find cells within search radius
    lat_mask = np.abs(lats - lat) <= search_radius
    lon_mask = np.abs(lons - lon) <= search_radius

    nearby_lats = lats[lat_mask]
    nearby_lons = lons[lon_mask]

    if len(nearby_lats) == 0 or len(nearby_lons) == 0:
        raise ValueError(
            f"No grid cells found within {search_radius} degrees of ({lat}, {lon})"
        )

    # Calculate distances from target point to cell CENTERS
    lon_grid, lat_grid = np.meshgrid(nearby_lons, nearby_lats)
    distances = np.sqrt((lat_grid - lat) ** 2 + (lon_grid - lon) ** 2)

    # Gaussian weights (sigma = half cell size gives natural falloff)
    weights = np.exp(-0.5 * (distances / sigma) ** 2)
    weights = weights / weights.sum()

    # Determine which variables to extract
    if variables is None:
        # Extract all numeric variables with the right dimensions
        variables = [
            v
            for v in ds.data_vars
            if set(ds[v].dims) >= {"decade", "lat", "lon"}
        ]

    results = {}
    for var in variables:
        if var not in ds:
            continue

        data_subset = ds[var].sel(lat=nearby_lats, lon=nearby_lons)
        var_results = {}

        for decade in ds.decade.values:
            decade_data = data_subset.sel(decade=decade).values

            # Handle NaN cells by excluding and re-normalizing weights
            valid_mask = ~np.isnan(decade_data)
            if valid_mask.any():
                valid_weights = weights[valid_mask]
                valid_weights = valid_weights / valid_weights.sum()
                var_results[int(decade)] = float(
                    np.sum(decade_data[valid_mask] * valid_weights)
                )
            else:
                var_results[int(decade)] = np.nan

        results[var] = var_results

    return results


def _calculate_cell_weights(
    polygon: Union[Polygon, MultiPolygon],
    lats: np.ndarray,
    lons: np.ndarray,
    cell_size: float = CELL_SIZE,
) -> np.ndarray:
    """Calculate fraction of each grid cell inside a polygon.

    Uses shapely intersection for accurate area calculation.

    Args:
        polygon: Shapely Polygon or MultiPolygon
        lats: 1D array of cell center latitudes
        lons: 1D array of cell center longitudes
        cell_size: Grid cell size in degrees

    Returns:
        2D array of weights (shape: lat x lon), values 0-1
    """
    half_cell = cell_size / 2.0
    weights = np.zeros((len(lats), len(lons)))

    # Get polygon bounds for quick filtering
    minx, miny, maxx, maxy = polygon.bounds

    for i, lat_val in enumerate(lats):
        # Skip rows outside polygon bounds
        if lat_val + half_cell < miny or lat_val - half_cell > maxy:
            continue

        for j, lon_val in enumerate(lons):
            # Skip columns outside polygon bounds
            if lon_val + half_cell < minx or lon_val - half_cell > maxx:
                continue

            # Create cell polygon from center
            cell = box(
                lon_val - half_cell,
                lat_val - half_cell,
                lon_val + half_cell,
                lat_val + half_cell,
            )

            # Calculate intersection area as fraction of cell
            try:
                intersection = polygon.intersection(cell)
                if not intersection.is_empty:
                    weights[i, j] = intersection.area / cell.area
            except Exception:
                # Skip invalid geometry intersections
                continue

    return weights


def extract_by_polygon(
    ds: xr.Dataset,
    polygon: Union[Polygon, MultiPolygon],
    variables: Optional[List[str]] = None,
    cell_size: float = CELL_SIZE,
) -> Dict[str, Dict[int, float]]:
    """Extract values within a polygon using area-weighted averaging.

    Calculates the fraction of each grid cell that falls inside the polygon
    and uses these fractions as weights for averaging.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        polygon: Shapely Polygon or MultiPolygon defining the region
        variables: List of variables to extract (default: all numeric variables)
        cell_size: Grid cell size in degrees (default: 0.5)

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Raises:
        ValueError: If polygon doesn't intersect any grid cells

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-123, 45), (-122, 45), (-122, 46), (-123, 46)])
        >>> result = extract_by_polygon(ds, poly)
    """
    lats = ds.lat.values
    lons = ds.lon.values

    # Calculate area weights
    weights = _calculate_cell_weights(polygon, lats, lons, cell_size)

    if weights.sum() == 0:
        raise ValueError("Polygon does not intersect any grid cells")

    weights_norm = weights / weights.sum()

    # Determine which variables to extract
    if variables is None:
        variables = [
            v
            for v in ds.data_vars
            if set(ds[v].dims) >= {"decade", "lat", "lon"}
        ]

    results = {}
    for var in variables:
        if var not in ds:
            continue

        var_results = {}
        for decade in ds.decade.values:
            data = ds[var].sel(decade=decade).values

            # Handle NaN cells by excluding and re-normalizing
            valid_mask = ~np.isnan(data) & (weights > 0)
            if valid_mask.any():
                valid_weights = weights_norm[valid_mask]
                valid_weights = valid_weights / valid_weights.sum()
                var_results[int(decade)] = float(
                    np.sum(data[valid_mask] * valid_weights)
                )
            else:
                var_results[int(decade)] = np.nan

        results[var] = var_results

    return results


def extract_by_region(
    ds: xr.Dataset,
    region_name: str,
    region_type: str = "country",
    variables: Optional[List[str]] = None,
    scale: str = "50m",
) -> Dict[str, Dict[int, float]]:
    """Extract values for a named geopolitical region.

    Loads the region boundary from Natural Earth and delegates to
    extract_by_polygon for area-weighted averaging.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        region_name: Name of the region (fuzzy matched)
        region_type: Type of region ("country" or "state")
        variables: List of variables to extract (default: all numeric variables)
        scale: Natural Earth scale ("50m" or "10m")

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Raises:
        ValueError: If region not found in Natural Earth data

    Example:
        >>> result = extract_by_region(ds, "France")
        >>> result = extract_by_region(ds, "California", region_type="state")
    """
    from .natural_earth import get_region_geometry

    geometry = get_region_geometry(region_name, region_type, scale)
    if geometry is None:
        raise ValueError(f"Region not found: {region_name} (type={region_type})")

    return extract_by_polygon(ds, geometry, variables)


def extract_all_variables(
    ds: xr.Dataset,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    polygon: Optional[Union[Polygon, MultiPolygon]] = None,
    region_name: Optional[str] = None,
    region_type: str = "country",
) -> Dict[str, Dict[int, float]]:
    """Convenience function to extract all variables using appropriate method.

    Exactly one of (lat/lon), polygon, or region_name must be provided.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        lat: Latitude for point extraction
        lon: Longitude for point extraction
        polygon: Shapely polygon for area extraction
        region_name: Region name for Natural Earth extraction
        region_type: Region type if using region_name

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Raises:
        ValueError: If invalid combination of arguments provided
    """
    # Count how many location methods are specified
    has_point = lat is not None and lon is not None
    has_polygon = polygon is not None
    has_region = region_name is not None

    method_count = sum([has_point, has_polygon, has_region])

    if method_count == 0:
        raise ValueError(
            "Must provide one of: (lat, lon), polygon, or region_name"
        )
    if method_count > 1:
        raise ValueError(
            "Provide only one of: (lat, lon), polygon, or region_name"
        )

    if has_point:
        return extract_by_point(ds, lat, lon)
    elif has_polygon:
        return extract_by_polygon(ds, polygon)
    else:
        return extract_by_region(ds, region_name, region_type)


def load_processed_dataset(
    processed_dir: Path,
    folder_pattern: str,
    scenario: str,
) -> xr.Dataset:
    """Load a processed NetCDF file by folder pattern and scenario.

    Args:
        processed_dir: Base directory for processed data
        folder_pattern: Substring to match in folder name (e.g., "cwood")
        scenario: Scenario code (e.g., "ssp126", "rcp26")

    Returns:
        xarray Dataset

    Raises:
        FileNotFoundError: If no matching file found
    """
    # Find matching directory
    matches = [
        d for d in processed_dir.iterdir()
        if d.is_dir() and folder_pattern in d.name
    ]

    if not matches:
        raise FileNotFoundError(
            f"No directory matching '{folder_pattern}' in {processed_dir}"
        )

    data_dir = matches[0]

    # Find scenario file
    file_pattern = f"*{scenario}*_processed.nc"
    files = list(data_dir.glob(file_pattern))

    if not files:
        raise FileNotFoundError(
            f"No file matching '{file_pattern}' in {data_dir}"
        )

    return xr.open_dataset(files[0])
