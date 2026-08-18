"""Spatial extraction utilities for ISIMIP processed data.

Extract climate data values by:
- Point (with Gaussian distance-weighted averaging)
- Polygon (with area-weighted averaging)
- Region (Natural Earth boundaries)

All functions work with processed NetCDF files having (decade, lat, lon) dimensions.
Observational (lat, lon) layers are lifted onto a single-period decade axis by
`as_period_dataset()` before extraction -- see its docstring.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import xarray as xr
from shapely.geometry import Point, Polygon, MultiPolygon, box


# ---------------------------------------------------------------------------------------
# Grid resolution -- INFERRED from a dataset's own coordinates, never assumed
# ---------------------------------------------------------------------------------------
#
# 0.5 deg remains the DEFAULT and the preferred grid: every layer shipped before 2026-08-14
# is 0.5 deg and none of their numbers may move. So `CELL_SIZE`, `radius = 0.5` and
# `sigma = 0.25` stay exactly what they were, and `grid_cell_size()` returns 0.5 for those
# layers, which routes them down an identical code path.
#
# What changes is only that a layer on ANOTHER regular grid (the 0.25 deg CaMa-Flood
# inundation layers, ingested 2026-08-14) no longer silently inherits 0.5 deg geometry:
# a 0.5 deg search radius on a 0.25 deg grid reaches 4x as many cells, which would blend a
# site with neighbours ~55 km away while claiming ~28 km resolution.
#
# The geometry is therefore expressed in CELLS. At 0.5 deg:
#     radius = 1.0 cell = 0.5 deg   (unchanged)
#     sigma  = 0.5 cell = 0.25 deg  (unchanged)
# ASSET-CATALOG.md "Spatial averaging -- the complete picture" documents what that produces
# and remains authoritative for the 0.5 deg case.
#
# NOT changed here, deliberately: distance is still measured in DEGREE space with no
# cos(lat) term, and polygon weights are still planar intersection fractions. Both are
# known approximations (ASSET-CATALOG.md), and correcting either would move shipped site
# values. That is a separate, declared decision -- not a side effect of resolution work.

#: Standard ISIMIP grid cell size, and the default when a dataset cannot be inspected.
CELL_SIZE = 0.5

#: Search radius and Gaussian sigma expressed as multiples of the cell size. At 0.5 deg
#: these reproduce the historical literals 0.5 and 0.25 exactly.
SEARCH_RADIUS_CELLS = 1.0
SIGMA_CELLS = 0.5

#: A grid within this tolerance of 0.5 deg takes the legacy path bit-for-bit.
_LEGACY_TOL = 1e-6


def grid_cell_size(obj) -> float:
    """Cell size in degrees, inferred from a dataset's / array's own lat-lon coordinates.

    Accepts anything carrying `lat` and `lon` coordinates (Dataset, DataArray). Falls back
    to CELL_SIZE when the coordinates are missing or degenerate, so callers that pass a
    bare array keep working.

    Raises ValueError on an irregular grid: silently averaging an irregular spacing would
    produce a search geometry that is wrong everywhere rather than obviously broken.
    """
    try:
        lat = np.asarray(obj["lat"].values, dtype=float)
        lon = np.asarray(obj["lon"].values, dtype=float)
    except (KeyError, TypeError, IndexError):
        return CELL_SIZE
    if lat.size < 2 or lon.size < 2:
        return CELL_SIZE
    dlat = np.abs(np.diff(lat))
    dlon = np.abs(np.diff(lon))
    for name, d in (("lat", dlat), ("lon", dlon)):
        if not np.allclose(d, d[0], rtol=0, atol=1e-6):
            raise ValueError(
                f"{name} spacing is irregular (min {d.min():.6f}, max {d.max():.6f}); "
                "extraction geometry is only defined on a regular grid."
            )
    size = float(dlat[0])
    if abs(size - float(dlon[0])) > 1e-6:
        raise ValueError(
            f"non-square cells: lat spacing {size:.6f} vs lon spacing {float(dlon[0]):.6f}"
        )
    # Snap a grid that is 0.5 deg to within floating-point noise onto the exact legacy
    # constant, so inference can never perturb a shipped layer's geometry.
    return CELL_SIZE if abs(size - CELL_SIZE) < _LEGACY_TOL else size


def search_radius_for(obj) -> float:
    """Point-extraction search radius in degrees for this object's grid (0.5 at 0.5 deg)."""
    return SEARCH_RADIUS_CELLS * grid_cell_size(obj)


def sigma_for(obj) -> float:
    """Gaussian sigma in degrees for this object's grid (0.25 at 0.5 deg)."""
    return SIGMA_CELLS * grid_cell_size(obj)

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


def as_period_dataset(ds: xr.Dataset) -> xr.Dataset:
    """Give an observational layer the single `decade` slice the extractor expects.

    Every projected layer is (decade, lat, lon). An `observational-historical-v1` layer
    is (lat, lon): it summarises ONE observed window, so it has no decade axis and one
    is deliberately not written to disk (faking a time axis on a file that has no time
    information is the kind of drift this repo forbids elsewhere).

    Rather than branch the whole extraction path, the axis is added AT READ TIME as a
    length-1 dimension labelled with the window's FIRST YEAR, taken from the file's own
    `temporal_window` attribute. Layers already carrying `decade` are returned untouched,
    so no projected layer changes behaviour.

    THE LABEL IS A PERIOD START, NOT A DECADE. For the tornado layers it is 1950,
    meaning "the 1950-2025 record", NOT "the 1950s". The pairing is what disambiguates
    it: these rows carry scenario `observed`, which no projected layer uses, and
    layers.csv carries the layer's `temporal_window` verbatim. Do not compute a decadal
    change from a layer whose scenario is `observed` -- there is only one period.
    """
    if "decade" in ds.dims:
        return ds
    window = str(ds.attrs.get("temporal_window", "")).strip()
    try:
        start = int(window.split("-")[0])
    except (ValueError, IndexError):
        raise ValueError(
            f"Layer has no decade dimension and no parseable `temporal_window` attribute "
            f"(got {window!r}); cannot place its values on a period axis."
        )
    return ds.expand_dims(decade=[start])


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
    search_radius: Optional[float] = None,
    sigma: Optional[float] = None,
) -> Dict[str, Dict[int, float]]:
    """Extract values at a point using Gaussian distance-weighted averaging.

    Weights nearby cell CENTRES by exp(-0.5 (d/sigma)^2) in DEGREE space, normalizes, sums.
    NaN cells are dropped and the surviving weights renormalized.

    Measured behaviour on the standard 0.5 deg grid with the defaults below -- see
    ASSET-CATALOG.md "Spatial averaging" for the full record, including the compounding with
    per-layer processing-time smoothing:

    - It is a 4-CELL BLEND, not a 3x3 smoother. The 3x3 stencil requires the site to sit
      exactly on a cell centre; 0.0005 deg off collapses it to 2x2, and 100% of 20,000
      random sites got 2x2. Footprint is 1 deg x 1 deg, ~111 km N-S.
    - The blend ranges from a flat 0.25/0.25/0.25/0.25 mean (site on a grid vertex) to
      0.78/0.10/0.10/0.01 (site at a cell centre), so sub-cell position changes the kernel's
      character, not just its numbers.
    - radius = 2*sigma, so the Gaussian is TRUNCATED at 13.5% of peak rather than tapered.
    - Longitude separation WRAPS the antimeridian (fixed 2026-08-12; before that 180E and
      180W returned different values). Longitude is NOT cos(lat)-scaled, so the ground
      footprint stretches with latitude -- 2.6x wider E-W than N-S at 67 deg N.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        lat: Target latitude (decimal degrees)
        lon: Target longitude (decimal degrees, -180 to 180)
        variables: List of variables to extract (default: all numeric variables)
        search_radius: Search radius in degrees. None (default) infers it from `ds`'s own
            grid as SEARCH_RADIUS_CELLS x cell size, which is 0.5 on a 0.5 deg grid --
            identical to the previous hardcoded default. Pass a number to override.
        sigma: Gaussian sigma in degrees. None (default) infers SIGMA_CELLS x cell size,
            which is 0.25 on a 0.5 deg grid -- identical to the previous default.

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Example:
        >>> ds = xr.open_dataset("processed.nc")
        >>> result = extract_by_point(ds, lat=45.5, lon=-122.7)
        >>> print(result["median"][2050])  # Value for 2050 decade
    """
    lon = normalize_longitude(lon)

    # Geometry follows THIS dataset's grid. On the 0.5 deg grid every shipped layer uses,
    # these resolve to exactly 0.5 and 0.25 -- the previous hardcoded defaults.
    if search_radius is None:
        search_radius = search_radius_for(ds)
    if sigma is None:
        sigma = sigma_for(ds)

    lats = ds.lat.values
    lons = ds.lon.values

    # Longitude separation WRAPS at the antimeridian. Without the wrap the search window
    # is one-sided at the seam and 180 deg E returns a different value from 180 deg W --
    # the same meridian. Measured on burntarea ssp585 at 17S, 2090s: 0.775 vs 0.962.
    lat_delta = np.abs(lats - lat)
    lon_delta = np.abs(lons - lon)
    lon_delta = np.minimum(lon_delta, 360.0 - lon_delta)

    lat_mask = lat_delta <= search_radius
    lon_mask = lon_delta <= search_radius

    nearby_lats = lats[lat_mask]
    nearby_lons = lons[lon_mask]

    if len(nearby_lats) == 0 or len(nearby_lons) == 0:
        raise ValueError(
            f"No grid cells found within {search_radius} degrees of ({lat}, {lon})"
        )

    # Distances from the target point to cell CENTRES, in degrees. Longitude uses the
    # wrapped separation computed above rather than a raw difference.
    lon_delta_grid, lat_grid = np.meshgrid(lon_delta[lon_mask], nearby_lats)
    distances = np.sqrt((lat_grid - lat) ** 2 + lon_delta_grid ** 2)

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
    cell_size: Optional[float] = None,
) -> Dict[str, Dict[int, float]]:
    """Extract values within a polygon using area-weighted averaging.

    Calculates the fraction of each grid cell that falls inside the polygon
    and uses these fractions as weights for averaging.

    Args:
        ds: xarray Dataset with (decade, lat, lon) dimensions
        polygon: Shapely Polygon or MultiPolygon defining the region
        variables: List of variables to extract (default: all numeric variables)
        cell_size: Grid cell size in degrees. None (default) infers it from `ds`'s own
            coordinates -- 0.5 on every layer shipped before 2026-08-14, so the cell
            boundaries and intersection fractions are unchanged for those.

    Returns:
        Dict mapping variable names to {decade: value} dicts

    Raises:
        ValueError: If polygon doesn't intersect any grid cells

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-123, 45), (-122, 45), (-122, 46), (-123, 46)])
        >>> result = extract_by_polygon(ds, poly)
    """
    if cell_size is None:
        cell_size = grid_cell_size(ds)

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
