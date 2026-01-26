"""Natural Earth and GADM boundary loading utilities.

Load and cache country/state/province/county boundaries from:
- Natural Earth (countries, states/provinces)
- GADM (counties/districts and finer admin levels)

Data sources:
- Natural Earth: https://www.naturalearthdata.com/
- GADM: https://gadm.org/
"""

from difflib import SequenceMatcher
from pathlib import Path
from typing import List, Optional, Tuple, Union

import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon


# Default cache directories
DEFAULT_CACHE_DIR = Path("data/naturalearth")
GADM_CACHE_DIR = Path("data/gadm")

# ISO3 country codes for GADM downloads
COUNTRY_ISO3 = {
    "united states": "USA",
    "usa": "USA",
    "us": "USA",
    "new zealand": "NZL",
    "nz": "NZL",
    "canada": "CAN",
    "australia": "AUS",
    "united kingdom": "GBR",
    "uk": "GBR",
    "germany": "DEU",
    "france": "FRA",
    "brazil": "BRA",
    "china": "CHN",
    "india": "IND",
    "japan": "JPN",
    "mexico": "MEX",
    "indonesia": "IDN",
    "russia": "RUS",
    "south africa": "ZAF",
}

# GADM column names for different admin levels
GADM_NAME_COLUMNS = {
    1: ["NAME_1", "VARNAME_1", "NL_NAME_1"],  # States/Provinces
    2: ["NAME_2", "VARNAME_2", "NL_NAME_2"],  # Counties/Districts
    3: ["NAME_3", "VARNAME_3", "NL_NAME_3"],  # Sub-districts
}

# Natural Earth S3 URLs for cultural vectors
NE_URLS = {
    "countries_110m": "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip",
    "countries_50m": "https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_0_countries.zip",
    "countries_10m": "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_countries.zip",
    "states_50m": "https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_1_states_provinces.zip",
    "states_10m": "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip",
}

# Column names to search for region names (in priority order)
COUNTRY_NAME_COLUMNS = ["NAME", "ADMIN", "NAME_LONG", "FORMAL_EN", "SOVEREIGNT"]
STATE_NAME_COLUMNS = ["name", "NAME", "name_en", "admin", "gn_name"]


def load_naturalearth(
    scale: str = "50m",
    layer: str = "countries",
    cache_dir: Optional[Path] = None,
) -> gpd.GeoDataFrame:
    """Load Natural Earth data with parquet caching for fast reloads.

    Downloads data on first use and caches as parquet for fast subsequent loads.

    Args:
        scale: Map scale - "110m" (coarse), "50m" (medium), "10m" (detailed)
        layer: Layer type - "countries" or "states"
        cache_dir: Cache directory (default: data/naturalearth/)

    Returns:
        GeoDataFrame with geometry and attribute columns

    Raises:
        ValueError: If invalid scale/layer combination
    """
    cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    key = f"{layer}_{scale}"
    cache_file = cache_dir / f"ne_{scale}_{layer}.parquet"

    if cache_file.exists():
        return gpd.read_parquet(cache_file)

    url = NE_URLS.get(key)
    if url is None:
        available = [k for k in NE_URLS.keys() if layer in k]
        raise ValueError(
            f"Unknown layer/scale: {key}. "
            f"Available for {layer}: {available}"
        )

    print(f"Downloading Natural Earth {scale} {layer}...")
    print(f"  URL: {url}")

    gdf = gpd.read_file(url)

    # Cache as parquet for faster reloads
    gdf.to_parquet(cache_file)
    print(f"  Cached to: {cache_file}")

    return gdf


def fuzzy_match_region(
    query: str,
    gdf: gpd.GeoDataFrame,
    name_columns: Optional[List[str]] = None,
    threshold: float = 0.6,
) -> Tuple[Optional[Union[Polygon, MultiPolygon]], Optional[str], float]:
    """Find a region by fuzzy name matching.

    Searches multiple name columns for the best match.

    Args:
        query: Search string (e.g., "United States", "California", "Brasil")
        gdf: GeoDataFrame with region geometries
        name_columns: Columns to search (default: auto-detect)
        threshold: Minimum similarity score (0-1, default: 0.6)

    Returns:
        Tuple of (geometry, matched_name, score) or (None, None, 0.0) if no match
    """
    if name_columns is None:
        # Auto-detect based on available columns
        name_columns = [c for c in COUNTRY_NAME_COLUMNS + STATE_NAME_COLUMNS if c in gdf.columns]

    query_lower = query.lower().strip()
    best_match_idx = None
    best_score = 0.0
    best_name = None

    for col in name_columns:
        if col not in gdf.columns:
            continue

        for idx, name in gdf[col].items():
            if name is None or str(name).strip() == "":
                continue

            name_str = str(name)
            name_lower = name_str.lower()

            # Exact match (case-insensitive)
            if name_lower == query_lower:
                return gdf.loc[idx, "geometry"], name_str, 1.0

            # Fuzzy match using SequenceMatcher
            score = SequenceMatcher(None, query_lower, name_lower).ratio()

            # Also check if query is contained in name (partial match)
            if query_lower in name_lower:
                score = max(score, 0.8)

            if score > best_score and score >= threshold:
                best_score = score
                best_match_idx = idx
                best_name = name_str

    if best_match_idx is not None:
        return gdf.loc[best_match_idx, "geometry"], best_name, best_score

    return None, None, 0.0


def get_region_geometry(
    name: str,
    region_type: str = "country",
    scale: str = "50m",
) -> Optional[Union[Polygon, MultiPolygon]]:
    """Get geometry for a named region with fuzzy matching.

    Args:
        name: Region name (fuzzy matched)
        region_type: "country" or "state"
        scale: Natural Earth scale ("50m" or "10m")

    Returns:
        Shapely geometry or None if not found

    Example:
        >>> geometry = get_region_geometry("United States")
        >>> geometry = get_region_geometry("California", region_type="state")
        >>> geometry = get_region_geometry("Brasil")  # Fuzzy matches "Brazil"
    """
    layer = "countries" if region_type == "country" else "states"

    # For states, 110m scale is not available
    if layer == "states" and scale == "110m":
        scale = "50m"

    name_columns = (
        COUNTRY_NAME_COLUMNS if region_type == "country" else STATE_NAME_COLUMNS
    )

    gdf = load_naturalearth(scale=scale, layer=layer)
    geometry, matched_name, score = fuzzy_match_region(gdf=gdf, query=name, name_columns=name_columns)

    if geometry is not None and matched_name != name:
        print(f"  Matched '{name}' to '{matched_name}' (score: {score:.2f})")

    return geometry


def search_regions(
    query: str,
    region_type: str = "country",
    scale: str = "50m",
    limit: int = 10,
) -> List[str]:
    """Search for regions matching a query string.

    Useful for interactive exploration when exact name is unknown.

    Args:
        query: Search string (partial match)
        region_type: "country" or "state"
        scale: Natural Earth scale
        limit: Maximum number of results

    Returns:
        List of matching region names (sorted alphabetically)

    Example:
        >>> search_regions("cal")
        ['California', 'Calabria', ...]
    """
    layer = "countries" if region_type == "country" else "states"

    # For states, 110m scale is not available
    if layer == "states" and scale == "110m":
        scale = "50m"

    gdf = load_naturalearth(scale=scale, layer=layer)

    # Search in the primary name column
    name_col = "NAME" if "NAME" in gdf.columns else "name"
    if name_col not in gdf.columns:
        name_col = gdf.columns[0]  # Fallback

    matches = []
    query_lower = query.lower()

    for name in gdf[name_col].dropna():
        name_str = str(name)
        if query_lower in name_str.lower():
            matches.append(name_str)

    return sorted(set(matches))[:limit]


def list_all_regions(
    region_type: str = "country",
    scale: str = "50m",
) -> List[str]:
    """List all available region names.

    Args:
        region_type: "country" or "state"
        scale: Natural Earth scale

    Returns:
        Sorted list of all region names
    """
    layer = "countries" if region_type == "country" else "states"

    if layer == "states" and scale == "110m":
        scale = "50m"

    gdf = load_naturalearth(scale=scale, layer=layer)

    name_col = "NAME" if "NAME" in gdf.columns else "name"
    if name_col not in gdf.columns:
        name_col = gdf.columns[0]

    return sorted(gdf[name_col].dropna().astype(str).unique().tolist())


def get_region_bounds(
    name: str,
    region_type: str = "country",
    scale: str = "50m",
) -> Optional[Tuple[float, float, float, float]]:
    """Get bounding box for a named region.

    Args:
        name: Region name
        region_type: "country" or "state"
        scale: Natural Earth scale

    Returns:
        Tuple of (minx, miny, maxx, maxy) or None if not found
    """
    geometry = get_region_geometry(name, region_type, scale)
    if geometry is None:
        return None
    return geometry.bounds


def get_region_centroid(
    name: str,
    region_type: str = "country",
    scale: str = "50m",
) -> Optional[Tuple[float, float]]:
    """Get centroid (center point) for a named region.

    Args:
        name: Region name
        region_type: "country" or "state"
        scale: Natural Earth scale

    Returns:
        Tuple of (lon, lat) or None if not found
    """
    geometry = get_region_geometry(name, region_type, scale)
    if geometry is None:
        return None
    centroid = geometry.centroid
    return (centroid.x, centroid.y)


# =============================================================================
# GADM (Global Administrative Areas Database) Support
# =============================================================================


def get_country_iso3(country: str) -> Optional[str]:
    """Get ISO3 country code for GADM download.

    Args:
        country: Country name (case-insensitive)

    Returns:
        ISO3 code or None if not found
    """
    country_lower = country.lower().strip()

    # Direct lookup
    if country_lower in COUNTRY_ISO3:
        return COUNTRY_ISO3[country_lower]

    # Fuzzy match
    for key, code in COUNTRY_ISO3.items():
        if country_lower in key or key in country_lower:
            return code

    return None


def load_gadm(
    country: str,
    admin_level: int = 2,
    cache_dir: Optional[Path] = None,
) -> gpd.GeoDataFrame:
    """Load GADM administrative boundaries for a country.

    Downloads GADM GeoPackage on first use and caches locally.
    GADM provides boundaries down to admin level 2-4 (varies by country).

    Args:
        country: Country name or ISO3 code
        admin_level: Administrative level (1=state, 2=county, 3=sub-county)
        cache_dir: Cache directory (default: data/gadm/)

    Returns:
        GeoDataFrame with administrative boundaries

    Raises:
        ValueError: If country not found or admin level not available
    """
    cache_dir = Path(cache_dir) if cache_dir else GADM_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Get ISO3 code
    iso3 = country.upper() if len(country) == 3 else get_country_iso3(country)
    if iso3 is None:
        raise ValueError(
            f"Unknown country: {country}. "
            f"Known countries: {list(COUNTRY_ISO3.keys())}"
        )

    # Check for cached parquet
    cache_file = cache_dir / f"gadm41_{iso3}_level{admin_level}.parquet"
    if cache_file.exists():
        return gpd.read_parquet(cache_file)

    # Download from GADM
    gpkg_cache = cache_dir / f"gadm41_{iso3}.gpkg"

    if not gpkg_cache.exists():
        url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_{iso3}.gpkg"
        print(f"Downloading GADM data for {iso3}...")
        print(f"  URL: {url}")
        print("  (This may take a few minutes for large countries)")

        import urllib.request
        urllib.request.urlretrieve(url, gpkg_cache)
        print(f"  Downloaded to: {gpkg_cache}")

    # Load the requested admin level
    layer_name = f"ADM_ADM_{admin_level}"
    try:
        gdf = gpd.read_file(gpkg_cache, layer=layer_name)
    except Exception:
        # Try alternative layer naming
        try:
            gdf = gpd.read_file(gpkg_cache, layer=f"gadm41_{iso3}_{admin_level}")
        except Exception as e:
            raise ValueError(
                f"Admin level {admin_level} not available for {iso3}. "
                f"Error: {e}"
            )

    # Cache as parquet for faster reloads
    gdf.to_parquet(cache_file)
    print(f"  Cached level {admin_level} to: {cache_file}")

    return gdf


def search_gadm_regions(
    query: str,
    country: str,
    admin_level: int = 2,
    limit: int = 10,
) -> List[str]:
    """Search for regions in GADM data.

    Args:
        query: Search string (partial match)
        country: Country name or ISO3 code
        admin_level: Administrative level to search
        limit: Maximum number of results

    Returns:
        List of matching region names

    Example:
        >>> search_gadm_regions("shasta", "USA", admin_level=2)
        ['Shasta']
    """
    gdf = load_gadm(country, admin_level)

    # Get name columns for this level
    name_cols = GADM_NAME_COLUMNS.get(admin_level, [f"NAME_{admin_level}"])

    matches = []
    query_lower = query.lower()

    for col in name_cols:
        if col not in gdf.columns:
            continue
        for name in gdf[col].dropna():
            name_str = str(name)
            if query_lower in name_str.lower():
                matches.append(name_str)

    return sorted(set(matches))[:limit]


def get_gadm_geometry(
    name: str,
    country: str,
    admin_level: int = 2,
    parent_name: Optional[str] = None,
) -> Optional[Union[Polygon, MultiPolygon]]:
    """Get geometry for a GADM administrative region.

    Args:
        name: Region name (fuzzy matched)
        country: Country name or ISO3 code
        admin_level: Administrative level (1=state, 2=county)
        parent_name: Parent region name for disambiguation (e.g., state name)

    Returns:
        Shapely geometry or None if not found

    Example:
        >>> geometry = get_gadm_geometry("Shasta", "USA", admin_level=2, parent_name="California")
    """
    gdf = load_gadm(country, admin_level)

    # Filter by parent if provided
    if parent_name and admin_level > 1:
        parent_col = f"NAME_{admin_level - 1}"
        if parent_col in gdf.columns:
            parent_lower = parent_name.lower()
            mask = gdf[parent_col].str.lower().str.contains(parent_lower, na=False)
            if mask.any():
                gdf = gdf[mask]

    # Get name columns for this level
    name_cols = GADM_NAME_COLUMNS.get(admin_level, [f"NAME_{admin_level}"])

    geometry, matched_name, score = fuzzy_match_region(
        query=name,
        gdf=gdf,
        name_columns=name_cols,
    )

    if geometry is not None and matched_name.lower() != name.lower():
        print(f"  Matched '{name}' to '{matched_name}' (score: {score:.2f})")

    return geometry


def list_gadm_regions(
    country: str,
    admin_level: int = 2,
    parent_name: Optional[str] = None,
) -> List[str]:
    """List all regions at a given admin level.

    Args:
        country: Country name or ISO3 code
        admin_level: Administrative level
        parent_name: Filter by parent region

    Returns:
        Sorted list of region names
    """
    gdf = load_gadm(country, admin_level)

    # Filter by parent if provided
    if parent_name and admin_level > 1:
        parent_col = f"NAME_{admin_level - 1}"
        if parent_col in gdf.columns:
            parent_lower = parent_name.lower()
            mask = gdf[parent_col].str.lower().str.contains(parent_lower, na=False)
            if mask.any():
                gdf = gdf[mask]

    name_col = f"NAME_{admin_level}"
    if name_col not in gdf.columns:
        name_col = gdf.columns[0]

    return sorted(gdf[name_col].dropna().astype(str).unique().tolist())
