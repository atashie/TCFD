#!/usr/bin/env python
"""Extract region polygon boundaries and update ex-locations-timber.csv.

Extracts WKT polygon boundaries from GADM/Natural Earth for region rows
that have the Region column populated but no Polygon data.
"""

import sys
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
from shapely.geometry import Point, Polygon, MultiPolygon

from scripts.utils.natural_earth import get_gadm_geometry, get_region_geometry


def simplify_geometry(geom, tolerance=0.01):
    """Simplify geometry to reduce WKT size.

    Args:
        geom: Shapely geometry
        tolerance: Simplification tolerance in degrees (~1km at equator)

    Returns:
        Simplified geometry
    """
    if geom is None:
        return None
    return geom.simplify(tolerance, preserve_topology=True)


def geometry_to_wkt(geom, reference_point=None):
    """Convert Shapely geometry to WKT string.

    For MultiPolygon:
    - If reference_point provided, selects polygon containing that point
    - Otherwise, selects the largest polygon

    Args:
        geom: Shapely Polygon or MultiPolygon
        reference_point: Optional (lon, lat) tuple for point containment check

    Returns:
        WKT string
    """
    if geom is None:
        return None

    # If MultiPolygon, select appropriate polygon
    if isinstance(geom, MultiPolygon):
        if reference_point:
            pt = Point(reference_point)
            for poly in geom.geoms:
                if poly.contains(pt):
                    geom = poly
                    break
            else:
                # Fallback to largest if point not contained in any polygon
                print(f"    Warning: reference point not contained, using largest polygon")
                geom = max(geom.geoms, key=lambda p: p.area)
        else:
            # Default: use largest polygon
            geom = max(geom.geoms, key=lambda p: p.area)

    return geom.wkt


def extract_region_polygon(region_name: str) -> str:
    """Extract polygon WKT for a region name.

    Args:
        region_name: Region identifier (e.g., "Shasta County, CA", "Louisiana")

    Returns:
        WKT polygon string, or None if not found
    """
    region_name = region_name.strip()
    geometry = None
    reference_point = None  # Optional point for MultiPolygon selection

    # US County format: "Name County, ST"
    if "County" in region_name and "," in region_name:
        parts = region_name.split(",")
        county_name = parts[0].strip().replace(" County", "")
        state_abbr = parts[1].strip()

        # Map state abbreviations to full names
        state_map = {
            "CA": "California",
            "OR": "Oregon",
            "WA": "Washington",
            "LA": "Louisiana",
            "TX": "Texas",
            "FL": "Florida",
        }
        state_name = state_map.get(state_abbr, state_abbr)

        print(f"  Extracting GADM county: {county_name}, {state_name}")
        geometry = get_gadm_geometry(
            name=county_name,
            country="USA",
            admin_level=2,
            parent_name=state_name
        )

    # Handle misspelled Louisiana
    elif region_name.lower() in ("lousiana", "louisiana"):
        print(f"  Extracting Natural Earth state: Louisiana")
        geometry = get_region_geometry("Louisiana", region_type="state")

    # New Zealand special case - use reference point to select North Island
    elif "new zealand" in region_name.lower():
        print(f"  Extracting Natural Earth country: New Zealand (selecting North Island)")
        geometry = get_region_geometry("New Zealand", region_type="country")
        # Reference point in Central North Island to select correct polygon
        reference_point = (176.1, -38.8)

    # Try state first, then country
    else:
        print(f"  Trying Natural Earth state: {region_name}")
        geometry = get_region_geometry(region_name, region_type="state")
        if geometry is None:
            print(f"  Trying Natural Earth country: {region_name}")
            geometry = get_region_geometry(region_name, region_type="country")

    if geometry is None:
        print(f"  WARNING: Could not find geometry for '{region_name}'")
        return None

    # Simplify to reduce file size
    geometry = simplify_geometry(geometry, tolerance=0.005)

    return geometry_to_wkt(geometry, reference_point=reference_point)


def main():
    """Main function to update CSV with region polygons."""
    # Paths
    project_root = Path(__file__).parent.parent
    csv_path = project_root / "location-analyses" / "ex-locations-timber.csv"

    print(f"Loading: {csv_path}")

    # Load CSV (tab-separated)
    df = pd.read_csv(csv_path, sep="\t")
    print(f"Loaded {len(df)} rows")
    print(f"Columns: {list(df.columns)}")

    # Find region rows that need polygon extraction
    updates = 0
    for idx, row in df.iterrows():
        region = row.get("Region", "")
        polygon = row.get("Polygon", "")

        # Skip if no region or already has polygon
        if pd.isna(region) or not str(region).strip():
            continue
        if pd.notna(polygon) and str(polygon).strip():
            continue

        location_name = row.get("Location Name", "")
        print(f"\nProcessing row {idx}: {location_name}")
        print(f"  Region: {region}")

        # Extract polygon
        wkt_polygon = extract_region_polygon(str(region))

        if wkt_polygon:
            df.at[idx, "Polygon"] = wkt_polygon
            updates += 1
            print(f"  SUCCESS: Added polygon ({len(wkt_polygon)} chars)")
        else:
            print(f"  FAILED: No polygon extracted")

    if updates > 0:
        # Save updated CSV
        df.to_csv(csv_path, sep="\t", index=False)
        print(f"\nUpdated {updates} rows in {csv_path}")
    else:
        print("\nNo updates needed")

    # Show summary
    print("\n=== Summary ===")
    for idx, row in df.iterrows():
        location = row.get("Location Name", "")
        has_polygon = pd.notna(row.get("Polygon")) and str(row.get("Polygon", "")).strip()
        has_region = pd.notna(row.get("Region")) and str(row.get("Region", "")).strip()
        status = "POLYGON" if has_polygon else ("REGION" if has_region else "POINT/OTHER")
        print(f"  {location}: {status}")


if __name__ == "__main__":
    main()
