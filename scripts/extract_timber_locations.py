#!/usr/bin/env python
"""Extract timber data for all locations in ex-locations-timber.csv.

Processes locations with species-specific datasets:
- Evergreen -> evgndltr cveg (SSP scenarios)
- Loblolly Pine -> needleleaf-evergreen-tree-temperate cveg (RCP scenarios)
- Deciduous Birch -> tebrsu cveg (RCP scenarios)
"""

import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr
from shapely import wkt
from shapely.geometry import Polygon, MultiPolygon

# Add project root to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.spatial_extract import (
    extract_by_point,
    extract_by_polygon,
    apply_percentile_inversion,
    normalize_longitude,
)
from scripts.utils.export_formatter import (
    LocationMetadata,
    HazardMapping,
    ExtractionResult,
    format_extraction_results,
    EXPORT_COLUMNS,
)
from scripts.utils.natural_earth import get_gadm_geometry, get_region_geometry


# Dataset mapping by species
SPECIES_DATASET_MAPPING = {
    "Evergreen": {
        "folder": "evgndltr_npp-cveg",
        "variable_pattern": "cveg-evgndltr",
        "scenarios": ["ssp126", "ssp370", "ssp585"],
        "hazard_measure": "Evergreen Needleleaf Vegetation Carbon (kg/m2)",
    },
    "Loblolly Pine": {
        "folder": "loblolly-temperate_cveg-needleleaf-evergreen-tree-temperate_annual",
        "variable_pattern": "cveg-needleleaf-evergreen-tree-temperate",
        "scenarios": ["rcp26", "rcp60", "rcp85"],
        "hazard_measure": "Temperate Needleleaf Evergreen Vegetation Carbon (kg/m2)",
    },
    "Deciduos Birch": {  # Note: misspelled in source CSV
        "folder": "tebrsu_cveg_annual",
        "variable_pattern": "cveg-tebrsu",
        "scenarios": ["rcp26", "rcp60", "rcp85"],
        "hazard_measure": "Temperate Broadleaf Summer Deciduous Vegetation Carbon (kg/m2)",
    },
}

# Scenario recoding to human-readable labels
SCENARIO_LABELS = {
    "ssp126": "Low Emissions",
    "ssp370": "Middle of the Road",
    "ssp585": "High Emissions",
    "rcp26": "Low Emissions",
    "rcp60": "Middle of the Road",
    "rcp85": "High Emissions",
}


def normalize_polygon_coords(wkt_str: str) -> Polygon:
    """Parse WKT polygon and normalize longitude coordinates.

    Adds 360 to any longitude < -180 to bring into -180 to 180 range.

    Args:
        wkt_str: WKT POLYGON string

    Returns:
        Shapely Polygon with normalized coordinates
    """
    # Parse WKT
    geom = wkt.loads(wkt_str)

    # Normalize coordinates
    if isinstance(geom, Polygon):
        exterior_coords = []
        for x, y in geom.exterior.coords:
            # Normalize longitude
            while x < -180:
                x += 360
            while x > 180:
                x -= 360
            exterior_coords.append((x, y))
        return Polygon(exterior_coords)
    elif isinstance(geom, MultiPolygon):
        normalized_polys = []
        for poly in geom.geoms:
            exterior_coords = []
            for x, y in poly.exterior.coords:
                while x < -180:
                    x += 360
                while x > 180:
                    x -= 360
                exterior_coords.append((x, y))
            normalized_polys.append(Polygon(exterior_coords))
        return MultiPolygon(normalized_polys)

    return geom


def load_processed_dataset(
    processed_dir: Path,
    folder_name: str,
    variable_pattern: str,
    scenario: str,
) -> xr.Dataset:
    """Load a processed NetCDF file.

    Args:
        processed_dir: Base directory for processed data
        folder_name: Name of the dataset folder
        variable_pattern: Variable name pattern for file matching
        scenario: Scenario code (e.g., "ssp126", "rcp26")

    Returns:
        xarray Dataset
    """
    folder_path = processed_dir / folder_name

    # Find matching file
    file_pattern = f"{variable_pattern}_{scenario}_processed.nc"
    files = list(folder_path.glob(file_pattern))

    if not files:
        raise FileNotFoundError(
            f"No file matching '{file_pattern}' in {folder_path}"
        )

    return xr.open_dataset(files[0])


def extract_for_location(
    ds: xr.Dataset,
    location_name: str,
    lat: Optional[float],
    lon: Optional[float],
    polygon_wkt: Optional[str],
    region_name: Optional[str],
    species: str,
    scenario: str,
    hazard_measure: str,
) -> List[ExtractionResult]:
    """Extract data for a single location using all applicable methods.

    Returns one ExtractionResult per extraction type.
    """
    results = []

    # Check if percentile inversion was already applied in the NetCDF
    # If percentile_direction attribute exists and is "higher_is_better",
    # the percentiles are already inverted - don't invert again
    already_inverted = ds.attrs.get("percentile_direction") == "higher_is_better"

    # Create hazard mapping
    hazard = HazardMapping(
        folder_pattern=species.lower(),
        hazard="Timber",
        hazard_measure=hazard_measure,
        percentile_direction="higher_is_better",  # cveg is higher = better
    )

    # Recode scenario to human-readable label
    scenario_label = SCENARIO_LABELS.get(scenario, scenario)

    # Point extraction
    if pd.notna(lat) and pd.notna(lon):
        try:
            lon_norm = normalize_longitude(lon)
            data = extract_by_point(ds, lat, lon_norm)
            # Apply percentile inversion for "higher_is_better" only if not already done
            if not already_inverted:
                data = apply_percentile_inversion(data, "higher_is_better")

            meta = LocationMetadata(
                location=f"{location_name} (point)",
                lat=lat,
                lon=lon_norm,
                asset_type="Timber",
                sub_asset_unit=species,
            )
            results.append(ExtractionResult(
                location_metadata=meta,
                hazard_mapping=hazard,
                scenario=scenario_label,
                data=data,
            ))
        except Exception as e:
            print(f"  Warning: Point extraction failed for {location_name}: {e}")

    # Polygon extraction
    if pd.notna(polygon_wkt) and polygon_wkt.strip():
        try:
            polygon = normalize_polygon_coords(polygon_wkt)
            data = extract_by_polygon(ds, polygon)
            if not already_inverted:
                data = apply_percentile_inversion(data, "higher_is_better")

            # Get centroid for lat/lon
            centroid = polygon.centroid

            meta = LocationMetadata(
                location=f"{location_name} (poly)",
                centroid_lat=centroid.y,
                centroid_lon=centroid.x,
                asset_type="Timber",
                sub_asset_unit=species,
            )
            results.append(ExtractionResult(
                location_metadata=meta,
                hazard_mapping=hazard,
                scenario=scenario_label,
                data=data,
            ))
        except Exception as e:
            print(f"  Warning: Polygon extraction failed for {location_name}: {e}")

    # Region extraction
    if pd.notna(region_name) and region_name.strip():
        try:
            geometry = None
            clean_region = region_name.strip()

            # Parse region format
            if "County" in region_name and "," in region_name:
                # e.g., "Shasta County, CA" -> GADM county extraction
                parts = region_name.split(",")
                county_part = parts[0].strip().replace(" County", "")
                state_abbr = parts[1].strip() if len(parts) > 1 else None

                # Map state abbreviations to full names
                state_map = {"CA": "California", "OR": "Oregon", "WA": "Washington"}
                state_name = state_map.get(state_abbr, state_abbr)

                print(f"    Using GADM for {county_part} County, {state_name}")
                geometry = get_gadm_geometry(
                    county_part, "USA", admin_level=2, parent_name=state_name
                )
                clean_region = f"{county_part} County, {state_name}"

            elif region_name.strip().lower() == "lousiana":
                # Fix misspelling - use Natural Earth for US states
                clean_region = "Louisiana"
                geometry = get_region_geometry("Louisiana", "state")

            elif "Louisiana" in region_name:
                # US state
                geometry = get_region_geometry("Louisiana", "state")
                clean_region = "Louisiana"

            elif "New Zealand" in region_name:
                # For NZ regions not in Natural Earth, use country-level
                # Central North Island is not an admin region
                print(f"    Note: '{region_name}' not an admin region, using NZ country")
                geometry = get_region_geometry("New Zealand", "country")
                clean_region = "New Zealand"

            else:
                # Try Natural Earth state first, then country
                geometry = get_region_geometry(clean_region, "state")
                if geometry is None:
                    geometry = get_region_geometry(clean_region, "country")

            if geometry is None:
                raise ValueError(f"Region not found: {region_name}")

            data = extract_by_polygon(ds, geometry)
            if not already_inverted:
                data = apply_percentile_inversion(data, "higher_is_better")

            meta = LocationMetadata(
                location=f"{location_name} (region)",
                region=clean_region,
                asset_type="Timber",
                sub_asset_unit=species,
            )
            results.append(ExtractionResult(
                location_metadata=meta,
                hazard_mapping=hazard,
                scenario=scenario_label,
                data=data,
            ))
        except Exception as e:
            print(f"  Warning: Region extraction failed for {location_name} ({region_name}): {e}")

    return results


def main():
    """Main extraction workflow."""
    # Paths
    project_root = Path(__file__).parent.parent
    locations_csv = project_root / "location-analyses" / "ex-locations-timber.csv"
    processed_dir = project_root / "data" / "processed"
    exports_dir = project_root / "data" / "exports"

    # Create exports directory if needed
    exports_dir.mkdir(parents=True, exist_ok=True)

    # Load locations
    print(f"Loading locations from {locations_csv}")
    df = pd.read_csv(locations_csv, sep="\t")  # Tab-separated
    print(f"Found {len(df)} location rows")
    print(f"Columns: {list(df.columns)}")

    # Collect all extraction results
    all_results: List[ExtractionResult] = []

    # Process each location row
    for idx, row in df.iterrows():
        location_name = row.get("Location Name", "")
        if not location_name or pd.isna(location_name):
            continue

        # Get base location name (remove extraction type suffix if present)
        base_name = re.sub(r'\s*\((point|poly|region)\)\s*$', '', location_name)

        species = row.get("Species of Interest", "")
        if pd.isna(species) or not species.strip():
            continue

        species = species.strip()

        # Get dataset config for this species
        if species not in SPECIES_DATASET_MAPPING:
            print(f"Warning: Unknown species '{species}' for {base_name}, skipping")
            continue

        config = SPECIES_DATASET_MAPPING[species]

        print(f"\nProcessing: {base_name} ({species})")

        # Get location details
        lat = row.get("Latitude")
        lon = row.get("Longitude")
        polygon_wkt = row.get("Polygon")
        region_name = row.get("Region")

        # Process each scenario
        for scenario in config["scenarios"]:
            print(f"  Scenario: {scenario}")

            try:
                # Load dataset
                ds = load_processed_dataset(
                    processed_dir,
                    config["folder"],
                    config["variable_pattern"],
                    scenario,
                )

                # Extract for this location
                results = extract_for_location(
                    ds=ds,
                    location_name=base_name,
                    lat=lat,
                    lon=lon,
                    polygon_wkt=polygon_wkt,
                    region_name=region_name,
                    species=species,
                    scenario=scenario,
                    hazard_measure=config["hazard_measure"],
                )

                all_results.extend(results)
                print(f"    Extracted {len(results)} result(s)")

                ds.close()

            except Exception as e:
                print(f"  Error processing scenario {scenario}: {e}")

    # Format and export results
    print(f"\nFormatting {len(all_results)} extraction results...")
    output_df = format_extraction_results(all_results)

    # Generate output filename with date
    date_str = datetime.now().strftime("%Y%m%d")
    output_path = exports_dir / f"timber_locations_{date_str}.csv"

    # Export
    output_df.to_csv(output_path, index=False)
    print(f"\nExported {len(output_df)} rows to {output_path}")

    # Print summary
    print("\nSummary:")
    print(f"  Locations: {output_df['Location'].nunique()}")
    print(f"  Scenarios: {output_df['Scenario'].unique().tolist()}")
    print(f"  Decades: {sorted(output_df['Decade'].unique().tolist())}")
    print(f"  Total rows: {len(output_df)}")


if __name__ == "__main__":
    main()
