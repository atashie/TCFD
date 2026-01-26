"""Export formatter for TCFD CSV output.

Converts extracted climate data to the standardized Export-Key.csv format
with all 28 columns including hazard scores, trends, and metadata.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd


# Relative Hazard Score thresholds (percentile -> label)
RELATIVE_HAZARD_THRESHOLDS = [
    (1, 20, "Very Low", 1),
    (21, 40, "Low", 2),
    (41, 60, "Medium", 3),
    (61, 80, "High", 4),
    (81, 100, "Very High", 5),
]

# Standard column order for export (matches Export-Key.csv)
EXPORT_COLUMNS = [
    "Location",
    "Region",
    "Subregion",
    "Lat",
    "Lon",
    "Hazard",
    "Hazard_Measure",
    "Decade",
    "Scenario",
    "Raw_Hazard_Value",
    "Percentile_Score",
    "Relative_Hazard_Score",
    "Decadal_Trend_Strength",
    "Decadal_Trend_Significance",
    "Long_Term_Trend_Strength",
    "Long_Term_Trend_Significance",
    "Relative_Hazard_Score_Number",
    "Trend_Aggregated_For_Looker",
    "Advanced_Data_Measures",
    "Advanced_Data_Measures_Units",
    "Raw_Hazard_Value_25th",
    "Raw_Hazard_Value_75th",
    "Asset_Type",
    "Sub-Asset_Unit",
    "Country",
    "State",
    "City",
    "Rating",
]

# NetCDF variable to Export column mapping
VARIABLE_MAPPING = {
    "median": "Raw_Hazard_Value",
    "percentile": "Percentile_Score",
    "trend": "Decadal_Trend_Strength",
    "trend_pvalue": "Decadal_Trend_Significance",
    "lower_ci": "Raw_Hazard_Value_25th",
    "upper_ci": "Raw_Hazard_Value_75th",
}


@dataclass
class LocationMetadata:
    """Metadata for a single extraction location."""

    location: str
    lat: Optional[float] = None
    lon: Optional[float] = None
    region: str = "NA"
    subregion: str = "NA"
    country: str = "NA"
    state: str = "NA"
    city: str = "NA"
    asset_type: str = "NA"
    sub_asset_unit: str = "NA"

    # For region/polygon extractions, centroid is calculated
    centroid_lat: Optional[float] = None
    centroid_lon: Optional[float] = None

    def get_lat(self) -> Optional[float]:
        """Get latitude (user-provided or centroid)."""
        return self.lat if self.lat is not None else self.centroid_lat

    def get_lon(self) -> Optional[float]:
        """Get longitude (user-provided or centroid)."""
        return self.lon if self.lon is not None else self.centroid_lon


@dataclass
class HazardMapping:
    """Mapping from processed data folder to hazard classification."""

    folder_pattern: str
    hazard: str
    hazard_measure: str
    units: str = ""
    advanced_data_measures: str = "NA"
    advanced_data_measures_units: str = "NA"
    percentile_direction: str = "higher_is_worse"  # or "higher_is_better"


@dataclass
class ExtractionResult:
    """Result of a single extraction operation."""

    location_metadata: LocationMetadata
    hazard_mapping: HazardMapping
    scenario: str
    data: Dict[str, Dict[int, float]] = field(default_factory=dict)
    # data format: {variable_name: {decade: value}}


def get_relative_hazard_score(percentile: float) -> tuple:
    """Convert percentile to relative hazard score.

    Thresholds (continuous, no gaps):
      0-20  -> Very Low (1)
      20-40 -> Low (2)
      40-60 -> Medium (3)
      60-80 -> High (4)
      80-100 -> Very High (5)

    Args:
        percentile: Percentile score (0-100)

    Returns:
        Tuple of (label, number) e.g., ("Medium", 3)
    """
    if pd.isna(percentile):
        return ("NA", np.nan)

    percentile = float(percentile)

    if percentile <= 20:
        return ("Very Low", 1)
    elif percentile <= 40:
        return ("Low", 2)
    elif percentile <= 60:
        return ("Medium", 3)
    elif percentile <= 80:
        return ("High", 4)
    else:
        return ("Very High", 5)


def calculate_trend_aggregated(
    trend_strength: float,
    trend_significance: float,
    long_term_strength: float,
    long_term_significance: float,
    significance_threshold: float = 0.05,
    decadal_threshold: float = 0.5,
) -> float:
    """Calculate aggregated trend for Looker visualization.

    Returns trend strength only if both decadal and long-term trends
    are statistically significant and in the same direction.

    Args:
        trend_strength: Decadal trend strength
        trend_significance: Decadal trend p-value
        long_term_strength: Long-term trend strength
        long_term_significance: Long-term trend p-value
        significance_threshold: P-value threshold for long-term significance
        decadal_threshold: P-value threshold for decadal significance

    Returns:
        Trend strength if significant, else 0
    """
    if any(pd.isna([trend_strength, trend_significance, long_term_strength, long_term_significance])):
        return 0.0

    # Check significance thresholds
    if trend_significance >= decadal_threshold:
        return 0.0
    if long_term_significance >= significance_threshold:
        return 0.0

    # Check direction consistency
    if (trend_strength > 0) != (long_term_strength > 0):
        return 0.0

    return float(trend_strength)


def format_single_extraction(
    result: ExtractionResult,
) -> pd.DataFrame:
    """Format a single extraction result to DataFrame rows.

    Each decade becomes a separate row.

    Args:
        result: ExtractionResult with extracted data

    Returns:
        DataFrame with one row per decade
    """
    rows = []
    meta = result.location_metadata
    hazard = result.hazard_mapping

    # Get decades from any variable
    decades = []
    for var_data in result.data.values():
        decades = sorted(var_data.keys())
        break

    if not decades:
        return pd.DataFrame(columns=EXPORT_COLUMNS)

    # Get long-term values from last decade (2090)
    last_decade = max(decades)
    long_term_trend = result.data.get("trend", {}).get(last_decade, np.nan)
    long_term_pvalue = result.data.get("trend_pvalue", {}).get(last_decade, np.nan)

    for decade in decades:
        # Extract values for this decade
        raw_value = result.data.get("median", {}).get(decade, np.nan)
        percentile = result.data.get("percentile", {}).get(decade, np.nan)
        trend = result.data.get("trend", {}).get(decade, np.nan)
        trend_pvalue = result.data.get("trend_pvalue", {}).get(decade, np.nan)
        lower_ci = result.data.get("lower_ci", {}).get(decade, np.nan)
        upper_ci = result.data.get("upper_ci", {}).get(decade, np.nan)

        # Calculate derived values
        rel_hazard_label, rel_hazard_num = get_relative_hazard_score(percentile)
        trend_agg = calculate_trend_aggregated(
            trend, trend_pvalue, long_term_trend, long_term_pvalue
        )

        row = {
            "Location": meta.location,
            "Region": meta.region,
            "Subregion": meta.subregion,
            "Lat": meta.get_lat(),
            "Lon": meta.get_lon(),
            "Hazard": hazard.hazard,
            "Hazard_Measure": hazard.hazard_measure,
            "Decade": decade,
            "Scenario": result.scenario,
            "Raw_Hazard_Value": raw_value,
            "Percentile_Score": percentile,
            "Relative_Hazard_Score": rel_hazard_label,
            "Decadal_Trend_Strength": trend,
            "Decadal_Trend_Significance": trend_pvalue,
            "Long_Term_Trend_Strength": long_term_trend,
            "Long_Term_Trend_Significance": long_term_pvalue,
            "Relative_Hazard_Score_Number": rel_hazard_num,
            "Trend_Aggregated_For_Looker": trend_agg,
            "Advanced_Data_Measures": hazard.advanced_data_measures,
            "Advanced_Data_Measures_Units": hazard.advanced_data_measures_units,
            "Raw_Hazard_Value_25th": lower_ci,
            "Raw_Hazard_Value_75th": upper_ci,
            "Asset_Type": meta.asset_type,
            "Sub-Asset_Unit": meta.sub_asset_unit,
            "Country": meta.country,
            "State": meta.state,
            "City": meta.city,
            "Rating": "NA",  # To be calculated across all locations
        }
        rows.append(row)

    return pd.DataFrame(rows, columns=EXPORT_COLUMNS)


def format_extraction_results(
    results: List[ExtractionResult],
) -> pd.DataFrame:
    """Format multiple extraction results to a single DataFrame.

    Args:
        results: List of ExtractionResult objects

    Returns:
        Combined DataFrame with all results
    """
    if not results:
        return pd.DataFrame(columns=EXPORT_COLUMNS)

    dfs = [format_single_extraction(r) for r in results]
    combined = pd.concat(dfs, ignore_index=True)

    # Ensure column order matches Export-Key.csv
    return combined[EXPORT_COLUMNS]


def export_to_csv(
    results: List[ExtractionResult],
    output_path: Union[str, Path],
    include_index: bool = False,
) -> Path:
    """Export extraction results to CSV file.

    Args:
        results: List of ExtractionResult objects
        output_path: Path for output CSV file
        include_index: Whether to include row index in CSV

    Returns:
        Path to the created CSV file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = format_extraction_results(results)
    df.to_csv(output_path, index=include_index)

    print(f"Exported {len(df)} rows to {output_path}")
    return output_path


def load_locations_csv(
    csv_path: Union[str, Path],
) -> List[LocationMetadata]:
    """Load location definitions from CSV file.

    Expected columns (all optional except Location):
    - Location: Required, user-defined name
    - Type: "point", "polygon", or "region"
    - Lat, Lon: For point extractions
    - Polygon_Coords: For polygon extractions (JSON list of [lon, lat] pairs)
    - Region_Name, Region_Type: For region extractions
    - Region, Subregion, Country, State, City: Optional metadata
    - Asset_Type, Sub-Asset_Unit: Optional asset classification

    Args:
        csv_path: Path to locations CSV

    Returns:
        List of LocationMetadata objects
    """
    df = pd.read_csv(csv_path)

    locations = []
    for _, row in df.iterrows():
        meta = LocationMetadata(
            location=row.get("Location", "Unknown"),
            lat=row.get("Lat"),
            lon=row.get("Lon"),
            region=row.get("Region", "NA"),
            subregion=row.get("Subregion", "NA"),
            country=row.get("Country", "NA"),
            state=row.get("State", "NA"),
            city=row.get("City", "NA"),
            asset_type=row.get("Asset_Type", "NA"),
            sub_asset_unit=row.get("Sub-Asset_Unit", "NA"),
        )
        locations.append(meta)

    return locations


def load_hazard_mapping_csv(
    csv_path: Union[str, Path],
) -> List[HazardMapping]:
    """Load hazard mapping definitions from CSV file.

    Expected columns:
    - folder_pattern: Substring to match in processed data folder
    - Hazard: Hazard category name
    - Hazard_Measure: Specific measure name
    - Units: (Optional) units for hazard measure
    - Advanced_Data_Measures: (Optional) threshold description
    - Advanced_Data_Measures_Units: (Optional) threshold units

    Args:
        csv_path: Path to hazard mapping CSV

    Returns:
        List of HazardMapping objects
    """
    df = pd.read_csv(csv_path)

    mappings = []
    for _, row in df.iterrows():
        mapping = HazardMapping(
            folder_pattern=row["folder_pattern"],
            hazard=row["Hazard"],
            hazard_measure=row["Hazard_Measure"],
            units=row.get("Units", ""),
            advanced_data_measures=row.get("Advanced_Data_Measures", "NA"),
            advanced_data_measures_units=row.get("Advanced_Data_Measures_Units", "NA"),
        )
        mappings.append(mapping)

    return mappings


def create_empty_export_df() -> pd.DataFrame:
    """Create empty DataFrame with all export columns.

    Useful for initializing before appending results.

    Returns:
        Empty DataFrame with Export-Key.csv columns
    """
    return pd.DataFrame(columns=EXPORT_COLUMNS)


def validate_export_df(df: pd.DataFrame) -> List[str]:
    """Validate that DataFrame has all required columns.

    Args:
        df: DataFrame to validate

    Returns:
        List of missing column names (empty if valid)
    """
    missing = [col for col in EXPORT_COLUMNS if col not in df.columns]
    return missing
