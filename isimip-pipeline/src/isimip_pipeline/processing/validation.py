"""Data validation module for ISIMIP climate datasets.

Performs comprehensive quality checks on climate data including:
- Fill value detection and characterization
- Missing data pattern analysis
- Spatial gap detection
- Temporal gap detection
- Statistical outlier detection
- Quality report generation

Designed to validate multi-model climate datasets before processing.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
from enum import Enum
from datetime import datetime

import numpy as np
import xarray as xr
from scipy import stats


class DataQualityLevel(Enum):
    """Data quality assessment levels."""
    EXCELLENT = "excellent"  # <1% issues
    GOOD = "good"  # 1-5% issues
    ACCEPTABLE = "acceptable"  # 5-10% issues
    POOR = "poor"  # 10-25% issues
    UNUSABLE = "unusable"  # >25% issues


@dataclass
class DataQualityIssue:
    """Represents a data quality issue."""
    issue_type: str
    severity: str  # "info", "warning", "error"
    description: str
    affected_fraction: float = 0.0  # Fraction of data affected (0-1)
    location: Optional[Dict[str, Any]] = None  # Optional location info


@dataclass
class FillValueInfo:
    """Information about detected fill values."""
    detected_fill_values: List[float]
    fill_value_count: int
    fill_value_percentage: float
    fill_value_locations: Dict[float, int]  # Value -> count mapping


@dataclass
class ValidationReport:
    """Complete validation report for a dataset."""
    variable: str
    timestamp: str
    n_cells: int
    n_time_steps: int
    total_points: int
    valid_points: int
    issues: List[DataQualityIssue] = field(default_factory=list)
    fill_value_info: Optional[FillValueInfo] = None
    gap_info: Dict[str, Any] = field(default_factory=dict)
    outlier_info: Dict[str, Any] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary."""
        return {
            "variable": self.variable,
            "timestamp": self.timestamp,
            "n_cells": self.n_cells,
            "n_time_steps": self.n_time_steps,
            "total_points": self.total_points,
            "valid_points": self.valid_points,
            "n_issues": len(self.issues),
            "summary": self.summary,
        }


def detect_fill_values(
    ds: xr.Dataset,
    variable: str,
    common_fill_values: Optional[List[float]] = None,
) -> Optional[FillValueInfo]:
    """Detect fill values in dataset variable.

    Args:
        ds: xarray Dataset.
        variable: Variable name to check.
        common_fill_values: List of common fill values to check.

    Returns:
        FillValueInfo with detected fill values, or None if none found.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found in dataset")

    if common_fill_values is None:
        common_fill_values = [
            1.00000002004088e+20,  # ISIMIP standard
            1.00000002e+20,
            -9999,
            -9999.0,
            -32768,
            32767,
            np.nan,
        ]

    data = ds[variable].values

    # Check variable attributes for fill value information
    var_attrs = ds[variable].attrs
    attr_fill_values = []

    if "_FillValue" in var_attrs:
        attr_fill_values.append(float(var_attrs["_FillValue"]))
    if "missing_value" in var_attrs:
        attr_fill_values.append(float(var_attrs["missing_value"]))

    # Combine lists
    all_fill_values = list(set(common_fill_values + attr_fill_values))

    detected_fills = {}
    total_fill_count = 0

    # Check for each potential fill value
    for fv in all_fill_values:
        if np.isnan(fv):
            # Check for NaN
            count = int(np.sum(np.isnan(data)))
        else:
            # Check for numeric value (with tolerance for floating point)
            if isinstance(fv, float):
                count = int(np.sum(np.abs(data - fv) < 0.1))
            else:
                count = int(np.sum(data == fv))

        if count > 0:
            detected_fills[fv] = count
            total_fill_count += count

    if not detected_fills:
        return None

    # Calculate percentage
    total_points = data.size
    fill_percentage = (total_fill_count / total_points) * 100

    return FillValueInfo(
        detected_fill_values=list(detected_fills.keys()),
        fill_value_count=total_fill_count,
        fill_value_percentage=fill_percentage,
        fill_value_locations=detected_fills,
    )


def find_missing_data_patterns(
    ds: xr.Dataset,
    variable: str,
) -> Dict[str, Any]:
    """Identify patterns in missing data.

    Args:
        ds: xarray Dataset.
        variable: Variable name to analyze.

    Returns:
        Dictionary with missing data patterns.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found")

    data = ds[variable].values

    # Count missing values (NaN, inf, etc.)
    is_missing = ~np.isfinite(data)
    n_missing = int(np.sum(is_missing))
    total_points = data.size
    missing_percentage = (n_missing / total_points) * 100

    # Analyze spatial patterns
    if len(data.shape) == 3:
        # (time, lat, lon)
        temporal_axis = 0
        spatial_shape = data.shape[1:]

        # Missing per grid cell (across time)
        missing_per_cell = np.sum(is_missing, axis=temporal_axis)
        completely_missing_cells = int(np.sum(missing_per_cell == data.shape[temporal_axis]))
        partially_missing_cells = int(np.sum((missing_per_cell > 0) & (missing_per_cell < data.shape[temporal_axis])))

    else:
        completely_missing_cells = 0
        partially_missing_cells = 0

    return {
        "missing_percentage": missing_percentage,
        "n_missing": n_missing,
        "total_points": total_points,
        "complete_cells": int(np.sum(missing_per_cell == 0)) if len(data.shape) == 3 else None,
        "partially_missing_cells": partially_missing_cells,
        "completely_missing_cells": completely_missing_cells,
    }


def detect_spatial_gaps(
    ds: xr.Dataset,
    variable: str,
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """Detect spatial gaps in data.

    Args:
        ds: xarray Dataset.
        variable: Variable name.
        threshold: Fraction threshold for considering cell as missing (0-1).

    Returns:
        Dictionary with spatial gap information.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found")

    data = ds[variable].values

    if len(data.shape) != 3:
        raise ValueError("Expected 3D data (time, lat, lon)")

    # Calculate missing fraction per grid cell
    is_missing = ~np.isfinite(data)
    missing_fraction_per_cell = np.mean(is_missing, axis=0)  # (lat, lon)

    # Count cells by missing fraction
    n_completely_missing = int(np.sum(missing_fraction_per_cell == 1.0))
    n_partially_missing = int(np.sum((missing_fraction_per_cell >= threshold) & (missing_fraction_per_cell < 1.0)))
    n_mostly_complete = int(np.sum(missing_fraction_per_cell < threshold))

    return {
        "n_missing_cells": n_completely_missing,
        "n_partially_missing_cells": n_partially_missing,
        "n_mostly_complete_cells": n_mostly_complete,
        "mean_missing_fraction": float(np.mean(missing_fraction_per_cell)),
        "max_missing_fraction": float(np.max(missing_fraction_per_cell)),
    }


def detect_temporal_gaps(
    ds: xr.Dataset,
    variable: str,
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """Detect temporal gaps in data.

    Args:
        ds: xarray Dataset.
        variable: Variable name.
        threshold: Fraction threshold for considering time step as gapped (0-1).

    Returns:
        Dictionary with temporal gap information.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found")

    data = ds[variable].values

    if len(data.shape) != 3:
        raise ValueError("Expected 3D data (time, lat, lon)")

    # Calculate missing fraction per time step
    is_missing = ~np.isfinite(data)
    missing_fraction_per_time = np.mean(is_missing, axis=(1, 2))  # (time,)

    # Identify gap periods
    gap_mask = missing_fraction_per_time >= threshold
    gap_starts = np.where(np.diff(np.concatenate([[False], gap_mask, [False]]).astype(int)) == 1)[0]
    gap_ends = np.where(np.diff(np.concatenate([[False], gap_mask, [False]]).astype(int)) == -1)[0]

    gap_periods = []
    for start, end in zip(gap_starts, gap_ends):
        gap_periods.append({
            "start": int(start),
            "end": int(end - 1),
            "duration": int(end - start),
        })

    return {
        "n_gap_periods": len(gap_periods),
        "gap_periods": gap_periods,
        "mean_missing_fraction_per_time": float(np.mean(missing_fraction_per_time)),
    }


def detect_outliers(
    ds: xr.Dataset,
    variable: str,
    method: str = "iqr",
    threshold: float = 3,
) -> Optional[Dict[str, Any]]:
    """Detect statistical outliers in data.

    Args:
        ds: xarray Dataset.
        variable: Variable name.
        method: "iqr" (interquartile range) or "zscore".
        threshold: Outlier threshold (IQR multiplier or z-score cutoff).

    Returns:
        Dictionary with outlier information, or None if no outliers.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found")

    data = ds[variable].values.flatten()

    # Remove NaN and inf
    valid_data = data[np.isfinite(data)]

    if len(valid_data) == 0:
        return None

    if method == "iqr":
        q1 = np.percentile(valid_data, 25)
        q3 = np.percentile(valid_data, 75)
        iqr = q3 - q1

        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr

        outlier_mask = (data < lower_bound) | (data > upper_bound)

    elif method == "zscore":
        mean = np.mean(valid_data)
        std = np.std(valid_data)

        if std == 0:
            return {"n_outliers": 0}

        z_scores = np.abs((data - mean) / std)
        outlier_mask = z_scores > threshold

    else:
        raise ValueError(f"Unknown method: {method}")

    n_outliers = int(np.sum(outlier_mask))

    if n_outliers == 0:
        return {"n_outliers": 0}

    outlier_values = data[outlier_mask]
    outlier_percentage = (n_outliers / len(data)) * 100
    outlier_locations = np.where(outlier_mask)[0]  # Flat indices of outliers

    return {
        "n_outliers": n_outliers,
        "outlier_percentage": outlier_percentage,
        "outlier_values": outlier_values,
        "outlier_locations": outlier_locations,
        "min_outlier": float(np.min(outlier_values)),
        "max_outlier": float(np.max(outlier_values)),
    }


def generate_validation_report(
    ds: xr.Dataset,
    variable: str,
    include_fill_values: bool = True,
    include_gaps: bool = True,
    include_outliers: bool = True,
) -> ValidationReport:
    """Generate comprehensive validation report.

    Args:
        ds: xarray Dataset.
        variable: Variable name to validate.
        include_fill_values: Whether to check for fill values.
        include_gaps: Whether to check for spatial/temporal gaps.
        include_outliers: Whether to check for outliers.

    Returns:
        ValidationReport with all checks.
    """
    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found")

    data = ds[variable].values
    is_valid = np.isfinite(data)
    n_valid = int(np.sum(is_valid))
    total_points = data.size

    # Get grid dimensions
    if len(data.shape) == 3:
        n_time, n_lat, n_lon = data.shape
        n_cells = n_lat * n_lon
    else:
        n_time = len(data) if data.size > 0 else 0
        n_cells = 1
        n_lat, n_lon = 1, 1

    issues = []
    fill_value_info = None
    gap_info = {}
    outlier_info = {}
    summary = {}

    # Check fill values
    if include_fill_values:
        fill_value_info = detect_fill_values(ds, variable)
        if fill_value_info:
            issues.append(DataQualityIssue(
                issue_type="fill_values",
                severity="warning" if fill_value_info.fill_value_percentage < 10 else "error",
                description=f"Detected {fill_value_info.fill_value_count} fill values "
                           f"({fill_value_info.fill_value_percentage:.1f}%)",
                affected_fraction=fill_value_info.fill_value_percentage / 100,
            ))

    # Check gaps
    if include_gaps:
        spatial_gaps = detect_spatial_gaps(ds, variable)
        temporal_gaps = detect_temporal_gaps(ds, variable)
        gap_info = {"spatial": spatial_gaps, "temporal": temporal_gaps}

        if spatial_gaps["n_missing_cells"] > 0:
            issues.append(DataQualityIssue(
                issue_type="spatial_gaps",
                severity="warning",
                description=f"Found {spatial_gaps['n_missing_cells']} grid cells with complete data loss",
                affected_fraction=spatial_gaps["n_missing_cells"] / n_cells if n_cells > 0 else 0,
            ))

        if temporal_gaps["n_gap_periods"] > 0:
            issues.append(DataQualityIssue(
                issue_type="temporal_gaps",
                severity="warning",
                description=f"Found {temporal_gaps['n_gap_periods']} temporal gap periods",
                affected_fraction=temporal_gaps["mean_missing_fraction_per_time"],
            ))

    # Check outliers
    if include_outliers:
        outlier_info = detect_outliers(ds, variable) or {}
        if outlier_info.get("n_outliers", 0) > 0:
            issues.append(DataQualityIssue(
                issue_type="outliers",
                severity="info",
                description=f"Found {outlier_info['n_outliers']} statistical outliers "
                           f"({outlier_info.get('outlier_percentage', 0):.1f}%)",
                affected_fraction=outlier_info.get("outlier_percentage", 0) / 100,
            ))

    # Generate summary
    data_loss_percentage = ((total_points - n_valid) / total_points) * 100

    if data_loss_percentage > 25:
        data_quality = DataQualityLevel.UNUSABLE
    elif data_loss_percentage > 10:
        data_quality = DataQualityLevel.POOR
    elif data_loss_percentage > 5:
        data_quality = DataQualityLevel.ACCEPTABLE
    elif data_loss_percentage > 1:
        data_quality = DataQualityLevel.GOOD
    else:
        data_quality = DataQualityLevel.EXCELLENT

    summary = {
        "data_quality": data_quality.value,
        "data_loss_percentage": data_loss_percentage,
        "valid_percentage": (n_valid / total_points) * 100,
        "n_issues": len(issues),
        "check_timestamp": datetime.now().isoformat(),
    }

    return ValidationReport(
        variable=variable,
        timestamp=datetime.now().isoformat(),
        n_cells=n_cells,
        n_time_steps=n_time,
        total_points=total_points,
        valid_points=n_valid,
        issues=issues,
        fill_value_info=fill_value_info,
        gap_info=gap_info,
        outlier_info=outlier_info,
        summary=summary,
    )


def validate_dataset(
    ds: xr.Dataset,
    variable: str,
    raise_on_error: bool = False,
) -> ValidationReport:
    """Validate dataset with default comprehensive checks.

    Args:
        ds: xarray Dataset.
        variable: Variable name to validate.
        raise_on_error: Whether to raise exception on major issues.

    Returns:
        ValidationReport.

    Raises:
        ValueError: If data quality is poor and raise_on_error=True.
    """
    report = generate_validation_report(ds, variable)

    if raise_on_error and report.summary["data_quality"] in ["poor", "unusable"]:
        raise ValueError(
            f"Data quality too poor for variable '{variable}': "
            f"{report.summary['data_loss_percentage']:.1f}% data loss"
        )

    return report
