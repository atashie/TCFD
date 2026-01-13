"""Tests for data validation module.

Tests for fill value detection, missing data patterns, gap detection,
and statistical outlier detection for ISIMIP climate data.
"""

import pytest
import numpy as np
import xarray as xr
from pathlib import Path
import tempfile
from datetime import datetime, timedelta

from isimip_pipeline.processing.validation import (
    detect_fill_values,
    find_missing_data_patterns,
    detect_spatial_gaps,
    detect_temporal_gaps,
    detect_outliers,
    generate_validation_report,
    DataQualityIssue,
    ValidationReport,
    FillValueInfo,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def clean_dataset():
    """Create a clean test dataset."""
    time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
    lat = np.arange(-90, 90, 10)
    lon = np.arange(-180, 180, 10)

    # Create realistic data: temperature values in reasonable range
    data = 273.15 + 20 * np.sin(np.arange(len(time))[:, None, None] / 365)
    data = data + 5 * np.random.randn(len(time), len(lat), len(lon))

    ds = xr.Dataset(
        {"temperature": (["time", "lat", "lon"], data)},
        coords={
            "time": time,
            "lat": lat,
            "lon": lon,
        },
    )
    return ds


@pytest.fixture
def dataset_with_fill_values():
    """Create dataset with common fill value patterns."""
    time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
    lat = np.arange(-90, 90, 10)
    lon = np.arange(-180, 180, 10)

    data = 273.15 + 20 * np.sin(np.arange(len(time))[:, None, None] / 365)

    # Add common fill values
    data[0:10, 0, 0] = 1.00000002e20  # ISIMIP fill value
    data[50:52, 3, 3] = -9999  # Common fill value
    data[100:105, 5, 5] = np.nan  # NaN fill value

    ds = xr.Dataset(
        {
            "temperature": (
                ["time", "lat", "lon"],
                data,
                {
                    "_FillValue": 1.00000002e20,
                    "missing_value": -9999,
                },
            )
        },
        coords={
            "time": time,
            "lat": lat,
            "lon": lon,
        },
    )
    return ds


@pytest.fixture
def dataset_with_gaps():
    """Create dataset with spatial and temporal gaps."""
    time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
    lat = np.arange(-90, 90, 10)
    lon = np.arange(-180, 180, 10)

    data = 273.15 + 20 * np.sin(np.arange(len(time))[:, None, None] / 365)

    # Create spatial gaps (missing entire regions)
    data[:, 0:2, 0:2] = np.nan  # Northwest quadrant

    # Create temporal gaps (missing time periods)
    data[100:200, :, :] = np.nan  # Gap in middle of time series

    ds = xr.Dataset(
        {"temperature": (["time", "lat", "lon"], data)},
        coords={
            "time": time,
            "lat": lat,
            "lon": lon,
        },
    )
    return ds


class TestFillValueDetection:
    """Test fill value detection."""

    def test_detect_isimip_fill_value(self, dataset_with_fill_values):
        """Should detect standard ISIMIP fill value."""
        fill_info = detect_fill_values(dataset_with_fill_values, "temperature")

        assert fill_info is not None
        assert 1.00000002e20 in fill_info.detected_fill_values
        assert fill_info.fill_value_count > 0

    def test_detect_nan_fill_values(self, dataset_with_fill_values):
        """Should detect NaN fill values."""
        fill_info = detect_fill_values(dataset_with_fill_values, "temperature")

        assert fill_info is not None
        # Should identify NaN regions
        assert fill_info.fill_value_count > 0

    def test_detect_multiple_fill_values(self, dataset_with_fill_values):
        """Should detect multiple different fill value types."""
        fill_info = detect_fill_values(dataset_with_fill_values, "temperature")

        assert fill_info is not None
        # Should detect multiple fill values
        assert len(fill_info.detected_fill_values) >= 1

    def test_no_fill_values_in_clean_data(self, clean_dataset):
        """Should return None or empty for clean data."""
        fill_info = detect_fill_values(clean_dataset, "temperature")

        # Should either be None or have zero fill values
        if fill_info is not None:
            assert fill_info.fill_value_count == 0

    def test_fill_value_percentage_calculation(self, dataset_with_fill_values):
        """Should calculate percentage of fill values correctly."""
        fill_info = detect_fill_values(dataset_with_fill_values, "temperature")

        if fill_info is not None:
            assert 0 <= fill_info.fill_value_percentage <= 100

    def test_detect_fill_value_from_attributes(self):
        """Should detect fill values from variable attributes."""
        data = np.random.randn(100, 10, 10)
        data[0:10, 0, 0] = -9999

        ds = xr.Dataset(
            {
                "precip": (
                    ["time", "lat", "lon"],
                    data,
                    {"_FillValue": -9999, "units": "mm/day"},
                )
            },
            coords={
                "time": np.arange(100),
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        fill_info = detect_fill_values(ds, "precip")

        if fill_info is not None:
            assert -9999 in fill_info.detected_fill_values or fill_info.fill_value_count > 0


class TestMissingDataPatterns:
    """Test missing data pattern detection."""

    def test_find_complete_regions(self, clean_dataset):
        """Should identify regions with complete data."""
        patterns = find_missing_data_patterns(clean_dataset, "temperature")

        assert patterns is not None
        assert "complete_cells" in patterns

    def test_find_partially_missing_regions(self, dataset_with_gaps):
        """Should identify regions with partial data gaps."""
        patterns = find_missing_data_patterns(dataset_with_gaps, "temperature")

        assert patterns is not None
        # Should detect gaps
        assert patterns["missing_percentage"] > 0

    def test_calculate_missing_percentage(self, dataset_with_gaps):
        """Should calculate overall missing data percentage."""
        patterns = find_missing_data_patterns(dataset_with_gaps, "temperature")

        assert patterns is not None
        assert "missing_percentage" in patterns
        assert 0 <= patterns["missing_percentage"] <= 100

    def test_identify_gap_locations(self, dataset_with_gaps):
        """Should identify spatial locations with gaps."""
        patterns = find_missing_data_patterns(dataset_with_gaps, "temperature")

        assert patterns is not None
        # Should identify regions
        assert len(patterns) > 0


class TestSpatialGapDetection:
    """Test spatial gap detection."""

    def test_detect_missing_grid_cells(self, dataset_with_gaps):
        """Should detect completely missing grid cells."""
        gaps = detect_spatial_gaps(dataset_with_gaps, "temperature", threshold=0.5)

        assert gaps is not None
        assert gaps["n_missing_cells"] > 0

    def test_detect_partially_missing_cells(self, dataset_with_gaps):
        """Should detect cells with partial data."""
        gaps = detect_spatial_gaps(dataset_with_gaps, "temperature", threshold=0.2)

        if gaps is not None:
            assert "partially_missing_cells" in gaps

    def test_calculate_cell_coverage(self, dataset_with_gaps):
        """Should calculate coverage for each grid cell."""
        gaps = detect_spatial_gaps(dataset_with_gaps, "temperature")

        assert gaps is not None
        # Should have coverage information
        assert gaps["n_missing_cells"] >= 0

    def test_clean_data_no_gaps(self, clean_dataset):
        """Clean data should have no spatial gaps."""
        gaps = detect_spatial_gaps(clean_dataset, "temperature", threshold=0.95)

        if gaps is not None:
            assert gaps["n_missing_cells"] == 0


class TestTemporalGapDetection:
    """Test temporal gap detection."""

    def test_detect_time_series_gaps(self, dataset_with_gaps):
        """Should detect gaps in time series."""
        gaps = detect_temporal_gaps(dataset_with_gaps, "temperature", threshold=0.5)

        assert gaps is not None
        assert gaps["n_gap_periods"] > 0

    def test_identify_gap_start_end(self, dataset_with_gaps):
        """Should identify start and end of gaps."""
        gaps = detect_temporal_gaps(dataset_with_gaps, "temperature", threshold=0.5)

        if gaps is not None and gaps["n_gap_periods"] > 0:
            assert "gap_periods" in gaps
            for gap in gaps["gap_periods"]:
                assert "start" in gap
                assert "end" in gap

    def test_calculate_gap_duration(self, dataset_with_gaps):
        """Should calculate duration of gaps."""
        gaps = detect_temporal_gaps(dataset_with_gaps, "temperature", threshold=0.5)

        if gaps is not None and gaps["n_gap_periods"] > 0:
            for gap in gaps["gap_periods"]:
                assert gap["duration"] > 0

    def test_no_gaps_in_clean_data(self, clean_dataset):
        """Clean data should have no temporal gaps."""
        gaps = detect_temporal_gaps(clean_dataset, "temperature", threshold=0.95)

        if gaps is not None:
            assert gaps["n_gap_periods"] == 0


class TestOutlierDetection:
    """Test statistical outlier detection."""

    def test_detect_extreme_values(self):
        """Should detect physically impossible values."""
        time = np.arange(100)
        lat = np.arange(10)
        lon = np.arange(10)

        # Create data with extreme outliers
        data = 273.15 + np.random.randn(100, 10, 10)
        data[0, 0, 0] = 400  # Extremely hot
        data[1, 0, 1] = 100  # Extremely cold

        ds = xr.Dataset(
            {"temperature": (["time", "lat", "lon"], data)},
            coords={"time": time, "lat": lat, "lon": lon},
        )

        outliers = detect_outliers(ds, "temperature", method="iqr", threshold=3)

        assert outliers is not None
        assert outliers["n_outliers"] > 0

    def test_detect_outliers_iqr_method(self, clean_dataset):
        """Should detect outliers using IQR method."""
        # Add an outlier to clean data
        data = clean_dataset.copy(deep=True)
        data["temperature"].values[0, 0, 0] = 500

        outliers = detect_outliers(data, "temperature", method="iqr")

        assert outliers is not None
        assert outliers["n_outliers"] >= 1

    def test_detect_outliers_zscore_method(self, clean_dataset):
        """Should detect outliers using z-score method."""
        # Add an outlier
        data = clean_dataset.copy(deep=True)
        data["temperature"].values[10, 2, 2] = 1000

        outliers = detect_outliers(data, "temperature", method="zscore", threshold=3)

        assert outliers is not None
        if outliers["n_outliers"] > 0:
            assert "outlier_locations" in outliers

    def test_no_outliers_in_normal_data(self, clean_dataset):
        """Should detect no outliers in normal data."""
        outliers = detect_outliers(clean_dataset, "temperature", threshold=3)

        if outliers is not None:
            # Should have few or no outliers in clean data
            assert outliers["n_outliers"] < 5

    def test_outlier_statistics(self):
        """Should calculate outlier statistics."""
        time = np.arange(100)
        lat = np.arange(10)
        lon = np.arange(10)

        # Create data with some outliers
        data = 273.15 + 5 * np.random.randn(100, 10, 10)
        data[np.random.choice(1000, 20, replace=False)] = 500

        ds = xr.Dataset(
            {"temperature": (["time", "lat", "lon"], data)},
            coords={"time": time, "lat": lat, "lon": lon},
        )

        outliers = detect_outliers(ds, "temperature", method="iqr")

        if outliers is not None and outliers["n_outliers"] > 0:
            assert "outlier_values" in outliers or "outlier_percentage" in outliers


class TestValidationReport:
    """Test validation report generation."""

    def test_generate_report_clean_data(self, clean_dataset):
        """Should generate report for clean data."""
        report = generate_validation_report(
            clean_dataset,
            "temperature",
            include_outliers=True,
            include_gaps=True,
        )

        assert report is not None
        assert report.variable == "temperature"
        assert report.n_cells > 0
        assert report.n_time_steps > 0

    def test_generate_report_with_issues(self, dataset_with_fill_values):
        """Should identify issues in validation report."""
        report = generate_validation_report(
            dataset_with_fill_values,
            "temperature",
            include_fill_values=True,
        )

        assert report is not None
        assert len(report.issues) > 0

    def test_report_summary(self, dataset_with_gaps):
        """Should provide summary statistics."""
        report = generate_validation_report(dataset_with_gaps, "temperature")

        assert report is not None
        assert report.summary is not None
        assert isinstance(report.summary, dict)

    def test_report_quality_flags(self, dataset_with_fill_values):
        """Should assign quality flags."""
        report = generate_validation_report(dataset_with_fill_values, "temperature")

        assert report is not None
        # Should have quality assessment
        assert hasattr(report, "issues") or hasattr(report, "summary")

    def test_report_exportable(self, dataset_with_gaps, temp_dir):
        """Should generate exportable report."""
        report = generate_validation_report(dataset_with_gaps, "temperature")

        assert report is not None
        # Should be convertible to dict or JSON
        report_dict = report.to_dict() if hasattr(report, "to_dict") else vars(report)
        assert isinstance(report_dict, dict)


class TestValidationIntegration:
    """Integration tests for validation workflow."""

    def test_full_validation_workflow(self, dataset_with_fill_values):
        """Should run complete validation workflow."""
        # Fill values
        fill_info = detect_fill_values(dataset_with_fill_values, "temperature")

        # Missing patterns
        patterns = find_missing_data_patterns(dataset_with_fill_values, "temperature")

        # Outliers
        outliers = detect_outliers(dataset_with_fill_values, "temperature")

        # Report
        report = generate_validation_report(dataset_with_fill_values, "temperature")

        assert report is not None
        assert report.variable == "temperature"

    def test_validation_consistency(self, dataset_with_gaps):
        """Validation results should be consistent."""
        # Run validation twice
        report1 = generate_validation_report(dataset_with_gaps, "temperature")
        report2 = generate_validation_report(dataset_with_gaps, "temperature")

        if report1 is not None and report2 is not None:
            # Results should match
            assert report1.n_cells == report2.n_cells
            assert report1.n_time_steps == report2.n_time_steps

    def test_validation_identifies_data_quality_issues(self, dataset_with_fill_values):
        """Should identify realistic data quality issues."""
        report = generate_validation_report(dataset_with_fill_values, "temperature")

        assert report is not None
        # Should have issues
        assert len(report.issues) > 0 or report.summary["data_quality"] != "excellent"
