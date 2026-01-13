"""Tests for search result table display."""

import pytest
import json
import tempfile
from pathlib import Path

from isimip_pipeline.search.result_table import (
    ResultTable,
    format_results,
    export_selection,
    group_by_attributes,
)
from isimip_pipeline.search.isimip_query import DatasetInfo


@pytest.fixture
def sample_datasets():
    """Create sample DatasetInfo objects for testing."""
    return [
        DatasetInfo(
            id="1",
            name="gfdl_ssp126_led_annual.nc",
            url="https://files.isimip.org/1.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp126",
            variable="led",
            model="gfdl-esm4",
            timestep="annual",
        ),
        DatasetInfo(
            id="2",
            name="gfdl_ssp370_led_annual.nc",
            url="https://files.isimip.org/2.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp370",
            variable="led",
            model="gfdl-esm4",
            timestep="annual",
        ),
        DatasetInfo(
            id="3",
            name="ukesm_ssp126_led_monthly.nc",
            url="https://files.isimip.org/3.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp126",
            variable="led",
            model="ukesm1-0-ll",
            timestep="monthly",
        ),
    ]


class TestResultTable:
    """Test ResultTable class."""

    def test_result_table_initializes(self, sample_datasets):
        """ResultTable should initialize with datasets."""
        table = ResultTable(sample_datasets)
        assert table is not None
        assert len(table.datasets) == 3

    def test_result_table_generates_rich_table(self, sample_datasets):
        """ResultTable should generate Rich Table object."""
        table = ResultTable(sample_datasets)
        rich_table = table.to_rich_table()

        # Should have a table with rows
        assert rich_table is not None
        assert rich_table.row_count == 3

    def test_result_table_summary(self, sample_datasets):
        """ResultTable should provide summary statistics."""
        table = ResultTable(sample_datasets)
        summary = table.get_summary()

        assert summary["total_datasets"] == 3
        assert "ssp126" in summary["scenarios"]
        assert "ssp370" in summary["scenarios"]


class TestGroupByAttributes:
    """Test grouping datasets by attributes."""

    def test_group_by_scenario(self, sample_datasets):
        """Should group datasets by climate scenario."""
        grouped = group_by_attributes(sample_datasets, "climate_scenario")

        assert "ssp126" in grouped
        assert "ssp370" in grouped
        assert len(grouped["ssp126"]) == 2
        assert len(grouped["ssp370"]) == 1

    def test_group_by_model(self, sample_datasets):
        """Should group datasets by model."""
        grouped = group_by_attributes(sample_datasets, "model")

        assert "gfdl-esm4" in grouped
        assert "ukesm1-0-ll" in grouped
        assert len(grouped["gfdl-esm4"]) == 2

    def test_group_by_timestep(self, sample_datasets):
        """Should group datasets by timestep."""
        grouped = group_by_attributes(sample_datasets, "timestep")

        assert "annual" in grouped
        assert "monthly" in grouped


class TestExportSelection:
    """Test exporting dataset selection to JSON."""

    def test_export_to_json(self, sample_datasets):
        """Should export datasets to JSON file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            output_path = Path(f.name)

        try:
            export_selection(sample_datasets, output_path)

            # Verify file contents
            with open(output_path) as f:
                data = json.load(f)

            assert len(data["datasets"]) == 3
            assert data["datasets"][0]["id"] == "1"
            assert data["datasets"][0]["url"] == "https://files.isimip.org/1.nc"
        finally:
            output_path.unlink()

    def test_export_includes_metadata(self, sample_datasets):
        """Export should include metadata about the selection."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            output_path = Path(f.name)

        try:
            export_selection(sample_datasets, output_path)

            with open(output_path) as f:
                data = json.load(f)

            assert "total" in data
            assert data["total"] == 3
        finally:
            output_path.unlink()


class TestFormatResults:
    """Test the format_results convenience function."""

    def test_format_results_returns_string(self, sample_datasets):
        """format_results should return formatted string."""
        output = format_results(sample_datasets)

        assert isinstance(output, str)
        assert "ssp126" in output or "led" in output

    def test_format_results_empty_list(self):
        """format_results should handle empty list."""
        output = format_results([])

        assert "no" in output.lower() or "0" in output
