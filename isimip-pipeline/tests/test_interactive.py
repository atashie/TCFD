"""Tests for interactive workflow module."""

import json
import tempfile
from pathlib import Path
from datetime import datetime
from unittest.mock import Mock, patch, MagicMock

import pytest

from isimip_pipeline.processing_log import DatasetEntry, ProcessingLog, save_processing_log
from isimip_pipeline.search.isimip_query import DatasetInfo
from isimip_pipeline.interactive import (
    save_selection_metadata,
    load_selection_metadata,
    parse_descriptive_name_from_folder,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_log():
    """Create a sample processing log."""
    log = ProcessingLog()

    entry = DatasetEntry(
        descriptive_name="drought-severity",
        variable="led",
        timestep="monthly",
        created_date=datetime(2026, 1, 10, 14, 30, 0),
        output_path="./outputs/drought-severity_led-monthly",
        file_count=15,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126", "ssp370", "ssp585"],
        gcm_models=["gfdl-esm4", "ukesm1-0-ll"],
        lsm_models=["clm5.0", "lpj-guess"],
        simulation_round="ISIMIP3b",
        query="drought exposure metrics",
    )
    log.add_entry(entry)

    return log


@pytest.fixture
def sample_datasets():
    """Create sample ISIMIP datasets."""
    return [
        DatasetInfo(
            id="1",
            name="led_clm45_gfdl-esm2m_historical.nc",
            url="https://example.com/led_clm45_gfdl-esm2m_historical.nc",
            simulation_round="ISIMIP2b",
            climate_scenario="historical",
            variable="led",
            model="clm45",
            timestep="monthly",
            size=1000000,
        ),
        DatasetInfo(
            id="2",
            name="led_clm45_gfdl-esm2m_rcp26.nc",
            url="https://example.com/led_clm45_gfdl-esm2m_rcp26.nc",
            simulation_round="ISIMIP2b",
            climate_scenario="rcp26",
            variable="led",
            model="clm45",
            timestep="monthly",
            size=1000000,
        ),
        DatasetInfo(
            id="3",
            name="led_clm45_gfdl-esm2m_annual.nc",
            url="https://example.com/led_clm45_gfdl-esm2m_annual.nc",
            simulation_round="ISIMIP2b",
            climate_scenario="rcp26",
            variable="led",
            model="clm45",
            timestep="annual",
            size=800000,
        ),
    ]


class TestSaveSelectionMetadata:
    """Tests for save_selection_metadata function."""

    def test_save_selection_creates_json(self, temp_dir, sample_datasets):
        """Save selection creates JSON file with metadata."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="drought exposure",
            descriptive_name="drought-test",
        )

        selection_file = output_dir / "selection.json"
        assert selection_file.exists()

    def test_save_selection_has_required_fields(self, temp_dir, sample_datasets):
        """Saved selection has all required fields."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="drought",
            descriptive_name="drought-test",
        )

        selection_file = output_dir / "selection.json"
        with open(selection_file) as f:
            data = json.load(f)

        assert "exported_at" in data
        assert "query" in data
        assert "descriptive_name" in data
        assert "total" in data
        assert "datasets" in data

    def test_save_selection_preserves_dataset_info(self, temp_dir, sample_datasets):
        """Saved selection preserves all dataset information."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="test",
            descriptive_name="test",
        )

        selection_file = output_dir / "selection.json"
        with open(selection_file) as f:
            data = json.load(f)

        assert len(data["datasets"]) == len(sample_datasets)
        assert data["datasets"][0]["variable"] == "led"
        assert data["datasets"][0]["timestep"] == "monthly"

    def test_save_selection_includes_query_metadata(self, temp_dir, sample_datasets):
        """Saved selection includes original query metadata."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        query = "drought exposure metrics"
        desc_name = "drought-severity"

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query=query,
            descriptive_name=desc_name,
        )

        selection_file = output_dir / "selection.json"
        with open(selection_file) as f:
            data = json.load(f)

        assert data["query"] == query
        assert data["descriptive_name"] == desc_name


class TestLoadSelectionMetadata:
    """Tests for load_selection_metadata function."""

    def test_load_selection_from_file(self, temp_dir, sample_datasets):
        """Load selection from saved file."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="drought",
            descriptive_name="drought-test",
        )

        loaded = load_selection_metadata(output_dir)

        assert loaded["query"] == "drought"
        assert loaded["descriptive_name"] == "drought-test"
        assert len(loaded["datasets"]) == len(sample_datasets)

    def test_load_missing_selection_raises_error(self, temp_dir):
        """Loading from missing file raises error."""
        output_dir = temp_dir / "nonexistent"
        output_dir.mkdir()

        with pytest.raises(FileNotFoundError):
            load_selection_metadata(output_dir)

    def test_load_selection_preserves_types(self, temp_dir, sample_datasets):
        """Loaded selection preserves correct data types."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="test",
            descriptive_name="test",
        )

        loaded = load_selection_metadata(output_dir)

        assert isinstance(loaded["total"], int)
        assert isinstance(loaded["datasets"], list)
        assert isinstance(loaded["exported_at"], str)


class TestParseDescriptiveNameFromFolder:
    """Tests for parse_descriptive_name_from_folder function."""

    def test_parse_simple_name(self):
        """Parse simple folder name."""
        name = parse_descriptive_name_from_folder("drought-severity_led-monthly")
        assert name == "drought-severity"

    def test_parse_name_with_suffix(self):
        """Parse folder name with -N suffix."""
        name = parse_descriptive_name_from_folder("drought-severity-2_led-monthly")
        assert name == "drought-severity-2"

    def test_parse_name_handles_unknown(self):
        """Parse handles unknown format gracefully."""
        # If format doesn't match, return whole name minus variable-timestep
        name = parse_descriptive_name_from_folder("something_unknown")
        assert name == "something"

    def test_parse_name_extracts_before_variable(self):
        """Extract name before variable code."""
        name = parse_descriptive_name_from_folder("fire-risk_burntarea-annual")
        assert name == "fire-risk"

    def test_parse_complex_name(self):
        """Parse complex folder names."""
        name = parse_descriptive_name_from_folder("long-descriptive-name-2_led-monthly")
        assert name == "long-descriptive-name-2"


class TestSelectionMetadataIntegration:
    """Integration tests for selection metadata functions."""

    def test_save_and_load_roundtrip(self, temp_dir, sample_datasets):
        """Save and load metadata roundtrips correctly."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        original = {
            "query": "drought metrics",
            "descriptive_name": "drought-severity",
            "datasets": sample_datasets,
        }

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query=original["query"],
            descriptive_name=original["descriptive_name"],
        )

        loaded = load_selection_metadata(output_dir)

        assert loaded["query"] == original["query"]
        assert loaded["descriptive_name"] == original["descriptive_name"]
        assert len(loaded["datasets"]) == len(sample_datasets)

    def test_parse_name_from_saved_folder(self, temp_dir, sample_datasets):
        """Parse folder name from saved metadata workflow."""
        folder_name = "drought-severity_led-monthly"
        output_dir = temp_dir / folder_name
        output_dir.mkdir()

        save_selection_metadata(
            output_dir,
            sample_datasets,
            query="drought",
            descriptive_name="drought-severity",
        )

        parsed_name = parse_descriptive_name_from_folder(folder_name)
        loaded = load_selection_metadata(output_dir)

        assert parsed_name == "drought-severity"
        assert loaded["descriptive_name"] == parsed_name


class TestInteractiveWorkflowHelpers:
    """Tests for interactive workflow helper functions."""

    def test_selection_metadata_preserves_all_dataset_fields(self, temp_dir):
        """All dataset fields are preserved in selection metadata."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        dataset = DatasetInfo(
            id="123",
            name="test.nc",
            url="https://example.com/test.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp126",
            variable="led",
            model="gfdl-esm4",
            timestep="monthly",
            size=5000000,
        )

        save_selection_metadata(
            output_dir,
            [dataset],
            query="test",
            descriptive_name="test",
        )

        loaded = load_selection_metadata(output_dir)
        loaded_ds = loaded["datasets"][0]

        assert loaded_ds["id"] == "123"
        assert loaded_ds["variable"] == "led"
        assert loaded_ds["timestep"] == "monthly"
        assert loaded_ds["climate_scenario"] == "ssp126"
        assert loaded_ds["url"] == "https://example.com/test.nc"

    def test_multiple_datasets_in_selection(self, temp_dir):
        """Handle multiple datasets in selection."""
        output_dir = temp_dir / "test_output"
        output_dir.mkdir()

        datasets = [
            DatasetInfo(
                id=str(i),
                name=f"dataset_{i}.nc",
                url=f"https://example.com/dataset_{i}.nc",
                simulation_round="ISIMIP3b",
                climate_scenario=f"ssp{i}",
                variable="led",
                model="gfdl-esm4",
                timestep="monthly",
                size=1000000,
            )
            for i in range(1, 4)
        ]

        save_selection_metadata(
            output_dir,
            datasets,
            query="test",
            descriptive_name="test",
        )

        loaded = load_selection_metadata(output_dir)

        assert len(loaded["datasets"]) == 3
        assert loaded["total"] == 3
