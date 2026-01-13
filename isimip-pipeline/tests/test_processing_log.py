"""Tests for processing_log module."""

import tempfile
from pathlib import Path
from datetime import datetime
from typing import Optional

import pytest
import yaml

from isimip_pipeline.processing_log import (
    DatasetEntry,
    ProcessingLog,
    load_processing_log,
    save_processing_log,
    add_dataset_entry,
    find_duplicate,
    search_log,
)


@pytest.fixture
def temp_log_file():
    """Create a temporary log file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        path = Path(f.name)
    yield path
    # Cleanup
    if path.exists():
        path.unlink()


@pytest.fixture
def sample_entry() -> DatasetEntry:
    """Create a sample dataset entry."""
    return DatasetEntry(
        descriptive_name="drought-severity",
        variable="led",
        timestep="monthly",
        created_date=datetime(2026, 1, 13, 10, 30, 0),
        output_path="./outputs/drought-severity_led-monthly",
        file_count=15,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126", "ssp370", "ssp585"],
        gcm_models=["gfdl-esm4", "ukesm1-0-ll"],
        lsm_models=["clm5.0", "lpj-guess"],
        simulation_round="ISIMIP3b",
        query="drought exposure metrics",
    )


class TestDatasetEntry:
    """Tests for DatasetEntry dataclass."""

    def test_create_entry(self, sample_entry):
        """Entry can be created with all fields."""
        assert sample_entry.descriptive_name == "drought-severity"
        assert sample_entry.variable == "led"
        assert sample_entry.timestep == "monthly"
        assert sample_entry.file_count == 15

    def test_entry_to_dict(self, sample_entry):
        """Entry can be converted to dict."""
        data = sample_entry.to_dict()
        assert data["variable"] == "led"
        assert data["timestep"] == "monthly"
        assert data["file_count"] == 15

    def test_entry_from_dict(self, sample_entry):
        """Entry can be created from dict."""
        data = sample_entry.to_dict()
        restored = DatasetEntry.from_dict(data)
        assert restored.variable == sample_entry.variable
        assert restored.timestep == sample_entry.timestep


class TestProcessingLog:
    """Tests for ProcessingLog class."""

    def test_create_empty_log(self):
        """Create empty processing log."""
        log = ProcessingLog()
        assert log.datasets == []
        assert log.last_updated is not None

    def test_add_entry_to_log(self, sample_entry):
        """Add entry to log."""
        log = ProcessingLog()
        log.add_entry(sample_entry)
        assert len(log.datasets) == 1
        assert log.datasets[0].variable == "led"

    def test_find_by_variable_timestep(self, sample_entry):
        """Find duplicate by variable+timestep."""
        log = ProcessingLog()
        log.add_entry(sample_entry)

        # Should find duplicate
        found = log.find_by_variable_timestep("led", "monthly")
        assert found is not None
        assert found.descriptive_name == "drought-severity"

        # Should not find different timestep
        found = log.find_by_variable_timestep("led", "annual")
        assert found is None

        # Should not find different variable
        found = log.find_by_variable_timestep("burntarea", "monthly")
        assert found is None

    def test_search_by_query(self, sample_entry):
        """Search log by query string."""
        log = ProcessingLog()
        log.add_entry(sample_entry)

        # Search by descriptive name
        results = log.search("drought")
        assert len(results) == 1

        # Search by variable
        results = log.search("led")
        assert len(results) == 1

        # No match
        results = log.search("nonexistent")
        assert len(results) == 0

    def test_search_case_insensitive(self, sample_entry):
        """Search is case insensitive."""
        log = ProcessingLog()
        log.add_entry(sample_entry)

        results = log.search("DROUGHT")
        assert len(results) == 1

    def test_multiple_entries(self):
        """Log can hold multiple entries."""
        log = ProcessingLog()

        entry1 = DatasetEntry(
            descriptive_name="drought-1",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path="./outputs/drought-1_led-monthly",
            file_count=15,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126"],
            gcm_models=["gfdl-esm4"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="drought",
        )

        entry2 = DatasetEntry(
            descriptive_name="drought-2",
            variable="led",
            timestep="annual",
            created_date=datetime.now(),
            output_path="./outputs/drought-2_led-annual",
            file_count=10,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126"],
            gcm_models=["gfdl-esm4"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="drought",
        )

        log.add_entry(entry1)
        log.add_entry(entry2)

        assert len(log.datasets) == 2
        # Should find first entry
        found = log.find_by_variable_timestep("led", "monthly")
        assert found.descriptive_name == "drought-1"


class TestLoadSaveLog:
    """Tests for load/save functionality."""

    def test_save_log_to_file(self, temp_log_file, sample_entry):
        """Save log to YAML file."""
        log = ProcessingLog()
        log.add_entry(sample_entry)

        save_processing_log(log, temp_log_file)

        assert temp_log_file.exists()
        with open(temp_log_file) as f:
            data = yaml.safe_load(f)
        assert len(data["datasets"]) == 1
        assert data["datasets"][0]["variable"] == "led"

    def test_load_log_from_file(self, temp_log_file, sample_entry):
        """Load log from YAML file."""
        log = ProcessingLog()
        log.add_entry(sample_entry)
        save_processing_log(log, temp_log_file)

        # Load it back
        loaded_log = load_processing_log(temp_log_file)
        assert len(loaded_log.datasets) == 1
        assert loaded_log.datasets[0].variable == "led"
        assert loaded_log.datasets[0].descriptive_name == "drought-severity"

    def test_load_nonexistent_file_creates_empty(self, temp_log_file):
        """Loading nonexistent file returns empty log."""
        # Remove the file if it exists
        if temp_log_file.exists():
            temp_log_file.unlink()

        log = load_processing_log(temp_log_file)
        assert len(log.datasets) == 0

    def test_roundtrip(self, temp_log_file):
        """Data roundtrips correctly through save/load."""
        original = ProcessingLog()
        entry1 = DatasetEntry(
            descriptive_name="test-1",
            variable="led",
            timestep="monthly",
            created_date=datetime(2026, 1, 1, 12, 0, 0),
            output_path="./outputs/test-1_led-monthly",
            file_count=20,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126", "ssp370"],
            gcm_models=["gfdl-esm4", "ukesm1"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="test",
        )
        original.add_entry(entry1)

        # Save and load
        save_processing_log(original, temp_log_file)
        loaded = load_processing_log(temp_log_file)

        # Verify all fields
        assert len(loaded.datasets) == 1
        e = loaded.datasets[0]
        assert e.descriptive_name == "test-1"
        assert e.variable == "led"
        assert e.timestep == "monthly"
        assert e.file_count == 20
        assert e.time_periods == ["2006-2100"]
        assert e.climate_scenarios == ["ssp126", "ssp370"]


class TestDuplicateFunctions:
    """Tests for module-level convenience functions."""

    def test_find_duplicate_in_log(self, temp_log_file, sample_entry):
        """find_duplicate convenience function works."""
        log = ProcessingLog()
        log.add_entry(sample_entry)
        save_processing_log(log, temp_log_file)

        # Load and check
        loaded = load_processing_log(temp_log_file)
        found = find_duplicate(loaded, "led", "monthly")
        assert found is not None

    def test_add_entry_to_log_file(self, temp_log_file, sample_entry):
        """add_dataset_entry adds to log and saves."""
        log = ProcessingLog()
        updated_log = add_dataset_entry(log, sample_entry)

        assert len(updated_log.datasets) == 1
        assert updated_log.datasets[0].variable == "led"

    def test_search_log_function(self, temp_log_file, sample_entry):
        """search_log convenience function works."""
        log = ProcessingLog()
        log.add_entry(sample_entry)
        save_processing_log(log, temp_log_file)

        loaded = load_processing_log(temp_log_file)
        results = search_log(loaded, "drought")
        assert len(results) == 1


class TestYAMLFormat:
    """Tests for YAML format correctness."""

    def test_yaml_structure(self, temp_log_file, sample_entry):
        """YAML file has expected structure."""
        log = ProcessingLog()
        log.add_entry(sample_entry)
        save_processing_log(log, temp_log_file)

        with open(temp_log_file) as f:
            data = yaml.safe_load(f)

        # Top-level keys
        assert "datasets" in data
        assert "last_updated" in data

        # Dataset fields
        ds = data["datasets"][0]
        assert ds["variable"] == "led"
        assert ds["timestep"] == "monthly"
        assert ds["descriptive_name"] == "drought-severity"
        assert isinstance(ds["climate_scenarios"], list)

    def test_empty_log_yaml(self, temp_log_file):
        """Empty log creates valid YAML."""
        log = ProcessingLog()
        save_processing_log(log, temp_log_file)

        with open(temp_log_file) as f:
            data = yaml.safe_load(f)

        assert "datasets" in data
        assert len(data["datasets"]) == 0
