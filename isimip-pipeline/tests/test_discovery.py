"""Tests for discovery module (local dataset search)."""

import tempfile
from pathlib import Path
from datetime import datetime

import pytest

from isimip_pipeline.processing_log import DatasetEntry, ProcessingLog, save_processing_log
from isimip_pipeline.discovery import (
    find_local_datasets,
    get_dataset_summary,
    verify_dataset_integrity,
)


@pytest.fixture
def sample_log():
    """Create a sample processing log with multiple entries."""
    log = ProcessingLog()

    # Entry 1: drought-monthly
    entry1 = DatasetEntry(
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
    log.add_entry(entry1)

    # Entry 2: drought-annual
    entry2 = DatasetEntry(
        descriptive_name="drought-annual",
        variable="led",
        timestep="annual",
        created_date=datetime(2026, 1, 11, 10, 0, 0),
        output_path="./outputs/drought-annual_led-annual",
        file_count=10,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126", "ssp585"],
        gcm_models=["gfdl-esm4"],
        lsm_models=["clm5.0"],
        simulation_round="ISIMIP3b",
        query="drought annual",
    )
    log.add_entry(entry2)

    # Entry 3: fire-annual
    entry3 = DatasetEntry(
        descriptive_name="fire-risk",
        variable="burntarea",
        timestep="annual",
        created_date=datetime(2026, 1, 12, 9, 15, 0),
        output_path="./outputs/fire-risk_burntarea-annual",
        file_count=8,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126", "ssp370"],
        gcm_models=["gfdl-esm4", "ipsl-cm6a-lr"],
        lsm_models=["caraib", "lpj-guess"],
        simulation_round="ISIMIP3b",
        query="fire burnt area",
    )
    log.add_entry(entry3)

    return log


class TestFindLocalDatasets:
    """Tests for find_local_datasets function."""

    def test_find_all_datasets(self, sample_log):
        """Find all datasets when no filters applied."""
        results = find_local_datasets(sample_log, query=None)
        assert len(results) == 3

    def test_find_by_query(self, sample_log):
        """Find datasets by search query."""
        results = find_local_datasets(sample_log, query="drought")
        assert len(results) == 2
        assert all("drought" in e.descriptive_name.lower() for e in results)

    def test_find_by_variable(self, sample_log):
        """Find datasets by variable name."""
        results = find_local_datasets(sample_log, variable="led")
        assert len(results) == 2
        assert all(e.variable == "led" for e in results)

    def test_find_by_timestep(self, sample_log):
        """Find datasets by timestep."""
        results = find_local_datasets(sample_log, timestep="annual")
        assert len(results) == 2
        assert all(e.timestep == "annual" for e in results)

    def test_find_by_scenario(self, sample_log):
        """Find datasets by climate scenario."""
        results = find_local_datasets(sample_log, scenario="ssp370")
        assert len(results) == 2
        # Both should have ssp370
        assert all("ssp370" in e.climate_scenarios for e in results)

    def test_find_by_multiple_filters(self, sample_log):
        """Find with multiple filters applied (AND logic)."""
        results = find_local_datasets(
            sample_log,
            variable="led",
            timestep="monthly",
        )
        assert len(results) == 1
        assert results[0].descriptive_name == "drought-severity"

    def test_find_by_variable_and_scenario(self, sample_log):
        """Find by variable and scenario."""
        results = find_local_datasets(
            sample_log,
            variable="led",
            scenario="ssp585",
        )
        assert len(results) == 2
        assert all(e.variable == "led" for e in results)

    def test_find_no_matches(self, sample_log):
        """Return empty list when no matches."""
        results = find_local_datasets(sample_log, query="nonexistent")
        assert len(results) == 0

    def test_find_query_case_insensitive(self, sample_log):
        """Query search is case insensitive."""
        results = find_local_datasets(sample_log, query="DROUGHT")
        assert len(results) == 2

    def test_find_by_scenario_no_match(self, sample_log):
        """No results if scenario not in dataset."""
        results = find_local_datasets(sample_log, scenario="ssp999")
        assert len(results) == 0


class TestGetDatasetSummary:
    """Tests for get_dataset_summary function."""

    def test_summary_contains_expected_fields(self, sample_log):
        """Summary has all expected metadata fields."""
        entry = sample_log.datasets[0]
        summary = get_dataset_summary(entry)

        assert summary["descriptive_name"] == "drought-severity"
        assert summary["variable"] == "led"
        assert summary["timestep"] == "monthly"
        assert summary["file_count"] == 15
        assert summary["simulation_round"] == "ISIMIP3b"

    def test_summary_formats_date(self, sample_log):
        """Summary includes formatted date."""
        entry = sample_log.datasets[0]
        summary = get_dataset_summary(entry)

        assert "created_date" in summary
        assert "2026-01-10" in summary["created_date"]

    def test_summary_lists_models(self, sample_log):
        """Summary includes GCM and LSM models."""
        entry = sample_log.datasets[0]
        summary = get_dataset_summary(entry)

        assert "gcm_models" in summary
        assert len(summary["gcm_models"]) == 2
        assert "lsm_models" in summary
        assert len(summary["lsm_models"]) == 2

    def test_summary_lists_scenarios(self, sample_log):
        """Summary includes all climate scenarios."""
        entry = sample_log.datasets[0]
        summary = get_dataset_summary(entry)

        assert "climate_scenarios" in summary
        assert len(summary["climate_scenarios"]) == 3
        assert "ssp126" in summary["climate_scenarios"]

    def test_summary_all_entries(self, sample_log):
        """Can get summary for all entries."""
        for entry in sample_log.datasets:
            summary = get_dataset_summary(entry)
            assert "descriptive_name" in summary
            assert "variable" in summary


class TestVerifyDatasetIntegrity:
    """Tests for verify_dataset_integrity function."""

    def test_verify_missing_directory(self):
        """Returns False if output directory doesn't exist."""
        entry = DatasetEntry(
            descriptive_name="nonexistent",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path="./nonexistent/path",
            file_count=0,
        )

        result = verify_dataset_integrity(entry)
        assert result is False

    def test_verify_existing_directory(self, tmp_path):
        """Returns True if output directory exists."""
        entry = DatasetEntry(
            descriptive_name="test",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path=str(tmp_path),
            file_count=0,
        )

        result = verify_dataset_integrity(entry)
        assert result is True

    def test_verify_with_processed_file(self, tmp_path):
        """Returns True if processed file exists."""
        processed_dir = tmp_path / "processed"
        processed_dir.mkdir()
        (processed_dir / "led_processed.nc").touch()

        entry = DatasetEntry(
            descriptive_name="test",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path=str(tmp_path),
            file_count=0,
        )

        result = verify_dataset_integrity(entry)
        assert result is True

    def test_verify_missing_processed_file(self, tmp_path):
        """Returns False if processed file missing."""
        processed_dir = tmp_path / "processed"
        processed_dir.mkdir()
        # No files created

        entry = DatasetEntry(
            descriptive_name="test",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path=str(tmp_path),
            file_count=1,  # Expect at least 1 file
        )

        result = verify_dataset_integrity(entry)
        assert result is False


class TestDiscoveryIntegration:
    """Integration tests for discovery functions."""

    def test_find_and_summarize(self, sample_log):
        """Can find datasets and summarize results."""
        results = find_local_datasets(sample_log, variable="led")
        assert len(results) == 2

        for entry in results:
            summary = get_dataset_summary(entry)
            assert summary["variable"] == "led"

    def test_find_multiple_scenarios(self, sample_log):
        """Can filter by multiple scenarios."""
        # led has ssp126, ssp370, ssp585
        # burntarea has ssp126, ssp370
        # Both have ssp126 and ssp370

        results_370 = find_local_datasets(sample_log, scenario="ssp370")
        assert len(results_370) == 2

        results_585 = find_local_datasets(sample_log, scenario="ssp585")
        # Both drought-severity and drought-annual have ssp585
        assert len(results_585) == 2
        assert all(e.variable == "led" for e in results_585)
