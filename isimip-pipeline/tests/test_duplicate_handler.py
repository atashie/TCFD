"""Tests for duplicate_handler module."""

import re
from enum import Enum
from datetime import datetime

import pytest

from isimip_pipeline.processing_log import DatasetEntry, ProcessingLog
from isimip_pipeline.duplicate_handler import (
    DuplicateAction,
    check_for_duplicate,
    generate_unique_name,
)


@pytest.fixture
def sample_log():
    """Create a processing log with sample entries."""
    log = ProcessingLog()

    entry1 = DatasetEntry(
        descriptive_name="drought-severity",
        variable="led",
        timestep="monthly",
        created_date=datetime(2026, 1, 10, 14, 30, 0),
        output_path="./outputs/drought-severity_led-monthly",
        file_count=15,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126", "ssp370"],
        gcm_models=["gfdl-esm4"],
        lsm_models=["clm5.0"],
        simulation_round="ISIMIP3b",
        query="drought metrics",
    )
    log.add_entry(entry1)

    entry2 = DatasetEntry(
        descriptive_name="fire-risk",
        variable="burntarea",
        timestep="annual",
        created_date=datetime(2026, 1, 11, 10, 0, 0),
        output_path="./outputs/fire-risk_burntarea-annual",
        file_count=10,
        time_periods=["2006-2100"],
        climate_scenarios=["ssp126"],
        gcm_models=["gfdl-esm4"],
        lsm_models=["clm5.0"],
        simulation_round="ISIMIP3b",
        query="fire metrics",
    )
    log.add_entry(entry2)

    return log


class TestDuplicateAction:
    """Tests for DuplicateAction enum."""

    def test_duplicate_action_values(self):
        """DuplicateAction enum has expected values."""
        assert hasattr(DuplicateAction, "SKIP")
        assert hasattr(DuplicateAction, "NEW_FOLDER")
        assert hasattr(DuplicateAction, "OVERWRITE")
        assert hasattr(DuplicateAction, "ABORT")

    def test_duplicate_action_string_values(self):
        """DuplicateAction values are correct strings."""
        assert DuplicateAction.SKIP.value == "skip"
        assert DuplicateAction.NEW_FOLDER.value == "new"
        assert DuplicateAction.OVERWRITE.value == "overwrite"
        assert DuplicateAction.ABORT.value == "abort"


class TestCheckForDuplicate:
    """Tests for check_for_duplicate function."""

    def test_find_duplicate_variable_timestep(self, sample_log):
        """Finds duplicate based on variable+timestep."""
        duplicate = check_for_duplicate(sample_log, "led", "monthly")
        assert duplicate is not None
        assert duplicate.descriptive_name == "drought-severity"

    def test_no_duplicate_different_timestep(self, sample_log):
        """No duplicate for different timestep with same variable."""
        duplicate = check_for_duplicate(sample_log, "led", "annual")
        assert duplicate is None

    def test_no_duplicate_different_variable(self, sample_log):
        """No duplicate for different variable."""
        duplicate = check_for_duplicate(sample_log, "potevap", "monthly")
        assert duplicate is None

    def test_no_duplicate_empty_log(self):
        """Returns None for empty log."""
        log = ProcessingLog()
        duplicate = check_for_duplicate(log, "led", "monthly")
        assert duplicate is None

    def test_multiple_entries_returns_first_match(self, sample_log):
        """Returns first match when multiple could exist."""
        # Add another led-monthly entry
        entry3 = DatasetEntry(
            descriptive_name="drought-metrics-2",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path="./outputs/drought-2_led-monthly",
            file_count=10,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126"],
            gcm_models=["gfdl-esm4"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="drought",
        )
        sample_log.add_entry(entry3)

        duplicate = check_for_duplicate(sample_log, "led", "monthly")
        assert duplicate is not None
        # Should get the first one
        assert duplicate.descriptive_name == "drought-severity"


class TestGenerateUniqueName:
    """Tests for generate_unique_name function."""

    def test_generate_unique_name_appends_suffix(self, sample_log):
        """Generates unique name by appending -2 suffix."""
        unique = generate_unique_name("drought", "led", "monthly", sample_log)
        assert unique.endswith("_led-monthly")
        assert "drought" in unique
        assert "-2" in unique

    def test_generate_unique_name_different_variable(self, sample_log):
        """Generates unique name for different variable (no existing)."""
        unique = generate_unique_name("new-drought", "led", "annual", sample_log)
        # Should not have -2 since no duplicate
        assert unique.endswith("_led-annual")
        assert "new-drought" in unique

    def test_generated_name_not_in_log(self, sample_log):
        """Generated name doesn't conflict with existing entries."""
        unique = generate_unique_name("drought", "led", "monthly", sample_log)
        # Parse the generated name
        parts = unique.rsplit("_", 1)
        name_part = parts[0]

        # Check it's not already in log
        duplicates = sample_log.search(name_part)
        # Should only find if it's the original
        if duplicates:
            assert all(d.descriptive_name.startswith("drought") for d in duplicates)

    def test_multiple_conflicts_increments(self, sample_log):
        """Increments number if multiple conflicts exist."""
        # Add drought-severity-2
        entry = DatasetEntry(
            descriptive_name="drought-severity-2",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path="./outputs/drought-severity-2_led-monthly",
            file_count=10,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126"],
            gcm_models=["gfdl-esm4"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="drought",
        )
        sample_log.add_entry(entry)

        unique = generate_unique_name("drought-severity", "led", "monthly", sample_log)
        # Should suggest -3 now
        assert "-3_" in unique or unique.endswith("-3")

    def test_format_follows_naming_convention(self, sample_log):
        """Generated name follows {name}_{variable}-{timestep} format."""
        unique = generate_unique_name("test", "led", "monthly", sample_log)
        # Pattern: name_variable-timestep or name-N_variable-timestep
        assert re.match(r".*_\w+-\w+$", unique)

    def test_cleans_descriptive_name(self, sample_log):
        """Cleans special characters from descriptive name."""
        unique = generate_unique_name("Test Drought!", "led", "monthly", sample_log)
        # Should be lowercase and no special chars
        name_part = unique.split("_")[0]
        assert name_part.islower()
        assert "!" not in name_part
        assert " " not in name_part


class TestDuplicateHandling:
    """Integration tests for duplicate handling."""

    def test_detect_and_generate_alternative(self, sample_log):
        """Workflow: detect duplicate and generate alternative name."""
        # Check for duplicate
        dup = check_for_duplicate(sample_log, "led", "monthly")
        assert dup is not None

        # Generate alternative
        alt_name = generate_unique_name(dup.descriptive_name, "led", "monthly", sample_log)
        assert alt_name != f"{dup.descriptive_name}_led-monthly"
        assert "_led-monthly" in alt_name

    def test_duplicate_with_existing_numbered_variant(self, sample_log):
        """Handle case where numbered variant already exists."""
        # Add drought-severity-2_led-monthly
        sample_log.add_entry(DatasetEntry(
            descriptive_name="drought-severity-2",
            variable="led",
            timestep="monthly",
            created_date=datetime.now(),
            output_path="./outputs/drought-severity-2_led-monthly",
            file_count=10,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126"],
            gcm_models=["gfdl-esm4"],
            lsm_models=["clm5.0"],
            simulation_round="ISIMIP3b",
            query="drought",
        ))

        # Try to generate name for drought-severity_led-monthly
        dup = check_for_duplicate(sample_log, "led", "monthly")
        assert dup.descriptive_name == "drought-severity"

        alt = generate_unique_name(dup.descriptive_name, "led", "monthly", sample_log)
        # Should suggest -3
        assert "-3_" in alt or alt.endswith("-3")
