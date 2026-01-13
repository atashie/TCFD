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
    group_by_variable_timestep,
    display_grouped_results,
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


class TestGroupByVariableTimestep:
    """Test grouping datasets by variable and timestep."""

    def test_group_by_variable_timestep_basic(self, sample_datasets):
        """Should group datasets by (variable, timestep) tuple."""
        grouped = group_by_variable_timestep(sample_datasets)

        # Should have two groups: led-annual and led-monthly
        assert len(grouped) == 2
        assert ("led", "annual") in grouped
        assert ("led", "monthly") in grouped

    def test_group_by_variable_timestep_counts(self, sample_datasets):
        """Should count files correctly in each group."""
        grouped = group_by_variable_timestep(sample_datasets)

        # led-annual: 2 files (ssp126, ssp370)
        # led-monthly: 1 file (ssp126)
        assert grouped[("led", "annual")]["file_count"] == 2
        assert grouped[("led", "monthly")]["file_count"] == 1

    def test_group_by_variable_timestep_scenarios(self, sample_datasets):
        """Should extract unique scenarios from each group."""
        grouped = group_by_variable_timestep(sample_datasets)

        # led-annual has ssp126 and ssp370
        scenarios_annual = grouped[("led", "annual")]["scenarios"]
        assert "ssp126" in scenarios_annual
        assert "ssp370" in scenarios_annual

        # led-monthly has only ssp126
        scenarios_monthly = grouped[("led", "monthly")]["scenarios"]
        assert "ssp126" in scenarios_monthly

    def test_group_by_variable_timestep_models(self, sample_datasets):
        """Should extract unique models from each group."""
        grouped = group_by_variable_timestep(sample_datasets)

        # led-annual has gfdl-esm4
        models_annual = grouped[("led", "annual")]["models"]
        assert "gfdl-esm4" in models_annual

        # led-monthly has ukesm1-0-ll
        models_monthly = grouped[("led", "monthly")]["models"]
        assert "ukesm1-0-ll" in models_monthly

    def test_group_by_variable_timestep_simulation_rounds(self, sample_datasets):
        """Should extract simulation rounds from each group."""
        grouped = group_by_variable_timestep(sample_datasets)

        # All datasets are ISIMIP3b
        assert "ISIMIP3b" in grouped[("led", "annual")]["simulation_rounds"]
        assert "ISIMIP3b" in grouped[("led", "monthly")]["simulation_rounds"]

    def test_group_by_variable_timestep_empty_list(self):
        """Should handle empty dataset list."""
        grouped = group_by_variable_timestep([])

        assert len(grouped) == 0

    def test_group_by_variable_timestep_single_dataset(self):
        """Should handle single dataset."""
        dataset = DatasetInfo(
            id="1",
            name="test.nc",
            url="https://example.com/test.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp126",
            variable="led",
            model="gfdl-esm4",
            timestep="monthly",
        )

        grouped = group_by_variable_timestep([dataset])

        assert len(grouped) == 1
        assert ("led", "monthly") in grouped
        assert grouped[("led", "monthly")]["file_count"] == 1

    def test_group_by_variable_timestep_multiple_variables(self):
        """Should correctly group multiple different variables."""
        datasets = [
            DatasetInfo(
                id="1",
                name="led_annual.nc",
                url="https://example.com/led_annual.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                variable="led",
                model="gfdl-esm4",
                timestep="annual",
            ),
            DatasetInfo(
                id="2",
                name="burntarea_monthly.nc",
                url="https://example.com/burntarea_monthly.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                variable="burntarea",
                model="gfdl-esm4",
                timestep="monthly",
            ),
            DatasetInfo(
                id="3",
                name="led_monthly.nc",
                url="https://example.com/led_monthly.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                variable="led",
                model="gfdl-esm4",
                timestep="monthly",
            ),
        ]

        grouped = group_by_variable_timestep(datasets)

        assert len(grouped) == 3
        assert ("led", "annual") in grouped
        assert ("burntarea", "monthly") in grouped
        assert ("led", "monthly") in grouped


class TestDisplayGroupedResults:
    """Test displaying grouped results to terminal."""

    def test_display_grouped_results_with_data(self, sample_datasets, capsys):
        """Should display grouped results without errors."""
        from rich.console import Console

        grouped = group_by_variable_timestep(sample_datasets)
        console = Console()

        # Should not raise an error
        try:
            display_grouped_results(grouped, console)
        except Exception as e:
            pytest.fail(f"display_grouped_results raised {type(e).__name__}: {e}")

    def test_display_grouped_results_empty(self, capsys):
        """Should handle empty grouped results."""
        from rich.console import Console

        console = Console()

        # Should not raise an error with empty dict
        try:
            display_grouped_results({}, console)
        except Exception as e:
            pytest.fail(f"display_grouped_results raised {type(e).__name__}: {e}")

    def test_display_grouped_results_single_group(self, capsys):
        """Should display single group correctly."""
        from rich.console import Console

        dataset = DatasetInfo(
            id="1",
            name="test.nc",
            url="https://example.com/test.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp126",
            variable="led",
            model="gfdl-esm4",
            timestep="monthly",
        )

        grouped = group_by_variable_timestep([dataset])
        console = Console()

        try:
            display_grouped_results(grouped, console)
        except Exception as e:
            pytest.fail(f"display_grouped_results raised {type(e).__name__}: {e}")

    def test_display_grouped_results_multiple_groups(self, sample_datasets, capsys):
        """Should display multiple groups correctly."""
        from rich.console import Console

        grouped = group_by_variable_timestep(sample_datasets)
        console = Console()

        try:
            display_grouped_results(grouped, console)
        except Exception as e:
            pytest.fail(f"display_grouped_results raised {type(e).__name__}: {e}")
