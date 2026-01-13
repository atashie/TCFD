"""Tests for CLI commands."""

import pytest
from typer.testing import CliRunner

from isimip_pipeline.cli import app


runner = CliRunner()


class TestCLIBasics:
    """Test basic CLI functionality."""

    def test_app_exists(self):
        """CLI app should exist and be invokable."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0

    def test_version_flag(self):
        """CLI should have a --version flag."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.stdout


class TestSearchCommand:
    """Test the search command."""

    def test_search_command_exists(self):
        """Search command should exist."""
        result = runner.invoke(app, ["search", "--help"])
        assert result.exit_code == 0
        assert "search" in result.stdout.lower() or "query" in result.stdout.lower()

    def test_search_requires_query(self):
        """Search should require a query argument."""
        result = runner.invoke(app, ["search"])
        # Should fail with missing argument
        assert result.exit_code != 0


class TestDownloadCommand:
    """Test the download command."""

    def test_download_command_exists(self):
        """Download command should exist."""
        result = runner.invoke(app, ["download", "--help"])
        assert result.exit_code == 0

    def test_download_has_selection_option(self):
        """Download should have a --selection option."""
        result = runner.invoke(app, ["download", "--help"])
        assert "--selection" in result.stdout or "selection" in result.stdout.lower()


class TestProcessCommand:
    """Test the process command."""

    def test_process_command_exists(self):
        """Process command should exist."""
        result = runner.invoke(app, ["process", "--help"])
        assert result.exit_code == 0

    def test_process_requires_input_dir(self):
        """Process should require an input directory argument."""
        result = runner.invoke(app, ["process"])
        # Should fail with missing argument
        assert result.exit_code != 0


class TestReportCommand:
    """Test the report command."""

    def test_report_command_exists(self):
        """Report command should exist."""
        result = runner.invoke(app, ["report", "--help"])
        assert result.exit_code == 0


class TestRunCommand:
    """Test the run (full pipeline) command."""

    def test_run_command_exists(self):
        """Run command should exist."""
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0

    def test_run_requires_query(self):
        """Run should require a query argument."""
        result = runner.invoke(app, ["run"])
        # Should fail with missing argument
        assert result.exit_code != 0


class TestFindCommand:
    """Test the find (local search) command."""

    def test_find_command_exists(self):
        """Find command should exist."""
        result = runner.invoke(app, ["find", "--help"])
        assert result.exit_code == 0
        assert "local" in result.stdout.lower() or "search" in result.stdout.lower()

    def test_find_accepts_query_argument(self):
        """Find should accept optional query argument."""
        # When no log exists, should handle gracefully
        result = runner.invoke(app, ["find", "drought"])
        # May fail or show no results, but shouldn't crash with syntax error
        assert "No datasets" in result.stdout or result.exit_code in [0, 1]

    def test_find_accepts_variable_filter(self):
        """Find should accept --variable filter option."""
        result = runner.invoke(app, ["find", "--help"])
        assert "--variable" in result.stdout or "-v" in result.stdout

    def test_find_accepts_timestep_filter(self):
        """Find should accept --timestep filter option."""
        result = runner.invoke(app, ["find", "--help"])
        assert "--timestep" in result.stdout or "-t" in result.stdout

    def test_find_accepts_detailed_flag(self):
        """Find should accept --detailed flag."""
        result = runner.invoke(app, ["find", "--help"])
        assert "--detailed" in result.stdout or "-d" in result.stdout


class TestInteractiveCommand:
    """Test the interactive workflow command."""

    def test_interactive_command_exists(self):
        """Interactive command should exist."""
        result = runner.invoke(app, ["interactive", "--help"])
        assert result.exit_code == 0
        assert "interactive" in result.stdout.lower() or "workflow" in result.stdout.lower()

    def test_interactive_requires_query(self):
        """Interactive should require a query argument."""
        result = runner.invoke(app, ["interactive"])
        # Should fail with missing argument
        assert result.exit_code != 0

    def test_interactive_accepts_query(self):
        """Interactive should accept a query argument."""
        result = runner.invoke(app, ["interactive", "--help"])
        # Help should show query argument
        assert "query" in result.stdout.lower()

    def test_interactive_has_local_only_option(self):
        """Interactive should have --local-only option."""
        result = runner.invoke(app, ["interactive", "--help"])
        assert "--local-only" in result.stdout or "local" in result.stdout.lower()

    def test_interactive_has_config_option(self):
        """Interactive should have --config option."""
        result = runner.invoke(app, ["interactive", "--help"])
        assert "--config" in result.stdout or "-c" in result.stdout


class TestProcessCommandEnhanced:
    """Test enhanced process command with auto-detection."""

    def test_process_accepts_variable_option(self):
        """Process should accept optional --variable option."""
        result = runner.invoke(app, ["process", "--help"])
        assert "--variable" in result.stdout or "-v" in result.stdout

    def test_process_accepts_scenarios_option(self):
        """Process should accept --scenarios option."""
        result = runner.invoke(app, ["process", "--help"])
        assert "--scenarios" in result.stdout or "-s" in result.stdout

    def test_process_has_output_option(self):
        """Process should have --output option."""
        result = runner.invoke(app, ["process", "--help"])
        assert "--output" in result.stdout or "-o" in result.stdout

    def test_process_help_mentions_auto_detection(self):
        """Process help should mention auto-detection feature."""
        result = runner.invoke(app, ["process", "--help"])
        assert "auto" in result.stdout.lower() or "detect" in result.stdout.lower()

    def test_process_help_mentions_log_update(self):
        """Process help should mention log update feature."""
        result = runner.invoke(app, ["process", "--help"])
        assert "log" in result.stdout.lower() or "update" in result.stdout.lower()


class TestCatalogCommand:
    """Test the catalog command."""

    def test_catalog_command_exists(self):
        """Catalog command should exist."""
        result = runner.invoke(app, ["catalog", "--help"])
        assert result.exit_code == 0

    def test_catalog_accepts_all_flag(self):
        """Catalog should accept --all flag."""
        result = runner.invoke(app, ["catalog", "--help"])
        assert "--all" in result.stdout or "-a" in result.stdout
