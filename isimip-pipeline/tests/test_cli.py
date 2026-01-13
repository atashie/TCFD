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
