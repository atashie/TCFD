"""Tests for configuration loading."""

import pytest
from pathlib import Path
import tempfile
import yaml

from isimip_pipeline.config import Config, load_config, DEFAULT_CONFIG


class TestConfigDefaults:
    """Test that Config provides sensible defaults."""

    def test_default_download_dir(self):
        """Config should have a default download directory."""
        config = Config()
        assert config.paths.download_dir == Path("./data/raw")

    def test_default_processed_dir(self):
        """Config should have a default processed directory."""
        config = Config()
        assert config.paths.processed_dir == Path("./data/processed")

    def test_default_reports_dir(self):
        """Config should have a default reports directory."""
        config = Config()
        assert config.paths.reports_dir == Path("./reports")

    def test_default_smoothing_bandwidth(self):
        """Config should have default processing parameters."""
        config = Config()
        assert config.processing.smoothing_bandwidth == 15

    def test_default_percentile_bins(self):
        """Config should have 100 percentile bins by default."""
        config = Config()
        assert config.processing.percentile_bins == 100


class TestConfigFromFile:
    """Test loading configuration from YAML file."""

    def test_load_config_overrides_defaults(self):
        """Values in config file should override defaults."""
        config_data = {
            "paths": {
                "download_dir": "/custom/download",
                "processed_dir": "/custom/processed",
            },
            "processing": {
                "smoothing_bandwidth": 20,
            },
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            config = load_config(config_path)
            assert config.paths.download_dir == Path("/custom/download")
            assert config.paths.processed_dir == Path("/custom/processed")
            assert config.processing.smoothing_bandwidth == 20
            # Default should still apply for unspecified values
            assert config.paths.reports_dir == Path("./reports")
        finally:
            config_path.unlink()

    def test_load_config_missing_file_returns_defaults(self):
        """Missing config file should return defaults, not error."""
        config = load_config(Path("/nonexistent/config.yaml"))
        assert config.paths.download_dir == Path("./data/raw")

    def test_load_config_with_api_keys(self):
        """Config should load API credentials."""
        config_data = {
            "api": {
                "you_api_key": "test-key-123",
                "you_agent_id": "agent-456",
            }
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            config = load_config(config_path)
            assert config.api.you_api_key == "test-key-123"
            assert config.api.you_agent_id == "agent-456"
        finally:
            config_path.unlink()


class TestConfigValidation:
    """Test configuration validation."""

    def test_invalid_trend_method_raises(self):
        """Invalid trend method should raise ValueError."""
        config_data = {
            "processing": {
                "trend_method": "invalid_method",
            }
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="trend_method"):
                load_config(config_path)
        finally:
            config_path.unlink()

    def test_valid_trend_methods(self):
        """Both theil_sen and ols should be valid trend methods."""
        for method in ["theil_sen", "ols"]:
            config_data = {"processing": {"trend_method": method}}

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False
            ) as f:
                yaml.dump(config_data, f)
                config_path = Path(f.name)

            try:
                config = load_config(config_path)
                assert config.processing.trend_method == method
            finally:
                config_path.unlink()
