"""Tests for ISIMIP metrics catalog."""

import pytest
import tempfile
from pathlib import Path
from datetime import datetime

from isimip_pipeline.catalog import (
    ISIMIPCatalog,
    VariableInfo,
    load_catalog,
    save_catalog,
)
from isimip_pipeline.search.isimip_query import DatasetInfo


class TestVariableInfo:
    """Test VariableInfo data structure."""

    def test_variable_info_creates(self):
        """VariableInfo should store variable metadata."""
        info = VariableInfo(
            name="led",
            long_name="Land Exposed to Drought",
            unit="fraction",
        )
        assert info.name == "led"
        assert info.long_name == "Land Exposed to Drought"

    def test_variable_info_tracks_scenarios(self):
        """VariableInfo should track available scenarios."""
        info = VariableInfo(name="led")
        info.scenarios.add("ssp126")
        info.scenarios.add("ssp585")
        assert "ssp126" in info.scenarios
        assert "ssp585" in info.scenarios

    def test_variable_info_tracks_models(self):
        """VariableInfo should track climate models."""
        info = VariableInfo(name="led")
        info.models.add("gfdl-esm4")
        info.models.add("ukesm1-0-ll")
        assert len(info.models) == 2

    def test_variable_info_tracks_file_count(self):
        """VariableInfo should track file count."""
        info = VariableInfo(name="led")
        info.file_count = 42
        assert info.file_count == 42


class TestISIMIPCatalog:
    """Test ISIMIPCatalog class."""

    def test_catalog_initializes_empty(self):
        """Catalog should initialize with empty variables."""
        catalog = ISIMIPCatalog()
        assert len(catalog.variables) == 0

    def test_catalog_add_variable(self):
        """Catalog should allow adding variables."""
        catalog = ISIMIPCatalog()
        catalog.add_variable("led", long_name="Land Exposed to Drought")
        assert "led" in catalog.variables
        assert catalog.variables["led"].long_name == "Land Exposed to Drought"

    def test_catalog_update_from_datasets(self):
        """Catalog should update from DatasetInfo list."""
        catalog = ISIMIPCatalog()
        datasets = [
            DatasetInfo(
                id="1",
                name="model1_ssp126_led.nc",
                url="https://example.com/1.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                model="gfdl-esm4",
                variable="led",
            ),
            DatasetInfo(
                id="2",
                name="model2_ssp585_led.nc",
                url="https://example.com/2.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp585",
                model="ukesm1-0-ll",
                variable="led",
            ),
        ]
        catalog.update_from_datasets(datasets)

        assert "led" in catalog.variables
        assert "ssp126" in catalog.variables["led"].scenarios
        assert "ssp585" in catalog.variables["led"].scenarios
        assert "gfdl-esm4" in catalog.variables["led"].models
        assert catalog.variables["led"].file_count == 2

    def test_catalog_merge_preserves_existing(self):
        """Updating catalog should merge with existing data."""
        catalog = ISIMIPCatalog()
        catalog.add_variable("led")
        catalog.variables["led"].scenarios.add("historical")
        catalog.variables["led"].file_count = 10

        datasets = [
            DatasetInfo(
                id="1",
                name="test.nc",
                url="https://example.com/1.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                variable="led",
            ),
        ]
        catalog.update_from_datasets(datasets)

        # Should have both old and new scenarios
        assert "historical" in catalog.variables["led"].scenarios
        assert "ssp126" in catalog.variables["led"].scenarios
        # File count should be cumulative
        assert catalog.variables["led"].file_count == 11

    def test_catalog_tracks_last_updated(self):
        """Catalog should track when it was last updated."""
        catalog = ISIMIPCatalog()
        catalog.update_from_datasets([])
        assert catalog.last_updated is not None

    def test_catalog_get_summary(self):
        """Catalog should provide summary statistics."""
        catalog = ISIMIPCatalog()
        catalog.add_variable("led", long_name="Drought")
        catalog.variables["led"].scenarios = {"ssp126", "ssp585"}
        catalog.variables["led"].file_count = 20

        catalog.add_variable("burntarea", long_name="Burnt Area")
        catalog.variables["burntarea"].scenarios = {"ssp126"}
        catalog.variables["burntarea"].file_count = 10

        summary = catalog.get_summary()
        assert summary["total_variables"] == 2
        assert summary["total_files"] == 30
        assert "ssp126" in summary["all_scenarios"]
        assert "ssp585" in summary["all_scenarios"]


class TestCatalogPersistence:
    """Test catalog save/load functionality."""

    def test_save_catalog_creates_file(self):
        """save_catalog should create YAML file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "catalog.yaml"
            catalog = ISIMIPCatalog()
            catalog.add_variable("led")

            save_catalog(catalog, path)
            assert path.exists()

    def test_load_catalog_reads_file(self):
        """load_catalog should read YAML file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "catalog.yaml"
            catalog = ISIMIPCatalog()
            catalog.add_variable("led", long_name="Drought Exposure")
            catalog.variables["led"].scenarios = {"ssp126", "ssp585"}
            catalog.variables["led"].models = {"gfdl-esm4"}
            catalog.variables["led"].file_count = 15
            save_catalog(catalog, path)

            loaded = load_catalog(path)
            assert "led" in loaded.variables
            assert loaded.variables["led"].long_name == "Drought Exposure"
            assert "ssp126" in loaded.variables["led"].scenarios
            assert loaded.variables["led"].file_count == 15

    def test_load_catalog_returns_empty_if_missing(self):
        """load_catalog should return empty catalog if file missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "nonexistent.yaml"
            catalog = load_catalog(path)
            assert len(catalog.variables) == 0

    def test_catalog_roundtrip(self):
        """Catalog should survive save/load roundtrip."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "catalog.yaml"

            original = ISIMIPCatalog()
            original.add_variable("led", long_name="Drought", unit="fraction")
            original.variables["led"].scenarios = {"historical", "ssp126", "ssp585"}
            original.variables["led"].models = {"gfdl-esm4", "ukesm1-0-ll"}
            original.variables["led"].simulation_rounds = {"ISIMIP3b"}
            original.variables["led"].file_count = 42

            original.add_variable("burntarea", long_name="Burnt Area", unit="fraction")
            original.variables["burntarea"].scenarios = {"ssp370"}
            original.variables["burntarea"].file_count = 10

            save_catalog(original, path)
            loaded = load_catalog(path)

            assert len(loaded.variables) == 2
            assert loaded.variables["led"].file_count == 42
            assert "ukesm1-0-ll" in loaded.variables["led"].models
            assert loaded.variables["burntarea"].long_name == "Burnt Area"


class TestCatalogDefaultPath:
    """Test default catalog path handling."""

    def test_default_catalog_path(self):
        """Catalog should have default path in ~/.isimip-pipeline/."""
        from isimip_pipeline.catalog import get_default_catalog_path
        path = get_default_catalog_path()
        assert ".isimip-pipeline" in str(path)
        assert path.name == "catalog.yaml"
