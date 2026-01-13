"""Tests for the data processor module."""

import pytest
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

from isimip_pipeline.processing.processor import (
    DataProcessor,
    find_netcdf_files,
    group_files_by_variable,
    load_and_aggregate,
)


@pytest.fixture
def sample_netcdf_dir(tmp_path):
    """Create sample NetCDF files for testing."""
    # Create sample data
    lat = np.arange(-90, 91, 30.0)
    lon = np.arange(-180, 180, 30.0)
    time = np.arange(0, 12)  # 12 months

    # Create two sample files for different models
    for model in ["gfdl-esm4", "ukesm1-0-ll"]:
        data = np.random.random((len(lon), len(lat), len(time))).astype(np.float32)

        ds = xr.Dataset(
            {"burntarea": (["lon", "lat", "time"], data)},
            coords={"lon": lon, "lat": lat, "time": time},
        )
        ds.attrs["climate_scenario"] = "ssp585"
        ds.attrs["simulation_round"] = "ISIMIP3b"

        filename = f"test_{model}_ssp585_burntarea_monthly.nc"
        ds.to_netcdf(tmp_path / filename)

    return tmp_path


class TestFindNetCDFFiles:
    """Test NetCDF file discovery."""

    def test_finds_nc_files(self, sample_netcdf_dir):
        """Should find .nc files in directory."""
        files = find_netcdf_files(sample_netcdf_dir)

        assert len(files) == 2
        assert all(f.suffix == ".nc" for f in files)

    def test_finds_nc4_files(self, tmp_path):
        """Should also find .nc4 files."""
        # Create a .nc4 file
        ds = xr.Dataset({"test": (["x"], [1, 2, 3])})
        ds.to_netcdf(tmp_path / "test.nc4")

        files = find_netcdf_files(tmp_path)

        assert len(files) == 1
        assert files[0].suffix == ".nc4"

    def test_returns_empty_for_no_files(self, tmp_path):
        """Should return empty list if no NetCDF files."""
        files = find_netcdf_files(tmp_path)

        assert files == []


class TestGroupFilesByVariable:
    """Test file grouping functionality."""

    def test_groups_by_variable_name(self, sample_netcdf_dir):
        """Should group files by variable name in filename."""
        files = find_netcdf_files(sample_netcdf_dir)
        groups = group_files_by_variable(files)

        assert "burntarea" in groups
        assert len(groups["burntarea"]) == 2


class TestLoadAndAggregate:
    """Test data loading and aggregation."""

    def test_loads_netcdf_files(self, sample_netcdf_dir):
        """Should load NetCDF files into xarray."""
        files = find_netcdf_files(sample_netcdf_dir)

        ds = load_and_aggregate(files, variable="burntarea")

        assert ds is not None
        assert "burntarea" in ds.data_vars

    def test_aggregates_monthly_to_yearly(self, sample_netcdf_dir):
        """Should aggregate monthly data to yearly."""
        files = find_netcdf_files(sample_netcdf_dir)

        ds = load_and_aggregate(files, variable="burntarea", aggregate="yearly")

        # Original has 12 months, aggregated should have fewer time steps
        assert ds is not None


class TestDataProcessor:
    """Test DataProcessor class."""

    def test_processor_initializes(self):
        """DataProcessor should initialize with config."""
        processor = DataProcessor(bandwidth=15)

        assert processor is not None
        assert processor.bandwidth == 15

    def test_processor_has_process_method(self):
        """DataProcessor should have process method."""
        processor = DataProcessor()

        assert hasattr(processor, "process")

    def test_process_returns_dataset(self, sample_netcdf_dir):
        """Process should return xarray Dataset with features."""
        processor = DataProcessor()

        result = processor.process(
            input_dir=sample_netcdf_dir,
            variable="burntarea",
        )

        assert isinstance(result, xr.Dataset)

    def test_process_extracts_features(self, sample_netcdf_dir):
        """Process should extract all 6 feature types."""
        processor = DataProcessor()

        result = processor.process(
            input_dir=sample_netcdf_dir,
            variable="burntarea",
        )

        # Should have value_class dimension with 6 types
        assert "value_class" in result.dims
        assert len(result.value_class) == 6
