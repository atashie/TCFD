"""Tests for processing output module."""

import pytest
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

from isimip_pipeline.processing.output import (
    OutputWriter,
    create_output_dataset,
    write_netcdf,
)


@pytest.fixture
def sample_features():
    """Create sample feature data for testing."""
    return {
        10: {
            "smoothed_median": 0.05,
            "percentile": 50,
            "trend": 0.001,
            "significance": 0.03,
            "lower_bound": 0.03,
            "upper_bound": 0.07,
        },
        20: {
            "smoothed_median": 0.08,
            "percentile": 65,
            "trend": 0.002,
            "significance": 0.01,
            "lower_bound": 0.06,
            "upper_bound": 0.10,
        },
    }


@pytest.fixture
def sample_lat_lon():
    """Create sample lat/lon arrays."""
    lat = np.arange(-90, 91, 1.0)
    lon = np.arange(-180, 180, 1.0)
    return lat, lon


class TestCreateOutputDataset:
    """Test output dataset creation."""

    def test_creates_xarray_dataset(self, sample_features, sample_lat_lon):
        """Should create xarray Dataset."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10, 20],
            scenarios=["ssp126"],
            variable="burntarea",
        )

        assert isinstance(ds, xr.Dataset)

    def test_has_correct_dimensions(self, sample_lat_lon):
        """Dataset should have lat, lon, decade, scenario, value_class dims."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10, 20, 30],
            scenarios=["ssp126", "ssp585"],
            variable="burntarea",
        )

        assert "lat" in ds.dims
        assert "lon" in ds.dims
        assert "decade" in ds.dims
        assert "scenario" in ds.dims
        assert "value_class" in ds.dims

    def test_value_class_has_6_types(self, sample_lat_lon):
        """Dataset should have 6 value classes."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10],
            scenarios=["ssp126"],
            variable="burntarea",
        )

        assert len(ds.value_class) == 6


class TestOutputWriter:
    """Test OutputWriter class."""

    def test_writer_initializes(self):
        """OutputWriter should initialize with output path."""
        writer = OutputWriter(output_dir=Path("./output"))

        assert writer is not None
        assert writer.output_dir == Path("./output")

    def test_writer_has_write_method(self):
        """OutputWriter should have write method."""
        writer = OutputWriter(output_dir=Path("./output"))

        assert hasattr(writer, "write")

    def test_writer_sets_attributes(self, sample_lat_lon):
        """Writer should set CF-compliant attributes."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10],
            scenarios=["ssp126"],
            variable="burntarea",
        )

        # Check for CF metadata
        assert "long_name" in ds.lat.attrs or "units" in ds.lat.attrs
        assert "long_name" in ds.lon.attrs or "units" in ds.lon.attrs


class TestWriteNetCDF:
    """Test NetCDF writing function."""

    def test_writes_netcdf_file(self, sample_lat_lon):
        """Should write NetCDF file to disk."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10],
            scenarios=["ssp126"],
            variable="test_var",
        )

        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_netcdf(ds, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        finally:
            output_path.unlink()

    def test_written_file_is_readable(self, sample_lat_lon):
        """Written file should be readable with xarray."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10, 20],
            scenarios=["ssp126"],
            variable="test_var",
        )

        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_netcdf(ds, output_path)

            # Re-read and verify
            ds_read = xr.open_dataset(output_path)
            assert "lat" in ds_read.dims
            assert "lon" in ds_read.dims
            ds_read.close()
        finally:
            output_path.unlink()

    def test_compression_enabled(self, sample_lat_lon):
        """Written file should use compression."""
        lat, lon = sample_lat_lon

        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=[10],
            scenarios=["ssp126"],
            variable="test_var",
        )

        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_netcdf(ds, output_path, compression=True)
            # File should exist and be smaller than uncompressed
            assert output_path.exists()
        finally:
            output_path.unlink()
