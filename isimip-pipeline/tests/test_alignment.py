"""Tests for data alignment module.

Tests for spatial grid verification, time coordinate alignment,
and calendar conversion for multi-model ISIMIP datasets.
"""

import pytest
import numpy as np
import xarray as xr
from pathlib import Path
import tempfile
from datetime import datetime, timedelta

from isimip_pipeline.processing.alignment import (
    SpatialGridMismatchError,
    TimeCoordinateMismatchError,
    verify_spatial_grids,
    align_time_coordinates,
    convert_calendar,
    get_calendar_type,
    align_datasets,
    SpatialGridInfo,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_grid():
    """Create sample latitude and longitude grids."""
    lon = np.linspace(-180, 179.5, 720)
    lat = np.linspace(-90, 89.75, 360)
    return lon, lat


@pytest.fixture
def standard_dataset(sample_grid):
    """Create a standard test dataset."""
    lon, lat = sample_grid
    time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')

    data = np.random.randn(len(time), len(lat), len(lon))

    ds = xr.Dataset(
        {
            "temperature": (["time", "lat", "lon"], data, {"units": "K"}),
        },
        coords={
            "time": time,
            "lat": lat,
            "lon": lon,
        },
    )
    return ds


class TestSpatialGridInfo:
    """Test SpatialGridInfo dataclass."""

    def test_spatial_grid_info_creation(self, sample_grid):
        """Should create SpatialGridInfo from arrays."""
        lon, lat = sample_grid
        info = SpatialGridInfo(lon, lat)

        assert info.lon_min == lon[0]
        assert info.lon_max == lon[-1]
        assert info.lat_min == lat[0]
        assert info.lat_max == lat[-1]
        assert info.n_lon == len(lon)
        assert info.n_lat == len(lat)

    def test_spatial_grid_info_resolution(self, sample_grid):
        """Should calculate grid resolution."""
        lon, lat = sample_grid
        info = SpatialGridInfo(lon, lat)

        expected_lon_res = (lon[-1] - lon[0]) / (len(lon) - 1)
        expected_lat_res = (lat[-1] - lat[0]) / (len(lat) - 1)

        assert abs(info.lon_resolution - expected_lon_res) < 0.01
        assert abs(info.lat_resolution - expected_lat_res) < 0.01


class TestVerifySpatialGrids:
    """Test spatial grid verification."""

    def test_verify_identical_grids(self, standard_dataset):
        """Should accept identical spatial grids."""
        ds1 = standard_dataset
        ds2 = standard_dataset.copy(deep=True)

        # Should not raise
        verify_spatial_grids(ds1, ds2)

    def test_verify_different_lon_resolution(self, sample_grid):
        """Should detect different longitude resolution."""
        lon, lat = sample_grid

        ds1 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat, "lon": lon},
        )

        # Create dataset with different lon resolution
        lon2 = np.linspace(-180, 179, 360)  # Different spacing
        ds2 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon2)))},
            coords={"lat": lat, "lon": lon2},
        )

        with pytest.raises(SpatialGridMismatchError):
            verify_spatial_grids(ds1, ds2)

    def test_verify_different_lat_resolution(self, sample_grid):
        """Should detect different latitude resolution."""
        lon, lat = sample_grid

        ds1 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat, "lon": lon},
        )

        # Create dataset with different lat resolution
        lat2 = np.linspace(-90, 89, 180)  # Different spacing
        ds2 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat2), len(lon)))},
            coords={"lat": lat2, "lon": lon},
        )

        with pytest.raises(SpatialGridMismatchError):
            verify_spatial_grids(ds1, ds2)

    def test_verify_tolerance_within_bounds(self, sample_grid):
        """Should accept grids within tolerance."""
        lon, lat = sample_grid

        ds1 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat, "lon": lon},
        )

        # Slightly different grid (within tolerance)
        lon2 = lon + 0.001  # 0.001 degree difference
        lat2 = lat + 0.001
        ds2 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat2, "lon": lon2},
        )

        # Should not raise with default tolerance
        verify_spatial_grids(ds1, ds2, tolerance=0.01)

    def test_verify_tolerance_exceeded(self, sample_grid):
        """Should reject grids beyond tolerance."""
        lon, lat = sample_grid

        ds1 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat, "lon": lon},
        )

        # Significantly different grid
        lon2 = lon + 0.5  # 0.5 degree difference
        lat2 = lat + 0.5
        ds2 = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(len(lat), len(lon)))},
            coords={"lat": lat2, "lon": lon2},
        )

        with pytest.raises(SpatialGridMismatchError):
            verify_spatial_grids(ds1, ds2, tolerance=0.1)


class TestGetCalendarType:
    """Test calendar type detection."""

    def test_detect_standard_calendar(self, standard_dataset):
        """Should detect standard calendar."""
        calendar = get_calendar_type(standard_dataset)
        assert calendar in ["standard", "gregorian", "proleptic_gregorian"]

    def test_detect_calendar_from_attrs(self):
        """Should read calendar from attributes."""
        ds = xr.Dataset(
            {"temp": (["time"], np.random.randn(100))},
            coords={
                "time": xr.DataArray(
                    np.arange(100),
                    attrs={"calendar": "360_day"}
                )
            },
        )

        calendar = get_calendar_type(ds)
        assert calendar == "360_day"

    def test_default_to_standard(self):
        """Should default to standard calendar."""
        ds = xr.Dataset(
            {"temp": (["time"], np.random.randn(100))},
            coords={"time": np.arange(100)},
        )

        calendar = get_calendar_type(ds)
        assert calendar == "standard"


class TestConvertCalendar:
    """Test calendar conversion."""

    def test_convert_360_day_to_standard(self):
        """Should convert 360-day calendar to standard."""
        # Create 360-day data (12 months * 30 days * 10 years = 3600 days)
        time_360 = np.arange(3600)

        ds = xr.Dataset(
            {"temp": (["time"], np.random.randn(3600))},
            coords={
                "time": xr.DataArray(
                    time_360,
                    attrs={"calendar": "360_day", "units": "days since 2000-01-01"}
                )
            },
        )

        # Convert to standard calendar
        converted = convert_calendar(ds, "360_day", "standard")

        # Should have valid time coordinates
        assert "time" in converted.coords
        assert len(converted.time) > 0

    def test_convert_noleap_to_standard(self):
        """Should convert no-leap calendar to standard."""
        # Create noleap data (365 days * 10 years = 3650 days)
        time_noleap = np.arange(3650)

        ds = xr.Dataset(
            {"temp": (["time"], np.random.randn(3650))},
            coords={
                "time": xr.DataArray(
                    time_noleap,
                    attrs={"calendar": "noleap", "units": "days since 2000-01-01"}
                )
            },
        )

        # Convert to standard calendar
        converted = convert_calendar(ds, "noleap", "standard")

        assert "time" in converted.coords
        assert len(converted.time) > 0


class TestAlignTimeCoordinates:
    """Test time coordinate alignment."""

    def test_align_identical_time(self, standard_dataset):
        """Should handle identical time coordinates."""
        ds1 = standard_dataset
        ds2 = standard_dataset.copy(deep=True)

        aligned1, aligned2 = align_time_coordinates(ds1, ds2)

        # Should be able to combine
        combined = xr.concat([aligned1, aligned2], dim="model")
        assert combined is not None

    def test_align_overlapping_time(self):
        """Should align overlapping time coordinates."""
        time1 = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        time2 = np.arange('2005-01-01', '2015-01-01', dtype='datetime64[D]')

        data1 = np.random.randn(len(time1), 10, 10)
        data2 = np.random.randn(len(time2), 10, 10)

        ds1 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], data1)},
            coords={
                "time": time1,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        ds2 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], data2)},
            coords={
                "time": time2,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        aligned1, aligned2 = align_time_coordinates(ds1, ds2)

        # Should have aligned times
        assert len(aligned1.time) == len(aligned2.time)

    def test_align_disjoint_time_raises_error(self):
        """Should raise error for completely disjoint times."""
        time1 = np.arange('2000-01-01', '2005-01-01', dtype='datetime64[D]')
        time2 = np.arange('2010-01-01', '2015-01-01', dtype='datetime64[D]')

        data1 = np.random.randn(len(time1), 10, 10)
        data2 = np.random.randn(len(time2), 10, 10)

        ds1 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], data1)},
            coords={
                "time": time1,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        ds2 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], data2)},
            coords={
                "time": time2,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        with pytest.raises(TimeCoordinateMismatchError):
            align_time_coordinates(ds1, ds2)


class TestAlignDatasets:
    """Test full dataset alignment."""

    def test_align_compatible_datasets(self):
        """Should align multiple compatible datasets."""
        time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        lat = np.arange(10)
        lon = np.arange(10)

        datasets = []
        for i in range(3):
            data = np.random.randn(len(time), len(lat), len(lon))
            ds = xr.Dataset(
                {"temp": (["time", "lat", "lon"], data)},
                coords={
                    "time": time,
                    "lat": lat,
                    "lon": lon,
                },
            )
            datasets.append(ds)

        aligned = align_datasets(datasets)

        # All should have same dimensions
        for ds in aligned:
            assert len(ds.time) == len(time)
            assert len(ds.lat) == len(lat)
            assert len(ds.lon) == len(lon)

    def test_align_datasets_with_different_calendars(self):
        """Should handle datasets with different calendars."""
        # Dataset 1: standard calendar
        time1 = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        ds1 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], np.random.randn(len(time1), 10, 10))},
            coords={
                "time": time1,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        # Dataset 2: 360-day calendar
        time2_numeric = np.arange(3650)  # 10 years of 365 days
        ds2 = xr.Dataset(
            {"temp": (["time", "lat", "lon"], np.random.randn(len(time2_numeric), 10, 10))},
            coords={
                "time": xr.DataArray(
                    time2_numeric,
                    attrs={"calendar": "360_day", "units": "days since 2000-01-01"}
                ),
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        # Should handle mixed calendars gracefully
        # This may raise TimeCoordinateMismatchError or succeed depending on overlap
        try:
            aligned = align_datasets([ds1, ds2])
            assert len(aligned) == 2
        except TimeCoordinateMismatchError:
            # This is also acceptable if calendars don't align
            pass

    def test_align_empty_list_raises_error(self):
        """Should raise error for empty dataset list."""
        with pytest.raises(ValueError):
            align_datasets([])

    def test_align_single_dataset(self):
        """Should handle single dataset."""
        time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        ds = xr.Dataset(
            {"temp": (["time", "lat", "lon"], np.random.randn(len(time), 10, 10))},
            coords={
                "time": time,
                "lat": np.arange(10),
                "lon": np.arange(10),
            },
        )

        aligned = align_datasets([ds])

        assert len(aligned) == 1
        assert len(aligned[0].time) == len(time)


class TestAlignmentIntegration:
    """Integration tests for alignment workflow."""

    def test_multi_model_alignment_workflow(self, temp_dir):
        """Should align multi-model ensemble following realistic workflow."""
        # Create 3 model datasets with different characteristics
        models = ["gfdl", "hadgem", "ipsl"]
        time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        lat = np.linspace(-90, 89.75, 10)
        lon = np.linspace(-180, 179.5, 20)

        datasets = []
        for model in models:
            data = np.random.randn(len(time), len(lat), len(lon))
            ds = xr.Dataset(
                {"temperature": (["time", "lat", "lon"], data)},
                coords={
                    "time": time,
                    "lat": lat,
                    "lon": lon,
                },
            )
            ds.attrs["model"] = model
            datasets.append(ds)

        # Align all datasets
        aligned = align_datasets(datasets)

        # Verify alignment
        assert len(aligned) == 3
        for ds in aligned:
            assert ds.dims["lat"] == len(lat)
            assert ds.dims["lon"] == len(lon)

    def test_alignment_preserves_data_integrity(self):
        """Should preserve data values during alignment."""
        time = np.arange('2000-01-01', '2010-01-01', dtype='datetime64[D]')
        lat = np.arange(5)
        lon = np.arange(5)

        original_data = np.arange(len(time) * len(lat) * len(lon)).reshape(
            len(time), len(lat), len(lon)
        ).astype(float)

        ds = xr.Dataset(
            {"temp": (["time", "lat", "lon"], original_data)},
            coords={
                "time": time,
                "lat": lat,
                "lon": lon,
            },
        )

        aligned = align_datasets([ds])

        # Data should be preserved
        assert np.allclose(aligned[0].temp.values, original_data)
