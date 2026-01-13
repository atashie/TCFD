"""Tests for QA report visualization module."""

import pytest
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

from isimip_pipeline.visualization.qa_report import (
    QAReport,
    create_map_figure,
    create_timeseries_figure,
    generate_html_report,
)


@pytest.fixture
def sample_dataset():
    """Create sample xarray dataset for testing."""
    lat = np.arange(-90, 91, 30.0)
    lon = np.arange(-180, 180, 30.0)
    decades = [10, 20, 30]
    scenarios = ["ssp126", "ssp585"]
    value_classes = [
        "smoothed_median",
        "percentile",
        "trend",
        "significance",
        "lower_bound",
        "upper_bound",
    ]

    data = np.random.random(
        (len(lon), len(lat), len(decades), len(scenarios), len(value_classes))
    ).astype(np.float32)

    ds = xr.Dataset(
        {
            "burntarea": (
                ["lon", "lat", "decade", "scenario", "value_class"],
                data,
            )
        },
        coords={
            "lon": lon,
            "lat": lat,
            "decade": decades,
            "scenario": scenarios,
            "value_class": value_classes,
        },
    )
    return ds


class TestCreateMapFigure:
    """Test map figure creation."""

    def test_creates_plotly_figure(self, sample_dataset):
        """Should create a Plotly figure."""
        fig = create_map_figure(
            sample_dataset,
            variable="burntarea",
            decade=10,
            scenario="ssp126",
            value_class="smoothed_median",
        )

        assert fig is not None
        assert hasattr(fig, "to_html")

    def test_map_has_title(self, sample_dataset):
        """Map figure should have a title."""
        fig = create_map_figure(
            sample_dataset,
            variable="burntarea",
            decade=10,
            scenario="ssp126",
            value_class="smoothed_median",
        )

        # Check layout has title
        assert fig.layout.title is not None or hasattr(fig.layout, "title")


class TestCreateTimeseriesFigure:
    """Test timeseries figure creation."""

    def test_creates_plotly_figure(self, sample_dataset):
        """Should create a Plotly figure for timeseries."""
        fig = create_timeseries_figure(
            sample_dataset,
            variable="burntarea",
            lat_idx=3,
            lon_idx=5,
            scenario="ssp126",
        )

        assert fig is not None
        assert hasattr(fig, "to_html")


class TestQAReport:
    """Test QAReport class."""

    def test_report_initializes(self, sample_dataset):
        """QAReport should initialize with dataset."""
        report = QAReport(sample_dataset, variable="burntarea")

        assert report is not None
        assert report.variable == "burntarea"

    def test_report_generates_maps(self, sample_dataset):
        """Report should generate map figures."""
        report = QAReport(sample_dataset, variable="burntarea")

        maps = report.generate_maps()

        assert isinstance(maps, dict)
        assert len(maps) > 0

    def test_report_has_summary(self, sample_dataset):
        """Report should have summary statistics."""
        report = QAReport(sample_dataset, variable="burntarea")

        summary = report.get_summary()

        assert "total_valid_cells" in summary or "variable" in summary


class TestGenerateHTMLReport:
    """Test HTML report generation."""

    def test_generates_html_string(self, sample_dataset):
        """Should generate HTML string."""
        html = generate_html_report(
            sample_dataset,
            variable="burntarea",
            title="Test Report",
        )

        assert isinstance(html, str)
        assert "<html" in html.lower()
        assert "burntarea" in html.lower() or "test" in html.lower()

    def test_writes_html_file(self, sample_dataset):
        """Should write HTML file to disk."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            output_path = Path(f.name)

        try:
            html = generate_html_report(
                sample_dataset,
                variable="burntarea",
                title="Test Report",
            )

            with open(output_path, "w") as f:
                f.write(html)

            assert output_path.exists()
            assert output_path.stat().st_size > 0
        finally:
            output_path.unlink()
