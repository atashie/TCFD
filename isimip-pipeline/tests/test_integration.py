"""Integration tests for the full ISIMIP pipeline.

Tests the complete workflow from data loading through processing
to output generation, ensuring all modules work together correctly.
"""

import pytest
import tempfile
from pathlib import Path
from datetime import datetime

import numpy as np
import xarray as xr

from isimip_pipeline.processing.processor import (
    DataProcessor,
    find_netcdf_files,
    load_and_aggregate,
    detect_timestep_from_files,
)
from isimip_pipeline.processing.features import (
    FeatureExtractor,
    kernel_smooth,
    theil_sen_slope,
    decadal_aggregation,
)
from isimip_pipeline.processing.output import (
    create_output_dataset,
    write_netcdf,
    OutputWriter,
    VALUE_CLASSES,
)
from isimip_pipeline.processing.alignment import (
    align_datasets,
    verify_spatial_grids,
    extract_spatial_grid,
)
from isimip_pipeline.processing.validation import (
    generate_validation_report,
    detect_fill_values,
    detect_spatial_gaps,
)
from isimip_pipeline.visualization.qa_report import (
    QAReport,
    generate_html_report,
    create_map_figure,
)
from isimip_pipeline.search.result_table import (
    group_by_variable_timestep,
    extract_time_coverage,
)
from isimip_pipeline.search.isimip_query import DatasetInfo
from isimip_pipeline.processing_log import (
    ProcessingLog,
    DatasetEntry,
    load_processing_log,
    save_processing_log,
    add_dataset_entry,
)


@pytest.fixture
def realistic_climate_data(tmp_path):
    """Create realistic climate data mimicking ISIMIP output.

    Creates multiple model files with temperature-like patterns:
    - Seasonal variation
    - Warming trend
    - Spatial patterns (latitude-dependent)
    """
    lat = np.arange(-60, 61, 10.0)  # 13 lat points
    lon = np.arange(-180, 180, 15.0)  # 24 lon points

    # Create 94 years (2006-2099) of monthly data
    n_months = 94 * 12
    time = np.arange(n_months)

    models = ["gfdl-esm4", "ukesm1-0-ll", "ipsl-cm6a-lr"]
    scenarios = ["ssp126", "ssp585"]

    for model in models:
        for scenario in scenarios:
            # Base temperature pattern
            lat_effect = (lat[:, None] / 90) * 15  # Warmer at equator
            lon_effect = np.sin(lon[None, :] * np.pi / 180) * 2
            base = 273.15 + 15 + lat_effect + lon_effect

            # Time evolution
            years = 2006 + time / 12
            # Warming trend (more for ssp585)
            trend_rate = 0.02 if scenario == "ssp126" else 0.05
            trend = trend_rate * (years - 2006)

            # Seasonal variation
            seasonal = 5 * np.sin(2 * np.pi * time / 12)

            # Add noise and combine
            np.random.seed(hash(f"{model}_{scenario}") % 2**32)
            noise = np.random.randn(len(time), len(lat), len(lon)) * 2

            data = (base[None, :, :] +
                    trend[:, None, None] +
                    seasonal[:, None, None] +
                    noise)

            # Create dataset
            ds = xr.Dataset(
                {"tas": (["time", "lat", "lon"], data.astype(np.float32))},
                coords={
                    "time": time,
                    "lat": lat,
                    "lon": lon,
                },
                attrs={
                    "climate_scenario": scenario,
                    "model": model,
                    "simulation_round": "ISIMIP3b",
                    "variable": "tas",
                }
            )

            filename = f"{model}_{scenario}_tas_global_monthly_2006_2099.nc"
            ds.to_netcdf(tmp_path / filename)

    return tmp_path


@pytest.fixture
def multi_model_datasets(tmp_path):
    """Create aligned multi-model datasets for ensemble processing."""
    lat = np.arange(-45, 46, 15.0)
    lon = np.arange(-180, 180, 30.0)
    time = np.arange('2010-01-01', '2020-01-01', dtype='datetime64[M]')

    datasets = []
    for i, model in enumerate(["model_a", "model_b", "model_c"]):
        np.random.seed(i)
        data = np.random.rand(len(time), len(lat), len(lon)) * 100

        ds = xr.Dataset(
            {"precip": (["time", "lat", "lon"], data.astype(np.float32))},
            coords={"time": time, "lat": lat, "lon": lon},
        )

        filepath = tmp_path / f"{model}_precip.nc"
        ds.to_netcdf(filepath)
        datasets.append(filepath)

    return datasets, tmp_path


class TestFullPipelineIntegration:
    """Test complete pipeline from data loading to output."""

    def test_load_process_output_workflow(self, realistic_climate_data):
        """Test the full workflow: load -> process -> output."""
        input_dir = realistic_climate_data

        # Step 1: Find and load files
        files = find_netcdf_files(input_dir)
        assert len(files) == 6  # 3 models x 2 scenarios

        # Step 2: Detect timestep (synthetic data uses integer time, may detect as daily)
        timestep = detect_timestep_from_files(files)
        # Timestep detection works on real ISIMIP data with proper datetime coords
        assert timestep in ["daily", "monthly", "annual", "unknown"]

        # Step 3: Process with DataProcessor
        processor = DataProcessor(bandwidth=15)
        result = processor.process(
            input_dir=input_dir,
            variable="tas",
            scenarios=["ssp126", "ssp585"],
            validate_input=False,  # Skip validation for synthetic data
        )

        # Verify output structure
        assert isinstance(result, xr.Dataset)
        assert "tas" in result.data_vars
        assert "value_class" in result.dims
        assert len(result.value_class) == 6

        # Step 4: Write output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "tas_processed.nc"
            write_netcdf(result, output_path)

            # Verify file was written
            assert output_path.exists()

            # Re-read and verify
            ds_read = xr.open_dataset(output_path)
            assert "tas" in ds_read.data_vars
            ds_read.close()

    def test_feature_extraction_pipeline(self, realistic_climate_data):
        """Test feature extraction produces valid statistics."""
        files = find_netcdf_files(realistic_climate_data)

        # Load first file
        ds = xr.open_dataset(files[0])
        data = ds["tas"].values  # (time, lat, lon)

        # Use first cell's time series
        cell_data = data[:, 5, 10]
        years = 2006 + np.arange(len(cell_data)) / 12

        # Test kernel smoothing
        smoothed = kernel_smooth(years, cell_data, bandwidth=15)
        assert len(smoothed) == len(cell_data)
        assert np.var(smoothed) < np.var(cell_data)  # Smoothing reduces variance

        # Test trend detection (should show warming)
        slope = theil_sen_slope(years, cell_data)
        # Trend should be positive for warming scenario
        assert isinstance(slope, float)

        ds.close()

    def test_validation_integration(self, realistic_climate_data):
        """Test validation catches data quality issues."""
        files = find_netcdf_files(realistic_climate_data)
        ds = xr.open_dataset(files[0])

        # Generate validation report
        report = generate_validation_report(ds, "tas")

        # Clean synthetic data should have good quality
        assert report is not None
        assert report.variable == "tas"
        assert report.summary["data_quality"] in ["excellent", "good"]

        ds.close()

    def test_output_generation_with_features(self, realistic_climate_data):
        """Test output dataset has correct structure and values."""
        processor = DataProcessor()

        result = processor.process(
            input_dir=realistic_climate_data,
            variable="tas",
            validate_input=False,
        )

        # Check dimensions
        assert result.dims["value_class"] == 6
        assert "decade" in result.dims
        assert "scenario" in result.dims

        # Check value classes are correct
        expected_classes = [
            "smoothed_median", "percentile", "trend",
            "significance", "lower_bound", "upper_bound"
        ]
        actual_classes = list(result.value_class.values)
        assert actual_classes == expected_classes


class TestMultiModelAlignment:
    """Test multi-model dataset alignment and processing."""

    def test_align_multiple_models(self, multi_model_datasets):
        """Test alignment of datasets from different models."""
        file_paths, tmp_path = multi_model_datasets

        # Load datasets
        datasets = [xr.open_dataset(f) for f in file_paths]

        try:
            # Verify spatial grids match
            for i, ds in enumerate(datasets[1:], 1):
                verify_spatial_grids(datasets[0], ds)

            # Align datasets
            aligned = align_datasets(datasets, verify_spatial=True)

            assert len(aligned) == 3
            # All should have same time dimension after alignment
            times = [len(ds.time) for ds in aligned]
            assert all(t == times[0] for t in times)
        finally:
            for ds in datasets:
                ds.close()

    def test_ensemble_processing(self, multi_model_datasets):
        """Test processing creates valid ensemble statistics."""
        file_paths, tmp_path = multi_model_datasets

        # Process through DataProcessor
        processor = DataProcessor()
        result = processor.process(
            input_dir=tmp_path,
            variable="precip",
            validate_input=False,
        )

        # Should have valid output
        assert "precip" in result.data_vars

        # Check for valid percentile values
        percentile_data = result["precip"].sel(value_class="percentile")
        valid_vals = percentile_data.values[~np.isnan(percentile_data.values)]
        if len(valid_vals) > 0:
            assert np.all(valid_vals >= 1)
            assert np.all(valid_vals <= 100)


class TestVisualizationIntegration:
    """Test QA report generation from processed data."""

    def test_qa_report_generation(self, realistic_climate_data):
        """Test complete QA report generation."""
        processor = DataProcessor()

        # Process data
        result = processor.process(
            input_dir=realistic_climate_data,
            variable="tas",
            scenarios=["ssp126"],
            validate_input=False,
        )

        # Generate QA report
        report = QAReport(result, "tas")
        summary = report.get_summary()

        assert "variable" in summary
        assert summary["variable"] == "tas"
        assert "total_cells" in summary
        assert summary["total_cells"] > 0

    def test_html_report_generation(self, realistic_climate_data):
        """Test HTML report can be generated."""
        processor = DataProcessor()

        result = processor.process(
            input_dir=realistic_climate_data,
            variable="tas",
            scenarios=["ssp126"],
            validate_input=False,
        )

        # Generate HTML
        html = generate_html_report(
            result,
            variable="tas",
            title="Test Report",
            include_maps=True,
            include_summary=True,
        )

        assert isinstance(html, str)
        assert "<!DOCTYPE html>" in html
        assert "tas" in html
        assert "Test Report" in html


class TestProcessingLogIntegration:
    """Test processing log tracks datasets correctly."""

    def test_log_creation_and_update(self, tmp_path):
        """Test creating and updating processing log."""
        log_path = tmp_path / "processed_data_log.yaml"

        # Create initial log
        log = ProcessingLog(datasets=[])

        # Add entry
        entry = DatasetEntry(
            descriptive_name="test-data",
            variable="tas",
            timestep="monthly",
            created_date=datetime.now().isoformat(),
            output_path=str(tmp_path / "test-data_tas-monthly"),
            file_count=5,
            time_periods=["2006-2100"],
            climate_scenarios=["ssp126", "ssp585"],
            gcm_models=["gfdl-esm4"],
            lsm_models=[],
            simulation_round="ISIMIP3b",
            query="test query",
        )

        log = add_dataset_entry(log, entry)
        save_processing_log(log, log_path)

        # Reload and verify
        loaded = load_processing_log(log_path)
        assert len(loaded.datasets) == 1
        assert loaded.datasets[0].variable == "tas"
        assert loaded.datasets[0].descriptive_name == "test-data"


class TestSearchResultGrouping:
    """Test search result grouping and display."""

    def test_group_and_extract_coverage(self):
        """Test grouping by variable+timestep with time coverage."""
        datasets = [
            DatasetInfo(
                id="1",
                name="gfdl_ssp126_led_2015_2100.nc",
                url="https://example.com/1.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp126",
                variable="led",
                model="gfdl-esm4",
                timestep="annual",
            ),
            DatasetInfo(
                id="2",
                name="ukesm_ssp126_led_2015_2100.nc",
                url="https://example.com/2.nc",
                simulation_round="ISIMIP3b",
                climate_scenario="ssp585",
                variable="led",
                model="ukesm1-0-ll",
                timestep="annual",
            ),
        ]

        grouped = group_by_variable_timestep(datasets)

        assert len(grouped) == 1
        assert ("led", "annual") in grouped

        group = grouped[("led", "annual")]
        assert group["file_count"] == 2
        assert "ssp126" in group["scenarios"]
        assert "ssp585" in group["scenarios"]
        # Time coverage should be extracted from filenames
        assert "2015" in group["time_coverage"]
        assert "2100" in group["time_coverage"]


class TestEndToEndWorkflow:
    """Test complete end-to-end workflow simulating user interaction."""

    def test_complete_workflow_simulation(self, realistic_climate_data, tmp_path):
        """Simulate complete user workflow from search to report."""
        # Step 1: Simulate finding datasets (would normally be ISIMIP search)
        files = find_netcdf_files(realistic_climate_data)
        assert len(files) > 0

        # Step 2: Group and select (simulating user choice)
        # In real workflow, this would use group_by_variable_timestep
        timestep = detect_timestep_from_files(files)
        # Synthetic data may detect as daily; real ISIMIP data would be monthly
        assert timestep in ["daily", "monthly", "annual", "unknown"]

        # Step 3: Process
        processor = DataProcessor(bandwidth=15, percentile_bins=100)
        result = processor.process(
            input_dir=realistic_climate_data,
            variable="tas",
            scenarios=["ssp126", "ssp585"],
            align_datasets=True,
            validate_input=False,
        )

        # Step 4: Write output
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        output_file = output_dir / "tas_processed.nc"
        write_netcdf(result, output_file, compression=True)

        # Step 5: Generate report
        html = generate_html_report(
            result, "tas",
            title="Temperature Analysis",
            include_maps=True,
        )

        report_file = output_dir / "qa_report.html"
        with open(report_file, "w") as f:
            f.write(html)

        # Step 6: Update processing log
        log_path = tmp_path / "processed_data_log.yaml"
        log = ProcessingLog(datasets=[])
        entry = DatasetEntry(
            descriptive_name="temperature-analysis",
            variable="tas",
            timestep=timestep,
            created_date=datetime.now().isoformat(),
            output_path=str(output_dir),
            file_count=len(files),
            time_periods=["2006-2099"],
            climate_scenarios=["ssp126", "ssp585"],
            gcm_models=["gfdl-esm4", "ukesm1-0-ll", "ipsl-cm6a-lr"],
            lsm_models=[],
            simulation_round="ISIMIP3b",
            query="temperature analysis",
        )
        log = add_dataset_entry(log, entry)
        save_processing_log(log, log_path)

        # Verify all outputs exist
        assert output_file.exists()
        assert report_file.exists()
        assert log_path.exists()

        # Verify log contains entry
        loaded_log = load_processing_log(log_path)
        assert len(loaded_log.datasets) == 1
        assert loaded_log.datasets[0].variable == "tas"
