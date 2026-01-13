"""Tests for download functionality."""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, AsyncMock
import asyncio

from isimip_pipeline.download.downloader import (
    Downloader,
    DownloadResult,
    DownloadStatus,
    download_files,
)
from isimip_pipeline.search.isimip_query import DatasetInfo


class TestDownloadResult:
    """Test DownloadResult data structure."""

    def test_download_result_success(self):
        """DownloadResult should represent successful download."""
        result = DownloadResult(
            url="https://example.com/test.nc",
            path=Path("/data/test.nc"),
            status=DownloadStatus.SUCCESS,
            size=1024,
        )
        assert result.status == DownloadStatus.SUCCESS
        assert result.size == 1024
        assert result.error is None

    def test_download_result_failure(self):
        """DownloadResult should represent failed download."""
        result = DownloadResult(
            url="https://example.com/test.nc",
            path=None,
            status=DownloadStatus.FAILED,
            error="Connection timeout",
        )
        assert result.status == DownloadStatus.FAILED
        assert result.error == "Connection timeout"

    def test_download_result_skipped(self):
        """DownloadResult should represent skipped download."""
        result = DownloadResult(
            url="https://example.com/test.nc",
            path=Path("/data/test.nc"),
            status=DownloadStatus.SKIPPED,
        )
        assert result.status == DownloadStatus.SKIPPED


class TestDownloader:
    """Test Downloader class."""

    def test_downloader_initializes(self):
        """Downloader should initialize with default settings."""
        downloader = Downloader()
        assert downloader is not None
        assert downloader.max_concurrent > 0

    def test_downloader_accepts_output_dir(self):
        """Downloader should accept custom output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(output_dir=Path(tmpdir))
            assert downloader.output_dir == Path(tmpdir)

    def test_downloader_accepts_max_concurrent(self):
        """Downloader should accept custom concurrency limit."""
        downloader = Downloader(max_concurrent=2)
        assert downloader.max_concurrent == 2


class TestDownloaderSkipExisting:
    """Test that downloader skips existing files."""

    def test_skip_existing_file(self):
        """Downloader should skip files that already exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create existing file
            existing = Path(tmpdir) / "existing.nc"
            existing.write_text("existing data")

            downloader = Downloader(output_dir=Path(tmpdir))

            # Check if file should be skipped
            should_skip = downloader.should_skip(
                "https://example.com/existing.nc", existing
            )
            assert should_skip is True

    def test_dont_skip_nonexistent_file(self):
        """Downloader should not skip files that don't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(output_dir=Path(tmpdir))
            missing = Path(tmpdir) / "missing.nc"

            should_skip = downloader.should_skip(
                "https://example.com/missing.nc", missing
            )
            assert should_skip is False


class TestDownloaderFromDatasetInfo:
    """Test downloading from DatasetInfo objects."""

    def test_extract_urls_from_dataset_info(self):
        """Downloader should extract URLs from DatasetInfo list."""
        datasets = [
            DatasetInfo(
                id="1",
                name="file1.nc",
                url="https://files.isimip.org/file1.nc",
                simulation_round="ISIMIP3b",
            ),
            DatasetInfo(
                id="2",
                name="file2.nc",
                url="https://files.isimip.org/file2.nc",
                simulation_round="ISIMIP3b",
            ),
        ]

        downloader = Downloader()
        urls = downloader.get_urls_from_datasets(datasets)

        assert len(urls) == 2
        assert "https://files.isimip.org/file1.nc" in urls
        assert "https://files.isimip.org/file2.nc" in urls


class TestDownloaderTempDirectory:
    """Test temp directory functionality."""

    def test_creates_temp_directory(self):
        """Downloader should create temp directory when use_temp=True."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(use_temp=True, temp_base=Path(tmpdir))
            assert downloader.temp_dir is not None
            assert downloader.temp_dir.exists()
            assert downloader.temp_dir.name.startswith("isimip_")
            # Cleanup
            downloader.cleanup()

    def test_temp_directory_in_specified_base(self):
        """Temp directory should be created in specified base."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(use_temp=True, temp_base=Path(tmpdir))
            assert downloader.temp_dir.parent == Path(tmpdir)
            downloader.cleanup()

    def test_cleanup_removes_temp_directory(self):
        """Cleanup should remove temp directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(use_temp=True, temp_base=Path(tmpdir))
            temp_path = downloader.temp_dir
            assert temp_path.exists()

            result = downloader.cleanup()
            assert result is True
            assert not temp_path.exists()

    def test_context_manager_cleans_up(self):
        """Using as context manager should cleanup on exit."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with Downloader(use_temp=True, temp_base=Path(tmpdir)) as downloader:
                temp_path = downloader.temp_dir
                assert temp_path.exists()
            # After context exit, temp should be gone
            assert not temp_path.exists()

    def test_no_temp_when_use_temp_false(self):
        """No temp directory when use_temp=False."""
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = Downloader(output_dir=Path(tmpdir), use_temp=False)
            assert downloader.temp_dir is None
            assert downloader.output_dir == Path(tmpdir)


class TestDownloadFilesFunction:
    """Test the convenience download_files function."""

    @patch("isimip_pipeline.download.downloader.Downloader")
    def test_download_files_creates_downloader(self, mock_downloader_class):
        """download_files should create a Downloader instance."""
        mock_downloader = Mock()
        mock_downloader.download_all = AsyncMock(return_value=[])
        mock_downloader_class.return_value = mock_downloader

        with tempfile.TemporaryDirectory() as tmpdir:
            results = download_files(
                urls=["https://example.com/test.nc"],
                output_dir=Path(tmpdir),
            )

        mock_downloader_class.assert_called_once()
