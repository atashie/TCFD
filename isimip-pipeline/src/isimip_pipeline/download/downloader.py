"""Download functionality for ISIMIP datasets."""

import asyncio
import shutil
import tempfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Callable
from urllib.parse import urlparse

import httpx
from rich.progress import Progress, TaskID

from isimip_pipeline.search.isimip_query import DatasetInfo

# Default temp directory base for downloads
DEFAULT_TEMP_BASE = Path.home() / "Downloads"


class DownloadStatus(Enum):
    """Status of a download operation."""

    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"
    IN_PROGRESS = "in_progress"


@dataclass
class DownloadResult:
    """Result of a download operation.

    Attributes:
        url: The URL that was downloaded.
        path: Local path where file was saved (None if failed).
        status: Status of the download.
        size: Size of downloaded file in bytes.
        error: Error message if download failed.
    """

    url: str
    path: Optional[Path]
    status: DownloadStatus
    size: Optional[int] = None
    error: Optional[str] = None


class Downloader:
    """Handles downloading files from ISIMIP repository.

    Supports concurrent downloads, progress tracking, resume, and temp directory management.
    """

    def __init__(
        self,
        output_dir: Optional[Path] = None,
        max_concurrent: int = 4,
        chunk_size: int = 8192,
        timeout: int = 600,
        use_temp: bool = False,
        temp_base: Optional[Path] = None,
    ):
        """Initialize the downloader.

        Args:
            output_dir: Directory to save downloaded files.
            max_concurrent: Maximum concurrent downloads.
            chunk_size: Size of chunks for streaming download.
            timeout: Request timeout in seconds.
            use_temp: If True, create a temp directory for downloads.
            temp_base: Base directory for temp folder (default: ~/Downloads).
        """
        self.max_concurrent = max_concurrent
        self.chunk_size = chunk_size
        self.timeout = timeout
        self.use_temp = use_temp
        self._temp_dir: Optional[Path] = None

        if use_temp:
            # Create temp directory in specified base or default
            base = temp_base or DEFAULT_TEMP_BASE
            base.mkdir(parents=True, exist_ok=True)
            self._temp_dir = Path(tempfile.mkdtemp(prefix="isimip_", dir=base))
            self.output_dir = self._temp_dir
        else:
            self.output_dir = output_dir or Path("./data/raw")

    @property
    def temp_dir(self) -> Optional[Path]:
        """Get the temp directory path if using temp storage."""
        return self._temp_dir

    def cleanup(self) -> bool:
        """Remove temp directory and all downloaded files.

        Returns:
            True if cleanup was successful, False otherwise.
        """
        if self._temp_dir and self._temp_dir.exists():
            try:
                shutil.rmtree(self._temp_dir)
                self._temp_dir = None
                return True
            except Exception:
                return False
        return True

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - cleanup temp files."""
        if self.use_temp:
            self.cleanup()
        return False

    def get_filename_from_url(self, url: str) -> str:
        """Extract filename from URL.

        Args:
            url: URL to extract filename from.

        Returns:
            Filename portion of the URL.
        """
        parsed = urlparse(url)
        return Path(parsed.path).name

    def get_output_path(self, url: str) -> Path:
        """Get the output path for a URL.

        Args:
            url: URL to get output path for.

        Returns:
            Full path where file will be saved.
        """
        filename = self.get_filename_from_url(url)
        return self.output_dir / filename

    def should_skip(self, url: str, path: Path) -> bool:
        """Check if download should be skipped.

        Args:
            url: URL to download.
            path: Local path where file would be saved.

        Returns:
            True if file exists and should be skipped.
        """
        return path.exists()

    def get_urls_from_datasets(self, datasets: List[DatasetInfo]) -> List[str]:
        """Extract URLs from DatasetInfo objects.

        Args:
            datasets: List of DatasetInfo objects.

        Returns:
            List of download URLs.
        """
        return [ds.url for ds in datasets if ds.url]

    async def download_one(
        self,
        client: httpx.AsyncClient,
        url: str,
        progress: Optional[Progress] = None,
        task_id: Optional[TaskID] = None,
    ) -> DownloadResult:
        """Download a single file.

        Args:
            client: HTTP client to use.
            url: URL to download.
            progress: Optional progress bar.
            task_id: Optional task ID for progress bar.

        Returns:
            DownloadResult with status and details.
        """
        output_path = self.get_output_path(url)

        # Check if should skip
        if self.should_skip(url, output_path):
            return DownloadResult(
                url=url,
                path=output_path,
                status=DownloadStatus.SKIPPED,
                size=output_path.stat().st_size,
            )

        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            async with client.stream("GET", url) as response:
                response.raise_for_status()

                total = int(response.headers.get("content-length", 0))

                if progress and task_id is not None:
                    progress.update(task_id, total=total)

                with open(output_path, "wb") as f:
                    downloaded = 0
                    async for chunk in response.aiter_bytes(self.chunk_size):
                        f.write(chunk)
                        downloaded += len(chunk)
                        if progress and task_id is not None:
                            progress.update(task_id, completed=downloaded)

            return DownloadResult(
                url=url,
                path=output_path,
                status=DownloadStatus.SUCCESS,
                size=downloaded,
            )

        except httpx.HTTPStatusError as e:
            return DownloadResult(
                url=url,
                path=None,
                status=DownloadStatus.FAILED,
                error=f"HTTP {e.response.status_code}: {e.response.reason_phrase}",
            )
        except httpx.RequestError as e:
            return DownloadResult(
                url=url,
                path=None,
                status=DownloadStatus.FAILED,
                error=str(e),
            )

    async def download_all(
        self,
        urls: List[str],
        progress_callback: Optional[Callable[[str, int, int], None]] = None,
    ) -> List[DownloadResult]:
        """Download multiple files concurrently.

        Args:
            urls: List of URLs to download.
            progress_callback: Optional callback for progress updates.

        Returns:
            List of DownloadResult objects.
        """
        semaphore = asyncio.Semaphore(self.max_concurrent)

        async def limited_download(client: httpx.AsyncClient, url: str) -> DownloadResult:
            async with semaphore:
                return await self.download_one(client, url)

        async with httpx.AsyncClient(timeout=self.timeout) as client:
            tasks = [limited_download(client, url) for url in urls]
            results = await asyncio.gather(*tasks)

        return list(results)

    def download_sync(self, urls: List[str]) -> List[DownloadResult]:
        """Synchronous wrapper for download_all.

        Args:
            urls: List of URLs to download.

        Returns:
            List of DownloadResult objects.
        """
        return asyncio.run(self.download_all(urls))


def download_files(
    urls: List[str],
    output_dir: Optional[Path] = None,
    max_concurrent: int = 4,
    use_temp: bool = False,
    temp_base: Optional[Path] = None,
) -> List[DownloadResult]:
    """Convenience function to download files.

    Args:
        urls: List of URLs to download.
        output_dir: Directory to save files.
        max_concurrent: Maximum concurrent downloads.
        use_temp: If True, use a temp directory.
        temp_base: Base directory for temp folder.

    Returns:
        List of DownloadResult objects.
    """
    downloader = Downloader(
        output_dir=output_dir,
        max_concurrent=max_concurrent,
        use_temp=use_temp,
        temp_base=temp_base,
    )
    return downloader.download_sync(urls)


def create_temp_downloader(
    temp_base: Optional[Path] = None,
    max_concurrent: int = 4,
) -> Downloader:
    """Create a downloader that uses temp directory.

    The temp directory will be in ~/Downloads by default.
    Use as context manager for automatic cleanup:

        with create_temp_downloader() as downloader:
            results = downloader.download_sync(urls)
            # process files in downloader.output_dir
        # temp files automatically cleaned up

    Args:
        temp_base: Base directory for temp folder (default: ~/Downloads).
        max_concurrent: Maximum concurrent downloads.

    Returns:
        Downloader configured with temp directory.
    """
    return Downloader(
        use_temp=True,
        temp_base=temp_base,
        max_concurrent=max_concurrent,
    )
