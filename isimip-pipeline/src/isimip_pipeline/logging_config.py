"""Logging configuration for the ISIMIP pipeline.

Provides structured logging with configurable verbosity levels and output formats.
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

# Default log format
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
SIMPLE_FORMAT = "%(levelname)s: %(message)s"
DETAILED_FORMAT = (
    "%(asctime)s - %(name)s - %(levelname)s - "
    "[%(filename)s:%(lineno)d] - %(message)s"
)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_style: str = "default",
    quiet: bool = False,
) -> logging.Logger:
    """Configure logging for the ISIMIP pipeline.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        log_file: Optional path to write log file.
        format_style: Format style ('default', 'simple', 'detailed').
        quiet: If True, suppress console output (only log to file).

    Returns:
        Configured root logger for the pipeline.
    """
    # Select format
    format_map = {
        "default": DEFAULT_FORMAT,
        "simple": SIMPLE_FORMAT,
        "detailed": DETAILED_FORMAT,
    }
    log_format = format_map.get(format_style, DEFAULT_FORMAT)

    # Get or create pipeline logger
    logger = logging.getLogger("isimip_pipeline")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Clear existing handlers
    logger.handlers = []

    # Console handler (unless quiet mode)
    if not quiet:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter(SIMPLE_FORMAT))
        logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)  # Log everything to file
        file_handler.setFormatter(logging.Formatter(DETAILED_FORMAT))
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a specific module.

    Args:
        name: Module name (e.g., 'processing', 'download').

    Returns:
        Logger instance for the module.
    """
    return logging.getLogger(f"isimip_pipeline.{name}")


def log_processing_start(
    logger: logging.Logger,
    input_dir: Path,
    variable: str,
    n_files: int,
) -> None:
    """Log the start of a processing job.

    Args:
        logger: Logger instance.
        input_dir: Input directory path.
        variable: Variable being processed.
        n_files: Number of input files.
    """
    logger.info(
        f"Starting processing: variable={variable}, "
        f"files={n_files}, input={input_dir}"
    )


def log_processing_complete(
    logger: logging.Logger,
    output_path: Path,
    duration_seconds: float,
    success: bool = True,
) -> None:
    """Log completion of a processing job.

    Args:
        logger: Logger instance.
        output_path: Output file path.
        duration_seconds: Processing duration in seconds.
        success: Whether processing succeeded.
    """
    status = "SUCCESS" if success else "FAILED"
    logger.info(
        f"Processing {status}: output={output_path}, "
        f"duration={duration_seconds:.1f}s"
    )


def log_download_progress(
    logger: logging.Logger,
    url: str,
    status: str,
    bytes_downloaded: int = 0,
) -> None:
    """Log download progress.

    Args:
        logger: Logger instance.
        url: URL being downloaded.
        status: Download status (started, progress, complete, failed).
        bytes_downloaded: Number of bytes downloaded so far.
    """
    filename = url.split("/")[-1]
    if bytes_downloaded:
        mb = bytes_downloaded / (1024 * 1024)
        logger.debug(f"Download {status}: {filename} ({mb:.1f} MB)")
    else:
        logger.debug(f"Download {status}: {filename}")


def log_validation_result(
    logger: logging.Logger,
    variable: str,
    quality: str,
    issues: int,
) -> None:
    """Log validation results.

    Args:
        logger: Logger instance.
        variable: Variable that was validated.
        quality: Quality level (excellent, good, acceptable, poor, unusable).
        issues: Number of issues found.
    """
    level = logging.WARNING if quality in ("poor", "unusable") else logging.INFO
    logger.log(
        level,
        f"Validation: variable={variable}, quality={quality}, issues={issues}"
    )


def log_api_call(
    logger: logging.Logger,
    endpoint: str,
    params: dict,
    success: bool = True,
    response_time_ms: float = 0,
) -> None:
    """Log API calls for debugging.

    Args:
        logger: Logger instance.
        endpoint: API endpoint called.
        params: Request parameters.
        success: Whether call succeeded.
        response_time_ms: Response time in milliseconds.
    """
    status = "OK" if success else "FAILED"
    logger.debug(
        f"API {status}: {endpoint} ({response_time_ms:.0f}ms) params={params}"
    )


class LogContext:
    """Context manager for logging operation timing.

    Usage:
        with LogContext(logger, "Processing data"):
            # ... do processing ...
    """

    def __init__(self, logger: logging.Logger, operation: str):
        """Initialize log context.

        Args:
            logger: Logger instance.
            operation: Description of the operation being logged.
        """
        self.logger = logger
        self.operation = operation
        self.start_time = None

    def __enter__(self):
        self.start_time = datetime.now()
        self.logger.info(f"Starting: {self.operation}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = (datetime.now() - self.start_time).total_seconds()
        if exc_type is None:
            self.logger.info(f"Completed: {self.operation} ({duration:.1f}s)")
        else:
            self.logger.error(
                f"Failed: {self.operation} ({duration:.1f}s) - {exc_val}"
            )
        return False  # Don't suppress exceptions
