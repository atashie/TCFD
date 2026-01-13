"""Custom exceptions and error handling for the ISIMIP pipeline.

Provides user-friendly error messages with actionable guidance for common issues.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, List


class ErrorCategory(Enum):
    """Categories of errors for grouping related issues."""

    FILE_NOT_FOUND = "file_not_found"
    CONFIGURATION = "configuration"
    NETWORK = "network"
    DATA_QUALITY = "data_quality"
    PROCESSING = "processing"
    VALIDATION = "validation"


@dataclass
class ErrorInfo:
    """Structured error information with actionable guidance.

    Attributes:
        message: Brief description of the error
        category: Type of error for categorization
        suggestions: List of actionable steps to resolve the issue
        details: Optional technical details for debugging
    """

    message: str
    category: ErrorCategory
    suggestions: List[str]
    details: Optional[str] = None

    def format(self) -> str:
        """Format error for display."""
        lines = [f"Error: {self.message}"]

        if self.suggestions:
            lines.append("\nSuggested actions:")
            for i, suggestion in enumerate(self.suggestions, 1):
                lines.append(f"  {i}. {suggestion}")

        if self.details:
            lines.append(f"\nDetails: {self.details}")

        return "\n".join(lines)


class PipelineError(Exception):
    """Base exception for pipeline errors with actionable guidance."""

    def __init__(self, info: ErrorInfo):
        self.info = info
        super().__init__(info.format())


class ConfigurationError(PipelineError):
    """Configuration-related errors."""
    pass


class FileNotFoundError(PipelineError):
    """File or directory not found errors."""
    pass


class NetworkError(PipelineError):
    """Network and API communication errors."""
    pass


class DataQualityError(PipelineError):
    """Data quality and validation errors."""
    pass


class ProcessingError(PipelineError):
    """Data processing errors."""
    pass


# =============================================================================
# Error Factory Functions - Create structured errors with guidance
# =============================================================================

def config_not_found_error(path: Path) -> ErrorInfo:
    """Create error info for missing configuration file."""
    return ErrorInfo(
        message=f"Configuration file not found: {path}",
        category=ErrorCategory.CONFIGURATION,
        suggestions=[
            f"Create config file at: {path}",
            "Copy from config.example.yaml and customize",
            "Run: isimip-pipeline --help for configuration details",
            "Use --config flag to specify alternate location",
        ],
    )


def config_invalid_error(path: Path, detail: str) -> ErrorInfo:
    """Create error info for invalid configuration."""
    return ErrorInfo(
        message=f"Invalid configuration in: {path}",
        category=ErrorCategory.CONFIGURATION,
        suggestions=[
            "Check YAML syntax (indentation, colons)",
            "Verify all required fields are present",
            "Compare against config.example.yaml template",
        ],
        details=detail,
    )


def api_key_missing_error() -> ErrorInfo:
    """Create error info for missing API key."""
    return ErrorInfo(
        message="API key not configured",
        category=ErrorCategory.CONFIGURATION,
        suggestions=[
            "Add 'you_api_key' to ~/.isimip-pipeline/config.yaml",
            "Use --no-llm flag to skip LLM query parsing",
            "Pipeline will use keyword fallback for query parsing",
        ],
    )


def directory_not_found_error(path: Path, purpose: str) -> ErrorInfo:
    """Create error info for missing directory."""
    return ErrorInfo(
        message=f"{purpose} directory not found: {path}",
        category=ErrorCategory.FILE_NOT_FOUND,
        suggestions=[
            f"Create the directory: mkdir -p {path}",
            "Check the path is correct",
            "Use absolute path if relative path is ambiguous",
        ],
    )


def no_netcdf_files_error(directory: Path) -> ErrorInfo:
    """Create error info when no NetCDF files found."""
    return ErrorInfo(
        message=f"No NetCDF files found in: {directory}",
        category=ErrorCategory.FILE_NOT_FOUND,
        suggestions=[
            "Verify files have .nc or .nc4 extension",
            "Check directory path is correct",
            "Run download command first if files not yet downloaded",
            "Example: isimip-pipeline download -s selection.json",
        ],
    )


def selection_file_not_found_error(path: Path) -> ErrorInfo:
    """Create error info for missing selection file."""
    return ErrorInfo(
        message=f"Selection file not found: {path}",
        category=ErrorCategory.FILE_NOT_FOUND,
        suggestions=[
            "Run search command to create selection file",
            "Example: isimip-pipeline search 'query' -o selection.json",
            "Or use interactive workflow: isimip-pipeline interactive 'query'",
            "Check the file path is correct",
        ],
    )


def processing_log_not_found_error() -> ErrorInfo:
    """Create error info when processing log doesn't exist."""
    return ErrorInfo(
        message="No processed datasets found in local log",
        category=ErrorCategory.FILE_NOT_FOUND,
        suggestions=[
            "Run 'isimip-pipeline run' to process datasets",
            "Or use interactive workflow: isimip-pipeline interactive",
            "Log will be created at: outputs/processed_data_log.yaml",
        ],
    )


def variable_not_found_error(variable: str, available: List[str]) -> ErrorInfo:
    """Create error info when specified variable not in files."""
    return ErrorInfo(
        message=f"Variable '{variable}' not found in input files",
        category=ErrorCategory.PROCESSING,
        suggestions=[
            f"Available variables: {', '.join(available)}",
            "Check spelling of variable name",
            "Omit --variable flag to auto-detect from files",
        ],
    )


def isimip_search_error(detail: str) -> ErrorInfo:
    """Create error info for ISIMIP API search failures."""
    return ErrorInfo(
        message="Failed to search ISIMIP repository",
        category=ErrorCategory.NETWORK,
        suggestions=[
            "Check internet connection",
            "ISIMIP server may be temporarily unavailable",
            "Try again in a few minutes",
            "Use --no-llm flag if LLM parsing is failing",
        ],
        details=detail,
    )


def no_search_results_error(query: str) -> ErrorInfo:
    """Create error info when search returns no results."""
    return ErrorInfo(
        message=f"No datasets found for query: '{query}'",
        category=ErrorCategory.PROCESSING,
        suggestions=[
            "Try broader search terms",
            "Check variable code spelling (e.g., 'led' not 'LED')",
            "Use --refresh flag to search for variable definitions",
            "Try: isimip-pipeline catalog --search <term>",
        ],
    )


def download_failed_error(url: str, detail: str) -> ErrorInfo:
    """Create error info for download failures."""
    return ErrorInfo(
        message=f"Download failed: {url}",
        category=ErrorCategory.NETWORK,
        suggestions=[
            "Check internet connection",
            "File may have moved or been removed",
            "Try running download again (will resume where it left off)",
            "Check ISIMIP repository directly if issue persists",
        ],
        details=detail,
    )


def data_quality_error(variable: str, quality: str, loss_pct: float) -> ErrorInfo:
    """Create error info for poor data quality."""
    return ErrorInfo(
        message=f"Data quality issue for '{variable}': {quality} ({loss_pct:.1f}% data loss)",
        category=ErrorCategory.DATA_QUALITY,
        suggestions=[
            "Review validation report for details",
            "Consider using a different data source",
            "Check for spatial/temporal gaps in input data",
            "Use --skip-validation flag to proceed anyway (not recommended)",
        ],
    )


def spatial_grid_mismatch_error(expected: str, actual: str) -> ErrorInfo:
    """Create error info for spatial grid mismatches."""
    return ErrorInfo(
        message="Spatial grid mismatch between datasets",
        category=ErrorCategory.PROCESSING,
        suggestions=[
            "Ensure all input files use same grid resolution",
            "Check coordinate names (lat/lon vs latitude/longitude)",
            "Consider regridding datasets to common grid",
            "Use alignment module with regrid=True option",
        ],
        details=f"Expected: {expected}, Found: {actual}",
    )


def calendar_conversion_error(source_cal: str, target_cal: str) -> ErrorInfo:
    """Create error info for calendar conversion issues."""
    return ErrorInfo(
        message=f"Calendar conversion failed: {source_cal} -> {target_cal}",
        category=ErrorCategory.PROCESSING,
        suggestions=[
            "Check source calendar type in dataset metadata",
            "Supported calendars: standard, 360_day, noleap",
            "Some date ranges may not convert between calendars",
            "Try processing datasets with same calendar separately",
        ],
    )


def insufficient_data_error(variable: str, n_timesteps: int, required: int) -> ErrorInfo:
    """Create error info when insufficient data for processing."""
    return ErrorInfo(
        message=f"Insufficient data for '{variable}': {n_timesteps} timesteps (need {required}+)",
        category=ErrorCategory.PROCESSING,
        suggestions=[
            "Download more data files for this variable",
            "Check if all years are present in input data",
            "Processing requires multi-decade time series",
            "Use a longer time period in search query",
        ],
    )


# =============================================================================
# Rich Console Formatting
# =============================================================================

def format_error_rich(info: ErrorInfo) -> str:
    """Format error info for Rich console output."""
    from rich.panel import Panel
    from rich.text import Text

    lines = []
    lines.append(f"[bold red]Error:[/bold red] {info.message}")

    if info.suggestions:
        lines.append("\n[bold yellow]Suggested actions:[/bold yellow]")
        for i, suggestion in enumerate(info.suggestions, 1):
            lines.append(f"  [dim]{i}.[/dim] {suggestion}")

    if info.details:
        lines.append(f"\n[dim]Details: {info.details}[/dim]")

    return "\n".join(lines)


def print_error(info: ErrorInfo, console=None) -> None:
    """Print formatted error to console."""
    if console is None:
        from rich.console import Console
        console = Console()

    console.print(format_error_rich(info))
