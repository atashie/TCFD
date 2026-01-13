"""Configuration loading and management for ISIMIP pipeline."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import yaml


@dataclass
class PathsConfig:
    """Path configuration settings."""

    download_dir: Path = field(default_factory=lambda: Path("./data/raw"))
    processed_dir: Path = field(default_factory=lambda: Path("./data/processed"))
    reports_dir: Path = field(default_factory=lambda: Path("./reports"))


@dataclass
class ApiConfig:
    """API configuration settings."""

    you_api_key: Optional[str] = None
    you_agent_id: Optional[str] = None
    isimip_timeout: int = 30


@dataclass
class ProcessingConfig:
    """Processing configuration settings."""

    smoothing_bandwidth: int = 15
    trend_method: str = "theil_sen"
    significance_method: str = "spearman"
    percentile_bins: int = 100
    decades: list = field(default_factory=lambda: [10, 20, 30, 40, 50, 60, 70, 80, 90])


@dataclass
class DownloadConfig:
    """Download configuration settings."""

    max_concurrent: int = 4
    chunk_size: int = 8192
    retry_attempts: int = 3


@dataclass
class VisualizationConfig:
    """Visualization configuration settings."""

    map_projection: str = "equirectangular"
    color_scale: str = "viridis"


@dataclass
class Config:
    """Main configuration container."""

    paths: PathsConfig = field(default_factory=PathsConfig)
    api: ApiConfig = field(default_factory=ApiConfig)
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    download: DownloadConfig = field(default_factory=DownloadConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)


# Valid values for validation
VALID_TREND_METHODS = {"theil_sen", "ols"}
VALID_SIGNIFICANCE_METHODS = {"spearman", "pearson"}


# Default configuration instance
DEFAULT_CONFIG = Config()


def load_config(config_path: Path) -> Config:
    """Load configuration from YAML file, falling back to defaults.

    Args:
        config_path: Path to YAML configuration file.

    Returns:
        Config object with values from file merged with defaults.

    Raises:
        ValueError: If configuration values are invalid.
    """
    config = Config()

    if not config_path.exists():
        return config

    with open(config_path, "r") as f:
        data = yaml.safe_load(f) or {}

    # Load paths
    if "paths" in data:
        paths_data = data["paths"]
        if "download_dir" in paths_data:
            config.paths.download_dir = Path(paths_data["download_dir"])
        if "processed_dir" in paths_data:
            config.paths.processed_dir = Path(paths_data["processed_dir"])
        if "reports_dir" in paths_data:
            config.paths.reports_dir = Path(paths_data["reports_dir"])

    # Load API config
    if "api" in data:
        api_data = data["api"]
        if "you_api_key" in api_data:
            config.api.you_api_key = api_data["you_api_key"]
        if "you_agent_id" in api_data:
            config.api.you_agent_id = api_data["you_agent_id"]
        if "isimip_timeout" in api_data:
            config.api.isimip_timeout = api_data["isimip_timeout"]

    # Load processing config with validation
    if "processing" in data:
        proc_data = data["processing"]
        if "smoothing_bandwidth" in proc_data:
            config.processing.smoothing_bandwidth = proc_data["smoothing_bandwidth"]
        if "trend_method" in proc_data:
            method = proc_data["trend_method"]
            if method not in VALID_TREND_METHODS:
                raise ValueError(
                    f"Invalid trend_method '{method}'. "
                    f"Must be one of: {VALID_TREND_METHODS}"
                )
            config.processing.trend_method = method
        if "significance_method" in proc_data:
            config.processing.significance_method = proc_data["significance_method"]
        if "percentile_bins" in proc_data:
            config.processing.percentile_bins = proc_data["percentile_bins"]
        if "decades" in proc_data:
            config.processing.decades = proc_data["decades"]

    # Load download config
    if "download" in data:
        dl_data = data["download"]
        if "max_concurrent" in dl_data:
            config.download.max_concurrent = dl_data["max_concurrent"]
        if "chunk_size" in dl_data:
            config.download.chunk_size = dl_data["chunk_size"]
        if "retry_attempts" in dl_data:
            config.download.retry_attempts = dl_data["retry_attempts"]

    # Load visualization config
    if "visualization" in data:
        viz_data = data["visualization"]
        if "map_projection" in viz_data:
            config.visualization.map_projection = viz_data["map_projection"]
        if "color_scale" in viz_data:
            config.visualization.color_scale = viz_data["color_scale"]

    return config
