"""Processing log management for tracking processed datasets.

Maintains a YAML-based log of all processed datasets with comprehensive metadata.
"""

from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any

import yaml


@dataclass
class DatasetEntry:
    """Information about a processed dataset.

    Attributes:
        descriptive_name: User-friendly name for the dataset
        variable: Variable code (e.g., 'led', 'burntarea')
        timestep: Time resolution (e.g., 'monthly', 'annual', 'daily')
        created_date: When the dataset was processed
        output_path: Path to output directory
        file_count: Number of NetCDF files
        time_periods: List of time period ranges (e.g., ['2006-2100'])
        climate_scenarios: Climate scenarios included (e.g., ['ssp126', 'ssp370'])
        gcm_models: Global Climate Models used
        lsm_models: Land Surface Models used
        simulation_round: ISIMIP simulation round (e.g., 'ISIMIP3b')
        query: Original search query
    """

    descriptive_name: str
    variable: str
    timestep: str
    created_date: datetime
    output_path: str
    file_count: int
    time_periods: List[str] = field(default_factory=list)
    climate_scenarios: List[str] = field(default_factory=list)
    gcm_models: List[str] = field(default_factory=list)
    lsm_models: List[str] = field(default_factory=list)
    simulation_round: str = ""
    query: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert entry to dictionary for YAML serialization."""
        data = asdict(self)
        # Convert datetime to ISO format string
        if isinstance(data["created_date"], datetime):
            data["created_date"] = data["created_date"].isoformat()
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DatasetEntry":
        """Create entry from dictionary (e.g., from YAML)."""
        data = data.copy()
        # Convert ISO format string back to datetime
        if isinstance(data.get("created_date"), str):
            data["created_date"] = datetime.fromisoformat(data["created_date"])
        return cls(**data)


@dataclass
class ProcessingLog:
    """Container for processing log entries.

    Manages a list of processed datasets with search and duplicate detection.
    """

    datasets: List[DatasetEntry] = field(default_factory=list)
    last_updated: datetime = field(default_factory=datetime.now)

    def add_entry(self, entry: DatasetEntry) -> None:
        """Add a new dataset entry to the log."""
        self.datasets.append(entry)
        self.last_updated = datetime.now()

    def find_by_variable_timestep(
        self, variable: str, timestep: str
    ) -> Optional[DatasetEntry]:
        """Find dataset by variable+timestep combination.

        Args:
            variable: Variable code
            timestep: Time resolution

        Returns:
            First matching entry or None
        """
        for entry in self.datasets:
            if entry.variable == variable and entry.timestep == timestep:
                return entry
        return None

    def search(self, query: str) -> List[DatasetEntry]:
        """Search log by query string.

        Searches in descriptive_name, variable, and query fields.
        Case-insensitive.

        Args:
            query: Search string

        Returns:
            List of matching entries
        """
        query_lower = query.lower()
        results = []

        for entry in self.datasets:
            # Search in descriptive name
            if query_lower in entry.descriptive_name.lower():
                results.append(entry)
                continue
            # Search in variable
            if query_lower in entry.variable.lower():
                results.append(entry)
                continue
            # Search in original query
            if query_lower in entry.query.lower():
                results.append(entry)

        return results

    def to_dict(self) -> Dict[str, Any]:
        """Convert log to dictionary for YAML serialization."""
        return {
            "datasets": [e.to_dict() for e in self.datasets],
            "last_updated": self.last_updated.isoformat(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ProcessingLog":
        """Create log from dictionary (e.g., from YAML)."""
        entries = [DatasetEntry.from_dict(e) for e in data.get("datasets", [])]
        last_updated = data.get("last_updated")
        if isinstance(last_updated, str):
            last_updated = datetime.fromisoformat(last_updated)
        else:
            last_updated = datetime.now()

        log = cls(datasets=entries, last_updated=last_updated)
        return log


def load_processing_log(path: Path) -> ProcessingLog:
    """Load processing log from YAML file.

    Creates empty log if file doesn't exist.

    Args:
        path: Path to YAML log file

    Returns:
        ProcessingLog instance
    """
    path = Path(path)

    if not path.exists():
        return ProcessingLog()

    with open(path, "r") as f:
        data = yaml.safe_load(f)

    if data is None:
        return ProcessingLog()

    return ProcessingLog.from_dict(data)


def save_processing_log(log: ProcessingLog, path: Path) -> None:
    """Save processing log to YAML file.

    Performs atomic write to prevent corruption.

    Args:
        log: ProcessingLog to save
        path: Path to write YAML file
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Write to temporary file first for atomic operation
    temp_path = path.with_suffix(".tmp")

    try:
        with open(temp_path, "w") as f:
            yaml.safe_dump(log.to_dict(), f, default_flow_style=False)

        # Atomic rename
        temp_path.replace(path)
    except Exception:
        # Clean up temp file on error
        if temp_path.exists():
            temp_path.unlink()
        raise


def find_duplicate(
    log: ProcessingLog, variable: str, timestep: str
) -> Optional[DatasetEntry]:
    """Convenience function to find duplicate in log.

    Args:
        log: ProcessingLog to search
        variable: Variable code
        timestep: Time resolution

    Returns:
        Matching entry or None
    """
    return log.find_by_variable_timestep(variable, timestep)


def add_dataset_entry(log: ProcessingLog, entry: DatasetEntry) -> ProcessingLog:
    """Convenience function to add entry to log.

    Args:
        log: ProcessingLog to update
        entry: DatasetEntry to add

    Returns:
        Updated ProcessingLog
    """
    log.add_entry(entry)
    return log


def search_log(log: ProcessingLog, query: str) -> List[DatasetEntry]:
    """Convenience function to search log.

    Args:
        log: ProcessingLog to search
        query: Search string

    Returns:
        List of matching entries
    """
    return log.search(query)
