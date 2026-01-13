"""Persistent ISIMIP metrics catalog.

Tracks available variables, scenarios, models, and file counts
across search and download operations.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Any

import yaml


def get_default_catalog_path() -> Path:
    """Get default path for catalog file.

    Returns:
        Path to ~/.isimip-pipeline/catalog.yaml
    """
    return Path.home() / ".isimip-pipeline" / "catalog.yaml"


@dataclass
class VariableInfo:
    """Information about an ISIMIP variable.

    Attributes:
        name: Variable identifier (e.g., 'led', 'burntarea').
        long_name: Human-readable name.
        unit: Unit of measurement.
        scenarios: Set of climate scenarios with data for this variable.
        models: Set of climate models with data for this variable.
        simulation_rounds: Set of ISIMIP simulation rounds.
        file_count: Total number of files seen for this variable.
        last_seen: When this variable was last encountered.
    """

    name: str
    long_name: Optional[str] = None
    unit: Optional[str] = None
    scenarios: Set[str] = field(default_factory=set)
    models: Set[str] = field(default_factory=set)
    simulation_rounds: Set[str] = field(default_factory=set)
    file_count: int = 0
    last_seen: Optional[datetime] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "long_name": self.long_name,
            "unit": self.unit,
            "scenarios": sorted(self.scenarios),
            "models": sorted(self.models),
            "simulation_rounds": sorted(self.simulation_rounds),
            "file_count": self.file_count,
            "last_seen": self.last_seen.isoformat() if self.last_seen else None,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "VariableInfo":
        """Create from dictionary."""
        last_seen = None
        if data.get("last_seen"):
            last_seen = datetime.fromisoformat(data["last_seen"])

        return cls(
            name=data["name"],
            long_name=data.get("long_name"),
            unit=data.get("unit"),
            scenarios=set(data.get("scenarios", [])),
            models=set(data.get("models", [])),
            simulation_rounds=set(data.get("simulation_rounds", [])),
            file_count=data.get("file_count", 0),
            last_seen=last_seen,
        )


class ISIMIPCatalog:
    """Catalog of ISIMIP variables and metadata.

    Maintains a persistent record of variables, scenarios, and models
    encountered during search and download operations.
    """

    def __init__(self):
        """Initialize empty catalog."""
        self.variables: Dict[str, VariableInfo] = {}
        self.last_updated: Optional[datetime] = None

    def add_variable(
        self,
        name: str,
        long_name: Optional[str] = None,
        unit: Optional[str] = None,
    ) -> VariableInfo:
        """Add or update a variable in the catalog.

        Args:
            name: Variable identifier.
            long_name: Human-readable name.
            unit: Unit of measurement.

        Returns:
            The VariableInfo object.
        """
        if name not in self.variables:
            self.variables[name] = VariableInfo(
                name=name,
                long_name=long_name,
                unit=unit,
            )
        else:
            # Update metadata if provided
            if long_name:
                self.variables[name].long_name = long_name
            if unit:
                self.variables[name].unit = unit

        return self.variables[name]

    def update_from_datasets(self, datasets: List[Any]) -> None:
        """Update catalog from list of DatasetInfo objects.

        Args:
            datasets: List of DatasetInfo from search results.
        """
        for ds in datasets:
            variable = getattr(ds, "variable", None)
            if not variable:
                # Try to extract from name
                continue

            # Add or get variable
            var_info = self.add_variable(variable)

            # Update scenarios
            scenario = getattr(ds, "climate_scenario", None)
            if scenario:
                var_info.scenarios.add(scenario)

            # Update models
            model = getattr(ds, "model", None)
            if model:
                var_info.models.add(model)

            # Update simulation rounds
            sim_round = getattr(ds, "simulation_round", None)
            if sim_round:
                var_info.simulation_rounds.add(sim_round)

            # Increment file count
            var_info.file_count += 1
            var_info.last_seen = datetime.now()

        self.last_updated = datetime.now()

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics for the catalog.

        Returns:
            Dictionary with summary statistics.
        """
        all_scenarios: Set[str] = set()
        all_models: Set[str] = set()
        total_files = 0

        for var_info in self.variables.values():
            all_scenarios.update(var_info.scenarios)
            all_models.update(var_info.models)
            total_files += var_info.file_count

        return {
            "total_variables": len(self.variables),
            "total_files": total_files,
            "all_scenarios": sorted(all_scenarios),
            "all_models": sorted(all_models),
            "last_updated": self.last_updated.isoformat() if self.last_updated else None,
        }

    def to_dict(self) -> Dict[str, Any]:
        """Convert catalog to dictionary for serialization."""
        return {
            "last_updated": self.last_updated.isoformat() if self.last_updated else None,
            "variables": {
                name: var.to_dict() for name, var in self.variables.items()
            },
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ISIMIPCatalog":
        """Create catalog from dictionary."""
        catalog = cls()

        if data.get("last_updated"):
            catalog.last_updated = datetime.fromisoformat(data["last_updated"])

        for name, var_data in data.get("variables", {}).items():
            catalog.variables[name] = VariableInfo.from_dict(var_data)

        return catalog


def save_catalog(catalog: ISIMIPCatalog, path: Optional[Path] = None) -> Path:
    """Save catalog to YAML file.

    Args:
        catalog: Catalog to save.
        path: Path to save to (default: ~/.isimip-pipeline/catalog.yaml).

    Returns:
        Path where catalog was saved.
    """
    if path is None:
        path = get_default_catalog_path()

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        yaml.safe_dump(catalog.to_dict(), f, default_flow_style=False, sort_keys=False)

    return path


def load_catalog(path: Optional[Path] = None) -> ISIMIPCatalog:
    """Load catalog from YAML file.

    Args:
        path: Path to load from (default: ~/.isimip-pipeline/catalog.yaml).

    Returns:
        Loaded catalog, or empty catalog if file doesn't exist.
    """
    if path is None:
        path = get_default_catalog_path()

    path = Path(path)

    if not path.exists():
        return ISIMIPCatalog()

    with open(path) as f:
        data = yaml.safe_load(f) or {}

    return ISIMIPCatalog.from_dict(data)


def update_catalog_from_datasets(
    datasets: List[Any],
    path: Optional[Path] = None,
) -> ISIMIPCatalog:
    """Load catalog, update with datasets, and save.

    Convenience function for updating the persistent catalog.

    Args:
        datasets: List of DatasetInfo objects.
        path: Path to catalog file.

    Returns:
        Updated catalog.
    """
    catalog = load_catalog(path)
    catalog.update_from_datasets(datasets)
    save_catalog(catalog, path)
    return catalog
