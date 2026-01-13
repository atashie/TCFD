"""Interactive workflow for guided dataset discovery and selection.

Orchestrates the 6-step workflow:
1. Local search - search pre-computed datasets
2. Remote search - if not found locally, search ISIMIP
3. User selection - pick dataset group to process
4. Duplicate check - verify we don't already have it
5-6. Download/process - handled by separate commands
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Any, Optional

from isimip_pipeline.search.isimip_query import DatasetInfo


def save_selection_metadata(
    output_dir: Path,
    datasets: List[DatasetInfo],
    query: str = "",
    descriptive_name: str = "",
) -> None:
    """Save dataset selection to JSON file with metadata.

    Saves selection metadata along with dataset information for later reference
    and to support the workflow steps.

    Args:
        output_dir: Directory to save selection.json
        datasets: List of DatasetInfo objects selected
        query: Original search query
        descriptive_name: User-provided descriptive name
    """
    from datetime import datetime

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    selection_data = {
        "exported_at": datetime.now().isoformat(),
        "query": query,
        "descriptive_name": descriptive_name,
        "total": len(datasets),
        "datasets": [
            {
                "id": ds.id,
                "name": ds.name,
                "url": ds.url,
                "simulation_round": ds.simulation_round,
                "climate_scenario": ds.climate_scenario,
                "variable": ds.variable,
                "model": ds.model,
                "timestep": ds.timestep,
                "size": ds.size,
            }
            for ds in datasets
        ],
    }

    selection_file = output_dir / "selection.json"
    with open(selection_file, "w") as f:
        json.dump(selection_data, f, indent=2)


def load_selection_metadata(output_dir: Path) -> Dict[str, Any]:
    """Load dataset selection from JSON file.

    Args:
        output_dir: Directory containing selection.json

    Returns:
        Dict with selection metadata and datasets

    Raises:
        FileNotFoundError: If selection.json doesn't exist
    """
    output_dir = Path(output_dir)
    selection_file = output_dir / "selection.json"

    if not selection_file.exists():
        raise FileNotFoundError(f"Selection file not found: {selection_file}")

    with open(selection_file, "r") as f:
        return json.load(f)


def parse_descriptive_name_from_folder(folder_name: str) -> str:
    """Extract descriptive name from folder structure.

    Parses folder name following pattern {name}_{variable}-{timestep}
    and extracts the {name} part.

    Args:
        folder_name: Folder name (e.g., "drought-severity_led-monthly")

    Returns:
        Descriptive name part (e.g., "drought-severity")
    """
    # Split on last underscore to separate name from variable-timestep
    if "_" in folder_name:
        parts = folder_name.rsplit("_", 1)
        return parts[0]

    # Fallback: if no underscore, return whole name
    return folder_name


def save_all_available_datasets(
    output_dir: Path,
    datasets: List[DatasetInfo],
    query: str = "",
) -> None:
    """Save all available datasets from search (not just selection).

    Saves the complete set of datasets found during search for reference.

    Args:
        output_dir: Directory to save all_available_datasets.json
        datasets: Complete list of datasets from search
        query: Original search query
    """
    from datetime import datetime

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_data = {
        "exported_at": datetime.now().isoformat(),
        "query": query,
        "total": len(datasets),
        "datasets": [
            {
                "id": ds.id,
                "name": ds.name,
                "url": ds.url,
                "simulation_round": ds.simulation_round,
                "climate_scenario": ds.climate_scenario,
                "variable": ds.variable,
                "model": ds.model,
                "timestep": ds.timestep,
                "size": ds.size,
            }
            for ds in datasets
        ],
    }

    all_file = output_dir / "all_available_datasets.json"
    with open(all_file, "w") as f:
        json.dump(all_data, f, indent=2)


def extract_variable_timestep_from_datasets(
    datasets: List[DatasetInfo],
) -> tuple:
    """Extract variable and timestep from dataset list.

    Gets the first variable and timestep found in the dataset list.

    Args:
        datasets: List of DatasetInfo objects

    Returns:
        Tuple of (variable, timestep), or ("unknown", "unknown") if not found
    """
    if not datasets:
        return ("unknown", "unknown")

    # Get first variable found
    variable = next(
        (ds.variable for ds in datasets if ds.variable),
        "unknown",
    )

    # Get first timestep found
    timestep = next(
        (ds.timestep for ds in datasets if ds.timestep),
        "unknown",
    )

    return (variable, timestep)
