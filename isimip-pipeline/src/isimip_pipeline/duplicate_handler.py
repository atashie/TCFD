"""Duplicate detection and resolution for datasets.

Handles detection of duplicate datasets based on variable+timestep
and provides utilities for resolving conflicts.
"""

import re
from enum import Enum
from typing import Optional

from isimip_pipeline.processing_log import DatasetEntry, ProcessingLog


class DuplicateAction(Enum):
    """Actions to take when duplicate is detected."""

    SKIP = "skip"               # Skip download, use existing dataset
    NEW_FOLDER = "new"          # Create new folder with different name
    OVERWRITE = "overwrite"     # Replace existing dataset
    ABORT = "abort"             # Cancel operation


def check_for_duplicate(
    log: ProcessingLog,
    variable: str,
    timestep: str,
) -> Optional[DatasetEntry]:
    """Check if dataset with same variable+timestep exists.

    Args:
        log: ProcessingLog to search
        variable: Variable code (e.g., 'led')
        timestep: Timestep (e.g., 'monthly')

    Returns:
        Existing DatasetEntry if found, None otherwise
    """
    return log.find_by_variable_timestep(variable, timestep)


def generate_unique_name(
    base_name: str,
    variable: str,
    timestep: str,
    log: ProcessingLog,
) -> str:
    """Generate unique folder name by appending suffix if needed.

    Creates a folder name following the pattern {name}_{variable}-{timestep}.
    If a duplicate exists, appends -2, -3, etc. until a unique name is found.

    Args:
        base_name: Descriptive name to base unique name on
        variable: Variable code
        timestep: Timestep
        log: ProcessingLog to check for conflicts

    Returns:
        Unique folder name (does not include .../outputs/ prefix)
    """
    # Clean base name
    clean_name = re.sub(r'[^a-z0-9-]', '', base_name.lower())
    clean_name = re.sub(r'-+', '-', clean_name).strip('-')

    # Build folder name template
    base_folder = f"{clean_name}_{variable}-{timestep}"

    # Check if any entry with this variable+timestep exists
    potential_dup = log.find_by_variable_timestep(variable, timestep)
    if potential_dup is None:
        # No conflict - safe to use base name
        return base_folder

    # Duplicate exists - need to generate unique suffix
    # Find the highest numbered variant that exists
    max_num = 2
    for entry in log.datasets:
        if entry.variable == variable and entry.timestep == timestep:
            # Check if name matches pattern: {name}-N
            pattern = re.escape(clean_name) + r"-(\d+)"
            match = re.search(pattern, entry.descriptive_name)
            if match:
                num = int(match.group(1))
                max_num = max(max_num, num + 1)

    # Return name with suffix
    numbered_name = f"{clean_name}-{max_num}"
    return f"{numbered_name}_{variable}-{timestep}"


def build_output_folder_name(
    descriptive_name: str,
    variable: str,
    timestep: str,
) -> str:
    """Build standardized folder name.

    Creates a folder name following the pattern {name}_{variable}-{timestep}
    with proper cleaning and normalization.

    Args:
        descriptive_name: User-friendly dataset name
        variable: Variable code
        timestep: Timestep

    Returns:
        Folder name (does not include .../outputs/ prefix)
    """
    # Clean descriptive name
    clean_name = re.sub(r'[^a-z0-9-]', '', descriptive_name.lower())
    clean_name = re.sub(r'-+', '-', clean_name).strip('-')

    return f"{clean_name}_{variable}-{timestep}"
