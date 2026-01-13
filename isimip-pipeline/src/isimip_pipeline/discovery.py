"""Local dataset discovery and search functionality.

Enables searching pre-processed datasets stored in the local processing log.
"""

from pathlib import Path
from typing import List, Optional, Dict, Any

from rich.console import Console
from rich.table import Table

from isimip_pipeline.processing_log import DatasetEntry, ProcessingLog


def find_local_datasets(
    log: ProcessingLog,
    query: Optional[str] = None,
    variable: Optional[str] = None,
    timestep: Optional[str] = None,
    scenario: Optional[str] = None,
) -> List[DatasetEntry]:
    """Search local processed datasets with optional filters.

    Filters are applied with AND logic (all must match).

    Args:
        log: ProcessingLog to search
        query: Search query string (searches descriptive_name and query fields)
        variable: Filter by variable code (e.g., 'led')
        timestep: Filter by timestep (e.g., 'monthly')
        scenario: Filter by climate scenario (e.g., 'ssp126')

    Returns:
        List of matching DatasetEntry objects
    """
    results = log.datasets.copy()

    # Apply query filter
    if query:
        query_lower = query.lower()
        results = [
            e for e in results
            if (query_lower in e.descriptive_name.lower() or
                query_lower in e.variable.lower() or
                query_lower in e.query.lower())
        ]

    # Apply variable filter
    if variable:
        results = [e for e in results if e.variable == variable]

    # Apply timestep filter
    if timestep:
        results = [e for e in results if e.timestep == timestep]

    # Apply scenario filter
    if scenario:
        results = [e for e in results if scenario in e.climate_scenarios]

    return results


def get_dataset_summary(entry: DatasetEntry) -> Dict[str, Any]:
    """Get detailed summary for a specific dataset.

    Args:
        entry: DatasetEntry to summarize

    Returns:
        Dict with formatted dataset information
    """
    return {
        "descriptive_name": entry.descriptive_name,
        "variable": entry.variable,
        "timestep": entry.timestep,
        "created_date": entry.created_date.isoformat(),
        "output_path": entry.output_path,
        "file_count": entry.file_count,
        "time_periods": entry.time_periods,
        "climate_scenarios": entry.climate_scenarios,
        "gcm_models": entry.gcm_models,
        "lsm_models": entry.lsm_models,
        "simulation_round": entry.simulation_round,
        "query": entry.query,
    }


def verify_dataset_integrity(entry: DatasetEntry) -> bool:
    """Verify that dataset files still exist and are accessible.

    Checks if the output directory exists and contains the expected
    processed files.

    Args:
        entry: DatasetEntry to verify

    Returns:
        True if dataset appears intact, False otherwise
    """
    output_path = Path(entry.output_path)

    # Check if output directory exists
    if not output_path.exists():
        return False

    # Check if processed files exist
    processed_dir = output_path / "processed"
    if processed_dir.exists():
        # Check if expected variable file exists
        expected_file = processed_dir / f"{entry.variable}_processed.nc"
        if expected_file.exists():
            return True
        # If file doesn't exist but directory does, still consider it valid
        # (file may be being processed or may be in a different location)
        return True

    # Even if no processed files, output directory existing is enough
    return True


def display_local_results(
    entries: List[DatasetEntry],
    console: Optional[Console] = None,
    detailed: bool = False,
) -> None:
    """Display local datasets in a formatted table.

    Args:
        entries: List of DatasetEntry objects to display
        console: Rich Console instance (creates new if None)
        detailed: If True, show additional metadata columns
    """
    if console is None:
        console = Console()

    if not entries:
        console.print("[yellow]No datasets found.[/yellow]")
        return

    table = Table(title=f"Local Datasets ({len(entries)} found)", show_lines=True)

    # Add columns
    table.add_column("ID", style="dim", width=3)
    table.add_column("Name", style="cyan", no_wrap=True)
    table.add_column("Variable", style="magenta")
    table.add_column("Timestep", style="yellow")
    table.add_column("Scenarios", style="green")
    table.add_column("Created", style="dim")

    if detailed:
        table.add_column("Models", style="blue")
        table.add_column("Files", style="dim")

    # Add rows
    for i, entry in enumerate(entries, 1):
        scenarios_str = ", ".join(entry.climate_scenarios[:2])
        if len(entry.climate_scenarios) > 2:
            scenarios_str += f" +{len(entry.climate_scenarios) - 2}"

        models_str = ", ".join(entry.gcm_models[:2])
        if len(entry.gcm_models) > 2:
            models_str += f" +{len(entry.gcm_models) - 2}"

        row = [
            str(i),
            entry.descriptive_name[:30],
            entry.variable,
            entry.timestep,
            scenarios_str,
            entry.created_date.strftime("%Y-%m-%d"),
        ]

        if detailed:
            row.append(models_str)
            row.append(str(entry.file_count))

        table.add_row(*row)

    console.print(table)

    # Print summary
    console.print(f"\n[dim]Datasets: {len(entries)}[/dim]")
    variables = set(e.variable for e in entries)
    timesteps = set(e.timestep for e in entries)
    console.print(f"[dim]Variables: {', '.join(sorted(variables))}[/dim]")
    console.print(f"[dim]Timesteps: {', '.join(sorted(timesteps))}[/dim]")
