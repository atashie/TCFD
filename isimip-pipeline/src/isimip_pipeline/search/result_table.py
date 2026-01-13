"""Result table display and export for ISIMIP search results."""

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

from rich.table import Table
from rich.console import Console

from isimip_pipeline.search.isimip_query import DatasetInfo


class ResultTable:
    """Formats and displays search results as a table.

    Provides methods to display results in terminal and export to JSON.
    """

    def __init__(self, datasets: List[DatasetInfo]):
        """Initialize with list of datasets.

        Args:
            datasets: List of DatasetInfo objects from search.
        """
        self.datasets = datasets

    def to_rich_table(self, title: str = "Search Results") -> Table:
        """Convert datasets to Rich Table for terminal display.

        Args:
            title: Table title.

        Returns:
            Rich Table object.
        """
        table = Table(title=title, show_lines=True)

        # Add columns
        table.add_column("#", style="dim", width=4)
        table.add_column("Name", style="cyan", no_wrap=True)
        table.add_column("Round", style="green")
        table.add_column("Scenario", style="yellow")
        table.add_column("Variable", style="magenta")
        table.add_column("Model", style="blue")
        table.add_column("Timestep", style="dim")

        # Add rows
        for i, ds in enumerate(self.datasets, 1):
            table.add_row(
                str(i),
                ds.name[:40] + "..." if len(ds.name) > 40 else ds.name,
                ds.simulation_round or "-",
                ds.climate_scenario or "-",
                ds.variable or "-",
                ds.model or "-",
                ds.timestep or "-",
            )

        return table

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics for the datasets.

        Returns:
            Dict with summary information.
        """
        scenarios = set()
        variables = set()
        models = set()
        rounds = set()
        timesteps = set()

        for ds in self.datasets:
            if ds.climate_scenario:
                scenarios.add(ds.climate_scenario)
            if ds.variable:
                variables.add(ds.variable)
            if ds.model:
                models.add(ds.model)
            if ds.simulation_round:
                rounds.add(ds.simulation_round)
            if ds.timestep:
                timesteps.add(ds.timestep)

        return {
            "total_datasets": len(self.datasets),
            "scenarios": sorted(scenarios),
            "variables": sorted(variables),
            "models": sorted(models),
            "simulation_rounds": sorted(rounds),
            "timesteps": sorted(timesteps),
        }

    def print_summary(self, console: Optional[Console] = None):
        """Print summary to console.

        Args:
            console: Rich Console instance (creates new if None).
        """
        if console is None:
            console = Console()

        summary = self.get_summary()

        console.print(f"\n[bold]Found {summary['total_datasets']} datasets[/bold]")

        if summary["scenarios"]:
            console.print(f"  Scenarios: {', '.join(summary['scenarios'])}")
        if summary["variables"]:
            console.print(f"  Variables: {', '.join(summary['variables'])}")
        if summary["models"]:
            console.print(f"  Models: {', '.join(summary['models'])}")
        if summary["timesteps"]:
            console.print(f"  Timesteps: {', '.join(summary['timesteps'])}")

    def display(self, console: Optional[Console] = None):
        """Display table and summary to console.

        Args:
            console: Rich Console instance (creates new if None).
        """
        if console is None:
            console = Console()

        if not self.datasets:
            console.print("[yellow]No datasets found.[/yellow]")
            return

        table = self.to_rich_table()
        console.print(table)
        self.print_summary(console)


def group_by_attributes(
    datasets: List[DatasetInfo],
    attribute: str,
) -> Dict[str, List[DatasetInfo]]:
    """Group datasets by a specific attribute.

    Args:
        datasets: List of DatasetInfo objects.
        attribute: Attribute name to group by (e.g., "climate_scenario").

    Returns:
        Dict mapping attribute values to lists of datasets.
    """
    grouped = defaultdict(list)

    for ds in datasets:
        value = getattr(ds, attribute, None)
        if value is not None:
            grouped[value].append(ds)
        else:
            grouped["unknown"].append(ds)

    return dict(grouped)


def export_selection(
    datasets: List[DatasetInfo],
    output_path: Path,
    query: str = "",
) -> None:
    """Export dataset selection to JSON file.

    Args:
        datasets: List of DatasetInfo objects to export.
        output_path: Path to output JSON file.
        query: Original search query (for metadata).
    """
    export_data = {
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

    with open(output_path, "w") as f:
        json.dump(export_data, f, indent=2)


def format_results(datasets: List[DatasetInfo]) -> str:
    """Format results as a string for display.

    Args:
        datasets: List of DatasetInfo objects.

    Returns:
        Formatted string representation.
    """
    if not datasets:
        return "No datasets found matching your query."

    lines = [f"Found {len(datasets)} datasets:\n"]

    for i, ds in enumerate(datasets[:20], 1):  # Limit to first 20
        lines.append(
            f"  {i}. {ds.name}"
            f" [{ds.climate_scenario or '-'}]"
            f" [{ds.variable or '-'}]"
        )

    if len(datasets) > 20:
        lines.append(f"\n  ... and {len(datasets) - 20} more")

    return "\n".join(lines)


def group_by_variable_timestep(
    datasets: List[DatasetInfo],
) -> Dict[tuple, Dict[str, Any]]:
    """Group datasets by (variable, timestep) with metadata.

    Creates a grouped structure that shows the different variable+timestep
    combinations available in the search results, with summary information
    about each group.

    Args:
        datasets: List of DatasetInfo objects.

    Returns:
        Dict mapping (variable, timestep) tuples to group metadata dicts.
        Each group dict contains:
        - 'datasets': List[DatasetInfo]
        - 'count': int (number of files)
        - 'scenarios': List[str]
        - 'models': List[str]
        - 'simulation_rounds': List[str]
        - 'time_coverage': str (e.g., '2006-2100')
    """
    grouped = defaultdict(lambda: {
        'datasets': [],
        'scenarios': set(),
        'models': set(),
        'simulation_rounds': set(),
    })

    for ds in datasets:
        var = ds.variable or 'unknown'
        ts = ds.timestep or 'unknown'
        key = (var, ts)

        grouped[key]['datasets'].append(ds)
        if ds.climate_scenario:
            grouped[key]['scenarios'].add(ds.climate_scenario)
        if ds.model:
            grouped[key]['models'].add(ds.model)
        if ds.simulation_round:
            grouped[key]['simulation_rounds'].add(ds.simulation_round)

    # Convert sets to sorted lists and add count
    result = {}
    for key, group in grouped.items():
        result[key] = {
            'datasets': group['datasets'],
            'count': len(group['datasets']),
            'scenarios': sorted(group['scenarios']),
            'models': sorted(group['models']),
            'simulation_rounds': sorted(group['simulation_rounds']),
            'time_coverage': '2006-2100',  # Placeholder; could be extracted
        }

    return result


def display_grouped_results(
    grouped: Dict[tuple, Dict[str, Any]],
    console: Optional[Console] = None,
) -> None:
    """Display grouped results with summary per group.

    Shows each variable+timestep combination as a group with details
    about scenarios, models, and file count.

    Args:
        grouped: Dict from group_by_variable_timestep.
        console: Rich Console instance (creates new if None).
    """
    if console is None:
        console = Console()

    if not grouped:
        console.print("[yellow]No dataset groups found.[/yellow]")
        return

    console.print(f"\n[bold]Found {len(grouped)} dataset group(s):[/bold]\n")

    for idx, ((variable, timestep), group) in enumerate(sorted(grouped.items()), 1):
        console.print(f"[cyan]{idx}. {variable} ({timestep})[/cyan]")
        console.print(f"   Files: {group['count']}")

        if group['scenarios']:
            scenarios_str = ', '.join(group['scenarios'][:5])
            if len(group['scenarios']) > 5:
                scenarios_str += f", ... ({len(group['scenarios'])} total)"
            console.print(f"   Scenarios: {scenarios_str}")

        if group['models']:
            models_str = ', '.join(group['models'][:3])
            if len(group['models']) > 3:
                models_str += f", ... ({len(group['models'])} total)"
            console.print(f"   Models: {models_str}")

        if group['simulation_rounds']:
            console.print(f"   Rounds: {', '.join(group['simulation_rounds'])}")

        if group['time_coverage']:
            console.print(f"   Coverage: {group['time_coverage']}")

        console.print()
