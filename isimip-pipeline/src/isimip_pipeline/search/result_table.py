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
