"""CLI entry points for ISIMIP pipeline."""

import json
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel

from isimip_pipeline import __version__
from isimip_pipeline.config import load_config, Config
from isimip_pipeline.processing_log import load_processing_log
from isimip_pipeline.discovery import find_local_datasets, display_local_results
from isimip_pipeline.duplicate_handler import (
    build_output_folder_name,
    check_for_duplicate,
)
from isimip_pipeline.interactive import (
    save_selection_metadata,
    save_all_available_datasets,
    extract_variable_timestep_from_datasets,
)
from isimip_pipeline.search.result_table import (
    group_by_variable_timestep,
    display_grouped_results,
)

app = typer.Typer(
    name="isimip-pipeline",
    help="Automated ISIMIP climate data discovery, download, and processing.",
    add_completion=False,
)

console = Console()


def get_config(config_path: Optional[Path]) -> Config:
    """Load config from path or default location."""
    if config_path:
        return load_config(config_path)

    # Try default locations
    default_paths = [
        Path.home() / ".isimip-pipeline" / "config.yaml",
        Path("config.yaml"),
    ]

    for path in default_paths:
        if path.exists():
            return load_config(path)

    return Config()


def version_callback(value: bool):
    """Print version and exit."""
    if value:
        console.print(f"isimip-pipeline version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
):
    """ISIMIP Pipeline - Climate data automation tool."""
    pass


@app.command()
def search(
    query: str = typer.Argument(..., help="Natural language search query"),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Save results to JSON file"
    ),
    use_llm: bool = typer.Option(
        True, "--llm/--no-llm", help="Use LLM for query parsing"
    ),
    refresh: bool = typer.Option(
        False, "--refresh", "-r", help="Force web search for variable definitions"
    ),
    limit: int = typer.Option(
        50, "--limit", "-l", help="Maximum number of results to show"
    ),
):
    """Search ISIMIP repository using natural language query.

    Examples:
        isimip-pipeline search "drought exposure metrics"
        isimip-pipeline search "wildfire burned area ssp585"
        isimip-pipeline search "led" --no-llm
        isimip-pipeline search "permafrost variables" --refresh
    """
    from isimip_pipeline.search.llm_parser import parse_natural_query
    from isimip_pipeline.search.isimip_query import ISIMIPQuery, SearchFilters
    from isimip_pipeline.search.result_table import ResultTable, export_selection

    cfg = get_config(config)

    console.print(Panel(f"[bold]Searching for:[/bold] {query}", title="ISIMIP Search"))

    # If refresh flag is set, prepend instruction to query
    search_query = query
    if refresh:
        search_query = f"search for: {query}"
        console.print("[dim]Refresh mode: forcing web search for definitions...[/dim]")

    # Parse query
    if use_llm and cfg.api.you_api_key:
        console.print("[dim]Parsing query with LLM...[/dim]")
        parsed = parse_natural_query(
            search_query,
            api_key=cfg.api.you_api_key,
            agent_id=cfg.api.you_agent_id or "",
        )
        filters = parsed.filters
        if parsed.explanation:
            console.print(f"[dim]Interpretation: {parsed.explanation}[/dim]")
        # Show if web search was used
        if parsed.raw_response and parsed.raw_response.get("web_search_used"):
            reason = parsed.raw_response.get("web_search_reason", "")
            console.print(f"[dim]Web search used: {reason}[/dim]")
    else:
        # Use keyword fallback
        from isimip_pipeline.search.llm_parser import LLMParser
        parser = LLMParser(api_key="", agent_id="")
        parsed = parser.keyword_fallback(query)
        filters = parsed.filters
        console.print(f"[dim]Extracted filters: {parsed.explanation}[/dim]")

    # Search ISIMIP
    console.print("[dim]Querying ISIMIP repository...[/dim]")

    try:
        isimip = ISIMIPQuery(timeout=cfg.api.isimip_timeout)

        # If no filters extracted, try direct text search
        if not any([filters.variable, filters.simulation_round, filters.climate_scenario]):
            datasets = isimip.search_by_query(query)
        else:
            datasets = isimip.search(filters)

        # Limit results
        if len(datasets) > limit:
            console.print(f"[yellow]Showing first {limit} of {len(datasets)} results[/yellow]")
            datasets = datasets[:limit]

        # Display results
        table = ResultTable(datasets)
        table.display(console)

        # Export if requested
        if output:
            export_selection(datasets, output, query=query)
            console.print(f"\n[green]Results saved to:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error searching ISIMIP: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def find(
    query: Optional[str] = typer.Argument(
        None, help="Search query (optional)"
    ),
    variable: Optional[str] = typer.Option(
        None, "--variable", "-v", help="Filter by variable (e.g., 'led')"
    ),
    timestep: Optional[str] = typer.Option(
        None, "--timestep", "-t", help="Filter by timestep (e.g., 'monthly')"
    ),
    scenario: Optional[str] = typer.Option(
        None, "--scenario", "-s", help="Filter by scenario (e.g., 'ssp126')"
    ),
    detailed: bool = typer.Option(
        False, "--detailed", "-d", help="Show detailed information"
    ),
):
    """Search locally processed datasets.

    Searches your local database of processed datasets stored in
    outputs/processed_data_log.yaml

    Examples:
        isimip-pipeline find                    # List all datasets
        isimip-pipeline find drought            # Search by query
        isimip-pipeline find -v led -t monthly  # Filter by variable and timestep
        isimip-pipeline find drought -d         # Detailed view
    """
    # Load processing log
    log_path = Path("./outputs/processed_data_log.yaml")

    try:
        log = load_processing_log(log_path)

        if not log.datasets:
            console.print("[yellow]No datasets found in local log.[/yellow]")
            console.print(
                "[dim]Run 'isimip-pipeline run' or 'isimip-pipeline interactive' "
                "to process and store datasets.[/dim]"
            )
            return

        # Search with filters
        results = find_local_datasets(
            log,
            query=query,
            variable=variable,
            timestep=timestep,
            scenario=scenario,
        )

        if not results:
            console.print("[yellow]No datasets match your search criteria.[/yellow]")
            return

        # Display results
        display_local_results(results, console, detailed=detailed)

    except Exception as e:
        console.print(f"[red]Error searching local datasets: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def interactive(
    query: str = typer.Argument(..., help="Natural language search query"),
    local_only: bool = typer.Option(
        False, "--local-only", help="Only search local datasets"
    ),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
):
    """Run interactive dataset discovery and selection workflow.

    Guides you through the workflow:
    1. Search local datasets
    2. Search ISIMIP if not found locally
    3. Select dataset group
    4. Check for duplicates

    After selection, you'll be guided to download and process the data.

    Examples:
        isimip-pipeline interactive "drought exposure"
        isimip-pipeline interactive "wildfire" --local-only
    """
    from isimip_pipeline.search.llm_parser import parse_natural_query, LLMParser
    from isimip_pipeline.search.isimip_query import ISIMIPQuery
    from rich.prompt import Prompt

    cfg = get_config(config)

    console.print(Panel(
        f"[bold]Interactive Dataset Workflow[/bold]\nQuery: {query}",
        title="ISIMIP Pipeline"
    ))

    # =========================================================================
    # STEP 1: LOCAL SEARCH
    # =========================================================================
    console.print("\n[bold cyan]Step 1/4: Searching local datasets...[/bold cyan]")

    log_path = Path("./outputs/processed_data_log.yaml")
    log = load_processing_log(log_path)

    local_results = find_local_datasets(log, query=query)

    if local_results:
        console.print(f"[green]Found {len(local_results)} local dataset(s)[/green]")
        display_local_results(local_results, console, detailed=True)

        use_existing = Prompt.ask(
            "[bold]Use one of these existing datasets?[/bold]",
            choices=["y", "n"],
            default="n",
        )

        if use_existing.lower() == "y":
            console.print("[green]Using existing dataset[/green]")
            if local_results:
                console.print(f"[dim]Dataset: {local_results[0].output_path}[/dim]")
            return

    # =========================================================================
    # STEP 2: REMOTE SEARCH
    # =========================================================================
    if not local_only:
        console.print("\n[bold cyan]Step 2/4: Searching ISIMIP repository...[/bold cyan]")

        try:
            # Parse query
            if cfg.api.you_api_key:
                parsed = parse_natural_query(
                    query,
                    api_key=cfg.api.you_api_key,
                    agent_id=cfg.api.you_agent_id or "",
                )
            else:
                parser = LLMParser(api_key="", agent_id="")
                parsed = parser.keyword_fallback(query)

            filters = parsed.filters
            if parsed.explanation:
                console.print(f"[dim]Interpretation: {parsed.explanation}[/dim]")

            # Search ISIMIP
            isimip = ISIMIPQuery(timeout=cfg.api.isimip_timeout)

            if not any([filters.variable, filters.simulation_round, filters.climate_scenario]):
                datasets = isimip.search_by_query(query)
            else:
                datasets = isimip.search(filters)

            if not datasets:
                console.print("[red]No datasets found in ISIMIP[/red]")
                raise typer.Exit(1)

            console.print(f"[green]Found {len(datasets)} datasets[/green]")

            # Group by variable+timestep
            grouped = group_by_variable_timestep(datasets)

            if not grouped:
                console.print("[red]No dataset groups found[/red]")
                raise typer.Exit(1)

            # =========================================================================
            # STEP 3: USER SELECTION
            # =========================================================================
            console.print("\n[bold cyan]Step 3/4: Select dataset group[/bold cyan]")

            display_grouped_results(grouped, console)

            # Let user pick a group
            group_choice = Prompt.ask(
                "[bold]Select group number[/bold]",
                choices=[str(i) for i in range(1, len(grouped) + 1)],
            )

            group_idx = int(group_choice) - 1
            selected_group_key = sorted(grouped.keys())[group_idx]
            selected_group = grouped[selected_group_key][
                "datasets"
            ]
            variable, timestep = selected_group_key

            console.print(
                f"[green]Selected: {variable} ({timestep})[/green]"
            )
            console.print(f"[dim]Files: {len(selected_group)}[/dim]")

            # Get descriptive name
            default_name = query.lower().replace(" ", "-")[:30]
            descriptive_name = Prompt.ask(
                "[bold]Enter descriptive name[/bold]",
                default=default_name,
            )

            # =========================================================================
            # STEP 4: DUPLICATE CHECK
            # =========================================================================
            console.print("\n[bold cyan]Step 4/4: Checking for duplicates...[/bold cyan]")

            log = load_processing_log(log_path)
            existing = check_for_duplicate(log, variable, timestep)

            if existing:
                console.print(
                    f"[yellow]⚠️  Duplicate detected[/yellow]\n"
                    f"Dataset '{variable}-{timestep}' already exists:\n"
                    f"  Name: {existing.descriptive_name}\n"
                    f"  Created: {existing.created_date.strftime('%Y-%m-%d')}\n"
                    f"  Path: {existing.output_path}"
                )

                action = Prompt.ask(
                    "[bold]What would you like to do?[/bold]",
                    choices=["skip", "new", "overwrite", "abort"],
                    default="skip",
                )

                if action == "skip":
                    console.print("[dim]Using existing dataset[/dim]")
                    return
                elif action == "abort":
                    console.print("[yellow]Operation cancelled[/yellow]")
                    return
                elif action == "overwrite":
                    from isimip_pipeline.duplicate_handler import generate_unique_name
                    folder_name = build_output_folder_name(
                        descriptive_name, variable, timestep
                    )
                else:  # new
                    from isimip_pipeline.duplicate_handler import generate_unique_name
                    folder_name = generate_unique_name(
                        descriptive_name, variable, timestep, log
                    )
            else:
                console.print("[green]No duplicates found[/green]")
                folder_name = build_output_folder_name(
                    descriptive_name, variable, timestep
                )

            # Save selection and all datasets
            output_dir = Path("./outputs") / folder_name
            output_dir.mkdir(parents=True, exist_ok=True)

            save_selection_metadata(
                output_dir,
                selected_group,
                query=query,
                descriptive_name=descriptive_name,
            )
            save_all_available_datasets(output_dir, datasets, query=query)

            console.print(f"\n[bold green]Ready to process![/bold green]")
            console.print(f"Output directory: [cyan]{output_dir}[/cyan]")
            console.print(
                f"\n[bold]Next steps:[/bold]\n"
                f"1. Download: isimip-pipeline download -s {output_dir}/selection.json\n"
                f"2. Process: isimip-pipeline process {output_dir}/raw"
            )

        except Exception as e:
            console.print(f"[red]Error during interactive workflow: {e}[/red]")
            raise typer.Exit(1)
    else:
        console.print("[yellow]Local-only mode, no remote search[/yellow]")


@app.command()
def download(
    selection: Path = typer.Option(
        ..., "--selection", "-s", help="JSON file with dataset selection"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output directory for downloads"
    ),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
    max_concurrent: int = typer.Option(
        4, "--concurrent", help="Maximum concurrent downloads"
    ),
):
    """Download selected datasets from ISIMIP repository.

    Examples:
        isimip-pipeline download --selection results.json
        isimip-pipeline download -s results.json -o ./data/raw
    """
    from isimip_pipeline.download.downloader import Downloader, DownloadStatus
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

    cfg = get_config(config)

    # Load selection
    if not selection.exists():
        console.print(f"[red]Selection file not found:[/red] {selection}")
        raise typer.Exit(1)

    with open(selection) as f:
        data = json.load(f)

    datasets = data.get("datasets", [])
    urls = [ds["url"] for ds in datasets if ds.get("url")]

    if not urls:
        console.print("[yellow]No URLs found in selection file[/yellow]")
        raise typer.Exit(1)

    console.print(Panel(
        f"[bold]Downloading {len(urls)} files[/bold]",
        title="ISIMIP Download"
    ))

    # Set up downloader
    output_dir = output or cfg.paths.download_dir
    downloader = Downloader(
        output_dir=output_dir,
        max_concurrent=max_concurrent,
    )

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download with progress
    console.print(f"[dim]Output directory: {output_dir}[/dim]")

    results = downloader.download_sync(urls)

    # Report results
    success = sum(1 for r in results if r.status == DownloadStatus.SUCCESS)
    skipped = sum(1 for r in results if r.status == DownloadStatus.SKIPPED)
    failed = sum(1 for r in results if r.status == DownloadStatus.FAILED)

    console.print(f"\n[green]Downloaded:[/green] {success}")
    console.print(f"[yellow]Skipped (existing):[/yellow] {skipped}")
    if failed:
        console.print(f"[red]Failed:[/red] {failed}")
        for r in results:
            if r.status == DownloadStatus.FAILED:
                console.print(f"  [red]- {r.url}: {r.error}[/red]")


@app.command()
def process(
    input_dir: Path = typer.Argument(..., help="Directory with raw NetCDF files"),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output directory for processed files"
    ),
    variable: Optional[str] = typer.Option(
        None, "--variable", "-v", help="Variable to extract (auto-detected if not specified)"
    ),
    scenarios: Optional[str] = typer.Option(
        None, "--scenarios", "-s", help="Comma-separated scenarios (e.g., ssp126,ssp585)"
    ),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
):
    """Process downloaded NetCDF files to extract features.

    Extracts 6 feature types per decade:
    1. Smoothed median value
    2. Percentile rank (1-100)
    3. Decadal trend (Theil-Sen slope)
    4. Trend significance (Spearman p-value)
    5. Lower confidence bound
    6. Upper confidence bound

    Auto-detects variable, timestep, and descriptive name from folder structure.
    Updates processing log after successful completion.

    Examples:
        isimip-pipeline process ./outputs/drought-severity_led-monthly/raw
        isimip-pipeline process ./data/raw --output ./data/processed
        isimip-pipeline process ./data/raw -v burntarea -s ssp126,ssp585
    """
    from datetime import datetime
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.prompt import Prompt
    from isimip_pipeline.processing.processor import (
        DataProcessor,
        find_netcdf_files,
        group_files_by_variable,
        detect_timestep_from_files,
    )
    from isimip_pipeline.processing.output import write_netcdf
    from isimip_pipeline.interactive import parse_descriptive_name_from_folder, load_selection_metadata

    cfg = get_config(config)

    console.print(Panel(
        f"[bold]Processing files from:[/bold] {input_dir}",
        title="ISIMIP Process"
    ))

    # Validate input directory
    if not input_dir.exists():
        console.print(f"[red]Input directory not found:[/red] {input_dir}")
        raise typer.Exit(1)

    # Find NetCDF files
    files = find_netcdf_files(input_dir)
    if not files:
        console.print(f"[red]No NetCDF files found in:[/red] {input_dir}")
        raise typer.Exit(1)

    console.print(f"[dim]Found {len(files)} NetCDF files[/dim]")

    # Group files by variable
    groups = group_files_by_variable(files)
    console.print(f"[dim]Variables detected: {', '.join(groups.keys())}[/dim]")

    # Auto-detect variable if not specified
    if variable:
        if variable not in groups:
            console.print(f"[red]Variable '{variable}' not found in files[/red]")
            console.print(f"[dim]Available: {', '.join(groups.keys())}[/dim]")
            raise typer.Exit(1)
        variables_to_process = [variable]
        detected_variable = variable
    else:
        variables_to_process = list(groups.keys())
        detected_variable = variables_to_process[0] if variables_to_process else "unknown"

    # Auto-detect timestep from files
    detected_timestep = detect_timestep_from_files(files)
    console.print(f"[dim]Detected timestep: {detected_timestep}[/dim]")

    # Auto-detect descriptive name from folder structure
    input_parent = input_dir.parent
    folder_name = input_parent.name
    detected_name = parse_descriptive_name_from_folder(folder_name)
    console.print(f"[dim]Detected name: {detected_name}[/dim]")

    # Try to load selection metadata for query info
    detected_query = ""
    try:
        selection = load_selection_metadata(input_parent)
        detected_query = selection.get("query", "")
    except FileNotFoundError:
        pass

    # Parse scenarios
    scenario_list = None
    if scenarios:
        scenario_list = [s.strip() for s in scenarios.split(",")]
        console.print(f"[dim]Scenarios: {', '.join(scenario_list)}[/dim]")

    # Display detected metadata for user confirmation
    console.print(f"\n[bold]Detected metadata:[/bold]")
    console.print(f"  Variable: [cyan]{detected_variable}[/cyan]")
    console.print(f"  Timestep: [cyan]{detected_timestep}[/cyan]")
    console.print(f"  Name: [cyan]{detected_name}[/cyan]")
    console.print(f"  Files: [cyan]{len(files)}[/cyan]")

    # Confirm with user before processing
    proceed = Prompt.ask(
        "[bold]Proceed with these settings?[/bold]",
        choices=["y", "n"],
        default="y",
    )

    if proceed.lower() != "y":
        console.print("[yellow]Processing cancelled[/yellow]")
        raise typer.Exit(0)

    # Set up output directory
    output_dir = output or cfg.paths.processed_dir
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize processor with config
    processor = DataProcessor(
        bandwidth=cfg.processing.smoothing_bandwidth,
        percentile_bins=cfg.processing.percentile_bins,
    )

    # Track processed files for logging
    processed_files = []
    processed_successfully = True

    # Process each variable
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        for var in variables_to_process:
            task = progress.add_task(f"Processing {var}...", total=None)

            try:
                result = processor.process(
                    input_dir=input_dir,
                    variable=var,
                    scenarios=scenario_list,
                )

                # Save output
                output_path = output_dir / f"{var}_processed.nc"
                write_netcdf(result, output_path)
                processed_files.append(str(output_path))

                progress.update(task, description=f"[green]Completed {var}[/green]")
                console.print(f"  [green]Saved:[/green] {output_path}")

            except Exception as e:
                processed_successfully = False
                progress.update(task, description=f"[red]Failed {var}[/red]")
                console.print(f"  [red]Error processing {var}: {e}[/red]")

    if not processed_successfully:
        console.print(f"\n[yellow]Processing completed with errors[/yellow]")
    else:
        console.print(f"\n[green]Processing complete![/green]")

    console.print(f"[dim]Output directory: {output_dir}[/dim]")

    # Update processing log if processing was successful
    if processed_successfully and processed_files:
        console.print("\n[dim]Updating processing log...[/dim]")

        try:
            from isimip_pipeline.processing_log import (
                load_processing_log,
                save_processing_log,
                DatasetEntry,
            )

            log_path = Path("./outputs/processed_data_log.yaml")
            log = load_processing_log(log_path)

            # Extract scenarios from processed data or use defaults
            climate_scenarios = scenario_list or ["unknown"]

            # Create entry for the log
            entry = DatasetEntry(
                descriptive_name=detected_name,
                variable=detected_variable,
                timestep=detected_timestep,
                created_date=datetime.now(),
                output_path=str(input_parent),
                file_count=len(files),
                time_periods=["2006-2100"],  # Could be extracted from data
                climate_scenarios=climate_scenarios,
                gcm_models=["unknown"],  # Could be extracted from filenames
                lsm_models=["unknown"],  # Could be extracted from filenames
                simulation_round="unknown",  # Could be extracted from data
                query=detected_query,
            )

            log.add_entry(entry)
            save_processing_log(log, log_path)

            console.print("[green]Processing log updated[/green]")

        except Exception as e:
            console.print(f"[yellow]Warning: Could not update processing log: {e}[/yellow]")


@app.command()
def report(
    processed_dir: Path = typer.Argument(
        ..., help="Directory with processed NetCDF files"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output path for HTML report"
    ),
    title: str = typer.Option(
        "ISIMIP QA Report", "--title", "-t", help="Report title"
    ),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
):
    """Generate interactive QA report from processed data.

    Creates an HTML report with:
    - Global heatmaps for each decade/scenario
    - Summary statistics
    - Data coverage information

    Examples:
        isimip-pipeline report ./data/processed
        isimip-pipeline report ./data/processed -o ./reports/qa.html
    """
    import xarray as xr
    from isimip_pipeline.visualization.qa_report import generate_html_report

    cfg = get_config(config)

    console.print(Panel(
        f"[bold]Generating report from:[/bold] {processed_dir}",
        title="ISIMIP Report"
    ))

    # Validate input directory
    processed_dir = Path(processed_dir)
    if not processed_dir.exists():
        console.print(f"[red]Directory not found:[/red] {processed_dir}")
        raise typer.Exit(1)

    # Find processed NetCDF files
    nc_files = list(processed_dir.glob("*_processed.nc"))
    if not nc_files:
        console.print(f"[red]No processed NetCDF files found in:[/red] {processed_dir}")
        console.print("[dim]Expected files matching *_processed.nc[/dim]")
        raise typer.Exit(1)

    console.print(f"[dim]Found {len(nc_files)} processed files[/dim]")

    # Set up output path
    reports_dir = cfg.paths.reports_dir
    if output:
        output_path = Path(output)
    else:
        output_path = reports_dir / "qa_report.html"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Generate report for each file
    all_html_parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        f"<title>{title}</title>",
        '<meta charset="utf-8">',
        '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>',
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 20px; }",
        "h1 { color: #333; }",
        "h2 { color: #666; margin-top: 30px; }",
        ".variable-section { margin-bottom: 40px; border-bottom: 1px solid #ddd; padding-bottom: 20px; }",
        "</style>",
        "</head>",
        "<body>",
        f"<h1>{title}</h1>",
    ]

    for nc_file in nc_files:
        variable = nc_file.stem.replace("_processed", "")
        console.print(f"[dim]Processing {variable}...[/dim]")

        try:
            ds = xr.open_dataset(nc_file)

            # Generate HTML for this variable
            html = generate_html_report(
                ds,
                variable=variable,
                title=f"{variable} Analysis",
                include_maps=True,
                include_summary=True,
            )

            # Extract body content (skip full HTML wrapper)
            body_start = html.find("<body>") + 6
            body_end = html.find("</body>")
            body_content = html[body_start:body_end]

            all_html_parts.append(f'<div class="variable-section">')
            all_html_parts.append(body_content)
            all_html_parts.append("</div>")

            ds.close()

        except Exception as e:
            console.print(f"  [red]Error processing {variable}: {e}[/red]")

    all_html_parts.extend(["</body>", "</html>"])

    # Write combined report
    with open(output_path, "w") as f:
        f.write("\n".join(all_html_parts))

    console.print(f"\n[green]Report generated:[/green] {output_path}")


@app.command()
def catalog(
    show_all: bool = typer.Option(
        False, "--all", "-a", help="Show all variable details"
    ),
    variable: Optional[str] = typer.Option(
        None, "--variable", "-v", help="Show details for specific variable"
    ),
    reset: bool = typer.Option(
        False, "--reset", help="Reset catalog to empty"
    ),
):
    """View or manage the persistent ISIMIP metrics catalog.

    The catalog tracks all variables, scenarios, and models encountered
    during search operations, building a local knowledge base over time.

    Examples:
        isimip-pipeline catalog           # Show summary
        isimip-pipeline catalog --all     # Show all variables
        isimip-pipeline catalog -v led    # Show details for 'led'
        isimip-pipeline catalog --reset   # Clear catalog
    """
    from rich.table import Table
    from isimip_pipeline.catalog import load_catalog, save_catalog, ISIMIPCatalog

    catalog_data = load_catalog()

    if reset:
        save_catalog(ISIMIPCatalog())
        console.print("[yellow]Catalog reset to empty[/yellow]")
        return

    if len(catalog_data.variables) == 0:
        console.print("[yellow]Catalog is empty. Run a search to populate it.[/yellow]")
        console.print("[dim]Example: isimip-pipeline run \"drought exposure\" --limit 5[/dim]")
        return

    summary = catalog_data.get_summary()

    if variable:
        # Show specific variable
        if variable not in catalog_data.variables:
            console.print(f"[red]Variable '{variable}' not found in catalog[/red]")
            console.print(f"[dim]Available: {', '.join(catalog_data.variables.keys())}[/dim]")
            return

        var_info = catalog_data.variables[variable]
        console.print(f"\n[bold]{variable}[/bold]")
        if var_info.long_name:
            console.print(f"  Name: {var_info.long_name}")
        if var_info.unit:
            console.print(f"  Unit: {var_info.unit}")
        console.print(f"  Files: {var_info.file_count}")
        console.print(f"  Scenarios: {', '.join(sorted(var_info.scenarios)) or 'None'}")
        console.print(f"  Models: {', '.join(sorted(var_info.models)) or 'None'}")
        console.print(f"  Simulation rounds: {', '.join(sorted(var_info.simulation_rounds)) or 'None'}")
        if var_info.last_seen:
            console.print(f"  Last seen: {var_info.last_seen.strftime('%Y-%m-%d %H:%M')}")
        return

    # Show summary
    console.print("\n[bold]ISIMIP Metrics Catalog[/bold]")
    console.print(f"  Variables tracked: {summary['total_variables']}")
    console.print(f"  Total files seen: {summary['total_files']}")
    console.print(f"  Scenarios: {', '.join(summary['all_scenarios']) or 'None'}")
    if summary['last_updated']:
        console.print(f"  Last updated: {summary['last_updated']}")

    if show_all:
        # Show table of all variables
        console.print("")
        table = Table(title="Variables")
        table.add_column("Variable", style="cyan")
        table.add_column("Long Name")
        table.add_column("Files", justify="right")
        table.add_column("Scenarios")
        table.add_column("Models")

        for name, var_info in sorted(catalog_data.variables.items()):
            table.add_row(
                name,
                var_info.long_name or "",
                str(var_info.file_count),
                ", ".join(sorted(var_info.scenarios)[:3]) + ("..." if len(var_info.scenarios) > 3 else ""),
                ", ".join(sorted(var_info.models)[:2]) + ("..." if len(var_info.models) > 2 else ""),
            )

        console.print(table)
    else:
        console.print(f"\n[dim]Variables: {', '.join(sorted(catalog_data.variables.keys()))}[/dim]")
        console.print("[dim]Use --all for details or -v <variable> for specific info[/dim]")


@app.command()
def run(
    query: str = typer.Argument(..., help="Natural language search query"),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Descriptive name for output folder (e.g., 'drought-severity')"
    ),
    scenarios: Optional[str] = typer.Option(
        None, "--scenarios", "-s", help="Comma-separated scenarios (e.g., ssp126,ssp370,ssp585)"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Override output directory (ignores --name)"
    ),
    limit: int = typer.Option(
        20, "--limit", "-l", help="Maximum datasets to download"
    ),
    skip_download: bool = typer.Option(
        False, "--skip-download", help="Skip download step (use existing files)"
    ),
    skip_process: bool = typer.Option(
        False, "--skip-process", help="Skip processing step"
    ),
    skip_report: bool = typer.Option(
        False, "--skip-report", help="Skip report generation"
    ),
    keep_raw: bool = typer.Option(
        False, "--keep-raw", help="Keep raw downloaded files (default: delete after processing)"
    ),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Path to config file"
    ),
):
    """Run complete pipeline: search, download, process, and report.

    This command executes all pipeline steps in sequence:
    1. Search ISIMIP repository for matching datasets
    2. Download NetCDF files to temp folder (~/Downloads/isimip_*)
    3. Process files to extract statistical features
    4. Generate an interactive HTML QA report
    5. Clean up temp files (unless --keep-raw is specified)

    Output is saved to ./outputs/{name}_{variable}/ by default.
    Use --name to specify a descriptive name (e.g., 'drought-severity').

    Examples:
        isimip-pipeline run "drought exposure" --name drought-severity
        isimip-pipeline run "wildfire burnt area" -n fire-risk
        isimip-pipeline run "drought" --scenarios ssp126,ssp585
        isimip-pipeline run "flood" -o ./custom_output  # Override path
        isimip-pipeline run "drought" --keep-raw  # Preserve raw files
    """
    import tempfile
    from datetime import datetime
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

    from isimip_pipeline.search.llm_parser import parse_natural_query, LLMParser
    from isimip_pipeline.search.isimip_query import ISIMIPQuery
    from isimip_pipeline.search.result_table import ResultTable, export_selection
    from isimip_pipeline.download.downloader import Downloader, DownloadStatus
    from isimip_pipeline.processing.processor import DataProcessor, find_netcdf_files, group_files_by_variable
    from isimip_pipeline.processing.output import write_netcdf
    from isimip_pipeline.visualization.qa_report import generate_html_report
    import xarray as xr

    cfg = get_config(config)

    # Output directory will be determined after search (needs variable name)
    # For now, just show the query
    console.print(Panel(
        f"[bold]Running full pipeline[/bold]\n"
        f"Query: {query}",
        title="ISIMIP Pipeline"
    ))

    # Parse scenarios
    scenario_list = None
    if scenarios:
        scenario_list = [s.strip() for s in scenarios.split(",")]

    # =========================================================================
    # STEP 1: SEARCH
    # =========================================================================
    console.print("\n[bold cyan]Step 1/4: Searching ISIMIP repository...[/bold cyan]")

    # Parse query
    if cfg.api.you_api_key:
        parsed = parse_natural_query(
            query,
            api_key=cfg.api.you_api_key,
            agent_id=cfg.api.you_agent_id or "",
        )
    else:
        parser = LLMParser(api_key="", agent_id="")
        parsed = parser.keyword_fallback(query)

    filters = parsed.filters
    if parsed.explanation:
        console.print(f"[dim]Interpretation: {parsed.explanation}[/dim]")

    # Apply scenario filter if specified
    if scenario_list and len(scenario_list) == 1:
        filters.climate_scenario = scenario_list[0]

    # Search ISIMIP
    try:
        isimip = ISIMIPQuery(timeout=cfg.api.isimip_timeout)

        if not any([filters.variable, filters.simulation_round, filters.climate_scenario]):
            datasets = isimip.search_by_query(query)
        else:
            datasets = isimip.search(filters)

        # Filter by scenarios if multiple specified
        if scenario_list and len(scenario_list) > 1:
            datasets = [d for d in datasets if d.climate_scenario in scenario_list]

        if not datasets:
            console.print("[red]No datasets found matching query[/red]")
            raise typer.Exit(1)

        # Extract variable and timestep from datasets
        variables_found = set(ds.variable for ds in datasets if ds.variable)
        variable_name = list(variables_found)[0] if variables_found else "unknown"

        timesteps_found = set(ds.timestep for ds in datasets if ds.timestep)
        timestep_name = list(timesteps_found)[0] if timesteps_found else "unknown"

        # Set up output directory structure: ./outputs/{name}_{variable}-{timestep}/
        if output:
            # User specified custom output path
            base_dir = Path(output)
        else:
            # Use default naming convention with timestep
            descriptive_name = name if name else query.lower().replace(" ", "-")[:30]
            # Use build_output_folder_name for consistent naming
            folder_name = build_output_folder_name(descriptive_name, variable_name, timestep_name)
            base_dir = Path("./outputs") / folder_name

        raw_dir = base_dir / "raw"
        processed_dir = base_dir / "processed"
        reports_dir = base_dir / "reports"

        console.print(f"[dim]Output directory: {base_dir}[/dim]")

        # Save ALL available datasets before limiting
        base_dir.mkdir(parents=True, exist_ok=True)
        all_datasets_file = base_dir / "all_available_datasets.json"
        export_selection(datasets, all_datasets_file, query=query)
        console.print(f"[dim]Saved all {len(datasets)} available datasets to {all_datasets_file.name}[/dim]")

        # Update persistent ISIMIP catalog
        from isimip_pipeline.catalog import update_catalog_from_datasets
        catalog = update_catalog_from_datasets(datasets)
        console.print(f"[dim]Updated ISIMIP catalog ({len(catalog.variables)} variables tracked)[/dim]")

        # Limit results for download
        if len(datasets) > limit:
            console.print(f"[yellow]Limiting to first {limit} of {len(datasets)} datasets for download[/yellow]")
            datasets = datasets[:limit]

        # Display summary
        table = ResultTable(datasets)
        summary = table.get_summary()
        console.print(f"[green]Selected {len(datasets)} datasets for download[/green]")
        console.print(f"[dim]Variables: {', '.join(summary['variables'])}[/dim]")
        console.print(f"[dim]Scenarios: {', '.join(summary['scenarios'])}[/dim]")

        # Save selection for download
        selection_file = base_dir / "selection.json"
        export_selection(datasets, selection_file, query=query)

    except Exception as e:
        console.print(f"[red]Search failed: {e}[/red]")
        raise typer.Exit(1)

    # =========================================================================
    # STEP 2: DOWNLOAD
    # =========================================================================
    # Track downloader for cleanup
    downloader: Optional[Downloader] = None

    if not skip_download:
        console.print("\n[bold cyan]Step 2/4: Downloading datasets...[/bold cyan]")

        urls = [ds.url for ds in datasets if ds.url]
        if not urls:
            console.print("[red]No download URLs found[/red]")
            raise typer.Exit(1)

        # Use temp directory in Downloads folder (will be cleaned up after processing)
        temp_base = Path.home() / "Downloads"
        downloader = Downloader(
            use_temp=not keep_raw,
            temp_base=temp_base,
            output_dir=raw_dir if keep_raw else None,
            max_concurrent=cfg.download.max_concurrent,
        )

        # Update raw_dir to point to actual download location
        raw_dir = downloader.output_dir

        if downloader.use_temp:
            console.print(f"[dim]Downloading {len(urls)} files to temp: {raw_dir}[/dim]")
            console.print(f"[dim](Raw files will be deleted after processing)[/dim]")
        else:
            raw_dir.mkdir(parents=True, exist_ok=True)
            console.print(f"[dim]Downloading {len(urls)} files to {raw_dir}[/dim]")

        results = downloader.download_sync(urls)

        success = sum(1 for r in results if r.status == DownloadStatus.SUCCESS)
        skipped = sum(1 for r in results if r.status == DownloadStatus.SKIPPED)
        failed = sum(1 for r in results if r.status == DownloadStatus.FAILED)

        console.print(f"[green]Downloaded: {success}[/green]")
        if skipped:
            console.print(f"[yellow]Skipped (existing): {skipped}[/yellow]")
        if failed:
            console.print(f"[red]Failed: {failed}[/red]")
            for r in results:
                if r.status == DownloadStatus.FAILED:
                    console.print(f"  [red]- {r.error}[/red]")
    else:
        console.print("\n[bold cyan]Step 2/4: Skipping download (--skip-download)[/bold cyan]")
        if not raw_dir.exists():
            console.print(f"[red]Raw directory not found: {raw_dir}[/red]")
            raise typer.Exit(1)

    # =========================================================================
    # STEP 3: PROCESS
    # =========================================================================
    if not skip_process:
        console.print("\n[bold cyan]Step 3/4: Processing data...[/bold cyan]")

        files = find_netcdf_files(raw_dir)
        if not files:
            console.print(f"[red]No NetCDF files found in {raw_dir}[/red]")
            raise typer.Exit(1)

        groups = group_files_by_variable(files)
        console.print(f"[dim]Found {len(files)} files, variables: {', '.join(groups.keys())}[/dim]")

        processed_dir.mkdir(parents=True, exist_ok=True)
        processor = DataProcessor(
            bandwidth=cfg.processing.smoothing_bandwidth,
            percentile_bins=cfg.processing.percentile_bins,
        )

        for var in groups.keys():
            console.print(f"[dim]Processing {var}...[/dim]")
            try:
                result = processor.process(
                    input_dir=raw_dir,
                    variable=var,
                    scenarios=scenario_list,
                )

                output_path = processed_dir / f"{var}_processed.nc"
                write_netcdf(result, output_path)
                console.print(f"  [green]Saved: {output_path.name}[/green]")

            except Exception as e:
                console.print(f"  [red]Error: {e}[/red]")
    else:
        console.print("\n[bold cyan]Step 3/4: Skipping processing (--skip-process)[/bold cyan]")

    # =========================================================================
    # STEP 4: REPORT
    # =========================================================================
    if not skip_report:
        console.print("\n[bold cyan]Step 4/4: Generating report...[/bold cyan]")

        nc_files = list(processed_dir.glob("*_processed.nc")) if processed_dir.exists() else []

        if not nc_files:
            console.print("[yellow]No processed files found, skipping report[/yellow]")
        else:
            reports_dir.mkdir(parents=True, exist_ok=True)
            report_path = reports_dir / "qa_report.html"

            all_html_parts = [
                "<!DOCTYPE html>",
                "<html>",
                "<head>",
                f"<title>ISIMIP Pipeline Report - {query}</title>",
                '<meta charset="utf-8">',
                '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>',
                "<style>",
                "body { font-family: Arial, sans-serif; margin: 20px; }",
                "h1 { color: #333; }",
                "h2 { color: #666; margin-top: 30px; }",
                ".summary { background: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px; }",
                ".variable-section { margin-bottom: 40px; border-bottom: 1px solid #ddd; padding-bottom: 20px; }",
                "</style>",
                "</head>",
                "<body>",
                f"<h1>ISIMIP Pipeline Report</h1>",
                f'<div class="summary">',
                f"<p><strong>Query:</strong> {query}</p>",
                f"<p><strong>Datasets:</strong> {len(datasets)}</p>",
                f"<p><strong>Generated:</strong> {datetime.now().isoformat()}</p>",
                "</div>",
            ]

            for nc_file in nc_files:
                variable = nc_file.stem.replace("_processed", "")
                console.print(f"[dim]Adding {variable} to report...[/dim]")

                try:
                    ds = xr.open_dataset(nc_file)
                    html = generate_html_report(ds, variable=variable, title=f"{variable} Analysis")

                    body_start = html.find("<body>") + 6
                    body_end = html.find("</body>")
                    body_content = html[body_start:body_end]

                    all_html_parts.append(f'<div class="variable-section">')
                    all_html_parts.append(body_content)
                    all_html_parts.append("</div>")
                    ds.close()

                except Exception as e:
                    console.print(f"  [red]Error: {e}[/red]")

            all_html_parts.extend(["</body>", "</html>"])

            with open(report_path, "w") as f:
                f.write("\n".join(all_html_parts))

            console.print(f"[green]Report saved: {report_path}[/green]")
    else:
        console.print("\n[bold cyan]Step 4/4: Skipping report (--skip-report)[/bold cyan]")

    # =========================================================================
    # CLEANUP
    # =========================================================================
    if downloader and downloader.use_temp and not keep_raw:
        console.print("\n[dim]Cleaning up temporary files...[/dim]")
        if downloader.cleanup():
            console.print("[dim]Temp files removed[/dim]")
        else:
            console.print("[yellow]Warning: Could not remove temp files[/yellow]")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    console.print("\n" + "=" * 60)
    console.print("[bold green]Pipeline complete![/bold green]")
    console.print(f"[dim]Output directory: {base_dir}[/dim]")
    console.print(f"[dim]  - selection.json: Dataset metadata[/dim]")
    if not skip_download and keep_raw:
        console.print(f"[dim]  - raw/: Downloaded NetCDF files[/dim]")
    if not skip_process:
        console.print(f"[dim]  - processed/: Processed feature files[/dim]")
    if not skip_report:
        console.print(f"[dim]  - reports/: HTML QA report[/dim]")


if __name__ == "__main__":
    app()
