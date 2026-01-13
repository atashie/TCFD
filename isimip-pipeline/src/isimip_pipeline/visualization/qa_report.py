"""QA Report generation for processed ISIMIP data.

Creates interactive HTML reports with Plotly visualizations:
- Global maps for each decade/scenario/value_class
- Time series showing decadal progression
- Difference maps for future vs baseline
- Summary statistics
"""

from typing import Dict, List, Any, Optional

import numpy as np
import xarray as xr
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


def create_map_figure(
    ds: xr.Dataset,
    variable: str,
    decade: int,
    scenario: str,
    value_class: str,
    title: Optional[str] = None,
) -> go.Figure:
    """Create a global geographic map figure for a specific slice.

    Args:
        ds: xarray Dataset with processed data.
        variable: Variable name in dataset.
        decade: Decade to display (10, 20, ..., 90).
        scenario: Scenario to display (e.g., "ssp126").
        value_class: Value class to display.
        title: Optional figure title.

    Returns:
        Plotly Figure object.
    """
    # Extract the data slice
    data = ds[variable].sel(
        decade=decade,
        scenario=scenario,
        value_class=value_class,
    )

    # Get coordinates - data shape is (lon, lat)
    lon_vals = ds.lon.values
    lat_vals = ds.lat.values
    z_data = data.values  # Shape: (lon, lat)

    # Create mesh grid for coordinates
    lon_grid, lat_grid = np.meshgrid(lon_vals, lat_vals, indexing='ij')

    # Flatten arrays
    lon_flat = lon_grid.flatten()
    lat_flat = lat_grid.flatten()
    z_flat = z_data.flatten()

    # Remove NaN values
    valid_mask = ~np.isnan(z_flat)
    lon_valid = lon_flat[valid_mask]
    lat_valid = lat_flat[valid_mask]
    z_valid = z_flat[valid_mask]

    # Handle empty data
    if len(z_valid) == 0:
        fig = go.Figure()
        fig.update_layout(
            title=f"No data available for {variable} - {value_class}",
            height=300,
        )
        return fig

    # Convert to Python lists for Plotly JSON serialization
    lon_valid = lon_valid.tolist()
    lat_valid = lat_valid.tolist()
    z_valid = z_valid.tolist()

    # Choose colorscale based on value_class
    if value_class == "trend":
        colorscale = "RdBu_r"  # Diverging for trends
        zmid = 0
    elif value_class == "significance":
        colorscale = "Reds_r"  # Low p-values = significant
        zmid = None
    elif value_class == "percentile":
        colorscale = "Viridis"
        zmid = None
    else:
        colorscale = "Viridis"
        zmid = None

    # Create geographic scatter plot
    fig = go.Figure()

    fig.add_trace(
        go.Scattergeo(
            lon=lon_valid,
            lat=lat_valid,
            mode='markers',
            marker=dict(
                size=3,
                color=z_valid,
                colorscale=colorscale,
                showscale=True,
                colorbar=dict(title=value_class),
                cmin=np.nanpercentile(z_valid, 2) if z_valid else 0,
                cmax=np.nanpercentile(z_valid, 98) if z_valid else 1,
            ),
            hovertemplate=(
                "Lon: %{lon:.1f}<br>"
                "Lat: %{lat:.1f}<br>"
                f"{value_class}: " + "%{marker.color:.4f}<extra></extra>"
            ),
        )
    )

    # Set layout with geographic projection
    default_title = f"{variable} - {value_class} ({scenario}, decade {decade})"
    fig.update_layout(
        title=title or default_title,
        height=500,
        geo=dict(
            showland=True,
            landcolor="rgb(243, 243, 243)",
            showocean=True,
            oceancolor="rgb(204, 229, 255)",
            showcoastlines=True,
            coastlinecolor="rgb(100, 100, 100)",
            coastlinewidth=0.5,
            showlakes=True,
            lakecolor="rgb(204, 229, 255)",
            showcountries=True,
            countrycolor="rgb(180, 180, 180)",
            countrywidth=0.3,
            projection_type="equirectangular",
            lataxis_range=[-90, 90],
            lonaxis_range=[-180, 180],
        ),
        margin=dict(l=0, r=0, t=50, b=0),
    )

    return fig


def create_timeseries_figure(
    ds: xr.Dataset,
    variable: str,
    lat_idx: int,
    lon_idx: int,
    scenario: str,
    value_class: str = "smoothed_median",
) -> go.Figure:
    """Create a timeseries figure showing decadal progression.

    Args:
        ds: xarray Dataset with processed data.
        variable: Variable name in dataset.
        lat_idx: Latitude index.
        lon_idx: Longitude index.
        scenario: Scenario to display.
        value_class: Value class to display.

    Returns:
        Plotly Figure object.
    """
    # Extract data at location
    data = ds[variable].isel(lat=lat_idx, lon=lon_idx).sel(
        scenario=scenario,
        value_class=value_class,
    )

    decades = ds.decade.values
    values = data.values

    # Get location coordinates
    lat = float(ds.lat.isel(lat=lat_idx).values)
    lon = float(ds.lon.isel(lon=lon_idx).values)

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=decades,
            y=values,
            mode="lines+markers",
            name=scenario,
            line=dict(width=2),
            marker=dict(size=8),
        )
    )

    fig.update_layout(
        title=f"{variable} at ({lat:.1f}°, {lon:.1f}°) - {scenario}",
        xaxis_title="Decade",
        yaxis_title=value_class,
        height=400,
    )

    return fig


def generate_html_report(
    ds: xr.Dataset,
    variable: str,
    title: str = "ISIMIP QA Report",
    include_maps: bool = True,
    include_summary: bool = True,
) -> str:
    """Generate complete HTML report.

    Report structure:
    1. Summary statistics
    2. Baseline vs End-of-Century (2010s vs 2090s) for all value classes
    3. Scenario Comparison (2090s median for all scenarios)

    Args:
        ds: xarray Dataset with processed data.
        variable: Variable name in dataset.
        title: Report title.
        include_maps: Whether to include map visualizations.
        include_summary: Whether to include summary statistics.

    Returns:
        HTML string.
    """
    report = QAReport(ds, variable)

    # Friendly names for value classes
    value_class_names = {
        "smoothed_median": "Median Value",
        "percentile": "Percentile Rank",
        "trend": "Trend (Theil-Sen Slope)",
        "significance": "Trend Significance (p-value)",
        "lower_bound": "Lower Bound (25th percentile)",
        "upper_bound": "Upper Bound (75th percentile)",
    }

    html_parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        f"<title>{title}</title>",
        '<meta charset="utf-8">',
        '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>',
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 20px; }",
        "h1 { color: #333; }",
        "h2 { color: #555; margin-top: 40px; border-bottom: 2px solid #ddd; padding-bottom: 10px; }",
        "h3 { color: #666; margin-top: 30px; }",
        ".summary { background: #f5f5f5; padding: 15px; border-radius: 5px; }",
        ".map-container { margin: 20px 0; }",
        ".section { margin-bottom: 50px; }",
        ".comparison-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }",
        "@media (max-width: 1200px) { .comparison-grid { grid-template-columns: 1fr; } }",
        "</style>",
        "</head>",
        "<body>",
        f"<h1>{title}</h1>",
    ]

    # Add summary section
    if include_summary:
        summary = report.get_summary()
        html_parts.append('<div class="summary">')
        html_parts.append("<h2>Summary</h2>")
        html_parts.append("<ul>")
        for key, value in summary.items():
            html_parts.append(f"<li><strong>{key}:</strong> {value}</li>")
        html_parts.append("</ul>")
        html_parts.append("</div>")

    # Add maps using structured report
    if include_maps:
        sections = report.generate_structured_report()
        scenarios = [str(s) for s in ds.scenario.values]
        value_classes = [str(vc) for vc in ds.value_class.values]

        # Get metadata from structured report
        metadata = sections.get("metadata", {})
        baseline_decade = metadata.get("baseline_decade", 10)
        endcentury_decade = metadata.get("endcentury_decade", 90)
        primary_scenario = metadata.get("primary_scenario", scenarios[0])

        # Section 1: Baseline vs End-of-Century Comparison
        html_parts.append('<div class="section">')
        html_parts.append(f"<h2>Baseline vs End-of-Century Comparison ({primary_scenario})</h2>")
        html_parts.append(f"<p>Comparing {baseline_decade}0s (baseline) with {endcentury_decade}0s (end-of-century) projections for all statistical measures.</p>")

        for vc in value_classes:
            vc_name = value_class_names.get(vc, vc)
            html_parts.append(f"<h3>{vc_name}</h3>")
            html_parts.append('<div class="comparison-grid">')

            # Baseline decade map
            fig_baseline = sections["baseline_comparison"].get(f"{baseline_decade}0s_{vc}")
            if fig_baseline:
                html_parts.append('<div class="map-container">')
                html_parts.append(fig_baseline.to_html(full_html=False, include_plotlyjs=False))
                html_parts.append("</div>")

            # End-of-century decade map
            fig_endcentury = sections["baseline_comparison"].get(f"{endcentury_decade}0s_{vc}")
            if fig_endcentury:
                html_parts.append('<div class="map-container">')
                html_parts.append(fig_endcentury.to_html(full_html=False, include_plotlyjs=False))
                html_parts.append("</div>")

            html_parts.append("</div>")  # Close comparison-grid

        html_parts.append("</div>")  # Close section

        # Section 2: Scenario Comparison (only if multiple scenarios)
        if len(scenarios) > 1:
            html_parts.append('<div class="section">')
            html_parts.append(f"<h2>Scenario Comparison ({endcentury_decade}0s Median Values)</h2>")
            html_parts.append("<p>Comparing end-of-century median values across different climate scenarios.</p>")

            for scenario in scenarios:
                fig = sections["scenario_comparison"].get(f"{endcentury_decade}0s_{scenario}")
                if fig:
                    html_parts.append('<div class="map-container">')
                    html_parts.append(fig.to_html(full_html=False, include_plotlyjs=False))
                    html_parts.append("</div>")

            html_parts.append("</div>")  # Close section

    html_parts.extend(["</body>", "</html>"])

    return "\n".join(html_parts)


class QAReport:
    """Generate QA reports for processed ISIMIP data.

    Creates interactive visualizations for quality assurance
    and data exploration.
    """

    def __init__(
        self,
        ds: xr.Dataset,
        variable: str,
    ):
        """Initialize QA report generator.

        Args:
            ds: xarray Dataset with processed data.
            variable: Variable name to analyze.
        """
        self.ds = ds
        self.variable = variable

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics for the dataset.

        Returns:
            Dict with summary statistics.
        """
        data = self.ds[self.variable]

        # Calculate basic statistics
        valid_count = int(np.sum(~np.isnan(data.values)))
        total_count = int(data.size)

        summary = {
            "variable": self.variable,
            "total_cells": total_count,
            "total_valid_cells": valid_count,
            "valid_fraction": f"{valid_count / total_count:.1%}",
            "decades": list(self.ds.decade.values),
            "scenarios": list(self.ds.scenario.values),
        }

        # Add value statistics for smoothed_median
        if "smoothed_median" in self.ds.value_class.values:
            sm_data = data.sel(value_class="smoothed_median")
            valid_values = sm_data.values[~np.isnan(sm_data.values)]
            if len(valid_values) > 0:
                summary["smoothed_median_mean"] = f"{np.mean(valid_values):.4f}"
                summary["smoothed_median_std"] = f"{np.std(valid_values):.4f}"

        return summary

    def generate_maps(
        self,
        decades: Optional[List[int]] = None,
        scenarios: Optional[List[str]] = None,
        value_classes: Optional[List[str]] = None,
    ) -> Dict[str, go.Figure]:
        """Generate map figures for specified combinations.

        Args:
            decades: List of decades (default: all).
            scenarios: List of scenarios (default: all).
            value_classes: List of value classes (default: ["smoothed_median"]).

        Returns:
            Dict mapping description to Figure.
        """
        if decades is None:
            decades = [int(d) for d in self.ds.decade.values]
        if scenarios is None:
            scenarios = [str(s) for s in self.ds.scenario.values]
        if value_classes is None:
            value_classes = ["smoothed_median"]

        maps = {}

        for decade in decades:
            for scenario in scenarios:
                for vc in value_classes:
                    key = f"{scenario}_{decade}s_{vc}"
                    fig = create_map_figure(
                        self.ds,
                        self.variable,
                        decade=decade,
                        scenario=scenario,
                        value_class=vc,
                    )
                    maps[key] = fig

        return maps

    def generate_structured_report(self) -> Dict[str, Dict[str, go.Figure]]:
        """Generate maps organized by report sections.

        Returns a structured dict with two sections:
        1. 'baseline_comparison': 2010s vs 2090s for all value classes (first scenario)
        2. 'scenario_comparison': 2090s median for all scenarios

        Returns:
            Dict with section names mapping to dicts of figures.
        """
        # Friendly names for value classes
        value_class_names = {
            "smoothed_median": "Median Value",
            "percentile": "Percentile Rank",
            "trend": "Trend (Theil-Sen Slope)",
            "significance": "Trend Significance (p-value)",
            "lower_bound": "Lower Bound (25th percentile)",
            "upper_bound": "Upper Bound (75th percentile)",
        }

        scenarios = [str(s) for s in self.ds.scenario.values]
        value_classes = [str(vc) for vc in self.ds.value_class.values]
        decades = [int(d) for d in self.ds.decade.values]
        primary_scenario = scenarios[0]  # Use first scenario for baseline comparison

        # Determine baseline and end-of-century decades
        baseline_decade = decades[0] if decades else 10  # First available decade
        endcentury_decade = decades[-1] if decades else 90  # Last available decade

        sections = {
            "baseline_comparison": {},
            "scenario_comparison": {},
            "metadata": {
                "baseline_decade": baseline_decade,
                "endcentury_decade": endcentury_decade,
                "primary_scenario": primary_scenario,
            }
        }

        # Section 1: Baseline vs End-of-Century - all value classes
        for vc in value_classes:
            vc_name = value_class_names.get(vc, vc)

            # Baseline decade
            fig_baseline = create_map_figure(
                self.ds,
                self.variable,
                decade=baseline_decade,
                scenario=primary_scenario,
                value_class=vc,
                title=f"{vc_name} - {baseline_decade}0s ({primary_scenario})",
            )
            sections["baseline_comparison"][f"{baseline_decade}0s_{vc}"] = fig_baseline

            # End-of-century decade
            fig_endcentury = create_map_figure(
                self.ds,
                self.variable,
                decade=endcentury_decade,
                scenario=primary_scenario,
                value_class=vc,
                title=f"{vc_name} - {endcentury_decade}0s ({primary_scenario})",
            )
            sections["baseline_comparison"][f"{endcentury_decade}0s_{vc}"] = fig_endcentury

        # Section 2: Scenario comparison - end-of-century median for all scenarios
        for scenario in scenarios:
            fig = create_map_figure(
                self.ds,
                self.variable,
                decade=endcentury_decade,
                scenario=scenario,
                value_class="smoothed_median",
                title=f"Median Value - {endcentury_decade}0s ({scenario})",
            )
            sections["scenario_comparison"][f"{endcentury_decade}0s_{scenario}"] = fig

        return sections

    def generate_timeseries(
        self,
        locations: List[tuple],
        scenario: str = "ssp126",
    ) -> Dict[str, go.Figure]:
        """Generate timeseries figures for specified locations.

        Args:
            locations: List of (lat_idx, lon_idx) tuples.
            scenario: Scenario to display.

        Returns:
            Dict mapping location to Figure.
        """
        timeseries = {}

        for lat_idx, lon_idx in locations:
            lat = float(self.ds.lat.isel(lat=lat_idx).values)
            lon = float(self.ds.lon.isel(lon=lon_idx).values)
            key = f"({lat:.1f}, {lon:.1f})"

            fig = create_timeseries_figure(
                self.ds,
                self.variable,
                lat_idx=lat_idx,
                lon_idx=lon_idx,
                scenario=scenario,
            )
            timeseries[key] = fig

        return timeseries
