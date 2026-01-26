"""Generate HTML visualization reports for extraction results.

Creates interactive time series, comparison charts, risk heatmaps,
and data tables from CSV extraction output.

Usage:
    python scripts/generate_extraction_report.py <csv_path> <output_dir> <report_name>

Example:
    python scripts/generate_extraction_report.py data/exports/timber_extraction_test.csv reports/extractions timber-test
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# Scenario colors (consistent with generate_maps.py)
SCENARIO_COLORS = {
    "SSP1-2.6": "#27ae60",  # Green (low emissions)
    "SSP3-7.0": "#f39c12",  # Orange (medium emissions)
    "SSP5-8.5": "#e74c3c",  # Red (high emissions)
    "ssp126": "#27ae60",
    "ssp370": "#f39c12",
    "ssp585": "#e74c3c",
    "RCP2.6": "#27ae60",
    "RCP6.0": "#f39c12",
    "RCP8.5": "#e74c3c",
    "rcp26": "#27ae60",
    "rcp60": "#f39c12",
    "rcp85": "#e74c3c",
}

# Risk score colors (1=very low to 5=very high)
RISK_COLORS = {
    1: "#27ae60",  # Very Low - Green
    2: "#2ecc71",  # Low - Light green
    3: "#f1c40f",  # Medium - Yellow
    4: "#e67e22",  # High - Orange
    5: "#e74c3c",  # Very High - Red
}

# CSS styling (consistent with generate_maps.py)
CSS_STYLES = """
<style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        background: #f5f6fa;
        color: #2c3e50;
        line-height: 1.6;
    }
    .header {
        background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        color: white;
        padding: 20px 30px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    }
    .header h1 { font-size: 1.8em; font-weight: 600; margin-bottom: 5px; }
    .header .subtitle { opacity: 0.8; font-size: 0.95em; }
    .nav {
        background: white;
        padding: 15px 30px;
        border-bottom: 1px solid #e0e0e0;
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
        align-items: center;
    }
    .nav-section {
        display: flex;
        gap: 8px;
        align-items: center;
        margin-right: 20px;
    }
    .nav-label {
        font-weight: 600;
        color: #7f8c8d;
        font-size: 0.85em;
        text-transform: uppercase;
    }
    .nav a {
        text-decoration: none;
        padding: 8px 16px;
        border-radius: 6px;
        font-size: 0.9em;
        transition: all 0.2s;
        color: #2c3e50;
        background: #ecf0f1;
    }
    .nav a:hover { background: #3498db; color: white; }
    .nav a.active { background: #3498db; color: white; font-weight: 600; }
    .content { padding: 30px; max-width: 1400px; margin: 0 auto; }
    .chart-container {
        background: white;
        border-radius: 12px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        padding: 20px;
        margin-bottom: 20px;
    }
    .chart-title {
        font-size: 1.2em;
        font-weight: 600;
        margin-bottom: 15px;
        color: #2c3e50;
    }
    .grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        gap: 20px;
    }
    .card {
        background: white;
        border-radius: 12px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        padding: 20px;
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 20px rgba(0,0,0,0.1);
    }
    .card a { text-decoration: none; color: inherit; display: block; }
    .card-title { font-weight: 600; font-size: 1.1em; margin-bottom: 8px; }
    .card-subtitle { color: #7f8c8d; font-size: 0.9em; }
    table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.9em;
    }
    th, td {
        padding: 12px 15px;
        text-align: left;
        border-bottom: 1px solid #e0e0e0;
    }
    th {
        background: #f8f9fa;
        font-weight: 600;
        color: #2c3e50;
        cursor: pointer;
    }
    th:hover { background: #ecf0f1; }
    tr:hover { background: #f8f9fa; }
    .risk-badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.85em;
        font-weight: 600;
        color: white;
    }
    .risk-1 { background: #27ae60; }
    .risk-2 { background: #2ecc71; }
    .risk-3 { background: #f1c40f; color: #2c3e50; }
    .risk-4 { background: #e67e22; }
    .risk-5 { background: #e74c3c; }
    .summary-stats {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
        gap: 15px;
        margin-bottom: 20px;
    }
    .stat-box {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        text-align: center;
    }
    .stat-value { font-size: 1.8em; font-weight: 700; color: #2c3e50; }
    .stat-label { font-size: 0.85em; color: #7f8c8d; margin-top: 5px; }
</style>
"""

# JavaScript for sortable tables
JS_SORTABLE = """
<script>
function sortTable(table, col, reverse) {
    var tb = table.tBodies[0];
    var tr = Array.prototype.slice.call(tb.rows, 0);
    reverse = -((+reverse) || -1);
    tr = tr.sort(function (a, b) {
        var aVal = a.cells[col].textContent.trim();
        var bVal = b.cells[col].textContent.trim();
        var aNum = parseFloat(aVal);
        var bNum = parseFloat(bVal);
        if (!isNaN(aNum) && !isNaN(bNum)) {
            return reverse * (aNum - bNum);
        }
        return reverse * aVal.localeCompare(bVal);
    });
    for(var i = 0; i < tr.length; ++i) tb.appendChild(tr[i]);
}

function makeSortable(table) {
    var th = table.tHead.rows[0].cells;
    for (var i = 0; i < th.length; i++) {
        (function(i) {
            var dir = 1;
            th[i].addEventListener('click', function() {
                sortTable(table, i, dir = 1 - dir);
            });
        })(i);
    }
}

document.addEventListener('DOMContentLoaded', function() {
    var tables = document.querySelectorAll('table.sortable');
    tables.forEach(function(table) { makeSortable(table); });
});
</script>
"""


class ExtractionReportGenerator:
    """Generate HTML visualization reports for extraction results."""

    def __init__(self, csv_path: Path, output_dir: Path):
        """Initialize generator.

        Args:
            csv_path: Path to extraction CSV file
            output_dir: Base directory for output reports
        """
        self.csv_path = Path(csv_path)
        self.output_dir = Path(output_dir)
        self.df: Optional[pd.DataFrame] = None
        self.locations: List[str] = []
        self.scenarios: List[str] = []
        self.decades: List[int] = []
        self.hazard_measures: List[str] = []

    def load_data(self) -> None:
        """Load and preprocess extraction data."""
        self.df = pd.read_csv(self.csv_path)

        # Extract unique values
        self.locations = sorted(self.df["Location"].unique())
        self.scenarios = sorted(self.df["Scenario"].unique())
        self.decades = sorted(self.df["Decade"].unique())
        self.hazard_measures = sorted(self.df["Hazard_Measure"].unique())

        print(f"Loaded {len(self.df)} rows")
        print(f"  Locations: {len(self.locations)}")
        print(f"  Scenarios: {len(self.scenarios)}")
        print(f"  Decades: {len(self.decades)}")
        print(f"  Hazard measures: {len(self.hazard_measures)}")

    def get_scenario_color(self, scenario: str) -> str:
        """Get color for a scenario."""
        return SCENARIO_COLORS.get(scenario, "#95a5a6")

    def _build_nav(
        self,
        current_page: str,
        report_dir: Path,
    ) -> str:
        """Build navigation HTML."""
        nav_items = []

        # Index link
        active = "active" if current_page == "index" else ""
        nav_items.append(f'<a href="index.html" class="{active}">Overview</a>')

        # Location time series
        nav_items.append('<span class="nav-label">Locations:</span>')
        for loc in self.locations:
            safe_loc = loc.replace(" ", "_").lower()
            active = "active" if current_page == f"ts_{safe_loc}" else ""
            nav_items.append(
                f'<a href="{safe_loc}_timeseries.html" class="{active}">{loc}</a>'
            )

        # Comparison pages
        nav_items.append('<span class="nav-label">Compare:</span>')
        for decade in [2050, 2090]:
            if decade in self.decades:
                active = "active" if current_page == f"comp_{decade}" else ""
                nav_items.append(
                    f'<a href="comparison_{decade}.html" class="{active}">{decade}s</a>'
                )

        # Other views
        nav_items.append('<span class="nav-label">Views:</span>')
        active = "active" if current_page == "heatmap" else ""
        nav_items.append(f'<a href="risk_heatmap.html" class="{active}">Risk Heatmap</a>')
        active = "active" if current_page == "table" else ""
        nav_items.append(f'<a href="data_table.html" class="{active}">Data Table</a>')

        return f'<nav class="nav">{" ".join(nav_items)}</nav>'

    def _build_page(
        self,
        title: str,
        subtitle: str,
        content: str,
        current_page: str,
        report_dir: Path,
    ) -> str:
        """Build complete HTML page."""
        nav = self._build_nav(current_page, report_dir)

        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    {CSS_STYLES}
    {JS_SORTABLE}
</head>
<body>
    <header class="header">
        <h1>{title}</h1>
        <div class="subtitle">{subtitle}</div>
    </header>
    {nav}
    <main class="content">
        {content}
    </main>
</body>
</html>
"""

    def generate_timeseries_page(
        self,
        location: str,
        report_dir: Path,
    ) -> None:
        """Generate time series page for a location."""
        loc_df = self.df[self.df["Location"] == location]

        # Create figure with subplots for each metric
        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                "Raw Hazard Value",
                "Percentile Score",
                "Decadal Trend",
                "Confidence Interval",
            ),
            vertical_spacing=0.12,
            horizontal_spacing=0.1,
        )

        for scenario in self.scenarios:
            scenario_df = loc_df[loc_df["Scenario"] == scenario].sort_values("Decade")
            color = self.get_scenario_color(scenario)

            # Raw value
            fig.add_trace(
                go.Scatter(
                    x=scenario_df["Decade"],
                    y=scenario_df["Raw_Hazard_Value"],
                    mode="lines+markers",
                    name=scenario,
                    line=dict(color=color, width=2),
                    marker=dict(size=8),
                    legendgroup=scenario,
                    showlegend=True,
                ),
                row=1,
                col=1,
            )

            # Percentile
            fig.add_trace(
                go.Scatter(
                    x=scenario_df["Decade"],
                    y=scenario_df["Percentile_Score"],
                    mode="lines+markers",
                    name=scenario,
                    line=dict(color=color, width=2),
                    marker=dict(size=8),
                    legendgroup=scenario,
                    showlegend=False,
                ),
                row=1,
                col=2,
            )

            # Trend
            fig.add_trace(
                go.Scatter(
                    x=scenario_df["Decade"],
                    y=scenario_df["Decadal_Trend_Strength"],
                    mode="lines+markers",
                    name=scenario,
                    line=dict(color=color, width=2),
                    marker=dict(size=8),
                    legendgroup=scenario,
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

            # Confidence interval
            if "Raw_Hazard_Value_25th" in scenario_df.columns:
                # Upper bound
                fig.add_trace(
                    go.Scatter(
                        x=scenario_df["Decade"],
                        y=scenario_df["Raw_Hazard_Value_75th"],
                        mode="lines",
                        name=f"{scenario} 75th",
                        line=dict(color=color, width=1, dash="dash"),
                        legendgroup=scenario,
                        showlegend=False,
                    ),
                    row=2,
                    col=2,
                )
                # Lower bound
                fig.add_trace(
                    go.Scatter(
                        x=scenario_df["Decade"],
                        y=scenario_df["Raw_Hazard_Value_25th"],
                        mode="lines",
                        name=f"{scenario} 25th",
                        line=dict(color=color, width=1, dash="dash"),
                        fill="tonexty",
                        fillcolor=f"rgba{tuple(list(int(color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4)) + [0.2])}",
                        legendgroup=scenario,
                        showlegend=False,
                    ),
                    row=2,
                    col=2,
                )
                # Median line
                fig.add_trace(
                    go.Scatter(
                        x=scenario_df["Decade"],
                        y=scenario_df["Raw_Hazard_Value"],
                        mode="lines",
                        name=scenario,
                        line=dict(color=color, width=2),
                        legendgroup=scenario,
                        showlegend=False,
                    ),
                    row=2,
                    col=2,
                )

        fig.update_layout(
            height=700,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
            margin=dict(l=60, r=40, t=80, b=60),
        )

        fig.update_xaxes(title_text="Decade", row=2, col=1)
        fig.update_xaxes(title_text="Decade", row=2, col=2)
        fig.update_yaxes(title_text="Value", row=1, col=1)
        fig.update_yaxes(title_text="Percentile (0-100)", row=1, col=2)
        fig.update_yaxes(title_text="Trend Strength", row=2, col=1)
        fig.update_yaxes(title_text="Value (with CI)", row=2, col=2)

        # Build content
        chart_html = fig.to_html(full_html=False, include_plotlyjs="cdn")

        # Summary stats
        latest = loc_df[loc_df["Decade"] == max(self.decades)]
        stats_html = '<div class="summary-stats">'
        for scenario in self.scenarios:
            s_data = latest[latest["Scenario"] == scenario]
            if not s_data.empty:
                val = s_data["Raw_Hazard_Value"].iloc[0]
                pct = s_data["Percentile_Score"].iloc[0]
                risk = s_data["Relative_Hazard_Score"].iloc[0]
                risk_num = int(s_data["Relative_Hazard_Score_Number"].iloc[0]) if pd.notna(s_data["Relative_Hazard_Score_Number"].iloc[0]) else 3
                stats_html += f"""
                <div class="stat-box">
                    <div class="stat-value" style="color: {self.get_scenario_color(scenario)}">{val:.3f}</div>
                    <div class="stat-label">{scenario}</div>
                    <div class="risk-badge risk-{risk_num}">{risk}</div>
                </div>
                """
        stats_html += "</div>"

        content = f"""
        <h2 class="chart-title">Climate Projections for {location}</h2>
        <p>Latest decade ({max(self.decades)}s) summary:</p>
        {stats_html}
        <div class="chart-container">
            {chart_html}
        </div>
        """

        # Get hazard info
        hazard = loc_df["Hazard"].iloc[0] if "Hazard" in loc_df.columns else "Climate"
        measure = loc_df["Hazard_Measure"].iloc[0] if "Hazard_Measure" in loc_df.columns else "Metric"

        safe_loc = location.replace(" ", "_").lower()
        html = self._build_page(
            title=f"{location} - Time Series",
            subtitle=f"{hazard}: {measure}",
            content=content,
            current_page=f"ts_{safe_loc}",
            report_dir=report_dir,
        )

        output_path = report_dir / f"{safe_loc}_timeseries.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"  Generated: {output_path.name}")

    def generate_comparison_page(
        self,
        decade: int,
        report_dir: Path,
    ) -> None:
        """Generate location comparison page for a specific decade."""
        decade_df = self.df[self.df["Decade"] == decade]

        # Create grouped bar chart
        fig = go.Figure()

        for scenario in self.scenarios:
            scenario_df = decade_df[decade_df["Scenario"] == scenario]
            if scenario_df.empty:
                continue

            # Group by location and take first value (handles multiple hazard measures)
            scenario_agg = scenario_df.groupby("Location").agg({
                "Raw_Hazard_Value": "first",
                "Percentile_Score": "first",
            }).reset_index()

            # Sort by location order (only include locations that exist in this scenario)
            loc_order = [loc for loc in self.locations if loc in scenario_agg["Location"].values]
            scenario_agg = scenario_agg.set_index("Location").loc[loc_order].reset_index()

            fig.add_trace(
                go.Bar(
                    name=scenario,
                    x=scenario_agg["Location"],
                    y=scenario_agg["Raw_Hazard_Value"],
                    marker_color=self.get_scenario_color(scenario),
                    text=scenario_agg["Raw_Hazard_Value"].round(3),
                    textposition="auto",
                )
            )

        fig.update_layout(
            barmode="group",
            height=500,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
            xaxis_title="Location",
            yaxis_title="Raw Hazard Value",
            margin=dict(l=60, r=40, t=80, b=100),
        )
        fig.update_xaxes(tickangle=45)

        # Percentile comparison
        fig2 = go.Figure()

        for scenario in self.scenarios:
            scenario_df = decade_df[decade_df["Scenario"] == scenario]
            if scenario_df.empty:
                continue

            scenario_agg = scenario_df.groupby("Location").agg({
                "Raw_Hazard_Value": "first",
                "Percentile_Score": "first",
            }).reset_index()

            loc_order = [loc for loc in self.locations if loc in scenario_agg["Location"].values]
            scenario_agg = scenario_agg.set_index("Location").loc[loc_order].reset_index()

            fig2.add_trace(
                go.Bar(
                    name=scenario,
                    x=scenario_agg["Location"],
                    y=scenario_agg["Percentile_Score"],
                    marker_color=self.get_scenario_color(scenario),
                )
            )

        fig2.update_layout(
            barmode="group",
            height=500,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
            xaxis_title="Location",
            yaxis_title="Percentile Score (0-100)",
            margin=dict(l=60, r=40, t=80, b=100),
        )
        fig2.update_xaxes(tickangle=45)

        chart1_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
        chart2_html = fig2.to_html(full_html=False, include_plotlyjs=False)

        content = f"""
        <div class="chart-container">
            <h3 class="chart-title">Raw Hazard Value - {decade}s</h3>
            {chart1_html}
        </div>
        <div class="chart-container">
            <h3 class="chart-title">Percentile Score - {decade}s</h3>
            {chart2_html}
        </div>
        """

        html = self._build_page(
            title=f"Location Comparison - {decade}s",
            subtitle=f"Comparing {len(self.locations)} locations across {len(self.scenarios)} scenarios",
            content=content,
            current_page=f"comp_{decade}",
            report_dir=report_dir,
        )

        output_path = report_dir / f"comparison_{decade}.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"  Generated: {output_path.name}")

    def generate_risk_heatmap(self, report_dir: Path) -> None:
        """Generate risk score heatmap page."""
        # Create heatmap data for each scenario
        content_parts = []

        for scenario in self.scenarios:
            scenario_df = self.df[self.df["Scenario"] == scenario]

            # Pivot to get locations x decades
            pivot = scenario_df.pivot_table(
                index="Location",
                columns="Decade",
                values="Relative_Hazard_Score_Number",
                aggfunc="first",
            )

            # Reorder locations
            pivot = pivot.reindex(self.locations)

            # Create heatmap
            fig = go.Figure(
                data=go.Heatmap(
                    z=pivot.values,
                    x=[str(d) for d in pivot.columns],
                    y=pivot.index,
                    colorscale=[
                        [0, "#27ae60"],
                        [0.25, "#2ecc71"],
                        [0.5, "#f1c40f"],
                        [0.75, "#e67e22"],
                        [1, "#e74c3c"],
                    ],
                    zmin=1,
                    zmax=5,
                    text=pivot.values,
                    texttemplate="%{text:.0f}",
                    textfont={"size": 12, "color": "white"},
                    hovertemplate="Location: %{y}<br>Decade: %{x}<br>Risk Score: %{z}<extra></extra>",
                    colorbar=dict(
                        title="Risk Score",
                        tickvals=[1, 2, 3, 4, 5],
                        ticktext=["Very Low", "Low", "Medium", "High", "Very High"],
                    ),
                )
            )

            fig.update_layout(
                height=max(400, len(self.locations) * 40),
                xaxis_title="Decade",
                yaxis_title="Location",
                margin=dict(l=150, r=100, t=50, b=60),
            )

            chart_html = fig.to_html(full_html=False, include_plotlyjs="cdn" if scenario == self.scenarios[0] else False)

            content_parts.append(f"""
            <div class="chart-container">
                <h3 class="chart-title">{scenario}</h3>
                {chart_html}
            </div>
            """)

        content = "\n".join(content_parts)

        # Add legend
        legend_html = """
        <div class="chart-container">
            <h3 class="chart-title">Risk Score Legend</h3>
            <div style="display: flex; gap: 20px; flex-wrap: wrap; margin-top: 10px;">
                <span class="risk-badge risk-1">1 = Very Low</span>
                <span class="risk-badge risk-2">2 = Low</span>
                <span class="risk-badge risk-3">3 = Medium</span>
                <span class="risk-badge risk-4">4 = High</span>
                <span class="risk-badge risk-5">5 = Very High</span>
            </div>
        </div>
        """

        html = self._build_page(
            title="Risk Score Heatmap",
            subtitle="Relative hazard scores by location and decade",
            content=legend_html + content,
            current_page="heatmap",
            report_dir=report_dir,
        )

        output_path = report_dir / "risk_heatmap.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"  Generated: {output_path.name}")

    def generate_data_table(self, report_dir: Path) -> None:
        """Generate sortable data table page."""
        # Select columns to display
        display_cols = [
            "Location",
            "Decade",
            "Scenario",
            "Raw_Hazard_Value",
            "Percentile_Score",
            "Relative_Hazard_Score",
            "Decadal_Trend_Strength",
            "Raw_Hazard_Value_25th",
            "Raw_Hazard_Value_75th",
        ]

        # Filter to existing columns
        display_cols = [c for c in display_cols if c in self.df.columns]

        # Build table HTML
        table_html = '<table class="sortable"><thead><tr>'
        for col in display_cols:
            table_html += f"<th>{col.replace('_', ' ')}</th>"
        table_html += "</tr></thead><tbody>"

        for _, row in self.df[display_cols].iterrows():
            table_html += "<tr>"
            for col in display_cols:
                val = row[col]
                if col == "Relative_Hazard_Score" and pd.notna(val):
                    risk_num = row.get("Relative_Hazard_Score_Number", 3)
                    if pd.notna(risk_num):
                        table_html += f'<td><span class="risk-badge risk-{int(risk_num)}">{val}</span></td>'
                    else:
                        table_html += f"<td>{val}</td>"
                elif isinstance(val, float) and pd.notna(val):
                    table_html += f"<td>{val:.4f}</td>"
                else:
                    table_html += f"<td>{val}</td>"
            table_html += "</tr>"

        table_html += "</tbody></table>"

        content = f"""
        <div class="chart-container">
            <h3 class="chart-title">Full Data Table</h3>
            <p style="margin-bottom: 15px; color: #7f8c8d;">Click column headers to sort. {len(self.df)} rows total.</p>
            {table_html}
        </div>
        """

        html = self._build_page(
            title="Data Table",
            subtitle="Complete extraction results",
            content=content,
            current_page="table",
            report_dir=report_dir,
        )

        output_path = report_dir / "data_table.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"  Generated: {output_path.name}")

    def generate_index(self, report_dir: Path, name: str) -> None:
        """Generate index page with overview and navigation cards."""
        # Summary statistics
        stats_html = f"""
        <div class="summary-stats">
            <div class="stat-box">
                <div class="stat-value">{len(self.locations)}</div>
                <div class="stat-label">Locations</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{len(self.scenarios)}</div>
                <div class="stat-label">Scenarios</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{len(self.decades)}</div>
                <div class="stat-label">Decades</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{len(self.df)}</div>
                <div class="stat-label">Total Rows</div>
            </div>
        </div>
        """

        # Location cards
        location_cards = '<h2 style="margin: 30px 0 15px 0;">Location Reports</h2><div class="grid">'
        for loc in self.locations:
            safe_loc = loc.replace(" ", "_").lower()
            loc_df = self.df[self.df["Location"] == loc]
            region = loc_df["Region"].iloc[0] if "Region" in loc_df.columns else "N/A"

            # Get latest risk score
            latest = loc_df[loc_df["Decade"] == max(self.decades)]
            avg_risk = latest["Relative_Hazard_Score_Number"].mean()
            risk_text = "N/A"
            risk_class = 3
            if pd.notna(avg_risk):
                risk_class = int(round(avg_risk))
                risk_labels = {1: "Very Low", 2: "Low", 3: "Medium", 4: "High", 5: "Very High"}
                risk_text = risk_labels.get(risk_class, "Medium")

            location_cards += f"""
            <div class="card">
                <a href="{safe_loc}_timeseries.html">
                    <div class="card-title">{loc}</div>
                    <div class="card-subtitle">Region: {region}</div>
                    <div style="margin-top: 10px;">
                        <span class="risk-badge risk-{risk_class}">{risk_text} Risk ({max(self.decades)}s)</span>
                    </div>
                </a>
            </div>
            """
        location_cards += "</div>"

        # View cards
        view_cards = '<h2 style="margin: 30px 0 15px 0;">Analysis Views</h2><div class="grid">'

        # Comparison cards
        for decade in [2050, 2090]:
            if decade in self.decades:
                view_cards += f"""
                <div class="card">
                    <a href="comparison_{decade}.html">
                        <div class="card-title">Location Comparison - {decade}s</div>
                        <div class="card-subtitle">Compare all locations at {decade}s</div>
                    </a>
                </div>
                """

        view_cards += """
        <div class="card">
            <a href="risk_heatmap.html">
                <div class="card-title">Risk Heatmap</div>
                <div class="card-subtitle">Visual risk matrix by location and decade</div>
            </a>
        </div>
        <div class="card">
            <a href="data_table.html">
                <div class="card-title">Data Table</div>
                <div class="card-subtitle">Full sortable data table</div>
            </a>
        </div>
        """
        view_cards += "</div>"

        content = stats_html + location_cards + view_cards

        html = self._build_page(
            title=f"Extraction Report: {name}",
            subtitle=f"Generated from {self.csv_path.name}",
            content=content,
            current_page="index",
            report_dir=report_dir,
        )

        output_path = report_dir / "index.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"  Generated: {output_path.name}")

    def generate_all(self, name: str) -> Path:
        """Generate all report pages.

        Args:
            name: Report name (used for directory and title)

        Returns:
            Path to the report directory
        """
        self.load_data()

        report_dir = self.output_dir / name
        report_dir.mkdir(parents=True, exist_ok=True)

        print(f"\nGenerating reports to: {report_dir}")

        # Generate each view
        print("\nGenerating time series pages...")
        for location in self.locations:
            self.generate_timeseries_page(location, report_dir)

        print("\nGenerating comparison pages...")
        for decade in [2050, 2090]:
            if decade in self.decades:
                self.generate_comparison_page(decade, report_dir)

        print("\nGenerating heatmap...")
        self.generate_risk_heatmap(report_dir)

        print("\nGenerating data table...")
        self.generate_data_table(report_dir)

        print("\nGenerating index...")
        self.generate_index(report_dir, name)

        print(f"\nReport complete: {report_dir / 'index.html'}")
        return report_dir


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate HTML visualization reports for extraction results"
    )
    parser.add_argument("csv_path", help="Path to extraction CSV file")
    parser.add_argument("output_dir", help="Base directory for output reports")
    parser.add_argument("name", help="Report name (used for directory and title)")

    args = parser.parse_args()

    generator = ExtractionReportGenerator(
        csv_path=Path(args.csv_path),
        output_dir=Path(args.output_dir),
    )
    generator.generate_all(args.name)


if __name__ == "__main__":
    main()
