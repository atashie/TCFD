"""Generate synoptic and diagnostic maps for processed climate data.

Creates individual HTML files for each map collection, comparing 2020s (current)
to 2090s (end-of-century) values across scenarios and metrics.

Following the process-metrics skill requirements:
- Individual HTML files per collection (avoid browser crashes)
- 2020s vs 2090s comparison structure
- 6σ anomaly detection (flag only, do not alter data)
"""

import numpy as np
import xarray as xr
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import json


# Configuration
DECADES = {
    "current": 2020,
    "future": 2090
}
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
SCENARIO_LABELS = {
    "ssp126": "SSP1-2.6 (Low Emissions)",
    "ssp370": "SSP3-7.0 (Intermediate)",
    "ssp585": "SSP5-8.5 (High Emissions)"
}
ANOMALY_SIGMA = 6  # Flag values > 6σ from 2020s mean

# Color scales by metric type
COLORSCALES = {
    "median": "Viridis",
    "percentile": "RdYlBu",
    "trend": "RdBu_r",
    "lower_ci": "Viridis",
    "upper_ci": "Viridis",
    "change": "RdBu",
    "anomaly": "Reds"
}

# Metric descriptions
METRIC_DESCRIPTIONS = {
    "median": "Ensemble Median Value",
    "percentile": "Percentile Rank (vs 2020s baseline)",
    "trend": "Decadal Trend",
    "lower_ci": "Lower Confidence Interval (25th percentile)",
    "upper_ci": "Upper Confidence Interval (75th percentile)",
    "change": "Absolute Change (2090s - 2020s)",
    "anomaly": "Anomaly Detection (>6σ from 2020s mean)"
}

# Colorbar labels by metric type (templates with {long_name} and {units} placeholders)
COLORBAR_LABELS = {
    "median": "{long_name} [{units}]",
    "percentile": "Percentile rank [1-100]",
    "trend": "Trend [{units} decade⁻¹]",
    "lower_ci": "{long_name} [{units}]",
    "upper_ci": "{long_name} [{units}]",
    "change": "Change [{units}]",
    "anomaly": "{long_name} [{units}]"
}


def log(msg: str):
    """Print with timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def create_map_figure(
    lons: np.ndarray,
    lats: np.ndarray,
    values: np.ndarray,
    title: str,
    colorscale: str = "Viridis",
    showscale: bool = True,
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    anomaly_mask: Optional[np.ndarray] = None,
    colorbar_title: str = "Value"
) -> go.Figure:
    """Create a Plotly geographic map figure.

    Args:
        lons: 1D array of longitudes
        lats: 1D array of latitudes
        values: 2D array of values (lat x lon)
        title: Map title
        colorscale: Plotly colorscale name
        showscale: Whether to show colorbar
        cmin: Minimum color value
        cmax: Maximum color value
        anomaly_mask: Optional 2D boolean mask for anomalies
        colorbar_title: Title for the colorbar (e.g., "Groundwater runoff [kg m⁻² s⁻¹]")

    Returns:
        Plotly Figure object
    """
    # Create meshgrid and flatten
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    lon_flat = lon_grid.flatten()
    lat_flat = lat_grid.flatten()
    val_flat = values.flatten()

    # Remove NaN values
    valid_mask = ~np.isnan(val_flat)
    lon_valid = lon_flat[valid_mask]
    lat_valid = lat_flat[valid_mask]
    val_valid = val_flat[valid_mask]

    # Auto-calculate color range if not specified
    if cmin is None:
        cmin = np.percentile(val_valid, 2) if len(val_valid) > 0 else 0
    if cmax is None:
        cmax = np.percentile(val_valid, 98) if len(val_valid) > 0 else 1

    fig = go.Figure()

    # Main data scatter
    fig.add_trace(go.Scattergeo(
        lon=lon_valid.tolist(),
        lat=lat_valid.tolist(),
        mode='markers',
        marker=dict(
            size=2,
            color=val_valid.tolist(),
            colorscale=colorscale,
            cmin=cmin,
            cmax=cmax,
            colorbar=dict(
                title=colorbar_title,
                exponentformat="power",  # Use ×10⁻⁶ style instead of μ
                showexponent="all"       # Show exponent on all tick labels
            ) if showscale else None,
            showscale=showscale
        ),
        hovertemplate="Lon: %{lon:.1f}<br>Lat: %{lat:.1f}<br>Value: %{marker.color:.3e}<extra></extra>"
    ))

    # Add anomaly markers if provided
    if anomaly_mask is not None:
        anom_flat = anomaly_mask.flatten()
        anom_valid = anom_flat[valid_mask]
        if np.any(anom_valid):
            anom_lons = lon_valid[anom_valid]
            anom_lats = lat_valid[anom_valid]
            fig.add_trace(go.Scattergeo(
                lon=anom_lons.tolist(),
                lat=anom_lats.tolist(),
                mode='markers',
                marker=dict(
                    size=6,
                    color='red',
                    symbol='x',
                    line=dict(width=1, color='darkred')
                ),
                name=f'Anomalies (n={len(anom_lons)})',
                hovertemplate="ANOMALY<br>Lon: %{lon:.1f}<br>Lat: %{lat:.1f}<extra></extra>"
            ))

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=14)),
        geo=dict(
            projection_type='equirectangular',
            showland=True,
            landcolor='rgb(243, 243, 243)',
            showocean=True,
            oceancolor='rgb(204, 229, 255)',
            showcoastlines=True,
            coastlinecolor='rgb(100, 100, 100)',
            coastlinewidth=0.5,
            showlakes=True,
            lakecolor='rgb(204, 229, 255)',
            showcountries=True,
            countrycolor='rgb(180, 180, 180)',
            lataxis=dict(range=[-90, 90]),
            lonaxis=dict(range=[-180, 180])
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        height=300,
        showlegend=bool(anomaly_mask is not None and np.any(anomaly_mask))
    )

    return fig


def generate_html_header(variable: str, metric: str, scenario: str) -> str:
    """Generate HTML header with navigation for per-scenario files.

    Args:
        variable: Variable code (e.g., 'qg')
        metric: Metric type (e.g., 'median', 'percentile', 'trend')
        scenario: Scenario code (e.g., 'ssp126', 'ssp370', 'ssp585')
    """
    # Build metric navigation (same metric, all scenarios)
    metric_nav = []
    for scen in SCENARIOS:
        active = "active" if scen == scenario else ""
        label = SCENARIO_LABELS.get(scen, scen).split()[0]  # e.g., "SSP1-2.6"
        metric_nav.append(f'<a href="{variable}_{metric}_{scen}.html" class="{active}">{label}</a>')

    # Build cross-metric navigation (same scenario, other metrics)
    cross_nav = []
    for m in ["median", "percentile", "trend", "confidence", "change", "anomaly"]:
        active = "active" if m == metric else ""
        cross_nav.append(f'<a href="{variable}_{m}_{scenario}.html" class="{active}">{m.title()}</a>')

    scenario_label = SCENARIO_LABELS.get(scenario, scenario)

    return f"""<!DOCTYPE html>
<html>
<head>
    <title>{variable.upper()} - {METRIC_DESCRIPTIONS.get(metric, metric)} - {scenario_label}</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; margin: -20px -20px 20px -20px; }}
        .header h1 {{ margin: 0 0 5px 0; }}
        .header .scenario {{ font-size: 18px; color: #3498db; margin-bottom: 5px; }}
        .header .subtitle {{ opacity: 0.8; font-size: 14px; }}
        .nav {{ background: #34495e; padding: 10px 20px; margin: -20px -20px 10px -20px; }}
        .nav a {{ color: white; text-decoration: none; padding: 8px 16px; margin-right: 5px; border-radius: 4px; }}
        .nav a:hover {{ background: #4a6278; }}
        .nav a.active {{ background: #3498db; }}
        .nav-section {{ margin-bottom: 5px; }}
        .nav-label {{ color: #95a5a6; font-size: 12px; margin-right: 10px; }}
        .comparison-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 30px; }}
        .map-container {{ background: white; border-radius: 8px; padding: 15px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .map-label {{ font-size: 14px; font-weight: bold; color: #2c3e50; text-align: center; margin-bottom: 10px; }}
        .footer {{ text-align: center; color: #666; padding: 20px; font-size: 12px; }}
        .stats-box {{ background: #ecf0f1; padding: 15px; border-radius: 4px; margin-bottom: 20px; }}
        .stats-box h3 {{ margin-top: 0; }}
        @media (max-width: 1200px) {{ .comparison-grid {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{variable.upper()} - {METRIC_DESCRIPTIONS.get(metric, metric)}</h1>
        <div class="scenario">{scenario_label}</div>
        <div class="subtitle">2020s (Current) vs 2090s (End of Century)</div>
    </div>
    <div class="nav">
        <div class="nav-section">
            <span class="nav-label">Scenario:</span>
            {' '.join(metric_nav)}
        </div>
        <div class="nav-section">
            <span class="nav-label">Metric:</span>
            {' '.join(cross_nav)}
            <a href="index.html">Index</a>
        </div>
    </div>
"""


def generate_html_footer(timestamp: str) -> str:
    """Generate HTML footer."""
    return f"""
    <div class="footer">
        Generated: {timestamp}<br>
        Data source: ISIMIP3b | Processing: process-metrics skill
    </div>
</body>
</html>
"""


class MapCollectionGenerator:
    """Generate individual HTML map collections for processed data."""

    def __init__(self, processed_dir: Path, output_dir: Path):
        self.processed_dir = Path(processed_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.data: Dict[str, xr.Dataset] = {}
        self.baseline_stats: Dict[str, Tuple[float, float]] = {}  # (mean, std) per scenario
        self.variable_units: str = ""
        self.variable_long_name: str = ""

    def load_data(self, variable: str):
        """Load all scenario data for a variable."""
        log(f"Loading data for {variable}...")
        self.data = {}

        for scenario in SCENARIOS:
            fpath = self.processed_dir / f"{variable}_{scenario}_processed.nc"
            if fpath.exists():
                self.data[scenario] = xr.open_dataset(fpath)
                log(f"  Loaded {scenario}: {fpath.name}")
            else:
                log(f"  WARNING: {fpath.name} not found")

        if not self.data:
            raise ValueError(f"No data files found for {variable}")

        # Extract metadata from first loaded dataset
        first_ds = list(self.data.values())[0]
        self.variable_units = first_ds.attrs.get("units", "")
        self.variable_long_name = first_ds.attrs.get("long_name", variable)
        log(f"  Metadata: {self.variable_long_name} [{self.variable_units}]")

        # Calculate baseline statistics for anomaly detection
        self._calculate_baseline_stats(variable)

    def _calculate_baseline_stats(self, variable: str):
        """Calculate 2020s mean and std for anomaly detection."""
        log("Calculating 2020s baseline statistics...")

        for scenario, ds in self.data.items():
            if DECADES["current"] in ds.decade.values:
                data_2020s = ds["median"].sel(decade=DECADES["current"]).values
                valid_data = data_2020s[~np.isnan(data_2020s)]
                if len(valid_data) > 0:
                    mean_val = float(np.mean(valid_data))
                    std_val = float(np.std(valid_data))
                    self.baseline_stats[scenario] = (mean_val, std_val)
                    log(f"  {scenario}: mean={mean_val:.3e}, std={std_val:.3e}")

    def generate_all_collections(self, variable: str):
        """Generate all map collections for a variable."""
        self.load_data(variable)

        # Create variable subdirectory
        var_dir = self.output_dir / variable
        var_dir.mkdir(exist_ok=True)

        # Generate each collection
        for metric in ["median", "percentile", "trend"]:
            self.generate_metric_comparison(variable, metric, var_dir)

        self.generate_confidence_comparison(variable, var_dir)
        self.generate_change_maps(variable, var_dir)
        self.generate_anomaly_maps(variable, var_dir)
        self.generate_index(variable, var_dir)

        # Close datasets
        for ds in self.data.values():
            ds.close()

    def generate_metric_comparison(self, variable: str, metric: str, output_dir: Path):
        """Generate 2020s vs 2090s comparison for a single metric (per-scenario files)."""
        log(f"Generating {metric} comparison maps...")

        # Get coordinate arrays from first dataset
        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        # Calculate consistent color range across all scenarios
        all_values = []
        for scenario, ds in self.data.items():
            if metric in ds.data_vars:
                for decade in [DECADES["current"], DECADES["future"]]:
                    if decade in ds.decade.values:
                        vals = ds[metric].sel(decade=decade).values
                        valid = vals[~np.isnan(vals)]
                        all_values.extend(valid.tolist())

        if all_values:
            cmin = np.percentile(all_values, 2)
            cmax = np.percentile(all_values, 98)
        else:
            cmin, cmax = 0, 1

        # Generate separate file for each scenario
        for scenario in SCENARIOS:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            if metric not in ds.data_vars:
                continue

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, metric, scenario)

            html += '<div class="comparison-grid">\n'

            for decade_label, decade in [("2020s (Current)", DECADES["current"]),
                                          ("2090s (Future)", DECADES["future"])]:
                if decade not in ds.decade.values:
                    html += f'<div class="map-container"><p>No data for {decade_label}</p></div>\n'
                    continue

                values = ds[metric].sel(decade=decade).values
                title = decade_label

                # Format colorbar label with variable metadata
                colorbar_label = COLORBAR_LABELS.get(metric, "Value").format(
                    long_name=self.variable_long_name,
                    units=self.variable_units
                )

                fig = create_map_figure(
                    lons, lats, values, title,
                    colorscale=COLORSCALES.get(metric, "Viridis"),
                    cmin=cmin, cmax=cmax,
                    colorbar_title=colorbar_label
                )

                html += '<div class="map-container">\n'
                html += f'<div class="map-label">{decade_label}</div>\n'
                html += fig.to_html(full_html=False, include_plotlyjs=False)
                html += '</div>\n'

            html += '</div>\n'
            html += generate_html_footer(timestamp)

            # Write file per scenario
            output_path = output_dir / f"{variable}_{metric}_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

    def generate_confidence_comparison(self, variable: str, output_dir: Path):
        """Generate confidence interval comparison maps (per-scenario files)."""
        log("Generating confidence interval comparison maps...")

        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        # Generate separate file for each scenario
        for scenario in SCENARIOS:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "confidence", scenario)

            for ci_metric, ci_label in [("lower_ci", "Lower CI (25th percentile)"),
                                        ("upper_ci", "Upper CI (75th percentile)")]:
                if ci_metric not in ds.data_vars:
                    continue

                html += f'<h3 style="color: #2c3e50; margin-top: 20px;">{ci_label}</h3>\n'
                html += '<div class="comparison-grid">\n'

                for decade_label, decade in [("2020s (Current)", DECADES["current"]),
                                              ("2090s (Future)", DECADES["future"])]:
                    if decade not in ds.decade.values:
                        continue

                    values = ds[ci_metric].sel(decade=decade).values
                    title = decade_label

                    # Format colorbar label for confidence interval
                    colorbar_label = COLORBAR_LABELS.get(ci_metric, "Value").format(
                        long_name=self.variable_long_name,
                        units=self.variable_units
                    )

                    fig = create_map_figure(lons, lats, values, title, colorscale="Viridis",
                                           colorbar_title=colorbar_label)

                    html += '<div class="map-container">\n'
                    html += f'<div class="map-label">{decade_label}</div>\n'
                    html += fig.to_html(full_html=False, include_plotlyjs=False)
                    html += '</div>\n'

                html += '</div>\n'

            html += generate_html_footer(timestamp)

            output_path = output_dir / f"{variable}_confidence_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

    def generate_change_maps(self, variable: str, output_dir: Path):
        """Generate absolute change maps (2090s - 2020s) per scenario."""
        log("Generating change maps...")

        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        # Calculate all changes for consistent color scaling across scenarios
        all_changes = []
        for scenario, ds in self.data.items():
            if DECADES["current"] in ds.decade.values and DECADES["future"] in ds.decade.values:
                val_2020 = ds["median"].sel(decade=DECADES["current"]).values
                val_2090 = ds["median"].sel(decade=DECADES["future"]).values
                change = val_2090 - val_2020
                valid = change[~np.isnan(change)]
                all_changes.extend(valid.tolist())

        if all_changes:
            max_abs = np.percentile(np.abs(all_changes), 98)
            cmin, cmax = -max_abs, max_abs
        else:
            cmin, cmax = -1, 1

        # Generate separate file for each scenario
        for scenario in SCENARIOS:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            if DECADES["current"] not in ds.decade.values or DECADES["future"] not in ds.decade.values:
                continue

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "change", scenario)

            html += '<div class="stats-box"><h3>Absolute Change: 2090s minus 2020s</h3>'
            html += '<p>Positive values (red) indicate increase, negative values (blue) indicate decrease.</p></div>\n'

            val_2020 = ds["median"].sel(decade=DECADES["current"]).values
            val_2090 = ds["median"].sel(decade=DECADES["future"]).values
            change = val_2090 - val_2020

            title = "Change (2090s - 2020s)"

            # Format colorbar label for change map
            colorbar_label = COLORBAR_LABELS.get("change", "Change").format(
                units=self.variable_units
            )

            fig = create_map_figure(
                lons, lats, change, title,
                colorscale="RdBu",
                cmin=cmin, cmax=cmax,
                colorbar_title=colorbar_label
            )

            # Single centered map for change
            html += '<div style="max-width: 800px; margin: 0 auto;">\n'
            html += '<div class="map-container">\n'
            html += f'<div class="map-label">{title}</div>\n'
            html += fig.to_html(full_html=False, include_plotlyjs=False)
            html += '</div>\n'
            html += '</div>\n'

            html += generate_html_footer(timestamp)

            output_path = output_dir / f"{variable}_change_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

    def generate_anomaly_maps(self, variable: str, output_dir: Path):
        """Generate anomaly detection maps (6-sigma threshold) per scenario."""
        log(f"Generating anomaly maps (threshold: {ANOMALY_SIGMA}-sigma)...")

        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        anomaly_summary = {}

        # Generate separate file for each scenario
        for scenario in SCENARIOS:
            if scenario not in self.data or scenario not in self.baseline_stats:
                continue

            ds = self.data[scenario]
            mean_2020s, std_2020s = self.baseline_stats[scenario]

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "anomaly", scenario)

            html += f'<div class="stats-box"><h3>Anomaly Detection: >{ANOMALY_SIGMA}σ from 2020s Mean</h3>'
            html += '<p>Red X markers indicate grid cells where the value deviates more than '
            html += f'{ANOMALY_SIGMA} standard deviations from the 2020s global mean.</p>'
            html += '<p><strong>Note:</strong> Anomalies are flagged for inspection only - data has NOT been altered.</p></div>\n'

            html += '<div class="comparison-grid">\n'

            scenario_anomalies = {}

            for decade_label, decade in [("2020s (Current)", DECADES["current"]),
                                          ("2090s (Future)", DECADES["future"])]:
                if decade not in ds.decade.values:
                    continue

                values = ds["median"].sel(decade=decade).values

                # Calculate anomaly mask
                anomaly_mask = np.abs(values - mean_2020s) > (ANOMALY_SIGMA * std_2020s)
                n_anomalies = int(np.sum(anomaly_mask & ~np.isnan(values)))

                decade_key = decade_label.split()[0]  # "2020s" or "2090s"
                scenario_anomalies[decade_key] = n_anomalies

                title = f"{decade_label} - {n_anomalies} anomalies"

                # Format colorbar label for anomaly map
                colorbar_label = COLORBAR_LABELS.get("anomaly", "Value").format(
                    long_name=self.variable_long_name,
                    units=self.variable_units
                )

                fig = create_map_figure(
                    lons, lats, values, title,
                    colorscale="Viridis",
                    anomaly_mask=anomaly_mask,
                    colorbar_title=colorbar_label
                )

                html += '<div class="map-container">\n'
                html += f'<div class="map-label">{decade_label}</div>\n'
                html += fig.to_html(full_html=False, include_plotlyjs=False)
                html += '</div>\n'

            html += '</div>\n'

            # Summary statistics
            html += '<div class="stats-box">\n'
            html += f'<strong>Baseline (2020s):</strong> mean = {mean_2020s:.3e}, std = {std_2020s:.3e}<br>\n'
            html += f'<strong>Threshold:</strong> |value - mean| > {ANOMALY_SIGMA} × {std_2020s:.3e} = {ANOMALY_SIGMA * std_2020s:.3e}<br>\n'
            for decade_key, count in scenario_anomalies.items():
                html += f'<strong>{decade_key}:</strong> {count} cells flagged<br>\n'
            html += '</div>\n'

            html += generate_html_footer(timestamp)

            output_path = output_dir / f"{variable}_anomaly_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

            anomaly_summary[scenario] = scenario_anomalies

        # Save anomaly summary as JSON
        summary_path = output_dir / f"{variable}_anomaly_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(anomaly_summary, f, indent=2)
        log(f"  Saved: {summary_path.name}")

    def generate_index(self, variable: str, output_dir: Path):
        """Generate index page with grid layout (metrics × scenarios)."""
        log("Generating index page...")

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Define metrics for the grid
        metrics = [
            ("median", "Median Values", "Ensemble median values (2020s vs 2090s)"),
            ("percentile", "Percentile Ranks", "Percentile ranks vs 2020s baseline"),
            ("trend", "Trends", "Decadal trend analysis"),
            ("confidence", "Confidence Intervals", "Lower (25th) and upper (75th) bounds"),
            ("change", "Change Maps", "Absolute change (2090s - 2020s)"),
            ("anomaly", "Anomaly Detection", f"Values >{ANOMALY_SIGMA}σ from 2020s mean"),
        ]

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>{variable.upper()} - Climate Projection Maps</title>
    <meta charset="utf-8">
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; margin: -20px -20px 20px -20px; }}
        .header h1 {{ margin: 0 0 10px 0; }}
        .header .subtitle {{ opacity: 0.8; }}
        table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px;
                 box-shadow: 0 2px 4px rgba(0,0,0,0.1); overflow: hidden; }}
        th {{ background: #34495e; color: white; padding: 15px; text-align: center; }}
        th.metric {{ background: #2c3e50; text-align: left; width: 200px; }}
        td {{ padding: 12px; text-align: center; border-bottom: 1px solid #ecf0f1; }}
        td.metric-name {{ text-align: left; background: #f8f9fa; font-weight: bold; color: #2c3e50; }}
        td.metric-desc {{ text-align: left; background: #f8f9fa; color: #666; font-size: 12px; }}
        a.btn {{ display: inline-block; background: #3498db; color: white; padding: 8px 16px;
                 border-radius: 4px; text-decoration: none; font-size: 13px; }}
        a.btn:hover {{ background: #2980b9; }}
        a.btn.ssp126 {{ background: #27ae60; }}
        a.btn.ssp126:hover {{ background: #219a52; }}
        a.btn.ssp370 {{ background: #f39c12; }}
        a.btn.ssp370:hover {{ background: #d68910; }}
        a.btn.ssp585 {{ background: #e74c3c; }}
        a.btn.ssp585:hover {{ background: #c0392b; }}
        .footer {{ text-align: center; color: #666; padding: 20px; font-size: 12px; }}
        .legend {{ margin: 20px 0; padding: 15px; background: white; border-radius: 8px;
                   box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .legend h3 {{ margin-top: 0; color: #2c3e50; }}
        .legend-item {{ display: inline-block; margin-right: 20px; }}
        .legend-color {{ display: inline-block; width: 20px; height: 20px; border-radius: 4px;
                         vertical-align: middle; margin-right: 5px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{variable.upper()} - Climate Projection Maps</h1>
        <div class="subtitle">2020s (Current) vs 2090s (End of Century) | {self.variable_long_name}</div>
    </div>

    <div class="legend">
        <h3>Scenario Legend</h3>
        <div class="legend-item">
            <span class="legend-color" style="background: #27ae60;"></span>
            <strong>SSP1-2.6</strong> - Low emissions (sustainable development)
        </div>
        <div class="legend-item">
            <span class="legend-color" style="background: #f39c12;"></span>
            <strong>SSP3-7.0</strong> - Intermediate emissions (regional rivalry)
        </div>
        <div class="legend-item">
            <span class="legend-color" style="background: #e74c3c;"></span>
            <strong>SSP5-8.5</strong> - High emissions (fossil-fueled development)
        </div>
    </div>

    <table>
        <tr>
            <th class="metric">Metric</th>
            <th>SSP1-2.6<br><small>(Low)</small></th>
            <th>SSP3-7.0<br><small>(Intermediate)</small></th>
            <th>SSP5-8.5<br><small>(High)</small></th>
        </tr>
"""

        for metric_key, metric_name, metric_desc in metrics:
            html += f"""        <tr>
            <td class="metric-name">{metric_name}<br><span style="font-weight:normal;font-size:11px;color:#888;">{metric_desc}</span></td>
            <td><a href="{variable}_{metric_key}_ssp126.html" class="btn ssp126">View</a></td>
            <td><a href="{variable}_{metric_key}_ssp370.html" class="btn ssp370">View</a></td>
            <td><a href="{variable}_{metric_key}_ssp585.html" class="btn ssp585">View</a></td>
        </tr>
"""

        html += f"""    </table>

    <div class="footer">
        Generated: {timestamp}<br>
        Data source: ISIMIP3b | Units: {self.variable_units}
    </div>
</body>
</html>
"""

        output_path = output_dir / "index.html"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        log(f"  Saved: {output_path.name}")


def main():
    """Main entry point."""
    # Get project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    processed_dir = project_root / "data" / "processed"
    output_dir = project_root / "reports" / "maps"

    log("=" * 60)
    log("Generating Synoptic and Diagnostic Maps")
    log("=" * 60)

    generator = MapCollectionGenerator(processed_dir, output_dir)

    # Generate maps for qg variable
    generator.generate_all_collections("qg")

    log("\n" + "=" * 60)
    log("Map generation complete!")
    log(f"Output directory: {output_dir / 'qg'}")
    log("=" * 60)


if __name__ == "__main__":
    main()
