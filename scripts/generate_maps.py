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

# SSP scenarios (ISIMIP3b)
SSP_SCENARIO_LABELS = {
    "ssp126": "SSP1-2.6 (Low Emissions)",
    "ssp370": "SSP3-7.0 (Intermediate)",
    "ssp585": "SSP5-8.5 (High Emissions)"
}
SSP_SCENARIO_COLORS = {
    "ssp126": "#27ae60",
    "ssp370": "#f39c12",
    "ssp585": "#e74c3c",
}

# RCP scenarios (ISIMIP2b)
RCP_SCENARIO_LABELS = {
    "rcp26": "RCP2.6 (Low Emissions)",
    "rcp45": "RCP4.5 (Intermediate-Low)",
    "rcp60": "RCP6.0 (Intermediate)",
    "rcp85": "RCP8.5 (High Emissions)",
    "picontrol": "Pre-Industrial Control",
    "historical": "Historical",
}
RCP_SCENARIO_COLORS = {
    "rcp26": "#27ae60",
    "rcp45": "#2ecc71",  # Light green - between rcp26 and rcp60
    "rcp60": "#f39c12",
    "rcp85": "#e74c3c",
    "picontrol": "#3498db",
    "historical": "#95a5a6",
}

# Legacy globals (for backward compatibility, overridden by auto-detection)
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
SCENARIO_LABELS = SSP_SCENARIO_LABELS

ANOMALY_SIGMA = 6  # Flag values > 6σ from 2020s mean

# Value class mapping for pipeline output format (from isimip-pipeline/processing/output.py)
# Maps value_class indices to generate_maps.py variable names
VALUE_CLASS_MAP = {
    0: "median",        # smoothed_median in pipeline
    1: "percentile",    # percentile rank
    2: "trend",         # LEGACY single trend -- pre-OUTPUT-SPEC layers only
    3: "significance",  # LEGACY, never populated -- pre-OUTPUT-SPEC layers only
    4: "lower_ci",      # lower_bound (Q25)
    5: "upper_ci",      # upper_bound (Q75)
}

# Metrics under the current contract (OUTPUT-SPEC.md). `trend` is retired in favour of
# the two slopes, which fail in OPPOSITE regimes and must both be mapped: sen collapses
# to 0 on zero-inflated fields, ols absorbs member level offsets when coverage is uneven.
# `trend` is still accepted so pre-spec layers keep rendering.
SLOPE_METRICS = ["ols_slope", "sen_slope"]
CONTRACT_METRICS = ["median", "percentile"] + SLOPE_METRICS + ["lower_ci", "upper_ci"]

#: The dashboard tab strip. Deliberately NOT one tab per contract variable:
#:  * `ols_slope` / `sen_slope` / `change` are all decadal-change views and belong on one
#:    "Trend" tab so the two estimators can be read against each other (they fail in
#:    OPPOSITE regimes -- see OUTPUT-SPEC.md -- so comparing them IS the diagnostic).
#:  * `lower_ci` / `upper_ci` are exactly what the "Confidence" tab already plots side by
#:    side, so separate tabs were pure duplication.
#:  * `members` shows every LSM x GCM member on a shared scale (needs {var}_members.nc).
DASHBOARD_TABS = ["median", "percentile", "trend", "confidence", "anomaly", "members"]

#: Tabs that are one page for the whole layer rather than one page per scenario.
SCENARIO_INDEPENDENT_TABS = {"members"}

#: Decades compared on the Trend tab. NOT the 2020s baseline -- its slopes are NaN by
#: contract and its change-from-itself is identically 0, so a 2020s panel is always blank.
TREND_DECADES = [2030, 2090]

#: How the symmetric (zero-centred) limit for diverging panels is chosen.
#:   95   -> the 95th percentile of |value| (current setting)
#:   None -> true max|value|
#:
#: Measured on the wildfire layer: max|ols_slope| = 38.2 %/dec against a median of 0.049,
#: so a true-max scale left 99.67% of cells inside the middle 10% of the colour range and
#: the maps read as blank. The 95th percentile spends the colour range on the body of the
#: distribution instead.
#:
#: Values BEYOND the limit are not lost or blanked: Plotly clamps `marker.color` to the
#: cmin/cmax endpoints, so anything more extreme renders in the full max (red) or min
#: (blue) colour. Both ends stay symmetric about zero either way.
SYMMETRIC_LIMIT_PERCENTILE = 95

#: Payload control. Every map is an SVG Scattergeo: one DOM marker per land cell, and the
#: JSON carries lon/lat/value per point. At 0.5 deg that is ~70,849 points per panel, so a
#: 6-panel Trend page was 11.2 MB and the 22-panel Members page 37.6 MB.
#:  * COORD_DECIMALS / VALUE_SIGFIGS trim float text without any visible change -- a 0.5 deg
#:    grid is exact at 2 dp, and 4 significant figures is far beyond what a colour can show.
#:  * MEMBERS_GRID_STRIDE block-averages the Members tab only. That tab answers "do the
#:    members differ in level or distribution?", which does not need 0.5 deg; stride 2 (1 deg)
#:    cuts its points 4x. Set to 1 for full resolution.
COORD_DECIMALS = 2
VALUE_SIGFIGS = 4
MEMBERS_GRID_STRIDE = 2

# Non-projection scenarios to exclude from report generation
# These are used to enhance baseline robustness but not shown as separate projections
EXCLUDED_SCENARIOS = {"picontrol", "historical"}


def detect_scenario_type(scenarios: List[str]) -> Tuple[Dict[str, str], Dict[str, str], str]:
    """Detect whether scenarios are SSP or RCP and return appropriate config."""
    if any(s.startswith('ssp') for s in scenarios):
        labels = {s: SSP_SCENARIO_LABELS.get(s, s) for s in scenarios}
        colors = {s: SSP_SCENARIO_COLORS.get(s, "#3498db") for s in scenarios}
        return labels, colors, "ISIMIP3b"
    elif any(s.startswith('rcp') or s in ['picontrol', 'historical'] for s in scenarios):
        labels = {s: RCP_SCENARIO_LABELS.get(s, s) for s in scenarios}
        colors = {s: RCP_SCENARIO_COLORS.get(s, "#3498db") for s in scenarios}
        return labels, colors, "ISIMIP2b"
    else:
        labels = {s: s for s in scenarios}
        colors = {s: "#3498db" for s in scenarios}
        return labels, colors, "ISIMIP"

# Color scales by metric type
COLORSCALES = {
    "median": "Viridis",
    "percentile": "RdYlBu_r",  # Reversed: low=blue (good), high=red (bad)
    "trend": "RdBu",  # legacy pre-spec layers
    "ols_slope": "RdBu",  # Blue=positive (good for "more is better" variables)
    "sen_slope": "RdBu",
    "lower_ci": "Viridis",
    "upper_ci": "Viridis",
    "change": "RdBu",
    "anomaly": "Reds"
}

# Metric descriptions
METRIC_DESCRIPTIONS = {
    "median": "Ensemble Median Value",
    "percentile": "Percentile Rank (vs 2020s baseline)",
    "ols_slope": "Decadal Trend \u2014 OLS slope",
    "sen_slope": "Decadal Trend \u2014 Theil-Sen slope",
    "lower_ci": "Lower Confidence Interval (25th percentile)",
    "upper_ci": "Upper Confidence Interval (75th percentile)",
    "change": "Absolute Change (2090s - 2020s)",
    "anomaly": "Anomaly Detection (>6σ from 2020s mean)",
    "trend": "Decadal Trend — OLS & Theil-Sen slopes and change vs baseline",
    "members": "Per-member raw values (every LSM × GCM on a shared scale)"
}

# Tab strip labels (m.title() gives "Ols_Slope"-style names, which read badly)
TAB_TITLES = {
    "median": "Median",
    "percentile": "Percentile",
    "trend": "Trend",
    "confidence": "Confidence",
    "anomaly": "Anomaly",
    "members": "Members",
}

# Hover number formats. Percentiles are integers on [1,100]; scientific notation there is
# noise, not precision.
HOVER_FORMATS = {
    "percentile": ".0f",
    "n_members": ".0f",
    "n_models": ".0f",
}
DEFAULT_HOVER_FORMAT = ".3e"

# Colorbar labels by metric type (templates with {long_name} and {units} placeholders)
COLORBAR_LABELS = {
    "median": "{long_name} [{units}]",
    "percentile": "Percentile rank [1-100]",
    "trend": "Trend [{units} decade⁻¹]",
    "ols_slope": "OLS slope [{units} decade⁻¹]",
    "sen_slope": "Theil-Sen slope [{units} decade⁻¹]",
    "lower_ci": "{long_name} [{units}]",
    "upper_ci": "{long_name} [{units}]",
    "change": "Change [{units}]",
    "anomaly": "{long_name} [{units}]"
}


def _sigfig(a: np.ndarray, n: int = VALUE_SIGFIGS) -> np.ndarray:
    """Round to n significant figures (np.round takes only scalar decimals)."""
    out = np.array(a, dtype=np.float64, copy=True)
    nz = np.isfinite(out) & (out != 0)
    if np.any(nz):
        mag = np.floor(np.log10(np.abs(out[nz])))
        factor = np.power(10.0, (n - 1 - mag))
        out[nz] = np.round(out[nz] * factor) / factor
    return out


def block_mean(values: np.ndarray, stride: int) -> np.ndarray:
    """NaN-aware block mean over (stride x stride) cells; trims a ragged edge.

    Plain slicing (values[::2, ::2]) would SAMPLE every other cell and drop the rest,
    which on a sparse hazard silently deletes burning cells. Averaging keeps them.
    """
    if stride <= 1:
        return values
    ny, nx = values.shape
    ny -= ny % stride
    nx -= nx % stride
    v = values[:ny, :nx].reshape(ny // stride, stride, nx // stride, stride)
    with np.errstate(invalid="ignore"):
        import warnings as _w
        with _w.catch_warnings():
            _w.filterwarnings("ignore", message="Mean of empty slice")
            return np.nanmean(v, axis=(1, 3))


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
    colorbar_title: str = "Value",
    reversescale: bool = False,
    hover_format: str = ".3e",
    subtitle: Optional[str] = None
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
        lon=np.round(lon_valid, COORD_DECIMALS).tolist(),
        lat=np.round(lat_valid, COORD_DECIMALS).tolist(),
        mode='markers',
        marker=dict(
            size=2,
            color=_sigfig(val_valid).tolist(),
            colorscale=colorscale,
            reversescale=reversescale,
            cmin=cmin,
            cmax=cmax,
            # NO colorbar title: a long title (e.g. "Annual burnt area [%]") reserves more
            # horizontal space than the map itself. The text moves to the figure title
            # directly above the map instead.
            colorbar=dict(
                exponentformat="power",  # Use ×10⁻⁶ style instead of μ
                showexponent="all"       # Show exponent on all tick labels
            ) if showscale else None,
            showscale=showscale
        ),
        hovertemplate=("Lon: %{lon:.1f}<br>Lat: %{lat:.1f}<br>"
                       "Value: %{marker.color:" + hover_format + "}<extra></extra>")
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

    # The figure title carries what used to be the colorbar title, sitting directly over
    # the map; `title` (the panel/decade label) becomes the smaller second line.
    heading = colorbar_title if colorbar_title else title
    if subtitle:
        heading = f"{heading}<br><span style='font-size:11px;color:#666'>{subtitle}</span>"
    elif colorbar_title and title:
        heading = f"{heading}<br><span style='font-size:11px;color:#666'>{title}</span>"

    fig.update_layout(
        title=dict(text=heading, x=0.5, xanchor="center", font=dict(size=13)),
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
        margin=dict(l=0, r=0, t=55, b=0),
        height=320,
        showlegend=bool(anomaly_mask is not None and np.any(anomaly_mask))
    )

    return fig


def generate_html_header(
    variable: str,
    metric: str,
    scenario: str,
    scenarios: List[str],
    scenario_labels: Dict[str, str]
) -> str:
    """Generate HTML header with navigation for per-scenario files.

    Args:
        variable: Variable code (e.g., 'qg')
        metric: Metric type (e.g., 'median', 'percentile', 'trend')
        scenario: Scenario code (e.g., 'ssp126', 'rcp26')
        scenarios: List of all scenarios to include in navigation
        scenario_labels: Dict mapping scenario codes to display labels
    """
    # Build metric navigation (same metric, all scenarios). Scenario-independent tabs have
    # exactly one page for the whole layer, so per-scenario links would 404 -- say so
    # instead of emitting dead links.
    metric_nav = []
    if metric in SCENARIO_INDEPENDENT_TABS:
        metric_nav.append('<span style="color:#95a5a6">all scenarios (one page)</span>')
    else:
        for scen in scenarios:
            active = "active" if scen == scenario else ""
            label = scenario_labels.get(scen, scen).split()[0]  # "SSP1-2.6" / "RCP2.6"
            metric_nav.append(
                f'<a href="{variable}_{metric}_{scen}.html" class="{active}">{label}</a>')

    # Build cross-metric navigation (same scenario, other tabs)
    cross_nav = []
    for m in DASHBOARD_TABS:
        active = "active" if m == metric else ""
        href = (f"{variable}_{m}.html" if m in SCENARIO_INDEPENDENT_TABS
                else f"{variable}_{m}_{scenario}.html")
        cross_nav.append(f'<a href="{href}" class="{active}">{TAB_TITLES.get(m, m.title())}</a>')

    scenario_label = scenario_labels.get(scenario, scenario)

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


def generate_html_footer(timestamp: str, data_source: str = "ISIMIP3b") -> str:
    """Generate HTML footer."""
    return f"""
    <div class="footer">
        Generated: {timestamp}<br>
        Data source: {data_source} | Processing: process-metrics skill
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
        # Instance variables for scenario configuration (auto-detected from data)
        self.scenarios: List[str] = []
        self.scenario_labels: Dict[str, str] = {}
        self.scenario_colors: Dict[str, str] = {}
        self.data_source: str = "ISIMIP"

    def load_data(self, variable: str):
        """Load all scenario data for a variable.

        Supports two formats:
        1. Single-file format: {variable}_processed.nc with 'scenario' dimension
        2. Per-scenario files: {variable}_{scenario}_processed.nc
        """
        log(f"Loading data for {variable}...")
        self.data = {}

        # Try single-file format first
        single_file = self.processed_dir / f"{variable}_processed.nc"
        if single_file.exists():
            log(f"  Found single-file format: {single_file.name}")
            ds = xr.open_dataset(single_file)

            # Check if scenario dimension exists
            if "scenario" in ds.dims:
                scenarios = [str(s) for s in ds.scenario.values]
                log(f"  Scenarios in file: {scenarios}")

                # Split into per-scenario datasets for compatibility
                for scenario in scenarios:
                    self.data[scenario] = ds.sel(scenario=scenario)

                # Auto-detect scenario type
                self.scenario_labels, self.scenario_colors, self.data_source = detect_scenario_type(scenarios)
                self.scenarios = scenarios
                log(f"  Detected data source: {self.data_source}")
            else:
                # Single file without scenario dimension - treat as single scenario
                scenario = ds.attrs.get("scenario", "unknown")
                self.data[scenario] = ds
                self.scenarios = [scenario]
                self.scenario_labels = {scenario: scenario}
                self.scenario_colors = {scenario: "#3498db"}
                self.data_source = "ISIMIP"
        else:
            # Fall back to per-scenario file format
            log(f"  Single-file not found, trying per-scenario files...")

            # Dynamically discover ALL scenario files from filesystem
            # Pattern: {variable}_{scenario}_processed.nc
            pattern = f"{variable}_*_processed.nc"
            scenario_files = sorted(self.processed_dir.glob(pattern))

            for fpath in scenario_files:
                # Extract scenario name from filename
                # e.g., "b30cm_rcp45_processed.nc" -> "rcp45"
                parts = fpath.stem.replace("_processed", "").split("_")
                if len(parts) >= 2:
                    scenario = parts[-1]  # Last part before "_processed"
                    self.data[scenario] = xr.open_dataset(fpath)
                    log(f"  Loaded {scenario}: {fpath.name}")

            if self.data:
                scenarios = list(self.data.keys())
                self.scenario_labels, self.scenario_colors, self.data_source = detect_scenario_type(scenarios)
                self.scenarios = scenarios
                log(f"  Detected data source: {self.data_source}")

        if not self.data:
            raise ValueError(f"No data files found for {variable}")

        # Extract metadata from first loaded dataset BEFORE format conversion
        first_ds = list(self.data.values())[0]

        # Check if this is pipeline format (has value_class dimension)
        if "value_class" in first_ds.dims:
            log("  Detected pipeline format (value_class dimension)")
            # Extract metadata from the original variable before conversion
            if variable in first_ds.data_vars:
                var_attrs = first_ds[variable].attrs
                self.variable_units = var_attrs.get("units", first_ds.attrs.get("units", ""))
                self.variable_long_name = var_attrs.get("long_name", first_ds.attrs.get("long_name", variable))
            else:
                self.variable_units = first_ds.attrs.get("units", "")
                self.variable_long_name = first_ds.attrs.get("long_name", variable)
            # Convert pipeline format to QG-style format with separate variables
            self._convert_pipeline_format(variable)
        else:
            # QG format: metadata from dataset attributes
            self.variable_units = first_ds.attrs.get("units", "")
            self.variable_long_name = first_ds.attrs.get("long_name", variable)

        # Diverging maps (trend/change): for "higher is worse" hazards (e.g. drought,
        # mortality) reverse RdBu so positive = red = worsening. Absent the attribute
        # (older processed files), keep the default RdBu (blue = positive).
        self.higher_is_worse = (
            first_ds.attrs.get("percentile_direction", "") == "higher_is_worse"
        )

        log(f"  Metadata: {self.variable_long_name} [{self.variable_units}]")
        if self.higher_is_worse:
            log("  Direction: higher_is_worse -> trend/change use reversed RdBu (red = worse)")

        # Filter out non-projection scenarios (picontrol, historical)
        # These are used to enhance baseline robustness but not shown as separate projections
        excluded = [s for s in self.scenarios if s in EXCLUDED_SCENARIOS]
        if excluded:
            log(f"  Excluding non-projection scenarios: {excluded}")
            for scenario in excluded:
                if scenario in self.data:
                    del self.data[scenario]
            self.scenarios = [s for s in self.scenarios if s not in EXCLUDED_SCENARIOS]
            # Update labels and colors to only include remaining scenarios
            self.scenario_labels = {s: self.scenario_labels[s] for s in self.scenarios}
            self.scenario_colors = {s: self.scenario_colors[s] for s in self.scenarios}
            log(f"  Remaining projection scenarios: {self.scenarios}")

        if not self.scenarios:
            raise ValueError("No projection scenarios found after filtering")

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

    def _convert_pipeline_format(self, variable: str):
        """Convert pipeline format (value_class dimension) to QG format (separate variables).

        The pipeline stores metrics in a single variable with a value_class dimension:
          data[variable][lon, lat, decade, scenario, value_class]

        This method converts each scenario's data to have separate variables:
          data[scenario]['median'][decade, lat, lon]
          data[scenario]['percentile'][decade, lat, lon]
          etc.
        """
        log("  Converting pipeline format to separate-variable format...")

        for scenario in self.scenarios:
            if scenario not in self.data:
                continue

            scenario_ds = self.data[scenario]

            # Check if the variable exists with value_class dimension
            if variable not in scenario_ds.data_vars:
                log(f"    WARNING: Variable '{variable}' not found in {scenario} dataset")
                continue

            var_data = scenario_ds[variable]

            # Create new dataset with separate variables for each metric
            new_ds = xr.Dataset()

            # Extract each metric from value_class dimension
            for vc_idx, var_name in VALUE_CLASS_MAP.items():
                if vc_idx >= len(scenario_ds.value_class):
                    continue

                # Extract the metric data
                metric_data = var_data.isel(value_class=vc_idx)

                # Handle dimension order: pipeline outputs (lon, lat, decade)
                # but generate_maps.py expects (decade, lat, lon)
                if "lon" in metric_data.dims and "lat" in metric_data.dims and "decade" in metric_data.dims:
                    # Transpose to (decade, lat, lon) order
                    metric_data = metric_data.transpose("decade", "lat", "lon")

                new_ds[var_name] = metric_data

            # Copy coordinate attributes
            # Convert decade indices (10, 20, ..., 90) to years (2010, 2020, ..., 2090) if needed
            decade_values = scenario_ds.decade.values
            if max(decade_values) < 100:  # Decade indices, not years
                decade_values = decade_values + 2000  # Convert to actual years
                log(f"    Converted decade indices to years: {list(decade_values)}")

            new_ds = new_ds.assign_coords({
                "decade": decade_values,
                "lat": scenario_ds.lat,
                "lon": scenario_ds.lon,
            })

            # Copy global attributes
            new_ds.attrs = scenario_ds.attrs.copy()

            # Replace the scenario data with converted format
            self.data[scenario] = new_ds
            log(f"    {scenario}: converted {len(VALUE_CLASS_MAP)} metrics")

        log("  Format conversion complete")

    def generate_all_collections(self, variable: str):
        """Generate all map collections for a variable."""
        self.load_data(variable)

        # Create variable subdirectory
        var_dir = self.output_dir / variable
        var_dir.mkdir(exist_ok=True)

        # Generate each collection. Which slope fields exist depends on whether the layer
        # predates OUTPUT-SPEC.md (single `trend`) or follows it (`ols_slope`+`sen_slope`),
        # so render whatever is actually present rather than assuming either.
        available = set()
        for ds in self.data.values():
            available |= set(ds.data_vars)
        wanted = ["median", "percentile"]
        present = [m for m in wanted if m in available]
        if not present:
            log(f"  WARNING: none of {wanted} found in {variable}; nothing to map")
        for metric in present:
            self.generate_metric_comparison(variable, metric, var_dir)

        # Slopes and change now live together on one Trend tab; `change` and the
        # per-slope tabs are no longer emitted as standalone pages.
        self.generate_trend_page(variable, var_dir)
        self.generate_confidence_comparison(variable, var_dir)
        self.generate_anomaly_maps(variable, var_dir)
        self.generate_members_page(variable, var_dir)
        self.generate_index(variable, var_dir)

        # Close datasets
        for ds in self.data.values():
            ds.close()

    def generate_metric_comparison(self, variable: str, metric: str, output_dir: Path):
        """Generate 2020s vs 2090s comparison for a single metric (per-scenario files).

        For metrics without a decade dimension (like trend), generates a single map.
        """
        log(f"Generating {metric} comparison maps...")

        # Get coordinate arrays from first dataset
        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        # Check if metric has decade dimension
        has_decade_dim = False
        for scenario, ds in self.data.items():
            if metric in ds.data_vars:
                has_decade_dim = 'decade' in ds[metric].dims
                break

        # Calculate consistent color range across all scenarios
        all_values = []
        for scenario, ds in self.data.items():
            if metric in ds.data_vars:
                if has_decade_dim:
                    for decade in [DECADES["current"], DECADES["future"]]:
                        if decade in ds.decade.values:
                            vals = ds[metric].sel(decade=decade).values
                            valid = vals[~np.isnan(vals)]
                            all_values.extend(valid.tolist())
                else:
                    # Metric without decade dimension (e.g., trend)
                    vals = ds[metric].values
                    valid = vals[~np.isnan(vals)]
                    all_values.extend(valid.tolist())

        if all_values:
            if metric == "trend":
                # Trend maps use symmetric scaling centered on zero (white=no change)
                max_abs = np.percentile(np.abs(all_values), 98)
                cmin, cmax = -max_abs, max_abs
            else:
                cmin = np.percentile(all_values, 2)
                cmax = np.percentile(all_values, 98)
        else:
            cmin, cmax = (-1, 1) if metric == "trend" else (0, 1)

        # Generate separate file for each scenario
        for scenario in self.scenarios:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            if metric not in ds.data_vars:
                continue

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, metric, scenario, self.scenarios, self.scenario_labels)

            html += '<div class="comparison-grid">\n'

            # Check if metric has decade dimension
            metric_has_decade = 'decade' in ds[metric].dims

            if metric_has_decade:
                # Standard case: show 2020s vs 2090s comparison
                for decade_label, decade in [("2020s (Current)", DECADES["current"]),
                                              ("2090s (Future)", DECADES["future"])]:
                    if decade not in ds.decade.values:
                        html += f'<div class="map-container"><p>No data for {decade_label}</p></div>\n'
                        continue

                    values = ds[metric].sel(decade=decade).values
                    title = decade_label

                    colorbar_label = COLORBAR_LABELS.get(metric, "Value").format(
                        long_name=self.variable_long_name,
                        units=self.variable_units
                    )

                    fig = create_map_figure(
                        lons, lats, values, title,
                        colorscale=COLORSCALES.get(metric, "Viridis"),
                        cmin=cmin, cmax=cmax,
                        colorbar_title=colorbar_label,
                        reversescale=(metric in ("trend",) + tuple(SLOPE_METRICS)
                                      and self.higher_is_worse),
                        hover_format=HOVER_FORMATS.get(metric, DEFAULT_HOVER_FORMAT),
                        subtitle=decade_label,
                    )

                    html += '<div class="map-container">\n'
                    html += f'<div class="map-label">{decade_label}</div>\n'
                    html += fig.to_html(full_html=False, include_plotlyjs=False)
                    html += '</div>\n'
            else:
                # No decade dimension (e.g., trend): show single map
                values = ds[metric].values
                title = METRIC_DESCRIPTIONS.get(metric, metric.title())

                colorbar_label = COLORBAR_LABELS.get(metric, "Value").format(
                    long_name=self.variable_long_name,
                    units=self.variable_units
                )

                fig = create_map_figure(
                    lons, lats, values, title,
                    colorscale=COLORSCALES.get(metric, "Viridis"),
                    cmin=cmin, cmax=cmax,
                    colorbar_title=colorbar_label,
                    reversescale=(metric in ("trend",) + tuple(SLOPE_METRICS)
                                  and self.higher_is_worse),
                    hover_format=HOVER_FORMATS.get(metric, DEFAULT_HOVER_FORMAT),
                )

                html += '<div class="map-container" style="grid-column: span 2;">\n'
                html += f'<div class="map-label">{title} (2020s-2090s)</div>\n'
                html += fig.to_html(full_html=False, include_plotlyjs=False)
                html += '</div>\n'

            html += '</div>\n'
            html += generate_html_footer(timestamp, self.data_source)

            # Write file per scenario
            output_path = output_dir / f"{variable}_{metric}_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

    def _symmetric_limit(self, metric: str, decades: List[int]) -> float:
        """Largest |value| for `metric` over `decades`, across ALL scenarios.

        Diverging panels must be centred on zero with equal blue/red extent, so the limit
        is a single symmetric magnitude rather than a (cmin, cmax) pair. It is pooled
        across scenarios so the same colour means the same rate on every page.
        """
        peak = 0.0
        for ds in self.data.values():
            if metric not in ds.data_vars:
                continue
            for d in decades:
                if d not in ds.decade.values:
                    continue
                v = ds[metric].sel(decade=d).values
                v = np.abs(v[np.isfinite(v)])
                if v.size:
                    peak = max(peak, float(np.max(v) if SYMMETRIC_LIMIT_PERCENTILE is None
                                            else np.percentile(v, SYMMETRIC_LIMIT_PERCENTILE)))
        return peak if peak > 0 else 1.0

    def _change_limit(self, base_decade: int, decades: List[int]) -> float:
        peak = 0.0
        for ds in self.data.values():
            if "median" not in ds.data_vars or base_decade not in ds.decade.values:
                continue
            b = ds["median"].sel(decade=base_decade).values
            for d in decades:
                if d not in ds.decade.values:
                    continue
                diff = ds["median"].sel(decade=d).values - b
                diff = np.abs(diff[np.isfinite(diff)])
                if diff.size:
                    peak = max(peak, float(np.max(diff) if SYMMETRIC_LIMIT_PERCENTILE is None
                                           else np.percentile(diff, SYMMETRIC_LIMIT_PERCENTILE)))
        return peak if peak > 0 else 1.0

    def _clamped_fraction(self, metric: str, decades: List[int], lim: float,
                          base_decade: Optional[int] = None) -> float:
        """Fraction of finite cells whose |value| exceeds the colour limit.

        Those cells are not lost -- Plotly clamps them to the endpoint colour -- but the
        reader deserves to know how much of the map is saturated.
        """
        beyond = total = 0
        for ds in self.data.values():
            if base_decade is not None:
                if "median" not in ds.data_vars or base_decade not in ds.decade.values:
                    continue
                b = ds["median"].sel(decade=base_decade).values
            elif metric not in ds.data_vars:
                continue
            for d in decades:
                if d not in ds.decade.values:
                    continue
                v = (ds["median"].sel(decade=d).values - b if base_decade is not None
                     else ds[metric].sel(decade=d).values)
                v = np.abs(v[np.isfinite(v)])
                beyond += int(np.sum(v > lim))
                total += v.size
        return (beyond / total) if total else 0.0

    def generate_trend_page(self, variable: str, output_dir: Path):
        """One Trend tab per scenario: both slopes AND change-vs-baseline.

        Design notes:
          * Decades are TREND_DECADES (2030s, 2090s), never the 2020s baseline -- its
            slopes are NaN by contract and its change-from-itself is identically 0.
          * Every panel is diverging, centred on 0, spanning +/- max|value| so blue and
            red cover equal magnitude. The limit is shared across scenarios per metric.
          * Showing OLS and Theil-Sen together is the point: they fail in opposite
            regimes, so their disagreement is the signal that a trend is not robust.
        """
        log("Generating trend maps (slopes + change)...")

        first_ds = list(self.data.values())[0]
        lons, lats = first_ds.lon.values, first_ds.lat.values
        base = DECADES["current"]

        available = set()
        for ds in self.data.values():
            available |= set(ds.data_vars)
        slope_metrics = [m for m in SLOPE_METRICS if m in available]
        if not slope_metrics and "trend" in available:
            slope_metrics = ["trend"]        # pre-OUTPUT-SPEC layers

        decades = [d for d in TREND_DECADES if d in first_ds.decade.values]
        limits = {m: self._symmetric_limit(m, decades) for m in slope_metrics}
        change_lim = self._change_limit(base, decades)
        sat = {m: self._clamped_fraction(m, decades, limits[m]) for m in slope_metrics}
        sat["change"] = self._clamped_fraction("median", decades, change_lim,
                                               base_decade=base)
        scale_note = ("95th pct of |value|" if SYMMETRIC_LIMIT_PERCENTILE is not None
                      else "true max|value|")

        # Diverging scales read blue=low / red=high once reversed; for a hazard where more
        # is worse we want red = increase.
        rev = bool(self.higher_is_worse)

        for scenario in self.scenarios:
            if scenario not in self.data:
                continue
            ds = self.data[scenario]
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "trend", scenario,
                                        self.scenarios, self.scenario_labels)

            html += ('<div class="stats-box"><h3>How to read this tab</h3>'
                     '<p style="margin:0;font-size:13px">All panels are centred on zero '
                     f'(white) and span &plusmn;max|value| so blue and red cover equal '
                     'magnitude. <b>OLS</b> and <b>Theil-Sen</b> fail in opposite regimes '
                     '(Theil-Sen collapses to exactly 0 on zero-inflated fields; OLS '
                     'absorbs member level offsets when coverage is uneven), so '
                     '<b>disagreement between them means the trend is not robust</b>. '
                     'The 2020s baseline is omitted: its slopes are NaN by contract and '
                     'its change from itself is zero. Colour limits are the '
                     f'<b>{scale_note}</b>; more extreme cells are clamped to the end '
                     'colour (never blanked) &mdash; the clamped share is noted on each '
                     'panel.</p></div>\n')

            for metric in slope_metrics:
                label = COLORBAR_LABELS.get(metric, "Slope").format(
                    long_name=self.variable_long_name, units=self.variable_units)
                html += f'<h2 style="color:#2c3e50">{METRIC_DESCRIPTIONS.get(metric, metric)}</h2>\n'
                html += '<div class="comparison-grid">\n'
                for d in decades:
                    if metric not in ds.data_vars or d not in ds.decade.values:
                        continue
                    values = ds[metric].sel(decade=d).values
                    lim = limits[metric]
                    fig = create_map_figure(
                        lons, lats, values, f"{d}s",
                        colorscale=COLORSCALES.get(metric, "RdBu"),
                        cmin=-lim, cmax=lim,
                        colorbar_title=label,
                        reversescale=rev,
                        hover_format=HOVER_FORMATS.get(metric, DEFAULT_HOVER_FORMAT),
                        subtitle=(f"{d}s &nbsp;|&nbsp; scale &plusmn;{lim:.3g} "
                                  f"({scale_note}); {100*sat[metric]:.1f}% clamped"),
                    )
                    html += '<div class="map-container">\n'
                    html += f'<div class="map-label">{d}s</div>\n'
                    html += fig.to_html(full_html=False, include_plotlyjs=False)
                    html += '</div>\n'
                html += '</div>\n'

            # Change vs the shared baseline, formerly its own tab.
            if "median" in ds.data_vars and base in ds.decade.values:
                clabel = COLORBAR_LABELS.get("change", "Change").format(
                    long_name=self.variable_long_name, units=self.variable_units)
                html += (f'<h2 style="color:#2c3e50">Absolute change vs {base}s baseline'
                         f'</h2>\n<div class="comparison-grid">\n')
                b = ds["median"].sel(decade=base).values
                for d in decades:
                    if d not in ds.decade.values:
                        continue
                    diff = ds["median"].sel(decade=d).values - b
                    fig = create_map_figure(
                        lons, lats, diff, f"{d}s - {base}s",
                        colorscale=COLORSCALES.get("change", "RdBu"),
                        cmin=-change_lim, cmax=change_lim,
                        colorbar_title=clabel,
                        reversescale=rev,
                        hover_format=DEFAULT_HOVER_FORMAT,
                        subtitle=(f"{d}s &minus; {base}s &nbsp;|&nbsp; scale "
                                  f"&plusmn;{change_lim:.3g} ({scale_note}); "
                                  f"{100*sat['change']:.1f}% clamped"),
                    )
                    html += '<div class="map-container">\n'
                    html += f'<div class="map-label">{d}s &minus; {base}s</div>\n'
                    html += fig.to_html(full_html=False, include_plotlyjs=False)
                    html += '</div>\n'
                html += '</div>\n'

            html += generate_html_footer(timestamp, self.data_source)
            output_path = output_dir / f"{variable}_trend_{scenario}.html"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log(f"  Saved: {output_path.name}")

    def generate_members_page(self, variable: str, output_dir: Path):
        """Per-member (LSM x GCM) raw-value maps on ONE shared colour scale.

        Reads an optional `{variable}_members.nc` sitting beside the processed scenario
        files, with dims (member, lat, lon). Processors emit it as a diagnostic; when it
        is absent the tab is skipped rather than failing, so pre-existing layers still
        render. A shared scale is the whole point -- it makes a member that runs hot, or
        one whose spatial distribution is unlike its siblings, visible at a glance.
        """
        log("Generating per-member maps...")
        path = self.processed_dir / f"{variable}_members.nc"
        if not path.exists():
            log(f"  SKIP: {path.name} not found "
                f"(processor did not emit per-member diagnostics)")
            return

        ds = xr.open_dataset(path)
        var = "value" if "value" in ds.data_vars else list(ds.data_vars)[0]
        members = [str(m) for m in ds["member"].values]
        lons, lats = ds.lon.values, ds.lat.values
        vals = ds[var].values                       # (member, lat, lon)

        # This tab is 22 panels on one page -- at full 0.5 deg it was 37.6 MB. Block-mean
        # to a coarser grid; the question it answers (do members differ in level or
        # distribution?) does not need per-cell resolution.
        stride = MEMBERS_GRID_STRIDE
        if stride > 1:
            vals = np.stack([block_mean(vals[i], stride) for i in range(vals.shape[0])])
            ny = lats.size - lats.size % stride
            nx = lons.size - lons.size % stride
            lats = lats[:ny].reshape(-1, stride).mean(axis=1)
            lons = lons[:nx].reshape(-1, stride).mean(axis=1)
            log(f"  block-averaged to stride {stride} -> "
                f"{vals.shape[1]}x{vals.shape[2]} grid ({lats.size}x{lons.size} coords)")

        finite = vals[np.isfinite(vals)]
        cmin = float(np.percentile(finite, 2)) if finite.size else 0.0
        cmax = float(np.percentile(finite, 98)) if finite.size else 1.0

        label = COLORBAR_LABELS.get("median", "{long_name} [{units}]").format(
            long_name=self.variable_long_name, units=self.variable_units)

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        scen0 = self.scenarios[0] if self.scenarios else ""
        html = generate_html_header(variable, "members", scen0,
                                    self.scenarios, self.scenario_labels)
        html += (f'<div class="stats-box"><h3>{len(members)} ensemble members '
                 f'&mdash; shared colour scale [{cmin:.3g}, {cmax:.3g}]</h3>'
                 '<p style="margin:0;font-size:13px">Every LSM &times; GCM member on the '
                 'SAME scale, so a member that runs systematically hot/cold, or whose '
                 'spatial distribution differs from its siblings, is visible immediately. '
                 f'Field: {ds.attrs.get("member_field", "per-member mean")}.</p></div>\n')

        # Quantified companion to the eyeball test.
        html += ('<table style="width:100%;border-collapse:collapse;background:white;'
                 'margin-bottom:20px;font-size:13px">'
                 '<tr style="background:#34495e;color:white">'
                 '<th style="padding:8px;text-align:left">Member</th>'
                 '<th style="padding:8px">Land mean</th><th style="padding:8px">Median</th>'
                 '<th style="padding:8px">P95</th><th style="padding:8px">Max</th>'
                 '<th style="padding:8px">Zero %</th><th style="padding:8px">Cells</th></tr>\n')
        for i, m in enumerate(members):
            v = vals[i][np.isfinite(vals[i])]
            if not v.size:
                continue
            html += (f'<tr><td style="padding:6px;border-bottom:1px solid #ecf0f1">{m}</td>'
                     f'<td style="padding:6px;text-align:center">{np.mean(v):.4g}</td>'
                     f'<td style="padding:6px;text-align:center">{np.median(v):.4g}</td>'
                     f'<td style="padding:6px;text-align:center">{np.percentile(v,95):.4g}</td>'
                     f'<td style="padding:6px;text-align:center">{np.max(v):.4g}</td>'
                     f'<td style="padding:6px;text-align:center">{100*np.mean(v==0):.1f}</td>'
                     f'<td style="padding:6px;text-align:center">{v.size:,}</td></tr>\n')
        html += '</table>\n'

        html += '<div class="comparison-grid">\n'
        for i, m in enumerate(members):
            fig = create_map_figure(
                lons, lats, vals[i], m,
                colorscale=COLORSCALES.get("median", "Viridis"),
                cmin=cmin, cmax=cmax,
                colorbar_title=label,
                hover_format=DEFAULT_HOVER_FORMAT,
                subtitle=m,
            )
            html += '<div class="map-container">\n'
            html += f'<div class="map-label">{m}</div>\n'
            html += fig.to_html(full_html=False, include_plotlyjs=False)
            html += '</div>\n'
        html += '</div>\n'
        html += generate_html_footer(timestamp, self.data_source)

        output_path = output_dir / f"{variable}_members.html"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        ds.close()
        log(f"  Saved: {output_path.name} ({len(members)} members)")

    def generate_confidence_comparison(self, variable: str, output_dir: Path):
        """Generate confidence interval comparison maps (per-scenario files)."""
        log("Generating confidence interval comparison maps...")

        first_ds = list(self.data.values())[0]
        lons = first_ds.lon.values
        lats = first_ds.lat.values

        # Generate separate file for each scenario
        for scenario in self.scenarios:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "confidence", scenario, self.scenarios, self.scenario_labels)

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

            html += generate_html_footer(timestamp, self.data_source)

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
        for scenario in self.scenarios:
            if scenario not in self.data:
                continue

            ds = self.data[scenario]
            if DECADES["current"] not in ds.decade.values or DECADES["future"] not in ds.decade.values:
                continue

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "change", scenario, self.scenarios, self.scenario_labels)

            pos_color, neg_color = ("red", "blue") if self.higher_is_worse else ("blue", "red")
            html += '<div class="stats-box"><h3>Absolute Change: 2090s minus 2020s</h3>'
            html += (f'<p>Positive values ({pos_color}) indicate increase, '
                     f'negative values ({neg_color}) indicate decrease.</p></div>\n')

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
                colorbar_title=colorbar_label,
                reversescale=self.higher_is_worse
            )

            # Single centered map for change
            html += '<div style="max-width: 800px; margin: 0 auto;">\n'
            html += '<div class="map-container">\n'
            html += f'<div class="map-label">{title}</div>\n'
            html += fig.to_html(full_html=False, include_plotlyjs=False)
            html += '</div>\n'
            html += '</div>\n'

            html += generate_html_footer(timestamp, self.data_source)

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
        for scenario in self.scenarios:
            if scenario not in self.data or scenario not in self.baseline_stats:
                continue

            ds = self.data[scenario]
            mean_2020s, std_2020s = self.baseline_stats[scenario]

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            html = generate_html_header(variable, "anomaly", scenario, self.scenarios, self.scenario_labels)

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

            html += generate_html_footer(timestamp, self.data_source)

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
            ("trend", "Trend",
             "OLS &amp; Theil-Sen slopes and change vs baseline (2030s, 2090s)"),
            ("confidence", "Confidence Intervals", "Lower (25th) and upper (75th) bounds"),
            ("anomaly", "Anomaly Detection", f"Values >{ANOMALY_SIGMA}σ from 2020s mean"),
            ("members", "Members", "Every LSM × GCM member on a shared scale"),
        ]

        # Build dynamic CSS for scenario button colors
        btn_css = ""
        for scenario in self.scenarios:
            color = self.scenario_colors.get(scenario, "#3498db")
            # Darken color for hover
            btn_css += f"""        a.btn.{scenario} {{ background: {color}; }}
        a.btn.{scenario}:hover {{ filter: brightness(0.85); }}
"""

        # Build legend items dynamically
        legend_items = ""
        for scenario in self.scenarios:
            color = self.scenario_colors.get(scenario, "#3498db")
            label = self.scenario_labels.get(scenario, scenario)
            legend_items += f"""        <div class="legend-item">
            <span class="legend-color" style="background: {color};"></span>
            <strong>{label.split()[0]}</strong> - {' '.join(label.split()[1:]) if len(label.split()) > 1 else ''}
        </div>
"""

        # Build table header with scenario columns
        scenario_headers = ""
        for scenario in self.scenarios:
            label = self.scenario_labels.get(scenario, scenario)
            short_label = label.split()[0]  # e.g., "SSP1-2.6" or "RCP2.6"
            scenario_headers += f'            <th>{short_label}</th>\n'

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
        a.btn {{ display: inline-block; background: #3498db; color: white; padding: 8px 16px;
                 border-radius: 4px; text-decoration: none; font-size: 13px; }}
        a.btn:hover {{ background: #2980b9; }}
{btn_css}        .footer {{ text-align: center; color: #666; padding: 20px; font-size: 12px; }}
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
        <h3>Scenario Legend ({self.data_source})</h3>
{legend_items}    </div>

    <table>
        <tr>
            <th class="metric">Metric</th>
{scenario_headers}        </tr>
"""

        for metric_key, metric_name, metric_desc in metrics:
            # Build scenario cells dynamically
            scenario_cells = ""
            if metric_key in SCENARIO_INDEPENDENT_TABS:
                scenario_cells = (f'            <td colspan="{len(self.scenarios)}">'
                                  f'<a href="{variable}_{metric_key}.html" class="btn">'
                                  f'View (all members)</a></td>\n')
            else:
                for scenario in self.scenarios:
                    scenario_cells += f'            <td><a href="{variable}_{metric_key}_{scenario}.html" class="btn {scenario}">View</a></td>\n'

            html += f"""        <tr>
            <td class="metric-name">{metric_name}<br><span style="font-weight:normal;font-size:11px;color:#888;">{metric_desc}</span></td>
{scenario_cells}        </tr>
"""

        html += f"""    </table>

    <div class="footer">
        Generated: {timestamp}<br>
        Data source: {self.data_source} | Units: {self.variable_units}
    </div>
</body>
</html>
"""

        output_path = output_dir / "index.html"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        log(f"  Saved: {output_path.name}")


def main():
    """Main entry point.

    Usage:
        python generate_maps.py [variable] [processed_dir] [output_dir]

    Examples:
        python generate_maps.py                          # Default: qg from data/processed
        python generate_maps.py leh                      # Generate maps for leh variable
        python generate_maps.py leh ./outputs/processed  # Custom processed directory
        python generate_maps.py leh ./outputs/processed ./reports/maps  # Full custom paths
    """
    import sys

    # Get project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent

    # Parse command-line arguments
    variable = sys.argv[1] if len(sys.argv) > 1 else "qg"
    processed_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else project_root / "data" / "processed"
    output_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else project_root / "reports" / "maps"

    log("=" * 60)
    log("Generating Synoptic and Diagnostic Maps")
    log("=" * 60)
    log(f"Variable: {variable}")
    log(f"Processed data: {processed_dir}")
    log(f"Output directory: {output_dir}")
    log("=" * 60)

    generator = MapCollectionGenerator(processed_dir, output_dir)

    # Generate maps for specified variable
    generator.generate_all_collections(variable)

    log("\n" + "=" * 60)
    log("Map generation complete!")
    log(f"Output directory: {output_dir / variable}")
    log("=" * 60)


if __name__ == "__main__":
    main()
