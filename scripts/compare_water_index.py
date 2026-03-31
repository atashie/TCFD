"""Interactive HTML comparison: Old (RCP) vs New (SSP) water index files.

Generates a trend-focused HTML report with:
- Per-pixel Theil-Sen slope maps (Old | New | Difference) showing decadal trends
- Spatial Spearman R² of trend patterns between old and new
- Trend sign agreement maps
- Reduced raw-value comparison maps (Annual_Mean at key decades only)
- Fixed scatter plots with percentile clipping
- Seasonal cycle comparison panels

Memory strategy: loads one (value_type, decade, scenario) slice at a time.

Usage:
    python scripts/compare_water_index.py
    python scripts/compare_water_index.py --old old.nc --new new.nc --output report.html
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats as sp_stats

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_water_variables import get_variable_config, SHARED_CONFIG, VALUE_TYPE_NAMES
from utils.water_index_compare import (
    spatial_correlation, summary_statistics, coverage_comparison,
    detect_anomalous_differences, ks_test, rmsd, compare_seasonal_cycle,
)


# Scenario pair mapping (SSP -> RCP for comparison)
SCENARIO_PAIRS = [
    ("ssp126", "rcp26", "SSP1-2.6 vs RCP2.6"),
    ("ssp370", "rcp60", "SSP3-7.0 vs RCP6.0"),
    ("ssp585", "rcp85", "SSP5-8.5 vs RCP8.5"),
]

# Reference locations for seasonal cycle comparison
REFERENCE_LOCATIONS = [
    ("Amazon Basin", -3.0, -60.0),
    ("Sahara Desert", 23.0, 5.0),
    ("Central US", 40.0, -95.0),
    ("Ganges Delta", 23.0, 89.0),
    ("Greenland", 72.0, -40.0),
    ("Congo Basin", 0.0, 22.0),
    ("Western Europe", 48.0, 5.0),
    ("Siberia", 60.0, 90.0),
    ("Australia", -25.0, 135.0),
    ("Southeast Asia", 15.0, 105.0),
]

# Rendering subsample step (keep file size manageable)
SPATIAL_STEP = 4


def log(msg: str):
    print(msg, flush=True)


# Module-level variable name, set in main() before generate_report()
_VARIABLE = "tws"


def load_slice(ds: xr.Dataset, scenario: str, value_type: int, decade: int) -> np.ndarray:
    """Load a single 2D slice (lat, lon) from a water index file."""
    try:
        return ds[_VARIABLE].sel(scenario=scenario, value_type=value_type, decade=decade).values
    except Exception:
        return np.full((360, 720), np.nan)


def find_common_decades(old_ds: xr.Dataset, new_ds: xr.Dataset) -> List[int]:
    """Find decades present in both files."""
    old_decades = set(int(d) for d in old_ds.decade.values)
    new_decades = set(int(d) for d in new_ds.decade.values)
    return sorted(old_decades & new_decades)


# ---------------------------------------------------------------------------
# Trend computation
# ---------------------------------------------------------------------------

def compute_pixel_trends(
    ds: xr.Dataset, scenario: str, decades: List[int], vt: int = 12,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute per-pixel Theil-Sen slope and Spearman rho across decades.

    Args:
        ds: Opened dataset.
        scenario: Scenario name to select.
        decades: List of decade values to use for trend.
        vt: Value type index (default 12 = Annual_Mean).

    Returns:
        (slope, rho) — both (lat, lon) arrays.
        slope: Theil-Sen slope (value change per decade).
        rho: Spearman rank correlation coefficient of value vs time.
    """
    n = len(decades)
    x = np.array(decades, dtype=float)

    # Stack: (lat, lon, n_decades)
    stack = np.stack([load_slice(ds, scenario, vt, d) for d in decades], axis=-1)

    nlat, nlon, _ = stack.shape
    slope = np.full((nlat, nlon), np.nan)
    rho = np.full((nlat, nlon), np.nan)

    # Theil-Sen: median of all pairwise slopes
    pair_slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            dy = stack[:, :, j] - stack[:, :, i]
            pair_slopes.append(dy / dx * 10.0)  # per decade

    # (lat, lon, n_pairs)
    all_slopes = np.stack(pair_slopes, axis=-1)
    slope = np.nanmedian(all_slopes, axis=-1)

    # Spearman rho per pixel (vectorized via ranking)
    # Rank along time axis
    # Use scipy for a random sample first; for full grid do manual ranking
    ranks_x = sp_stats.rankdata(x)  # same for all pixels
    # Rank y along axis=-1 per pixel
    # argsort of argsort gives ranks
    order = np.argsort(np.argsort(stack, axis=-1), axis=-1).astype(float) + 1.0
    # Where stack is NaN, set rank to NaN
    nan_mask = np.isnan(stack)
    order[nan_mask] = np.nan

    # Spearman = Pearson on ranks
    n_valid = np.sum(~nan_mask, axis=-1)
    mean_rx = np.mean(ranks_x)
    mean_ry = np.nanmean(order, axis=-1, keepdims=True)

    dx_r = ranks_x[None, None, :] - mean_rx
    dy_r = order - mean_ry
    dy_r[nan_mask] = 0.0  # zero contribution from NaN positions

    num = np.nansum(dx_r * dy_r, axis=-1)
    denom_x = np.sqrt(np.sum(dx_r**2))  # same for all pixels
    # For dy, only use non-NaN
    denom_y = np.sqrt(np.nansum(dy_r**2, axis=-1))

    valid_denom = (denom_x * denom_y) > 0
    rho[valid_denom] = num[valid_denom] / (denom_x * denom_y[valid_denom])

    # Mask pixels with too few valid decades
    too_few = n_valid < 4
    slope[too_few] = np.nan
    rho[too_few] = np.nan

    return slope, rho


# ---------------------------------------------------------------------------
# Visualization helpers
# ---------------------------------------------------------------------------

def _subsample(arr2d: np.ndarray, step: int = SPATIAL_STEP) -> np.ndarray:
    sub = arr2d[::step, ::step]
    # Round to 4 significant figures instead of fixed decimal places,
    # so small values (e.g., 1e-6 kg/m2/s) are not truncated to zero.
    valid = sub[~np.isnan(sub)]
    if len(valid) == 0:
        return sub
    max_abs = np.max(np.abs(valid))
    if max_abs == 0:
        return sub
    # Number of decimal places to keep 4 significant figures
    decimals = max(0, int(4 - np.floor(np.log10(max_abs))) - 1)
    return np.round(sub, decimals)


def _to_list(arr2d: np.ndarray):
    """Convert 2D array to nested Python list with None for NaN (Plotly JSON compat)."""
    return np.where(np.isnan(arr2d), None, arr2d).tolist()


def create_three_panel_map(
    old_2d: np.ndarray,
    new_2d: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    title: str,
    stats: Dict,
    diff_label: str = "Difference (New - Old)",
    colorscale: str = "Viridis",
    diff_colorscale: str = "RdBu",
) -> go.Figure:
    """Create three-panel map: Old | New | Difference."""
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=["Old (RCP)", "New (SSP)", diff_label],
        horizontal_spacing=0.03,
    )

    # Shared colorscale bounds for Old/New (clip outliers)
    combined = np.concatenate([old_2d[~np.isnan(old_2d)], new_2d[~np.isnan(new_2d)]])
    if len(combined) > 0:
        vmin = float(np.percentile(combined, 2))
        vmax = float(np.percentile(combined, 98))
    else:
        vmin, vmax = 0, 1

    lat_sub = lats[::SPATIAL_STEP].tolist()
    lon_sub = lons[::SPATIAL_STEP].tolist()
    old_sub = _to_list(_subsample(old_2d))
    new_sub = _to_list(_subsample(new_2d))

    fig.add_trace(go.Heatmap(
        z=old_sub, x=lon_sub, y=lat_sub,
        colorscale=colorscale, zmin=vmin, zmax=vmax,
        colorbar=dict(x=0.3, len=0.8, title="Value"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>value: %{z:.4g}<extra>Old</extra>",
    ), row=1, col=1)

    fig.add_trace(go.Heatmap(
        z=new_sub, x=lon_sub, y=lat_sub,
        colorscale=colorscale, zmin=vmin, zmax=vmax,
        colorbar=dict(x=0.63, len=0.8, title="Value"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>value: %{z:.4g}<extra>New</extra>",
    ), row=1, col=2)

    # Difference panel (RdBu, centered on zero per GUARDRAILS rule 5)
    diff = new_2d - old_2d
    diff_valid = diff[~np.isnan(diff)]
    max_abs = float(np.percentile(np.abs(diff_valid), 98)) if len(diff_valid) > 0 else 1.0
    diff_sub = _to_list(_subsample(diff))

    fig.add_trace(go.Heatmap(
        z=diff_sub, x=lon_sub, y=lat_sub,
        colorscale=diff_colorscale, zmin=-max_abs, zmax=max_abs, zmid=0,
        colorbar=dict(x=0.98, len=0.8, title="Diff"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>diff: %{z:.4g}<extra>Diff</extra>",
    ), row=1, col=3)

    stats_text = (
        f"Pearson r={stats.get('pearson_r', 0):.3f} | "
        f"Spearman r={stats.get('spearman_r', 0):.3f} | "
        f"RMSD={stats.get('rmsd', 0):.4g} | "
        f"Coverage overlap={stats.get('overlap_pct', 0):.1f}%"
    )

    fig.update_layout(
        title=dict(text=f"{title}<br><sub>{stats_text}</sub>", x=0.5),
        height=400, width=1800,
        margin=dict(l=40, r=40, t=80, b=40),
    )
    return fig


def create_trend_map(
    old_trend: np.ndarray,
    new_trend: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    title: str,
    trend_corr: Dict,
) -> go.Figure:
    """Create three-panel trend map: Old slope | New slope | Difference.

    Uses RdBu diverging colorscale centered on zero for all three panels,
    since trends can be positive or negative.
    """
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=[
            "Old (RCP) — Sen Slope/decade",
            "New (SSP) — Sen Slope/decade",
            "Slope Difference (New - Old)",
        ],
        horizontal_spacing=0.03,
    )

    # Shared scale for old/new slope panels
    combined = np.concatenate([
        old_trend[~np.isnan(old_trend)],
        new_trend[~np.isnan(new_trend)],
    ])
    if len(combined) > 0:
        max_abs = float(np.percentile(np.abs(combined), 98))
    else:
        max_abs = 1.0

    lat_sub = lats[::SPATIAL_STEP].tolist()
    lon_sub = lons[::SPATIAL_STEP].tolist()

    fig.add_trace(go.Heatmap(
        z=_to_list(_subsample(old_trend)), x=lon_sub, y=lat_sub,
        colorscale="RdBu_r", zmin=-max_abs, zmax=max_abs, zmid=0,
        colorbar=dict(x=0.3, len=0.8, title="Slope"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>slope: %{z:.4g}<extra>Old</extra>",
    ), row=1, col=1)

    fig.add_trace(go.Heatmap(
        z=_to_list(_subsample(new_trend)), x=lon_sub, y=lat_sub,
        colorscale="RdBu_r", zmin=-max_abs, zmax=max_abs, zmid=0,
        colorbar=dict(x=0.63, len=0.8, title="Slope"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>slope: %{z:.4g}<extra>New</extra>",
    ), row=1, col=2)

    # Difference
    diff = new_trend - old_trend
    diff_valid = diff[~np.isnan(diff)]
    diff_max = float(np.percentile(np.abs(diff_valid), 98)) if len(diff_valid) > 0 else 1.0

    fig.add_trace(go.Heatmap(
        z=_to_list(_subsample(diff)), x=lon_sub, y=lat_sub,
        colorscale="RdBu_r", zmin=-diff_max, zmax=diff_max, zmid=0,
        colorbar=dict(x=0.98, len=0.8, title="Diff"),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>diff: %{z:.4g}<extra>Diff</extra>",
    ), row=1, col=3)

    stats_text = (
        f"Spatial Pearson r={trend_corr.get('pearson_r', 0):.3f} | "
        f"Spearman r={trend_corr.get('spearman_r', 0):.3f} | "
        f"Spearman R²={trend_corr.get('spearman_r', 0)**2:.3f}"
    )

    fig.update_layout(
        title=dict(text=f"{title}<br><sub>{stats_text}</sub>", x=0.5),
        height=400, width=1800,
        margin=dict(l=40, r=40, t=80, b=40),
    )
    return fig


def create_trend_agreement_map(
    old_trend: np.ndarray,
    new_trend: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    title: str,
) -> go.Figure:
    """Map showing where old and new trends agree in sign.

    Color coding:
    - Both positive (wetting): blue
    - Both negative (drying): red
    - Disagreement: grey
    - Insufficient data: white/NaN
    """
    mask = ~np.isnan(old_trend) & ~np.isnan(new_trend)
    agreement = np.full_like(old_trend, np.nan)

    both_pos = mask & (old_trend > 0) & (new_trend > 0)
    both_neg = mask & (old_trend < 0) & (new_trend < 0)
    disagree = mask & ~both_pos & ~both_neg

    agreement[both_pos] = 1.0    # Both wetting
    agreement[both_neg] = -1.0   # Both drying
    agreement[disagree] = 0.0    # Disagree

    n_agree = int(both_pos.sum() + both_neg.sum())
    n_disagree = int(disagree.sum())
    n_total = int(mask.sum())
    pct = n_agree / n_total * 100 if n_total > 0 else 0

    lat_sub = lats[::SPATIAL_STEP].tolist()
    lon_sub = lons[::SPATIAL_STEP].tolist()

    # Custom discrete colorscale
    colorscale = [
        [0.0, "#d62728"],    # -1: both drying (red)
        [0.25, "#d62728"],
        [0.25, "#999999"],   # 0: disagree (grey)
        [0.75, "#999999"],
        [0.75, "#1f77b4"],   # +1: both wetting (blue)
        [1.0, "#1f77b4"],
    ]

    fig = go.Figure(go.Heatmap(
        z=_to_list(_subsample(agreement)), x=lon_sub, y=lat_sub,
        colorscale=colorscale, zmin=-1.5, zmax=1.5,
        colorbar=dict(
            title="Agreement",
            tickvals=[-1, 0, 1],
            ticktext=["Both Drying", "Disagree", "Both Wetting"],
        ),
        hovertemplate="lat: %{y:.1f}<br>lon: %{x:.1f}<br>%{z}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=f"{title}<br><sub>Sign agreement: {n_agree:,}/{n_total:,} pixels ({pct:.1f}%) | "
                 f"Disagree: {n_disagree:,} ({100-pct:.1f}%)</sub>",
            x=0.5,
        ),
        height=450, width=900,
        margin=dict(l=40, r=40, t=80, b=40),
    )
    return fig


def create_trend_scatter(
    old_trend: np.ndarray,
    new_trend: np.ndarray,
    title: str,
) -> go.Figure:
    """Scatter plot of old vs new per-pixel trends."""
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Old Trend vs New Trend", "Trend Difference Distribution"],
    )

    mask = ~np.isnan(old_trend) & ~np.isnan(new_trend)
    if mask.sum() < 10:
        fig.update_layout(title=f"{title} - Insufficient data")
        return fig

    o, n_ = old_trend[mask], new_trend[mask]

    # Clip to 2nd-98th percentile for visualization
    o_lo, o_hi = np.percentile(o, 2), np.percentile(o, 98)
    n_lo, n_hi = np.percentile(n_, 2), np.percentile(n_, 98)
    clip_mask = (o >= o_lo) & (o <= o_hi) & (n_ >= n_lo) & (n_ <= n_hi)
    o_clip, n_clip = o[clip_mask], n_[clip_mask]

    # Subsample for rendering
    max_pts = 25_000
    if len(o_clip) > max_pts:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(o_clip), max_pts, replace=False)
        o_s, n_s = o_clip[idx], n_clip[idx]
    else:
        o_s, n_s = o_clip, n_clip

    # Scatter
    fig.add_trace(go.Scattergl(
        x=o_s.tolist(), y=n_s.tolist(), mode="markers",
        marker=dict(size=2, opacity=0.3, color="#3498db"),
        name="Pixels",
    ), row=1, col=1)

    # 1:1 line
    vmin = min(float(o_s.min()), float(n_s.min()))
    vmax = max(float(o_s.max()), float(n_s.max()))
    fig.add_trace(go.Scatter(
        x=[vmin, vmax], y=[vmin, vmax], mode="lines",
        line=dict(color="red", dash="dash"), name="1:1",
    ), row=1, col=1)

    # OLS fit line
    slope_fit, intercept_fit = np.polyfit(o_s, n_s, 1)
    fig.add_trace(go.Scatter(
        x=[vmin, vmax],
        y=[intercept_fit + slope_fit * vmin, intercept_fit + slope_fit * vmax],
        mode="lines", line=dict(color="green", dash="dot"),
        name=f"OLS (slope={slope_fit:.2f})",
    ), row=1, col=1)

    fig.update_xaxes(title_text="Old (RCP) Trend/decade", row=1, col=1)
    fig.update_yaxes(title_text="New (SSP) Trend/decade", row=1, col=1)

    # Histogram of trend differences
    diff = n_ - o
    diff_clip = diff[(diff >= np.percentile(diff, 1)) & (diff <= np.percentile(diff, 99))]
    fig.add_trace(go.Histogram(
        x=diff_clip.tolist(), nbinsx=100, name="Trend Diff",
        marker_color="#2ecc71", opacity=0.7,
    ), row=1, col=2)
    fig.update_xaxes(title_text="New Trend - Old Trend", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)

    corr = spatial_correlation(old_trend, new_trend)
    fig.update_layout(
        title=dict(
            text=f"{title}<br><sub>Pearson r={corr['pearson_r']:.3f}, "
                 f"Spearman r={corr['spearman_r']:.3f}, "
                 f"R²={corr['spearman_r']**2:.3f}, "
                 f"n={corr['n_valid']:,}</sub>",
            x=0.5,
        ),
        height=450, width=1200,
        margin=dict(l=60, r=40, t=80, b=40),
    )
    return fig


def create_scatter_histogram(
    old_2d: np.ndarray, new_2d: np.ndarray, title: str,
) -> go.Figure:
    """Scatter plot of old vs new raw values, with percentile clipping."""
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Old vs New (scatter)", "Difference Distribution"],
    )

    mask = ~np.isnan(old_2d) & ~np.isnan(new_2d)
    if mask.sum() < 10:
        fig.update_layout(title=f"{title} - Insufficient data")
        return fig

    o, n_ = old_2d[mask], new_2d[mask]

    # Clip to 2nd-98th percentile to avoid outlier compression
    o_lo, o_hi = np.percentile(o, 2), np.percentile(o, 98)
    n_lo, n_hi = np.percentile(n_, 2), np.percentile(n_, 98)
    clip_mask = (o >= o_lo) & (o <= o_hi) & (n_ >= n_lo) & (n_ <= n_hi)
    o_clip, n_clip = o[clip_mask], n_[clip_mask]

    max_pts = 25_000
    if len(o_clip) > max_pts:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(o_clip), max_pts, replace=False)
        o_s, n_s = o_clip[idx], n_clip[idx]
    else:
        o_s, n_s = o_clip, n_clip

    fig.add_trace(go.Scattergl(
        x=o_s.tolist(), y=n_s.tolist(), mode="markers",
        marker=dict(size=2, opacity=0.3, color="#3498db"),
        name="Cells",
    ), row=1, col=1)

    vmin = min(float(o_s.min()), float(n_s.min()))
    vmax = max(float(o_s.max()), float(n_s.max()))
    fig.add_trace(go.Scatter(
        x=[vmin, vmax], y=[vmin, vmax], mode="lines",
        line=dict(color="red", dash="dash"), name="1:1 line",
    ), row=1, col=1)
    fig.update_xaxes(title_text="Old (RCP)", row=1, col=1)
    fig.update_yaxes(title_text="New (SSP)", row=1, col=1)

    diff = n_ - o
    diff_clip = diff[(diff >= np.percentile(diff, 1)) & (diff <= np.percentile(diff, 99))]
    fig.add_trace(go.Histogram(
        x=diff_clip.tolist(), nbinsx=100, name="Difference",
        marker_color="#2ecc71", opacity=0.7,
    ), row=1, col=2)
    fig.update_xaxes(title_text="New - Old", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)

    corr = spatial_correlation(old_2d, new_2d)
    fig.update_layout(
        title=dict(
            text=f"{title}<br><sub>Pearson r={corr['pearson_r']:.3f}, "
                 f"Spearman r={corr['spearman_r']:.3f}, "
                 f"n={corr['n_valid']:,} (clipped to P2-P98)</sub>",
            x=0.5,
        ),
        height=400, width=1200,
        margin=dict(l=60, r=40, t=80, b=40),
    )
    return fig


def create_seasonal_comparison(
    old_ds: xr.Dataset, new_ds: xr.Dataset,
    old_scenario: str, new_scenario: str,
    decade: int,
) -> go.Figure:
    """Create seasonal cycle comparison for reference locations."""
    n_locs = len(REFERENCE_LOCATIONS)
    rows = (n_locs + 1) // 2
    fig = make_subplots(
        rows=rows, cols=2,
        subplot_titles=[name for name, _, _ in REFERENCE_LOCATIONS],
        vertical_spacing=0.08,
    )

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for idx, (name, lat, lon) in enumerate(REFERENCE_LOCATIONS):
        row = idx // 2 + 1
        col = idx % 2 + 1

        old_monthly = np.array([
            float(load_slice(old_ds, old_scenario, m, decade)
                  [np.argmin(np.abs(old_ds.lat.values - lat)),
                   np.argmin(np.abs(old_ds.lon.values - lon))])
            for m in range(12)
        ])
        new_monthly = np.array([
            float(load_slice(new_ds, new_scenario, m, decade)
                  [np.argmin(np.abs(new_ds.lat.values - lat)),
                   np.argmin(np.abs(new_ds.lon.values - lon))])
            for m in range(12)
        ])

        fig.add_trace(go.Scatter(
            x=month_names, y=old_monthly.tolist(), mode="lines+markers",
            name=f"Old ({old_scenario})", line=dict(color="#e74c3c"),
            legendgroup="old", showlegend=(idx == 0),
        ), row=row, col=col)
        fig.add_trace(go.Scatter(
            x=month_names, y=new_monthly.tolist(), mode="lines+markers",
            name=f"New ({new_scenario})", line=dict(color="#3498db"),
            legendgroup="new", showlegend=(idx == 0),
        ), row=row, col=col)

    fig.update_layout(
        title=f"Seasonal Cycle Comparison - Decade {decade}s",
        height=300 * rows, width=1200,
        margin=dict(l=60, r=40, t=80, b=40),
    )
    return fig


# ---------------------------------------------------------------------------
# Statistics table
# ---------------------------------------------------------------------------

def build_trend_stats_table(
    trend_results: List[Dict],
) -> str:
    """Build HTML table summarizing trend correlation across scenario pairs."""
    rows_html = []
    for r in trend_results:
        r2 = r["spearman_r"] ** 2
        rows_html.append(
            f"<tr><td>{r['label']}</td>"
            f"<td>{r['pearson_r']:.3f}</td>"
            f"<td>{r['spearman_r']:.3f}</td>"
            f"<td style='font-weight:bold'>{r2:.3f}</td>"
            f"<td>{r['sign_agree_pct']:.1f}%</td>"
            f"<td>{r['n_valid']:,}</td></tr>"
        )

    return (
        '<table border="1" cellpadding="6" cellspacing="0" '
        'style="border-collapse: collapse; font-size: 13px; margin: 10px 0;">'
        "<tr style='background:#f0f0f0'><th>Scenario Pair</th>"
        "<th>Pearson r</th><th>Spearman r</th><th>Spearman R²</th>"
        "<th>Trend Sign Agreement</th><th>Valid Pixels</th></tr>"
        + "\n".join(rows_html)
        + "</table>"
    )


def build_value_stats_table(
    old_ds: xr.Dataset, new_ds: xr.Dataset, decade: int,
) -> str:
    """Build compact statistics table for raw value comparison at one decade."""
    rows_html = []
    for ssp, rcp, label in SCENARIO_PAIRS:
        for vt in [12, 16]:  # Annual_Mean, Q50
            vt_name = VALUE_TYPE_NAMES[vt]
            old_2d = load_slice(old_ds, rcp, vt, decade)
            new_2d = load_slice(new_ds, ssp, vt, decade)

            corr = spatial_correlation(old_2d, new_2d)
            cov = coverage_comparison(old_2d, new_2d)
            rmsd_val = rmsd(old_2d, new_2d)

            rows_html.append(
                f"<tr><td>{label}</td><td>{vt_name}</td>"
                f"<td>{corr['pearson_r']:.3f}</td>"
                f"<td>{corr['spearman_r']:.3f}</td>"
                f"<td>{rmsd_val:.4g}</td>"
                f"<td>{cov['overlap_pct']:.1f}%</td></tr>"
            )

    return (
        '<table border="1" cellpadding="6" cellspacing="0" '
        'style="border-collapse: collapse; font-size: 13px; margin: 10px 0;">'
        f"<tr style='background:#f0f0f0'><th>Scenario Pair</th><th>Value Type</th>"
        f"<th>Pearson r</th><th>Spearman r</th><th>RMSD</th><th>Coverage Overlap</th></tr>"
        + "\n".join(rows_html)
        + "</table>"
    )


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(old_path: Path, new_path: Path, output_path: Path):
    """Generate the full HTML comparison report."""
    log(f"Old file: {old_path}")
    log(f"New file: {new_path}")

    old_ds = xr.open_dataset(old_path)
    new_ds = xr.open_dataset(new_path)

    common_decades = find_common_decades(old_ds, new_ds)
    log(f"Common decades: {common_decades}")

    # Use 2020-2090 for trend computation (skip 2010 — different data sources)
    trend_decades = [d for d in common_decades if d >= 2020]
    log(f"Trend decades: {trend_decades}")

    lats = new_ds.lat.values
    lons = new_ds.lon.values

    html_parts = []
    plotly_included = False

    def add_fig(fig):
        nonlocal plotly_included
        include_js = "cdn" if not plotly_included else False
        plotly_included = True
        html_parts.append(fig.to_html(full_html=False, include_plotlyjs=include_js))

    # --- HTML Header ---
    html_parts.append(f"""
    <html>
    <head>
        <title>Water Index Comparison: {_VARIABLE} (RCP vs SSP)</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; background: #fafafa; }}
            h1 {{ color: #2c3e50; }}
            h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; }}
            h3 {{ color: #555; }}
            .section {{ margin: 30px 0; padding: 20px; background: white;
                        border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            .note {{ color: #666; font-style: italic; margin: 10px 0; }}
            .highlight {{ background: #ffffcc; padding: 4px 8px; border-radius: 4px; }}
        </style>
    </head>
    <body>
    <h1>Water Index Comparison: {_VARIABLE.upper()} &mdash; RCP (Old) vs SSP (New)</h1>
    <p>Old file: <code>{old_path.name}</code> | New file: <code>{new_path.name}</code></p>
    <p>Common decades: {common_decades} | Trend decades: {trend_decades}</p>
    <p class="note">RCP and SSP scenarios are analogous, not equivalent. Differences are expected.
    Old RCP and new SSP use different impact model ensembles.</p>

    <h2>Table of Contents</h2>
    <ol>
        <li><a href="#trends">Trend Analysis (Sen Slopes &amp; Spearman R²)</a></li>
        <li><a href="#agreement">Trend Sign Agreement Maps</a></li>
        <li><a href="#trend-scatter">Trend Scatter Plots</a></li>
        <li><a href="#values">Raw Value Comparison (Annual Mean, 2050s)</a></li>
        <li><a href="#scatter">Value Distribution Scatter</a></li>
        <li><a href="#seasonal">Seasonal Cycle Comparison</a></li>
    </ol>
    """)

    # =======================================================================
    # Section 1: TREND ANALYSIS (main focus)
    # =======================================================================
    html_parts.append('<div class="section" id="trends">')
    html_parts.append("<h2>1. Trend Analysis &mdash; Theil-Sen Slopes &amp; Spatial Correlation</h2>")
    html_parts.append(
        '<p>Per-pixel Theil-Sen slope of Annual Mean (vt12) across decades '
        f'{trend_decades[0]}s-{trend_decades[-1]}s. Units: change per decade. '
        'Spatial Spearman R² measures how well old and new trend <em>patterns</em> match.</p>'
    )

    trend_results = []
    old_trends = {}
    new_trends = {}

    for ssp, rcp, label in SCENARIO_PAIRS:
        log(f"  Computing trends: {label}...")

        old_slope, old_rho = compute_pixel_trends(old_ds, rcp, trend_decades, vt=12)
        new_slope, new_rho = compute_pixel_trends(new_ds, ssp, trend_decades, vt=12)

        old_trends[(ssp, rcp)] = old_slope
        new_trends[(ssp, rcp)] = new_slope

        # Spatial correlation of trend maps
        corr = spatial_correlation(old_slope, new_slope)

        # Sign agreement
        mask = ~np.isnan(old_slope) & ~np.isnan(new_slope)
        agree = ((old_slope > 0) & (new_slope > 0)) | ((old_slope < 0) & (new_slope < 0))
        n_total = int(mask.sum())
        n_agree = int((mask & agree).sum())
        pct = n_agree / n_total * 100 if n_total > 0 else 0

        trend_results.append({
            "label": label,
            "pearson_r": corr["pearson_r"],
            "spearman_r": corr["spearman_r"],
            "sign_agree_pct": pct,
            "n_valid": corr["n_valid"],
        })

        # Three-panel trend map
        html_parts.append(f"<h3>{label}</h3>")
        fig = create_trend_map(old_slope, new_slope, lats, lons, label, corr)
        add_fig(fig)

    # Summary table
    html_parts.append("<h3>Trend Correlation Summary</h3>")
    html_parts.append(build_trend_stats_table(trend_results))
    html_parts.append("</div>")

    # =======================================================================
    # Section 2: TREND SIGN AGREEMENT
    # =======================================================================
    html_parts.append('<div class="section" id="agreement">')
    html_parts.append("<h2>2. Trend Sign Agreement Maps</h2>")
    html_parts.append(
        "<p>Where do old and new agree on the <em>direction</em> of change? "
        "Blue = both wetting, Red = both drying, Grey = disagree.</p>"
    )

    for ssp, rcp, label in SCENARIO_PAIRS:
        old_slope = old_trends[(ssp, rcp)]
        new_slope = new_trends[(ssp, rcp)]
        fig = create_trend_agreement_map(old_slope, new_slope, lats, lons, label)
        add_fig(fig)

    html_parts.append("</div>")

    # =======================================================================
    # Section 3: TREND SCATTER PLOTS
    # =======================================================================
    html_parts.append('<div class="section" id="trend-scatter">')
    html_parts.append("<h2>3. Trend Scatter Plots</h2>")
    html_parts.append(
        "<p>Old per-pixel trend vs New per-pixel trend. Points near the 1:1 line "
        "indicate matching trend magnitude; OLS fit shows systematic bias.</p>"
    )

    for ssp, rcp, label in SCENARIO_PAIRS:
        old_slope = old_trends[(ssp, rcp)]
        new_slope = new_trends[(ssp, rcp)]
        fig = create_trend_scatter(old_slope, new_slope, label)
        add_fig(fig)

    html_parts.append("</div>")

    # =======================================================================
    # Section 4: RAW VALUE COMPARISON (reduced — just Annual_Mean at 2050)
    # =======================================================================
    focus_decade = 2050 if 2050 in common_decades else common_decades[-1]

    html_parts.append('<div class="section" id="values">')
    html_parts.append(f"<h2>4. Raw Value Comparison &mdash; {focus_decade}s</h2>")
    html_parts.append(
        f"<p>Annual Mean (vt12) at {focus_decade}s. Note: scales differ due to normalization "
        "of new SSP data. Spatial patterns should be broadly similar.</p>"
    )

    html_parts.append(build_value_stats_table(old_ds, new_ds, focus_decade))

    for ssp, rcp, label in SCENARIO_PAIRS:
        vt = 12
        vt_name = VALUE_TYPE_NAMES[vt]

        old_2d = load_slice(old_ds, rcp, vt, focus_decade)
        new_2d = load_slice(new_ds, ssp, vt, focus_decade)

        corr = spatial_correlation(old_2d, new_2d)
        cov = coverage_comparison(old_2d, new_2d)
        rmsd_val = rmsd(old_2d, new_2d)
        stats_dict = {**corr, **cov, "rmsd": rmsd_val}

        title = f"{label} | {focus_decade}s | {vt_name}"
        fig = create_three_panel_map(old_2d, new_2d, lats, lons, title, stats_dict)
        add_fig(fig)

    html_parts.append("</div>")

    # =======================================================================
    # Section 5: SCATTER & DISTRIBUTION (fixed with clipping)
    # =======================================================================
    html_parts.append('<div class="section" id="scatter">')
    html_parts.append(f"<h2>5. Value Distribution Scatter &mdash; {focus_decade}s</h2>")
    html_parts.append(
        "<p>Per-pixel old vs new values (clipped to 2nd-98th percentile). "
        "Scale differences are expected due to normalization.</p>"
    )

    for ssp, rcp, label in SCENARIO_PAIRS:
        for vt in [12, 16]:  # Annual_Mean, Q50
            vt_name = VALUE_TYPE_NAMES[vt]
            old_2d = load_slice(old_ds, rcp, vt, focus_decade)
            new_2d = load_slice(new_ds, ssp, vt, focus_decade)
            fig = create_scatter_histogram(old_2d, new_2d, f"{label} | {focus_decade}s | {vt_name}")
            add_fig(fig)

    html_parts.append("</div>")

    # =======================================================================
    # Section 6: SEASONAL CYCLES
    # =======================================================================
    html_parts.append('<div class="section" id="seasonal">')
    html_parts.append("<h2>6. Seasonal Cycle Comparisons</h2>")

    for ssp, rcp, label in SCENARIO_PAIRS:
        decade = focus_decade
        log(f"  Seasonal cycles: {label} / {decade}s...")
        fig = create_seasonal_comparison(old_ds, new_ds, rcp, ssp, decade)
        html_parts.append(f"<h3>{label} &mdash; {decade}s</h3>")
        add_fig(fig)

    html_parts.append("</div>")

    # --- Close ---
    html_parts.append("</body></html>")
    old_ds.close()
    new_ds.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(html_parts))

    size_mb = output_path.stat().st_size / (1024 * 1024)
    log(f"\nReport saved: {output_path} ({size_mb:.1f} MB)")


def main():
    global _VARIABLE

    parser = argparse.ArgumentParser(
        description="Generate interactive HTML comparison of old (RCP) vs new (SSP) water index files"
    )
    parser.add_argument(
        "--variable", type=str, default="tws",
        help="Variable to compare (default: tws). Options: tws, rootmoist, dis, qr, potevap, precip",
    )
    parser.add_argument(
        "--old", type=Path, default=None,
        help="Old (RCP) water index file (auto-derived from variable if not specified)",
    )
    parser.add_argument(
        "--new", type=Path, default=None,
        help="New (SSP) water index file (auto-derived from variable if not specified)",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output HTML file (auto-derived from variable if not specified)",
    )
    args = parser.parse_args()

    # Derive paths from variable config
    var_config = get_variable_config(args.variable)
    _VARIABLE = var_config.name

    old_path = args.old or Path(rf"C:\Cai_data\WaterIndex\waterIndexUnderlyingData_{var_config.name}.nc")
    new_path = args.new or Path(rf"C:\Cai_data\WaterIndex\{var_config.output_filename}")
    output_path = args.output or (PROJECT_ROOT / "reports" / f"water_{var_config.name}_comparison.html")

    for fpath, label in [(old_path, "Old"), (new_path, "New")]:
        if not fpath.exists():
            log(f"ERROR: {label} file not found: {fpath}")
            sys.exit(1)

    generate_report(old_path, new_path, output_path)


if __name__ == "__main__":
    main()
