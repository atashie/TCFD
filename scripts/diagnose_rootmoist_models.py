"""Diagnostic analysis of raw rootmoist data per impact model.

Computes per-model distribution statistics (mean, SD, IQR, min, max, percentiles)
across all GCMs/scenarios, and generates boxplots + summary tables for visual
comparison of model baselines.

Two-level comparison:
  1. Model-level (miroc-integ-land vs web-dhm-sg) — key diagnostic for normalization decision
  2. Per-model x GCM breakdown — within-model variation

Usage:
    python scripts/diagnose_rootmoist_models.py
    python scripts/diagnose_rootmoist_models.py --raw-dir data/raw/water_rootmoist --output reports/rootmoist_model_diagnostics.html
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import xarray as xr
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).parent))
from config_water_variables import get_variable_config, SHARED_CONFIG

# ── Configuration ──────────────────────────────────────────────────────────

cfg = get_variable_config("rootmoist")
MODELS = cfg.models
GCMS = list(SHARED_CONFIG.gcms)
SCENARIOS = list(SHARED_CONFIG.scenarios)
VARIABLE = cfg.name

# Sample years for the diagnostic (use first decade 2015-2024 to capture baseline behavior)
SAMPLE_YEARS = list(range(2015, 2025))

# Subsample stride for spatial dimension (every Nth cell) to keep memory manageable
SPATIAL_STRIDE = 4  # sample ~1/16 of cells

MODEL_COLORS = {
    "miroc-integ-land": "#d62728",
    "web-dhm-sg": "#1f77b4",
}


def find_raw_file(raw_dir: Path, model: str, gcm: str, scenario: str) -> Path | None:
    """Find the raw NetCDF file for a given model/gcm/scenario."""
    subdir = raw_dir / model / f"{gcm}_{scenario}"
    if not subdir.exists():
        return None
    nc_files = list(subdir.glob("*.nc"))
    return nc_files[0] if nc_files else None


def load_annual_means(filepath: Path, years: list) -> np.ndarray | None:
    """Load raw monthly rootmoist, compute annual means for specified years.

    Returns 2D array of shape (n_years, n_sampled_cells) with annual mean values.
    """
    try:
        ds = xr.open_dataset(filepath)
        var_name = VARIABLE if VARIABLE in ds.data_vars else list(ds.data_vars)[0]
        da = ds[var_name]

        # Select years and subsample spatially
        times = da.time.values
        try:
            time_years = np.array([t.year for t in times])
        except AttributeError:
            import pandas as pd
            time_years = pd.DatetimeIndex(times).year.values

        annual_values = []
        for yr in years:
            yr_mask = time_years == yr
            if yr_mask.sum() == 0:
                continue
            time_indices = np.where(yr_mask)[0]
            monthly = da.isel(time=time_indices).values  # (months, lat, lon)
            monthly_sub = monthly[:, ::SPATIAL_STRIDE, ::SPATIAL_STRIDE]
            ann_mean = np.nanmean(monthly_sub, axis=0)  # (lat_sub, lon_sub)
            annual_values.append(ann_mean.ravel())

        ds.close()

        if not annual_values:
            return None
        return np.stack(annual_values, axis=0)  # (n_years, n_cells)

    except Exception as e:
        print(f"  ERROR loading {filepath.name}: {e}", flush=True)
        return None


def compute_stats(values: np.ndarray) -> dict:
    """Compute summary statistics from a flat array of valid values."""
    v = values[~np.isnan(values)]
    if len(v) == 0:
        return {"n": 0}
    return {
        "n": int(len(v)),
        "mean": float(np.mean(v)),
        "std": float(np.std(v)),
        "min": float(np.min(v)),
        "p01": float(np.percentile(v, 1)),
        "p05": float(np.percentile(v, 5)),
        "p25": float(np.percentile(v, 25)),
        "p50": float(np.percentile(v, 50)),
        "p75": float(np.percentile(v, 75)),
        "p95": float(np.percentile(v, 95)),
        "p99": float(np.percentile(v, 99)),
        "max": float(np.max(v)),
        "iqr": float(np.percentile(v, 75) - np.percentile(v, 25)),
    }


def main():
    parser = argparse.ArgumentParser(description="Diagnostic analysis of raw rootmoist model data")
    parser.add_argument("--raw-dir", type=Path,
                        default=Path("data/raw/water_rootmoist"),
                        help="Path to raw rootmoist data directory")
    parser.add_argument("--output", type=Path,
                        default=Path("reports/rootmoist_model_diagnostics.html"),
                        help="Output HTML report path")
    parser.add_argument("--json-stats", type=Path, default=None,
                        help="Also save stats as JSON")
    args = parser.parse_args()

    raw_dir = args.raw_dir
    if not raw_dir.is_absolute():
        raw_dir = Path(__file__).parent.parent / raw_dir

    output_path = args.output
    if not output_path.is_absolute():
        output_path = Path(__file__).parent.parent / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Raw data dir: {raw_dir}", flush=True)
    print(f"Models: {MODELS}", flush=True)
    print(f"Sampling years: {SAMPLE_YEARS[0]}-{SAMPLE_YEARS[-1]}", flush=True)
    print(f"Spatial stride: every {SPATIAL_STRIDE}th cell\n", flush=True)

    # ── Collect data per model ─────────────────────────────────────────────
    model_data = {}           # model -> flat array of all annual-mean values
    model_gcm_data = {}       # model -> {gcm -> flat array}
    model_stats = {}          # model -> stats dict
    model_gcm_stats = {}      # model -> {gcm -> stats dict}

    for model in MODELS:
        print(f"Processing {model}...", flush=True)
        all_values = []
        gcm_values = defaultdict(list)

        for gcm in GCMS:
            for scenario in SCENARIOS:
                fpath = find_raw_file(raw_dir, model, gcm, scenario)
                if fpath is None:
                    print(f"  MISSING: {model}/{gcm}_{scenario}", flush=True)
                    continue

                ann = load_annual_means(fpath, SAMPLE_YEARS)
                if ann is None:
                    continue

                flat = ann.ravel()
                valid = flat[~np.isnan(flat)]
                all_values.append(valid)
                gcm_values[gcm].append(valid)

        if all_values:
            pooled = np.concatenate(all_values)
            model_data[model] = pooled
            model_stats[model] = compute_stats(pooled)
            print(f"  {model}: {len(pooled):,} values, "
                  f"mean={model_stats[model]['mean']:.1f}, "
                  f"std={model_stats[model]['std']:.1f}, "
                  f"range=[{model_stats[model]['min']:.1f}, {model_stats[model]['max']:.1f}]",
                  flush=True)
        else:
            print(f"  {model}: NO DATA", flush=True)

        model_gcm_data[model] = {}
        model_gcm_stats[model] = {}
        for gcm in GCMS:
            if gcm in gcm_values and gcm_values[gcm]:
                pooled_gcm = np.concatenate(gcm_values[gcm])
                model_gcm_data[model][gcm] = pooled_gcm
                model_gcm_stats[model][gcm] = compute_stats(pooled_gcm)

    # ── Print summary table ────────────────────────────────────────────────
    print("\n" + "=" * 90)
    print(f"{'Model':<20} {'Mean':>10} {'SD':>10} {'Min':>10} {'P25':>10} "
          f"{'Median':>10} {'P75':>10} {'Max':>10} {'IQR':>10}")
    print("-" * 90)
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        print(f"{model:<20} {s['mean']:>10.1f} {s['std']:>10.1f} {s['min']:>10.1f} "
              f"{s['p25']:>10.1f} {s['p50']:>10.1f} {s['p75']:>10.1f} "
              f"{s['max']:>10.1f} {s['iqr']:>10.1f}")
    print("=" * 90)

    # ── Normalization recommendation ──────────────────────────────────────
    if len(model_stats) >= 2:
        medians = [model_stats[m]["p50"] for m in model_stats]
        iqrs = [model_stats[m]["iqr"] for m in model_stats]
        median_ratio = max(medians) / min(medians) if min(medians) > 0 else float("inf")
        iqr_ratio = max(iqrs) / min(iqrs) if min(iqrs) > 0 else float("inf")

        print(f"\n--- Normalization Decision ---")
        print(f"Median ratio (max/min): {median_ratio:.2f}")
        print(f"IQR ratio (max/min): {iqr_ratio:.2f}")
        if median_ratio > 2 or iqr_ratio > 2:
            print("RECOMMENDATION: Normalization STRONGLY recommended (ratio > 2x)")
        elif median_ratio > 1.5 or iqr_ratio > 1.5:
            print("RECOMMENDATION: Normalization recommended (ratio > 1.5x)")
        else:
            print("RECOMMENDATION: Normalization optional (models have similar baselines)")

    # ── Generate HTML report with Plotly ───────────────────────────────────

    # Figure 1: Boxplots per model (main comparison)
    fig1 = go.Figure()
    for model in MODELS:
        if model not in model_data:
            continue
        vals = model_data[model]
        if len(vals) > 200_000:
            rng = np.random.default_rng(42)
            vals = rng.choice(vals, 200_000, replace=False)
        fig1.add_trace(go.Box(
            y=vals.tolist(),
            name=model,
            marker_color=MODEL_COLORS.get(model, "#888"),
            boxmean='sd',
        ))
    fig1.update_layout(
        title="Rootmoist Annual Mean Distribution by Impact Model (2015-2024, all GCMs/SSPs pooled)",
        yaxis_title="Rootmoist (kg m⁻²)",
        xaxis_title="Impact Model",
        height=600,
        template="plotly_white",
    )

    # Figure 2: Boxplots per model x GCM
    fig2 = go.Figure()
    for model in MODELS:
        if model not in model_gcm_data:
            continue
        for gcm in GCMS:
            if gcm not in model_gcm_data[model]:
                continue
            vals = model_gcm_data[model][gcm]
            if len(vals) > 50_000:
                rng = np.random.default_rng(42)
                vals = rng.choice(vals, 50_000, replace=False)
            label = f"{model}<br>{gcm}"
            fig2.add_trace(go.Box(
                y=vals.tolist(),
                name=label,
                marker_color=MODEL_COLORS.get(model, "#888"),
                boxmean='sd',
            ))

    fig2.update_layout(
        title="Rootmoist Annual Mean Distribution by Impact Model x GCM (2015-2024)",
        yaxis_title="Rootmoist (kg m⁻²)",
        xaxis_title="Model x GCM",
        height=700,
        template="plotly_white",
        showlegend=False,
        xaxis_tickangle=-45,
    )

    # Figure 3: Overlapping histograms (log-y for tail comparison)
    fig3 = go.Figure()
    for model in MODELS:
        if model not in model_data:
            continue
        vals = model_data[model]
        lo, hi = np.percentile(vals, [0.5, 99.5])
        clipped = vals[(vals >= lo) & (vals <= hi)]
        fig3.add_trace(go.Histogram(
            x=clipped.tolist(),
            name=model,
            marker_color=MODEL_COLORS.get(model, "#888"),
            opacity=0.5,
            nbinsx=200,
        ))
    fig3.update_layout(
        title="Rootmoist Value Distribution Overlap (0.5th-99.5th percentile, all GCMs/SSPs pooled)",
        xaxis_title="Rootmoist (kg m⁻²)",
        yaxis_title="Count",
        yaxis_type="log",
        barmode="overlay",
        height=500,
        template="plotly_white",
    )

    # Figure 4: Percentile comparison table as heatmap
    pct_labels = ["P01", "P05", "P25", "Median", "P75", "P95", "P99"]
    pct_keys = ["p01", "p05", "p25", "p50", "p75", "p95", "p99"]
    z_data = []
    model_labels = []
    for model in MODELS:
        if model not in model_stats:
            continue
        model_labels.append(model)
        z_data.append([model_stats[model][k] for k in pct_keys])

    fig4 = go.Figure(data=go.Heatmap(
        z=z_data,
        x=pct_labels,
        y=model_labels,
        text=[[f"{v:.1f}" for v in row] for row in z_data],
        texttemplate="%{text}",
        colorscale="Viridis",
        colorbar_title="kg m⁻²",
    ))
    fig4.update_layout(
        title="Percentile Comparison Across Models",
        height=300,
        template="plotly_white",
    )

    # ── Build HTML ─────────────────────────────────────────────────────────
    stats_html = "<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse; font-family:monospace;'>\n"
    stats_html += ("<tr><th>Model</th><th>N cells</th><th>Mean</th><th>SD</th>"
                   "<th>Min</th><th>P01</th><th>P05</th><th>P25</th><th>Median</th>"
                   "<th>P75</th><th>P95</th><th>P99</th><th>Max</th><th>IQR</th></tr>\n")
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        stats_html += (f"<tr><td><b>{model}</b></td><td>{s['n']:,}</td>"
                       f"<td>{s['mean']:.1f}</td><td>{s['std']:.1f}</td>"
                       f"<td>{s['min']:.1f}</td><td>{s['p01']:.1f}</td>"
                       f"<td>{s['p05']:.1f}</td><td>{s['p25']:.1f}</td>"
                       f"<td>{s['p50']:.1f}</td><td>{s['p75']:.1f}</td>"
                       f"<td>{s['p95']:.1f}</td><td>{s['p99']:.1f}</td>"
                       f"<td>{s['max']:.1f}</td><td>{s['iqr']:.1f}</td></tr>\n")
    stats_html += "</table>\n"

    # Per-model x GCM table
    gcm_html = "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse; font-family:monospace; font-size:0.85em;'>\n"
    gcm_html += "<tr><th>Model</th><th>GCM</th><th>Mean</th><th>SD</th><th>Median</th><th>IQR</th><th>Min</th><th>Max</th></tr>\n"
    for model in MODELS:
        if model not in model_gcm_stats:
            continue
        for gcm in GCMS:
            if gcm not in model_gcm_stats[model]:
                continue
            s = model_gcm_stats[model][gcm]
            gcm_html += (f"<tr><td>{model}</td><td>{gcm}</td>"
                         f"<td>{s['mean']:.1f}</td><td>{s['std']:.1f}</td>"
                         f"<td>{s['p50']:.1f}</td><td>{s['iqr']:.1f}</td>"
                         f"<td>{s['min']:.1f}</td><td>{s['max']:.1f}</td></tr>\n")
    gcm_html += "</table>\n"

    # Normalization recommendation HTML
    norm_html = ""
    if len(model_stats) >= 2:
        medians = [model_stats[m]["p50"] for m in model_stats]
        iqrs = [model_stats[m]["iqr"] for m in model_stats]
        median_ratio = max(medians) / min(medians) if min(medians) > 0 else float("inf")
        iqr_ratio = max(iqrs) / min(iqrs) if min(iqrs) > 0 else float("inf")
        if median_ratio > 2 or iqr_ratio > 2:
            rec = "STRONGLY recommended"
            color = "#dc3545"
        elif median_ratio > 1.5 or iqr_ratio > 1.5:
            rec = "Recommended"
            color = "#ffc107"
        else:
            rec = "Optional (models have similar baselines)"
            color = "#28a745"

        norm_html = f"""
        <div class="note" style="border-left-color: {color};">
        <b>Normalization Decision:</b><br>
        Median ratio (max/min): {median_ratio:.2f}<br>
        IQR ratio (max/min): {iqr_ratio:.2f}<br>
        <b style="color: {color};">Recommendation: Normalization {rec}</b>
        </div>
        """

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Rootmoist Model Diagnostic Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
    body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }}
    h1 {{ color: #333; }}
    h2 {{ color: #555; margin-top: 40px; }}
    .note {{ background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; margin: 16px 0; }}
    table {{ margin: 16px 0; }}
    th {{ background: #f0f0f0; }}
</style>
</head><body>

<h1>Rootmoist Model Diagnostic Report</h1>
<p>Analysis of raw ISIMIP3b rootmoist data across {len(MODELS)} impact models ({', '.join(MODELS)}), 5 GCMs, 3 SSP scenarios.</p>
<p>Sample period: {SAMPLE_YEARS[0]}-{SAMPLE_YEARS[-1]} | Spatial subsampling: every {SPATIAL_STRIDE}th cell | Units: kg m&minus;2</p>

<div class="note">
<b>Purpose:</b> Identify systematic baseline offsets between impact models before ensemble aggregation.
If models have very different absolute rootmoist ranges, normalization (robust z-score) may be
needed before computing cross-model ensemble statistics.
</div>

{norm_html}

<h2>1. Summary Statistics by Model</h2>
{stats_html}

<h2>2. Boxplot Comparison: Models</h2>
<div id="fig1"></div>

<h2>3. Overlapping Histograms (log scale)</h2>
<div id="fig3"></div>

<h2>4. Percentile Heatmap</h2>
<div id="fig4"></div>

<h2>5. Boxplot: Model x GCM Breakdown</h2>
<div id="fig2"></div>

<h2>6. Detailed Statistics: Model x GCM</h2>
{gcm_html}

</body>
<script>
Plotly.newPlot('fig1', {json.dumps(fig1.to_dict()['data'], default=str)}, {json.dumps(fig1.to_dict()['layout'], default=str)});
Plotly.newPlot('fig2', {json.dumps(fig2.to_dict()['data'], default=str)}, {json.dumps(fig2.to_dict()['layout'], default=str)});
Plotly.newPlot('fig3', {json.dumps(fig3.to_dict()['data'], default=str)}, {json.dumps(fig3.to_dict()['layout'], default=str)});
Plotly.newPlot('fig4', {json.dumps(fig4.to_dict()['data'], default=str)}, {json.dumps(fig4.to_dict()['layout'], default=str)});
</script>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    print(f"\nHTML report saved: {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")

    if args.json_stats:
        json_path = args.json_stats
        if not json_path.is_absolute():
            json_path = Path(__file__).parent.parent / json_path
        combined = {"model_stats": model_stats, "model_gcm_stats": model_gcm_stats}
        json_path.write_text(json.dumps(combined, indent=2))
        print(f"JSON stats saved: {json_path}")


if __name__ == "__main__":
    main()
