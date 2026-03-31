"""Diagnostic analysis of raw streamflow/discharge (dis) data per impact model.

Computes per-model distribution statistics (mean, SD, IQR, min, max, percentiles)
across all GCMs/scenarios, and generates boxplots, histograms, seasonal cycle
comparison, and summary tables for visual comparison of model baselines.

Discharge values span many orders of magnitude (0 in deserts to 100,000+ m3/s
for major rivers), so log-scale visualizations are essential.

Usage:
    python scripts/diagnose_dis_models.py
    python scripts/diagnose_dis_models.py --raw-dir data/raw/water_dis --output reports/dis_model_diagnostics.html
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

# ── Configuration ──────────────────────────────────────────────────────────

MODELS = ["cwatm", "h08", "jules-w2", "miroc-integ-land", "watergap2-2e"]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

# Sample years for the diagnostic (use first decade 2015-2024 to capture baseline behavior)
SAMPLE_YEARS = list(range(2015, 2025))

# Subsample stride for spatial dimension (every Nth cell) to keep memory manageable
SPATIAL_STRIDE = 4  # sample ~1/16 of cells

MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

MODEL_COLORS = {
    "cwatm": "#1f77b4",
    "h08": "#ff7f0e",
    "jules-w2": "#2ca02c",
    "miroc-integ-land": "#d62728",
    "watergap2-2e": "#9467bd",
}


def find_raw_files(raw_dir: Path, model: str, gcm: str, scenario: str) -> list[Path]:
    """Find the raw NetCDF file(s) for a given model/gcm/scenario.

    Dis files may be split into two time segments (e.g., 2015-2060 and 2061-2100).
    Returns all matching files sorted by name.
    """
    subdir = raw_dir / model / f"{gcm}_{scenario}"
    if not subdir.exists():
        return []
    nc_files = sorted(subdir.glob("*.nc"))
    # Only keep monthly files (exclude daily files that may have been downloaded)
    nc_files = [f for f in nc_files if "daily" not in f.name]
    # Filter to prefer 2015soc over other social forcing variants
    soc_2015 = [f for f in nc_files if "_2015soc_" in f.name]
    soc_from_hist = [f for f in nc_files if "2015soc-from-histsoc" in f.name]
    if soc_2015:
        return soc_2015
    if soc_from_hist:
        return soc_from_hist
    return nc_files


def load_annual_means(filepaths: list[Path], years: list) -> np.ndarray | None:
    """Load raw monthly dis, compute annual means for specified years.

    Returns 2D array of shape (n_years, n_sampled_cells) with annual mean dis values in m3/s.
    """
    try:
        datasets = []
        for fpath in filepaths:
            datasets.append(xr.open_dataset(fpath))

        if len(datasets) == 1:
            ds = datasets[0]
        else:
            ds = xr.concat(datasets, dim="time")
            ds = ds.sortby("time")

        var_name = "dis" if "dis" in ds.data_vars else list(ds.data_vars)[0]
        da = ds[var_name]

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

        for d in datasets:
            d.close()

        if not annual_values:
            return None
        return np.stack(annual_values, axis=0)  # (n_years, n_cells)

    except Exception as e:
        print(f"  ERROR loading {filepaths[0].name}: {e}", flush=True)
        return None


def load_monthly_means(filepaths: list[Path], years: list) -> np.ndarray | None:
    """Load raw monthly dis, return per-month means.

    Returns array of shape (12,) with mean dis for each calendar month across
    the sampled years, spatially averaged over land cells.
    """
    try:
        datasets = []
        for fpath in filepaths:
            datasets.append(xr.open_dataset(fpath))

        if len(datasets) == 1:
            ds = datasets[0]
        else:
            ds = xr.concat(datasets, dim="time")
            ds = ds.sortby("time")

        var_name = "dis" if "dis" in ds.data_vars else list(ds.data_vars)[0]
        da = ds[var_name]

        times = da.time.values
        try:
            time_years = np.array([t.year for t in times])
            time_months = np.array([t.month for t in times])
        except AttributeError:
            import pandas as pd
            dti = pd.DatetimeIndex(times)
            time_years = dti.year.values
            time_months = dti.month.values

        monthly_spatial_means = [[] for _ in range(12)]
        for yr in years:
            for mon in range(1, 13):
                mask = (time_years == yr) & (time_months == mon)
                if mask.sum() == 0:
                    continue
                idx = np.where(mask)[0]
                data = da.isel(time=idx).values
                data_sub = data[:, ::SPATIAL_STRIDE, ::SPATIAL_STRIDE]
                spatial_mean = np.nanmean(data_sub)
                monthly_spatial_means[mon - 1].append(spatial_mean)

        for d in datasets:
            d.close()

        result = np.full(12, np.nan)
        for m in range(12):
            if monthly_spatial_means[m]:
                result[m] = np.mean(monthly_spatial_means[m])
        return result

    except Exception as e:
        print(f"  ERROR loading monthly means: {e}", flush=True)
        return None


def compute_stats(values: np.ndarray) -> dict:
    """Compute summary statistics from a flat array of valid values."""
    v = values[~np.isnan(values)]
    if len(v) == 0:
        return {"n": 0}
    n_negative = int(np.sum(v < 0))
    n_zero = int(np.sum(v == 0))
    # Log-scale stats for positive values (discharge spans many orders of magnitude)
    v_pos = v[v > 0]
    log_mean = float(np.mean(np.log10(v_pos))) if len(v_pos) > 0 else float("nan")
    return {
        "n": int(len(v)),
        "n_negative": n_negative,
        "pct_negative": float(n_negative / len(v) * 100),
        "n_zero": n_zero,
        "pct_zero": float(n_zero / len(v) * 100),
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
        "log10_mean": log_mean,
    }


def main():
    parser = argparse.ArgumentParser(description="Diagnostic analysis of raw dis model data")
    parser.add_argument("--raw-dir", type=Path,
                        default=Path("data/raw/water_dis"),
                        help="Path to raw dis data directory")
    parser.add_argument("--output", type=Path,
                        default=Path("reports/dis_model_diagnostics.html"),
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
    print(f"Sampling years: {SAMPLE_YEARS[0]}-{SAMPLE_YEARS[-1]}", flush=True)
    print(f"Spatial stride: every {SPATIAL_STRIDE}th cell", flush=True)
    print(f"Units: m3 s-1 (native)\n", flush=True)

    # ── Collect data per model ─────────────────────────────────────────────
    model_data = {}
    model_gcm_data = {}
    model_stats = {}
    model_gcm_stats = {}
    model_seasonal = {}

    for model in MODELS:
        print(f"Processing {model}...", flush=True)
        all_values = []
        gcm_values = defaultdict(list)
        seasonal_cycles = []

        for gcm in GCMS:
            for scenario in SCENARIOS:
                fpaths = find_raw_files(raw_dir, model, gcm, scenario)
                if not fpaths:
                    continue

                ann = load_annual_means(fpaths, SAMPLE_YEARS)
                if ann is None:
                    continue

                flat = ann.ravel()
                valid = flat[~np.isnan(flat)]
                all_values.append(valid)
                gcm_values[gcm].append(valid)

                seasonal = load_monthly_means(fpaths, SAMPLE_YEARS)
                if seasonal is not None:
                    seasonal_cycles.append(seasonal)

        if all_values:
            pooled = np.concatenate(all_values)
            model_data[model] = pooled
            model_stats[model] = compute_stats(pooled)
            s = model_stats[model]
            print(f"  {model}: {len(pooled):,} values, "
                  f"mean={s['mean']:.1f} m3/s, "
                  f"median={s['p50']:.1f} m3/s, "
                  f"range=[{s['min']:.2f}, {s['max']:.0f}], "
                  f"negative={s['pct_negative']:.1f}%, "
                  f"zero={s['pct_zero']:.1f}%",
                  flush=True)
        else:
            print(f"  {model}: NO DATA", flush=True)

        if seasonal_cycles:
            model_seasonal[model] = np.nanmean(seasonal_cycles, axis=0)

        model_gcm_data[model] = {}
        model_gcm_stats[model] = {}
        for gcm in GCMS:
            if gcm in gcm_values and gcm_values[gcm]:
                pooled_gcm = np.concatenate(gcm_values[gcm])
                model_gcm_data[model][gcm] = pooled_gcm
                model_gcm_stats[model][gcm] = compute_stats(pooled_gcm)

    # ── Print summary table ────────────────────────────────────────────────
    print("\n" + "=" * 120)
    print(f"{'Model':<20} {'Mean':>10} {'SD':>10} {'Min':>8} {'P25':>8} "
          f"{'Median':>8} {'P75':>10} {'Max':>12} {'IQR':>10} {'%Neg':>6} {'%Zero':>6} {'log10m':>7}")
    print("-" * 120)
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        print(f"{model:<20} {s['mean']:>10.1f} {s['std']:>10.1f} {s['min']:>8.2f} "
              f"{s['p25']:>8.2f} {s['p50']:>8.2f} {s['p75']:>10.1f} "
              f"{s['max']:>12.0f} {s['iqr']:>10.2f} {s['pct_negative']:>5.1f}% {s['pct_zero']:>5.1f}% {s['log10_mean']:>7.2f}")
    print("=" * 120)

    # Print normalization recommendation
    if len(model_stats) >= 2:
        means = [model_stats[m]["mean"] for m in MODELS if m in model_stats]
        medians = [model_stats[m]["p50"] for m in MODELS if m in model_stats]
        iqrs = [model_stats[m]["iqr"] for m in MODELS if m in model_stats]
        mean_ratio = max(means) / min(means) if min(means) > 0 else float("inf")
        median_ratio = max(medians) / min(medians) if min(medians) > 0 else float("inf")
        iqr_ratio = max(iqrs) / min(iqrs) if min(iqrs) > 0 else float("inf")
        # Also check log-scale alignment
        log_means = [model_stats[m]["log10_mean"] for m in MODELS if m in model_stats and not np.isnan(model_stats[m]["log10_mean"])]
        log_spread = max(log_means) - min(log_means) if log_means else float("inf")
        print(f"\nNormalization indicators:")
        print(f"  Mean ratio (max/min):   {mean_ratio:.2f}x")
        print(f"  Median ratio (max/min): {median_ratio:.2f}x")
        print(f"  IQR ratio (max/min):    {iqr_ratio:.2f}x")
        print(f"  Log10(mean) spread:     {log_spread:.2f} decades")
        if mean_ratio > 3 or median_ratio > 3:
            print("  >> RECOMMENDATION: Normalization likely needed (large scale divergence)")
        elif mean_ratio > 1.5 or median_ratio > 1.5:
            print("  >> RECOMMENDATION: Normalization may be beneficial (moderate divergence)")
        else:
            print("  >> RECOMMENDATION: Normalization probably unnecessary (models well-aligned)")

    # ── Generate HTML report with Plotly ───────────────────────────────────

    # Figure 1: Boxplots per model (log scale for discharge)
    fig1 = go.Figure()
    for model in MODELS:
        if model not in model_data:
            continue
        vals = model_data[model]
        # Use positive values only for log-scale box plot
        vals_pos = vals[vals > 0]
        if len(vals_pos) > 200_000:
            rng = np.random.default_rng(42)
            vals_pos = rng.choice(vals_pos, 200_000, replace=False)
        fig1.add_trace(go.Box(
            y=np.log10(vals_pos).tolist(),
            name=model,
            marker_color=MODEL_COLORS.get(model, "#888"),
            boxmean='sd',
        ))
    fig1.update_layout(
        title="Discharge Annual Mean Distribution by Model (log10 scale, 2015-2024)",
        yaxis_title="log10(Discharge, m³ s⁻¹)",
        xaxis_title="Impact Model",
        height=600,
        template="plotly_white",
    )

    # Figure 2: Boxplots per model x GCM (log scale)
    fig2 = go.Figure()
    for model in MODELS:
        if model not in model_gcm_data:
            continue
        for gcm in GCMS:
            if gcm not in model_gcm_data[model]:
                continue
            vals = model_gcm_data[model][gcm]
            vals_pos = vals[vals > 0]
            if len(vals_pos) > 50_000:
                rng = np.random.default_rng(42)
                vals_pos = rng.choice(vals_pos, 50_000, replace=False)
            label = f"{model}<br>{gcm}"
            fig2.add_trace(go.Box(
                y=np.log10(vals_pos).tolist(),
                name=label,
                marker_color=MODEL_COLORS.get(model, "#888"),
                boxmean='sd',
            ))
    fig2.update_layout(
        title="Discharge Annual Mean by Model x GCM (log10 scale, 2015-2024)",
        yaxis_title="log10(Discharge, m³ s⁻¹)",
        xaxis_title="Model x GCM",
        height=700,
        template="plotly_white",
        showlegend=False,
        xaxis_tickangle=-45,
    )

    # Figure 3: Overlapping histograms (log-x, log-y)
    fig3 = go.Figure()
    for model in MODELS:
        if model not in model_data:
            continue
        vals = model_data[model]
        vals_pos = vals[vals > 0]
        log_vals = np.log10(vals_pos)
        fig3.add_trace(go.Histogram(
            x=log_vals.tolist(),
            name=model,
            marker_color=MODEL_COLORS.get(model, "#888"),
            opacity=0.5,
            nbinsx=200,
        ))
    fig3.update_layout(
        title="Discharge Distribution Overlap (log10 scale, positive values only)",
        xaxis_title="log10(Discharge, m³ s⁻¹)",
        yaxis_title="Count",
        yaxis_type="log",
        barmode="overlay",
        height=500,
        template="plotly_white",
    )

    # Figure 4: Percentile heatmap
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
        colorbar_title="m³/s",
    ))
    fig4.update_layout(
        title="Percentile Comparison Across Models (m³ s⁻¹)",
        height=350,
        template="plotly_white",
    )

    # Figure 5: Seasonal cycle comparison
    fig5 = go.Figure()
    for model in MODELS:
        if model not in model_seasonal:
            continue
        cycle = model_seasonal[model]
        fig5.add_trace(go.Scatter(
            x=MONTH_NAMES,
            y=cycle.tolist(),
            mode="lines+markers",
            name=model,
            marker_color=MODEL_COLORS.get(model, "#888"),
            line=dict(width=2),
        ))
    fig5.update_layout(
        title="Mean Seasonal Cycle by Impact Model (2015-2024, global spatial mean, all GCMs/SSPs pooled)",
        xaxis_title="Month",
        yaxis_title="Discharge (m³ s⁻¹)",
        height=500,
        template="plotly_white",
    )

    # ── Build HTML ─────────────────────────────────────────────────────────
    stats_html = "<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse; font-family:monospace;'>\n"
    stats_html += ("<tr><th>Model</th><th>N</th><th>Mean</th><th>SD</th>"
                   "<th>Min</th><th>P01</th><th>P05</th><th>P25</th><th>Median</th>"
                   "<th>P75</th><th>P95</th><th>P99</th><th>Max</th><th>IQR</th>"
                   "<th>%Neg</th><th>%Zero</th><th>log10(mean)</th></tr>\n")
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        stats_html += (f"<tr><td><b>{model}</b></td><td>{s['n']:,}</td>"
                       f"<td>{s['mean']:.1f}</td><td>{s['std']:.1f}</td>"
                       f"<td>{s['min']:.2f}</td><td>{s['p01']:.2f}</td>"
                       f"<td>{s['p05']:.2f}</td><td>{s['p25']:.2f}</td>"
                       f"<td>{s['p50']:.2f}</td><td>{s['p75']:.1f}</td>"
                       f"<td>{s['p95']:.0f}</td><td>{s['p99']:.0f}</td>"
                       f"<td>{s['max']:.0f}</td><td>{s['iqr']:.2f}</td>"
                       f"<td>{s['pct_negative']:.1f}%</td><td>{s['pct_zero']:.1f}%</td>"
                       f"<td>{s['log10_mean']:.2f}</td></tr>\n")
    stats_html += "</table>\n"

    gcm_html = "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse; font-family:monospace; font-size:0.85em;'>\n"
    gcm_html += "<tr><th>Model</th><th>GCM</th><th>Mean</th><th>SD</th><th>Median</th><th>IQR</th><th>Min</th><th>Max</th><th>%Neg</th><th>%Zero</th></tr>\n"
    for model in MODELS:
        if model not in model_gcm_stats:
            continue
        for gcm in GCMS:
            if gcm not in model_gcm_stats[model]:
                continue
            s = model_gcm_stats[model][gcm]
            gcm_html += (f"<tr><td>{model}</td><td>{gcm}</td>"
                         f"<td>{s['mean']:.1f}</td><td>{s['std']:.1f}</td>"
                         f"<td>{s['p50']:.2f}</td><td>{s['iqr']:.2f}</td>"
                         f"<td>{s['min']:.2f}</td><td>{s['max']:.0f}</td>"
                         f"<td>{s['pct_negative']:.1f}%</td><td>{s['pct_zero']:.1f}%</td></tr>\n")
    gcm_html += "</table>\n"

    seasonal_html = ""
    if model_seasonal:
        seasonal_html = "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse; font-family:monospace;'>\n"
        seasonal_html += "<tr><th>Model</th>" + "".join(f"<th>{m}</th>" for m in MONTH_NAMES) + "<th>Annual</th></tr>\n"
        for model in MODELS:
            if model not in model_seasonal:
                continue
            cycle = model_seasonal[model]
            seasonal_html += f"<tr><td><b>{model}</b></td>"
            for v in cycle:
                seasonal_html += f"<td>{v:.1f}</td>"
            seasonal_html += f"<td><b>{np.nanmean(cycle):.1f}</b></td></tr>\n"
        seasonal_html += "</table>\n"

    # Normalization summary
    norm_html = ""
    if len(model_stats) >= 2:
        means = [model_stats[m]["mean"] for m in MODELS if m in model_stats]
        medians = [model_stats[m]["p50"] for m in MODELS if m in model_stats]
        iqrs = [model_stats[m]["iqr"] for m in MODELS if m in model_stats]
        mean_ratio = max(means) / min(means) if min(means) > 0 else float("inf")
        median_ratio = max(medians) / min(medians) if min(medians) > 0 else float("inf")
        iqr_ratio = max(iqrs) / min(iqrs) if min(iqrs) > 0 else float("inf")

        if mean_ratio > 3 or median_ratio > 3:
            box_class = "warn"
            recommendation = "Normalization likely needed — large scale divergence between models."
        elif mean_ratio > 1.5 or median_ratio > 1.5:
            box_class = "note"
            recommendation = "Normalization may be beneficial — moderate divergence between models."
        else:
            box_class = "good"
            recommendation = "Normalization probably unnecessary — models are well-aligned."

        norm_html = f"""
        <div class="{box_class}">
        <b>Normalization Assessment:</b><br>
        Mean ratio (max/min): {mean_ratio:.2f}x | Median ratio: {median_ratio:.2f}x | IQR ratio: {iqr_ratio:.2f}x<br>
        <b>{recommendation}</b>
        </div>
        """

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Streamflow/Discharge (dis) Model Diagnostic Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
    body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }}
    h1 {{ color: #333; }}
    h2 {{ color: #555; margin-top: 40px; }}
    .note {{ background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; margin: 16px 0; }}
    .warn {{ background: #f8d7da; padding: 12px; border-left: 4px solid #dc3545; margin: 16px 0; }}
    .good {{ background: #d4edda; padding: 12px; border-left: 4px solid #28a745; margin: 16px 0; }}
    table {{ margin: 16px 0; }}
    th {{ background: #f0f0f0; }}
</style>
</head><body>

<h1>Streamflow/Discharge (dis) Model Diagnostic Report</h1>
<p>Analysis of raw ISIMIP3b dis data across {len(MODELS)} impact models, {len(GCMS)} GCMs, {len(SCENARIOS)} SSP scenarios.</p>
<p>Sample period: {SAMPLE_YEARS[0]}-{SAMPLE_YEARS[-1]} | Spatial subsampling: every {SPATIAL_STRIDE}th cell | Units: m³ s⁻¹</p>

<div class="note">
<b>Purpose:</b> Identify systematic baseline offsets between impact models before ensemble aggregation.
Discharge spans many orders of magnitude (0 in deserts to 100,000+ m³/s for major rivers), so
log-scale comparisons are used alongside linear statistics. Note: some GCM gaps exist
(miroc-integ-land and watergap2-2e missing mri-esm2-0).
</div>

{norm_html}

<h2>1. Summary Statistics by Model</h2>
{stats_html}

<h2>2. Boxplot Comparison: Models (log10 scale)</h2>
<div id="fig1"></div>

<h2>3. Overlapping Histograms (log10 scale)</h2>
<div id="fig3"></div>

<h2>4. Percentile Heatmap</h2>
<div id="fig4"></div>

<h2>5. Seasonal Cycle Comparison</h2>
<p>Global spatial mean discharge by calendar month, averaged across all GCMs, SSPs, and sample years.</p>
<div id="fig5"></div>
{seasonal_html}

<h2>6. Boxplot: Model x GCM Breakdown (log10 scale)</h2>
<div id="fig2"></div>

<h2>7. Detailed Statistics: Model x GCM</h2>
{gcm_html}

</body>
<script>
Plotly.newPlot('fig1', {json.dumps(fig1.to_dict()['data'], default=str)}, {json.dumps(fig1.to_dict()['layout'], default=str)});
Plotly.newPlot('fig2', {json.dumps(fig2.to_dict()['data'], default=str)}, {json.dumps(fig2.to_dict()['layout'], default=str)});
Plotly.newPlot('fig3', {json.dumps(fig3.to_dict()['data'], default=str)}, {json.dumps(fig3.to_dict()['layout'], default=str)});
Plotly.newPlot('fig4', {json.dumps(fig4.to_dict()['data'], default=str)}, {json.dumps(fig4.to_dict()['layout'], default=str)});
Plotly.newPlot('fig5', {json.dumps(fig5.to_dict()['data'], default=str)}, {json.dumps(fig5.to_dict()['layout'], default=str)});
</script>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    print(f"\nHTML report saved: {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")

    if args.json_stats:
        json_path = args.json_stats
        if not json_path.is_absolute():
            json_path = Path(__file__).parent.parent / json_path
        combined = {
            "model_stats": model_stats,
            "model_gcm_stats": model_gcm_stats,
            "model_seasonal": {m: v.tolist() for m, v in model_seasonal.items()},
        }
        json_path.write_text(json.dumps(combined, indent=2))
        print(f"JSON stats saved: {json_path}")


if __name__ == "__main__":
    main()
