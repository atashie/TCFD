"""Diagnostic analysis of raw groundwater recharge (qr) data per impact model.

Computes per-model distribution statistics (mean, SD, IQR, min, max, percentiles)
across all GCMs/scenarios, and generates boxplots + summary tables for visual
comparison of model baselines. Checks for negative recharge values and spatial
coverage differences.

Usage:
    python scripts/diagnose_qr_models.py
    python scripts/diagnose_qr_models.py --raw-dir data/raw/water_qr --output reports/qr_model_diagnostics.html
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

MODELS = ["cwatm", "h08", "miroc-integ-land", "watergap2-2e"]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

# Sample years for the diagnostic (use first decade 2015-2024 to capture baseline behavior)
SAMPLE_YEARS = list(range(2015, 2025))

# Subsample stride for spatial dimension (every Nth cell) to keep memory manageable
SPATIAL_STRIDE = 4  # sample ~1/16 of cells

# Unit conversion: kg/m2/s -> mm/month
KG_M2_S_TO_MM_MONTH = 60 * 60 * 24 * 30.4375  # 2,629,800

MODEL_COLORS = {
    "cwatm": "#1f77b4",
    "h08": "#ff7f0e",
    "miroc-integ-land": "#d62728",
    "watergap2-2e": "#9467bd",
}


def find_raw_file(raw_dir: Path, model: str, gcm: str, scenario: str) -> Path | None:
    """Find the raw NetCDF file for a given model/gcm/scenario."""
    subdir = raw_dir / model / f"{gcm}_{scenario}"
    if not subdir.exists():
        return None
    nc_files = list(subdir.glob("*.nc"))
    return nc_files[0] if nc_files else None


def load_annual_means(filepath: Path, years: list) -> np.ndarray | None:
    """Load raw monthly qr, convert to mm/month, compute annual means for specified years.

    Returns 2D array of shape (n_years, n_sampled_cells) with annual mean qr values in mm/month.
    """
    try:
        ds = xr.open_dataset(filepath)
        var_name = "qr" if "qr" in ds.data_vars else list(ds.data_vars)[0]
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
            # Convert kg/m2/s -> mm/month
            monthly_mm = monthly * KG_M2_S_TO_MM_MONTH
            # Subsample spatially
            monthly_sub = monthly_mm[:, ::SPATIAL_STRIDE, ::SPATIAL_STRIDE]
            # Annual mean of monthly recharge rates
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
    n_negative = int(np.sum(v < 0))
    n_zero = int(np.sum(v == 0))
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
    }


def main():
    parser = argparse.ArgumentParser(description="Diagnostic analysis of raw qr model data")
    parser.add_argument("--raw-dir", type=Path,
                        default=Path("data/raw/water_qr"),
                        help="Path to raw qr data directory")
    parser.add_argument("--output", type=Path,
                        default=Path("reports/qr_model_diagnostics.html"),
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
    print(f"Unit conversion: kg/m2/s x {KG_M2_S_TO_MM_MONTH:.0f} -> mm/month\n", flush=True)

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
            s = model_stats[model]
            print(f"  {model}: {len(pooled):,} values, "
                  f"mean={s['mean']:.2f} mm/mo, "
                  f"median={s['p50']:.2f} mm/mo, "
                  f"range=[{s['min']:.2f}, {s['max']:.2f}], "
                  f"negative={s['pct_negative']:.1f}%, "
                  f"zero={s['pct_zero']:.1f}%",
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
    print("\n" + "=" * 110)
    print(f"{'Model':<20} {'Mean':>8} {'SD':>8} {'Min':>8} {'P25':>8} "
          f"{'Median':>8} {'P75':>8} {'Max':>10} {'IQR':>8} {'%Neg':>6} {'%Zero':>6}")
    print("-" * 110)
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        print(f"{model:<20} {s['mean']:>8.2f} {s['std']:>8.2f} {s['min']:>8.2f} "
              f"{s['p25']:>8.2f} {s['p50']:>8.2f} {s['p75']:>8.2f} "
              f"{s['max']:>10.1f} {s['iqr']:>8.2f} {s['pct_negative']:>5.1f}% {s['pct_zero']:>5.1f}%")
    print("=" * 110)

    # ── Generate HTML report with Plotly ───────────────────────────────────

    # Figure 1: Boxplots per model
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
        title="QR Annual Mean Distribution by Impact Model (2015-2024, all GCMs/SSPs pooled)",
        yaxis_title="Groundwater Recharge (mm month⁻¹)",
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
        title="QR Annual Mean Distribution by Impact Model x GCM (2015-2024)",
        yaxis_title="Groundwater Recharge (mm month⁻¹)",
        xaxis_title="Model x GCM",
        height=700,
        template="plotly_white",
        showlegend=False,
        xaxis_tickangle=-45,
    )

    # Figure 3: Overlapping histograms (log-y)
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
        title="QR Value Distribution Overlap (0.5th-99.5th percentile, all GCMs/SSPs pooled)",
        xaxis_title="Groundwater Recharge (mm month⁻¹)",
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
        text=[[f"{v:.2f}" for v in row] for row in z_data],
        texttemplate="%{text}",
        colorscale="Viridis",
        colorbar_title="mm/mo",
    ))
    fig4.update_layout(
        title="Percentile Comparison Across Models",
        height=300,
        template="plotly_white",
    )

    # ── Build HTML ─────────────────────────────────────────────────────────
    stats_html = "<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse; font-family:monospace;'>\n"
    stats_html += ("<tr><th>Model</th><th>N</th><th>Mean</th><th>SD</th>"
                   "<th>Min</th><th>P01</th><th>P05</th><th>P25</th><th>Median</th>"
                   "<th>P75</th><th>P95</th><th>P99</th><th>Max</th><th>IQR</th>"
                   "<th>%Neg</th><th>%Zero</th></tr>\n")
    for model in MODELS:
        if model not in model_stats:
            continue
        s = model_stats[model]
        stats_html += (f"<tr><td><b>{model}</b></td><td>{s['n']:,}</td>"
                       f"<td>{s['mean']:.2f}</td><td>{s['std']:.2f}</td>"
                       f"<td>{s['min']:.2f}</td><td>{s['p01']:.2f}</td>"
                       f"<td>{s['p05']:.2f}</td><td>{s['p25']:.2f}</td>"
                       f"<td>{s['p50']:.2f}</td><td>{s['p75']:.2f}</td>"
                       f"<td>{s['p95']:.2f}</td><td>{s['p99']:.2f}</td>"
                       f"<td>{s['max']:.1f}</td><td>{s['iqr']:.2f}</td>"
                       f"<td>{s['pct_negative']:.1f}%</td><td>{s['pct_zero']:.1f}%</td></tr>\n")
    stats_html += "</table>\n"

    gcm_html = "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse; font-family:monospace; font-size:0.85em;'>\n"
    gcm_html += "<tr><th>Model</th><th>GCM</th><th>Mean</th><th>SD</th><th>Median</th><th>IQR</th><th>Min</th><th>Max</th><th>%Neg</th></tr>\n"
    for model in MODELS:
        if model not in model_gcm_stats:
            continue
        for gcm in GCMS:
            if gcm not in model_gcm_stats[model]:
                continue
            s = model_gcm_stats[model][gcm]
            gcm_html += (f"<tr><td>{model}</td><td>{gcm}</td>"
                         f"<td>{s['mean']:.2f}</td><td>{s['std']:.2f}</td>"
                         f"<td>{s['p50']:.2f}</td><td>{s['iqr']:.2f}</td>"
                         f"<td>{s['min']:.2f}</td><td>{s['max']:.1f}</td>"
                         f"<td>{s['pct_negative']:.1f}%</td></tr>\n")
    gcm_html += "</table>\n"

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Groundwater Recharge (qr) Model Diagnostic Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
    body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }}
    h1 {{ color: #333; }}
    h2 {{ color: #555; margin-top: 40px; }}
    .note {{ background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; margin: 16px 0; }}
    .warn {{ background: #f8d7da; padding: 12px; border-left: 4px solid #dc3545; margin: 16px 0; }}
    table {{ margin: 16px 0; }}
    th {{ background: #f0f0f0; }}
</style>
</head><body>

<h1>Groundwater Recharge (qr) Model Diagnostic Report</h1>
<p>Analysis of raw ISIMIP3b qr data across {len(MODELS)} impact models, {len(GCMS)} GCMs, {len(SCENARIOS)} SSP scenarios.</p>
<p>Sample period: {SAMPLE_YEARS[0]}-{SAMPLE_YEARS[-1]} | Spatial subsampling: every {SPATIAL_STRIDE}th cell | Units: mm month<sup>-1</sup></p>
<p>Unit conversion: raw kg m<sup>-2</sup> s<sup>-1</sup> x {KG_M2_S_TO_MM_MONTH:,.0f} = mm month<sup>-1</sup></p>

<div class="note">
<b>Purpose:</b> Identify systematic baseline offsets between impact models before ensemble aggregation.
If models have very different absolute groundwater recharge ranges, normalization (robust z-score) may be
needed before computing cross-model ensemble statistics. Also check for negative recharge values
(net groundwater discharge to surface) which may be physically meaningful or model artifacts.
</div>

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
