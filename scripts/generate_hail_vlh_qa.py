#!/usr/bin/env python3
"""QA/QC page for the global very-large-hail observational layer.

WHY NOT generate_maps.py
    That tool is bound to the OUTPUT-SPEC decadal contract -- it selects on `ds.decade`,
    builds a Trend tab from the two slopes and a Members tab from `n_members`, none of
    which exist here. Same reasoning as `generate_tornado_qa.py`, which is the other
    renderer for `observational-historical-v1`. This is deliberately a third narrow tool
    rather than a branch through the one that renders 23 shipped layers.

WHY IT EXISTS
    `qa_reviewed_on` is null. The registry convention is that the date goes in only after
    a human has looked at the maps, and nobody can look without maps. This makes them. It
    does not perform the review.

WHAT TO LOOK FOR
    1. Does the global maximum sit north-west of Cordoba, Argentina, with the Uruguay /
       Paraguay / southern Brazil tri-border, the US Great Plains, north-east South Africa
       and southern Alberta/Saskatchewan behind it? Those are the paper's own hotspots. A
       maximum somewhere else means the ingest is wrong, whatever the contract check says.
    2. COASTLINE ARTEFACTS. The source excludes maritime cells more than 100 km offshore
       and this layer further masks to land. A rim of high values hugging every coast would
       be the mask boundary, not weather.
    3. Does the percentile panel look like a hazard field or like the land mask? The zero
       tier is 8.2% of covered cells; if the map reads as a binary land/no-land split then
       the tier boundary is doing the work the ranking should be doing.
    4. The relative-interval panel should be WIDE where the rate is low and narrow where it
       is high -- that is what a Poisson sampling interval does. If it is wide in the
       hotspots something is wrong with the posterior.
    5. Native-resolution zoom over South America: look for grid stripes or blocky steps.
       ERA5 is 0.25 deg natively so there should be none; visible lattice structure would
       mean the padding or the coordinate join went wrong.

DISPLAY RESOLUTION IS NOT THE FILE'S RESOLUTION
    The global panels are BLOCK-AVERAGED to 0.5 deg for the browser -- 720x1440 per panel
    is ~1M numbers serialised as text, which is several MB per panel. `block_mean` is used,
    never `values[::2, ::2]`: on a sparse hazard, striding deletes cells rather than
    averaging them. The zoom panel is at the file's native 0.25 deg precisely so the
    coarsening cannot hide a lattice artefact. The FILE is 0.25 deg throughout.

USAGE
    python3 scripts/generate_hail_vlh_qa.py
"""

from __future__ import annotations

import argparse
import html
import sys
import warnings
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.viz_common import FONT_STACK, SEQUENTIAL_BLUE, ORDINAL_RISK, plotly_colorscale  # noqa: E402

LAYER = "hail-vlh-archamo_ge5cm"          # the processed-layer folder
PROCESSED = Path("data/processed") / LAYER / "hail-vlh_observed_processed.nc"
# The maps folder uses the SHORT hazard name, matching reports/maps/{tornado,landslide,...}.
# It is deliberately not derived from LAYER: the processed folder carries the source and
# threshold in its name, the maps folder is the one a person browses.
OUT_DIR = Path("reports/maps/hail-vlh")

ZOOM = dict(name="South America (native 0.25 deg)", lat=(-45.0, -15.0), lon=(-75.0, -45.0))

REFERENCE_SITES = {
    "Cordoba, Argentina": (-31.4, -64.2),
    "Mendoza, Argentina": (-32.9, -68.8),
    "Johannesburg, South Africa": (-26.2, 28.0),
    "Oklahoma City, USA": (35.5, -97.5),
    "Calgary, Canada": (51.0, -114.1),
    "Milan, Italy": (45.5, 9.2),
    "Brisbane, Australia": (-27.5, 153.0),
    "Beijing, China": (39.9, 116.4),
}


def block_mean(a: np.ndarray, k: int = 2) -> np.ndarray:
    """Average k x k blocks, ignoring NaN. NOT slicing -- see the module docstring."""
    ny, nx = a.shape
    a = a[: ny // k * k, : nx // k * k].reshape(ny // k, k, nx // k, k)
    with np.errstate(invalid="ignore"):
        return np.nanmean(a, axis=(1, 3))


def coarse_coords(v: np.ndarray, k: int = 2) -> np.ndarray:
    n = v.size // k * k
    return v[:n].reshape(-1, k).mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--processed", type=Path, default=PROCESSED)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    # Empty blocks are the ocean; averaging nothing is the expected case, not a defect.
    warnings.filterwarnings("ignore", message="Mean of empty slice")

    if not args.processed.exists():
        print(f"layer not found: {args.processed}\nRun scripts/process_hail_vlh_battaglioli.py first.",
              file=sys.stderr)
        return 2
    ds = xr.open_dataset(args.processed)
    lat = ds["lat"].values
    lon = ds["lon"].values
    med = ds["median"].values
    pct = ds["percentile"].values
    lo = ds["lower_ci"].values
    hi = ds["upper_ci"].values
    cov = ds["coverage_mask"].values == 1
    n_ev = ds["n_events"].values

    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.where(med > 0, (hi - lo) / med, np.nan)

    clat = coarse_coords(lat)
    clon = coarse_coords(lon)
    # `n_events` is not panelled: it is the rate times 74 and adds a quarter of a million
    # numbers to the payload to show the same map twice.
    panels = [
        ("Expected events per box per year", block_mean(med), "events yr⁻¹", SEQUENTIAL_BLUE, None, 4),
        ("Percentile (two-tier, covered cells only)", block_mean(pct), "1–100", ORDINAL_RISK, None, 1),
        ("Relative interval width (upper−lower)/median", block_mean(rel), "ratio", SEQUENTIAL_BLUE, (0, 4), 2),
        ("Coverage mask (1 = land inside 57°S–84°N)", block_mean(cov.astype(float), 4), "share", SEQUENTIAL_BLUE, (0, 1), 2),
    ]

    fig = make_subplots(rows=3, cols=2,
                        subplot_titles=[p[0] for p in panels] + [ZOOM["name"], ""],
                        horizontal_spacing=0.06, vertical_spacing=0.07)
    for i, (title, z, unit, scale, zr, dp) in enumerate(panels):
        r, c = divmod(i, 2)
        xs = coarse_coords(lon, 4) if z.shape[1] != clon.size else clon
        ys = coarse_coords(lat, 4) if z.shape[0] != clat.size else clat
        fig.add_trace(go.Heatmap(
            z=np.round(z, dp), x=xs, y=ys, colorscale=plotly_colorscale(scale),
            zmin=None if zr is None else zr[0], zmax=None if zr is None else zr[1],
            showscale=True, colorbar=dict(len=0.26, y=0.86 - 0.335 * r, x=0.455 + 0.545 * c, thickness=9),
            hovertemplate=f"%{{y:.2f}}, %{{x:.2f}}<br>{unit}: %{{z:.4f}}<extra></extra>"),
            row=r + 1, col=c + 1)

    la0, la1 = ZOOM["lat"]
    lo0, lo1 = ZOOM["lon"]
    si = (lat >= la0) & (lat <= la1)
    sj = (lon >= lo0) & (lon <= lo1)
    fig.add_trace(go.Heatmap(
        z=np.round(med[np.ix_(si, sj)], 5), x=lon[sj], y=lat[si],
        colorscale=plotly_colorscale(SEQUENTIAL_BLUE), showscale=False,
        hovertemplate="%{y:.2f}, %{x:.2f}<br>events yr⁻¹: %{z:.4f}<extra></extra>"), row=3, col=1)

    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=False, zeroline=False)
    for c in (1, 2):
        for r in (1, 2):
            fig.update_xaxes(range=[-180, 180], row=r, col=c)
            fig.update_yaxes(range=[-58, 85], row=r, col=c)
    fig.update_xaxes(range=[lo0, lo1], row=3, col=1)
    fig.update_yaxes(range=[la0, la1], row=3, col=1)
    fig.update_annotations(font_size=12)
    fig.update_layout(height=1180, width=1500, font=dict(family=FONT_STACK, size=11),
                      margin=dict(t=60, b=30, l=40, r=40),
                      plot_bgcolor="white", paper_bgcolor="white", showlegend=False)

    rows = []
    for name, (la, lo_) in REFERENCE_SITES.items():
        i = int(np.argmin(np.abs(lat - la)))
        j = int(np.argmin(np.abs(lon - lo_)))
        rows.append(f"<tr><td>{html.escape(name)}</td><td>{med[i, j]:.4f}</td>"
                    f"<td>{lo[i, j]:.4f} – {hi[i, j]:.4f}</td><td>{pct[i, j]:.1f}</td></tr>")

    zero_share = 100.0 * np.mean(n_ev[cov] <= 0)
    stats = {
        "covered cells (land, 57°S–84°N)": f"{int(cov.sum()):,}",
        "zero tier (model gives no VLH)": f"{zero_share:.2f}% of covered",
        "median rate over covered cells": f"{np.nanmedian(med[cov]):.5f} events yr⁻¹",
        "99th percentile of rate": f"{np.nanpercentile(med[cov], 99):.4f} events yr⁻¹",
        "maximum rate": f"{np.nanmax(med):.4f} events yr⁻¹",
    }
    stat_rows = "".join(f"<tr><td>{html.escape(k)}</td><td>{html.escape(v)}</td></tr>" for k, v in stats.items())

    def attr(name: str) -> str:
        return html.escape(str(ds.attrs.get(name, "")).strip())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "hail-vlh-qa.html"
    body = fig.to_html(full_html=False, include_plotlyjs="cdn", default_width="1500px")
    out.write_text(f"""<!doctype html><meta charset="utf-8">
<title>Very large hail — QA review</title>
<style>
 body {{ font-family: {FONT_STACK}; margin: 28px 40px; color:#1a1a18; max-width: 1560px; }}
 h1 {{ font-size: 21px; margin-bottom: 2px; }} h2 {{ font-size: 15px; margin-top: 26px; }}
 .sub {{ color:#5a5a55; font-size:13px; margin-bottom:18px; }}
 table {{ border-collapse: collapse; font-size: 13px; margin: 8px 0 18px; }}
 td, th {{ border-bottom: 1px solid #e3e2de; padding: 4px 14px 4px 0; text-align: left; }}
 .caveat {{ background:#fbf6ec; border-left:3px solid #eda100; padding:10px 14px; margin:8px 0;
            font-size:13px; max-width: 1100px; }}
 li {{ font-size: 13.5px; margin-bottom: 5px; max-width: 1100px; }}
</style>
<h1>Very large hail (≥ 5 cm) — historical baseline frequency — QA review</h1>
<div class="sub">{attr('source_dataset')}<br>
Layer <code>{html.escape(LAYER)}</code> · contract <code>observational-historical-v1</code> ·
window {attr('temporal_window')} ({attr('n_years')} years) · {attr('qa_reviewed_on')}</div>

<div class="caveat"><b>Display is block-averaged to 0.5°; the file is 0.25°.</b> The zoom panel is
at native resolution so the coarsening cannot hide a lattice artefact.</div>

<h2>What to look for</h2>
<ol>
<li>Does the global maximum sit north-west of Córdoba, Argentina, with the Uruguay/Paraguay/southern
Brazil tri-border, the US Great Plains, north-east South Africa and southern Alberta behind it?
Those are the source paper's own hotspots.</li>
<li><b>Coastline artefacts.</b> The source drops maritime cells beyond 100 km offshore and this layer
masks to land. A rim of high values hugging every coast would be the mask, not weather.</li>
<li>Does the percentile panel read as a hazard field or as the land mask? The zero tier is
{zero_share:.1f}% of covered cells.</li>
<li>The relative-interval panel should be wide where the rate is low and narrow where it is high.
Wide intervals in the hotspots would mean the posterior is wrong.</li>
<li>Grid stripes or blocky steps in the native-resolution zoom would mean the coordinate join or the
padding went wrong.</li>
</ol>

<h2>Reference sites</h2>
<table><tr><th>Site</th><th>events yr⁻¹</th><th>interquartile</th><th>percentile</th></tr>
{''.join(rows)}</table>

<h2>Field summary</h2>
<table>{stat_rows}</table>

<h2>Caveats carried in the file</h2>
<div class="caveat"><b>Threshold.</b> {attr('threshold_caveat')}</div>
<div class="caveat"><b>Baseline window.</b> {attr('baseline_window_caveat')}</div>
<div class="caveat"><b>Resolution.</b> {attr('resolution_caveat')}</div>
<div class="caveat"><b>Calibration.</b> {attr('reporting_bias_caveat')}</div>
<div class="caveat"><b>Coverage.</b> {attr('coverage_caveat')}</div>
<div class="caveat"><b>No trend.</b> {attr('no_trend_rationale')}</div>
<div class="caveat"><b>Trend context, measured.</b> {attr('trend_context_measured')}</div>
<div class="caveat"><b>Interval meaning.</b> {attr('statistic_rationale')}</div>

<h2>Maps</h2>
{body}
""", encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
