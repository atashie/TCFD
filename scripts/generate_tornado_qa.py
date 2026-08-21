#!/usr/bin/env python3
"""QA maps for the CONUS observational tornado layers.

WHY THIS IS NOT generate_maps.py
    `generate_maps.py` is and stays the tool for gridded TCFD layers, and the two are
    deliberately not merged. But it is bound to the OUTPUT-SPEC decadal contract: it
    selects on `ds.decade`, builds a Trend tab from `ols_slope`/`sen_slope` and a
    Members tab from `n_members`, none of which exist on an
    `observational-historical-v1` layer. Teaching it a second contract would put a
    branch through the tool that renders 23 shipped layers, to serve 8 that have no
    decade axis at all. This is a narrow QA renderer for that one contract instead.

WHAT IT IS FOR
    `qa_reviewed_on` is null in every tornado layer, and the registry comment is
    explicit that the date goes in only after a human has actually looked. Nobody
    could look, because there were no maps. This makes the maps. It does not, and
    cannot, perform the review.

WHAT TO LOOK FOR (the review this page is meant to support)
    1. Does the Plains/Dixie maximum sit where tornado climatology puts it, and does
       it MOVE between rungs? On `all` it should include metro artefacts; on `f2plus`
       those should vanish while Oklahoma/Kansas/Alabama persist.
    2. Are there straight-line or round-coordinate artefacts? At 0.25 deg the
       pre-1976 geocoding quantization should be absorbed. Visible grid stripes mean
       it is not, and the resolution argument needs revisiting.
    3. Does the percentile panel look like a hazard field or like a coastline? The
       zero tier is 31-86% of CONUS depending on rung; a percentile map that is
       mostly flat is telling you the tier boundary is doing all the work.
    4. State-border discontinuities. Reporting practice varies by NWS office, and a
       visible edge along a state line is an observing artefact, never weather.

USAGE
    python3 scripts/generate_tornado_qa.py                    # both windows
    python3 scripts/generate_tornado_qa.py --window full
"""

from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent))
from utils import viz_common as vc  # noqa: E402

PROCESSED = Path("data/processed")
OUT_DIR = Path("reports/maps/tornado")     # QA/QC maps live in reports/maps/{hazard}/
RUNGS = ["all", "f1plus", "f2plus", "f3plus"]

METRICS = [
    ("median", "Crossing rate", "per 1e4 km2 per yr", "rate"),
    ("percentile", "Percentile (two-tier)", "1-100", "pct"),
    ("n_events", "Paths crossing cell", "count", "count"),
]


def assert_colorbars_clear(fig, cbar_y: float, thickness_px: int, top_margin_px: int) -> None:
    """Fail the build if any colourbar overlaps a map panel, or overflows the margin.

    This exists because the first version of this page shipped with the colourbars
    sitting on top of the first row of maps, which is invisible in the code and
    obvious only once a human opens the file. A geometry bug that only a reviewer can
    see is exactly the kind worth asserting.
    """
    tops = []
    for key, ax in fig.layout.to_plotly_json().items():
        if key.startswith("yaxis") and isinstance(ax, dict) and "domain" in ax:
            tops.append(ax["domain"][1])
    highest_panel = max(tops) if tops else 1.0
    if cbar_y <= highest_panel:
        raise SystemExit(
            f"colourbar y={cbar_y} sits at or below the top panel edge ({highest_panel}); "
            "it would be drawn over the maps."
        )
    plot_h = fig.layout.height - top_margin_px - fig.layout.margin.b
    needed = (cbar_y - 1.0) * plot_h + thickness_px + 40  # ticks + title headroom
    if needed > top_margin_px - 46:  # 46px reserved for the figure title
        raise SystemExit(
            f"colourbar needs {needed:.0f}px above the plot area but the top margin is "
            f"{top_margin_px}px; raise top_margin or lower CBAR_Y."
        )


def layer_dirs(window: str):
    out = []
    for rung in RUNGS:
        d = PROCESSED / f"tornado-spc_{rung}_{window}"
        # Scenario token is `observed` (observational-historical-v1 contract), renamed from
        # `historical` on 2026-08-18 minutes after this page was last rendered -- the stale
        # token left this renderer silently finding nothing (caught at QA review 2026-08-20).
        f = d / f"tornado-{rung}_observed_processed.nc"
        if f.exists():
            out.append((rung, f))
    return out


def build(window: str) -> Path | None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    layers = layer_dirs(window)
    if not layers:
        print(f"  no layers found for window={window!r}")
        return None

    # CONUS crop. The fields are NaN outside the domain, and on a go.Heatmap every global
    # NaN row still costs bytes in the base64 payload -- measured 153.7 MB/page uncropped
    # float64 vs ~4 MB cropped float32 (the same dtype/crop rules as the landslide page).
    with xr.open_dataset(layers[0][1]) as ds0:
        fin0 = np.isfinite(ds0[METRICS[0][0]].values)
        rows = np.where(fin0.any(axis=1))[0]
        cols = np.where(fin0.any(axis=0))[0]
        r0, r1 = max(rows.min() - 2, 0), min(rows.max() + 3, fin0.shape[0])
        c0, c1 = max(cols.min() - 2, 0), min(cols.max() + 3, fin0.shape[1])

    scale_rate = vc.plotly_colorscale(vc.SEQUENTIAL_RED)
    scale_pct = vc.plotly_colorscale(vc.ORDINAL_RISK)

    # Shared upper limit per metric so rungs are visually comparable. p99 rather than
    # max: a single extreme cell would otherwise flatten every other rung to one colour.
    caps = {}
    for key, *_ in METRICS:
        vals = []
        for _, f in layers:
            with xr.open_dataset(f) as ds:
                v = ds[key].values
                vals.append(v[np.isfinite(v)])
        allv = np.concatenate(vals) if vals else np.array([0.0])
        caps[key] = float(np.quantile(allv, 0.99)) or 1.0

    fig = make_subplots(
        rows=len(layers), cols=len(METRICS),
        subplot_titles=[f"{rung} — {label}" for rung, _ in layers for _, label, _, _ in METRICS],
        horizontal_spacing=0.06, vertical_spacing=0.05,
    )

    # COLOURBAR PLACEMENT.
    # The first version put the three colourbars at paper y=0.89, which is INSIDE the
    # first row of maps (row 1 spans y-domain [0.7875, 1.0]) -- they sat on top of the
    # `all` rung's panels and obscured exactly the maps a reviewer looks at first.
    # They now sit in the top margin, above the subplot titles at y=1.0, anchored from
    # their bottom edge so they grow upward into the margin rather than down into the
    # plot. Column centres are read from the axis domains rather than hardcoded, so
    # they stay correct if the spacing changes. assert_colorbars_clear() below fails
    # the build if anything ever drifts back over a panel.
    col_centre = []
    for c in range(1, len(METRICS) + 1):
        ax = "xaxis" if c == 1 else f"xaxis{c}"
        d = fig.layout[ax].domain
        col_centre.append((d[0] + d[1]) / 2)

    CBAR_Y = 1.045          # paper coords; 1.0 is the top of the plotting area
    CBAR_THICKNESS = 11

    ref_rows = []
    caveats = {}

    for r, (rung, f) in enumerate(layers, start=1):
        ds = xr.open_dataset(f)
        lat = ds.lat.values[r0:r1]
        lon = ds.lon.values[c0:c1]
        if not caveats:
            caveats = {k: ds.attrs.get(k, "") for k in
                       ("statistic", "statistic_rationale", "resolution_caveat",
                        "reporting_bias_caveat", "coverage_caveat", "no_trend_rationale",
                        "unrated_handling", "output_contract")}

        for c, (key, label, units, kind) in enumerate(METRICS, start=1):
            # float32, not round: plotly>=6 serialises numpy as base64 typed arrays, so
            # rounding changes the payload by exactly zero and dtype is the whole term.
            z = ds[key].values[r0:r1, c0:c1].astype(np.float32)
            fig.add_trace(
                go.Heatmap(
                    z=z, x=lon, y=lat,
                    colorscale=scale_pct if kind == "pct" else scale_rate,
                    zmin=1 if kind == "pct" else 0,
                    zmax=100 if kind == "pct" else caps[key],
                    showscale=(r == 1),
                    colorbar=dict(len=0.24, thickness=CBAR_THICKNESS,
                                  x=col_centre[c - 1], xanchor="center",
                                  y=CBAR_Y, yanchor="bottom",
                                  orientation="h",
                                  tickfont=dict(size=9),
                                  title=dict(text=f"{label} [{units}]", side="top",
                                             font=dict(size=10)))
                    if r == 1 else None,
                    hovertemplate=f"{label}: %{{z}}<br>%{{y:.2f}}N %{{x:.2f}}E<extra></extra>",
                ),
                row=r, col=c,
            )
            fig.update_yaxes(scaleanchor=f"x{'' if (r-1)*3+c == 1 else (r-1)*3+c}",
                             scaleratio=1.0, row=r, col=c)

        ref_rows.append((rung, ds.attrs.get("cells_with_at_least_one_crossing"),
                         ds.attrs.get("cells_in_mask"),
                         ds.attrs.get("zero_tier_share_of_mask"),
                         ds.attrs.get("magnitude_threshold")))
        ds.close()

    top_margin = 200
    fig.update_layout(
        height=300 * len(layers) + top_margin + 40, width=1500,
        title=dict(text=f"CONUS tornado hazard — QA — window {window}", x=0.01,
                   y=0.985, yref="container", yanchor="top",
                   font=dict(size=18, family=vc.FONT_STACK)),
        font=dict(family=vc.FONT_STACK, size=11),
        margin=dict(l=40, r=20, t=top_margin, b=40),
        paper_bgcolor="white", plot_bgcolor="#f0efec",
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)

    assert_colorbars_clear(fig, CBAR_Y, CBAR_THICKNESS, top_margin)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"tornado-qa-{window}.html"

    body = fig.to_html(full_html=False, include_plotlyjs="cdn", default_width="1500px")

    tbl = "".join(
        f"<tr><td><code>{html.escape(str(a))}</code></td><td>{html.escape(str(e))}</td>"
        f"<td>{html.escape(str(m))}</td><td>{html.escape(str(z))}%</td>"
        f"<td>{html.escape(str(t))}</td></tr>"
        for a, e, m, z, t in ref_rows
    )
    cav = "".join(
        f"<h3>{html.escape(k)}</h3><p>{html.escape(str(v))}</p>"
        for k, v in caveats.items() if str(v).strip()
    )

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Tornado QA — {html.escape(window)}</title>
<style>
 body{{font-family:{vc.FONT_STACK};margin:24px;max-width:1560px;color:#1a1a18;background:#fff}}
 h1{{font-size:22px}} h2{{font-size:17px;margin-top:32px;border-bottom:1px solid #ddd;padding-bottom:4px}}
 h3{{font-size:13px;margin-bottom:2px;color:#555;font-family:ui-monospace,monospace}}
 p{{margin-top:2px;line-height:1.5;max-width:95ch}}
 table{{border-collapse:collapse;margin:12px 0}} td,th{{border:1px solid #ddd;padding:5px 10px;font-size:13px;text-align:left}}
 th{{background:#f0efec}} .warn{{background:#fff4f4;border-left:4px solid #e34948;padding:10px 14px;margin:16px 0}}
</style></head><body>
<h1>CONUS tornado hazard — QA review — window {html.escape(window)}</h1>
<div class="warn"><strong>This page supports a review; it is not one.</strong>
<code>qa_reviewed_on</code> stays null until a human has read these maps and the caveats below.
A layer passing its contract check means it is <em>shaped</em> right, not that it means what its name says.</div>

<h2>What to look for</h2>
<ol>
<li>Does the maximum sit where tornado climatology puts it, and does it <em>move</em> between rungs?
    Metro artefacts should be visible on <code>all</code> and gone by <code>f2plus</code>.</li>
<li>Grid stripes or round-coordinate hot spots — the pre-1976 geocoding quantization should be
    absorbed at 0.25°. If it is visible, the resolution argument needs revisiting.</li>
<li>Does the percentile panel read as a hazard field or as a coastline? The zero tier is
    31–86% of CONUS depending on rung.</li>
<li>State-border discontinuities — reporting practice varies by NWS office. An edge along a
    state line is an observing artefact, never weather.</li>
</ol>

<h2>Rung summary</h2>
<table><tr><th>rung</th><th>occupied cells</th><th>cells in mask</th><th>zero tier</th><th>threshold</th></tr>
{tbl}</table>

<h2>Maps</h2>
{body}

<h2>Caveats carried in the files' own attributes</h2>
{cav}
</body></html>"""

    out.write_text(page)
    print(f"  wrote {out}  ({out.stat().st_size/1e6:.1f} MB)")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", default=None, help="'full', 'from1996', or omit for both")
    args = ap.parse_args()

    windows = [args.window] if args.window else ["full", "from1996"]
    print("Generating tornado QA pages")
    made = [build(w) for w in windows]
    made = [m for m in made if m]
    if not made:
        print("Nothing generated. Run scripts/process_tornado_spc.py first.")
        return 1
    print(f"\n{len(made)} page(s). qa_reviewed_on remains null until a human reads them.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
