#!/usr/bin/env python3
"""QA maps for the hail severity-response fields built by `build_hail_severity_fields.py`.

NOT a customer dashboard and NOT a `generate_maps.py` product: that tool renders layers that
satisfy OUTPUT-SPEC (six tabs keyed to `median`/`percentile`/slopes), and these fields have
no decade axis, no percentile and no trend, so it has nothing to render. This is the QA view
for judging whether the spatial pattern is physical before anyone decides what to do with it.

One row per quantity: the historical panel on a sequential scale, then the per-scenario
change on a diverging scale centred at zero with a symmetric limit, both from
`scripts/utils/viz_common.py` so the vocabulary matches the rest of the product.

Cells below the event floor are already masked in the input file and render as gaps. Gaps
are the honest rendering: the event population is satellite-detected storms, so an empty
cell means "not sampled", never "no hail".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.viz_common import (  # noqa: E402
    FONT_STACK,
    SEQUENTIAL_BLUE,
    diverging_colorscale,
    plotly_colorscale,
    symmetric_limit,
)

# Pooled first (the primary fields), then the storm-weighted pair, because the two
# weightings disagree in sign over several regions and the QA view exists to show that.
QUANTITIES = [
    ("mean_ke", "Mean KE per landed stone — pooled over stones", "J"),
    ("mean_diameter", "Mean diameter per landed stone — pooled over stones", "mm"),
    ("production", "Share of embryos producing hail at the ground", "fraction"),
    ("ake", "Mean KE per seeded embryo — blended (production x intensity)", "J"),
    ("storm_mean_ke", "Mean KE — storm-weighted (each storm counts once)", "J"),
    ("storm_mean_diameter", "Mean diameter — storm-weighted", "mm"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fields", type=Path, default=Path("reports/maps/hail-severity/hail_severity_10deg.nc"))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if not args.fields.exists():
        print(f"fields file not found: {args.fields}", file=sys.stderr)
        return 2
    ds = xr.open_dataset(args.fields)
    scenarios = [str(s) for s in ds["scenario"].values]
    out = args.out or args.fields.with_suffix(".html")

    ncols = 1 + len(scenarios)
    titles = []
    for _, label, unit in QUANTITIES:
        titles.append(f"{label} — historical ({unit})")
        titles.extend(f"{label} — {sc} change (%)" for sc in scenarios)
    fig = make_subplots(rows=len(QUANTITIES), cols=ncols, subplot_titles=titles,
                        horizontal_spacing=0.02, vertical_spacing=0.07)

    lat = ds["lat"].values
    lon = ds["lon"].values
    clamped_counts: list[tuple[str, float, int]] = []
    for r, (key, label, unit) in enumerate(QUANTITIES, start=1):
        hist = ds[f"{key}_historical"].isel(scenario=0).values
        fig.add_trace(go.Heatmap(z=hist, x=lon, y=lat, colorscale=plotly_colorscale(SEQUENTIAL_BLUE),
                                 showscale=False, hovertemplate=f"%{{y}}, %{{x}}<br>{unit}: %{{z:.2f}}<extra></extra>"),
                      row=r, col=1)
        chg = ds[f"{key}_change_pct"].values
        # symmetric_limit returns the clamped COUNT as well, and viz_common requires it to be
        # reported rather than discarded -- a clamped cell is drawn at the endpoint colour and
        # would otherwise read as "at the limit" when it is beyond it.
        limit, clamped = symmetric_limit(chg[np.isfinite(chg)])
        clamped_counts.append((label, limit, clamped))
        for c, sc in enumerate(scenarios, start=2):
            z = ds[f"{key}_change_pct"].sel(scenario=sc).values
            lo = ds[f"{key}_change_pct_member_min"].sel(scenario=sc).values
            hi = ds[f"{key}_change_pct_member_max"].sel(scenario=sc).values
            n = ds["n_events"].sel(scenario=sc).values
            fig.add_trace(go.Heatmap(
                z=z, x=lon, y=lat, zmid=0, zmin=-limit, zmax=limit,
                colorscale=diverging_colorscale("light"), showscale=(c == ncols and r == 1),
                customdata=np.dstack([lo, hi, n]),
                hovertemplate=("%{y}, %{x}<br>change: %{z:+.1f}%"
                               "<br>member range: %{customdata[0]:+.1f} to %{customdata[1]:+.1f}%"
                               "<br>events: %{customdata[2]:.0f}<extra></extra>")),
                row=r, col=c)

    fig.update_xaxes(showgrid=False, zeroline=False, range=[-180, 180])
    fig.update_yaxes(showgrid=False, zeroline=False, range=[-60, 75])
    fig.update_annotations(font_size=11)
    fig.update_layout(
        height=300 * len(QUANTITIES), width=460 * ncols,
        font=dict(family=FONT_STACK, size=11),
        title=dict(text="Hail severity response — conditional on a severe hailstorm occurring "
                        "(EC-Earth3, 20 members, median; NOT a frequency field)", font=dict(size=15)),
        margin=dict(t=110, b=40, l=40, r=40), plot_bgcolor="white", paper_bgcolor="white")

    clamp_note = "; ".join(f"{lab}: ±{lim:.0f}%, {n} cells clamped" for lab, lim, n in clamped_counts)
    note = (f"Source: {ds.attrs.get('source','')}<br>"
            f"{ds.attrs.get('event_conditioning','')}<br>"
            f"{ds.attrs.get('frequency_caveat','')}<br>"
            f"{ds.attrs.get('embryo_support','')}<br>"
            f"n_models: {ds.attrs.get('n_models','')}<br>"
            f"diverging limits (95th pct of |change|) — {clamp_note}")
    fig.add_annotation(text=note, xref="paper", yref="paper", x=0, y=1.045,
                       showarrow=False, align="left", font=dict(size=10, color="#5a5a55"))

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"wrote {out} ({out.stat().st_size/1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
