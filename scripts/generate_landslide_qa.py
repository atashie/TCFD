#!/usr/bin/env python3
"""QA maps for the global observational landslide hazard layer.

WHY THIS IS NOT generate_maps.py
    Same reason `generate_tornado_qa.py` is not: `generate_maps.py` is bound to the
    OUTPUT-SPEC decadal contract -- it selects on `ds.decade`, builds a Trend tab from
    `ols_slope`/`sen_slope` and a Members tab from `n_members`, none of which exist on
    an `observational-historical-v1` layer. This is a narrow renderer for that contract.

WHY IT IS NOT generate_tornado_qa.py EITHER
    That one is hardcoded to four magnitude rungs on a CONUS mask and renders three
    metrics that this layer does not have (`n_events`). This layer is one global panel
    set whose interesting structure is the relationship between four fields -- the
    conditional rate, the areal mean it is NOT ranked on, the hazard-bearing area
    fraction, and the percentile. Merging the two would put a branch through a working
    renderer to serve a different shape.

WHAT IT IS FOR
    `qa_reviewed_on` is null. The registry comment is explicit that the date goes in
    only after a human has actually looked. This makes the maps. It does not, and
    cannot, perform the review.

WHAT TO LOOK FOR (the review this page is meant to support)
    1. THE ONE THAT MATTERS: `median` and `percentile` are deliberately NOT ranked on
       the same quantity (Spearman 0.34). Compare panels 1 and 4. Cells that are dark
       in `percentile` and pale in `median` should be places that are mostly flat with
       a genuinely hazardous minority -- large river valleys, coastal plains backed by
       scarps. If that pattern instead looks arbitrary, the ranking choice is wrong.
    2. Do the maxima sit where landslide climatology puts them -- the Himalayan arc,
       the Andes, the Philippines/Indonesia, Japan, the Alps, the Apennines, the East
       African Rift, Central America? A maximum in a flat region is a defect.
    3. `hazard_area_fraction` should look like a terrain map, because it essentially is
       one. If it correlates with country boundaries instead, the source model's
       training inventory is showing through and the reporting-bias caveat understates
       the problem.
    4. Coastlines. The published mask is the union of the ISIMIP3b land mask (upsampled
       0.5 -> 0.25 deg, so blocky) with any cell carrying hazard. Look for a fringe of
       published cells sitting in open water, and for hazard-bearing land dropped at the
       coast. Neither is fatal at 0.25 deg; both bound how far the layer can be trusted
       near a shoreline.
    5. The 60 S / 72 N edges are the SOURCE extent, not our choice. They should be clean
       horizontal cuts. Anything ragged means the row indexing drifted.

USAGE
    python3 scripts/generate_landslide_qa.py
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

LAYER = Path("data/processed/landslide-arup_rf-median_hist/landslide-arup_observed_processed.nc")
OUT_DIR = Path("reports/maps/landslide")   # QA/QC maps live in reports/maps/{hazard}/

#: (variable, label, units, colour kind). Order is the review order in the docstring.
METRICS = [
    ("median", "Rate on hazard-bearing ground", "km-2 yr-1", "rate"),
    ("areal_mean_rate", "Areal mean rate (what percentile ranks on)", "km-2 yr-1", "rate"),
    ("hazard_area_fraction", "Hazard-bearing area fraction", "0-1", "frac"),
    ("percentile", "Percentile (two-tier)", "1-100", "pct"),
]

CAVEAT_ATTRS = (
    "output_contract", "statistic", "statistic_rationale", "percentile_basis",
    "resolution_caveat", "reporting_bias_caveat", "coverage_caveat",
    "no_trend_rationale", "no_scenario_rationale", "source_zero_is_ambiguous",
    "land_mask", "median_branch_measured", "source_licence", "source_licence_status",
    "attribution_required",
)

#: Sites the layer is checked at, and why each one is here. Reproduced on the page so a
#: reviewer sees the numbers next to the maps rather than having to trust a docstring.
REF_SITES = [
    ("Baguio PH", 16.41, 120.60, "typhoon-driven slides, expected near the global max"),
    ("Wenchuan CN", 31.00, 103.40, "steep Longmenshan front"),
    ("Medellin CO", 6.25, -75.56, "Andean city, frequent fatal slides"),
    ("Shimla IN", 31.10, 77.17, "Himalayan front"),
    ("Bergen NO", 60.39, 5.32, "steep and very wet, high latitude"),
    ("Cusco PE", -13.53, -71.97, "Andes"),
    ("Kathmandu NP", 27.71, 85.32, "Himalaya, but a broad flat valley floor"),
    ("Apennines IT", 43.00, 12.50, "most landslide-mapped terrain in Europe"),
    ("Freetown SL", 8.48, -13.23, "2017 Regent slide"),
    ("Rio de Janeiro BR", -22.91, -43.18, "coastal cell, mostly water -- resolution probe"),
    ("Amsterdam NL", 52.37, 4.90, "control: flat, must be tier 1"),
    ("Des Moines IA", 41.59, -93.62, "control: flat, must be tier 1"),
    ("Cairo EG", 30.04, 31.24, "control: arid, one escarpment pixel"),
    ("mid-Pacific", 0.00, -160.00, "control: ocean, must be NaN"),
]


def assert_colorbars_clear(fig, cbar_y: float, thickness_px: int, top_margin_px: int) -> None:
    """Fail the build if a colourbar would be drawn over a map panel.

    Lifted from generate_tornado_qa.py deliberately rather than imported: that module is
    a script, not a library, and a geometry bug only a reviewer can see is worth
    asserting in both places.
    """
    tops = [ax["domain"][1] for key, ax in fig.layout.to_plotly_json().items()
            if key.startswith("yaxis") and isinstance(ax, dict) and "domain" in ax]
    highest = max(tops) if tops else 1.0
    if cbar_y <= highest:
        raise SystemExit(f"colourbar y={cbar_y} is at or below the top panel edge ({highest})")
    plot_h = fig.layout.height - top_margin_px - fig.layout.margin.b
    needed = (cbar_y - 1.0) * plot_h + thickness_px + 40
    if needed > top_margin_px - 46:
        raise SystemExit(f"colourbar needs {needed:.0f}px, top margin is {top_margin_px}px")


def build() -> Path:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if not LAYER.exists():
        raise SystemExit(f"missing {LAYER} -- run scripts/process_landslide_arup.py first")

    ds = xr.open_dataset(LAYER)
    # Crop to the SOURCE extent plus a margin, for payload only. The file itself stays
    # global; this keeps the 72 N / 59.5 S edges visible (review point 5) while dropping
    # ~190 rows that are NaN in every panel and cost ~25% of the page.
    ds = ds.sel(lat=slice(-66.0, 76.0))
    lat, lon = ds.lat.values, ds.lon.values

    scale_rate = vc.plotly_colorscale(vc.SEQUENTIAL_RED)
    scale_pct = vc.plotly_colorscale(vc.ORDINAL_RISK)
    scale_frac = vc.plotly_colorscale(vc.SEQUENTIAL_BLUE)

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[f"{lbl}  [{u}]" for _, lbl, u, _ in METRICS],
        horizontal_spacing=0.07, vertical_spacing=0.10,
    )

    CBAR_Y, CBAR_THICK = 1.045, 11
    top_margin = 190

    #: PAYLOAD: THE TERM IS DTYPE, NOT DECIMAL PLACES. Measured 2026-08-19 -- plotly 6.9
    #: serialises numpy arrays as base64 TYPED ARRAYS (`bdata`), not as JSON text, so
    #: rounding the values changes the page size by nothing at all. The first version of
    #: this file claimed the opposite in a comment and rounded to fix it; the page went
    #: 50.4 -> 39.8 MB purely from the latitude crop above, and `"null"` appears once in
    #: the whole document. Cost is bytes-per-cell: 818k cells x 4 panels x 8 bytes
    #: (float64) x 4/3 (base64) = 35 MB. float32 halves it and is far more precision than
    #: a colour ramp can show. Rounding is kept only as a guard in case a future plotly
    #: falls back to text encoding.
    DECIMALS = {"median": 5, "areal_mean_rate": 6, "hazard_area_fraction": 3, "percentile": 1}

    for i, (key, label, units, kind) in enumerate(METRICS):
        r, c = i // 2 + 1, i % 2 + 1
        z = np.round(ds[key].values.astype(np.float32), DECIMALS[key]).astype(np.float32)

        if kind == "pct":
            scale, zmin, zmax = scale_pct, 1, 100
        elif kind == "frac":
            scale, zmin, zmax = scale_frac, 0, 1
        else:
            # p99 of the OCCUPIED values, not the max: one extreme cell (0.72) would
            # otherwise flatten the entire map to the bottom of the colour ramp.
            finite = z[np.isfinite(z) & (z > 0)]
            scale, zmin = scale_rate, 0
            zmax = float(np.quantile(finite, 0.99)) if finite.size else 1.0

        fig.add_trace(
            go.Heatmap(
                z=z, x=lon, y=lat, colorscale=scale, zmin=zmin, zmax=zmax,
                colorbar=dict(len=0.30, thickness=CBAR_THICK,
                              x=[0.21, 0.79][c - 1], xanchor="center",
                              y=CBAR_Y, yanchor="bottom", orientation="h",
                              tickfont=dict(size=9),
                              title=dict(text=f"{label} [{units}]", side="top",
                                         font=dict(size=10))) if r == 1 else None,
                showscale=(r == 1),
                hovertemplate=f"{label}: %{{z}}<br>%{{y:.2f}}N %{{x:.2f}}E<extra></extra>",
            ),
            row=r, col=c,
        )
        n = i + 1
        fig.update_yaxes(scaleanchor=f"x{'' if n == 1 else n}", scaleratio=1.0, row=r, col=c)

    fig.update_layout(
        height=430 * 2 + top_margin, width=1600,
        title=dict(text="Global rainfall-triggered landslide hazard (Arup 1980-2018) "
                        "— QA — 0.25 deg",
                   x=0.01, y=0.985, yref="container", yanchor="top",
                   font=dict(size=18, family=vc.FONT_STACK)),
        font=dict(family=vc.FONT_STACK, size=11),
        margin=dict(l=40, r=20, t=top_margin, b=40),
        paper_bgcolor="white", plot_bgcolor="#f0efec",
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    assert_colorbars_clear(fig, CBAR_Y, CBAR_THICK, top_margin)

    # --- reference-site table, read from the file at run time ---
    rows = []
    for name, la, lo, why in REF_SITES:
        s = ds.sel(lat=la, lon=lo, method="nearest")
        fmt = lambda v: "NaN" if not np.isfinite(v) else f"{v:.5f}"  # noqa: E731
        rows.append((name, why,
                     fmt(float(s["median"].values)),
                     fmt(float(s["areal_mean_rate"].values)),
                     fmt(float(s["hazard_area_fraction"].values)),
                     fmt(float(s["percentile"].values))))

    caveats = {k: ds.attrs.get(k, "") for k in CAVEAT_ATTRS}
    summary = {k: ds.attrs.get(k, "") for k in
               ("temporal_window", "n_years", "spatial_resolution", "native_resolution",
                "cells_published", "cells_occupied", "zero_tier_share_of_mask", "units")}
    ds.close()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "landslide-qa.html"
    body = fig.to_html(full_html=False, include_plotlyjs="cdn", default_width="1600px")

    esc = html.escape
    ref_tbl = "".join(
        f"<tr><td>{esc(n)}</td><td class='why'>{esc(w)}</td><td>{m}</td><td>{a}</td>"
        f"<td>{f}</td><td><strong>{p}</strong></td></tr>"
        for n, w, m, a, f, p in rows)
    sum_tbl = "".join(f"<tr><td><code>{esc(k)}</code></td><td>{esc(str(v))}</td></tr>"
                      for k, v in summary.items())
    cav = "".join(f"<h3>{esc(k)}</h3><p>{esc(str(v))}</p>"
                  for k, v in caveats.items() if str(v).strip())

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Landslide QA — Arup 1980-2018 @ 0.25°</title>
<style>
 body{{font-family:{vc.FONT_STACK};margin:24px;max-width:1660px;color:#1a1a18;background:#fff}}
 h1{{font-size:22px}} h2{{font-size:17px;margin-top:32px;border-bottom:1px solid #ddd;padding-bottom:4px}}
 h3{{font-size:13px;margin-bottom:2px;color:#555;font-family:ui-monospace,monospace}}
 p{{margin-top:2px;line-height:1.5;max-width:100ch}}
 table{{border-collapse:collapse;margin:12px 0}}
 td,th{{border:1px solid #ddd;padding:5px 10px;font-size:13px;text-align:left}}
 th{{background:#f0efec}} td.why{{color:#666;font-size:12px}}
 .warn{{background:#fff4f4;border-left:4px solid #e34948;padding:10px 14px;margin:16px 0}}
 .note{{background:#f4f7ff;border-left:4px solid #2a78d6;padding:10px 14px;margin:16px 0}}
</style></head><body>
<h1>Global rainfall-triggered landslide hazard — QA review</h1>

<div class="warn"><strong>This page supports a review; it is not one.</strong>
<code>qa_reviewed_on</code> stays null until a human has read these maps and the caveats below.
Passing the contract check means the file is <em>shaped</em> right, not that it means what its
name says — both withdrawn sugarcane layers passed every algebraic check.</div>

<div class="note"><strong>Read this before panel 1.</strong> <code>median</code> is the rate on
the cell's <em>hazard-bearing</em> ground; <code>percentile</code> ranks on
<code>areal_mean_rate</code> over the <em>whole</em> cell. Spearman between the two orderings is
0.34, so they are close to independent by design. <strong>A cell can show
<code>median</code>&nbsp;=&nbsp;0 at a high percentile</strong> — most of it is flat, the cell as
a whole still carries substantial landslide activity. That is the intended behaviour, not a bug.</div>

<h2>What to look for</h2>
<ol>
<li><strong>Panels 1 vs 4.</strong> Dark in <code>percentile</code>, pale in <code>median</code>
    should mean "mostly flat with a hazardous minority" — river valleys, plains backed by
    scarps. If the pattern looks arbitrary, the ranking choice is wrong.</li>
<li>Do the maxima sit where landslide climatology puts them — Himalayan arc, Andes,
    Philippines/Indonesia, Japan, Alps, Apennines, East African Rift, Central America?
    A maximum in flat terrain is a defect.</li>
<li><code>hazard_area_fraction</code> should look like a terrain map. If it tracks country
    boundaries instead, the source model's training inventory is showing through and the
    reporting-bias caveat understates the problem.</li>
<li>Coastlines: the mask is the ISIMIP3b land mask (upsampled, so blocky) ∪ any hazard-bearing
    cell. Look for published cells in open water, and hazard-bearing land dropped at the coast.</li>
<li>The 60°S / 72°N edges are the <em>source</em> extent and should be clean horizontal cuts.
    Ragged edges mean the row indexing drifted.</li>
</ol>

<h2>Maps</h2>
{body}

<h2>Reference sites — read from the file at run time</h2>
<table><tr><th>site</th><th>why it is here</th><th>median</th><th>areal mean</th>
<th>hazard frac</th><th>percentile</th></tr>{ref_tbl}</table>
<p>Rio de Janeiro is the resolution probe: a coastal cell that is mostly water and flat, so a
mid-range percentile there is the 28&nbsp;km support showing, not a modelling error.</p>

<h2>Layer summary</h2>
<table><tr><th>attribute</th><th>value</th></tr>{sum_tbl}</table>

<h2>Caveats and declarations carried in the file's own attributes</h2>
{cav}
</body></html>"""

    out.write_text(page)
    print(f"  wrote {out}  ({out.stat().st_size / 1e6:.1f} MB)")
    return out


def main() -> int:
    argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    ).parse_args()
    print("Generating landslide QA page")
    build()
    print("\nqa_reviewed_on stays null until a human reads it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
