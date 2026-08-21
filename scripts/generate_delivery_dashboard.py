#!/usr/bin/env python
"""Build the QA/QC dashboard for a customer delivery.

    python scripts/generate_delivery_dashboard.py deliveries/acme-timber/20260812

Writes `dashboard.html` into the delivery folder. One self-contained page (Plotly from CDN)
with four views over the same filtered slice:

    Map          one marker per site, metric toggle-able (value / percentile / trend)
    Summaries    bar rollups by hazard, by asset class, and a portfolio band histogram
    Time series  one location at a time, small multiples per hazard, scenarios as series
    Table        every filtered row, sortable, with a text search

WHY A SEPARATE SCRIPT FROM generate_maps.py
-------------------------------------------
`generate_maps.py` renders GRIDDED layers: ~70,849 SVG markers per panel, where marker
count is the binding constraint and the output is a multi-file collection. This renders a
POINT portfolio: hundreds of sites, where the whole delivery embeds as a few hundred KB of
JSON and interactivity is the point. Same conventions, opposite payload profile. Shared
vocabulary lives in `scripts/utils/viz_common.py`; the two scripts should not be merged.

SCOPING OF THE FILTER ROW
-------------------------
One filter row scopes everything, per the dataviz rule against per-chart filters. Two
deliberate exceptions, both because the view's whole purpose is to vary that dimension:

  * Time series ignores the DECADE filter -- it plots every decade.
  * Summaries ignore the HAZARD filter -- "impacts by hazard" is the chart.
  * The Values table (since 2026-08-21) answers to its OWN controls only -- per-column
    dropdowns plus its search box -- because bar filters silently intersecting with
    column filters produced empty results that read as broken filters.

All are stated on the page so a reader is never guessing what a chart is showing.
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from html import escape as html_escape
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.viz_common import (
    is_forcing_scenario,  # noqa: E402
    FONT_STACK, ORDINAL_RISK, ORDINAL_RISK_DARK, PLOTLY_CDN, RISK_BAND_ORDER,
    SCENARIO_TIER, SEQUENTIAL_BLUE, SEQUENTIAL_RED, SEQUENTIAL_RED_DARK,
    TIER_COLOR_DARK, TIER_COLOR_LIGHT, TIER_DASH, TIER_LABELS, TIER_ORDER, TIER_SYMBOL,
    check_tier_collisions, css_tokens, diverging_colorscale, plotly_colorscale, sigfig,
    tier_of,
)


def load_delivery(path: Path):
    need = ["locations.csv", "assets.csv", "layers.csv", "values.csv"]
    missing = [n for n in need if not (path / n).exists()]
    if missing:
        raise SystemExit(f"Not a delivery folder ({path}): missing {', '.join(missing)}")
    scores = (pd.read_csv(path / "climate_score.csv")
              if (path / "climate_score.csv").exists() else pd.DataFrame())
    return (
        pd.read_csv(path / "locations.csv"),
        pd.read_csv(path / "assets.csv"),
        pd.read_csv(path / "layers.csv"),
        pd.read_csv(path / "values.csv"),
        scores,
        json.loads((path / "manifest.json").read_text())
        if (path / "manifest.json").exists() else {},
    )


def build_payload(locations, assets, layers, values, scores):
    """Compact JSON for the browser.

    Values are keyed on (location, layer, scenario, decade), NOT on asset: extraction is
    purely coordinate-driven, so two assets at one site read identically for a given layer.
    Carrying them per asset would multiply the payload for no information. Asset identity
    stays in `assets` and is joined back for the asset-class rollups.
    """
    assets = assets.copy()
    assets["layer_ids"] = assets["layer_ids"].fillna("").astype(str)

    aid_to_loc = dict(zip(assets["asset_id"], assets["location_id"]))
    vals = values.copy()
    vals["location_id"] = vals["asset_id"].map(aid_to_loc)

    keep = ["location_id", "layer_id", "scenario", "decade", "value", "lower_ci",
            "upper_ci", "percentile", "ols_slope", "sen_slope", "slopes_agree",
            "n_members", "n_models", "data_status"]
    dedup = vals[keep].drop_duplicates(["location_id", "layer_id", "scenario", "decade"])

    def clean(x):
        if pd.isna(x):
            return None
        return float(sigfig(np.array([float(x)]))[0]) if isinstance(x, (int, float, np.floating)) else x

    rows = []
    for _, r in dedup.iterrows():
        agree = r["slopes_agree"]
        rows.append({
            "loc": r["location_id"], "lay": r["layer_id"], "scn": r["scenario"],
            "dec": int(r["decade"]),
            "v": clean(r["value"]), "lo": clean(r["lower_ci"]), "hi": clean(r["upper_ci"]),
            "p": clean(r["percentile"]),
            "ols": clean(r["ols_slope"]), "sen": clean(r["sen_slope"]),
            "ag": None if pd.isna(agree) else bool(agree),
            "nm": None if pd.isna(r["n_members"]) else int(r["n_members"]),
            "nmo": None if pd.isna(r["n_models"]) else int(r["n_models"]),
            "st": r["data_status"],
        })

    def s(x):
        return "" if pd.isna(x) else str(x)

    return {
        "locations": [
            {"id": r["location_id"], "name": r["name"], "lat": float(r["lat"]),
             "lon": float(r["lon"]), "country": s(r.get("country")),
             "state": s(r.get("state")), "city": s(r.get("city"))}
            for _, r in locations.iterrows()
        ],
        "assets": [
            {"id": r["asset_id"], "loc": r["location_id"], "type": s(r["asset_type"]),
             "sub": s(r.get("sub_asset_unit")),
             "layers": [x for x in s(r["layer_ids"]).split(";") if x]}
            for _, r in assets.iterrows()
        ],
        "layers": {
            r["layer_id"]: {
                "hazard": s(r["hazard"]), "measure": s(r["hazard_measure"]),
                "units": s(r.get("units")), "slope_units": s(r.get("slope_units")),
                "direction": s(r.get("percentile_direction")) or "higher_is_worse",
                "statistic": s(r.get("decadal_statistic")),
                "recommended": s(r.get("recommended_slope")) or "sen_slope",
                "qa": s(r.get("qa_reviewed_on")),
                "note": s(r.get("delivery_note")),
                "scenarios": [x for x in s(r.get("scenarios")).split(";") if x],
                "members": s(r.get("n_members_by_scenario")),
            }
            for _, r in layers.iterrows()
        },
        "values": rows,
        "tiers": {s: tier_of(s) for s in sorted(values["scenario"].unique())},
        # Climate Score stays at ASSET grain here. Rolling it up to a location or the
        # portfolio is the dashboard's job and depends on the live asset-type filter, so
        # pre-aggregating would freeze a filter choice into the payload.
        # v2 (2026-08-20): n counts FAMILIES present, ne the asset type's weighted
        # families, and h the family list -- rollups compare rows on identical h, since
        # full coverage is structurally unreachable (permafrost exists at few sites;
        # cyclone and sea level publish no high tier).
        "scores": [
            {"a": r["asset_id"], "t": r["scenario_tier"], "d": int(r["decade"]),
             "s": float(r["climate_score"]), "n": int(r["n_hazards"]),
             "ne": int(r["n_hazards_expected"]), "h": str(r.get("hazards") or "")}
            for _, r in scores.iterrows()
        ] if len(scores) else [],
    }


def esc(text) -> str:
    """HTML-escape a value before it is interpolated into markup.

    Customer names, location names, asset types and layer notes are untrusted text that
    reaches the page. Anything not passed through here can inject markup.
    """
    return html_escape(str(text), quote=True)


def render(delivery: Path, payload: dict, manifest: dict, warnings: list) -> Path:
    customer = esc(manifest.get("customer", delivery.parent.name))
    generated = manifest.get("generated_utc", "")
    try:
        generated = datetime.fromisoformat(generated).strftime("%Y-%m-%d %H:%M UTC")
    except ValueError:
        pass
    counts_line = (f'{len(payload["locations"])} locations · '
                   f'{len(payload["assets"])} assets · '
                   f'{len(payload["layers"])} layers')
    unreviewed = [lid for lid, m in payload["layers"].items() if m["qa"] == "NOT CONFIRMED"]

    banner = ""
    if unreviewed:
        banner = (
            '<div class="banner">'
            f'<strong>Not yet quality-reviewed: {esc(", ".join(unreviewed))}.</strong> '
            'Treat these numbers as provisional.'
            '</div>'
        )
    # Build warnings are OPERATOR-facing and are printed to the terminal at build time;
    # they no longer render into the page (user decision 2026-08-20 -- the banner read as
    # a defect list to customers while describing expected portfolio structure).

    tokens = css_tokens()
    # SAFE JSON EMBEDDING. json.dumps escapes quotes and backslashes but NOT `<`, `>` or
    # `&`, so a customer location literally named `Depot </script><script>...` terminates
    # this script element early -- breaking the page and executing whatever follows.
    # Verified 2026-08-13 by generating a delivery with that name. The \uXXXX forms are
    # valid JSON and decode to the same characters, so the data is unchanged.
    data_json = (json.dumps(payload, separators=(",", ":"))
                 .replace("<", "\\u003c").replace(">", "\\u003e")
                 .replace("&", "\\u0026").replace("\u2028", "\\u2028")
                 .replace("\u2029", "\\u2029"))

    # BUILD STAMP. A browser serving a cached dashboard is indistinguishable from a
    # correctly regenerated one by eye -- that ambiguity has already sent us chasing two
    # phantom bugs. The stamp is a hash of the payload AND the page logic, so any real
    # change moves it, and `--check-stamp` can confirm what is on disk from the terminal.
    build = hashlib.sha256((data_json + JS + tokens).encode()).hexdigest()[:8]
    built_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{customer} — Climate hazard QA dashboard</title>
<script src="{PLOTLY_CDN}"></script>
<style>
{tokens}
* {{ box-sizing: border-box; }}
body {{ margin:0; background: var(--plane); color: var(--text-primary);
       font-family: {FONT_STACK}; font-size: 14px; line-height: 1.5; }}
header {{ padding: 20px 24px 12px; position: relative; }}
header h1 {{ margin:0 0 4px; font-size: 20px; font-weight: 600; }}
header .sub {{ color: var(--text-secondary); font-size: 13px; }}
.faq {{ position:absolute; top:20px; right:24px; }}
.faq summary {{ list-style:none; cursor:pointer; font-size:13px; padding:6px 14px;
    border:1px solid var(--axis); border-radius:6px; background:var(--surface);
    color:var(--text-secondary); user-select:none; }}
.faq summary::-webkit-details-marker {{ display:none; }}
.faq[open] summary {{ background: var(--series-1); color:#fff; }}
.faq-panel {{ position:absolute; right:0; top:40px; width:min(480px, 90vw);
    max-height:70vh; overflow:auto; background:var(--surface);
    border:1px solid var(--border); border-radius:8px; padding:14px 16px;
    box-shadow:0 8px 30px rgba(0,0,0,0.18); z-index:40; font-size:13px; }}
.faq-panel h3 {{ margin:12px 0 4px; font-size:13px; color: var(--text-primary); }}
.faq-panel h3:first-child {{ margin-top:0; }}
.faq-panel p {{ margin:0 0 6px; color: var(--text-secondary); }}
.banner {{ margin: 8px 24px; padding: 10px 14px; border-radius: 6px;
          background: color-mix(in srgb, #fab219 18%, var(--surface));
          border: 1px solid var(--border); font-size: 13px; }}
.banner.warn {{ background: color-mix(in srgb, #ec835a 18%, var(--surface)); }}
/* Pinned so the controls stay reachable while reading the charts further down. Full-bleed
   rather than inset, so nothing scrolls through the gap at its sides. */
.filters {{ display:flex; flex-wrap:wrap; gap:14px; align-items:flex-end;
           padding: 12px 24px; background: var(--surface);
           border-bottom:1px solid var(--axis);
           position: sticky; top: 0; z-index: 30;
           box-shadow: 0 1px 6px rgba(0,0,0,0.06); }}
.filters label {{ display:block; font-size:11px; text-transform:uppercase;
                 letter-spacing:.04em; color: var(--text-muted); margin-bottom:4px; }}
select, input[type=search] {{ font:inherit; font-size:13px; padding:6px 8px;
    background: var(--surface); color: var(--text-primary);
    border:1px solid var(--axis); border-radius:5px; min-width:150px; }}
.grid {{ display:grid; grid-template-columns: 1fr 1fr; gap:16px; padding:0 24px 24px; }}
.card {{ background: var(--surface); border:1px solid var(--border); border-radius:8px;
        padding:16px; min-width:0; }}
.card.full {{ grid-column: 1 / -1; }}
.card h2 {{ margin:0 0 2px; font-size:15px; font-weight:600; }}
.card .note {{ color: var(--text-secondary); font-size:12px; margin-bottom:10px; }}
.toggle {{ display:inline-flex; gap:2px; margin-bottom:10px; border:1px solid var(--axis);
          border-radius:6px; overflow:hidden; flex-wrap:wrap; }}
.toggle button {{ font:inherit; font-size:12px; padding:6px 12px; border:0;
    background: var(--surface); color: var(--text-secondary); cursor:pointer; }}
.toggle button[aria-pressed=true] {{ background: var(--series-1); color:#fff; }}
button.showbtn {{ font:inherit; font-size:12px; padding:6px 12px; margin-right:14px;
    border:1px solid var(--axis); border-radius:6px; background:var(--surface);
    color:var(--text-secondary); cursor:pointer; }}
button.showbtn[aria-pressed=true] {{ background: var(--series-1); color:#fff; }}
table {{ border-collapse:collapse; width:100%; font-size:12px;
        font-variant-numeric: tabular-nums; }}
th, td {{ text-align:left; padding:5px 8px; border-bottom:1px solid var(--grid);
         white-space:nowrap; }}
th {{ position:sticky; top:0; background: var(--surface); cursor:pointer;
     color: var(--text-secondary); font-weight:600; z-index:1; }}
/* Per-column dropdown filters sit in a second header row, below the sort row. */
tr.colf th {{ top:29px; cursor:default; padding:2px 4px; }}
tr.colf select {{ min-width:0; width:100%; font-size:11px; padding:3px 4px; }}
/* Portfolio overview: the story is one number, so it is a stat tile, not a chart. */
.tiles {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(190px,1fr));
         gap:14px; }}
.tile .k {{ font-size:11px; text-transform:uppercase; letter-spacing:.04em;
           color: var(--text-muted); margin-bottom:6px; }}
.tile .v {{ font-size:34px; line-height:1.1; font-weight:600;
           color: var(--text-primary); }}   /* proportional figures, not tabular-nums */
.tile .s {{ font-size:12px; color: var(--text-secondary); margin-top:4px; }}
td.num {{ text-align:right; }}
.scroll {{ max-height:460px; overflow:auto; }}
.chip {{ display:inline-block; padding:1px 7px; border-radius:9px; font-size:11px;
        border:1px solid var(--border); color: var(--text-secondary); }}
footer {{ padding: 8px 24px 32px; color: var(--text-muted); font-size:12px; }}
@media (max-width: 1100px) {{ .grid {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<header>
  <h1>{customer} — climate hazard QA dashboard</h1>
  <div class="sub">Generated {generated} &nbsp;·&nbsp; {counts_line}</div>
  <details class="faq">
    <summary>FAQ</summary>
    <div class="faq-panel">
      <h3>What is a risk percentile?</h3>
      <p>How a site compares with the rest of the world on one hazard, from 1 (lowest
         risk) to 100 (highest). Every hazard is already oriented so 100 is always worst,
         including ones where "more" is good.</p>
      <h3>What is an emissions pathway?</h3>
      <p>The low / medium / high scenario tiers. They group the two climate-scenario
         families (RCP and SSP) from different model generations; the tiers are only
         approximately comparable, which any narrative built on them should say.</p>
      <h3>What is the Climate Score?</h3>
      <p>The average risk percentile across the hazard families relevant to an asset type
         (each family counts 0 or 1). A hazard with no data at a site is left out of that
         site's average, and the hover shows how many hazards were included of how many
         are weighted.</p>
      <h3>Why do some hazards look identical in every decade and pathway?</h3>
      <p>Tornado, landslide and hail come from measured historical records, not climate
         models. They enter the score as the same "current conditions" value everywhere,
         which also slightly damps how much the score moves between pathways.</p>
      <h3>What are OLS and Sen?</h3>
      <p>Two methods of estimating a trend per decade. Each hazard has a recommended one
         (in layers.csv); where the two disagree, treat the trend as uncertain.</p>
      <h3>Why is a panel blank, or flat at zero?</h3>
      <p>Blank with "no coverage" = the observing record does not cover that site
         (tornado outside the United States) — unknown, not safe. A flat dotted zero
         marked "absent" = the hazard physically does not occur there (permafrost outside
         permafrost regions).</p>
      <h3>Why do time comparisons mention "the same assets"?</h3>
      <p>Any change over time or gap between pathway lines is computed only over assets
         with identical hazard coverage in everything compared — otherwise a change in
         coverage would read as a change in risk.</p>
      <h3>What do the map colours mean?</h3>
      <p>Percentile views: green (low) → yellow → orange → red (high). The Climate Score
         map: blue (low) → white → red (high). Trend views: red = worsening, blue =
         improving, capped so a few extreme sites don't wash out the rest.</p>
      <h3>Why does the table stop at 800 rows?</h3>
      <p>To keep the page fast. The truncation is stated in the table, and the delivered
         CSVs always carry every row.</p>
      <h3>Where is the full documentation?</h3>
      <p>README.md, delivered alongside the CSVs, including how the score is built and
         the delivery's caveats.</p>
    </div>
  </details>
</header>
{banner}

<div class="filters">
  <div><label for="f-tier">Forcing tier</label><select id="f-tier"></select></div>
  <div><label for="f-hazard">Hazard</label><select id="f-hazard"></select></div>
  <div><label for="f-asset">Asset type</label><select id="f-asset"></select></div>
  <div><label for="f-decade">Decade</label><select id="f-decade"></select></div>
  <div style="color:var(--text-muted);font-size:12px;max-width:330px">
    These filters apply to every panel, except where a panel says otherwise.
  </div>
</div>

<div class="grid">
  <div class="card full">
    <h2>Portfolio overview — Climate Score</h2>
    <div class="note" id="tiles-note"></div>
    <div class="tiles" id="tiles"></div>
  </div>

  <div class="card full">
    <h2>Site map</h2>
    <div class="note" id="map-note"></div>
    <div class="toggle" id="metric-toggle"></div>
    <div id="map" style="height:520px"></div>
  </div>

  <div class="card"><h2>Average risk percentile by hazard</h2>
    <div class="note">All hazards, selected decade and pathway.</div>
    <div id="bar-hazard" style="height:300px"></div></div>

  <div class="card"><h2>Average risk percentile by asset type</h2>
    <div class="note">Average across each asset type's relevant hazards.</div>
    <div id="bar-asset" style="height:300px"></div></div>

  <div class="card"><h2>Portfolio risk distribution</h2>
    <div class="note" id="band-note"></div>
    <div id="bar-band" style="height:300px"></div></div>

  <div class="card"><h2>Climate Score over time</h2>
    <div style="margin-bottom:10px">
      <label for="score-sel" style="font-size:11px;text-transform:uppercase;
          letter-spacing:.04em;color:var(--text-muted)">Series</label>
      <select id="score-sel"></select>
    </div>
    <div class="note" id="score-note"></div>
    <div id="score-series" style="height:300px"></div></div>

  <div class="card full">
    <h2>Time series by location</h2>
    <div class="note">One location at a time: one chart per hazard, plus the Climate
        Score. Lines are pathways; all decades shown. Hidden by default.</div>
    <div style="margin-bottom:10px">
      <button id="series-show" class="showbtn" aria-pressed="false">Show plots</button>
      <label for="f-loc" style="font-size:11px;text-transform:uppercase;
          letter-spacing:.04em;color:var(--text-muted)">Location</label>
      <select id="f-loc"></select>
      <span class="toggle" id="series-toggle" style="margin-left:10px"></span>
    </div>
    <div id="series" style="height:420px; display:none"></div>
  </div>

  <div class="card full">
    <h2>Values</h2>
    <div class="note">Filter with the dropdowns and search (the page filters above don't
        apply here); click a heading to sort.</div>
    <input type="search" id="f-search" placeholder="Filter rows…" style="margin-bottom:10px">
    <div class="scroll"><table id="table"></table></div>
  </div>
</div>

<footer><span id="foot"></span>
  <span style="color:var(--text-muted)">· build
  <code style="font-variant-numeric:tabular-nums">{build}</code></span></footer>

<script>
const DATA = {data_json};
/* Blue = raw magnitude ("how much"). Red = risk ("how bad"), where 100 is worst by
   construction. Dark-mode variants are truncated so the extreme step stays visible. */
const SEQ_VALUE = {json.dumps(plotly_colorscale(SEQUENTIAL_BLUE))};
const SEQ_RISK = {json.dumps(plotly_colorscale(SEQUENTIAL_RED))};
const SEQ_RISK_DARK = {json.dumps(plotly_colorscale(SEQUENTIAL_RED_DARK))};
const DIV = {json.dumps(diverging_colorscale('light'))};
const DIV_DARK = {json.dumps(diverging_colorscale('dark'))};
const ORD = {json.dumps(ORDINAL_RISK)};
const ORD_DARK = {json.dumps(ORDINAL_RISK_DARK)};
const BANDS = {json.dumps(RISK_BAND_ORDER)};
const TIER_ORDER = {json.dumps(TIER_ORDER)};
const TIER_LABELS = {json.dumps(TIER_LABELS)};
const TIER_COLOR = {json.dumps(TIER_COLOR_LIGHT)};
const TIER_COLOR_DARK = {json.dumps(TIER_COLOR_DARK)};
/* Secondary encoding, applied unconditionally: blue and red are near-identical in
   luminance and yellow/red converge under deuteranopia, so hue alone cannot carry tier
   identity. See viz_common.TIER_COLOR_LIGHT for the measurements. */
const TIER_SYMBOL = {json.dumps(TIER_SYMBOL)};
const TIER_DASH = {json.dumps(TIER_DASH)};
{JS}
</script>
</body>
</html>
"""
    out = delivery / "dashboard.html"
    out.write_text(html)
    return out, build


JS = r"""
const $ = id => document.getElementById(id);
// Untrusted text -> markup. Location names, asset types and statuses come from a customer
// CSV and several of these are built with innerHTML, so every interpolated DATA value goes
// through this. (The payload itself is \u-escaped server-side so it cannot break out of
// the script element; this covers what we then write into the DOM.)
const esc = v => String(v ?? "").replace(/[&<>"\']/g,
  c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));

// Any uncaught page error becomes a visible banner instead of a silently dead control --
// added 2026-08-21 while chasing a table-filter fault that never surfaced an error we
// could be told about. If this banner ever appears, report its text verbatim.
function showError(msg) {
  let el = document.getElementById("err-banner");
  if (!el) {
    el = document.createElement("div");
    el.id = "err-banner";
    el.style.cssText = "position:fixed;bottom:0;left:0;right:0;z-index:99;" +
      "background:#7f1d1d;color:#fff;font:12px/1.5 ui-monospace,monospace;padding:8px 16px;";
    document.body.appendChild(el);
  }
  el.textContent = "Page error — please report this text: " + msg;
}
window.addEventListener("error",
  e => showError((e.message || "?") + " @ line " + (e.lineno || "?")));
window.addEventListener("unhandledrejection", e => showError(String(e.reason)));
const locById = Object.fromEntries(DATA.locations.map(l => [l.id, l]));
const fmt = (x, d=3) => x === null || x === undefined ? "—" :
    (Math.abs(x) >= 0.01 && Math.abs(x) < 1e5 ? Number(x).toFixed(d) : Number(x).toExponential(2));

// Which layers does a location actually carry, for a given asset-type filter?
function locLayers(assetType) {
  const m = {};
  for (const a of DATA.assets) {
    if (assetType !== "__all" && a.type !== assetType) continue;
    (m[a.loc] = m[a.loc] || new Set());
    a.layers.forEach(l => m[a.loc].add(l));
  }
  return m;
}

function ink(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
// Theme-dependent scales. The CSS tokens swap themselves; Plotly scales are JS values, so
// they have to be chosen at render time against the same signal.
function isDark() {
  const stamp = document.documentElement.getAttribute("data-theme");
  if (stamp) return stamp === "dark";
  return matchMedia("(prefers-color-scheme: dark)").matches;
}
const riskScale = () => isDark() ? SEQ_RISK_DARK : SEQ_RISK;
// Percentile ramp: green -> yellow -> orange -> red, continuous over 1-100 (user call,
// 2026-08-21). Used wherever a raw risk PERCENTILE is colour-coded (percentile map view
// and the percentile bar charts); the Climate Score map keeps its blue-white-red ramp.
const PCT_RAMP = [[0, "#1a9850"], [0.25, "#91cf60"], [0.5, "#fee08b"],
                  [0.75, "#f46d43"], [1, "#d73027"]];
const divScale  = () => isDark() ? DIV_DARK : DIV;
const ordScale  = () => isDark() ? ORD_DARK : ORD;
const tierColor = t => (isDark() ? TIER_COLOR_DARK : TIER_COLOR)[t]
                       || (isDark() ? TIER_COLOR_DARK : TIER_COLOR).medium;
// Every tier line carries hue + symbol + dash, so identity survives CVD and print.
function tierStyle(t) {
  return {
    line: {color: tierColor(t), width: 2, dash: TIER_DASH[t] || "solid"},
    marker: {size: 7, symbol: TIER_SYMBOL[t] || "circle",
             line: {width: 2, color: ink("--surface")}},
  };
}
function baseLayout(extra) {
  return Object.assign({
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
    font: {color: ink("--text-secondary"), size: 12,
           family: 'system-ui, -apple-system, "Segoe UI", sans-serif'},
    margin: {l: 60, r: 20, t: 24, b: 44},
    xaxis: {gridcolor: ink("--grid"), linecolor: ink("--axis"), zeroline: false},
    yaxis: {gridcolor: ink("--grid"), linecolor: ink("--axis"), zeroline: false},
    showlegend: false,
  }, extra || {});
}
const CFG = {responsive: true, displaylogo: false,
             modeBarButtonsToRemove: ["select2d","lasso2d","autoScale2d"]};

// Basemap follows the page theme (a fixed dark basemap was tried 2026-08-21 and reverted
// the same day, both user calls); the marker ring is fixed DARK GREY so points read as
// outlined against the light land.
const geoLayout = () => ({showland: true, landcolor: ink("--grid"), showcountries: true,
  countrycolor: ink("--surface"), showocean: true, oceancolor: ink("--plane"),
  showlakes: true, lakecolor: ink("--plane"), coastlinecolor: ink("--axis"),
  coastlinewidth: 0.5, projection: {type: "natural earth"}, bgcolor: "rgba(0,0,0,0)"});
const MARKER_RING = {width: 1.5, color: "#444444"};

// ---- filter state -------------------------------------------------------------------
// Climate Score leads: it is the portfolio-level answer, and every other metric is a
// component of it. The hazard filter does not apply while it is selected.
const state = {tier: "__all", hazard: null, asset: "__all", decade: null,
               metric: "climate_score", loc: null, seriesMetric: "percentile",
               scoreMetric: "__score", search: "", sortCol: null, sortDir: 1,
               colFilters: {}, seriesVisible: false};

// --- Climate Score helpers ------------------------------------------------------------
// v2 model (2026-08-20): a score row is the mean over the asset type's WEIGHTED hazard
// families PRESENT at the site, and full coverage is structurally unreachable (permafrost
// exists at few sites; cyclone and sea level publish no high tier). So rows are shown as
// delivered -- renormalized, with n of ne families disclosed -- and never filtered away
// for being "incomplete". What IS guarded is comparison: any rollup across decades first
// restricts to assets whose family set (r.h) is IDENTICAL in every compared cell, so a
// composition change can never read as a risk change.
const assetById = Object.fromEntries(DATA.assets.map(a => [a.id, a]));
function scoreRows({tier, decade} = {}) {
  return DATA.scores.filter(r => {
    const a = assetById[r.a];
    if (!a) return false;
    if (state.asset !== "__all" && a.type !== state.asset) return false;
    if (tier && tier !== "__all" && r.t !== tier) return false;
    if (decade !== undefined && decade !== null && r.d !== decade) return false;
    return true;
  });
}
const mean = xs => xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : null;

// BALANCED PANEL -- the single most important guard on this page.
//
// Layer scenario coverage is NOT uniform across tiers: `cyclone` and `sealevel-2b`
// publish no high-tier scenario, so the high-tier line structurally lacks those families
// for every asset. Averaging each cell over whatever happens to be there makes the lines
// describe different portfolios -- measured on an early delivery, the high-tier 2020s
// mean read 39.9 against 42.1 for low and medium, purely composition.
//
// The v2 guard: an asset joins a comparison panel only if, WITHIN EACH tier compared, its
// family set (r.h) is identical across every decade compared. Family sets may differ
// BETWEEN tiers -- that gap is structural and is stated in the chart note (tierGaps())
// rather than blanking the chart, because no asset could ever satisfy cross-tier equality.
function balancedAssets(tiers, decades, locId) {
  const cand = DATA.assets
    .filter(a => (state.asset === "__all" || a.type === state.asset)
                 && (!locId || a.loc === locId))
    .map(a => a.id);
  const cell = {};
  DATA.scores.forEach(r => { (cell[r.a] = cell[r.a] || {})[r.t + "|" + r.d] = r.h; });
  return cand.filter(id => {
    const m = cell[id];
    if (!m) return false;
    return tiers.every(t => {
      let fam = null;
      for (const d of decades) {
        const h = m[t + "|" + d];
        if (h === undefined) return false;
        if (fam === null) fam = h;
        else if (h !== fam) return false;
      }
      return true;
    });
  });
}

// Which layers are missing a tier entirely? That is a portfolio-level fact, not a bug.
function tierGaps() {
  const gaps = [];
  for (const [lid, m] of Object.entries(DATA.layers)) {
    const present = new Set(m.scenarios.map(s => DATA.tiers[s] || "medium"));
    const missing = TIER_ORDER.filter(t => !present.has(t));
    if (missing.length) gaps.push(`${lid} has no ${missing.join("/")} scenario`);
  }
  return gaps;
}

function tierOf(scn) { return DATA.tiers[scn] || "medium"; }

function initControls() {
  const layerIds = Object.keys(DATA.layers);
  const decades = [...new Set(DATA.values.map(v => v.dec))].sort((a,b)=>a-b);
  const assetTypes = [...new Set(DATA.assets.map(a => a.type))].sort();
  const tiers = [...new Set(DATA.values.map(v => tierOf(v.scn)))]
        .sort((a,b) => TIER_ORDER.indexOf(a) - TIER_ORDER.indexOf(b));

  fill($("f-tier"), [["__all","All tiers"]].concat(tiers.map(t => [t, TIER_LABELS[t]])));
  // Hazard dropdowns sort alphabetically by label, and the map initialises on the first
  // option (user call, 2026-08-21).
  const hazardOpts = layerIds
      .map(l => [l, DATA.layers[l].hazard + " (" + l + ")"])
      .sort((a, b) => a[1].localeCompare(b[1]));
  fill($("f-hazard"), hazardOpts);
  fill($("f-loc"), DATA.locations.map(l => [l.id, l.name]));

  state.hazard = hazardOpts[0][0];
  state.decade = decades.includes(2050) ? 2050 : decades[decades.length-1];
  state.loc = DATA.locations[0].id;
  $("f-hazard").value = state.hazard;
  $("f-loc").value = state.loc;
  syncFilterOptions();

  $("f-tier").onchange   = e => { state.tier = e.target.value; renderAll(); };
  $("f-hazard").onchange = e => { state.hazard = e.target.value;
                                  syncFilterOptions(); renderAll(); };
  $("f-asset").onchange  = e => { state.asset = e.target.value; renderAll(); };
  $("f-decade").onchange = e => { state.decade = +e.target.value; renderAll(); };
  $("f-loc").onchange    = e => { state.loc = e.target.value; renderSeries(); };
  // The panel grid draws one panel per standard-set layer and is the heaviest thing on
  // the page, so it renders only on demand (user call, 2026-08-21).
  $("series-show").onclick = () => {
    state.seriesVisible = !state.seriesVisible;
    const b = $("series-show");
    b.setAttribute("aria-pressed", String(state.seriesVisible));
    b.textContent = state.seriesVisible ? "Hide plots" : "Show plots";
    renderSeries();
  };
  $("f-search").oninput  = e => { state.search = e.target.value.toLowerCase();
                                  updateFilterOptions(-1); renderTable(); };

  buildToggle("metric-toggle",
    [["climate_score","Climate Score"],["value","Value"],["percentile","Percentile"],
     ["ols_slope","OLS trend"],["sen_slope","Sen trend"]],
    () => state.metric, v => { state.metric = v; renderMap(); });
  // Climate Score is a PANEL in this grid, not a mode -- it sits beside the hazard
  // percentiles it is the mean of. The toggle therefore only switches the hazard panels
  // between their risk percentile and their native-unit value.
  buildToggle("series-toggle", [["percentile","Percentile"],["value","Value"]],
    () => state.seriesMetric, v => { state.seriesMetric = v; renderSeries(); });

  // Climate Score, or the per-hazard percentile it is the mean of. Same axis meaning
  // (1-100 risk percentile) either way, so the two are directly comparable. A DROPDOWN,
  // not a toggle: with the v2 standard set there are 30+ layers and a button row cannot
  // list them (user call, 2026-08-20).
  fill($("score-sel"),
    [["__score","Climate Score (all hazards)"]].concat(hazardOpts));
  $("score-sel").value = state.scoreMetric;
  $("score-sel").onchange = e => { state.scoreMetric = e.target.value; renderScoreSeries(); };

  initTable();   // builds the table header + filter row ONCE; renderTable fills tbody
}

function fill(sel, pairs) {
  sel.innerHTML = pairs.map(([v,l]) =>
    `<option value="${esc(v)}">${esc(l)}</option>`).join("");
}

// CASCADING FILTERS (user call, 2026-08-21): the hazard selection repopulates the Decade
// and Asset-type options with only the choices that actually show data for that layer, so
// a reader is never offered a decade an observational layer does not have, or an asset
// type with no data anywhere on that layer. Invalid selections snap to a valid one.
function syncFilterOptions() {
  const lid = state.hazard;
  const has = DATA.values.filter(v => v.lay === lid && (v.p !== null || v.v !== null));
  const decs = [...new Set(has.map(v => v.dec))].sort((a,b) => a-b);
  const single = decs.length === 1 && decs[0] === 2020;
  fill($("f-decade"), decs.map(d => [d, single ? "Current" : d + "s"]));
  if (!decs.includes(state.decade))
    state.decade = decs.includes(2050) ? 2050 : decs[decs.length - 1];
  $("f-decade").value = state.decade;

  const locsWith = new Set(has.map(v => v.loc));
  const types = [...new Set(DATA.assets.filter(a => locsWith.has(a.loc)).map(a => a.type))]
      .sort();
  fill($("f-asset"), [["__all","All asset types"]].concat(types.map(t => [t, t])));
  if (state.asset !== "__all" && !types.includes(state.asset)) state.asset = "__all";
  $("f-asset").value = state.asset;
}
function buildToggle(id, opts, get, set) {
  const el = $(id);
  el.innerHTML = opts.map(([v,l]) =>
    `<button data-v="${esc(v)}" aria-pressed="${get()===v}">${esc(l)}</button>`).join("");
  el.querySelectorAll("button").forEach(b => b.onclick = () => {
    set(b.dataset.v);
    el.querySelectorAll("button").forEach(x =>
      x.setAttribute("aria-pressed", x.dataset.v === get()));
  });
}

// ---- row selection ------------------------------------------------------------------
// Scenario is filtered by TIER, so a mixed-round portfolio (cyclone is RCP-only, wildfire
// SSP-only) still answers "show me the high-forcing view" with each layer's own code.
function rows({hazard, decade, allDecades, allHazards} = {}) {
  const ll = locLayers(state.asset);
  return DATA.values.filter(v => {
    if (!allHazards && v.lay !== (hazard || state.hazard)) return false;
    if (!allDecades && v.dec !== (decade || state.decade)) return false;
    if (state.tier !== "__all" && tierOf(v.scn) !== state.tier) return false;
    const s = ll[v.loc];
    return s && s.has(v.lay);
  });
}

// ---- map ----------------------------------------------------------------------------
function renderMap() {
  if (state.metric === "climate_score") return renderScoreMap();
  const lay = DATA.layers[state.hazard];
  const isSlope = state.metric === "ols_slope" || state.metric === "sen_slope";
  const key = {value:"v", percentile:"p", ols_slope:"ols", sen_slope:"sen"}[state.metric];

  // One point per (location, scenario) is ambiguous on a map, so collapse to the tier's
  // representative scenario for this layer -- the worst-case reading when several match.
  const byLoc = {};
  for (const r of rows()) {
    const cur = byLoc[r.loc];
    if (!cur || (r[key] !== null && (cur[key] === null || r[key] > cur[key]))) byLoc[r.loc] = r;
  }
  const pts = Object.values(byLoc).filter(r => r[key] !== null);

  let cs, cmin, cmax, note, rev = false;
  if (isSlope) {
    const [L, clamped, rule] = pointLimit(pts.map(r => r[key]));
    cs = divScale(); cmin = -L; cmax = L;
    // Red always means WORSE. On a higher_is_better layer (stored carbon, productivity)
    // an increase is an improvement, so the scale reverses.
    rev = lay.direction === "higher_is_better";
    note = `Trend per decade. Red = worsening, blue = improving; colour capped at ` +
           `±${fmt(L)}` + (clamped ? ` (${clamped} site(s) at the cap)` : "") + ".";
  } else if (state.metric === "percentile") {
    cs = PCT_RAMP; cmin = 1; cmax = 100;
    note = "Risk percentile, 1 (lowest) – 100 (highest). Green = low, red = high.";
  } else {
    cs = SEQ_VALUE; cmin = null; cmax = null;
    note = `${lay.measure} [${lay.units || "—"}] · darker = larger`;
  }

  const nonRobust = pts.filter(r => r.ag === false).length;
  if (isSlope && nonRobust)
    note += ` ${nonRobust} of ${pts.length} sites: trend uncertain.`;
  $("map-note").textContent = note;

  const trace = {
    type: "scattergeo", mode: "markers",
    lat: pts.map(r => locById[r.loc].lat), lon: pts.map(r => locById[r.loc].lon),
    text: pts.map(r => {
      const l = locById[r.loc];
      return `<b>${esc(l.name)}</b><br>${esc(lay.hazard)} · ${esc(r.scn)} · ${r.dec}s` +
        `<br>value ${fmt(r.v)} ${esc(lay.units)}` +
        `<br>percentile ${r.p === null ? "—" : Math.round(r.p)}` +
        `<br>ols ${fmt(r.ols)} · sen ${fmt(r.sen)}` +
        `<br>slopes agree: ${r.ag === null ? "n/a (inactive)" : r.ag}` +
        `<br>members ${r.nm ?? "—"} · models ${r.nmo ?? "—"} · ${esc(r.st)}`;
    }),
    hovertemplate: "%{text}<extra></extra>",
    marker: {
      size: 13, color: pts.map(r => r[key]), colorscale: cs, reversescale: rev,
      cmin: cmin === null ? undefined : cmin, cmax: cmax === null ? undefined : cmax,
      line: MARKER_RING,
      colorbar: {thickness: 12, len: 0.7, outlinewidth: 0,
                 tickfont: {color: ink("--text-secondary")}},
    },
  };
  Plotly.react("map", [trace], baseLayout({
    margin: {l: 0, r: 0, t: 0, b: 0},
    geo: geoLayout(),
  }), CFG);
}

// Climate Score map: one marker per LOCATION, scored as the mean over that location's
// assets (after the asset-type filter). The hazard filter deliberately does not apply --
// the whole point of the score is that it spans hazards.
function renderScoreMap() {
  const rowsIn = scoreRows({tier: state.tier, decade: state.decade});
  const byLoc = {};
  for (const r of rowsIn) {
    const loc = assetById[r.a].loc;
    (byLoc[loc] = byLoc[loc] || []).push(r);
  }
  const pts = Object.entries(byLoc).map(([loc, rs]) => ({
    loc, s: mean(rs.map(x => x.s)), n: rs.length,
    fmin: Math.min(...rs.map(x => x.n)), fmax: Math.max(...rs.map(x => x.n)),
    ne: Math.max(...rs.map(x => x.ne)),
    tiers: [...new Set(rs.map(x => x.t))],
  }));

  $("map-note").textContent =
    `Overall climate risk per location, ${state.decade}s` +
    (state.tier === "__all" ? ", all pathways" :
       `, ${TIER_LABELS[state.tier].toLowerCase()} pathway`) +
    ". Blue = low, white = middle, red = high.";

  const trace = {
    type: "scattergeo", mode: "markers",
    lat: pts.map(p => locById[p.loc].lat), lon: pts.map(p => locById[p.loc].lon),
    text: pts.map(p => {
      const l = locById[p.loc];
      const fam = p.fmin === p.fmax ? String(p.fmax) : `${p.fmin}–${p.fmax}`;
      return `<b>${esc(l.name)}</b><br>Climate Score ${p.s.toFixed(1)}` +
        `<br>families ${fam} of ${p.ne} weighted · ${p.n} asset-row(s)` +
        `<br>tier(s): ${esc(p.tiers.join(", "))}`;
    }),
    hovertemplate: "%{text}<extra></extra>",
    marker: {
      // Diverging blue -> white -> red over the 1-100 score (user call, 2026-08-21):
      // more colour variability than a one-hue ramp. Theme-dependent variant, since the
      // basemap follows the page theme again.
      size: 15, color: pts.map(p => p.s), colorscale: divScale(), cmin: 1, cmax: 100,
      line: MARKER_RING,
      colorbar: {thickness: 12, len: 0.7, outlinewidth: 0,
                 tickfont: {color: ink("--text-secondary")}},
    },
  };
  Plotly.react("map", [trace], baseLayout({
    margin: {l: 0, r: 0, t: 0, b: 0},
    geo: geoLayout(),
  }), CFG);
}

// Portfolio overview. One number is a stat tile, not a chart.
function renderTiles() {
  const decades = [...new Set(DATA.scores.map(r => r.d))].sort((a,b)=>a-b);
  const cur = scoreRows({tier: state.tier, decade: state.decade});
  const score = mean(cur.map(r => r.s));

  // Baseline is the layer contract's shared 2020s panel where present, else the earliest
  // decade with complete coverage. The delta compares TWO decades, so it needs a common
  // asset set across both -- otherwise it reports a composition change as a risk change.
  const baseDec = decades.includes(2020) ? 2020
        : decades.find(d => scoreRows({tier: state.tier, decade: d}).length);
  const tiersForDelta = state.tier === "__all" ? TIER_ORDER.filter(
        t => DATA.scores.some(r => r.t === t)) : [state.tier];
  const deltaPanel = new Set(balancedAssets(tiersForDelta, [baseDec, state.decade]));
  const pick = d => mean(DATA.scores
      .filter(r => r.d === d && deltaPanel.has(r.a)
              && (state.tier === "__all" || r.t === state.tier))
      .map(r => r.s));
  const now = pick(state.decade), base = pick(baseDec);
  const delta = (now !== null && base !== null) ? now - base : null;

  // Worst location at the selected slice.
  const byLoc = {};
  cur.forEach(r => (byLoc[assetById[r.a].loc] = byLoc[assetById[r.a].loc] || []).push(r.s));
  const worst = Object.entries(byLoc).map(([l, xs]) => [l, mean(xs)])
        .sort((a,b) => b[1]-a[1])[0];

  const nAssets = new Set(cur.map(r => r.a)).size;
  const nHaz = cur.length ? Math.max(...cur.map(r => r.n)) : 0;

  const gaps = tierGaps();
  $("tiles-note").textContent =
    `Average risk percentile across the hazards relevant to each asset type. Uses the ` +
    `pathway, asset-type and decade filters.` +
    (gaps.length ? ` Some hazards lack data for some pathways.` : "");
  $("tiles").innerHTML = [
    tile("Portfolio Climate Score", score === null ? "—" : score.toFixed(1),
         `${state.decade}s · ` + (state.tier === "__all" ? "all tiers"
            : TIER_LABELS[state.tier].toLowerCase()) +
         ` · ${new Set(cur.map(r => r.a)).size} assets`),
    tile("Change vs " + baseDec + "s", delta === null ? "—" :
         (delta >= 0 ? "+" : "") + delta.toFixed(1),
         delta === null ? "not enough overlap to compare"
            : `same ${deltaPanel.size} assets in both decades · ` +
              (delta > 0 ? "rising risk" : delta < 0 ? "falling risk" : "flat")),
    tile("Highest-risk location", worst ? worst[1].toFixed(1) : "—",
         worst ? locById[worst[0]].name : ""),
    tile("Coverage", nAssets + " assets",
         nHaz + " hazards (most included at any site) · " + cur.length + " score rows"),
  ].join("");
}
function tile(k, v, s) {
  return `<div class="tile"><div class="k">${esc(k)}</div><div class="v">${esc(v)}</div>` +
         `<div class="s">${esc(s)}</div></div>`;
}

// Climate Score over time: ALL tiers on one axis, so the scenario spread is the story.
// Scoped by asset type only -- the forcing-tier and decade filters deliberately do not
// apply, because varying those two is exactly what this chart shows.
function renderScoreSeries() {
  if (state.scoreMetric !== "__score") return renderHazardScoreSeries();
  const tiers = TIER_ORDER.filter(t => DATA.scores.some(r => r.t === t));
  const allDecades = [...new Set(DATA.scores.map(r => r.d))].sort((a,b)=>a-b);

  // Keep the widest decade span on which SOME asset is complete in every tier, then
  // restrict to the assets complete across that whole grid.
  const decades = allDecades.filter(d => balancedAssets(tiers, [d]).length > 0);
  const panel = new Set(balancedAssets(tiers, decades));
  const candidates = DATA.assets.filter(
      a => state.asset === "__all" || a.type === state.asset).length;

  const traces = tiers.map(t => {
    const xs = [], ys = [];
    decades.forEach(d => {
      const rs = DATA.scores.filter(r => r.t === t && r.d === d && panel.has(r.a));
      if (rs.length) { xs.push(d); ys.push(mean(rs.map(r => r.s))); }
    });
    return {
      type: "scatter", mode: "lines+markers", name: TIER_LABELS[t],
      x: xs, y: ys,
      ...tierStyle(t),
      hovertemplate: `${TIER_LABELS[t]}<br>%{x}s: %{y:.1f}<extra></extra>`,
    };
  });

  const droppedDec = allDecades.filter(d => !decades.includes(d));
  const gaps = tierGaps();

  // The decade filter asks "is SOME asset complete here?" and the panel then asks "which
  // assets are complete EVERYWHERE?" -- those can intersect to nothing (asset A complete
  // only in 2020, asset B only in 2030). The chart is then blank for a reason no reader
  // could infer, so say it outright rather than rendering empty axes.
  if (!panel.size) {
    Plotly.purge("score-series");
    $("score-note").textContent =
      `No group of assets keeps the same hazard coverage across all decades, so a ` +
      `consistent line can't be drawn. Try narrowing the asset-type filter.`;
    $("score-series").innerHTML = "";
    return;
  }

  $("score-note").textContent =
    `Portfolio average over time, one line per pathway ` +
    `(${panel.size} of ${candidates} assets on a consistent basis).`;

  Plotly.react("score-series", traces, baseLayout({
    showlegend: true,
    legend: {orientation: "h", y: -0.2, x: 0, font: {color: ink("--text-secondary")}},
    margin: {l: 52, r: 18, t: 8, b: 58},
    xaxis: {tickformat: "d", gridcolor: ink("--grid"), linecolor: ink("--axis"),
            zeroline: false},
    yaxis: {title: {text: "Climate Score"}, gridcolor: ink("--grid"),
            linecolor: ink("--axis"), zeroline: false},
    shapes: [{type: "line", xref: "x", yref: "paper", x0: state.decade, x1: state.decade,
              y0: 0, y1: 1, line: {color: ink("--axis"), width: 1}}],
  }), CFG);
}

// The same chart for ONE hazard: portfolio mean risk percentile over time, per tier. This
// is the component the Climate Score averages, on the same 1-100 axis, so switching
// between them shows which hazard drives the aggregate.
//
// No balanced panel is needed here: a single hazard has one scenario set, so every asset
// carrying it is present in exactly the tiers that hazard publishes. A tier the hazard does
// not publish simply produces no line -- and that absence is stated rather than left blank.
function renderHazardScoreSeries() {
  const lid = state.scoreMetric;
  const lay = DATA.layers[lid];
  const assetsWith = DATA.assets.filter(a =>
      (state.asset === "__all" || a.type === state.asset) && a.layers.includes(lid));
  const locs = new Set(assetsWith.map(a => a.loc));
  const decades = [...new Set(DATA.values.filter(v => v.lay === lid).map(v => v.dec))]
      .sort((a,b)=>a-b);
  const tiers = TIER_ORDER.filter(t => lay.scenarios.some(s => (DATA.tiers[s]||"medium") === t));

  const traces = tiers.map(t => {
    const xs = [], ys = [];
    decades.forEach(d => {
      const rs = DATA.values.filter(v => v.lay === lid && v.dec === d && locs.has(v.loc)
                                    && (DATA.tiers[v.scn]||"medium") === t && v.p !== null);
      if (rs.length) { xs.push(d); ys.push(mean(rs.map(r => r.p))); }
    });
    return {type: "scatter", mode: "lines+markers", name: TIER_LABELS[t],
            x: xs, y: ys, ...tierStyle(t),
            hovertemplate: `${TIER_LABELS[t]}<br>%{x}s: %{y:.1f}<extra></extra>`};
  });

  const missing = TIER_ORDER.filter(t => !tiers.includes(t));
  $("score-note").textContent =
    `${lay.hazard}: average percentile across ${locs.size} site(s), per pathway.` +
    (missing.length ? ` No ${missing.join("/")} line (no data for that pathway).` : "");

  Plotly.react("score-series", traces, baseLayout({
    showlegend: true,
    legend: {orientation: "h", y: -0.2, x: 0, font: {color: ink("--text-secondary")}},
    margin: {l: 52, r: 18, t: 8, b: 58},
    xaxis: {tickformat: "d", gridcolor: ink("--grid"), linecolor: ink("--axis"),
            zeroline: false},
    yaxis: {title: {text: "Risk percentile"}, range: [0, 105],
            gridcolor: ink("--grid"), linecolor: ink("--axis"), zeroline: false},
    shapes: [{type: "line", xref: "x", yref: "paper", x0: state.decade, x1: state.decade,
              y0: 0, y1: 1, line: {color: ink("--axis"), width: 1}}],
  }), CFG);
}

// Small portfolios have no heavy tail, so a percentile limit would clamp real sites for
// no readability gain. Mirrors viz_common.point_symmetric_limit.
function pointLimit(vals) {
  const v = vals.filter(x => x !== null && isFinite(x)).map(Math.abs).sort((a,b)=>a-b);
  if (!v.length) return [1, 0, "no data"];
  if (v.length < 40) return [v[v.length-1] || 1, 0, "max |value|, no clamping"];
  const L = v[Math.floor(0.95 * (v.length - 1))] || v[v.length-1];
  return [L, v.filter(x => x > L).length, "95th pct of |value|"];
}

// ---- bar summaries ------------------------------------------------------------------
function meanBy(rowsIn, keyFn) {
  const acc = {};
  for (const r of rowsIn) {
    if (r.p === null) continue;
    const k = keyFn(r);
    if (k === null) continue;
    (acc[k] = acc[k] || []).push(r.p);
  }
  return Object.entries(acc).map(([k, xs]) =>
    [k, xs.reduce((a,b)=>a+b,0)/xs.length, xs.length]).sort((a,b) => b[1]-a[1]);
}

function barTrace(pairs, hover) {
  return {
    type: "bar", orientation: "h",
    y: pairs.map(p => p[0]), x: pairs.map(p => p[1]),
    // Bars encode a risk PERCENTILE, so they take the same green-to-red ramp as the
    // percentile map (user call, 2026-08-21).
    marker: {color: pairs.map(p => p[1]), colorscale: PCT_RAMP, cmin: 0, cmax: 100},
    text: pairs.map(p => p[1].toFixed(0)), textposition: "outside",
    customdata: pairs.map(p => p[2]),
    hovertemplate: hover + "<extra></extra>",
  };
}

function renderBars() {
  const all = rows({allHazards: true});

  Plotly.react("bar-hazard",
    [barTrace(meanBy(all, r => DATA.layers[r.lay].hazard),
      "%{y}<br>mean percentile %{x:.1f}<br>%{customdata} site-scenario rows")],
    baseLayout({margin:{l:150,r:40,t:8,b:36},
                xaxis:{range:[0,105], title:{text:"Mean risk percentile"},
                       gridcolor:ink("--grid"), zeroline:false}}), CFG);

  // An asset's exposure = its own layers only, so join through the asset table.
  const byType = [];
  for (const a of DATA.assets) {
    if (state.asset !== "__all" && a.type !== state.asset) continue;
    for (const r of all) {
      if (r.loc === a.loc && a.layers.includes(r.lay) && r.p !== null)
        byType.push({t: a.type, p: r.p});
    }
  }
  const acc = {};
  byType.forEach(x => (acc[x.t] = acc[x.t] || []).push(x.p));
  const pairs = Object.entries(acc).map(([k,xs]) =>
    [k, xs.reduce((a,b)=>a+b,0)/xs.length, xs.length]).sort((a,b)=>b[1]-a[1]);
  Plotly.react("bar-asset",
    [barTrace(pairs, "%{y}<br>mean percentile %{x:.1f}<br>%{customdata} asset-hazard rows")],
    baseLayout({margin:{l:150,r:40,t:8,b:36},
                xaxis:{range:[0,105], title:{text:"Mean risk percentile"},
                       gridcolor:ink("--grid"), zeroline:false}}), CFG);

  // Portfolio: ordered risk bands -> the ORDINAL ramp (ordered categories, not nominal).
  const counts = Object.fromEntries(BANDS.map(b => [b, 0]));
  let n = 0;
  for (const x of byType) { const b = band(x.p); if (b) { counts[b]++; n++; } }
  $("band-note").textContent =
    `${n} asset-hazard combinations in the ${state.decade}s, grouped into five risk ` +
    `bands by their percentile.`;
  Plotly.react("bar-band", [{
    type: "bar", x: BANDS, y: BANDS.map(b => counts[b]),
    marker: {color: ordScale()},
    hovertemplate: "%{x}<br>%{y} combinations<extra></extra>",
    text: BANDS.map(b => counts[b]), textposition: "outside",
  }], baseLayout({margin:{l:50,r:20,t:8,b:40},
                  yaxis:{title:{text:"Asset-hazard combinations"},
                         gridcolor:ink("--grid"), zeroline:false}}), CFG);
}

function band(p) {
  if (p === null || !isFinite(p)) return null;
  if (p <= 20) return BANDS[0];
  if (p <= 40) return BANDS[1];
  if (p <= 60) return BANDS[2];
  if (p <= 80) return BANDS[3];
  return BANDS[4];
}

// ---- time series --------------------------------------------------------------------
// Small multiples for one location: the Climate Score panel FIRST, then one panel per
// hazard. The Score is the mean of the hazard percentiles beside it, so putting them in one
// grid lets a reader see which hazard drives the aggregate without switching views.
//
// The metric toggle (Percentile / Value) applies to the HAZARD panels only. The Climate
// Score is a 1-100 score by construction and has no native-unit analogue, so its panel is
// unchanged by the toggle and its axis title says so.
//
// ONE render path, ONE Plotly lifecycle. The previous design swapped between two functions
// that both drew into #series, and clearing innerHTML without Plotly.purge left stale
// internal state bound to the div -- Plotly.react then failed silently and the panel went
// blank. Never clear that node by hand; go through resetSeries().
function resetSeries(placeholder) {
  Plotly.purge("series");                 // drop Plotly's state BEFORE touching the DOM
  $("series").innerHTML = placeholder || "";
}

function locationScoreTraces() {
  const tiers = TIER_ORDER.filter(t => DATA.scores.some(r => r.t === t));
  const allDecades = [...new Set(DATA.scores.map(r => r.d))].sort((a,b)=>a-b);
  // Balanced panel, as everywhere a Climate Score is rolled up: an asset that drops out of
  // one tier would otherwise make that tier's line a different asset mix.
  const decades = allDecades.filter(d => balancedAssets(tiers, [d], state.loc).length > 0);
  const panel = new Set(balancedAssets(tiers, decades, state.loc));
  const candidates = DATA.assets.filter(a => a.loc === state.loc
      && (state.asset === "__all" || a.type === state.asset)).length;

  const traces = tiers.map(t => {
    const xs = [], ys = [];
    decades.forEach(d => {
      const rs = DATA.scores.filter(r => panel.has(r.a) && r.t === t && r.d === d);
      if (rs.length) { xs.push(d); ys.push(mean(rs.map(r => r.s))); }
    });
    return {type: "scatter", mode: "lines+markers", name: TIER_LABELS[t],
            x: xs, y: ys, ...tierStyle(t), legendgroup: t,
            hovertemplate: `${TIER_LABELS[t]}<br>%{x}s: %{y:.1f}<extra></extra>`};
  });
  return {traces, panelSize: panel.size, candidates};
}

function renderSeries() {
  const box = $("series");
  if (!state.seriesVisible) {
    Plotly.purge("series");            // free the traces when hidden, not just hide them
    box.innerHTML = "";
    box.style.display = "none";
    return;
  }
  box.style.display = "";
  const ll = locLayers(state.asset);
  const have = ll[state.loc] ? [...ll[state.loc]] : [];
  if (!have.length) {
    resetSeries('<div style="color:var(--text-muted);padding:24px">No hazards for this ' +
                'location under the current asset-type filter.</div>');
    return;
  }
  resetSeries();

  const key = state.seriesMetric === "value" ? "v" : "p";
  const score = locationScoreTraces();
  // An empty balanced panel means the Score cannot be shown on a common basis here; the
  // hazard panels are still valid, so drop only the Score panel rather than the grid.
  const hasScore = score.panelSize > 0 && score.traces.some(t => t.x.length);
  // Panel 0 is the Climate Score when it can be computed; hazards follow, sorted by
  // hazard label so related panels sit together and a reader can scan for one.
  have.sort((a, b) => (DATA.layers[a].hazard + a).localeCompare(DATA.layers[b].hazard + b));
  const layerDecs = {};
  have.forEach(lid => {
    layerDecs[lid] = [...new Set(DATA.values.filter(v => v.lay === lid).map(v => v.dec))]
        .sort((a,b) => a-b);
  });
  const panels = (hasScore ? [{kind: "score"}] : [])
      .concat(have.map(lid => ({kind: "hazard", lid,
                                single: layerDecs[lid].length === 1})));

  // With the v2 standard set this grid holds 30+ panels. Cramming them into a fixed-height
  // box made every panel unreadable (user call, 2026-08-20): the height scales with the
  // row count -- ~290px per panel, with enough inter-row gap that one row's x labels
  // never collide with the next row's title (user call, 2026-08-21) -- and the page
  // simply scrolls.
  const cols = panels.length <= 2 ? panels.length : (panels.length <= 4 ? 2 : 3);
  const rowsN = Math.ceil(panels.length / cols);
  $("series").style.height = (rowsN * 290 + 140) + "px";
  const traces = [], annos = [];
  let legendDone = false;

  panels.forEach((pan, i) => {
    const ax = i === 0 ? "" : (i + 1);
    if (pan.kind === "score") {
      score.traces.forEach(tr => {
        traces.push(Object.assign({}, tr, {
          xaxis: "x" + ax, yaxis: "y" + ax, showlegend: !legendDone}));
      });
      legendDone = true;
      annos.push({text: "Climate Score" +
          (score.panelSize < score.candidates
            ? ` · ${score.panelSize}/${score.candidates} assets` : ""),
        xref: "x" + ax + " domain", yref: "y" + ax + " domain",
        x: 0, y: 1.14, showarrow: false, xanchor: "left",
        font: {size: 12, color: ink("--text-primary")}});
      return;
    }
    const lay = DATA.layers[pan.lid];
    const here = DATA.values.filter(v => v.loc === state.loc && v.lay === pan.lid);
    const anyVal = here.some(v => v[key] !== null);
    const statuses = [...new Set(here.map(v => v.st))];
    // A hazard absent at the site (off its layer's footprint) plots as a flat zero series
    // rather than an empty panel (user call, 2026-08-21) -- but labelled as ABSENT, so a
    // shown-as-0 absence never reads as a measured zero. Offshore/off-grid stays an
    // explicit "no data" note: there the value is unknown, not absent.
    const offMask = !anyVal && statuses.length === 1 && statuses[0] === "OFF_LAYER_MASK";
    // An OBSERVED layer with no data here (tornado outside CONUS) means UNOBSERVED, not
    // absent -- its panel stays blank rather than drawing a zero (user call, 2026-08-21).
    const obsNoCover = offMask && (lay.scenarios || []).includes("observed");
    const hover = pan.single ? "Current" : "%{x}s";

    // ALL scenarios, ALWAYS -- the forcing-tier filter never reaches a time series.
    const scns = [...new Set(DATA.values.filter(v => v.lay === pan.lid).map(v => v.scn))]
        .sort((a,b) => TIER_ORDER.indexOf(tierOf(a)) - TIER_ORDER.indexOf(tierOf(b)));
    if (!offMask) scns.forEach(scn => {
      const pts = DATA.values
        .filter(v => v.loc === state.loc && v.lay === pan.lid && v.scn === scn)
        .sort((a,b) => a.dec - b.dec);
      traces.push({
        type: "scatter", mode: pan.single ? "markers" : "lines+markers", name: scn,
        x: pts.map(p => p.dec), y: pts.map(p => p[key]),
        xaxis: "x" + ax, yaxis: "y" + ax,
        // Colour follows the forcing TIER, so a scenario keeps its hue across panels even
        // though RCP and SSP layers use different codes.
        ...tierStyle(tierOf(scn)),
        legendgroup: tierOf(scn), showlegend: !legendDone,
        hovertemplate: `${esc(scn)}<br>${hover}: %{y:.4g}<extra></extra>`,
      });
    });
    if (offMask && !obsNoCover) {
      const decs = layerDecs[pan.lid];
      traces.push({
        type: "scatter", mode: decs.length > 1 ? "lines+markers" : "markers",
        x: decs, y: decs.map(() => 0), xaxis: "x" + ax, yaxis: "y" + ax,
        showlegend: false,
        line: {color: ink("--text-muted"), width: 1.5, dash: "dot"},
        marker: {size: 6, symbol: "circle-open", color: ink("--text-muted")},
        hovertemplate: "hazard absent at this site (off footprint) — shown as 0<extra></extra>",
      });
    }
    legendDone = true;
    annos.push({text: lay.hazard + (state.seriesMetric === "value"
                    ? (lay.units ? ` [${lay.units}]` : "") : " [percentile]") +
                    (obsNoCover ? " — no coverage at this site"
                     : offMask ? " — absent here, shown as 0" : ""),
                xref: "x" + ax + " domain", yref: "y" + ax + " domain",
                x: 0, y: 1.14, showarrow: false, xanchor: "left",
                font: {size: 12, color: offMask ? ink("--text-muted")
                                               : ink("--text-primary")}});

    // An offshore/off-grid site still produces a panel with axes and no line, which looks
    // exactly like a rendering bug -- it has already been reported as one. Say why.
    if (!anyVal && !offMask) {
      annos.push({text: `No data at this site<br>(${esc(statuses.join(", ")) || "no rows"})`,
                  xref: "x" + ax + " domain", yref: "y" + ax + " domain",
                  x: 0.5, y: 0.5, showarrow: false, xanchor: "center", align: "center",
                  font: {size: 12, color: ink("--text-muted")}});
    }
  });

  const layout = baseLayout({
    grid: {rows: rowsN, columns: cols, pattern: "independent",
           roworder: "top to bottom", xgap: 0.14, ygap: rowsN > 4 ? 0.26 : 0.34},
    annotations: annos, showlegend: true,
    legend: {orientation: "h", y: -0.14, x: 0, font: {color: ink("--text-secondary")}},
    margin: {l: 56, r: 16, t: 34, b: 62},
    // No selected-decade marker lines here (removed 2026-08-21, user call): with
    // yref "paper" they ran the full height of the grid, striking through every row
    // and its titles.
  });
  panels.forEach((pan, i) => {
    const ax = i === 0 ? "" : (i + 1);
    layout["xaxis" + ax] = {gridcolor: ink("--grid"), linecolor: ink("--axis"),
                            zeroline: false, tickformat: "d"};
    // A single-period (observational) hazard is one reading, not a series: one marker,
    // labelled "Current" rather than the window's start year (user call, 2026-08-21).
    // The true observed window is in layers.csv.
    if (pan.kind === "hazard" && pan.single) {
      layout["xaxis" + ax].tickvals = layerDecs[pan.lid];
      layout["xaxis" + ax].ticktext = ["Current"];
      layout["xaxis" + ax].range = [layerDecs[pan.lid][0] - 6, layerDecs[pan.lid][0] + 6];
    }
    layout["yaxis" + ax] = {gridcolor: ink("--grid"), linecolor: ink("--axis"),
                            zeroline: false};
    // Score and percentile share the 1-100 risk axis, so pin both for comparability.
    if (pan.kind === "score" || state.seriesMetric === "percentile")
      layout["yaxis" + ax].range = [0, 105];
  });
  Plotly.react("series", traces, layout, CFG);
}

// ---- table --------------------------------------------------------------------------
// SELF-CONTAINED SCOPING (2026-08-21). The Values table answers to its OWN controls
// only -- the per-column dropdowns and the search box; the top filter bar does not reach
// it. Earlier wirings scoped it by the bar's decade and tier, so a column-filter pick
// could silently intersect with a bar filter to zero rows and read as "the filter does
// not work". This is the third stated exception to the one-filter-row rule.
//
// STABLE HEADER (same date). The header -- label row and filter row -- is built ONCE in
// initTable() and never rebuilt: only <tbody> re-renders on a change. Rebuilding the
// whole table replaced the very <select> the user was interacting with while its change
// event was still dispatching, which is the likeliest culprit for the reported
// filter-resets, and it re-created every dropdown on every keystroke besides.
//
// [label, getter, numeric, filterable]. "Agree" and "Status" dropped 2026-08-21 (user
// calls); both remain in values.csv. Decade, Percentile and Members render as integers.
const COLS = [
  ["Location", r => locById[r.loc].name, false, true],
  ["Hazard", r => DATA.layers[r.lay].hazard, false, true],
  ["Layer", r => r.lay, false, true],
  ["Scenario", r => r.scn, false, true],
  ["Decade", r => r.dec, true, true],
  ["Value", r => r.v, true, false],
  ["25th", r => r.lo, true, false],
  ["75th", r => r.hi, true, false],
  ["Percentile", r => r.p, true, false],
  ["OLS", r => r.ols, true, false],
  ["Sen", r => r.sen, true, false],
  ["Members", r => r.nm, true, false],
];
const INT_COLS = new Set(["Decade", "Percentile", "Members"]);

// Rows passing the search plus every column filter EXCEPT skipIdx (pass -1 for all).
function rowsPassing(skipIdx) {
  let base = DATA.values;
  if (state.search) {
    base = base.filter(r => (locById[r.loc].name + " " + DATA.layers[r.lay].hazard + " " +
        r.lay + " " + r.scn).toLowerCase().includes(state.search));
  }
  return base.filter(r => COLS.every((c, i) => {
    if (i === skipIdx || !c[3]) return true;
    const f = state.colFilters[i];
    return f === undefined || f === "__all" || String(c[1](r)) === f;
  }));
}
const tableRows = () => rowsPassing(-1);

const optSort = (a, b) => {
  const na = Number(a), nb = Number(b);
  const an = isFinite(na), bn = isFinite(nb);
  if (an && bn) return na - nb;
  if (an !== bn) return an ? -1 : 1;
  return a.localeCompare(b);
};

// CASCADING COLUMN FILTERS (user call, 2026-08-21): picking in one dropdown narrows the
// OPTIONS of the others to values that still have rows (each list is computed ignoring
// only its own filter, so a choice never deletes its own alternatives). The select being
// changed is never rebuilt mid-event; a pick invalidated upstream resets to All.
function updateFilterOptions(changedIdx) {
  const t = $("table");
  COLS.forEach((c, i) => {
    if (!c[3] || i === changedIdx) return;
    const sel = t.querySelector(`select[data-f="${i}"]`);
    if (!sel) return;
    const vals = [...new Set(rowsPassing(i).map(r => String(c[1](r))))].sort(optSort);
    let cur = state.colFilters[i] ?? "__all";
    if (cur !== "__all" && !vals.includes(cur)) { cur = "__all"; state.colFilters[i] = cur; }
    sel.innerHTML = `<option value="__all">All</option>` +
      vals.map(v => `<option value="${esc(v)}">${esc(v)}</option>`).join("");
    sel.value = cur;
  });
}

function initTable() {
  const t = $("table");
  const head = "<tr>" + COLS.map((c, i) =>
      `<th class="sorth" data-i="${i}">${esc(c[0])}</th>`).join("") + "</tr>";
  const filt = "<tr class='colf'>" + COLS.map((c, i) => {
      if (!c[3]) return "<th></th>";
      // Options come from the FULL dataset, once. A static list can produce an honest
      // "0 rows" for an unusual combination, which beats options that churn under the
      // user's pointer.
      const vals = [...new Set(DATA.values.map(r => String(c[1](r))))].sort(optSort);
      return `<th><select data-f="${i}"><option value="__all">All</option>` +
        vals.map(v => `<option value="${esc(v)}">${esc(v)}</option>`).join("") +
        `</select></th>`;
    }).join("") + "</tr>";
  t.innerHTML = "<thead>" + head + filt + "</thead><tbody></tbody>";
  t.onclick = e => {
    const th = e.target.closest("th.sorth");
    if (!th || !t.contains(th)) return;
    const i = +th.dataset.i;
    state.sortDir = state.sortCol === i ? -state.sortDir : 1;
    state.sortCol = i;
    renderTable();
  };
  t.onchange = e => {
    const sel = e.target.closest("select[data-f]");
    if (!sel || !t.contains(sel)) return;
    const i = +sel.dataset.f;
    state.colFilters[i] = sel.value;
    updateFilterOptions(i);
    renderTable();
  };
}

function renderTable() {
  const t = $("table");
  if (!t.tHead) initTable();
  let data = tableRows();
  if (state.sortCol !== null) {
    const get = COLS[state.sortCol][1];
    data = data.slice().sort((a, b) => {
      const x = get(a), y = get(b);
      if (x === null) return 1; if (y === null) return -1;
      return (x > y ? 1 : x < y ? -1 : 0) * state.sortDir;
    });
  }
  // Sort indicators update IN PLACE on the label row; the filter row is never touched.
  t.tHead.rows[0].querySelectorAll("th").forEach((th, i) => {
    th.textContent = COLS[i][0] +
      (state.sortCol === i ? (state.sortDir > 0 ? " ▲" : " ▼") : "");
  });
  // ROW CAP: keeps the DOM small however broad the filters; truncation is stated in the
  // table itself, never silent.
  const CAP = 800;
  const shown = data.slice(0, CAP);
  const trunc = data.length > CAP
    ? `<tr><td colspan="${COLS.length}" style="color:var(--text-muted)">Showing the ` +
      `first ${CAP} of ${data.length} rows — narrow the filters or search to see the ` +
      `rest. The CSVs carry every row.</td></tr>`
    : "";
  t.tBodies[0].innerHTML = shown.map(r => "<tr>" + COLS.map(c => {
      const v = c[1](r);
      const txt = c[2]
        ? (typeof v === "number"
            ? (INT_COLS.has(c[0]) ? String(Math.round(v)) : fmt(v))
            : (v ?? "—"))
        : (v ?? "—");
      return `<td class="${c[2] ? "num" : ""}">${esc(txt)}</td>`;
    }).join("") + "</tr>").join("") + trunc;
}

function renderAll() {
  renderTiles(); renderMap(); renderBars(); renderScoreSeries();
  renderSeries(); renderTable();
}

initControls();
renderAll();
$("foot").textContent =
  "Full guidance on reading this delivery: README.md, alongside the CSVs. " +
  "Quick definitions: FAQ, top right.";
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", renderAll);
"""


def build_dashboard(delivery: Path, quiet: bool = False):
    """Build the dashboard for a finished delivery. Returns (path, build_id, warnings).

    The importable entry point: `generate_customer_delivery.py` calls this at the end of
    every run so a delivery is never shipped as CSVs alone. Keep it free of argparse and
    sys.exit so it composes.
    """
    locations, assets, layers, values, scores, manifest = load_delivery(delivery)
    payload = build_payload(locations, assets, layers, values, scores)

    warnings = check_tier_collisions(
        {lid: m["scenarios"] for lid, m in payload["layers"].items()}
    )
    # A layer with NO forcing scenario at all is observational, not a gap in a projected
    # layer. Saying it "has no high/low-tier scenario" invites someone to go looking for
    # the missing files; there are none to find, and the cross-tier panel is simply not a
    # question this layer answers.
    tiers_seen = {lid: {tier_of(s) for s in meta["scenarios"] if is_forcing_scenario(s)}
                  for lid, meta in payload["layers"].items()}
    observational = sorted(lid for lid, t in tiers_seen.items() if not t)
    union = set().union(*(t for t in tiers_seen.values() if t)) if tiers_seen else set()
    for lid, present in sorted(tiers_seen.items()):
        if lid in observational:
            continue
        gap = union - present
        if gap:
            warnings.append(
                f"{lid} has no {'/'.join(sorted(gap))}-tier scenario, so assets carrying "
                f"it drop out of the balanced panel for cross-tier comparisons")
    for lid in observational:
        warnings.append(
            f"{lid} is OBSERVATIONAL (scenario "
            f"{'/'.join(sorted(payload['layers'][lid]['scenarios']))}): it carries no "
            f"forcing pathway, so it is shown on its own and excluded from cross-tier "
            f"panels. Where its hazard family is weighted for an asset, it enters the "
            f"Climate Score as a CONSTANT across every decade and tier (user decision "
            f"2026-08-20), which damps the score's scenario and time contrast in "
            f"proportion to its weight share. This is not a missing scenario.")
    unknown = sorted(s for s in set(payload["tiers"]) - set(SCENARIO_TIER)
                     if is_forcing_scenario(s))
    if unknown:
        warnings.append(
            f"scenario code(s) not in SCENARIO_TIER, defaulted to 'medium': "
            f"{', '.join(unknown)}")
    if not len(scores):
        warnings.append("no climate_score.csv -- the Climate Score views will be absent")

    out, build = render(delivery, payload, manifest, warnings)
    if not quiet:
        print(f"  dashboard: {out.name}  ({out.stat().st_size // 1024} KB)  build {build}")
        for w in warnings:
            print(f"    WARNING: {w}")
    return out, build, warnings


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("delivery", type=Path, help="Delivery folder")
    p.add_argument("--check-stamp", action="store_true",
                   help="Print the build stamp of the dashboard already on disk and exit. "
                        "Compare it with the stamp shown in the page header: a mismatch "
                        "means the browser is serving a cached copy.")
    args = p.parse_args()

    if args.check_stamp:
        f = args.delivery / "dashboard.html"
        if not f.exists():
            print(f"No dashboard at {f}")
            return 1
        html = f.read_text()
        m = re.search(r"dashboard build <code[^>]*>([0-9a-f]{8})</code>", html)
        stamp = m.group(1) if m else "NONE (built before build stamps existed)"
        mtime = datetime.fromtimestamp(f.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        print(f"on disk: build {stamp}   written {mtime}   {f.stat().st_size // 1024} KB")
        print("If the page header shows a different build, hard-reload "
              "(Cmd-Shift-R) — the browser is serving a cached copy.")
        return 0

    out, build, _ = build_dashboard(args.delivery, quiet=True)
    size_kb = out.stat().st_size / 1024
    print(f"Wrote {out}  ({size_kb:.0f} KB)  build {build}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
