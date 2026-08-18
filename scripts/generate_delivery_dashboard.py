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

  * Time series ignores the DECADE filter -- it plots every decade, and marks the selected
    one. Filtering it by decade would leave a single point.
  * Summaries ignore the HAZARD filter -- "impacts by hazard" is the chart.

Both are stated on the page so a reader is never guessing what a chart is showing.
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
        "scores": [
            {"a": r["asset_id"], "t": r["scenario_tier"], "d": int(r["decade"]),
             "s": float(r["climate_score"]), "n": int(r["n_hazards"])}
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
    unreviewed = [lid for lid, m in payload["layers"].items() if m["qa"] == "NOT CONFIRMED"]

    banner = ""
    if unreviewed:
        banner = (
            '<div class="banner">'
            f'<strong>QA review not confirmed for {esc(", ".join(unreviewed))}.</strong> '
            'Passing the output contract means a file is shaped right, not that its input '
            'is about what its name says. Treat these numbers as provisional.'
            '</div>'
        )
    if warnings:
        banner += '<div class="banner warn"><strong>Build warnings:</strong> ' + \
                  esc("; ".join(warnings)) + '</div>'

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
header {{ padding: 20px 24px 12px; }}
header h1 {{ margin:0 0 4px; font-size: 20px; font-weight: 600; }}
header .sub {{ color: var(--text-secondary); font-size: 13px; }}
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
table {{ border-collapse:collapse; width:100%; font-size:12px;
        font-variant-numeric: tabular-nums; }}
th, td {{ text-align:left; padding:5px 8px; border-bottom:1px solid var(--grid);
         white-space:nowrap; }}
th {{ position:sticky; top:0; background: var(--surface); cursor:pointer;
     color: var(--text-secondary); font-weight:600; z-index:1; }}
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
  <div class="sub">Delivery {delivery.name} &nbsp;·&nbsp; generated {generated}
      &nbsp;·&nbsp; <span id="counts"></span></div>
  <div class="sub" style="margin-top:4px">
      dashboard build <code style="font-variant-numeric:tabular-nums">{build}</code>
      &nbsp;·&nbsp; {built_at}
      &nbsp;·&nbsp; <span style="color:var(--text-muted)">if this does not match
      <code>--check-stamp</code>, the page is cached — hard-reload</span></div>
</header>
{banner}

<div class="filters">
  <div><label for="f-tier">Forcing tier</label><select id="f-tier"></select></div>
  <div><label for="f-hazard">Hazard</label><select id="f-hazard"></select></div>
  <div><label for="f-asset">Asset type</label><select id="f-asset"></select></div>
  <div><label for="f-decade">Decade</label><select id="f-decade"></select></div>
  <div style="color:var(--text-muted);font-size:12px;max-width:330px">
    Scopes every view. Time series shows all decades; summaries show all hazards.
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

  <div class="card"><h2>Mean risk percentile by hazard</h2>
    <div class="note">All hazards, ignoring the hazard filter. Selected decade and tier.</div>
    <div id="bar-hazard" style="height:300px"></div></div>

  <div class="card"><h2>Mean risk percentile by asset class</h2>
    <div class="note">Averaged over each asset's own hazards.</div>
    <div id="bar-asset" style="height:300px"></div></div>

  <div class="card"><h2>Portfolio risk distribution</h2>
    <div class="note" id="band-note"></div>
    <div id="bar-band" style="height:300px"></div></div>

  <div class="card"><h2>Climate Score over time</h2>
    <div class="toggle" id="score-toggle"></div>
    <div class="note" id="score-note"></div>
    <div id="score-series" style="height:300px"></div></div>

  <div class="card full">
    <h2>Time series by location</h2>
    <div class="note">One location at a time. The <strong>Climate Score</strong> panel
        comes first, then one panel per hazard — the Score is the mean of those hazard
        percentiles, so they read together. The toggle switches the hazard panels only;
        the Score is a 1–100 score by construction.
        <strong>All decades and all scenarios are always shown</strong> — neither the
        decade nor the forcing-tier filter applies here, because seeing the full spread is
        the point of the chart.</div>
    <div style="margin-bottom:10px">
      <label for="f-loc" style="font-size:11px;text-transform:uppercase;
          letter-spacing:.04em;color:var(--text-muted)">Location</label>
      <select id="f-loc"></select>
      <span class="toggle" id="series-toggle" style="margin-left:10px"></span>
    </div>
    <div id="series" style="height:420px"></div>
  </div>

  <div class="card full">
    <h2>Values</h2>
    <div class="note">Every row in the filtered slice. Click a header to sort.</div>
    <input type="search" id="f-search" placeholder="Filter rows…" style="margin-bottom:10px">
    <div class="scroll"><table id="table"></table></div>
  </div>
</div>

<footer id="foot"></footer>

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

// ---- filter state -------------------------------------------------------------------
// Climate Score leads: it is the portfolio-level answer, and every other metric is a
// component of it. The hazard filter does not apply while it is selected.
const state = {tier: "__all", hazard: null, asset: "__all", decade: null,
               metric: "climate_score", loc: null, seriesMetric: "percentile",
               scoreMetric: "__score", search: "", sortCol: null, sortDir: 1};

// --- Climate Score helpers ------------------------------------------------------------
// An asset's score is complete only when every hazard in its catalog set contributed.
// ISIMIP3b layers have no 2010s panel, so a 2010s row can legitimately rest on ONE hazard
// while the 2020s rests on three -- plotting both on one axis would read as risk tripling
// when only the hazard set changed. Incomplete rows are excluded and counted, never
// silently averaged in.
const assetById = Object.fromEntries(DATA.assets.map(a => [a.id, a]));
function expectedHazards(assetId) {
  const a = assetById[assetId];
  return a ? a.layers.length : 0;
}
function scoreRows({tier, decade, completeOnly = true} = {}) {
  return DATA.scores.filter(r => {
    const a = assetById[r.a];
    if (!a) return false;
    if (state.asset !== "__all" && a.type !== state.asset) return false;
    if (tier && tier !== "__all" && r.t !== tier) return false;
    if (decade !== undefined && decade !== null && r.d !== decade) return false;
    if (completeOnly && r.n < expectedHazards(r.a)) return false;
    return true;
  });
}
const mean = xs => xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : null;

// BALANCED PANEL -- the single most important guard on this page.
//
// Layer scenario coverage is NOT uniform across tiers: `cyclone` publishes rcp26/rcp60 and
// no rcp85, so every cyclone-carrying asset is INCOMPLETE at the high tier. Averaging each
// tier over whatever assets it happens to have makes the tier lines describe different
// portfolios. Measured on the example delivery: the high-tier 2020s mean read 39.9 against
// 42.1 for low and medium -- entirely a composition artifact, since the shared 2020s
// baseline is bit-identical across scenarios and the three MUST be equal on a common basis.
//
// So any chart that compares across tiers or decades first restricts to the assets present
// and complete in EVERY cell of that grid, and states what it dropped.
function balancedAssets(tiers, decades, locId) {
  const cand = new Set(DATA.assets
    .filter(a => (state.asset === "__all" || a.type === state.asset)
                 && (!locId || a.loc === locId))
    .map(a => a.id));
  const have = new Set();
  DATA.scores.forEach(r => {
    if (r.n >= expectedHazards(r.a)) have.add(r.a + "|" + r.t + "|" + r.d);
  });
  return [...cand].filter(id =>
    tiers.every(t => decades.every(d => have.has(id + "|" + t + "|" + d))));
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
  fill($("f-hazard"), layerIds.map(l => [l, DATA.layers[l].hazard + " (" + l + ")"]));
  fill($("f-asset"), [["__all","All asset types"]].concat(assetTypes.map(a => [a,a])));
  fill($("f-decade"), decades.map(d => [d, d + "s"]));
  fill($("f-loc"), DATA.locations.map(l => [l.id, l.name]));

  state.hazard = layerIds[0];
  state.decade = decades.includes(2050) ? 2050 : decades[decades.length-1];
  state.loc = DATA.locations[0].id;
  $("f-hazard").value = state.hazard;
  $("f-decade").value = state.decade;
  $("f-loc").value = state.loc;

  $("f-tier").onchange   = e => { state.tier = e.target.value; renderAll(); };
  $("f-hazard").onchange = e => { state.hazard = e.target.value; renderAll(); };
  $("f-asset").onchange  = e => { state.asset = e.target.value; renderAll(); };
  $("f-decade").onchange = e => { state.decade = +e.target.value; renderAll(); };
  $("f-loc").onchange    = e => { state.loc = e.target.value; renderSeries(); };
  $("f-search").oninput  = e => { state.search = e.target.value.toLowerCase(); renderTable(); };

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
  // (1-100 risk percentile) either way, so the two are directly comparable.
  buildToggle("score-toggle",
    [["__score","Climate Score"]].concat(layerIds.map(l => [l, DATA.layers[l].hazard])),
    () => state.scoreMetric, v => { state.scoreMetric = v; renderScoreSeries(); });

  $("counts").textContent = DATA.locations.length + " locations · " +
      DATA.assets.length + " assets · " + Object.keys(DATA.layers).length + " layers · " +
      DATA.values.length + " site-hazard-scenario-decade rows";
}

function fill(sel, pairs) {
  sel.innerHTML = pairs.map(([v,l]) =>
    `<option value="${esc(v)}">${esc(l)}</option>`).join("");
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
    note = `Diverging scale centred on zero, ±${fmt(L)} (${rule})` +
           (clamped ? `, ${clamped} site(s) clamped` : "") +
           ` · ${lay.slope_units || "per decade"} · red = worsening` +
           (rev ? " (scale reversed: this layer is higher_is_better)" : "") +
           ` · recommended slope for this layer: ${lay.recommended}`;
  } else if (state.metric === "percentile") {
    cs = riskScale(); cmin = 1; cmax = 100;
    note = "Percentile 1–100 against the shared 2020s baseline. Already oriented for " +
           "risk — 100 is worst on every layer, including higher_is_better ones. Darker = worse.";
  } else {
    cs = SEQ_VALUE; cmin = null; cmax = null;
    note = `${lay.measure} [${lay.units || "—"}] · statistic: ${lay.statistic} · darker = larger`;
  }

  const nonRobust = pts.filter(r => r.ag === false).length;
  if (isSlope && nonRobust)
    note += ` · ${nonRobust}/${pts.length} site(s) have non-robust trends (slopes disagree)`;
  $("map-note").textContent = note + (lay.note ? " — " + lay.note : "");

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
      line: {width: 2, color: ink("--surface")},   // 2px surface ring, not a border
      colorbar: {thickness: 12, len: 0.7, outlinewidth: 0,
                 tickfont: {color: ink("--text-secondary")}},
    },
  };
  Plotly.react("map", [trace], baseLayout({
    margin: {l: 0, r: 0, t: 0, b: 0},
    geo: {showland: true, landcolor: ink("--grid"), showcountries: true,
          countrycolor: ink("--surface"), showocean: true, oceancolor: ink("--plane"),
          coastlinecolor: ink("--axis"), coastlinewidth: 0.5,
          projection: {type: "natural earth"}, bgcolor: "rgba(0,0,0,0)"},
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
    haz: Math.max(...rs.map(x => x.n)),
    tiers: [...new Set(rs.map(x => x.t))],
  }));

  const dropped = scoreRows({tier: state.tier, decade: state.decade, completeOnly: false})
      .length - rowsIn.length;
  $("map-note").textContent =
    `Mean risk percentile across all of a location's hazards, ${state.decade}s` +
    (state.tier === "__all" ? ", averaged over all forcing tiers" :
       ` under ${TIER_LABELS[state.tier].toLowerCase()}`) +
    ". Higher = higher aggregate physical risk. Darker = worse. " +
    "The hazard filter does not apply to this metric. " +
    "Keyed on forcing tier, not a native scenario code — no ISIMIP code spans both " +
    "rounds, so a per-code score would cover only a subset of hazards." +
    (dropped ? ` ${dropped} asset-row(s) excluded for incomplete hazard coverage.` : "");

  const trace = {
    type: "scattergeo", mode: "markers",
    lat: pts.map(p => locById[p.loc].lat), lon: pts.map(p => locById[p.loc].lon),
    text: pts.map(p => {
      const l = locById[p.loc];
      return `<b>${esc(l.name)}</b><br>Climate Score ${p.s.toFixed(1)}` +
        `<br>${p.haz} hazard(s) · ${p.n} asset-row(s)` +
        `<br>tier(s): ${esc(p.tiers.join(", "))}`;
    }),
    hovertemplate: "%{text}<extra></extra>",
    marker: {
      size: 15, color: pts.map(p => p.s), colorscale: riskScale(), cmin: 1, cmax: 100,
      line: {width: 2, color: ink("--surface")},
      colorbar: {thickness: 12, len: 0.7, outlinewidth: 0,
                 tickfont: {color: ink("--text-secondary")}},
    },
  };
  Plotly.react("map", [trace], baseLayout({
    margin: {l: 0, r: 0, t: 0, b: 0},
    geo: {showland: true, landcolor: ink("--grid"), showcountries: true,
          countrycolor: ink("--surface"), showocean: true, oceancolor: ink("--plane"),
          coastlinecolor: ink("--axis"), coastlinewidth: 0.5,
          projection: {type: "natural earth"}, bgcolor: "rgba(0,0,0,0)"},
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
    `Unweighted mean of risk percentile across each asset's hazards. Scoped by the ` +
    `asset-type, forcing-tier and decade filters; the hazard filter does not apply. ` +
    `Assets without full hazard coverage are excluded rather than scored on a subset.` +
    (gaps.length ? ` Note: ${gaps.join("; ")} — assets carrying those hazards cannot be ` +
      `scored at the missing tier at all.` : "");
  $("tiles").innerHTML = [
    tile("Portfolio Climate Score", score === null ? "—" : score.toFixed(1),
         `${state.decade}s · ` + (state.tier === "__all" ? "all tiers"
            : TIER_LABELS[state.tier].toLowerCase()) +
         ` · ${new Set(cur.map(r => r.a)).size} assets`),
    tile("Change vs " + baseDec + "s", delta === null ? "—" :
         (delta >= 0 ? "+" : "") + delta.toFixed(1),
         delta === null ? "insufficient common coverage"
            : `${deltaPanel.size} assets on a common basis · ` +
              (delta > 0 ? "rising risk" : delta < 0 ? "falling risk" : "flat")),
    tile("Highest-risk location", worst ? worst[1].toFixed(1) : "—",
         worst ? locById[worst[0]].name : ""),
    tile("Coverage", nAssets + " assets", nHaz + " hazards · " + cur.length + " score rows"),
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
    const perDecade = decades.map(d =>
      `${d}s: ${balancedAssets(tiers, [d]).length}`).join(", ");
    $("score-note").textContent =
      `No asset has complete hazard coverage in every tier across all decades shown, so ` +
      `there is no common basis to plot. Assets qualifying in individual decades — ` +
      `${perDecade || "none"}. Averaging those would compare different asset sets across ` +
      `time. Narrow the asset-type filter to a class with uniform coverage.` +
      (gaps.length ? ` Cause: ${gaps.join("; ")}.` : "");
    $("score-series").innerHTML = "";
    return;
  }

  $("score-note").textContent =
    `Portfolio mean over a BALANCED PANEL of ${panel.size} of ${candidates} asset(s) — ` +
    `those with complete hazard coverage in every tier and decade shown, so the three ` +
    `lines describe the same portfolio. Scoped by the asset-type filter only; the ` +
    `forcing-tier and decade filters do not apply, because varying them is what this ` +
    `chart shows.` +
    (droppedDec.length ? ` ${droppedDec.map(d=>d+"s").join(", ")} omitted (no asset has ` +
      `full coverage there — ISIMIP3b layers have no 2010s panel).` : "") +
    (panel.size < candidates && gaps.length
      ? ` ${candidates - panel.size} asset(s) excluded because ${gaps.join("; ")}.` : "");

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
    `${lay.hazard} — portfolio mean risk percentile across ${locs.size} site(s) ` +
    `carrying this hazard. Same 1–100 axis as the Climate Score, which is the mean of ` +
    `these across an asset's hazards. Scoped by the asset-type filter only.` +
    (missing.length ? ` No ${missing.join("/")}-tier line: this layer publishes no such ` +
      `scenario (${lay.scenarios.join(", ")}).` : "") +
    (lay.direction === "higher_is_better"
      ? " This layer is higher_is_better — the percentile is already inverted, so a rising"
      + " line still means rising risk." : "");

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
    marker: {color: tierColor('low')},  // one series -> one colour; never a ramp on nominal bars
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
    `${n} asset-hazard-scenario combinations in the ${state.decade}s, binned from ` +
    `percentile. Bands are a DISPLAY derivation — they are deliberately not stored in ` +
    `values.csv, which carries the percentile they come from.`;
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
  // Panel 0 is the Climate Score when it can be computed; hazards follow.
  const panels = (hasScore ? [{kind: "score"}] : [])
      .concat(have.map(lid => ({kind: "hazard", lid})));

  const cols = panels.length <= 2 ? panels.length : (panels.length <= 4 ? 2 : 3);
  const rowsN = Math.ceil(panels.length / cols);
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
    // ALL scenarios, ALWAYS -- the forcing-tier filter never reaches a time series.
    const scns = [...new Set(DATA.values.filter(v => v.lay === pan.lid).map(v => v.scn))]
        .sort((a,b) => TIER_ORDER.indexOf(tierOf(a)) - TIER_ORDER.indexOf(tierOf(b)));
    scns.forEach(scn => {
      const pts = DATA.values
        .filter(v => v.loc === state.loc && v.lay === pan.lid && v.scn === scn)
        .sort((a,b) => a.dec - b.dec);
      traces.push({
        type: "scatter", mode: "lines+markers", name: scn,
        x: pts.map(p => p.dec), y: pts.map(p => p[key]),
        xaxis: "x" + ax, yaxis: "y" + ax,
        // Colour follows the forcing TIER, so a scenario keeps its hue across panels even
        // though RCP and SSP layers use different codes.
        ...tierStyle(tierOf(scn)),
        legendgroup: tierOf(scn), showlegend: !legendDone,
        hovertemplate: `${esc(scn)}<br>%{x}s: %{y:.4g}<extra></extra>`,
      });
    });
    legendDone = true;
    annos.push({text: lay.hazard + (state.seriesMetric === "value"
                    ? (lay.units ? ` [${lay.units}]` : "") : " [percentile]"),
                xref: "x" + ax + " domain", yref: "y" + ax + " domain",
                x: 0, y: 1.14, showarrow: false, xanchor: "left",
                font: {size: 12, color: ink("--text-primary")}});

    // An off-mask or offshore site produces a panel with axes and no line, which looks
    // exactly like a rendering bug -- it has already been reported as one. Say why.
    const here = DATA.values.filter(v => v.loc === state.loc && v.lay === pan.lid);
    if (!here.some(v => v[key] !== null)) {
      const statuses = esc([...new Set(here.map(v => v.st))].join(", ")) || "no rows";
      annos.push({text: `No data at this site<br>(${statuses})`,
                  xref: "x" + ax + " domain", yref: "y" + ax + " domain",
                  x: 0.5, y: 0.5, showarrow: false, xanchor: "center", align: "center",
                  font: {size: 12, color: ink("--text-muted")}});
    }
  });

  const layout = baseLayout({
    grid: {rows: rowsN, columns: cols, pattern: "independent",
           roworder: "top to bottom", xgap: 0.14, ygap: 0.34},
    annotations: annos, showlegend: true,
    legend: {orientation: "h", y: -0.14, x: 0, font: {color: ink("--text-secondary")}},
    margin: {l: 56, r: 16, t: 34, b: 62},
    shapes: panels.map((_, i) => ({
      type: "line", xref: "x" + (i === 0 ? "" : i+1), yref: "paper",
      x0: state.decade, x1: state.decade, y0: 0, y1: 1,
      line: {color: ink("--axis"), width: 1},
    })),
  });
  panels.forEach((pan, i) => {
    const ax = i === 0 ? "" : (i + 1);
    layout["xaxis" + ax] = {gridcolor: ink("--grid"), linecolor: ink("--axis"),
                            zeroline: false, tickformat: "d"};
    layout["yaxis" + ax] = {gridcolor: ink("--grid"), linecolor: ink("--axis"),
                            zeroline: false};
    // Score and percentile share the 1-100 risk axis, so pin both for comparability.
    if (pan.kind === "score" || state.seriesMetric === "percentile")
      layout["yaxis" + ax].range = [0, 105];
  });
  Plotly.react("series", traces, layout, CFG);
}

// ---- table --------------------------------------------------------------------------
const COLS = [
  ["Location", r => locById[r.loc].name, false],
  ["Hazard", r => DATA.layers[r.lay].hazard, false],
  ["Scenario", r => r.scn, false],
  ["Decade", r => r.dec, true],
  ["Value", r => r.v, true],
  ["25th", r => r.lo, true],
  ["75th", r => r.hi, true],
  ["Percentile", r => r.p, true],
  ["OLS", r => r.ols, true],
  ["Sen", r => r.sen, true],
  ["Agree", r => r.ag === null ? "—" : String(r.ag), false],
  ["Members", r => r.nm, true],
  ["Status", r => r.st, false],
];

function renderTable() {
  let data = rows({allHazards: true});
  if (state.search) {
    data = data.filter(r => (locById[r.loc].name + " " + DATA.layers[r.lay].hazard + " " +
        r.scn + " " + r.st).toLowerCase().includes(state.search));
  }
  if (state.sortCol !== null) {
    const get = COLS[state.sortCol][1];
    data = data.slice().sort((a,b) => {
      const x = get(a), y = get(b);
      if (x === null) return 1; if (y === null) return -1;
      return (x > y ? 1 : x < y ? -1 : 0) * state.sortDir;
    });
  }
  const head = "<thead><tr>" + COLS.map((c,i) =>
      `<th data-i="${i}">${esc(c[0])}${state.sortCol===i ? (state.sortDir>0?" ▲":" ▼") : ""}</th>`
    ).join("") + "</tr></thead>";
  const body = "<tbody>" + data.map(r => "<tr>" + COLS.map(c => {
      const v = c[1](r);
      const txt = c[2] ? (typeof v === "number" ? fmt(v) : (v ?? "—")) : (v ?? "—");
      return `<td class="${c[2] ? "num" : ""}">${esc(txt)}</td>`;
    }).join("") + "</tr>").join("") + "</tbody>";
  const t = $("table");
  t.innerHTML = head + body;
  t.querySelectorAll("th").forEach(th => th.onclick = () => {
    const i = +th.dataset.i;
    state.sortDir = state.sortCol === i ? -state.sortDir : 1;
    state.sortCol = i;
    renderTable();
  });
}

function renderAll() {
  renderTiles(); renderMap(); renderBars(); renderScoreSeries();
  renderSeries(); renderTable();
}

initControls();
renderAll();
$("foot").textContent =
  "Percentiles are delivered already oriented for risk (100 = worst) — higher_is_better " +
  "layers were inverted at processing time and are NOT re-inverted here. Slopes are per " +
  "decade as stored. Read the slope named in each layer's recommended_slope; disagreement " +
  "between the two is the robustness signal, and there is no p-value under this contract.";
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
            f"forcing pathway, so it is shown on its own and excluded from the Climate "
            f"Score and from cross-tier panels. This is not a missing scenario.")
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
