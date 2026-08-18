"""Shared visualization vocabulary for TCFD dashboards.

Extracted 2026-08-12 because three scripts had independently reinvented the same
decisions and drifted:

    generate_maps.py        symmetric limit at the 95th percentile of |value|
    compare_water_index.py  the same idea at the 98th, written twice more inline
    extract_timber_locations.py  its own scenario -> "Low/Middle/High" recoding

None of those are per-script taste. A diverging scale that is not centred on zero
misreports the sign of a trend, and a scenario that reads "Low Emissions" in one
artifact and "rcp26" in another is the same scenario. They belong in one place.

WHAT THIS MODULE IS NOT: it is not a restatement of the dashboard conventions in the
`/isimip-process-visualize` skill (tab structure, payload limits, per-metric hover formats).
That skill stays authoritative for the gridded-layer dashboards. This module holds the
values and helpers those conventions operate on.

PALETTE PROVENANCE
------------------
The hexes below are copied VERBATIM from the `dataviz` skill's `references/palette.md`,
which is a validated instance -- its own documentation records the measured results:

  * categorical slots 1-3 clear the all-pairs gate in BOTH modes (worst pair CVD dE 9.2
    light / 9.4 dark, normal-vision 24.0 light / 20.9 dark). Three slots is the documented
    all-pairs cap; a fourth puts yellow beside orange and fails. Our scenario dimension has
    exactly three tiers, which is why it fits.
  * the ordinal ramp starts at step 250 (#86b6ef) because the documented light-mode floor
    is "no lighter than step 250" (2.06:1 on the light surface).

DO NOT SUBSTITUTE OTHER HEXES WITHOUT RUNNING `scripts/validate_palette.js` from that
skill. The validation belongs to these specific values, not to the general idea of
"blue, orange, aqua". (No JS runtime was available on 2026-08-12, which is precisely why
the values were taken verbatim rather than re-derived.)
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------------------
# Palette (verbatim from dataviz references/palette.md)
# ---------------------------------------------------------------------------------------

#: Categorical slots, in the fixed order the palette validated. Assign by ENTITY, never by
#: rank -- a filter that drops a series must not repaint the survivors.
CATEGORICAL_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a"]
CATEGORICAL_DARK = ["#3987e5", "#d95926", "#199e70"]

#: FORCING-TIER colours: blue / yellow / red for low / medium / high (user decision
#: 2026-08-12), the conventional climate-scenario reading. These are palette slots 1, 4 and
#: 8 -- NOT the validated adjacent triple (1,2,3), so the documented all-pairs result does
#: not transfer and could not be re-measured here (no JS runtime).
#:
#: What WAS measured: blue and red sit at near-identical WCAG luminance (1.12:1 apart on
#: light, 1.05:1 on dark), so there is effectively NO lightness fallback if hue fails; and
#: light-mode yellow is 2.11:1 on the surface, below the 3:1 bar, which the reference
#: palette already flags for this hue. Yellow and red also converge under deuteranopia.
#:
#: Therefore secondary encoding is applied UNCONDITIONALLY rather than treated as optional:
#: every tier carries a distinct marker symbol and line dash in addition to its hue, the
#: legend is always present, and the table view carries every value. Identity is never
#: colour-alone here.
TIER_COLOR_LIGHT = {"low": "#2a78d6", "medium": "#eda100", "high": "#e34948"}
TIER_COLOR_DARK = {"low": "#3987e5", "medium": "#c98500", "high": "#e66767"}

#: The secondary channels. Ordered low -> high so the encoding itself reads as a progression.
TIER_SYMBOL = {"low": "circle", "medium": "square", "high": "diamond"}
TIER_DASH = {"low": "solid", "medium": "dot", "high": "dash"}

#: Sequential single-hue blue ramp, light -> dark. Used for RAW MAGNITUDE (`value`), where
#: large is not inherently bad -- keeping it blue also distinguishes "how much" from "how
#: risky" at a glance.
SEQUENTIAL_BLUE = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]

#: Sequential single-hue RED ramp for RISK magnitude -- `percentile` and Climate Score,
#: where 100 is worst by construction and red carries that meaning without a legend lookup.
#: (User decision 2026-08-12.)
#:
#: NOT from the dataviz reference palette, which ships only a blue sequential ramp. Built
#: here and verified with what can be checked WITHOUT the JS validator, which is unavailable
#: on this machine: WCAG relative luminance decreases monotonically across every step, and
#: each step's contrast against the surface it renders on was measured. A single-hue
#: light->dark ramp is the sequential rule, so the CVD pair checks that govern CATEGORICAL
#: palettes do not apply -- but if a validator ever becomes available, re-run it here.
#: Anchored on the palette's own red (#e34948) and status-critical (#d03b3b).
SEQUENTIAL_RED = [
    "#fde3e0", "#fbcfcb", "#f9b8b3", "#f59d99", "#f0837f", "#e96562",
    "#e34948", "#d03b3b", "#bb3030", "#a52727", "#8c1f1f", "#741919", "#5c1313",
]

#: Dark-surface variant: the same ramp truncated before the near-black end, which would
#: otherwise vanish into #1a1a19. Darkest step measured 2.42:1 on the dark surface.
SEQUENTIAL_RED_DARK = [
    "#fde3e0", "#fbcfcb", "#f9b8b3", "#f59d99", "#f0837f", "#e96562",
    "#e34948", "#d03b3b", "#bb3030", "#a52727",
]

#: Ordinal ramp for ORDERED categories (the five risk bands). Red for the same reason as
#: above. Every step clears the documented 2:1 floor on the surface it renders on --
#: lightest light-mode step 2.01, darkest dark-mode step 2.97 -- which is why the two modes
#: use different steps rather than one list.
ORDINAL_RISK = ["#f59d99", "#e96562", "#d03b3b", "#a52727", "#741919"]
ORDINAL_RISK_DARK = ["#fbcfcb", "#f59d99", "#e96562", "#d03b3b", "#bb3030"]

#: Diverging poles for signed quantities (slopes). Warm/cool so they read as opposite,
#: with a NEUTRAL GRAY midpoint -- a hue at the midpoint would stop it reading as "nothing".
DIVERGING_LOW = "#2a78d6"      # blue
DIVERGING_MID_LIGHT = "#f0efec"
DIVERGING_MID_DARK = "#383835"
DIVERGING_HIGH = "#e34948"     # red

INK = {
    "light": {
        "surface": "#fcfcfb", "plane": "#f9f9f7", "primary": "#0b0b0b",
        "secondary": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
        "axis": "#c3c2b7", "border": "rgba(11,11,11,0.10)",
    },
    "dark": {
        "surface": "#1a1a19", "plane": "#0d0d0d", "primary": "#ffffff",
        "secondary": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
        "axis": "#383835", "border": "rgba(255,255,255,0.10)",
    },
}

FONT_STACK = 'system-ui, -apple-system, "Segoe UI", sans-serif'


def plotly_colorscale(hexes: Sequence[str]) -> List[List]:
    """Even-spaced Plotly colorscale from an ordered hex list."""
    n = len(hexes)
    return [[i / (n - 1), h] for i, h in enumerate(hexes)]


def diverging_colorscale(mode: str = "light") -> List[List]:
    mid = DIVERGING_MID_LIGHT if mode == "light" else DIVERGING_MID_DARK
    return [[0.0, DIVERGING_LOW], [0.5, mid], [1.0, DIVERGING_HIGH]]


# ---------------------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------------------

#: ISIMIP2b/RCP and ISIMIP3b/SSP scenarios mapped onto a shared forcing TIER.
#:
#: The delivery CSV carries native codes only, by user decision -- harmonization is the
#: report's job, and this module is the report side. The tier exists so a mixed-round
#: portfolio (cyclone is RCP-only, wildfire is SSP-only) can share one colour legend and one
#: filter. Colour is assigned by TIER so `rcp26` and `ssp126` are the same blue; the NATIVE
#: CODE is what gets displayed. RCP and SSP tiers are only approximately comparable and any
#: narrative must say so.
SCENARIO_TIER: Dict[str, str] = {
    "rcp26": "low", "ssp126": "low",
    "rcp45": "medium", "rcp60": "medium", "ssp245": "medium", "ssp370": "medium",
    "rcp85": "high", "ssp585": "high",
}

#: Scenario codes that are NOT forcing pathways at all. An observational layer summarises
#: a measured historical window; it has no radiative forcing, so it belongs to no tier and
#: must never be defaulted into one. Folding `observed` into "medium" would make the medium
#: pathway look worse than low and high for a reason that has nothing to do with forcing.
NON_FORCING_SCENARIOS = frozenset({"observed"})


def is_forcing_scenario(scenario: str) -> bool:
    """False for an observational scenario code, which has no forcing tier."""
    return str(scenario).lower() not in NON_FORCING_SCENARIOS


TIER_ORDER = ["low", "medium", "high"]
TIER_LABELS = {"low": "Low forcing", "medium": "Medium forcing", "high": "High forcing"}


def tier_of(scenario: str) -> str:
    """Forcing tier for a scenario code. Unknown codes fall to 'medium' but are reported.

    Check `is_forcing_scenario()` first for anything that aggregates across tiers: a
    non-forcing code has no meaningful answer here and the 'medium' fallback is a
    placeholder, not a classification.
    """
    return SCENARIO_TIER.get(scenario.lower(), "medium")


def tier_color(tier: str, mode: str = "light") -> str:
    table = TIER_COLOR_LIGHT if mode == "light" else TIER_COLOR_DARK
    return table.get(tier, table["medium"])


def check_tier_collisions(scenarios_by_layer: Dict[str, Sequence[str]]) -> List[str]:
    """Report layers where two scenarios collapse onto one tier.

    `rcp45` and `rcp60` both map to 'medium'. No shipped layer publishes both, but if one
    ever does, colouring by tier would draw two different scenarios in the same hue with no
    way to tell them apart. Surface it rather than rendering a quietly wrong legend.
    """
    problems = []
    for layer_id, scenarios in scenarios_by_layer.items():
        seen: Dict[str, List[str]] = {}
        for s in scenarios:
            seen.setdefault(tier_of(s), []).append(s)
        for tier, group in seen.items():
            if len(group) > 1:
                problems.append(
                    f"{layer_id}: {', '.join(group)} all map to tier '{tier}' -- "
                    f"they would render in one colour"
                )
    return problems


# ---------------------------------------------------------------------------------------
# Risk bands (display-only)
# ---------------------------------------------------------------------------------------

#: Percentile -> ordered band. DISPLAY-ONLY: it is derived from `percentile` and therefore
#: deliberately absent from the delivery CSV under the no-derived-columns rule. It lives
#: here because a rollup chart needs ordered classes; it is not a stored value class.
#:
#: Thresholds carried over from the retired `export_formatter.RELATIVE_HAZARD_THRESHOLDS`
#: so a band means the same thing it did in prior deliveries.
RISK_BANDS = [
    (0, 20, "Very Low"),
    (20, 40, "Low"),
    (40, 60, "Medium"),
    (60, 80, "High"),
    (80, 101, "Very High"),
]
RISK_BAND_ORDER = [b[2] for b in RISK_BANDS]


def risk_band(percentile: Optional[float]) -> Optional[str]:
    if percentile is None or not np.isfinite(percentile):
        return None
    for lo, hi, label in RISK_BANDS:
        if percentile <= hi and (percentile > lo or lo == 0):
            return label
    return RISK_BAND_ORDER[-1]


# ---------------------------------------------------------------------------------------
# Scale helpers
# ---------------------------------------------------------------------------------------

#: Percentile of |value| used for a diverging limit on GRIDDED data. On a 0.5 deg global
#: grid a handful of outliers set max|value| and leave ~99.7% of cells near-white, so the
#: body of the distribution needs the colour range (measured on the wildfire layer:
#: ols_slope limit 38.2 true-max vs 1.16 at the 95th).
GRID_SYMMETRIC_PERCENTILE = 95

#: Site portfolios are hundreds of points, not 70k, and have no comparable tail -- so a
#: percentile limit would clamp real sites for no readability gain. Below this count the
#: limit is true max and NOTHING is clamped.
POINT_SYMMETRIC_MIN_N = 40


def symmetric_limit(
    values: Iterable[float],
    percentile: Optional[float] = GRID_SYMMETRIC_PERCENTILE,
) -> Tuple[float, int]:
    """Limit L for a diverging scale centred on zero, plus the count clamped beyond it.

    The scale is always `[-L, +L]` so zero sits at the neutral midpoint. Values beyond L are
    clamped by Plotly to the endpoint colour, never blanked -- but the clamped COUNT must be
    reported alongside the chart, which is why it is returned rather than discarded.
    """
    v = np.asarray(list(values), dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 1.0, 0
    mag = np.abs(v)
    limit = float(np.max(mag)) if percentile is None else float(np.percentile(mag, percentile))
    if limit <= 0:
        limit = float(np.max(mag)) or 1.0
    return limit, int((mag > limit).sum())


def point_symmetric_limit(values: Iterable[float]) -> Tuple[float, int, str]:
    """Diverging limit for a SITE portfolio, with the rule used stated for the subtitle."""
    v = [x for x in values if x is not None and np.isfinite(x)]
    if len(v) >= POINT_SYMMETRIC_MIN_N:
        limit, clamped = symmetric_limit(v, GRID_SYMMETRIC_PERCENTILE)
        return limit, clamped, f"{GRID_SYMMETRIC_PERCENTILE}th pct of |value|"
    limit, clamped = symmetric_limit(v, None)
    return limit, clamped, "max |value| (no clamping)"


def sigfig(a: np.ndarray, n: int = 4) -> np.ndarray:
    """Round to n significant figures. Payload control -- finer than a colour can encode."""
    out = np.array(a, dtype=np.float64, copy=True)
    nz = np.isfinite(out) & (out != 0)
    if np.any(nz):
        mag = np.floor(np.log10(np.abs(out[nz])))
        out[nz] = np.round(out[nz] * np.power(10.0, n - 1 - mag)) / np.power(10.0, n - 1 - mag)
    return out


# ---------------------------------------------------------------------------------------
# HTML chrome
# ---------------------------------------------------------------------------------------

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.27.0.min.js"


def css_tokens() -> str:
    """Theme-aware CSS custom properties.

    Light values on bare `:root`; dark redefined under BOTH the OS media query and an
    explicit `[data-theme="dark"]` stamp, so a viewer's toggle wins in either direction.
    """
    def block(mode: str) -> str:
        ink = INK[mode]
        tiers = TIER_COLOR_LIGHT if mode == "light" else TIER_COLOR_DARK
        # Emit the TIER tokens the page actually paints with, not the generic categorical
        # slots. Carrying unused `--series-2/3` left the retired orange and aqua sitting in
        # the stylesheet of a page that no longer draws with them -- dead colour tokens are
        # what a future edit picks up by accident.
        return "\n".join(
            [f"  --surface: {ink['surface']};",
             f"  --plane: {ink['plane']};",
             f"  --text-primary: {ink['primary']};",
             f"  --text-secondary: {ink['secondary']};",
             f"  --text-muted: {ink['muted']};",
             f"  --grid: {ink['grid']};",
             f"  --axis: {ink['axis']};",
             f"  --border: {ink['border']};"]
            + [f"  --tier-{t}: {tiers[t]};" for t in TIER_ORDER]
            # `--series-1` is the UI accent (active toggle); it tracks the low tier so the
            # chrome and the charts stay one system.
            + [f"  --series-1: {tiers['low']};"]
        )

    return f""":root {{
  color-scheme: light;
{block('light')}
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    color-scheme: dark;
{block('dark')}
  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
{block('dark')}
}}"""
