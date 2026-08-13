"""Inline-SVG figure vocabulary for customer reports. No JavaScript, no plotly, no images.

WHY SVG AND NOT PLOTLY
----------------------
`generate_delivery_dashboard.py` is an INSTRUMENT -- it exists to be poked at, so an
interactive engine earns its 3.5 MB. A report is a DOCUMENT: it gets read once, printed,
attached to a filing, and re-read by an auditor two years later. Those are opposite
requirements, and this machine settles the argument anyway:

    no pandoc, no weasyprint, no wkhtmltopdf, no kaleido, no node, no headless Chrome
    (measured 2026-08-13; only Safari and cupsfilter are present)

so a static raster or a server-rendered PDF is not producible here at all. Inline SVG is:
~10 KB per figure, prints at printer resolution rather than screen resolution, survives with
JavaScript disabled or stripped by an email gateway, and is diffable in git -- which matters
because a figure in a disclosure document is an assertion, and an assertion should be
reviewable line by line.

The trade is that this module has to draw its own axes. That is the whole cost, and it is
bounded: a report needs four figure types, not forty.

COLOUR IS NOT CHOSEN HERE
-------------------------
Every colour is a CSS custom property emitted by `report_common.css_tokens()`, which derives
them from `viz_common`. This module never writes a hex literal. Two consequences worth
stating:

  - light/dark mode is handled by the page redefining the tokens; one SVG serves both, and
    print forces the light set.
  - continuous risk is QUANTIZED to the five `viz_common.RISK_BANDS` classes rather than
    drawn from the 13-step ramp. In an interactive chart a reader hovers to recover the
    exact value; in a static one they cannot, so a shade they must map back to a number by
    eye is a shade that will be misread. Five labelled classes are honest about the
    precision a printed swatch actually carries.

ACCESSIBILITY IS STRUCTURAL, NOT OPTIONAL
-----------------------------------------
Forcing tiers are blue/yellow/red, which was a user decision, and the measurement behind
`viz_common.TIER_COLOR_LIGHT` is that blue and red sit 1.12:1 apart in luminance -- so there
is effectively NO lightness fallback if hue perception fails. Every tier series therefore
carries a distinct marker shape AND a distinct dash pattern, unconditionally. Do not strip
them as clutter; they are the channel that survives deuteranopia and greyscale printing.

Every figure emits `role="img"` with a `<title>` and a `<desc>`, and every figure in a report
is paired with a data table by `report_common` -- the figure is never the only way to reach a
number.
"""

from __future__ import annotations

import hashlib
import html
from typing import Dict, List, Optional, Sequence, Tuple

from .viz_common import (
    RISK_BANDS,
    RISK_BAND_ORDER,
    TIER_DASH,
    TIER_LABELS,
    TIER_ORDER,
    TIER_SYMBOL,
    risk_band,
)


# ---------------------------------------------------------------------------------------
# Mark specs (dataviz references/marks-and-anatomy.md)
# ---------------------------------------------------------------------------------------

#: 2px surface-coloured gap between adjacent fills, so a stacked bar reads as segments
#: rather than as one smeared block. Applied as a stroke in the surface colour.
SEGMENT_GAP = 2
#: Rounded data-end on a bar. The end away from the baseline is rounded; the baseline end
#: stays square so the bar visibly starts at zero.
BAR_RADIUS = 4
LINE_WIDTH = 2
#: Markers are >= 8px so they remain hit-able and, more importantly here, so the SHAPE is
#: legible at print size -- an undersized diamond and an undersized square are the same
#: smudge, which would destroy the secondary encoding.
MARKER_SIZE = 9
GRID_WIDTH = 1

#: SVG user-space width. Height varies per figure. The page scales this with
#: `width: 100%; height: auto`, so these are proportions, not pixels.
FIG_WIDTH = 720

#: Dash arrays matching viz_common.TIER_DASH, in user units.
DASH_ARRAY = {"solid": "none", "dot": "2 4", "dash": "8 5"}


def _esc(text) -> str:
    """XML-escape. SVG is XML: a site named `Depot & Sons <NC>` breaks the document."""
    return html.escape(str(text), quote=True)


def fmt(value: Optional[float], places: int = 1) -> str:
    """Format a number for a label. Bare '-' for missing, never 'nan' or 'None'."""
    if value is None:
        return "–"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "–"
    if f != f:  # NaN
        return "–"
    if abs(f) >= 1000:
        return f"{f:,.0f}"
    if abs(f) >= 100 or places == 0:
        return f"{f:.0f}"
    return f"{f:.{places}f}"


def band_index(percentile: Optional[float]) -> Optional[int]:
    """0-4 index into RISK_BANDS, or None. The quantization the module docstring describes."""
    label = risk_band(percentile)
    if label is None:
        return None
    return RISK_BAND_ORDER.index(label)


def band_var(percentile: Optional[float]) -> str:
    idx = band_index(percentile)
    return "var(--ink-axis)" if idx is None else f"var(--risk-{idx + 1})"


#: Average glyph advance as a fraction of font-size for the UI stack, used to wrap labels.
#: An estimate is the only option -- SVG cannot measure text without a layout engine, and
#: there is no browser here. Erring narrow costs a little whitespace; erring wide clips the
#: label, which is the failure this replaces.
CHAR_WIDTH_EM = 0.55

SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _wrap(text: str, width_px: float, font_size: float, max_lines: int = 2) -> List[str]:
    """Wrap a label to fit a fixed column, breaking on spaces, ellipsising past max_lines.

    Site labels routinely run past 50 characters ("Gulf Platform Alpha — Warehouse (Offshore
    supply base)") and a fixed single-line column silently truncates them at the plot edge,
    which loses exactly the part that distinguishes one asset at a site from another.
    """
    max_chars = max(8, int(width_px / (font_size * CHAR_WIDTH_EM)))
    words, lines, cur = str(text).split(), [], ""
    for w in words:
        trial = f"{cur} {w}".strip()
        if len(trial) <= max_chars or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = w
            if len(lines) == max_lines:
                break
    if cur and len(lines) < max_lines:
        lines.append(cur)
    if not lines:
        return [""]
    consumed = len(" ".join(lines).split())
    if consumed < len(words):
        tail = lines[-1]
        lines[-1] = (tail[: max_chars - 1] + "…") if len(tail) >= max_chars - 1 else tail + " …"
    return lines


def _scale_for(limit: float) -> Tuple[float, str]:
    """Pick a display multiplier so small slopes stop rendering as 0.000.

    Returns (multiplier, suffix). Slopes here span three orders of magnitude between layers
    -- cyclone tops out around 9e-4 per decade while wildfire reaches 3e-1 -- so a fixed
    number of decimal places renders whole figures as columns of zeros. Each figure picks
    its own factor and states it in the axis label, which is the same convention a scientific
    axis uses.
    """
    if not limit or limit != limit or limit <= 0:
        return 1.0, ""
    import math

    if 0.1 <= limit < 10000:
        return 1.0, ""
    exp = -int(math.floor(math.log10(limit)))
    return 10.0 ** exp, f" ×10{str(-exp).translate(SUPERSCRIPT)}"


def _marker(shape: str, cx: float, cy: float, size: float, fill: str) -> str:
    """One marker in the tier's assigned SHAPE. A 2px surface ring keeps overlapping
    markers separable where two tiers cross."""
    r = size / 2
    ring = f'stroke="var(--ink-surface)" stroke-width="2"'
    if shape == "circle":
        return f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{fill}" {ring}/>'
    if shape == "square":
        return (
            f'<rect x="{cx - r:.2f}" y="{cy - r:.2f}" width="{size:.2f}" '
            f'height="{size:.2f}" fill="{fill}" {ring}/>'
        )
    # diamond
    pts = f"{cx:.2f},{cy - r:.2f} {cx + r:.2f},{cy:.2f} {cx:.2f},{cy + r:.2f} {cx - r:.2f},{cy:.2f}"
    return f'<polygon points="{pts}" fill="{fill}" {ring}/>'


def _frame(width: int, height: int, title: str, desc: str, body: str) -> str:
    return (
        f'<svg viewBox="0 0 {width} {height}" width="100%" height="auto" '
        f'role="img" aria-label="{_esc(title)}" '
        f'xmlns="http://www.w3.org/2000/svg" class="fig">'
        f"<title>{_esc(title)}</title><desc>{_esc(desc)}</desc>"
        f"{body}</svg>"
    )


def _text(x, y, s, *, anchor="start", cls="fig-label", size=12, weight=None) -> str:
    w = f' font-weight="{weight}"' if weight else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" class="{cls}" '
        f'font-size="{size}"{w}>{_esc(s)}</text>'
    )


# ---------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------


def ranked_bar(
    rows: Sequence[Tuple[str, Optional[float]]],
    *,
    title: str,
    axis_max: float = 100,
    axis_label: str = "Risk percentile (100 = worst)",
    label_width: int = 210,
    colour_by_band: bool = True,
) -> str:
    """Horizontal bars, one per row, ordered as given. The workhorse ranking figure.

    Every bar is directly labelled with its value -- there are at most a few dozen and the
    reader needs the number, not an impression of it. Colour encodes the risk BAND, so the
    bar length and the fill say the same thing twice; that redundancy is deliberate for
    greyscale printing, where only the length survives.

    A row with no value draws no bar and is annotated, rather than drawing a zero-length bar
    that reads as "assessed, and found to be nil".
    """
    label_size = 11.5
    wrapped = [_wrap(label, label_width - 8, label_size) for label, _ in rows]
    n_lines = max((len(w) for w in wrapped), default=1)
    line_h = 13.5
    row_h = max(26, n_lines * line_h + 8)
    gap, pad_top, pad_bottom = 7, 34, 30
    plot_left = label_width
    plot_right = FIG_WIDTH - 58
    plot_w = plot_right - plot_left
    height = pad_top + len(rows) * (row_h + gap) + pad_bottom

    out: List[str] = []
    out.append(_text(0, 16, title, cls="fig-title", size=13, weight="600"))

    # Recessive gridlines at 0/20/40/60/80/100, drawn behind the bars.
    for frac in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
        x = plot_left + frac * plot_w
        out.append(
            f'<line x1="{x:.2f}" y1="{pad_top - 6}" x2="{x:.2f}" '
            f'y2="{height - pad_bottom + 4}" stroke="var(--ink-grid)" '
            f'stroke-width="{GRID_WIDTH}"/>'
        )
        out.append(
            _text(x, height - pad_bottom + 18, fmt(frac * axis_max, 0),
                  anchor="middle", cls="fig-tick", size=11)
        )

    for i, (label, value) in enumerate(rows):
        y = pad_top + i * (row_h + gap)
        lines = wrapped[i]
        # Vertically centre the wrapped block against the bar.
        first = y + row_h / 2 - (len(lines) - 1) * line_h / 2 + 4
        for k, line in enumerate(lines):
            out.append(_text(0, first + k * line_h, line, cls="fig-label", size=label_size))
        if value is None or value != value:
            out.append(
                _text(plot_left + 4, y + row_h * 0.62, "not assessed at this site",
                      cls="fig-null", size=11)
            )
            continue
        frac = max(0.0, min(1.0, float(value) / axis_max))
        w = max(frac * plot_w, 2.0)
        fill = band_var(value) if colour_by_band else "var(--tier-medium)"
        # Rounded data-end only: rx on a <rect> rounds all four corners, so the baseline
        # end is squared off again with a small overlay. A bar that is rounded at the
        # baseline looks detached from zero.
        out.append(
            f'<rect x="{plot_left:.2f}" y="{y}" width="{w:.2f}" height="{row_h}" '
            f'rx="{BAR_RADIUS}" fill="{fill}"/>'
        )
        out.append(
            f'<rect x="{plot_left:.2f}" y="{y}" width="{min(BAR_RADIUS, w):.2f}" '
            f'height="{row_h}" fill="{fill}"/>'
        )
        out.append(
            _text(plot_left + w + 6, y + row_h * 0.62, fmt(value), cls="fig-value", size=12)
        )

    out.append(
        _text(plot_left, height - 4, axis_label, cls="fig-axis-title", size=11)
    )
    desc = "; ".join(f"{lab} {fmt(v)}" for lab, v in rows)
    return _frame(FIG_WIDTH, height, title, f"{axis_label}. {desc}", "".join(out))


def decade_line(
    series: Dict[str, Sequence[Tuple[int, Optional[float]]]],
    *,
    title: str,
    y_label: str,
    y_min: Optional[float] = None,
    y_max: Optional[float] = None,
    tier_keyed: bool = True,
) -> str:
    """One line per forcing tier across decades.

    `series` maps a tier key (`low`/`medium`/`high`) to (decade, value) pairs. Tier keying
    is what lets a mixed-round portfolio share one legend: `rcp26` and `ssp126` are one blue
    line, and the native codes are named in the caption rather than the legend.

    Missing points BREAK the line rather than interpolating across the gap. An interpolated
    segment through a decade a layer does not publish is a fabricated number in a figure.
    """
    pad_l, pad_r, pad_t, pad_b = 52, 96, 34, 44
    height = 300
    plot_w = FIG_WIDTH - pad_l - pad_r
    plot_h = height - pad_t - pad_b

    pts = [(d, v) for s in series.values() for d, v in s if v is not None and v == v]
    if not pts:
        return empty_figure(title, "No values available for this selection.")
    decades = sorted({d for d, _ in pts})
    vals = [v for _, v in pts]
    lo = y_min if y_min is not None else min(vals)
    hi = y_max if y_max is not None else max(vals)
    if hi == lo:
        hi, lo = hi + 1, lo - 1
    span_x = max(decades) - min(decades) or 1

    def X(d):
        return pad_l + (d - min(decades)) / span_x * plot_w

    def Y(v):
        return pad_t + plot_h - (float(v) - lo) / (hi - lo) * plot_h

    out: List[str] = [_text(0, 16, title, cls="fig-title", size=13, weight="600")]

    for k in range(5):
        frac = k / 4
        y = pad_t + plot_h - frac * plot_h
        out.append(
            f'<line x1="{pad_l}" y1="{y:.2f}" x2="{pad_l + plot_w}" y2="{y:.2f}" '
            f'stroke="var(--ink-grid)" stroke-width="{GRID_WIDTH}"/>'
        )
        out.append(
            _text(pad_l - 8, y + 4, fmt(lo + frac * (hi - lo)), anchor="end",
                  cls="fig-tick", size=11)
        )
    for d in decades:
        out.append(
            _text(X(d), pad_t + plot_h + 18, f"{d}s", anchor="middle", cls="fig-tick", size=11)
        )

    legend_y = pad_t + 4
    for tier in TIER_ORDER:
        if tier not in series:
            continue
        colour = f"var(--tier-{tier})" if tier_keyed else "var(--tier-medium)"
        dash = DASH_ARRAY[TIER_DASH[tier]]
        shape = TIER_SYMBOL[tier]
        ordered = sorted(series[tier], key=lambda p: p[0])

        # Break the path at every gap.
        run: List[str] = []
        for d, v in ordered:
            if v is None or v != v:
                if len(run) > 1:
                    out.append(
                        f'<polyline points="{" ".join(run)}" fill="none" stroke="{colour}" '
                        f'stroke-width="{LINE_WIDTH}" stroke-dasharray="{dash}" '
                        f'stroke-linecap="round" stroke-linejoin="round"/>'
                    )
                run = []
                continue
            run.append(f"{X(d):.2f},{Y(v):.2f}")
        if len(run) > 1:
            out.append(
                f'<polyline points="{" ".join(run)}" fill="none" stroke="{colour}" '
                f'stroke-width="{LINE_WIDTH}" stroke-dasharray="{dash}" '
                f'stroke-linecap="round" stroke-linejoin="round"/>'
            )
        for d, v in ordered:
            if v is not None and v == v:
                out.append(_marker(shape, X(d), Y(v), MARKER_SIZE, colour))

        # Direct label at the last finite point, so the legend is a convenience rather
        # than the only route to identity.
        finite = [(d, v) for d, v in ordered if v is not None and v == v]
        if finite:
            d, v = finite[-1]
            out.append(
                _text(X(d) + 10, Y(v) + 4, TIER_LABELS[tier], cls="fig-value", size=11)
            )
        out.append(_marker(shape, FIG_WIDTH - pad_r + 14, legend_y, MARKER_SIZE, colour))
        out.append(
            _text(FIG_WIDTH - pad_r + 24, legend_y + 4, TIER_LABELS[tier],
                  cls="fig-label", size=11)
        )
        legend_y += 18

    out.append(_text(0, height - 6, y_label, cls="fig-axis-title", size=11))
    desc = "; ".join(
        f"{TIER_LABELS[t]}: " + ", ".join(f"{d}s {fmt(v)}" for d, v in sorted(s))
        for t, s in series.items()
    )
    return _frame(FIG_WIDTH, height, title, f"{y_label}. {desc}", "".join(out))


def band_stack(
    counts: Dict[str, int],
    *,
    title: str,
    total_label: str = "assets",
) -> str:
    """One stacked bar showing how a population distributes across the five risk bands.

    Bands are ordered Very Low -> Very High always, never by size: the axis is ordinal and
    re-sorting it by count would destroy the reading. Segments carry a 2px surface gap.
    """
    total = sum(counts.values())
    if not total:
        return empty_figure(title, "No assets in this selection.")
    height, bar_y, bar_h = 118, 30, 34
    out: List[str] = [_text(0, 16, title, cls="fig-title", size=13, weight="600")]

    x = 0.0
    for i, label in enumerate(RISK_BAND_ORDER):
        n = counts.get(label, 0)
        if not n:
            continue
        w = n / total * FIG_WIDTH
        out.append(
            f'<rect x="{x:.2f}" y="{bar_y}" width="{max(w - SEGMENT_GAP, 1):.2f}" '
            f'height="{bar_h}" rx="{BAR_RADIUS}" fill="var(--risk-{i + 1})"/>'
        )
        # Label inside when it fits, below when it does not -- never omitted.
        if w > 62:
            out.append(
                _text(x + w / 2 - SEGMENT_GAP / 2, bar_y + bar_h * 0.62, str(n),
                      anchor="middle", cls="fig-value-on-fill", size=12, weight="600")
            )
        out.append(
            _text(x, bar_y + bar_h + 18, f"{label} ({n})", cls="fig-tick", size=11)
            if w > 96
            else ""
        )
        x += w

    out.append(
        _text(0, height - 8, f"{total} {total_label} · bands are ordinal, left to right",
              cls="fig-axis-title", size=11)
    )
    desc = ", ".join(f"{lab} {counts.get(lab, 0)}" for lab in RISK_BAND_ORDER)
    return _frame(FIG_WIDTH, height, title, f"{total} {total_label}: {desc}", "".join(out))


def trend_strip(
    rows: Sequence[Tuple[str, Optional[float], Optional[bool]]],
    *,
    title: str,
    units: str,
) -> str:
    """Signed slopes on a zero-centred diverging axis, one row per item.

    Two things this figure does that a plain bar chart would not:

      - the axis is SYMMETRIC about zero with a neutral midpoint, so a small positive and a
        small negative slope look equally small (GUARDRAILS §5).
      - `slopes_agree` is drawn as a hatched overlay rather than dropped. Under this
        contract there is no p-value; disagreement between the OLS and Theil-Sen estimators
        IS the robustness signal, and a bar with no such marking is asserting a trend the
        data does not support.
    """
    finite = [v for _, v, _ in rows if v is not None and v == v]
    if not finite:
        return empty_figure(title, "No slope values available for this selection.")
    limit = max(abs(v) for v in finite) or 1.0
    scale, scale_suffix = _scale_for(limit)

    label_width, label_size, line_h = 236, 11.5, 13.5
    wrapped = [_wrap(label, label_width - 8, label_size) for label, _v, _a in rows]
    n_lines = max((len(w) for w in wrapped), default=1)
    row_h = max(24, n_lines * line_h + 6)
    gap, pad_top, pad_bottom = 7, 34, 42
    plot_left, plot_right = label_width, FIG_WIDTH - 76
    mid = (plot_left + plot_right) / 2
    half = (plot_right - plot_left) / 2
    height = pad_top + len(rows) * (row_h + gap) + pad_bottom

    # SVG ids are document-global even inside inline <svg>, so two trend strips on one page
    # would declare the same id twice and every url(#...) would resolve to the first. The id
    # is therefore derived from the figure's own content: deterministic (the report build
    # stamp depends on it) and unique per distinct figure.
    pat = "hatch-" + hashlib.sha256(
        (title + units + repr(list(rows))).encode()
    ).hexdigest()[:8]

    out: List[str] = [
        _text(0, 16, title, cls="fig-title", size=13, weight="600"),
        f'<defs><pattern id="{pat}" width="6" height="6" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)">'
        '<rect width="6" height="6" fill="var(--ink-surface)"/>'
        '<line x1="0" y1="0" x2="0" y2="6" stroke="var(--ink-muted)" stroke-width="3"/>'
        "</pattern></defs>",
    ]
    out.append(
        f'<line x1="{mid:.2f}" y1="{pad_top - 6}" x2="{mid:.2f}" '
        f'y2="{height - pad_bottom + 4}" stroke="var(--ink-axis)" stroke-width="{GRID_WIDTH}"/>'
    )

    for i, (label, value, agree) in enumerate(rows):
        y = pad_top + i * (row_h + gap)
        lines = wrapped[i]
        first = y + row_h / 2 - (len(lines) - 1) * line_h / 2 + 4
        for k, line in enumerate(lines):
            out.append(_text(0, first + k * line_h, line, cls="fig-label", size=label_size))
        if value is None or value != value:
            out.append(_text(mid + 6, y + row_h * 0.65, "no trend published",
                             cls="fig-null", size=11))
            continue
        w = abs(float(value)) / limit * half
        x = mid if value >= 0 else mid - w
        colour = "var(--div-high)" if value >= 0 else "var(--div-low)"
        out.append(
            f'<rect x="{x:.2f}" y="{y}" width="{max(w, 1.5):.2f}" height="{row_h}" '
            f'rx="2" fill="{colour}"/>'
        )
        if agree is False:
            out.append(
                f'<rect x="{x:.2f}" y="{y}" width="{max(w, 1.5):.2f}" height="{row_h}" '
                f'rx="2" fill="url(#{pat})" opacity="0.55"/>'
            )
        anchor_x = mid + w + 6 if value >= 0 else mid - w - 6
        out.append(
            _text(anchor_x, y + row_h * 0.65, fmt(float(value) * scale, 2),
                  anchor="start" if value >= 0 else "end", cls="fig-value", size=11)
        )

    out.append(
        _text(0, height - 22,
              f"Slope in {units}{scale_suffix}; axis symmetric about zero at "
              f"±{fmt(limit * scale, 2)}",
              cls="fig-axis-title", size=11)
    )
    out.append(
        _text(0, height - 6,
              "Hatched = the two slope estimators disagree; treat that trend as not robust",
              cls="fig-axis-title", size=11)
    )
    desc = "; ".join(
        f"{lab} {fmt((v * scale) if v is not None and v == v else v, 2)}"
        f"{'' if a is not False else ' (not robust)'}"
        for lab, v, a in rows
    )
    return _frame(
        FIG_WIDTH, height, title, f"Slopes in {units}{scale_suffix}. {desc}", "".join(out)
    )


def empty_figure(title: str, reason: str) -> str:
    """A figure that cannot be drawn says why.

    An empty panel is indistinguishable from a rendering bug -- that ambiguity has already
    produced bug reports against the dashboard. A report is worse: nobody is there to ask.
    """
    height = 96
    body = (
        _text(0, 16, title, cls="fig-title", size=13, weight="600")
        + f'<rect x="0" y="26" width="{FIG_WIDTH}" height="52" fill="var(--ink-plane)" '
          f'stroke="var(--ink-grid)" stroke-width="1" rx="4"/>'
        + _text(FIG_WIDTH / 2, 57, reason, anchor="middle", cls="fig-null", size=12)
    )
    return _frame(FIG_WIDTH, height, title, reason, body)


def band_legend() -> str:
    """The five ordinal risk classes, with their percentile ranges. Always shown wherever a
    band colour is used, because a colour whose thresholds are not stated is not a
    measurement."""
    swatch, height = 16, 30
    out: List[str] = []
    x = 0.0
    for i, (lo, hi, label) in enumerate(RISK_BANDS):
        upper = 100 if hi > 100 else hi
        out.append(
            f'<rect x="{x:.2f}" y="6" width="{swatch}" height="{swatch}" rx="3" '
            f'fill="var(--risk-{i + 1})"/>'
        )
        out.append(_text(x + swatch + 6, 19, f"{label} ({lo}–{upper})",
                         cls="fig-label", size=11))
        x += swatch + 12 + 8.0 * len(f"{label} ({lo}-{upper})")
    return _frame(
        FIG_WIDTH, height, "Risk band legend",
        "Percentile bands: " + ", ".join(f"{l} {lo} to {min(h,100)}" for lo, h, l in RISK_BANDS),
        "".join(out),
    )
