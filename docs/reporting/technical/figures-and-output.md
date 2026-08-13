# Figures and output mechanics

## The environment decides the format

Measured on this machine, 2026-08-13:

| Tool | Present |
|---|---|
| pandoc, weasyprint, wkhtmltopdf, prince | no |
| kaleido, matplotlib, jinja2, markdown | no |
| node, npx, bun, deno | no |
| headless Chrome / Chromium | no |
| Safari, `cupsfilter` | yes |
| Python: plotly, pandas, xarray, yaml | yes |

So a server-rendered PDF is not producible here at all, and neither is a static raster of a
plotly chart. **HTML with inline SVG figures, printed to PDF from the browser, is not a
preference — it is the only route that works.** It is also the right one for a document:

- ~10 KB per figure against ~3.5 MB for the plotly runtime
- prints at printer resolution, not screen resolution
- survives JavaScript being disabled or stripped by a mail gateway
- diffable in git, which matters because a figure in a filing is an assertion

## Producing a PDF

```bash
open deliveries/<c>/<d>/report_compliance.html
# Safari: File ▸ Print ▸ PDF ▸ Save as PDF
```

The print stylesheet in `report_common.PAGE_CSS` sets a 16 mm page margin, forces the light
token set (a dark-surface figure prints as a black rectangle), repeats table headers across
pages, and marks figures, callouts, tables and rows `page-break-inside: avoid`. Navigation
carries `.no-print`.

**Pagination has never been verified here.** No browser can be driven on this machine. Any
session that generates a report and does not open it must report the report as *unreviewed*,
exactly as the layer workflow requires for maps.

## The figure vocabulary

`scripts/utils/report_figures.py`. Five figures, deliberately — a report needs four kinds of
chart, not forty.

| Function | Use for |
|---|---|
| `ranked_bar` | sites or hazards ordered by a 0–100 measure |
| `decade_line` | one line per forcing tier across decades |
| `band_stack` | how a population distributes across the five risk bands |
| `trend_strip` | signed slopes on a zero-centred diverging axis |
| `band_legend`, `empty_figure` | the band key; a figure that cannot be drawn |

Rules the module already enforces, which you should not undo:

- **No hex literal ever appears in an SVG.** Every colour is a CSS custom property emitted by
  `report_common.css_tokens()` from `viz_common`. One SVG serves light, dark and print.
- **Continuous risk is quantized to the five `RISK_BANDS` classes**, not drawn from the
  13-step ramp. A reader of a printed chart cannot hover to recover a value, so a shade they
  must map back to a number by eye is a shade that will be misread.
- **Tier series carry a distinct marker shape and dash unconditionally.** Blue and red sit
  1.12:1 apart in luminance, so there is effectively no lightness fallback if hue perception
  fails. These are the accessibility channel, not decoration.
- **Missing points break a line rather than interpolating.** An interpolated segment through
  a decade a layer does not publish is a fabricated number in a figure.
- **A figure that cannot be drawn says why.** An empty panel is indistinguishable from a
  rendering bug, and that ambiguity has already produced bug reports against the dashboard.
- **`trend_strip` hatches non-robust bars.** There is no p-value here; estimator disagreement
  is the robustness signal, and an unmarked bar asserts a trend the data does not support.
- **Pattern ids are content-derived** (`hatch-<8 hex>`). SVG ids are document-global even
  inside inline `<svg>`, so two identical ids on one page would make every `url(#…)` resolve
  to the first.

Every figure is paired with a data table by `report_common.figure_block()`. A figure
communicates shape; a filing needs the number.

## Why this is separate from `generate_maps.py` and the dashboard

Three renderers, three payload profiles, deliberately not merged:

| | payload | binding constraint |
|---|---|---|
| `generate_maps.py` | ~70k SVG markers per gridded panel | marker count |
| `generate_delivery_dashboard.py` | hundreds of points as embedded JSON | interactivity |
| `report_figures.py` | a few dozen values as inline SVG | print fidelity |

They share `scripts/utils/viz_common.py` — the validated palette, the scenario→tier colour
mapping, the symmetric-limit helper — and nothing else.

## Build stamps

Both reports carry an 8-hex content stamp in the header and support `--check-stamp`. A
cached page is indistinguishable from a fresh one by eye, and that ambiguity has produced two
phantom bug reports already. **When someone says a regenerated report looks wrong, check the
stamp before debugging the code.**

The stamp is matched on `class="stamp">build ([0-9a-f]{8})<`, not on a bare `build ` substring
— caveat prose legitimately contains the word "build", and a loose split once pulled a whole
sentence out as the stamp.

## Markdown

`report_common.markdown()` handles paragraphs, headings, bullet and numbered lists,
blockquotes, bold, italic, code and links. That is the entire grammar and it is intentional:
no markdown engine is installed, the narrative is written against this converter, and a
converter whose grammar fits on one screen cannot surprise a document filed with a regulator.

**Escaping happens before markup, never after.** A customer-supplied name containing markup
must render as text.
