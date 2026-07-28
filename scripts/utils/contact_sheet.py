"""Render a per-member contact sheet: one small global panel per ensemble member.

GUARDRAILS §11. Every distribution statistic used to value-check a member (min/max,
zero counts, unique values, units, calendar) is invariant under spatial
rearrangement, so none of them can see a spatial defect. This is the artifact that
can: a grid of thumbnails where a natively coarse member, a banded one, a
hemisphere flip or a broken land mask is obvious at a glance. A `~4°×5°` member once
passed the full tabular check twice plus 37 algebraic QA checks and shipped anyway,
because nobody had looked at an image of it.

Panels are **full 0.5° resolution, one pixel per grid cell, nearest-neighbour**.
Downsampling is deliberately avoided: it would smooth away the block structure the
sheet exists to reveal, which is the one thing it must not do.

Each panel is an embedded base64 PNG rather than a plotly trace. 17 members at
259,200 cells would be ~4.4M JSON numbers (tens of MB); as PNGs the whole sheet is
~1-2 MB, renders instantly, and stays pixel-exact. PNG encoding is pure stdlib
(zlib + struct), so this adds no dependency.
"""

import base64
import html
import struct
import zlib
from pathlib import Path

import numpy as np

# 32-point viridis, linearly interpolated to 256. Perceptually uniform, so a
# smooth field looks smooth and a blocky one looks blocky -- which is the point.
_VIRIDIS32 = [
    (68, 1, 84), (71, 13, 96), (72, 24, 106), (72, 35, 116), (71, 45, 123),
    (69, 55, 129), (66, 64, 134), (62, 73, 137), (59, 82, 139), (55, 91, 141),
    (51, 99, 141), (47, 107, 142), (44, 114, 142), (41, 122, 142), (38, 130, 142),
    (35, 137, 142), (33, 145, 140), (31, 152, 139), (31, 160, 136), (34, 167, 133),
    (40, 174, 128), (50, 182, 122), (63, 188, 115), (78, 195, 107), (94, 201, 98),
    (112, 207, 87), (132, 212, 75), (152, 216, 62), (173, 220, 48), (194, 223, 35),
    (216, 226, 25), (253, 231, 37),
]

# Diverging blue-white-red for signed fields (change/trend), white at zero.
_RDBU32 = [
    (5, 48, 97), (14, 66, 124), (24, 84, 148), (38, 101, 167), (57, 120, 182),
    (81, 139, 195), (108, 158, 207), (137, 177, 218), (166, 194, 226), (191, 209, 233),
    (211, 222, 239), (226, 233, 244), (238, 242, 248), (246, 248, 251), (251, 251, 252),
    (253, 253, 253), (253, 251, 249), (253, 245, 240), (252, 235, 226), (250, 222, 208),
    (246, 205, 187), (241, 185, 164), (233, 162, 140), (222, 138, 117), (208, 114, 96),
    (192, 90, 77), (174, 68, 61), (154, 48, 47), (132, 32, 37), (110, 20, 30),
    (86, 12, 24), (67, 0, 18),
]


def _lut(anchors):
    """Expand a 32-point colour list to a 256x3 uint8 lookup table."""
    a = np.asarray(anchors, dtype=float)
    xs = np.linspace(0, 255, len(a))
    out = np.empty((256, 3), dtype=np.uint8)
    for c in range(3):
        out[:, c] = np.interp(np.arange(256), xs, a[:, c]).round().astype(np.uint8)
    return out


_LUTS = {"viridis": _lut(_VIRIDIS32), "rdbu": _lut(_RDBU32)}


def _png(rgba: np.ndarray) -> bytes:
    """Encode an (h, w, 4) uint8 array as a PNG. Pure stdlib."""
    h, w = rgba.shape[:2]
    raw = b"".join(b"\x00" + rgba[y].tobytes() for y in range(h))  # filter byte 0

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 6))
            + chunk(b"IEND", b""))


def _panel_png(field: np.ndarray, vmin: float, vmax: float, cmap: str) -> str:
    """Colour-map a 2D field to a base64 PNG data URI; NaN renders transparent."""
    lut = _LUTS[cmap]
    finite = np.isfinite(field)
    span = (vmax - vmin) or 1.0
    idx = np.clip((field - vmin) / span, 0.0, 1.0)
    idx = np.where(finite, idx, 0.0)
    ints = (idx * 255).round().astype(np.uint8)
    rgba = np.zeros(field.shape + (4,), dtype=np.uint8)
    rgba[..., :3] = lut[ints]
    rgba[..., 3] = np.where(finite, 255, 0)          # ocean / off-mask = transparent
    return "data:image/png;base64," + base64.b64encode(_png(rgba)).decode("ascii")


def render_contact_sheet(panels, out_path, layer_id, decade, units="",
                         diverging=False, note="", columns=4):
    """Write a per-member contact-sheet HTML.

    Args:
        panels: ordered ``{member_label: 2D array (lat, lon)}``. Row order must
            match the source latitude vector; arrays are drawn as-is so a
            hemisphere flip stays visible rather than being silently corrected.
        out_path: destination ``.html`` path.
        layer_id: layer id, for the heading.
        decade: decade the panels represent, e.g. 2020.
        units: value units, shown in the shared-scale caption.
        diverging: use a zero-centred blue-white-red scale (for signed fields).
        note: optional extra caption line.
        columns: panels per row.

    Returns:
        The written ``Path``.
    """
    out_path = Path(out_path)
    if not panels:
        raise ValueError("no panels to render")

    pooled = np.concatenate([f[np.isfinite(f)].ravel() for f in panels.values()])
    if pooled.size == 0:
        raise ValueError("every panel is empty")

    # ONE shared colour scale across all panels. Per-panel autoscaling would hide
    # exactly what this sheet is for: a member whose magnitudes are far off the
    # others (the coarse member that shipped was also the most extreme-valued, and
    # per-panel scaling would have made it look unremarkable).
    if diverging:
        m = float(np.percentile(np.abs(pooled), 98)) or 1.0
        vmin, vmax, cmap = -m, m, "rdbu"
    else:
        vmin = float(np.percentile(pooled, 2))
        vmax = float(np.percentile(pooled, 98))
        cmap = "viridis"

    cards = []
    for label, field in panels.items():
        f = np.asarray(field, dtype=float)
        fin = f[np.isfinite(f)]
        stats = (f"median {np.median(fin):.3g} &middot; max {fin.max():.3g} "
                 f"&middot; {fin.size:,} cells" if fin.size else "no data")
        cards.append(
            f'<figure><img src="{_panel_png(f, vmin, vmax, cmap)}" '
            f'alt="{html.escape(label)}" loading="lazy">'
            f'<figcaption><b>{html.escape(label)}</b><br><span>{stats}</span>'
            f'</figcaption></figure>')

    ramp = "".join(
        f'<span style="background:rgb({r},{g},{b})"></span>'
        for r, g, b in _LUTS[cmap][::4])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(f"""<!doctype html>
<meta charset="utf-8"><title>Per-member contact sheet — {html.escape(layer_id)}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:1.5rem auto;
       max-width:1500px;padding:0 1rem;color:#222;background:#fff}}
 h1{{font-size:1.25rem;margin:0 0 .25rem}}
 .sub{{color:#555;font-size:.9rem;margin:0 0 1rem}}
 .grid{{display:grid;grid-template-columns:repeat({columns},minmax(0,1fr));gap:1rem}}
 figure{{margin:0}}
 img{{width:100%;height:auto;display:block;image-rendering:pixelated;
      background:#eaf2fb;border:1px solid #ddd}}
 figcaption{{font-size:.8rem;line-height:1.3;padding-top:.3rem}}
 figcaption span{{color:#666}}
 .scale{{display:flex;align-items:center;gap:.5rem;margin:.75rem 0 1.25rem;
         font-size:.85rem;color:#444}}
 .ramp{{display:flex;height:12px;width:260px;border:1px solid #ccc}}
 .ramp span{{flex:1}}
 .why{{background:#f7f7f9;border-left:3px solid #999;padding:.6rem .8rem;
       font-size:.85rem;color:#444;margin:0 0 1.25rem}}
 a{{color:#1a5fb4}}
</style>
<h1>Per-member contact sheet — {html.escape(layer_id)}</h1>
<p class="sub">{decade}s decadal mean, one panel per ensemble member.
Full 0.5&deg; resolution, one pixel per grid cell, nearest-neighbour &mdash; not
downsampled. Transparent = outside that member's land mask.
{html.escape(note)}</p>
<div class="scale"><span>shared scale</span>
  <span>{vmin:.3g}</span><div class="ramp">{ramp}</div><span>{vmax:.3g}</span>
  <span>{html.escape(units)}</span></div>
<p class="why"><b>What to look for:</b> block structure or seams (a member running
natively coarser than the grid it declares), banding, hard unrealistic edges, a
land mask that disagrees with its neighbours, hemisphere flips, and hot or cold
patches unrelated to geography. The scale is shared across panels, so a member
whose magnitudes sit far from the rest stands out too. Distribution statistics
cannot show any of this &mdash; see GUARDRAILS &sect;11.</p>
<div class="grid">{''.join(cards)}</div>
""")
    return out_path
