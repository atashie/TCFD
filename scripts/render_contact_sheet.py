"""Render a per-member contact sheet to PNG so it can actually be LOOKED AT.

The `/isimip-process-visualize` skill requires viewing one small global panel per member
before choosing statistics, because every statistic in a value-check table is invariant
under spatial rearrangement and therefore cannot see a spatial defect: a ~4x5 degree member
once passed the full table twice and 37 algebraic QA checks, and a human caught it by
looking. The plotly dashboard's Members tab covers this in a browser, but no browser can be
driven in this environment, so nothing there is ever seen.

This writes a PNG instead, which CAN be opened and inspected here.

There is no matplotlib and no Pillow in this venv, so the PNG is encoded directly from
numpy with zlib -- a plain 8-bit RGB, non-interlaced, one IDAT. That is the whole reason
this file exists rather than three lines of `plt.imsave`.

Look for: block structure (a member on a coarser native grid), seams and banding, hard
unrealistic edges, land-mask errors, hemisphere flips, and patches unrelated to geography.
Check PER MEMBER -- the pooled ensemble dilutes one bad member.

Usage:
    python scripts/render_contact_sheet.py {members.nc} {out.png} [--cols 5] [--clip-pct 99]
    python scripts/render_contact_sheet.py {processed.nc} {out.png} --var percentile
"""

import argparse
import struct
import sys
import zlib
from pathlib import Path

import numpy as np
import xarray as xr

#: Perceptually-ordered dark->bright ramp (viridis control points, linearly interpolated).
#: Bright = high, so a hot member stands out against its neighbours in the grid.
_VIRIDIS = np.array([
    (68, 1, 84), (72, 40, 120), (62, 74, 137), (49, 104, 142), (38, 130, 142),
    (31, 158, 137), (53, 183, 121), (109, 205, 89), (180, 222, 44), (253, 231, 37),
], dtype=np.float64)

NAN_RGB = (28, 28, 32)      # near-black so masked area reads as absent, not as low
LABEL_RGB = (255, 255, 255)

#: 5x7 bitmap font, enough to label members. Anything unmapped renders blank.
_FONT = {
    "A": ["01110", "10001", "10001", "11111", "10001", "10001", "10001"],
    "B": ["11110", "10001", "10001", "11110", "10001", "10001", "11110"],
    "C": ["01110", "10001", "10000", "10000", "10000", "10001", "01110"],
    "D": ["11110", "10001", "10001", "10001", "10001", "10001", "11110"],
    "E": ["11111", "10000", "10000", "11110", "10000", "10000", "11111"],
    "F": ["11111", "10000", "10000", "11110", "10000", "10000", "10000"],
    "G": ["01110", "10001", "10000", "10111", "10001", "10001", "01111"],
    "H": ["10001", "10001", "10001", "11111", "10001", "10001", "10001"],
    "I": ["11111", "00100", "00100", "00100", "00100", "00100", "11111"],
    "J": ["00111", "00010", "00010", "00010", "00010", "10010", "01100"],
    "K": ["10001", "10010", "10100", "11000", "10100", "10010", "10001"],
    "L": ["10000", "10000", "10000", "10000", "10000", "10000", "11111"],
    "M": ["10001", "11011", "10101", "10101", "10001", "10001", "10001"],
    "N": ["10001", "11001", "10101", "10011", "10001", "10001", "10001"],
    "O": ["01110", "10001", "10001", "10001", "10001", "10001", "01110"],
    "P": ["11110", "10001", "10001", "11110", "10000", "10000", "10000"],
    "Q": ["01110", "10001", "10001", "10001", "10101", "10010", "01101"],
    "R": ["11110", "10001", "10001", "11110", "10100", "10010", "10001"],
    "S": ["01111", "10000", "10000", "01110", "00001", "00001", "11110"],
    "T": ["11111", "00100", "00100", "00100", "00100", "00100", "00100"],
    "U": ["10001", "10001", "10001", "10001", "10001", "10001", "01110"],
    "V": ["10001", "10001", "10001", "10001", "10001", "01010", "00100"],
    "W": ["10001", "10001", "10001", "10101", "10101", "11011", "10001"],
    "X": ["10001", "01010", "00100", "00100", "00100", "01010", "10001"],
    "Y": ["10001", "01010", "00100", "00100", "00100", "00100", "00100"],
    "Z": ["11111", "00001", "00010", "00100", "01000", "10000", "11111"],
    "0": ["01110", "10001", "10011", "10101", "11001", "10001", "01110"],
    "1": ["00100", "01100", "00100", "00100", "00100", "00100", "01110"],
    "2": ["01110", "10001", "00001", "00110", "01000", "10000", "11111"],
    "3": ["11111", "00010", "00100", "00010", "00001", "10001", "01110"],
    "4": ["00010", "00110", "01010", "10010", "11111", "00010", "00010"],
    "5": ["11111", "10000", "11110", "00001", "00001", "10001", "01110"],
    "6": ["00110", "01000", "10000", "11110", "10001", "10001", "01110"],
    "7": ["11111", "00001", "00010", "00100", "01000", "01000", "01000"],
    "8": ["01110", "10001", "10001", "01110", "10001", "10001", "01110"],
    "9": ["01110", "10001", "10001", "01111", "00001", "00010", "01100"],
    "-": ["00000", "00000", "00000", "11111", "00000", "00000", "00000"],
    ".": ["00000", "00000", "00000", "00000", "00000", "01100", "01100"],
    " ": ["00000", "00000", "00000", "00000", "00000", "00000", "00000"],
}


def colormap(norm):
    """(H, W) in [0,1] (NaN allowed) -> (H, W, 3) uint8."""
    out = np.empty(norm.shape + (3,), np.uint8)
    bad = ~np.isfinite(norm)
    x = np.clip(np.nan_to_num(norm), 0.0, 1.0) * (len(_VIRIDIS) - 1)
    lo = np.floor(x).astype(int)
    hi = np.minimum(lo + 1, len(_VIRIDIS) - 1)
    t = (x - lo)[..., None]
    out[:] = (_VIRIDIS[lo] * (1 - t) + _VIRIDIS[hi] * t).astype(np.uint8)
    out[bad] = NAN_RGB
    return out


def draw_text(canvas, text, y, x, rgb=LABEL_RGB, scale=1):
    """Blit `text` into an (H, W, 3) uint8 canvas at row y, col x."""
    for ch in text.upper():
        glyph = _FONT.get(ch)
        if glyph is None:
            x += 6 * scale
            continue
        for r, row in enumerate(glyph):
            for c, bit in enumerate(row):
                if bit == "1":
                    y0, x0 = y + r * scale, x + c * scale
                    if 0 <= y0 < canvas.shape[0] - scale and 0 <= x0 < canvas.shape[1] - scale:
                        canvas[y0:y0 + scale, x0:x0 + scale] = rgb
        x += 6 * scale
    return x


def write_png(path, rgb):
    """Minimal 8-bit RGB PNG. No matplotlib, no Pillow in this venv."""
    h, w, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[i].tobytes() for i in range(h))

    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    png = (b"\x89PNG\r\n\x1a\n"
           + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 6))
           + chunk(b"IEND", b""))
    Path(path).write_bytes(png)
    return len(png)


def block_mean(a, k):
    """Downsample by k with a block MEAN. Never slice -- `a[::k, ::k]` samples every kth
    cell and silently deletes exposed cells on a sparse hazard."""
    h, w = a.shape
    a = a[:h - h % k, :w - w % k]
    return np.nanmean(a.reshape(a.shape[0] // k, k, a.shape[1] // k, k), axis=(1, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--var", default=None,
                    help="variable to render; default 'value' for a members file")
    ap.add_argument("--cols", type=int, default=5)
    ap.add_argument("--stride", type=int, default=2, help="block-mean downsample factor")
    ap.add_argument("--clip-pct", type=float, default=99.0,
                    help="shared colour limit = this percentile of finite values across "
                         "ALL panels; a shared scale is the point -- a per-panel scale "
                         "would hide exactly the hot member this is meant to expose")
    ap.add_argument("--decade", type=int, default=None,
                    help="for a processed file with a `decade` dim, which decade to render")
    args = ap.parse_args()

    ds = xr.open_dataset(args.src)
    var = args.var or ("value" if "value" in ds.data_vars else list(ds.data_vars)[0])
    da = ds[var]

    if "member" in da.dims:
        panels = [(str(m), da.sel(member=m).values) for m in ds["member"].values]
        title = f"{Path(args.src).name}  var={var}  ({len(panels)} members)"
    elif "decade" in da.dims:
        decs = ds["decade"].values
        if args.decade is not None:
            decs = [d for d in decs if int(d) == args.decade]
        panels = [(f"{int(d)}S", da.sel(decade=d).values) for d in decs]
        title = f"{Path(args.src).name}  var={var}"
    else:
        panels = [(var, da.values)]
        title = Path(args.src).name

    # ORIENT NORTH-UP BEFORE RENDERING. Array row 0 is whatever the file happens to store
    # first, and ISIMIP files disagree: the sea-level and GEBCO grids ascend (row 0 = -90),
    # the ISIMIP2b land mask descends. Drawing rows top-down without checking renders an
    # ascending file upside down -- which is worse here than anywhere else, because this
    # sheet exists to catch hemisphere flips and would instead manufacture one. Caught
    # 2026-08-14 on `sealevel-2b`: its dense Arctic coastline (7,141 cells at 60-90N,
    # against 0 in Antarctica) drew as a bright band along the BOTTOM edge and read as
    # Antarctic contamination.
    lat_name = next((c for c in ("lat", "latitude", "y") if c in ds.coords), None)
    flip = False
    if lat_name is not None:
        lv = np.asarray(ds[lat_name].values, float)
        flip = lv.size > 1 and lv[1] > lv[0]
    if flip:
        panels = [(n, np.asarray(v)[::-1, :]) for n, v in panels]
        print("  latitude ascends in this file -> flipped to render north-up")

    grids = [block_mean(np.asarray(v, float), args.stride) for _, v in panels]
    allv = np.concatenate([g[np.isfinite(g)] for g in grids])
    if allv.size == 0:
        print("ERROR: no finite values to render")
        return 2
    vmax = float(np.percentile(allv, args.clip_pct))
    vmin = float(allv.min())
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = float(allv.max()) or 1.0

    ph, pw = grids[0].shape
    cols = min(args.cols, len(grids))
    rows = (len(grids) + cols - 1) // cols
    pad, top = 4, 22
    lab = 10
    H = top + rows * (ph + lab + pad) + pad
    W = pad + cols * (pw + pad)
    canvas = np.zeros((H, W, 3), np.uint8)
    canvas[:] = (18, 18, 22)

    draw_text(canvas, title[:110], 6, pad)
    for i, ((name, _), g) in enumerate(zip(panels, grids)):
        r, c = divmod(i, cols)
        y = top + r * (ph + lab + pad) + lab
        x = pad + c * (pw + pad)
        norm = (g - vmin) / (vmax - vmin)
        canvas[y:y + ph, x:x + pw] = colormap(norm)
        draw_text(canvas, name[:pw // 6], y - lab + 1, x)

    n = write_png(args.out, canvas)
    print(f"wrote {args.out}  {W}x{H}px  {n/1024:.0f} KB")
    print(f"  {len(panels)} panels, shared scale [{vmin:.5g}, {vmax:.5g}] "
          f"(p{args.clip_pct} of all finite values), stride {args.stride}")
    print(f"  dark slate = NaN/masked, dark purple = low, yellow = high")
    return 0


if __name__ == "__main__":
    sys.exit(main())
