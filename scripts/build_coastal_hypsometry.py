"""Stage 2 of the coastal layer: per-0.5-degree coastal-land elevation histograms.

Reduces the GEBCO_2026 15-arcsec global terrain grid to one elevation HISTOGRAM per 0.5
degree cell, over land that is coastal and hydrologically connected to the ocean. Stage 3
then slides each histogram against the sea-level delta to get exposure and percentile,
which is why a histogram is stored rather than a summary: the whole hypsometric curve is
needed, and it is needed at every decade, so reducing it now would mean re-reading 7 GB
per decade.

MEMORY IS THE BINDING CONSTRAINT -- 16 GB, AND THE GRID IS 7.5 GB
-----------------------------------------------------------------
GEBCO_2026 is 43200 x 86400 int16 = 7.46 GB. Loading it whole leaves nothing for the
label array, so every pass here streams in latitude bands and nothing global is ever held
at 15 arcsec. The one global array that IS held is the 1-arcmin connectivity grid
(10800 x 21600), which is 233 MB as bool and 933 MB as int32 labels.

WHY CONNECTIVITY IS COMPUTED AT ALL
-----------------------------------
"Below sea level" is not "floods". The Qattara Depression sits at -133 m about 60 km from
the Mediterranean and the Dead Sea shore at -430 m about 70 km from it -- both inside any
100 km coastal buffer, and both would read as catastrophically inundated under a plain
elevation threshold. They are not connected to the sea, and neither is the Caspian, the
Turfan Basin or Death Valley. So the ocean is defined as the connected component of water,
not as "elevation <= 0", and a cell counts as exposed only if it is in that component.

The connectivity grid is min-reduced from 15 arcsec to 1 arcmin. MIN, not mean: a channel
narrower than the block still reads as water, so straits stay open. That array is used
ONLY to decide what is connected -- never for an elevation value.

Connectivity is evaluated once at the highest water level considered rather than per
decade. Sea-level change here is sub-metre, so the connected set barely moves, and
evaluating at the top of the range is the permissive end -- any cell that could connect,
does.

WHAT THIS CANNOT SEE, AND IT MATTERS
------------------------------------
At 15 arcsec a grid cell is ~460 m at the equator. Dikes, levees, seawalls and barriers
are narrower than that and are absent from the DEM, so a defended coast reads exactly like
an undefended one. The Netherlands is the standing example: its polders are below sea
level and its dikes are sub-pixel, so this method will show them connected and exposed.
That is a property of every bathtub-style elevation model and it belongs in the layer's
caveats, not in a footnote.

GEBCO is a DIGITAL SURFACE model -- it carries vegetation canopy and buildings. In coastal
zones SRTM-class surface models run +2.49 m (Australia) to +3.67 m (US) high, against a
sea-level signal here of +0.21 to +0.34 m. The bias is an order of magnitude larger than
the signal, which is why Stage 3 reports a CHANGE in exposed area computed from the same
DEM at both ends rather than an absolute inundation count.

Run:  .venv/bin/python3 scripts/build_coastal_hypsometry.py --inspect   # header only
      .venv/bin/python3 scripts/build_coastal_hypsometry.py
Out:  data/interim/coastal/ocean_connected_1arcmin.npz
      data/interim/coastal/hypsometry_halfdeg.nc
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import ndimage

GEBCO = Path("data/raw/gebco")
OUT = Path("data/interim/coastal")

FINE_ROWS, FINE_COLS = 43200, 86400      # 15 arcsec
COARSE_FACTOR = 4                        # -> 1 arcmin for connectivity
BAND_ROWS = 1440                         # 1440 x 86400 int16 = 249 MB per read

# Water level at which connectivity is evaluated: the permissive end of the range.
CONNECT_LEVEL_M = 1.0

# Histogram over elevation-above-present-sea-level. 0.05 m bins across [-30, +30] resolve
# a 0.2-0.8 m signal into 4-16 bins; everything outside is kept in the two tail counters
# so no land is silently dropped from the denominator.
BIN_LO, BIN_HI, BIN_W = -30.0, 30.0, 0.05
N_BINS = int(round((BIN_HI - BIN_LO) / BIN_W))

HALF_DEG = 0.5
NLAT_OUT, NLON_OUT = 360, 720
SUB = FINE_ROWS // NLAT_OUT              # 120 fine rows per 0.5 deg row
assert SUB * NLAT_OUT == FINE_ROWS and (FINE_COLS // NLON_OUT) * NLON_OUT == FINE_COLS

# How far inland to bother. The user's constraint is 100 km; applied at 0.5 deg
# granularity, which is ~55 km per cell, so two rings is ~110 km. Land inside a coastal
# 0.5 deg cell is within ~80 km of that cell's own coast anyway, so this is a compute
# filter rather than a scientific threshold, and it is deliberately the generous side.
INLAND_RINGS = 2


def find_grid() -> Path:
    cands = sorted(GEBCO.glob("*.nc"))
    if not cands:
        raise SystemExit(f"no GEBCO .nc found in {GEBCO} -- is the download finished?")
    return max(cands, key=lambda p: p.stat().st_size)


def inspect(path: Path):
    ds = xr.open_dataset(path, decode_times=False)
    print(f"file    : {path}  ({path.stat().st_size / 1e9:.2f} GB)")
    print(f"dims    : {dict(ds.sizes)}")
    print(f"vars    : {list(ds.data_vars)}")
    for v in ds.data_vars:
        print(f"  {v}: dtype={ds[v].dtype} attrs={ds[v].attrs}")
    for c in ds.coords:
        vals = ds[c].values
        print(f"coord {c}: {vals[:2]} ... {vals[-2:]}  (n={vals.size}, "
              f"{'ascending' if vals[1] > vals[0] else 'DESCENDING'})")
    print(f"global  : { {k: ds.attrs[k] for k in list(ds.attrs)[:6]} }")
    ds.close()


def elevation_var(ds) -> str:
    for name in ("elevation", "z", "Band1"):
        if name in ds.data_vars:
            return name
    raise SystemExit(f"cannot find an elevation variable among {list(ds.data_vars)}")


def pass_a_connectivity(path: Path) -> np.ndarray:
    """Ocean-connected water at 1 arcmin, as a bool array (n/4 x m/4)."""
    cache = OUT / "ocean_connected_1arcmin.npz"
    if cache.exists():
        print(f"pass A: reusing {cache}")
        return np.load(cache)["ocean"]

    ds = xr.open_dataset(path, decode_times=False, chunks=None)
    var = elevation_var(ds)
    nr, nc = ds[var].shape
    cr, cc = nr // COARSE_FACTOR, nc // COARSE_FACTOR
    coarse = np.empty((cr, cc), dtype="int16")

    print(f"pass A: min-reducing {nr}x{nc} -> {cr}x{cc} in bands of {BAND_ROWS}")
    for r0 in range(0, nr, BAND_ROWS):
        r1 = min(r0 + BAND_ROWS, nr)
        band = ds[var][r0:r1].values.astype("int16")
        b = band.reshape(band.shape[0] // COARSE_FACTOR, COARSE_FACTOR,
                         cc, COARSE_FACTOR)
        coarse[r0 // COARSE_FACTOR:r1 // COARSE_FACTOR] = b.min(axis=(1, 3))
        del band, b
        print(f"  rows {r0:6d}-{r1:6d}", end="\r")
    ds.close()
    print()

    water = coarse <= CONNECT_LEVEL_M
    del coarse
    lab, n = ndimage.label(water, structure=np.ones((3, 3), np.uint8))
    print(f"pass A: {n:,} water components")

    # Join across the antimeridian, then take the largest -- the world ocean.
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for r in range(lab.shape[0]):
        for dr in (-1, 0, 1):
            rr = r + dr
            if 0 <= rr < lab.shape[0] and lab[r, 0] and lab[rr, -1]:
                a, b = find(lab[r, 0]), find(lab[rr, -1])
                if a != b:
                    parent[b] = a
    roots = np.array([find(i) for i in range(n + 1)])
    roots[0] = 0
    sizes = np.bincount(roots[lab].ravel())
    sizes[0] = 0
    ocean_root = int(sizes.argmax())
    ocean = roots[lab] == ocean_root
    print(f"pass A: ocean component = {ocean.sum():,} of {water.sum():,} water cells "
          f"({100 * ocean.sum() / water.sum():.1f}%)")

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, ocean=ocean)
    return ocean


def pass_b_hypsometry(path: Path, ocean_c: np.ndarray):
    """Per-0.5-degree histogram of coastal-land elevation."""
    # Which 0.5 deg cells to process: those containing land, within INLAND_RINGS of a cell
    # containing ocean.
    f = ocean_c.shape[0] // NLAT_OUT                      # 1-arcmin rows per 0.5 deg
    g = ocean_c.shape[1] // NLON_OUT
    has_ocean = ocean_c.reshape(NLAT_OUT, f, NLON_OUT, g).any(axis=(1, 3))
    near = has_ocean.copy()
    for _ in range(INLAND_RINGS):
        near = (near
                | np.roll(near, 1, 0) | np.roll(near, -1, 0)
                | np.roll(near, 1, 1) | np.roll(near, -1, 1))
    has_land = ~ocean_c.reshape(NLAT_OUT, f, NLON_OUT, g).all(axis=(1, 3))
    todo = near & has_land
    print(f"pass B: {todo.sum():,} of {NLAT_OUT * NLON_OUT:,} half-degree cells to process")

    ds = xr.open_dataset(path, decode_times=False)
    var = elevation_var(ds)
    hist = np.zeros((NLAT_OUT, NLON_OUT, N_BINS), dtype="int32")
    below = np.zeros((NLAT_OUT, NLON_OUT), dtype="int32")
    above = np.zeros((NLAT_OUT, NLON_OUT), dtype="int32")
    n_land = np.zeros((NLAT_OUT, NLON_OUT), dtype="int32")

    subc = FINE_COLS // NLON_OUT
    for jr in range(NLAT_OUT):
        if not todo[jr].any():
            continue
        band = ds[var][jr * SUB:(jr + 1) * SUB].values.astype("float32")
        # Upsample the 1-arcmin ocean flag to 15 arcsec for this band.
        oc = np.repeat(np.repeat(ocean_c[jr * f:(jr + 1) * f], COARSE_FACTOR, 0),
                       COARSE_FACTOR, 1)
        land = ~oc
        for ic in np.flatnonzero(todo[jr]):
            blk = band[:, ic * subc:(ic + 1) * subc]
            lm = land[:, ic * subc:(ic + 1) * subc]
            vals = blk[lm]
            if vals.size == 0:
                continue
            n_land[jr, ic] = vals.size
            below[jr, ic] = int((vals < BIN_LO).sum())
            above[jr, ic] = int((vals >= BIN_HI).sum())
            inb = vals[(vals >= BIN_LO) & (vals < BIN_HI)]
            if inb.size:
                idx = ((inb - BIN_LO) / BIN_W).astype("int32")
                hist[jr, ic] = np.bincount(idx, minlength=N_BINS)
        del band, oc, land
        print(f"  row {jr + 1}/{NLAT_OUT}", end="\r")
    ds.close()
    print()

    edges = BIN_LO + BIN_W * np.arange(N_BINS + 1)
    out = xr.Dataset(
        {
            "hist": (("lat", "lon", "bin"), hist,
                     {"long_name": "count of land sub-cells per elevation bin"}),
            "n_below": (("lat", "lon"), below,
                        {"long_name": f"land sub-cells below {BIN_LO} m"}),
            "n_above": (("lat", "lon"), above,
                        {"long_name": f"land sub-cells at or above {BIN_HI} m"}),
            "n_land": (("lat", "lon"), n_land,
                       {"long_name": "total ocean-disconnected (land) sub-cells"}),
            "processed": (("lat", "lon"), todo.astype("int8"),
                          {"long_name": "cell was within the coastal buffer and had land"}),
        },
        coords={
            "lat": np.arange(-90 + HALF_DEG / 2, 90, HALF_DEG),
            "lon": np.arange(-180 + HALF_DEG / 2, 180, HALF_DEG),
            "bin_left": ("bin", edges[:-1]),
        },
        attrs={
            "dem": "GEBCO_2026 ice-surface elevation, 15 arcsec",
            "dem_datum": "GEBCO states its grids assume all source data referred to mean "
                         "sea level",
            "dem_type": "DIGITAL SURFACE MODEL -- includes vegetation canopy and buildings; "
                        "SRTM-class coastal bias is +2.49 to +3.67 m, an order of magnitude "
                        "larger than the sea-level signal",
            "connectivity": f"ocean = largest connected water component at 1 arcmin, "
                            f"evaluated at {CONNECT_LEVEL_M} m; min-reduced from 15 arcsec "
                            f"so narrow straits stay open",
            "defences_represented": "NO -- dikes, levees and seawalls are narrower than a "
                                    "grid cell and absent from the DEM",
            "inland_buffer": f"{INLAND_RINGS} half-degree rings from a cell containing ocean",
            "bin_low_m": BIN_LO, "bin_high_m": BIN_HI, "bin_width_m": BIN_W,
        },
    )
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "hypsometry_halfdeg.nc"
    out.to_netcdf(p, encoding={"hist": {"zlib": True, "complevel": 4}})
    print(f"pass B: wrote {p} ({p.stat().st_size / 1e6:.0f} MB)")
    print(f"  cells with land: {int((n_land > 0).sum()):,}; "
          f"total land sub-cells: {n_land.sum():,}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true", help="print the GEBCO header and stop")
    a = ap.parse_args()
    path = find_grid()
    if a.inspect:
        inspect(path)
        return
    ocean_c = pass_a_connectivity(path)
    pass_b_hypsometry(path, ocean_c)


if __name__ == "__main__":
    main()
