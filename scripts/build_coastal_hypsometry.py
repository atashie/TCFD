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

# Water level at which connectivity is evaluated. This is the LECZ bound itself, and not
# by coincidence: McGranahan, Balk & Anderson (2007) define the Low Elevation Coastal Zone
# as "the contiguous area along the coast that is less than 10 metres above sea level",
# hydrologically connected to the sea -- which is a connectivity criterion evaluated at
# 10 m. Setting it lower truncates the hypsometry instead of just excluding closed basins:
# at 1 m, land sitting at 3 m fell outside the connectable set and was counted as
# "unconnected" when it is merely higher, which emptied 95% of coastal land out of the
# histogram.
CONNECT_LEVEL_M = 10.0

# ONE CONNECTIVITY LEVEL CANNOT DO BOTH JOBS. Evaluated at 10 m, the Salton Sink connects
# to the Gulf of California across the Colorado delta's alluvial fan and the Danakil
# Depression connects to the Red Sea -- so the Imperial Valley, 100 km inland at -72 m,
# would read as maximum coastal risk. Evaluated at 2 m instead, the percentile grading
# collapses: coastal plain at 5 m falls outside the connected set and scores 1 as though it
# were a mountain. So the two questions are asked separately. EXPO_LEVEL_M gates what can
# actually be inundated and only has to clear the largest sea-level delta in the ensemble
# (measured: +1.328 m at rcp60); CONNECT_LEVEL_M gates LECZ pool membership and the graded
# percentile.
EXPO_LEVEL_M = 2.0

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


def find_grid() -> tuple[Path, Path]:
    """(elevation grid, type-identifier grid). Named explicitly, not picked by size."""
    elev, tid = GEBCO / "GEBCO_2026.nc", GEBCO / "GEBCO_2026_TID.nc"
    for f in (elev, tid):
        if not f.exists():
            raise SystemExit(f"missing {f} -- is the download finished?")
    return elev, tid


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


def pass_a_connectivity(elev_path: Path, tid_path: Path) -> np.ndarray:
    """Ocean-connected water at 1 arcmin, as a bool array (n/4 x m/4).

    "Connectable" is marine water OR land low enough to be reachable by the sea at
    CONNECT_LEVEL_M. Marine water comes from the TID grid rather than from an elevation
    threshold, which is the correction that matters: taking `elevation <= 1 m` as sea
    swallowed the Dutch polders -- below sea level, dikes narrower than a 460 m cell -- so
    the most exposed land on the planet was reclassified as ocean and left the denominator
    entirely. TID calls those cells Land (verified: 52.4N 4.9E -> TID 0) and the North Sea
    beside them 17.
    """
    cache = OUT / "ocean_connected_1arcmin.npz"
    if cache.exists():
        print(f"pass A: reusing {cache}")
        z = np.load(cache)
        return {"expo": z["expo"], "pool": z["pool"]}

    de = xr.open_dataset(elev_path, decode_times=False)
    dt = xr.open_dataset(tid_path, decode_times=False)
    var = elevation_var(de)
    nr, nc = de[var].shape
    cr, cc = nr // COARSE_FACTOR, nc // COARSE_FACTOR
    levels = {"expo": EXPO_LEVEL_M, "pool": CONNECT_LEVEL_M}
    water = {k: np.empty((cr, cc), dtype=bool) for k in levels}

    print(f"pass A: reducing {nr}x{nc} -> {cr}x{cc} in bands of {BAND_ROWS}")
    for r0 in range(0, nr, BAND_ROWS):
        r1 = min(r0 + BAND_ROWS, nr)
        e = de[var][r0:r1].values.astype("int16")
        t = dt["tid"][r0:r1].values.astype("int8")
        marine = t != 0
        del t
        for k, lv in levels.items():
            conn = marine | (e <= lv)
            # ANY connectable sub-cell makes the 1-arcmin cell connectable, so channels
            # narrower than the coarse cell stay open.
            water[k][r0 // COARSE_FACTOR:r1 // COARSE_FACTOR] = conn.reshape(
                conn.shape[0] // COARSE_FACTOR, COARSE_FACTOR, cc,
                COARSE_FACTOR).any(axis=(1, 3))
            del conn
        del e, marine
        print(f"  rows {r0:6d}-{r1:6d}", end="\r")
    de.close(); dt.close()
    print()
    out = {}
    for key in ("expo", "pool"):
        out[key] = _largest_component(water.pop(key), key)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, expo=out["expo"], pool=out["pool"])
    return out


def _largest_component(water: np.ndarray, label: str) -> np.ndarray:
    lab, n = ndimage.label(water, structure=np.ones((3, 3), np.uint8))
    print(f"pass A [{label}]: {n:,} water components")

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
    roots = np.array([find(i) for i in range(n + 1)], dtype=np.int32)
    roots[0] = 0
    lab = roots[lab]            # relabel in place-ish; drops the pre-union labels
    del parent
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    ocean_root = int(sizes.argmax())
    ocean = lab == ocean_root
    del lab
    print(f"pass A [{label}]: ocean component = {ocean.sum():,} of {water.sum():,} "
          f"water cells ({100 * ocean.sum() / water.sum():.1f}%)")
    return ocean


def pass_b_hypsometry(elev_path: Path, tid_path: Path, ocean_c: np.ndarray):
    """Per-0.5-degree histogram of coastal-land elevation.

    LAND COMES FROM TID, EXPOSURE FROM ELEVATION. Land is TID == 0, the publisher's own
    verdict, so below-sea-level land behind sub-pixel defences stays in the denominator
    instead of being flood-filled into the ocean.

    Land is split into FLOODABLE (hydrologically connected to the sea) and UNCONNECTED.
    Only floodable land enters the histogram, because land the sea cannot reach cannot be
    inundated by it however low it sits -- the Dead Sea shore at -430 m is not at coastal
    risk. Unconnected land is still counted, so it dilutes the exposure fraction and takes
    the safe end of the percentile rather than vanishing.
    """
    pool_c, expo_c = ocean_c["pool"], ocean_c["expo"]
    f = pool_c.shape[0] // NLAT_OUT
    g = pool_c.shape[1] // NLON_OUT
    has_ocean = expo_c.reshape(NLAT_OUT, f, NLON_OUT, g).any(axis=(1, 3))
    near = has_ocean.copy()
    for _ in range(INLAND_RINGS):
        near = (near
                | np.roll(near, 1, 0) | np.roll(near, -1, 0)
                | np.roll(near, 1, 1) | np.roll(near, -1, 1))
    todo = near
    print(f"pass B: {todo.sum():,} of {NLAT_OUT * NLON_OUT:,} half-degree cells in the "
          f"coastal buffer")

    de = xr.open_dataset(elev_path, decode_times=False)
    dt = xr.open_dataset(tid_path, decode_times=False)
    var = elevation_var(de)
    cell_ix, rows_hist, rows_below, rows_above, rows_nland, rows_unconn = [], [], [], [], [], []
    rows_ehist, rows_ebelow = [], []
    subc = FINE_COLS // NLON_OUT

    for jr in range(NLAT_OUT):
        if not todo[jr].any():
            continue
        band = de[var][jr * SUB:(jr + 1) * SUB].values.astype("float32")
        tband = dt["tid"][jr * SUB:(jr + 1) * SUB].values.astype("int8")
        conn = np.repeat(np.repeat(pool_c[jr * f:(jr + 1) * f], COARSE_FACTOR, 0),
                         COARSE_FACTOR, 1)
        econn = np.repeat(np.repeat(expo_c[jr * f:(jr + 1) * f], COARSE_FACTOR, 0),
                          COARSE_FACTOR, 1)
        is_land = tband == 0
        for ic in np.flatnonzero(todo[jr]):
            sl = slice(ic * subc, (ic + 1) * subc)
            lm = is_land[:, sl]
            if not lm.any():
                continue
            # Reachable = land inside the ocean-connected LECZ component. Everything else
            # -- land above the LECZ bound, and land below it sitting in a closed basin
            # like the Dead Sea or Qattara -- is SAFE from the sea and is counted but never
            # exposed. Both go in the denominator so the exposure fraction is a share of
            # the cell's real coastal land, not of a subset chosen by the method.
            def _hist(mask):
                v = band[:, sl][mask]
                hh = np.zeros(N_BINS, dtype="int32")
                nbel = 0
                if v.size:
                    nbel = int((v < BIN_LO).sum())
                    ib = v[(v >= BIN_LO) & (v < BIN_HI)]
                    if ib.size:
                        hh = np.bincount(((ib - BIN_LO) / BIN_W).astype("int32"),
                                         minlength=N_BINS).astype("int32")
                return hh, nbel

            h, n_bel = _hist(lm & conn[:, sl])            # LECZ pool -> graded percentile
            eh, e_bel = _hist(lm & econn[:, sl])          # floodable -> exposure only
            n_land_c = int(lm.sum())
            cell_ix.append(jr * NLON_OUT + ic)
            rows_hist.append(h)
            rows_below.append(n_bel)
            rows_above.append(n_land_c - int(h.sum()) - n_bel)   # safe land
            rows_nland.append(n_land_c)
            rows_unconn.append(int((lm & ~conn[:, sl]).sum()))
            rows_ehist.append(eh)
            rows_ebelow.append(e_bel)
        del band, tband, conn, econn, is_land
        print(f"  row {jr + 1}/{NLAT_OUT}", end="\r")
    de.close(); dt.close()
    print()

    edges = BIN_LO + BIN_W * np.arange(N_BINS + 1)
    hist = np.asarray(rows_hist, dtype="int32") if rows_hist else np.zeros((0, N_BINS), "int32")
    ehist = np.asarray(rows_ehist, dtype="int32") if rows_ehist else np.zeros((0, N_BINS), "int32")
    n_land = np.asarray(rows_nland, dtype="int32")
    unconn = np.asarray(rows_unconn, dtype="int32")
    out = xr.Dataset(
        {
            "hist": (("cell", "bin"), hist,
                     {"long_name": "count of FLOODABLE land sub-cells per elevation bin"}),
            "n_below": (("cell",), np.asarray(rows_below, dtype="int32"),
                        {"long_name": f"floodable land sub-cells below {BIN_LO} m"}),
            "n_above": (("cell",), np.asarray(rows_above, dtype="int32"),
                        {"long_name": "land sub-cells the sea cannot reach: above the LECZ "
                                      "bound, or below it inside a closed basin. Counted in "
                                      "the denominator, never exposed, percentile 1."}),
            "hist_expo": (("cell", "bin"), ehist,
                          {"long_name": "count of FLOODABLE land sub-cells per elevation "
                                        "bin -- connected to the sea at expo_level_m. Used "
                                        "for exposure only; the graded percentile uses "
                                        "`hist`."}),
            "n_below_expo": (("cell",), np.asarray(rows_ebelow, dtype="int32"),
                             {"long_name": f"floodable land sub-cells below {BIN_LO} m"}),
            "n_land": (("cell",), n_land,
                       {"long_name": "land sub-cells (TID == 0) in the cell"}),
            "n_unconnected": (("cell",), unconn,
                              {"long_name": "land sub-cells not hydrologically connected "
                                            "to the sea; counted but never exposed"}),
            "cell_flat_index": (("cell",), np.asarray(cell_ix, dtype="int32"),
                                {"long_name": f"flat index into the "
                                              f"{NLAT_OUT}x{NLON_OUT} half-degree grid"}),
        },
        coords={"bin_left": ("bin", edges[:-1])},
        attrs={
            "dem": "GEBCO_2026 ice-surface elevation, 15 arcsec; land base SRTM15+ v2.8",
            "dem_datum": "elevation standard_name is height_above_mean_sea_level",
            "dem_type": "DIGITAL SURFACE MODEL -- includes vegetation canopy and buildings; "
                        "SRTM-class coastal bias is +2.49 to +3.67 m, an order of magnitude "
                        "larger than the sea-level signal",
            "land_definition": "GEBCO_2026 TID == 0 ('Land'), NOT an elevation threshold",
            "connectivity": f"largest connected component of (TID != 0 OR elevation <= "
                            f"{CONNECT_LEVEL_M} m) at 1 arcmin, any-reduced from 15 arcsec "
                            f"so narrow straits stay open",
            "defences_represented": "NO -- dikes, levees and seawalls are narrower than a "
                                    "grid cell and absent from the DEM",
            "expo_level_m": EXPO_LEVEL_M,
            "expo_level_rationale": "connectivity for EXPOSURE is evaluated at 2 m, not at "
                                    "the 10 m LECZ bound: at 10 m the Salton Sink connects "
                                    "through the Colorado delta and the Danakil Depression "
                                    "through the Red Sea, so inland basins would read as "
                                    "maximum coastal risk. 2 m clears the largest measured "
                                    "sea-level delta (+1.328 m at rcp60).",
            "inland_buffer": f"{INLAND_RINGS} half-degree rings from a cell containing ocean",
            "bin_low_m": BIN_LO, "bin_high_m": BIN_HI, "bin_width_m": BIN_W,
        },
    )
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "hypsometry_halfdeg.nc"
    out.to_netcdf(p, encoding={"hist": {"zlib": True, "complevel": 4}})
    print(f"pass B: wrote {p} ({p.stat().st_size / 1e6:.0f} MB)")
    print(f"  cells with land: {int((n_land > 0).sum()):,}; land sub-cells: "
          f"{int(n_land.sum()):,}")
    print(f"  in LECZ pool (graded): {int(hist.sum() + np.asarray(rows_below).sum()):,}; "
          f"floodable (exposure):  {int(ehist.sum() + np.asarray(rows_ebelow).sum()):,}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true", help="print the GEBCO header and stop")
    a = ap.parse_args()
    elev_path, tid_path = find_grid()
    if a.inspect:
        inspect(elev_path)
        return
    ocean_c = pass_a_connectivity(elev_path, tid_path)
    pass_b_hypsometry(elev_path, tid_path, ocean_c)


if __name__ == "__main__":
    main()
