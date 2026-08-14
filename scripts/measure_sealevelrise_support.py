"""Where is the ISIMIP2b sea-level field actually defined, and can a coast reach it?

The first measurement found the field finite on only 62.5% of the grid, spanning BOTH
land and ocean, and 10 of 12 named coastal sites landing on empty cells. That makes the
support geometry a pipeline question, not a curiosity: if a coastal land cell has no
sea-level neighbour within a short distance, the delta applied to its elevation is
borrowed from somewhere else, and how far it was borrowed from is a caveat the customer
has to see.

Measures:
  1. the geography of the support (by latitude band, and against the land-sea mask)
  2. for every COASTAL land cell, the distance to the nearest defined sea-level cell
  3. the percentile axis in absolute terms (the first pass showed the DELTA is nearly
     percentile-independent, which is only interesting if the levels differ absolutely)
  4. the extreme values, which did not move monotonically with time

Run:  .venv/bin/python3 scripts/measure_sealevelrise_support.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

RAW = Path("data/raw/sealevelrise_2b")
SLR = RAW / "total_year_GFDL-ESM2M_rcp26_2006_2099.nc4"
MASK = RAW / "landseamask_generic.nc4"


def main():
    ds = xr.open_dataset(SLR, decode_times=False)
    da = ds["total"]
    lat = ds["lat"].values
    lon = ds["lon"].values
    pct = ds["percentile"].values

    field = da.isel(time=-1, percentile=3).values          # p50, 2099
    finite = np.isfinite(field)

    # THE MASK IS LATITUDE-FLIPPED RELATIVE TO THE SLR FILE. The sea-level file ascends
    # (-89.75 -> +89.75); ISIMIP2b_landseamask_generic descends (+89.75 -> -89.75).
    # Crossing them as raw arrays silently puts Antarctica's mask over the Arctic, and the
    # result looks plausible until you notice 411 "land" cells between 45N and 60N. Align
    # on the coordinate, never on array order. LSM is 1 over land and NaN over ocean --
    # there is no 0, so `> 0` alone would also have worked but `isfinite` is the real test.
    lsm = xr.open_dataset(MASK)["LSM"].squeeze().sortby("lat")
    assert np.allclose(lsm["lat"].values, lat), "mask lat does not align with SLR lat"
    assert np.allclose(lsm["lon"].values, lon), "mask lon does not align with SLR lon"
    land = np.isfinite(np.asarray(lsm.values))

    print("=" * 78)
    print("1. SUPPORT GEOGRAPHY")
    print("=" * 78)
    print(f"grid {field.shape}, finite {finite.sum():,} ({100*finite.mean():.1f}%)")
    print(f"\n{'lat band':>14}  {'cells':>7}  {'finite':>7}  {'%':>6}  "
          f"{'land':>7}  {'land&fin':>8}  {'land no-val':>11}")
    edges = [-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (lat >= lo) & (lat < hi)
        f = finite[sel]; l = land[sel]
        print(f"{lo:+4.0f}..{hi:+4.0f}      {f.size:7,}  {f.sum():7,}  "
              f"{100*f.mean():5.1f}%  {l.sum():7,}  {int((f & l).sum()):8,}  "
              f"{int((l & ~f).sum()):11,}")

    # Is the support a coastal band around the ocean, or something else?
    print("\n--- is the ocean support a coastal band? ---")
    ocean = ~land
    print(f"ocean cells {ocean.sum():,}; finite over ocean {int((ocean & finite).sum()):,} "
          f"({100*(ocean & finite).sum()/ocean.sum():.1f}%)")
    print(f"land  cells {land.sum():,}; finite over land  {int((land & finite).sum()):,} "
          f"({100*(land & finite).sum()/land.sum():.1f}%)")

    # --- 2. can each coastal land cell reach a value? ---------------------------------
    print("\n" + "=" * 78)
    print("2. COASTAL LAND CELLS -> DISTANCE TO NEAREST DEFINED SEA-LEVEL CELL")
    print("=" * 78)
    # A coastal land cell = land with at least one ocean cell in its 8-neighbourhood.
    nb = np.zeros_like(ocean, dtype=int)
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dj == 0 and di == 0:
                continue
            nb += np.roll(np.roll(ocean, dj, axis=0), di, axis=1).astype(int)
    coastal = land & (nb > 0)
    print(f"coastal land cells: {coastal.sum():,}")
    print(f"  of which already have a value : {int((coastal & finite).sum()):,} "
          f"({100*(coastal & finite).sum()/coastal.sum():.1f}%)")
    print(f"  of which need a lookup        : {int((coastal & ~finite).sum()):,}")

    jj, ii = np.where(finite)
    fin_lat = lat[jj]
    fin_lon = lon[ii]
    need_j, need_i = np.where(coastal & ~finite)
    dists = []
    # Great-circle distance, sampled -- exhaustive is O(n*m) and the distribution is what
    # matters, not each individual cell.
    rng = np.random.default_rng(0)
    sample = rng.choice(len(need_j), size=min(1500, len(need_j)), replace=False)
    for k in sample:
        la, lo_ = lat[need_j[k]], lon[need_i[k]]
        dlat = np.radians(fin_lat - la)
        dlon = np.radians(fin_lon - lo_)
        a = (np.sin(dlat / 2) ** 2
             + np.cos(np.radians(la)) * np.cos(np.radians(fin_lat)) * np.sin(dlon / 2) ** 2)
        dists.append(6371.0 * 2 * np.arcsin(np.sqrt(a)).min())
    dists = np.array(dists)
    print(f"\n  lookup distance over {len(sample)} sampled cells (km):")
    for q in (50, 75, 90, 95, 99, 100):
        print(f"    p{q:<3d} {np.percentile(dists, q):8.1f}")
    print(f"    mean {dists.mean():.1f}   >100 km: {100*(dists > 100).mean():.1f}%   "
          f">250 km: {100*(dists > 250).mean():.1f}%")

    # --- 3. percentile axis in absolute terms -----------------------------------------
    print("\n" + "=" * 78)
    print("3. PERCENTILE AXIS -- absolute values, 2090s mean")
    print("=" * 78)
    yrs = 1661 + ds["time"].values.astype(int)
    f9 = np.where((yrs >= 2090) & (yrs <= 2099))[0]
    b2 = np.where((yrs >= 2020) & (yrs <= 2029))[0]
    for i, p in enumerate(pct):
        v90 = da.isel(time=f9, percentile=i).mean("time").values
        v20 = da.isel(time=b2, percentile=i).mean("time").values
        print(f"  p{p:<5g} 2020s mean {np.nanmean(v20):+.4f}   2090s mean {np.nanmean(v90):+.4f}"
              f"   delta {np.nanmean(v90 - v20):+.4f} m")

    # --- 4. the extremes ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4. EXTREME VALUES -- 2006 max was +3.51 m, larger than 2099's")
    print("=" * 78)
    for t, label in ((0, "2006"), (93, "2099")):
        v = da.isel(time=t, percentile=3).values
        for fn, nm in ((np.nanargmax, "max"), (np.nanargmin, "min")):
            k = fn(np.where(np.isfinite(v), v, np.nan))
            j, i = np.unravel_index(k, v.shape)
            print(f"  {label} {nm}: {v[j, i]:+.3f} m at lat {lat[j]:+.2f} lon {lon[i]:+.2f} "
                  f"({'land' if land[j, i] else 'ocean'})")
    v = da.isel(percentile=3).values
    print(f"\n  |value| > 1 m anywhere in the record: "
          f"{int(np.nansum(np.abs(v) > 1)):,} cell-years of {int(np.isfinite(v).sum()):,}")

    ds.close()


if __name__ == "__main__":
    main()
