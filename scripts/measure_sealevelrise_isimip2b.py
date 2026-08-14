"""Measure the ISIMIP2b regional sea-level-rise input before framing anything on it.

GUARDRAILS 9: data nature is MEASURED, never inferred from a name, a long_name or a
sibling. This product declares `units = m` and a 7-level `percentile` axis whose
long_name is literally "generic", so the axis contents, the spatial support (ocean? land?
both?) and the magnitude of the signal are all unknown until read.

The questions that decide whether a coastal-exposure layer can be built on it at all:

  1. What are the 7 percentile levels, and is the axis ordered?
  2. Where is the field DEFINED -- if it is ocean-only, every coastal land cell needs a
     nearest-ocean lookup, and the cost of that lookup is a pipeline stage, not a detail.
  3. Is it an ANOMALY (starts near 0 in 2006) or an absolute height above a geoid? The
     delta-based design only works if we know which.
  4. How big is the 2020s -> 2090s delta, and does it vary regionally? A spatially flat
     field would mean the 0.5 deg grid carries no information a global scalar wouldn't.
  5. Does every coastal land cell actually have a defined sea-level neighbour?

Run:  .venv/bin/python3 scripts/measure_sealevelrise_isimip2b.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

RAW = Path("data/raw/sealevelrise_2b")
FILES = {
    "rcp26": RAW / "total_year_GFDL-ESM2M_rcp26_2006_2099.nc4",
    "rcp60": RAW / "total_year_GFDL-ESM2M_rcp60_2006_2099.nc4",
}
MASK = RAW / "landseamask_generic.nc4"

# Reference coasts, (lat, lon) in degrees east. Deliberately spread across the regional
# fingerprint pattern: SLR is NOT globally uniform, and a layer that cannot reproduce
# that is not carrying the information the 0.5 deg grid claims to.
SITES = {
    "Ganges-Brahmaputra delta": (22.5, 90.5),
    "Mekong delta": (10.0, 106.0),
    "Nile delta": (31.5, 31.0),
    "Netherlands coast": (52.5, 4.5),
    "Miami / SE Florida": (25.5, -80.5),
    "New Orleans / Mississippi": (29.5, -90.0),
    "Jakarta": (-6.0, 106.5),
    "Lagos": (6.5, 3.5),
    "Shanghai": (31.0, 121.5),
    "Alaska (Anchorage)": (61.0, -150.0),
    "Iceland (Reykjavik)": (64.0, -22.0),
    "Tuvalu / C Pacific": (-8.5, 179.0),
}


def wrap_lon(lon, ref):
    """Match the file's longitude convention (0..360 vs -180..180)."""
    return lon % 360 if ref.max() > 180 else ((lon + 180) % 360) - 180


def main():
    print("=" * 78)
    print("ISIMIP2b sealevelrise -- measured characteristics")
    print("=" * 78)

    ds = xr.open_dataset(FILES["rcp26"], decode_times=False)
    print(f"\nvariables: {list(ds.data_vars)}")
    var = "total"
    da = ds[var]
    print(f"dims       : {dict(da.sizes)}")
    print(f"attrs      : {da.attrs}")

    # --- 1. the percentile axis -------------------------------------------------------
    print("\n--- 1. percentile axis ---")
    pct = ds["percentile"].values
    print(f"values     : {pct}")
    print(f"attrs      : {ds['percentile'].attrs}")

    # --- 2. spatial support -----------------------------------------------------------
    print("\n--- 2. where is the field defined? ---")
    last = da.isel(time=-1)
    # Use a mid percentile level for the support test; NaN structure is level-independent
    # only if measured, so check that too.
    finite_by_level = [int(np.isfinite(last.isel(percentile=i).values).sum()) for i in range(len(pct))]
    print(f"finite cells per percentile level: {finite_by_level}")

    field = last.isel(percentile=len(pct) // 2).values
    finite = np.isfinite(field)
    total_cells = field.size
    print(f"finite cells: {finite.sum():,} of {total_cells:,} ({100*finite.sum()/total_cells:.1f}%)")

    lsm = xr.open_dataset(MASK)
    lsm_var = [v for v in lsm.data_vars if lsm[v].ndim >= 2][0]
    land = np.asarray(lsm[lsm_var].squeeze().values)
    land_bool = np.isfinite(land) & (land > 0)
    print(f"landseamask var={lsm_var!r} -> land cells: {land_bool.sum():,}")

    if land_bool.shape == finite.shape:
        print(f"  finite AND land : {int((finite & land_bool).sum()):,}")
        print(f"  finite AND ocean: {int((finite & ~land_bool).sum()):,}")
        print(f"  land WITHOUT a value: {int((land_bool & ~finite).sum()):,}")
    else:
        print(f"  !! shape mismatch mask{land_bool.shape} vs field{finite.shape} -- cannot cross")

    # --- 3. anomaly or absolute? ------------------------------------------------------
    print("\n--- 3. anomaly or absolute height? ---")
    first = da.isel(time=0, percentile=len(pct) // 2).values
    ok = np.isfinite(first)
    print(f"2006 field: min {np.nanmin(first):+.4f}  mean {np.nanmean(first):+.4f}  "
          f"max {np.nanmax(first):+.4f} m  (n={ok.sum():,})")
    print(f"2099 field: min {np.nanmin(field):+.4f}  mean {np.nanmean(field):+.4f}  "
          f"max {np.nanmax(field):+.4f} m")

    # --- 4. magnitude and regional structure ------------------------------------------
    print("\n--- 4. 2020s -> 2090s delta (per percentile level) ---")
    yrs = ds["time"].values
    # time is "years since 1661-1-1" on a 365_day calendar -> integer year offsets
    years = 1661 + yrs.astype(int)
    print(f"year range : {years.min()}-{years.max()} ({len(years)} steps)")
    b = (years >= 2020) & (years <= 2029)
    f = (years >= 2090) & (years <= 2099)
    print(f"baseline 2020s: {b.sum()} yrs, future 2090s: {f.sum()} yrs")

    for scen, path in FILES.items():
        d = xr.open_dataset(path, decode_times=False)[var]
        base = d.isel(time=np.where(b)[0]).mean("time")
        fut = d.isel(time=np.where(f)[0]).mean("time")
        delta = (fut - base)
        row = []
        for i in range(len(pct)):
            v = delta.isel(percentile=i).values
            row.append(f"p{pct[i]:g}={np.nanmean(v):+.3f}")
        print(f"  {scen} global-mean delta (m): {'  '.join(row)}")
        mid = delta.isel(percentile=len(pct) // 2).values
        print(f"    spatial spread at mid level: min {np.nanmin(mid):+.3f}  "
              f"max {np.nanmax(mid):+.3f}  sd {np.nanstd(mid):.3f} m")
        d.close()

    # --- 5. reference sites ------------------------------------------------------------
    print("\n--- 5. reference coasts, 2020s -> 2090s delta at the mid percentile (m) ---")
    lat_c, lon_c = ds["lat"].values, ds["lon"].values
    for scen, path in FILES.items():
        d = xr.open_dataset(path, decode_times=False)[var]
        base = d.isel(time=np.where(b)[0]).mean("time").isel(percentile=len(pct) // 2)
        fut = d.isel(time=np.where(f)[0]).mean("time").isel(percentile=len(pct) // 2)
        delta = (fut - base).values
        print(f"  [{scen}]")
        for name, (la, lo) in SITES.items():
            j = int(np.abs(lat_c - la).argmin())
            i = int(np.abs(lon_c - wrap_lon(lo, lon_c)).argmin())
            v = delta[j, i]
            # If the target cell is empty, report how far the nearest defined cell is --
            # that distance IS the coastal-lookup problem, quantified.
            if not np.isfinite(v):
                fin = np.isfinite(delta)
                jj, ii = np.where(fin)
                dist = np.hypot((jj - j) * 0.5, (ii - i) * 0.5 * np.cos(np.radians(la)))
                k = int(dist.argmin())
                print(f"    {name:28s} EMPTY at cell; nearest defined {dist[k]:.2f} deg away "
                      f"-> {delta[jj[k], ii[k]]:+.3f}")
            else:
                print(f"    {name:28s} {v:+.3f}")
        d.close()

    ds.close()
    lsm.close()


if __name__ == "__main__":
    main()
