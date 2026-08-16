"""GUARDRAILS 9 + 12 value check for ISIMIP3b bias-adjusted daily precipitation.

Downloads ONE chunk per probed member, measures what the file actually contains, and
deletes it. Nothing here is inferred from a variable name or a sidecar.

WHY THE UNIT MUST BE MEASURED, NOT READ. ISIMIP publishes `pr` as a FLUX in kg m-2 s-1;
a day count wants mm/day, and the conversion is x86400 (1 kg m-2 == 1 mm of water). The
five group-I GCMs publish sidecars with NO netcdf_header at all (verified 2026-08-14 for
tasmax, re-verified 2026-08-16 for pr), so their declared unit is unreadable without
opening the file -- and a declared unit has been wrong before in this repository
(`burntarea` declares % on a 0-1 fraction). So the scale is inferred from VALUES against
physical reality: a global land mean near 2-3 mm/day is right, near 3e-5 means the flux
was not converted, and near 2.6e6 means it was converted twice.

WHY PRECIPITATION IS NOT A SECOND TEMPERATURE LADDER. Daily temperature varies about
two-fold across the globe in Kelvin and a fixed threshold means roughly the same thing
everywhere. Daily precipitation climatology varies ~100x, so an absolute mm/day threshold
is a different question in Manaus than in Phoenix. This script therefore measures the
per-cell distribution as well as the mean -- specifically each cell's wet-day count and
its high quantiles -- because those decide whether an absolute ladder, a local-percentile
ladder, or both are defensible for pluvial flooding.
"""

from __future__ import annotations

import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BASE = "https://files.isimip.org/ISIMIP3b"
GROUP_I = "InputData/climate/atmosphere/bias-adjusted/global/daily"
EXTENDED = "SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily"
STAGING = Path("data/raw/_pr_probe")
MASK = Path("data/masks/ISIMIP3b_landseamask_no-ant.nc")

#: One per root, plus the two GCMs that publish `pr` but NOT tasmax -- those are exactly
#: the members a 14-GCM decision turns on, so they are checked, not assumed comparable.
PROBES = [
    ("GFDL-ESM4", GROUP_I, "gfdl-esm4", "r1i1p1f1"),
    ("CNRM-ESM2-1", EXTENDED, "cnrm-esm2-1", "r1i1p1f2"),
    ("CESM2-WACCM", EXTENDED, "cesm2-waccm", "r1i1p1f1"),
    ("IITM-ESM", EXTENDED, "iitm-esm", "r1i1p1f1"),
]
SPAN = "2021_2030"
SCEN = "ssp585"

#: (name, lat, lon, what a working precipitation field must show)
SITES = [
    ("Manaus",       -3.12,  -60.02, "humid tropics: very wet, >2000 mm/yr"),
    ("Cherrapunji",  25.30,   91.70, "wettest inhabited place on earth"),
    ("Singapore",     1.35,  103.82, "equatorial: wet all year"),
    ("Phoenix AZ",   33.45, -112.07, "hot desert: ~200 mm/yr, few wet days"),
    ("Cairo",        30.04,   31.24, "hyper-arid: near zero, almost no wet days"),
    ("London",       51.51,   -0.13, "temperate maritime: ~600 mm/yr, MANY light days"),
    ("Seattle",      47.61, -122.33, "wet temperate, winter-dominated"),
    ("Mumbai",       19.08,   72.88, "monsoon: extreme daily intensities"),
]


def http(url, timeout=3600):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "curl/8.4"}), timeout=timeout)


def land_mask(lats):
    with Dataset(MASK) as ds:
        raw = ds.variables["mask"][:]
        a = np.asarray(raw.filled(np.nan) if np.ma.isMaskedArray(raw) else raw,
                       "f8").squeeze()
        mlat = np.asarray(ds.variables["lat"][:])
    m = np.isfinite(a) & (a > 0.5)
    if np.allclose(lats, mlat[::-1]):
        m = m[::-1]
    elif not np.allclose(lats, mlat):
        raise SystemExit("mask/grid mismatch")
    return m


def main():
    STAGING.mkdir(parents=True, exist_ok=True)
    for gcm, root, tok, ens in PROBES:
        stem = f"{tok}_{ens}_w5e5_{SCEN}_pr_global_daily_{SPAN}"
        url = f"{BASE}/{root}/{SCEN}/{gcm}/{stem}.nc"
        dest = STAGING / f"{stem}.nc"
        print("=" * 78)
        print(f"{gcm}  ({'group-I' if root == GROUP_I else 'extended'})")
        print("=" * 78)
        if not dest.exists():
            t0 = time.time()
            with http(url) as r, open(dest, "wb") as fh:
                while True:
                    b = r.read(1 << 22)
                    if not b:
                        break
                    fh.write(b)
            print(f"  downloaded {dest.stat().st_size/1e9:.2f} GB in "
                  f"{time.time()-t0:.0f}s")
        try:
            ds = Dataset(dest)
            v = ds.variables["pr"]
            t = ds.variables["time"]
            cal = getattr(t, "calendar", None) or "standard"
            yrs = np.array([d.year for d in num2date(t[:], t.units, calendar=cal)])
            lats = np.asarray(ds.variables["lat"][:])
            lons = np.asarray(ds.variables["lon"][:])
            print(f"  declared units : {getattr(v, 'units', 'UNDECLARED')!r}")
            print(f"  long_name      : {getattr(v, 'long_name', 'NONE')!r}")
            print(f"  calendar       : {cal}   days/yr {len(yrs)/len(set(yrs)):.2f}")

            m = land_mask(lats)
            # one full year, land only
            sel = yrs == sorted(set(yrs))[1]
            a = np.asarray(v[np.flatnonzero(sel)], "float64")[:, m]
            ocean = np.asarray(v[0], "float64")
            print(f"  finite over the WHOLE grid? "
                  f"{np.isfinite(ocean).mean():.1%} -- isfinite is "
                  f"{'NOT a mask' if np.isfinite(ocean).mean() > 0.99 else 'a mask'}")

            mmday = a * 86400.0
            gm = np.nanmean(mmday)
            scale = ("kg m-2 s-1 -> mm/day via x86400 (MEASURED: land mean "
                     f"{gm:.2f} mm/day)" if 0.5 < gm < 10 else
                     f"UNEXPECTED land mean {gm:.4g} mm/day AFTER x86400 -- do not proceed")
            print(f"  scale          : {scale}")
            print(f"  raw min/max    : {np.nanmin(a):.3e} / {np.nanmax(a):.3e} kg m-2 s-1")
            print(f"  mm/day max     : {np.nanmax(mmday):.1f}")
            neg = int((a < 0).sum())
            print(f"  negatives      : {neg}" + ("  <-- INVALID" if neg else ""))
            for thr in (0.1, 1.0, 10.0, 20.0, 50.0):
                print(f"    days >= {thr:5.1f} mm : "
                      f"{(mmday >= thr).sum(axis=0).mean():6.1f} per cell-year")
            wet = mmday >= 1.0
            wd = wet.sum(axis=0)
            print(f"  wet days (>=1mm): min {wd.min()}  median {np.median(wd):.0f}  "
                  f"max {wd.max()}  cells with ZERO wet days: "
                  f"{int((wd == 0).sum()):,} of {int(m.sum()):,}")
            q = np.array([np.percentile(c[c >= 1.0], 95) if (c >= 1.0).sum() > 10
                          else np.nan for c in mmday.T[::97]])
            print(f"  per-cell p95 of WET days (mm/day), sampled: "
                  f"min {np.nanmin(q):.1f}  median {np.nanmedian(q):.1f}  "
                  f"max {np.nanmax(q):.1f}  -- a {np.nanmax(q)/max(np.nanmin(q),1e-9):.0f}x "
                  f"range is why an absolute threshold is not a fixed question")

            print("  reference sites (mm/yr, and days >=20mm):")
            for name, la, lo in [(s[0], s[1], s[2]) for s in SITES]:
                i = int(np.argmin(np.abs(lats - la)))
                j = int(np.argmin(np.abs(lons - (lo % 360 if lons.max() > 180 else lo))))
                if not m[i, j]:
                    print(f"    {name:<13} OCEAN/masked")
                    continue
                col = np.asarray(v[np.flatnonzero(sel), i, j], "float64") * 86400.0
                exp = [s[3] for s in SITES if s[0] == name][0]
                print(f"    {name:<13} {col.sum():7.0f} mm/yr  "
                      f"{int((col >= 20).sum()):3d} d>=20mm   {exp}")
            ds.close()
        finally:
            dest.unlink(missing_ok=True)
    print("\nprobe files deleted.")


if __name__ == "__main__":
    sys.exit(main())
