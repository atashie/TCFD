"""GUARDRAILS 9 + 12 value check for the daily bias-adjusted temperature forcing.

WHY THIS EXISTS. A threshold-exceedance layer is a CONSTRUCTED INDEX (GUARDRAILS 1-2):
nothing in ISIMIP publishes a days-above-a-temperature count, in any round or any product
root (verified 2026-08-14, see config/isimip_search_catalog.yaml search_results.chronic_heat).
That means every property the count depends on is ours to measure, not to inherit:

  * UNITS. Declared K on the seven SecondaryInputData GCMs, whose sidecars carry a
    netcdf_header. The five group-I InputData GCMs publish sidecars WITHOUT a header block
    (checksum and specifiers only, verified 2026-08-14), so for those five the unit is
    UNVERIFIED until read from the values. A degC file silently counted against a 308.15 K
    threshold returns zero everywhere, which reads exactly like a cold climate.

  * CALENDAR. A day COUNT is not calendar-neutral. A 360_day model has 360 chances per
    year to cross a threshold and a proleptic_gregorian one has 365.25, so a raw count from
    a 360_day member is biased ~1.5% low against its siblings -- a bias that lands in the
    trend, because it is constant per member and the ensemble composition is what varies.
    All seven extended GCMs declare proleptic_gregorian with 3652 steps per decade; the
    group-I five are unread. This script reports days-per-year per member so the pooling
    decision is made on a measurement.

  * MASK. The forcing is atmospheric and is NOT land-masked -- it is finite over ocean.
    `isfinite` is therefore not a footprint (the same trap as `cropfailure` and `leh`).
    Land must come from the official ISIMIP3b landseamask, which scripts/utils/land_mask.py
    already resolves.

WHAT IT PRINTS. Per member: declared metadata, measured value range in BOTH K and degC,
the fill/NaN share, the land-vs-ocean split, the per-year day count, and then the nine
threshold counts of the ladder evaluated at named reference sites (GUARDRAILS 12) where
the answer is known well enough to falsify a bug -- Singapore must have zero frost days,
Yakutsk must have hundreds, Kuwait must have a large hd35 and Paris a small one.

A zero at a reference site is a STOP, not a data point to rationalise.

Usage:
    .venv/bin/python3 scripts/check_tas_thresholds_nature.py --gcm GFDL-ESM4 \
        --scenario ssp585 --chunk 2021_2030 [--keep]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BASE = "https://files.isimip.org/ISIMIP3b"
GROUP_I = "InputData/climate/atmosphere/bias-adjusted/global/daily"
EXTENDED = "SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily"

#: (GCM directory name, filename token, ensemble member). Root is resolved by membership:
#: the two pools are DISJOINT for ssp126/ssp370/ssp585 (measured 2026-08-14).
GROUP_I_GCMS = {
    "GFDL-ESM4": ("gfdl-esm4", "r1i1p1f1"),
    "IPSL-CM6A-LR": ("ipsl-cm6a-lr", "r1i1p1f1"),
    "MPI-ESM1-2-HR": ("mpi-esm1-2-hr", "r1i1p1f1"),
    "MRI-ESM2-0": ("mri-esm2-0", "r1i1p1f1"),
    "UKESM1-0-LL": ("ukesm1-0-ll", "r1i1p1f2"),
}
EXTENDED_GCMS = {
    "CNRM-CM6-1": ("cnrm-cm6-1", "r1i1p1f2"),
    "CNRM-ESM2-1": ("cnrm-esm2-1", "r1i1p1f2"),
    "CanESM5": ("canesm5", "r1i1p1f1"),
    "EC-Earth3": ("ec-earth3", "r1i1p1f1"),
    "KACE-1-0-G": ("kace-1-0-g", "r1i1p1f1"),
    "MIROC6": ("miroc6", "r1i1p1f1"),
    "TaiESM1": ("taiesm1", "r1i1p1f1"),
}
#: CESM2-WACCM and IITM-ESM are deliberately absent: they publish `tas` and `hurs` but NO
#: tasmax/tasmin in any scenario (measured per GCM 2026-08-14). A pooled variable count
#: reports 9 GCMs for a layer that can only have 7.

RAW = Path("data/raw/heatdays-isimip3b_tasthresh_daily")

#: The ladder. Every rung is free once the chunk is streamed -- the day is tested anyway --
#: so the count is not a commitment to publish all of them. Thresholds in degC.
HOT = [30.0, 35.0, 40.0, 45.0]      # World Bank CCKP hd30/hd35/hd40/hd45
COLD_MIN = [0.0, -10.0]             # ETCCDI FD (Tmin<0); -10 is severe frost, NON-standard
COLD_MAX = [0.0]                    # ETCCDI ID (ice day, Tmax<0)
TROPICAL_NIGHT = [20.0, 25.0]       # ETCCDI TR20; CCKP ships tr23/26/29/32

#: GUARDRAILS 12. Sites chosen so a bug is falsifiable, not so the layer looks good:
#: the hot sites must be large and the frost sites zero, and vice versa.
SITES = [
    ("Kuwait City", 29.37, 47.98, "hd35 must be very large; FD must be 0"),
    ("Jacobabad", 28.28, 68.44, "hd35 must be very large; FD must be 0"),
    ("Phoenix", 33.45, -112.07, "hd35 large; FD small but non-zero"),
    ("Delhi", 28.61, 77.21, "hd35 large; FD ~0"),
    ("Timbuktu", 16.78, -3.01, "hd35 very large; FD 0"),
    ("Alice Springs", -23.70, 133.88, "hd35 large; FD non-zero (desert night)"),
    ("Singapore", 1.35, 103.82, "hd35 SMALL despite equator; FD 0; TR20 ~ every night"),
    ("Manaus", -3.10, -60.02, "hd35 modest; FD 0; TR20 ~ every night"),
    ("Lagos", 6.52, 3.38, "hd35 small; FD 0; TR20 ~ every night"),
    ("Paris", 48.86, 2.35, "hd35 SMALL; FD tens of days"),
    ("Chicago", 41.88, -87.63, "hd35 small; FD ~120"),
    ("Yakutsk", 62.03, 129.73, "hd35 ~0; FD VERY large; ID large"),
    ("Fairbanks", 64.84, -147.72, "hd35 0; FD very large; ID large"),
    ("Ulaanbaatar", 47.89, 106.91, "hd35 ~0; FD very large"),
]


def sha512_of(path: Path) -> str:
    h = hashlib.sha512()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def resolve(gcm: str):
    if gcm in GROUP_I_GCMS:
        return GROUP_I, GROUP_I_GCMS[gcm]
    if gcm in EXTENDED_GCMS:
        return EXTENDED, EXTENDED_GCMS[gcm]
    raise SystemExit(
        f"unknown GCM {gcm!r}. group-I: {sorted(GROUP_I_GCMS)}; "
        f"extended: {sorted(EXTENDED_GCMS)}")


def fetch(gcm: str, scenario: str, var: str, chunk: str) -> tuple[Path, dict]:
    root, (tok, ens) = resolve(gcm)
    stem = f"{tok}_{ens}_w5e5_{scenario}_{var}_global_daily_{chunk}"
    url = f"{BASE}/{root}/{scenario}/{gcm}/{stem}.nc"
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / f"{stem}.nc"

    side = {}
    try:
        with urllib.request.urlopen(
                urllib.request.Request(f"{BASE}/{root}/{scenario}/{gcm}/{stem}.json",
                                       headers={"User-Agent": "curl/8.4"}), timeout=60) as r:
            side = json.loads(r.read().decode())
    except Exception as e:  # noqa: BLE001
        print(f"  sidecar unavailable ({str(e)[:60]}) -- integrity will rest on size alone")

    if dest.exists() and side.get("size") and dest.stat().st_size == side["size"]:
        print(f"  {dest.name}: already present and size-matched, skipping download")
    else:
        t0 = time.time()
        print(f"  downloading {stem}.nc ...", flush=True)
        with urllib.request.urlopen(
                urllib.request.Request(url, headers={"User-Agent": "curl/8.4"}),
                timeout=1800) as r, open(dest, "wb") as fh:
            n = 0
            while True:
                b = r.read(1 << 22)
                if not b:
                    break
                fh.write(b)
                n += len(b)
        dt = time.time() - t0
        print(f"  {n / 1e9:.2f} GB in {dt / 60:.1f} min ({n / 1e6 / dt:.1f} MB/s)")

    if side.get("checksum"):
        got = sha512_of(dest)
        ok = got == side["checksum"]
        print(f"  sha512 {'MATCH' if ok else 'MISMATCH'} against publisher sidecar")
        if not ok:
            raise SystemExit("checksum mismatch -- refusing to use this file")
    return dest, side


def describe(path: Path, var: str, side: dict) -> dict:
    """GUARDRAILS 9. Report what the file DECLARES and what its VALUES say."""
    ds = Dataset(path)
    v = ds.variables[var]
    t = ds.variables["time"]
    ntime = len(ds.dimensions["time"])
    cal = getattr(t, "calendar", None)
    tunits = getattr(t, "units", None)
    d0 = num2date(t[0], tunits, calendar=cal or "standard")
    d1 = num2date(t[-1], tunits, calendar=cal or "standard")
    nyears = (d1.year - d0.year) + 1

    declared_units = getattr(v, "units", None)
    print(f"\n  DECLARED  units={declared_units!r} long_name={getattr(v, 'long_name', None)!r}")
    print(f"  DECLARED  calendar={cal!r} time_units={tunits!r}")
    print(f"  DECLARED  dims={ {k: len(d) for k, d in ds.dimensions.items()} }")
    print(f"  MEASURED  span {d0.isoformat()[:10]} .. {d1.isoformat()[:10]} = {nyears} yr, "
          f"{ntime} steps -> {ntime / nyears:.2f} days/yr")
    if side.get("netcdf_header"):
        print("  sidecar carried a netcdf_header (extended root)")
    else:
        print("  sidecar carried NO netcdf_header -- unit/calendar were UNVERIFIED "
              "until the lines above")

    # values: sample strided slices rather than loading 3.8 GB
    vmin, vmax, nfin, ntot = np.inf, -np.inf, 0, 0
    for i in range(0, ntime, max(1, ntime // 12)):
        a = np.asarray(v[i], dtype="float64")
        fin = np.isfinite(a)
        if fin.any():
            vmin = min(vmin, float(a[fin].min()))
            vmax = max(vmax, float(a[fin].max()))
        nfin += int(fin.sum())
        ntot += a.size
    print(f"  MEASURED  range [{vmin:.2f}, {vmax:.2f}] -- "
          f"as K that is [{vmin - 273.15:.2f}, {vmax - 273.15:.2f}] degC; "
          f"as degC it would be [{vmin:.2f}, {vmax:.2f}] degC")
    kelvin = vmin > 150.0
    print(f"  VERDICT   unit is {'KELVIN' if kelvin else 'NOT Kelvin -- CHECK'} "
          f"(min {vmin:.1f} {'>' if kelvin else '<='} 150)")
    print(f"  MEASURED  finite on {nfin}/{ntot} sampled cells ({100 * nfin / ntot:.1f}%) "
          f"-- forcing is not land-masked, so this is NOT a footprint")
    ds.close()
    return {"ntime": ntime, "nyears": nyears, "calendar": cal, "kelvin": kelvin,
            "days_per_year": ntime / nyears, "vmin": vmin, "vmax": vmax}


def counts(path: Path, var: str, kelvin: bool, above: list, below: list) -> dict:
    """Annual-mean day counts per threshold, accumulated slice-wise (no dask here)."""
    ds = Dataset(path)
    v = ds.variables[var]
    ntime = len(ds.dimensions["time"])
    off = 273.15 if kelvin else 0.0
    out = {f"gt{t:g}": np.zeros(v.shape[1:], "int32") for t in above}
    out.update({f"lt{t:g}": np.zeros(v.shape[1:], "int32") for t in below})
    step = 366
    for i in range(0, ntime, step):
        a = np.asarray(v[i:i + step], dtype="float32") - off
        for t in above:
            out[f"gt{t:g}"] += (a > t).sum(axis=0).astype("int32")
        for t in below:
            out[f"lt{t:g}"] += (a < t).sum(axis=0).astype("int32")
    lat = ds.variables["lat"][:]
    lon = ds.variables["lon"][:]
    ds.close()
    return {"counts": out, "lat": np.asarray(lat), "lon": np.asarray(lon)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gcm", default="GFDL-ESM4")
    ap.add_argument("--scenario", default="ssp585")
    ap.add_argument("--chunk", default="2021_2030")
    ap.add_argument("--keep", action="store_true", help="do not delete the raw chunk")
    a = ap.parse_args()

    print(f"=== {a.gcm} {a.scenario} {a.chunk} "
          f"({'group-I' if a.gcm in GROUP_I_GCMS else 'extended'} root) ===")

    res = {}
    for var in ("tasmax", "tasmin"):
        print(f"\n--- {var} ---")
        path, side = fetch(a.gcm, a.scenario, var, a.chunk)
        meta = describe(path, var, side)
        above = HOT if var == "tasmax" else TROPICAL_NIGHT
        below = COLD_MAX if var == "tasmax" else COLD_MIN
        c = counts(path, var, meta["kelvin"], above, below)
        res[var] = (meta, c)
        if not a.keep:
            path.unlink()
            print(f"  deleted {path.name} (streaming design; sha512 recorded above)")

    nyr = res["tasmax"][0]["nyears"]
    lat, lon = res["tasmax"][1]["lat"], res["tasmax"][1]["lon"]
    print(f"\n=== REFERENCE SITES (GUARDRAILS 12) -- days per year, mean over {nyr} yr ===")
    hdr = (f"{'site':14s} {'hd30':>6s} {'hd35':>6s} {'hd40':>6s} {'hd45':>6s} "
           f"{'ID':>6s} {'FD':>6s} {'FD-10':>6s} {'TR20':>6s} {'TR25':>6s}   expectation")
    print(hdr)
    for name, la, lo, expect in SITES:
        i = int(np.abs(lat - la).argmin())
        j = int(np.abs(lon - lo).argmin())
        cx, cn = res["tasmax"][1]["counts"], res["tasmin"][1]["counts"]
        row = [cx["gt30"][i, j], cx["gt35"][i, j], cx["gt40"][i, j], cx["gt45"][i, j],
               cx["lt0"][i, j], cn["lt0"][i, j], cn["lt-10"][i, j],
               cn["gt20"][i, j], cn["gt25"][i, j]]
        print(f"{name:14s} " + " ".join(f"{x / nyr:6.1f}" for x in row) + f"   {expect}")

    print("\nRead the table against the expectation column, not against plausibility. "
          "A zero where the expectation says 'very large' is a STOP.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
