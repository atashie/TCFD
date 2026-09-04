#!/usr/bin/env python
"""Member-resolved dis panels for the WRI basin runner's paired ledger.

The runner's review (WRI v2/reviews/2026-09-04, finding 9) established
that the ensemble-mean family files cannot satisfy the B3 member-pairing
contract (same impact model + GCM on supply and demand, statistics after
pairing). This builder produces the missing supply side for the PAIRABLE
roster — the two impact models whose demand stacks are complete from
their own variables on both bases:

    miroc-integ-land, watergap2-2e   x  5 GCMs  x  3 scenarios

Reduction is the FAMILY CONVENTION applied per member (mirrored from
process_water_variable.py, asserted here in comments):
  VT 0-11   mean of that calendar month's instances within the decade
  VT 12     annual = unweighted mean of VT 0-11
  VT 13/14  Q15 / Q25 = np.nanpercentile (linear) of the member's
            per-year annual means (each year = unweighted mean of its
            12 monthly values) within the decade
Decades 2010..2090 with the 2010s = 2015-2019 stub (config convention).
Values stay in the native m3/s; the runner converts to volumes.

Trust chain per raw file (the process_role1_demand pattern): sole 3-dim
data variable named 'dis', units 'm3 s-1', 1,032 monthly steps decoding
to 2015-01..2100-12 under the file's own calendar, family grid
orientation (lat 89.75 descending), land-cell count in [67000, 67500],
no negative discharge beyond -1e-6. VT12 == mean(VT0-11) is exact by
construction and asserted. sha256 of every raw file lands in provenance.

Output: data/processed/wri_runner/dis_member_panels_v1.nc
  dims (member 10, scenario 3, decade 9, value_type 15, lat, lon), f4.

Run:  cd ~/github/TCFD && .venv/bin/python scripts/build_dis_member_panels.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import netCDF4
import numpy as np

TCFD = Path(__file__).resolve().parents[1]
RAW = TCFD / "data" / "raw" / "water_dis"
OUT_DIR = TCFD / "data" / "processed" / "wri_runner"

MODELS = ("miroc-integ-land", "watergap2-2e")
GCMS = ("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0",
        "ukesm1-0-ll")
SCENARIOS = ("ssp126", "ssp370", "ssp585")
DECADES = tuple(range(2010, 2100, 10))
N_LAT, N_LON, N_VT = 360, 720, 15
VT_NAMES = {**{i: f"month_{i+1:02d}_mean" for i in range(12)},
            12: "annual_mean", 13: "annual_Q15", 14: "annual_Q25"}


def fail(msg: str) -> None:
    print(f"\nDIS MEMBER PANELS REFUSED: {msg}", file=sys.stderr)
    sys.exit(1)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def decade_years(decade: int) -> tuple:
    return (2015, 2019) if decade == 2010 else (decade, decade + 9)


def reduce_member(path: Path) -> tuple:
    """(panels[3? no — per file is one scenario] -> [9, 15, lat, lon],
    receipts dict) for one member-scenario raw file."""
    with netCDF4.Dataset(path) as ds:
        data_vars = [k for k, v in ds.variables.items()
                     if len(v.dimensions) == 3]
        if data_vars != ["dis"]:
            fail(f"{path.name}: data variables {data_vars} != ['dis']")
        v = ds.variables["dis"]
        if getattr(v, "units", "") != "m3 s-1":
            fail(f"{path.name}: units {getattr(v, 'units', '?')}")
        lat = np.asarray(ds.variables["lat"][:])
        if not np.allclose(lat, 89.75 - 0.5 * np.arange(N_LAT), atol=1e-6):
            fail(f"{path.name}: not family lat orientation")
        tvar = ds.variables["time"]
        if len(tvar) != 1032:
            fail(f"{path.name}: {len(tvar)} time steps != 1032")
        first = netCDF4.num2date(tvar[0], tvar.units,
                                 getattr(tvar, "calendar", "standard"))
        last = netCDF4.num2date(tvar[-1], tvar.units,
                                getattr(tvar, "calendar", "standard"))
        if (first.year, first.month) != (2015, 1) or \
                (last.year, last.month) != (2100, 12):
            fail(f"{path.name}: time axis {first}..{last}")
        arr = np.ma.filled(v[:].astype(np.float32), np.nan)

    land = np.isfinite(arr[0])
    n_land = int(land.sum())
    if not (67000 <= n_land <= 67500):
        fail(f"{path.name}: land-cell count {n_land}")
    vmin = float(np.nanmin(arr))
    if vmin < -1e-6:
        fail(f"{path.name}: negative discharge {vmin}")
    years = arr.reshape(86, 12, N_LAT, N_LON)     # 2015..2100
    out = np.full((len(DECADES), N_VT, N_LAT, N_LON), np.nan,
                  dtype=np.float32)
    for d_i, dec in enumerate(DECADES):
        y0, y1 = decade_years(dec)
        sl = years[y0 - 2015:y1 - 2015 + 1]       # (n_years, 12, lat, lon)
        with np.errstate(invalid="ignore"):
            monthly = np.nanmean(sl, axis=0)      # family VT 0-11 rule
            out[d_i, :12] = monthly
            out[d_i, 12] = np.nanmean(monthly, axis=0)   # VT12 = mean(0-11)
            annual = np.nanmean(sl, axis=1)       # per-year annual means
            out[d_i, 13] = np.nanpercentile(annual, 15, axis=0)
            out[d_i, 14] = np.nanpercentile(annual, 25, axis=0)
        dev = np.nanmax(np.abs(out[d_i, 12] - np.nanmean(out[d_i, :12],
                                                         axis=0)))
        if not (dev <= 1e-4 or np.isnan(dev)):
            fail(f"{path.name}: VT12 identity dev {dev}")
    receipts = {"n_land": n_land, "min_m3s": vmin,
                "global_mean_2020s_m3s": round(float(
                    np.nansum(out[1, 12])), 1)}
    return out, receipts


def main() -> None:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "dis_member_panels_v1.nc"
    members = [f"{m}/{g}" for m in MODELS for g in GCMS]

    onc = netCDF4.Dataset(out_path, "w")
    onc.createDimension("member", len(members))
    onc.createDimension("scenario", len(SCENARIOS))
    onc.createDimension("decade", len(DECADES))
    onc.createDimension("value_type", N_VT)
    onc.createDimension("lat", N_LAT)
    onc.createDimension("lon", N_LON)
    onc.createVariable("lat", "f8", ("lat",))[:] = \
        89.75 - 0.5 * np.arange(N_LAT)
    onc.createVariable("lon", "f8", ("lon",))[:] = \
        -179.75 + 0.5 * np.arange(N_LON)
    onc.createVariable("decade", "i4", ("decade",))[:] = DECADES
    v_m = onc.createVariable("member", str, ("member",))
    for i, m in enumerate(members):
        v_m[i] = m
    v_s = onc.createVariable("scenario", str, ("scenario",))
    for i, s in enumerate(SCENARIOS):
        v_s[i] = s
    v_vt = onc.createVariable("value_type", "i4", ("value_type",))
    v_vt[:] = np.arange(N_VT)
    v_vt.description = json.dumps(VT_NAMES)
    v_d = onc.createVariable("dis", "f4",
                             ("member", "scenario", "decade",
                              "value_type", "lat", "lon"),
                             zlib=True, complevel=4,
                             fill_value=np.float32(np.nan))
    v_d.units = "m3 s-1"

    hashes = {}
    receipts = {}
    for mi, member in enumerate(members):
        model, gcm = member.split("/")
        for si, scen in enumerate(SCENARIOS):
            d = RAW / model / f"{gcm}_{scen}"
            files = sorted(d.glob("*_dis_global_monthly_2015_2100.nc"))
            if len(files) != 1:
                fail(f"{d}: {len(files)} candidate files")
            path = files[0]
            panels, rec = reduce_member(path)
            v_d[mi, si] = panels
            hashes[path.name] = sha256_file(path)
            receipts[f"{member}/{scen}"] = rec
            print(f"  {member} {scen}: land {rec['n_land']:,} | "
                  f"2020s global-mean sum {rec['global_mean_2020s_m3s']:,}"
                  f" m3/s | {time.time() - t0:.0f}s")

    def git(*a):
        return subprocess.run(["git", "-C", str(TCFD), *a],
                              capture_output=True, text=True).stdout.strip()
    prov = {"generator": "scripts/build_dis_member_panels.py",
            "self_sha256": sha256_file(Path(__file__)),
            "git_commit": git("rev-parse", "HEAD"),
            "git_dirty": bool(git("status", "--porcelain")),
            "raw_sha256": hashes,
            "reduction": "family convention per member (mirrors "
                         "process_water_variable.py): VT0-11 monthly "
                         "means in decade; VT12 = mean(VT0-11); "
                         "VT13/14 = nanpercentile(15/25, linear) of "
                         "per-year annual means; 2010s = 2015-2019",
            "purpose": "paired-ledger supply side (WRI runner review "
                       "finding 9); pairable roster only",
            "run_utc": datetime.now(timezone.utc).isoformat(
                timespec="seconds")}
    onc.title = "Member-resolved dis panels (paired-ledger roster)"
    onc.provenance = json.dumps(prov)
    onc.close()
    (OUT_DIR / "dis_member_panels_v1.receipts.json").write_text(
        json.dumps({"receipts": receipts, "provenance": prov}, indent=2))
    print(f"\nBUILT {out_path}\nwall {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
