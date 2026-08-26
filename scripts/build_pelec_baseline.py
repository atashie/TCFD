#!/usr/bin/env python
"""Build the static 2020s thermoelectric water-demand baseline from ISIMIP3b.

User decision 2026-08-26: process and store pelecww / pelecuse as static values
for the 2020s, to serve as the baseline against which future scenario and
mitigation analyses are run.

WHAT THESE DATA ARE — and are not (measured, see
waterRiskIndex_beta/v2/enumeration/WS2_role3_assessment.md):

  ARE : a credible present-day SPATIAL MAP of potential thermoelectric water
        withdrawal and consumption. Contiguous-US withdrawal totals 186.3 km3/yr
        against USGS 2015's ~184 km3/yr — agreement within ~1%.
  NOT : a projection. Under 2015soc the field is CONSTANT: Jan 2015 == Jan 2100,
        the seasonal max/min ratio is exactly 1.000, and all 15 GCM x SSP files
        are byte-identical. ISIMIP3b publishes no ssp*soc variant for pelec, so
        no transient thermoelectric pathway exists at all.

Because of that, "the 2020s baseline" is the 2020-2029 mean, and the script
VERIFIES that this equals the all-period mean rather than assuming it. If a
future ISIMIP release makes pelec transient, that assertion fails loudly instead
of silently averaging away a real signal.

Outputs a single small NetCDF holding, per variable:
    <var>          flux, kg m-2 s-1  (native units, as published)
    <var>_km3yr    volume, km3/yr per cell (flux x cell area x 365.25 d)
The volume field makes subwatershed aggregation a plain sum when the BasinATLAS
publish fabric is wired.

Usage:
    python build_pelec_baseline.py [--dry-run]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

RAW = Path("/Users/arik/github/TCFD/data/raw/ws2_role3/pelec")
OUT_DIR = Path("/Users/arik/github/TCFD/data/processed/ws2_role3")
OUT = OUT_DIR / "pelec_baseline_2020s.nc"
VARS = ["pelecww", "pelecuse"]
SEC_PER_YEAR = 365.25 * 86400.0
BASELINE_YEARS = (2020, 2029)
FILE_START_YEAR = 2015  # monthly files run 2015-01 .. 2100-12
FILL = 1e19             # values >= this are the netCDF fill, not data


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def cell_area_m2(lat: np.ndarray, n_lon: int, dlat: float = 0.5, dlon: float = 0.5) -> np.ndarray:
    """Spherical-band cell area, (lat, lon) in m^2."""
    r = 6_371_000.0
    la = np.deg2rad(lat)
    band = (r ** 2) * np.deg2rad(dlon) * (
        np.sin(la + np.deg2rad(dlat / 2)) - np.sin(la - np.deg2rad(dlat / 2)))
    return np.repeat(np.abs(band)[:, None], n_lon, axis=1)


def load_clean(path: str, var: str, sl: slice | None = None) -> np.ndarray:
    nc = Dataset(path)
    try:
        a = np.asarray(nc.variables[var][sl] if sl is not None else nc.variables[var][:])
    finally:
        nc.close()
    a = a.astype("float64")
    a[a >= FILL] = np.nan
    return a


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    i0 = (BASELINE_YEARS[0] - FILE_START_YEAR) * 12
    i1 = (BASELINE_YEARS[1] - FILE_START_YEAR + 1) * 12
    print(f"baseline window: {BASELINE_YEARS[0]}-{BASELINE_YEARS[1]} "
          f"-> month index [{i0}:{i1}] ({i1 - i0} months)")

    fields, prov = {}, {}
    lat = lon = None

    for var in VARS:
        members = sorted(glob.glob(str(RAW / var / "*" / "*.nc")))
        if not members:
            raise SystemExit(f"no files found for {var} under {RAW / var}")
        print(f"\n{var}: {len(members)} member files")

        ref_path = members[0]
        nc = Dataset(ref_path)
        try:
            lat_v = np.asarray(nc.variables["lat"][:])
            lon_v = np.asarray(nc.variables["lon"][:])
            units = getattr(nc.variables[var], "units", "?")
            long_name = getattr(nc.variables[var], "long_name", "")
            n_time = len(nc.dimensions["time"])
        finally:
            nc.close()
        if units != "kg m-2 s-1":
            raise SystemExit(f"REFUSING: {var} units {units!r} != 'kg m-2 s-1'")
        if n_time != 1032:
            raise SystemExit(f"REFUSING: {var} has {n_time} months, expected 1032")
        if lat is None:
            lat, lon = lat_v, lon_v
        elif not (np.array_equal(lat, lat_v) and np.array_equal(lon, lon_v)):
            raise SystemExit(f"REFUSING: {var} grid differs from {VARS[0]}")

        # --- the invariance contract: every member identical to the reference ---
        ref_first = load_clean(ref_path, var, slice(0, 1))[0]
        differing = []
        for m in members[1:]:
            if not np.allclose(np.nan_to_num(load_clean(m, var, slice(0, 1))[0]),
                               np.nan_to_num(ref_first), equal_nan=True):
                differing.append(Path(m).parent.name)
        print(f"  cross-member identity: {len(members) - len(differing)}/{len(members)} identical"
              + (f"  DIFFERING: {differing}" if differing else ""))

        base = np.nanmean(load_clean(ref_path, var, slice(i0, i1)), axis=0)
        allp = np.nanmean(load_clean(ref_path, var), axis=0)
        static = bool(np.allclose(np.nan_to_num(base), np.nan_to_num(allp), equal_nan=True))
        print(f"  2020s mean == all-period mean: {static}")
        if not static:
            raise SystemExit(
                "REFUSING: pelec is no longer time-invariant — a real temporal signal "
                "exists and must not be collapsed to a static baseline. Re-assess.")
        if differing:
            raise SystemExit(
                "REFUSING: members are no longer identical — a GCM/SSP signal exists "
                "and the single-map assumption is void. Re-assess.")

        fields[var] = base
        prov[var] = {
            "n_member_files": len(members),
            "members_identical": True,
            "reference_file": Path(ref_path).name,
            "reference_sha256": sha256(ref_path),
            "units_native": units,
            "long_name": long_name,
            "time_invariant_verified": True,
        }

    area = cell_area_m2(lat, len(lon))
    for var in VARS:
        vol = np.nan_to_num(fields[var]) * area * SEC_PER_YEAR / 1e12  # km3/yr
        fields[var + "_km3yr"] = vol
        print(f"\n{var}: global total = {vol.sum():,.1f} km3/yr, "
              f"{int((np.nan_to_num(fields[var]) > 0).sum())} non-zero cells")
    ratio = fields["pelecuse_km3yr"].sum() / fields["pelecww_km3yr"].sum()
    print(f"consumption / withdrawal = {100 * ratio:.1f}%")

    if args.dry_run:
        print("\n[dry-run] nothing written.")
        return

    try:
        commit = subprocess.run(["git", "-C", "/Users/arik/github/TCFD", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True).stdout.strip() or "unknown"
    except Exception:
        commit = "unknown"

    ds = Dataset(str(OUT), "w", format="NETCDF4")
    try:
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        vlat = ds.createVariable("lat", "f8", ("lat",)); vlat[:] = lat
        vlat.units, vlat.long_name = "degrees_north", "latitude"
        vlon = ds.createVariable("lon", "f8", ("lon",)); vlon[:] = lon
        vlon.units, vlon.long_name = "degrees_east", "longitude"
        va = ds.createVariable("cell_area", "f8", ("lat", "lon"), zlib=True, complevel=4)
        va[:] = area; va.units, va.long_name = "m2", "grid cell area"

        for var in VARS:
            v = ds.createVariable(var, "f4", ("lat", "lon"), zlib=True, complevel=4,
                                  fill_value=np.nan)
            v[:] = fields[var]
            v.units = "kg m-2 s-1"
            v.long_name = prov[var]["long_name"]
            v.cell_methods = f"time: mean over {BASELINE_YEARS[0]}-{BASELINE_YEARS[1]}"
            vv = ds.createVariable(var + "_km3yr", "f4", ("lat", "lon"), zlib=True,
                                   complevel=4, fill_value=np.nan)
            vv[:] = fields[var + "_km3yr"]
            vv.units = "km3 yr-1"
            vv.long_name = prov[var]["long_name"] + " (volume per grid cell)"

        ds.title = "ISIMIP3b thermoelectric water demand — static 2020s baseline"
        ds.summary = (
            "Potential thermoelectric water withdrawal (pelecww) and consumption "
            "(pelecuse), WaterGAP2-2e, 2015soc, averaged over 2020-2029.")
        ds.WHAT_THIS_SHOWS = (
            "A present-day SPATIAL BASELINE of thermoelectric water demand. "
            "Magnitude validated: contiguous-US withdrawal 186.3 km3/yr vs USGS 2015 "
            "~184 km3/yr (within ~1%). Consumption is ~2.3% of withdrawal, the "
            "once-through-dominated ratio.")
        ds.WHAT_THIS_DOES_NOT_SHOW = (
            "NOT a projection and NOT an uncertainty range. The source field is "
            "constant in time (Jan 2015 == Jan 2100, seasonal max/min ratio 1.000) "
            "and all 15 GCM x SSP files are byte-identical, so there is no temporal "
            "trend, no seasonality, no scenario differentiation and no model spread. "
            "ISIMIP3b publishes no ssp*soc variant for pelec, so no transient "
            "thermoelectric pathway exists. Do not read change over time from this "
            "layer, and do not present it as scenario-dependent.")
        ds.INTENDED_USE = (
            "Baseline context for scenario / mitigation analysis: the existing "
            "thermoelectric demand a new or modified water user competes with. "
            "Forward-looking power or data-centre archetypes must come from our own "
            "withdrawal model, not from this layer.")
        ds.source = "ISIMIP3b OutputData/water_global, WaterGAP2-2e, 2015soc, default"
        ds.institution_of_source = "University of Frankfurt (WaterGAP2-2e)"
        ds.baseline_window = f"{BASELINE_YEARS[0]}-{BASELINE_YEARS[1]}"
        ds.processing_code_commit = commit
        ds.created = "2026-08-26"
        ds.provenance_json = json.dumps(prov)
    finally:
        ds.close()

    print(f"\nwrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")
    prov_path = OUT_DIR / "pelec_baseline_2020s_provenance.json"
    prov_path.write_text(json.dumps(
        {"output": OUT.name, "output_sha256": sha256(str(OUT)),
         "baseline_window": list(BASELINE_YEARS), "variables": prov,
         "global_totals_km3yr": {v: float(fields[v + "_km3yr"].sum()) for v in VARS},
         "code_commit": commit, "created": "2026-08-26"}, indent=2))
    print(f"wrote {prov_path.name}")


if __name__ == "__main__":
    main()
