#!/usr/bin/env python
"""Mask the SSP precip file to the land-surface-model land mask.

Rationale (user decision 2026-08-26, see
waterRiskIndex_beta/v2/enumeration/A3_precip_mask_finding.md):

precip comes from the W5E5/GCM bias-adjusted ATMOSPHERIC forcing and is therefore
valid globally (259,200 cells), while the five land-surface-model variables
(dis, qr, potevap, tws, rootmoist) exist only on land (67,421 cells). The engine's
IDW interpolation zeroes NaN neighbours and renormalises weights, so a globally
valid precip field gives coastal locations a land+ocean stencil for precip while
the water-balance variables use a land-only stencil at the same location. Legacy
(RCP-generation) precip was land-only, so leaving this would also be a silent
departure from every delivery to date.

This script makes precip's mask match the canonical LSM land mask.

Safety properties:
  - the unmasked original is PRESERVED (renamed to *_UNMASKED.nc), never deleted;
  - land values are verified BIT-IDENTICAL before/after — masking only removes
    ocean cells, it never perturbs a retained value;
  - the canonical mask is verified constant across every scenario x value_type x
    decade slice of the source variable before it is used;
  - --dry-run performs every check and writes nothing.

Usage:
    python mask_precip_to_land.py --dry-run
    python mask_precip_to_land.py
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

PROC = Path("/Users/arik/github/TCFD/data/processed/waterindex-ssp")
PRECIP = PROC / "waterIndexUnderlyingData_precip_ssp.nc"
UNMASKED = PROC / "waterIndexUnderlyingData_precip_ssp_UNMASKED.nc"
# Canonical land mask source. dis/qr/potevap/tws share an identical mask
# (67,421 cells); rootmoist is a strict subset (67,097, no extra cells).
MASK_SOURCE = PROC / "waterIndexUnderlyingData_dis_ssp.nc"
MASK_VAR = "dis"
EXPECTED_LAND_CELLS = 67_421


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def canonical_land_mask() -> np.ndarray:
    """Land mask from the LSM reference variable, verified slice-invariant."""
    nc = Dataset(str(MASK_SOURCE))
    try:
        arr = np.asarray(nc.variables[MASK_VAR][:])  # (lat, lon, scen, vt, dec)
    finally:
        nc.close()
    finite = np.isfinite(arr)
    # every scenario x value_type x decade slice must agree on which cells exist
    ref = finite[:, :, 0, 12, 0]
    collapsed = finite.reshape(finite.shape[0], finite.shape[1], -1)
    same = np.all(collapsed == ref[:, :, None], axis=2)
    n_bad = int((~same).sum())
    if n_bad:
        raise SystemExit(
            f"REFUSING: {MASK_VAR} mask is not slice-invariant "
            f"({n_bad} cells differ across scenario/value_type/decade) — "
            "the canonical land mask is not well defined.")
    n = int(ref.sum())
    if n != EXPECTED_LAND_CELLS:
        raise SystemExit(f"REFUSING: land mask has {n} cells, expected {EXPECTED_LAND_CELLS}")
    print(f"  canonical land mask from {MASK_SOURCE.name}: {n} cells, "
          f"slice-invariant across all {collapsed.shape[2]} slices")
    return ref


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if UNMASKED.exists():
        raise SystemExit(f"REFUSING: {UNMASKED.name} already exists — masking appears "
                         "to have been run already. Inspect before re-running.")
    if not PRECIP.exists():
        raise SystemExit(f"missing {PRECIP}")

    print(f"source : {PRECIP}")
    land = canonical_land_mask()

    nc = Dataset(str(PRECIP))
    try:
        pr = np.asarray(nc.variables["pr"][:])
        dims = nc.variables["pr"].dimensions
        gattrs = {a: nc.getncattr(a) for a in nc.ncattrs()}
        coords = {}
        for cv in ("lat", "lon", "scenario", "value_type", "decade"):
            v = nc.variables[cv]
            coords[cv] = (np.asarray(v[:]), v.dtype, v.dimensions,
                          {a: v.getncattr(a) for a in v.ncattrs()})
        var_attrs = {a: nc.variables["pr"].getncattr(a)
                     for a in nc.variables["pr"].ncattrs() if a != "_FillValue"}
    finally:
        nc.close()

    before = int(np.isfinite(pr[:, :, 0, 12, 0]).sum())
    print(f"  precip valid cells before: {before} "
          f"({100 * before / (pr.shape[0] * pr.shape[1]):.1f}% of grid)")
    if before == EXPECTED_LAND_CELLS:
        raise SystemExit("REFUSING: precip already appears land-masked — nothing to do.")

    masked = pr.copy()
    masked[~land, :, :, :] = np.nan
    after = int(np.isfinite(masked[:, :, 0, 12, 0]).sum())
    print(f"  precip valid cells after : {after} "
          f"({100 * after / (pr.shape[0] * pr.shape[1]):.1f}% of grid)")

    # --- verification: retained land values must be bit-identical ---
    a, b = pr[land], masked[land]
    both_nan = np.isnan(a) & np.isnan(b)
    identical = np.array_equal(a[~both_nan], b[~both_nan])
    ocean_clear = bool(np.all(np.isnan(masked[~land])))
    print(f"  land values bit-identical : {identical}")
    print(f"  all ocean cells now NaN   : {ocean_clear}")
    print(f"  land cell count == mask   : {after == EXPECTED_LAND_CELLS}")
    if not (identical and ocean_clear and after == EXPECTED_LAND_CELLS):
        raise SystemExit("REFUSING: post-mask verification failed — nothing written.")

    if args.dry_run:
        print("\n[dry-run] all checks passed; no files written.")
        return

    src_sha = sha256(PRECIP)
    print(f"\n  archiving original -> {UNMASKED.name} (sha256 {src_sha[:16]}…)")
    shutil.move(str(PRECIP), str(UNMASKED))

    try:
        commit = subprocess.run(["git", "-C", "/Users/arik/github/TCFD", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True).stdout.strip() or "unknown"
    except Exception:
        commit = "unknown"

    out = Dataset(str(PRECIP), "w", format="NETCDF4")
    try:
        sizes = dict(zip(dims, masked.shape))
        for d, n in sizes.items():
            out.createDimension(d, n)
        for cv, (vals, dt, vdims, vattrs) in coords.items():
            # _FillValue is immutable after creation — it must go through the
            # createVariable keyword, not setncattr.
            cfill = vattrs.pop("_FillValue", None)
            cvar = (out.createVariable(cv, dt, vdims) if cfill is None
                    else out.createVariable(cv, dt, vdims, fill_value=cfill))
            cvar[:] = vals
            for k, v in vattrs.items():
                cvar.setncattr(k, v)
        pv = out.createVariable("pr", "f4", dims, zlib=True, complevel=4,
                                shuffle=True, fill_value=np.nan)
        pv[:] = masked
        for k, v in var_attrs.items():
            pv.setncattr(k, v)
        for k, v in gattrs.items():
            out.setncattr(k, v)
        out.setncattr("land_mask_applied", "true")
        out.setncattr("land_mask_source", f"{MASK_SOURCE.name} ({MASK_VAR}), "
                                          f"{EXPECTED_LAND_CELLS} cells, slice-invariant")
        out.setncattr("land_mask_rationale",
                      "precip is global atmospheric forcing; masked to the LSM land mask so it "
                      "shares the IDW stencil used by dis/qr/potevap/tws at coastal sites "
                      "(matches legacy RCP-generation behaviour). NOTE: this does NOT make all "
                      "six variables identical -- rootmoist is 324 cells short of this mask, so "
                      "those cells still differ. See v2/enumeration/A3_precip_mask_finding.md")
        out.setncattr("land_mask_applied_date", "2026-08-26")
        out.setncattr("unmasked_source_file", UNMASKED.name)
        out.setncattr("unmasked_source_sha256", src_sha)
        out.setncattr("masking_code_commit", commit)
    finally:
        out.close()

    # --- re-read and re-verify what actually landed on disk ---
    nc = Dataset(str(PRECIP))
    try:
        rt = np.asarray(nc.variables["pr"][:])
    finally:
        nc.close()
    rt_cells = int(np.isfinite(rt[:, :, 0, 12, 0]).sum())
    ra, rb = pr[land], rt[land]
    bn = np.isnan(ra) & np.isnan(rb)
    ok = np.array_equal(ra[~bn], rb[~bn]) and rt_cells == EXPECTED_LAND_CELLS
    print(f"  re-read from disk: {rt_cells} land cells, land values identical to source: {ok}")
    if not ok:
        raise SystemExit("FATAL: on-disk verification failed — restore from "
                         f"{UNMASKED.name}")
    print(f"  wrote {PRECIP.name} ({PRECIP.stat().st_size / 1e6:.1f} MB, "
          f"sha256 {sha256(PRECIP)[:16]}…)")
    print("\nDONE — precip now shares the LSM land mask.")


if __name__ == "__main__":
    main()
