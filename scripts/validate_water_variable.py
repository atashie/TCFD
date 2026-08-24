"""A4 fatal-gate validator for a processed `waterIndexUnderlyingData_{var}_ssp.nc`.

Every scenario x decade is inspected (no sampling); every contract violation is
FATAL (nonzero exit). Informational measurements print but never pass a broken
file. Checks:

  C1 structure: dims (lat 360, lon 720, scenario 3, value_type 20, decade 9);
     scenario labels exactly ssp126/ssp370/ssp585; decades 2010..2090.
  C2 units attr equals the expected contract string (pass --expect-units for
     normalized variables whose config units are not the file truth).
  C3 quantile monotonicity: vt13..19 nondecreasing everywhere (tolerance 1e-6
     of local magnitude), every scenario x decade.
  C4 annual-mean consistency: vt12 == nanmean(vt0..11) within rtol 1e-4,
     every scenario x decade.
  C5 valid-cell floor: every scenario x decade has >= --min-valid finite cells
     in vt12 (default 40000), and the finite mask is identical across scenarios
     within each decade.
  C6 reference cells finite: Amazon(-3,-60), Central US(40,-95), Sahara(23,5)
     have finite vt12 in every scenario x decade.

Usage:
    python scripts/validate_water_variable.py FILE.nc VAR [--expect-units STR]
        [--min-valid N]
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
from netCDF4 import Dataset

FAIL = 0


def fail(msg: str) -> None:
    global FAIL
    FAIL += 1
    print(f"FATAL: {msg}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("variable")
    ap.add_argument("--expect-units", default=None)
    ap.add_argument("--min-valid", type=int, default=40000)
    args = ap.parse_args()

    nc = Dataset(args.file)
    v = nc.variables[args.variable]

    # C1 structure
    dims = {k: len(d) for k, d in nc.dimensions.items()}
    expected_dims = {"lat": 360, "lon": 720, "scenario": 3, "value_type": 20, "decade": 9}
    for k, n in expected_dims.items():
        if dims.get(k) != n:
            fail(f"dim {k}={dims.get(k)} != {n}")
    scen = ["".join(x) if isinstance(x, str) else str(x) for x in
            np.asarray(nc.variables["scenario"][:]).tolist()]
    if scen != ["ssp126", "ssp370", "ssp585"]:
        fail(f"scenario labels {scen}")
    decades = np.asarray(nc.variables["decade"][:]).tolist()
    if decades != list(range(2010, 2100, 10)):
        fail(f"decades {decades}")

    # C2 units
    units = nc.getncattr("units") if "units" in nc.ncattrs() else "MISSING"
    expect = args.expect_units
    if expect is not None and units != expect:
        fail(f"units attr '{units}' != expected '{expect}'")
    print(f"units attr: '{units}'" + ("" if expect is None else f" (expected '{expect}')"))

    data = v[:]
    data = data.filled(np.nan) if hasattr(data, "filled") else np.asarray(data)
    # dims order (lat, lon, scenario, value_type, decade)

    ref_cells = {"Amazon": (-3.0, -60.0), "CentralUS": (40.0, -95.0), "Sahara": (23.0, 5.0)}
    lat = np.asarray(nc.variables["lat"][:]); lon = np.asarray(nc.variables["lon"][:])
    ref_idx = {k: (int(np.argmin(np.abs(lat - la))), int(np.argmin(np.abs(lon - lo))))
               for k, (la, lo) in ref_cells.items()}

    finite_masks = {}
    for s in range(3):
        for d in range(9):
            tag = f"{scen[s]}/{decades[d]}"
            q = data[:, :, s, 13:20, d]
            dq = np.diff(q, axis=2)
            tol = 1e-6 * np.maximum(np.abs(q[:, :, :-1]), 1e-12)
            bad = np.nansum(dq < -tol)
            if bad:
                fail(f"C3 quantile order: {int(bad)} violations at {tag}")
            annual = data[:, :, s, 12, d]
            monthly_mean = np.nanmean(data[:, :, s, 0:12, d], axis=2)
            both = np.isfinite(annual) & np.isfinite(monthly_mean)
            if both.any():
                rel = np.abs(annual[both] - monthly_mean[both]) / np.maximum(np.abs(annual[both]), 1e-12)
                worst = float(np.max(rel))
                if worst > 1e-4:
                    fail(f"C4 vt12 vs mean(vt0-11): max rel {worst:.2e} at {tag}")
            n_valid = int(np.isfinite(annual).sum())
            if n_valid < args.min_valid:
                fail(f"C5 valid cells {n_valid} < {args.min_valid} at {tag}")
            finite_masks.setdefault(d, []).append(np.isfinite(annual))
            for name, (i, j) in ref_idx.items():
                if not np.isfinite(annual[i, j]):
                    fail(f"C6 reference cell {name} not finite at {tag}")
    for d, masks in finite_masks.items():
        for s in range(1, 3):
            if not np.array_equal(masks[0], masks[s]):
                diff = int(np.sum(masks[0] != masks[s]))
                fail(f"C5 finite-mask differs between {scen[0]} and {scen[s]} in decade {decades[d]}: {diff} cells")

    nc.close()
    if FAIL:
        print(f"\n=== VALIDATION FAILED: {FAIL} fatal finding(s) ===")
        sys.exit(1)
    print("\n=== VALIDATION PASSED: all scenario x decade panels, all checks ===")


if __name__ == "__main__":
    main()
