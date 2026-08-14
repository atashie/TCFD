#!/usr/bin/env python3
"""Backfill `spatial_resolution_degrees` onto every processed layer that lacks it.

    python scripts/backfill_grid_attribute.py --dry-run
    python scripts/backfill_grid_attribute.py --run

WHY THIS EXISTS
---------------
Until 2026-08-14 every shipped layer was 0.5 deg, so resolution was an assumption rather
than a recorded fact and nothing validated it -- the first 0.25 deg layer passed
`test_shared_baseline.py` cleanly while contradicting OUTPUT-SPEC.md. That test now checks
the grid, which means existing files need the attribute before the check can be required.

The value is MEASURED FROM THE FILE'S OWN COORDINATES, never assumed to be 0.5: the point of
the exercise is to stop asserting the grid. A file whose spacing is irregular or non-square
is reported and skipped rather than stamped with a number that would then look verified.

THE DATA IS NOT TOUCHED. Every contract variable is digested before and after and asserted
identical. Note that the FILE checksum does change -- writing a global attribute rewrites
the header -- so only the variable arrays can be, and are, asserted equal.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import netCDF4

PROCESSED = Path("data/processed")
ATTR = "spatial_resolution_degrees"


def array_digest(path: Path) -> str:
    h = hashlib.sha256()
    with netCDF4.Dataset(path) as ds:
        for name in sorted(ds.variables):
            h.update(name.encode())
            h.update(np.asarray(ds[name][:]).tobytes())
    return h.hexdigest()


def measure(path: Path):
    """(cell_size, problem). cell_size is None when the grid is not regular and square."""
    with netCDF4.Dataset(path) as ds:
        if "lat" not in ds.variables or "lon" not in ds.variables:
            return None, "no lat/lon coordinates"
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
    if lat.size < 2 or lon.size < 2:
        return None, "degenerate coordinates"
    dlat, dlon = np.abs(np.diff(lat)), np.abs(np.diff(lon))
    if not (np.allclose(dlat, dlat[0], rtol=0, atol=1e-6)
            and np.allclose(dlon, dlon[0], rtol=0, atol=1e-6)):
        return None, "irregular spacing"
    if abs(float(dlat[0]) - float(dlon[0])) > 1e-6:
        return None, f"non-square cells ({dlat[0]:.6f} vs {dlon[0]:.6f})"
    return round(float(dlat[0]), 6), ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true", help="write (default is a dry run)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    files = sorted(PROCESSED.glob("*/*_processed.nc"))
    if not files:
        print(f"no processed layers under {PROCESSED}", file=sys.stderr)
        return 1

    todo, already, bad = [], 0, []
    for f in files:
        with netCDF4.Dataset(f) as ds:
            present = ATTR in ds.ncattrs()
        cell, problem = measure(f)
        if problem:
            bad.append((f, problem))
            continue
        if present:
            already += 1
            continue
        todo.append((f, cell))

    by_res: dict[float, int] = {}
    for _f, cell in todo:
        by_res[cell] = by_res.get(cell, 0) + 1
    print(f"{len(files)} processed files: {len(todo)} to stamp, {already} already carry "
          f"{ATTR}, {len(bad)} unusable")
    for cell, n in sorted(by_res.items()):
        print(f"  {cell:>8.4f} deg : {n} file(s)")
    for f, problem in bad:
        print(f"  SKIP {f.parent.name}/{f.name}: {problem}")

    if not args.run:
        print("\nDRY RUN -- re-run with --run to write.")
        return 0

    for f, cell in todo:
        before = array_digest(f)
        with netCDF4.Dataset(f, "a") as ds:
            ds.setncattr(ATTR, float(cell))
        after = array_digest(f)
        if before != after:
            print(f"  ABORT: {f} variable data changed during a header write", file=sys.stderr)
            return 2
        print(f"  stamped {f.parent.name}/{f.name} = {cell} deg (arrays verified identical)")
    print(f"\nstamped {len(todo)} file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
