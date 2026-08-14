#!/usr/bin/env python3
"""Regression tests for grid-resolution handling in extraction and the delivery domain.

    python scripts/test_extraction_resolution.py

WHY THIS EXISTS
---------------
On 2026-08-14 the product gained its first layer on a grid other than 0.5 deg: the ISIMIP3b
CaMa-Flood inundation layers at 0.25 deg. Two things had to be true at once, and neither is
self-evident from reading the diff:

  1. THE 0.5 DEG PATH MUST NOT MOVE. Every layer shipped before that date is 0.5 deg and
     their delivered site values are already in customer hands. Extraction geometry used to
     be the literals `search_radius=0.5`, `sigma=0.25`, `cell_size=0.5`; it is now inferred
     from each dataset's own coordinates. Inference that returned 0.4999999 instead of 0.5
     would silently change boundary inclusion and every Gaussian weight. Test 1 pins this.

  2. A MIXED-RESOLUTION REGISTRY MUST NOT COLLAPSE THE DOMAIN. `_domain_mask()` used to
     union every layer's finite mask with `|`, which is an xarray binary op and therefore
     ALIGNS ON COORDINATES. 0.5 deg centres sit at +/-(k+0.25) and 0.25 deg centres at
     +/-(k+0.125), so they share NO coordinate values and the union collapsed to shape
     (0, 0) -- returning OUTSIDE_DOMAIN for every customer site in every delivery, with no
     exception raised. Test 2 pins the fix and test 3 reproduces the original trap so the
     old behaviour cannot quietly return.

These are cheap and deterministic; run them alongside test_shared_baseline.py whenever
extraction, the delivery domain, or a layer's grid changes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.delivery import Domain, _layer_geometry, _point_in_domain, _uniform_or_none  # noqa: E402
from utils.spatial_extract import (  # noqa: E402
    extract_by_point, extract_by_polygon, grid_cell_size, search_radius_for, sigma_for,
)

#: A 0.5 deg layer to test against. Any shipped layer works; this one is small and stable.
HALF_DEG_LAYER = "data/processed/drought-isimip3b_driedarea_annual/driedarea_ssp585_processed.nc"

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def synthetic(cell: float, value: bool = True) -> xr.DataArray:
    lat = np.arange(90 - cell / 2, -90, -cell)
    lon = np.arange(-180 + cell / 2, 180, cell)
    return xr.DataArray(np.full((lat.size, lon.size), value, bool),
                        coords={"lat": lat, "lon": lon}, dims=["lat", "lon"])


def test_half_degree_unchanged() -> None:
    print("\n1. The 0.5 deg path is bit-identical through inferred geometry")
    path = Path(HALF_DEG_LAYER)
    if not path.exists():
        check("0.5 deg reference layer present", False, f"{path} missing -- cannot verify")
        return
    ds = xr.open_dataset(path)
    check("cell size inferred as exactly 0.5", grid_cell_size(ds) == 0.5, str(grid_cell_size(ds)))
    check("search radius resolves to the historical 0.5", search_radius_for(ds) == 0.5)
    check("sigma resolves to the historical 0.25", sigma_for(ds) == 0.25)

    rng = np.random.default_rng(0)
    sites = [(float(rng.uniform(-60, 70)), float(rng.uniform(-180, 180))) for _ in range(400)]
    sites += [(0.0, 180.0), (0.0, -180.0), (89.9, 12.3), (-89.9, -45.0)]  # seam and poles
    compared = diffs = 0
    for lat, lon in sites:
        try:
            old = extract_by_point(ds, lat, lon, variables=["median"],
                                   search_radius=0.5, sigma=0.25)
            new = extract_by_point(ds, lat, lon, variables=["median"])
        except ValueError:
            continue
        for dec, a in old["median"].items():
            b = new["median"][dec]
            compared += 1
            if not ((np.isnan(a) and np.isnan(b)) or a == b):
                diffs += 1
    check("point extraction identical", diffs == 0, f"{compared} values, {diffs} differ")

    from shapely.geometry import box
    poly = box(-123.0, 45.0, -120.0, 47.0)
    old_p = extract_by_polygon(ds, poly, variables=["median"], cell_size=0.5)["median"]
    new_p = extract_by_polygon(ds, poly, variables=["median"])["median"]
    pdiff = sum(0 if ((np.isnan(old_p[d]) and np.isnan(new_p[d])) or old_p[d] == new_p[d])
                else 1 for d in old_p)
    check("polygon extraction identical", pdiff == 0, f"{len(old_p)} values, {pdiff} differ")
    ds.close()


def test_mixed_resolution_domain() -> None:
    print("\n2. A mixed-resolution registry still yields a usable domain")
    dom = Domain(masks=[("half", synthetic(0.5)), ("quarter", synthetic(0.25))],
                 consulted=["half", "quarter"])
    inside = [(45.25, -122.75), (0.125, 0.125), (-33.5, 151.25)]
    check("sites resolve INSIDE the domain",
          all(_point_in_domain(dom, la, lo) for la, lo in inside),
          f"{sum(_point_in_domain(dom, la, lo) for la, lo in inside)}/{len(inside)}")

    empty = Domain(masks=[("half", synthetic(0.5, value=False)),
                          ("quarter", synthetic(0.25, value=False))],
                   consulted=["half", "quarter"])
    check("an all-false domain still reports OUTSIDE",
          not _point_in_domain(empty, 45.25, -122.75))

    geom = _layer_geometry(dom)
    check("per-layer geometry follows each grid",
          geom["half"]["search_radius_degrees"] == 0.5
          and geom["quarter"]["search_radius_degrees"] == 0.25,
          f"half {geom['half']['search_radius_degrees']}, "
          f"quarter {geom['quarter']['search_radius_degrees']}")
    check("mixed delivery reports NO single radius",
          _uniform_or_none(geom, "search_radius_degrees") is None)
    check("uniform delivery still reports one radius",
          _uniform_or_none(_layer_geometry(
              Domain(masks=[("a", synthetic(0.5)), ("b", synthetic(0.5))],
                     consulted=["a", "b"])), "search_radius_degrees") == 0.5)


def test_the_original_trap() -> None:
    print("\n3. The coordinate-alignment trap that caused this is reproduced, not forgotten")
    half, quarter = synthetic(0.5), synthetic(0.25)
    shared = len(set(half.lat.values) & set(quarter.lat.values))
    union = half | quarter
    check("0.5 and 0.25 deg grids share no coordinates", shared == 0, f"{shared} shared")
    check("their xarray union is EMPTY -- never union masks across grids",
          union.size == 0, f"union shape {tuple(union.shape)}")


def main() -> int:
    print("=" * 74)
    print("Extraction / domain resolution regression")
    print("=" * 74)
    test_half_degree_unchanged()
    test_mixed_resolution_domain()
    test_the_original_trap()
    print("\n" + "=" * 74)
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
