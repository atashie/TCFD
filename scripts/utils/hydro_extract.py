"""Point extraction for river-network layers — snap to the nearest river cell IN-BASIN.

WHY THIS EXISTS.  `delivery.py` extracts a point by Gaussian-blending the four
neighbouring cells (`EXTRACT_SIGMA=0.25`, `EXTRACT_SEARCH_RADIUS=0.5`), dropping NaN
neighbours and renormalising the weights.  That is correct for a field that varies
smoothly over ~55 km — temperature, soil moisture, burned area.  It is **wrong for a
routed river variable**, in two ways that both produce a confidently incorrect number:

  1. Neighbouring cells can lie in DIFFERENT CATCHMENTS.  A site 20 km from the Nile but
     on the far side of the divide would receive a blend containing Nile water it has no
     access to.
  2. With a sparse river mask, dropping NaN neighbours and renormalising means the blend
     collapses onto whichever neighbour happens to be a river — silently snapping the
     asset to a river it may not share a basin with.

The rule here replaces the blend entirely.  There is NO interpolation: an asset takes the
value of exactly one cell, chosen as the nearest river cell **within the asset's own
basin**, and the choice is recorded so a reader can audit it.

WHAT THE DELIVERED NUMBER MEANS.  Not "water stress at this site".  It is "water stress on
the river this site would draw from, measured `snap_km` away, on the same drainage".  That
distinction must reach the customer; `snap_km` and `snapped_lat/lon` exist to make it
concrete rather than rhetorical.

Refusal is a feature.  An asset with no river in its basin, or whose nearest in-basin river
is beyond `max_snap_km`, returns a status and no value.  Returning a distant river's ratio
would be worse than returning nothing.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

#: Beyond this, the river is not the asset's water and we refuse rather than reach.
#: One 0.5 deg cell is ~55 km N-S, so this permits roughly one cell of reach plus a
#: margin for a site sitting near a cell corner. It is a DECLARED product parameter,
#: not a physical constant -- raise it only with a reason recorded in the manifest.
DEFAULT_MAX_SNAP_KM = 50.0

STATUS_OK = "OK"
STATUS_OUTSIDE_DOMAIN = "OUTSIDE_ROUTING_DOMAIN"   # outside DDM30 (beyond 83.75N/55.75S)
STATUS_NO_BASIN = "NO_BASIN"                       # ocean, or an unrouted cell
STATUS_NO_RIVER_IN_BASIN = "NO_RIVER_IN_BASIN"     # basin holds no cell above the mask
STATUS_SNAP_TOO_FAR = "SNAP_TOO_FAR"               # nearest in-basin river exceeds the cap

EARTH_R_KM = 6371.0088


@dataclass(frozen=True)
class SnapResult:
    status: str
    row: Optional[int] = None
    col: Optional[int] = None
    snapped_lat: Optional[float] = None
    snapped_lon: Optional[float] = None
    snap_km: Optional[float] = None
    basin_id: Optional[int] = None
    upstream_cells: Optional[int] = None

    @property
    def ok(self) -> bool:
        return self.status == STATUS_OK


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance. Index distance will not do: 0.5 deg of longitude is ~55 km
    at the equator and ~19 km at 70N, so a nearest-neighbour search in index space picks
    the wrong cell at high latitude."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    dl = (dl + np.pi) % (2 * np.pi) - np.pi          # wrap the dateline
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * EARTH_R_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


class RiverSnapper:
    """Holds the routing grid and a fixed river mask; snaps points to in-basin rivers.

    `river_mask` is the SAME fixed mask the layer uses (derived once from the shared
    2020s baseline climatology). Deriving a different mask here would let an asset snap
    to a cell the layer itself masked out.
    """

    def __init__(self, basins: np.ndarray, lat: np.ndarray, lon: np.ndarray,
                 river_mask: np.ndarray, upstream_cells: Optional[np.ndarray] = None,
                 max_snap_km: float = DEFAULT_MAX_SNAP_KM):
        if basins.shape != river_mask.shape:
            raise ValueError(f"basins {basins.shape} != river_mask {river_mask.shape}")
        self.basins = basins
        self.lat = np.asarray(lat, float)
        self.lon = np.asarray(lon, float)
        self.river_mask = np.asarray(river_mask, bool)
        self.upstream_cells = upstream_cells
        self.max_snap_km = float(max_snap_km)

        # Precompute per-basin river-cell index lists once; a delivery snaps many assets.
        self._by_basin: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        rr, cc = np.nonzero(self.river_mask)
        if rr.size:
            bid = self.basins[rr, cc]
            good = np.isfinite(bid)
            rr, cc, bid = rr[good], cc[good], bid[good].astype(np.int64)
            order = np.argsort(bid, kind="stable")
            rr, cc, bid = rr[order], cc[order], bid[order]
            uniq, start = np.unique(bid, return_index=True)
            stop = np.append(start[1:], bid.size)
            for b, s, e in zip(uniq, start, stop):
                self._by_basin[int(b)] = (rr[s:e], cc[s:e])

    def _locate(self, lat: float, lon: float) -> Optional[tuple[int, int]]:
        lon = ((float(lon) + 180.0) % 360.0) - 180.0
        if not (self.lat.min() - 0.25 <= lat <= self.lat.max() + 0.25):
            return None
        r = int(np.argmin(np.abs(self.lat - lat)))
        c = int(np.argmin(np.abs(((self.lon - lon + 180) % 360) - 180)))
        return r, c

    def snap(self, lat: float, lon: float) -> SnapResult:
        loc = self._locate(lat, lon)
        if loc is None:
            return SnapResult(STATUS_OUTSIDE_DOMAIN)
        r, c = loc
        bid = self.basins[r, c]
        if not np.isfinite(bid):
            return SnapResult(STATUS_NO_BASIN)
        bid = int(bid)
        cand = self._by_basin.get(bid)
        if cand is None:
            return SnapResult(STATUS_NO_RIVER_IN_BASIN, basin_id=bid)
        rr, cc = cand
        d = haversine_km(lat, lon, self.lat[rr], self.lon[cc])
        k = int(np.argmin(d))
        if d[k] > self.max_snap_km:
            return SnapResult(STATUS_SNAP_TOO_FAR, basin_id=bid, snap_km=float(d[k]))
        up = None
        if self.upstream_cells is not None:
            up = int(self.upstream_cells[rr[k], cc[k]])
        return SnapResult(
            STATUS_OK, row=int(rr[k]), col=int(cc[k]),
            snapped_lat=float(self.lat[rr[k]]), snapped_lon=float(self.lon[cc[k]]),
            snap_km=float(d[k]), basin_id=bid, upstream_cells=up,
        )
