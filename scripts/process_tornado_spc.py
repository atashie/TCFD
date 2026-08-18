#!/usr/bin/env python3
"""Process the NOAA SPC tornado database into a CONUS historical tornado hazard layer.

PRODUCT: TCFD/CDP. This layer does NOT satisfy OUTPUT-SPEC.md's decadal contract and
must not be read as if it did -- see "OUTPUT CONTRACT" below. Options review and the
ISIMIP-absence receipt: docs/tornado-data-options-2026-08-18.md.

===========================================================================
WHAT IS PUBLISHED, AND WHY IT IS A RATE RATHER THAN A COUNT
===========================================================================
The obvious design -- annual tornado count per cell, pooled over years, reported as
median with an inter-quartile range -- was measured before it was written, and it is
DEGENERATE. Within-cell quantiles of the annual count, over occupied cells only
(cells that have EVER recorded a tornado, the most favourable denominator available):

    resolution  block      q25   median   q75
    0.25 deg    annual       0        0     0
    0.25 deg    5-year       0        0     1
    0.25 deg    decadal      0        0     1
    0.50 deg    annual       0        0     0
    0.50 deg    decadal      0        2     5

91.6% of (cell, year) pairs are zero at 0.25 deg; 97.3% at 0.1 deg. So median,
lower_ci and upper_ci would have published as identically (0, 0, 0) over the entire
map -- three of four variables carrying no information, in a file that would still
pass a structural contract check. This is the documented pooled_mean_zero_inflated
failure (CLAUDE.md: the median branch "erased 93% of exposed land on one layer") in
the most extreme form yet measured in this repository.

The estimand for a point process is therefore a RATE, not a pooled count:

    median     = tornado crossings per 10^4 km^2 per year, over the whole window
    lower_ci   = 25th percentile of the rate's posterior
    upper_ci   = 75th percentile of the rate's posterior

All three are percentiles of the SAME Gamma(k + 1/2, T) posterior (Jeffreys prior on
a Poisson rate). That is what the user's "quantiles 0.25 and 0.75" has to mean once
the annual series is known to be degenerate: quartiles OF THE RATE ESTIMATE, not of
the count series. It answers the question a customer actually has -- is 3 events in
76 years distinguishable from 6 -- which an empirical IQR cannot, being 0 either way.

CORRECTED 2026-08-18 after external review. The first build published `median` as the
MLE k/T next to Gamma posterior quartiles, then CLIPPED the lower quartile down to it
so the triple would order. That clip covered a genuine incoherence -- an MLE on the
boundary at k=0 cannot be bracketed by quantiles of a continuous posterior -- and it
made the published interval something other than the Gamma quartile interval it was
documented as. Both are now percentiles of one posterior, ordering holds without
clipping, and `--estimator mle` reproduces the old behaviour for comparison.

The zero tier of the percentile is therefore keyed on the OBSERVED COUNT rather than
on the rate being exactly zero -- see two_tier_percentile.

===========================================================================
COUNTING: TRACKS, NOT TOUCHDOWNS
===========================================================================
The hazard question is "does a tornado pass over this location", not "does one
first touch down here". Each tornado therefore increments EVERY cell its damage
path crosses by 1 (not by a fraction -- each such cell did experience a tornado).
Consequence: the sum of n_events over cells EXCEEDS the number of tornadoes, by
design, and n_events is a per-cell crossing count rather than a partition of the
record.

Measured on this database: median path length 1.0 mi, so most tornadoes are
genuinely sub-cell and the choice barely moves the map at 0.25 deg -- 4.2% of
tracks exceed one cell (vs 14.3% at 0.1 deg, max 234.7 mi). But long-track
tornadoes are disproportionately the violent ones, and touchdown-only assignment
displaces their risk systematically upstream. 35.5% of records carry no usable end
point and fall back to touchdown; `--geometry touchdown` forces that everywhere for
comparison.

===========================================================================
RESOLUTION: 0.25 DEG, AND WHY NOT FINER
===========================================================================
0.25 deg nests exactly inside the 0.5 deg ISIMIP grid every other layer uses -- but
the GRID nests, the VALUES DO NOT AGGREGATE. Measured 2026-08-18: averaging the four
0.25 deg children area-weighted onto their 0.5 deg parent reads 7.9% high at the
median cell (p90 +23%), and total path-crossings inflate 1.123x from 0.5 to 0.25 deg,
because one track crossing several children is counted once at the coarser
resolution and several times at the finer. Do not resample this layer between
resolutions; recompute it. 0.1 deg was rejected on three
independent measurements, any one of which is disqualifying:

  * It is below the source's own positional accuracy. 25.6% of 1950-1975 records
    sit on an exact 0.1 deg latitude multiple, with 1,373 distinct latitudes for
    16,830 events -- a county-centroid geocode, not a survey position. A 0.1 deg
    grid RESOLVES that quantization into the map as stripes on round coordinates.
    (2000-2025: 6.0% and 20,139 distinct values -- the artefact is era-dependent,
    so it would also fake a spatial trend.)
  * 46.1% of occupied cells would hold exactly one tornado in 76 years (79.0% for
    F2+). That maps the record, not the hazard.
  * 14.3% of tracks exceed one cell, so track geometry stops being optional while
    a third of records have none.

THIS LAYER SCREENS; IT DOES NOT SITE. A damage path is ~100 m wide against a
~28 km cell, so `resolution_caveat` is set and both reports must disclose it. The
cell value is a track-INTERSECTION frequency divided by cell area, which is not the
same thing as an areal tornado intensity (see above). The probability that a given
building is struck is far smaller and is NOT the footprint-to-cell-area ratio: it
depends on swath width, path length and orientation, none of which this layer
models. Do not convert a cell value into a strike probability.

===========================================================================
NO TREND, AND THE REASON MATTERS
===========================================================================
No ols_slope and no sen_slope are emitted. NOT because 76 years cannot support a
trend -- it comfortably can -- but because the trend would measure the OBSERVING
SYSTEM. Measured on this file, per year by decade:

    period        all/yr    F0/yr    F2+/yr   unrated/yr
    1950-1959      479.3    105.3     188.8          0.0
    1990-1999     1213.7    737.1     149.4          0.0
    2020-2025     1349.0    481.0     153.5        204.3

Reports rose 2.8x; weak reports rose 7.7x; strong (F2+) reports FELL 19% and have
no modern trend. Radar, spotter networks, population and cameras explain the
growth; climate does not. Recording "not computable" instead of this would invite
someone to notice it is computable and add a slope.

===========================================================================
THE MAGNITUDE THRESHOLD IS A LADDER, NOT A PICK
===========================================================================
Restricting to F2+ is the standard mitigation for the reporting bias above, and it
costs 82% of the records; at 0.25 deg it leaves 38.5% of occupied cells resting on
a single event. The fix for the confound and the cure for the thinness pull
opposite ways, so every rung is built in one pass and the comparison is the
deliverable: `all`, `f1plus`, `f2plus`, `f3plus`.

Unrated reports (mag = -9) are counted in `all` and EXCLUDED from every F-rung,
because they cannot be verified against a threshold. They now run ~204/yr (~15% of
the modern record), so the F-rungs understate recent years -- `n_unrated` is
written to the file so the size of that effect is auditable rather than asserted.

===========================================================================
COVERAGE MASK
===========================================================================
SPC observes the United States only. Cells in Mexico, Canada and the ocean are
zero because they are UNOBSERVED, and scoring them would rank reporting systems
rather than hazard -- the same defect as putting Bangladesh at the 1st percentile
on a global grid. The percentile is therefore ranked over CONUS land cells only,
and everything outside the mask is NaN, never 0.

===========================================================================
OUTPUT CONTRACT
===========================================================================
This is `observational-historical-v1`, NOT the OUTPUT-SPEC.md decadal contract.
Deliberately absent: the `decade` dimension, `ols_slope`, `sen_slope`, `n_members`,
`n_models`. n_members/n_models are NOT faked to 1 -- there is no ensemble here, and
a 1 would read as a thin ensemble rather than as a different kind of product.
scripts/test_shared_baseline.py will reject this file; that is correct behaviour
and needs a verifier decision before delivery, not a silent workaround.

USAGE
    python3 scripts/process_tornado_spc.py                       # full record 1950-2025
    python3 scripts/process_tornado_spc.py --start-year 1996     # modern observing era
    python3 scripts/process_tornado_spc.py --geometry touchdown  # sensitivity
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import gamma as gamma_dist
from scipy.stats import rankdata

RAW_DIR = Path("data/raw/tornado-spc")
SPC_CSV = RAW_DIR / "spc_tornadoes_1950_2025.csv"
NE_GEOJSON = RAW_DIR / "ne_50m_admin_0_countries.geojson"
OUT_ROOT = Path("data/processed")

# THE GRID IS GLOBAL; THE DATA IS CONUS. These are different things and conflating
# them broke a delivery.
#
# The first build published a regional 24-50N / 66-125W grid. That is defensible as a
# file and wrong as a product: `spatial_extract.extract_by_point` RAISES when a site
# falls outside the grid extent, so a Honolulu or Frankfurt warehouse crashed the whole
# delivery instead of being reported as unmodelled. That guard is correct and must stay
# strict -- an off-grid point on a global layer means bad coordinates -- so the fix is to
# publish on the global grid every other layer uses and let the MASK carry the meaning.
#
# Result: any site on Earth resolves. Inside CONUS it gets a value; outside it gets NaN,
# which the delivery reports as OFF_LAYER_MASK / OUTSIDE_DOMAIN -- "this layer does not
# model your site" -- and never as a low score. Cost is ~0 on disk: the array is
# overwhelmingly NaN and compresses away.
GRID_LAT_MIN, GRID_LAT_MAX = -90.0, 90.0
GRID_LON_MIN, GRID_LON_MAX = -180.0, 180.0

# Where SPC actually observes. Used to filter records and to bound the (expensive)
# point-in-polygon mask test; AK and HI fall outside it by construction, which is
# intended -- they are separate tornado regimes with far thinner records.
LAT_MIN, LAT_MAX = 24.0, 50.0
LON_MIN, LON_MAX = -125.0, -66.0

EARTH_RADIUS_KM = 6371.0

RUNGS = {
    "all": (None, "All reported tornadoes, any rating, including unrated"),
    "f1plus": (1, "(E)F1 and stronger"),
    "f2plus": (2, "(E)F2 and stronger -- the reporting-stable subset"),
    "f3plus": (3, "(E)F3 and stronger -- violent tornadoes"),
}


# --------------------------------------------------------------------------
# grid
# --------------------------------------------------------------------------
def build_grid(res: float):
    lat_edges = np.arange(GRID_LAT_MIN, GRID_LAT_MAX + res / 2, res)
    lon_edges = np.arange(GRID_LON_MIN, GRID_LON_MAX + res / 2, res)
    lat_centres = (lat_edges[:-1] + lat_edges[1:]) / 2
    lon_centres = (lon_edges[:-1] + lon_edges[1:]) / 2
    return lat_edges, lon_edges, lat_centres, lon_centres


def cell_area_km2(lat_edges: np.ndarray, lon_edges: np.ndarray) -> np.ndarray:
    """Spherical cell area. Area varies by ~40% across this box, so the rate must be
    area-normalised or Texas outranks Iowa purely by being further south."""
    dlon = np.deg2rad(lon_edges[1] - lon_edges[0])
    sin_edges = np.sin(np.deg2rad(lat_edges))
    band = EARTH_RADIUS_KM**2 * dlon * (sin_edges[1:] - sin_edges[:-1])  # per lat band
    return np.repeat(band[:, None], len(lon_edges) - 1, axis=1)


# --------------------------------------------------------------------------
# CONUS mask
# --------------------------------------------------------------------------
def build_conus_mask(lat_centres, lon_centres) -> np.ndarray:
    from shapely.geometry import shape, Point
    from shapely.prepared import prep
    from shapely.ops import unary_union

    with NE_GEOJSON.open() as fh:
        gj = json.load(fh)

    geoms = []
    for feat in gj["features"]:
        props = feat.get("properties", {})
        names = {str(props.get(k, "")).lower() for k in
                 ("NAME", "ADMIN", "NAME_LONG", "SOVEREIGNT", "GEOUNIT")}
        if "united states of america" in names or "united states" in names:
            geoms.append(shape(feat["geometry"]))
    if not geoms:
        raise SystemExit("FAILED: no USA feature found in Natural Earth GeoJSON -- "
                         "property names may have changed; inspect before assuming absence.")

    usa = prep(unary_union(geoms))
    mask = np.zeros((len(lat_centres), len(lon_centres)), dtype=bool)
    # Point-in-polygon over a global 0.25 deg grid is ~1M tests; restrict to the CONUS
    # bounding box first, which is where the only True cells can be anyway.
    lat_idx = np.where((lat_centres > LAT_MIN) & (lat_centres < LAT_MAX))[0]
    lon_idx = np.where((lon_centres > LON_MIN) & (lon_centres < LON_MAX))[0]
    for i in lat_idx:
        la = float(lat_centres[i])
        for j in lon_idx:
            if usa.contains(Point(float(lon_centres[j]), la)):
                mask[i, j] = True
    return mask


# --------------------------------------------------------------------------
# counting
# --------------------------------------------------------------------------
def track_cells(row, res: float, nlat: int, nlon: int, geometry: str):
    """Return the set of (i, j) cells this tornado's damage path crosses.

    Sampled densely along the segment (res/16). This is DENSE SAMPLING, NOT EXACT
    RASTER TRAVERSAL: a line can clip an arbitrarily small corner of a cell between
    two samples, so a corner-grazing cell can still be missed. An exact supercover
    (Amanatides-Woo) traversal would close that gap. The residual is confined to
    corner clips of near-diagonal tracks and is small against a median path length
    of 1.0 mi, but the earlier claim that "no crossed cell is skipped" was false and
    is withdrawn. An unmoving or endpoint-less record contributes only its touchdown
    cell.
    """
    slat, slon = row[0], row[1]
    elat, elon = row[2], row[3]
    # SPC encodes a missing endpoint as literal 0, but nothing guarantees a future
    # release will not carry NaN -- and a NaN span makes int(np.ceil(...)) raise.
    # Verified 2026-08-18: zero non-finite coordinates in the current file.
    if not (np.isfinite(slat) and np.isfinite(slon)):
        return set()
    if not (np.isfinite(elat) and np.isfinite(elon)):
        elat = elon = 0.0

    def to_ij(la, lo):
        i = int((la - GRID_LAT_MIN) // res)
        j = int((lo - GRID_LON_MIN) // res)
        if 0 <= i < nlat and 0 <= j < nlon:
            return (i, j)
        return None

    start = to_ij(slat, slon)
    if geometry == "touchdown" or elat == 0 or elon == 0:
        return {start} if start else set()

    span = max(abs(elat - slat), abs(elon - slon))
    if span == 0:
        return {start} if start else set()

    n = int(np.ceil(span / (res / 16))) + 1
    lats = np.linspace(slat, elat, n)
    lons = np.linspace(slon, elon, n)
    out = set()
    for la, lo in zip(lats, lons):
        ij = to_ij(la, lo)
        if ij:
            out.add(ij)
    return out


def count_grid(df: pd.DataFrame, res: float, nlat: int, nlon: int, geometry: str):
    counts = np.zeros((nlat, nlon), dtype=np.int64)
    arr = df[["slat", "slon", "elat", "elon"]].to_numpy()
    fallback = 0
    for row in arr:
        if row[2] == 0 or row[3] == 0:
            fallback += 1
        for (i, j) in track_cells(row, res, nlat, nlon, geometry):
            counts[i, j] += 1
    return counts, fallback


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def rate_and_interval(counts: np.ndarray, area_km2: np.ndarray, n_years: int,
                      estimator: str = "posterior"):
    """Central rate estimate and its 25th/75th percentiles.

    Units: crossings per 10^4 km^2 per year.

    TWO ESTIMATORS, AND WHY THE DEFAULT CHANGED (2026-08-18, external review).
    The first version published `median` = the MLE k/T alongside Gamma(k+1/2, T)
    posterior quartiles, and then CLIPPED the lower quartile down to the MLE. That
    clip was a fudge covering a real incoherence: k/T is a frequentist point estimate
    that sits on the boundary at k=0, while every proper continuous posterior quantile
    is strictly positive, so the triple could not be ordered without forcing it. The
    resulting interval was then NOT the Gamma quartile interval it was documented as.

      "posterior" (default) -- median, lower_ci, upper_ci are the 50th, 25th and 75th
          percentiles of the SAME Gamma(k+1/2, T) posterior. Coherent by construction,
          ordered without clipping, and honest at k=0: a cell with no report in 76
          years does not have zero tornado risk, it has a rate small enough that no
          report was likely. Consequence: `median` is small-but-positive where nothing
          was recorded, so THE ZERO TIER IS KEYED ON THE OBSERVED COUNT, not on the
          rate being exactly 0 (see two_tier_percentile).
      "mle" -- the original k/T with clipped bounds, retained only so the change is
          reproducible and comparable. It is documented in the output as a hybrid,
          because that is what it is.

    THE JEFFREYS PRIOR IS DEFENSIBLE, THE AREA NORMALISATION IS THE WEAK LINK.
    p(lambda) ~ lambda^-1/2 is standard for a Poisson rate, and the Poisson form
    survives track counting -- thinning a Poisson process of tornado objects to those
    intersecting one cell is still Poisson. Measured 2026-08-18, the variance-to-mean
    ratio of annual per-cell counts is 0.99 at the median cell (p90 1.70; 16.4% of
    cells above 1.5 for all reports, 6.2% for F2+), so the marginal Poisson assumption
    holds well typically and understates the interval in the clustered tail.
    What does NOT hold is reading the result as an area-normalised intensity -- see
    `rate_is_resolution_dependent` in the output attributes.
    """
    exposure = (area_km2 / 1e4) * n_years  # T, in 10^4 km^2 * years
    alpha = counts + 0.5

    if estimator == "posterior":
        med = gamma_dist.ppf(0.50, a=alpha, scale=1.0 / exposure)
        lo = gamma_dist.ppf(0.25, a=alpha, scale=1.0 / exposure)
        hi = gamma_dist.ppf(0.75, a=alpha, scale=1.0 / exposure)
    elif estimator == "mle":
        med = counts / exposure
        lo = np.minimum(gamma_dist.ppf(0.25, a=alpha, scale=1.0 / exposure), med)
        hi = gamma_dist.ppf(0.75, a=alpha, scale=1.0 / exposure)
    else:
        raise ValueError(f"unknown estimator {estimator!r}")

    assert np.all(lo <= med) and np.all(med <= hi), "CI ordering invariant violated"
    return med, lo, hi


def two_tier_percentile(rate: np.ndarray, counts: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Cells with NO observed crossing -> 1. The rest ranked 2..100 among themselves.

    THE TIER IS KEYED ON THE OBSERVED COUNT, NOT ON THE RATE BEING ZERO. Under the
    posterior estimator the rate is small-but-positive everywhere, so `rate == 0`
    would select nothing; and the tier is meant to separate "nothing was ever
    recorded here" from "something was", which is a statement about the observation,
    not about the estimator. Keying on counts makes the two estimators produce the
    same tiering, which is what lets them be compared at all.

    NOT A PERCENTILE OF ALL CELLS, and the name overstates it. It is a rank score
    over the OCCUPIED population: percentile 50 is the median among cells that have
    ever recorded a tornado, not median risk, and the zero tier is 31-86% of CONUS
    depending on rung. Ties take average ranks, so with heavy ties the lowest
    occupied cells need not land exactly on 2 (measured: 2.005 on `all`, 2.000 on
    the F-rungs). That inversion is why this is a must-disclose caveat.
    """
    pct = np.full(rate.shape, np.nan)
    inmask = mask & np.isfinite(rate)
    zero = inmask & (counts <= 0)
    nonzero = inmask & (counts > 0)

    pct[zero] = 1.0
    vals = rate[nonzero]
    if vals.size == 0:
        return pct
    if vals.size == 1:
        pct[nonzero] = 100.0
        return pct
    ranks = rankdata(vals, method="average")
    pct[nonzero] = 2.0 + 98.0 * (ranks - 1) / (vals.size - 1)
    return pct


# --------------------------------------------------------------------------
# reference sites -- check the layer where the thing exists, before and after
# --------------------------------------------------------------------------
REFERENCE_SITES = [
    ("Moore, OK", 35.34, -97.49, "high"),
    ("Joplin, MO", 37.08, -94.51, "high"),
    ("Tuscaloosa, AL", 33.21, -87.57, "high"),
    ("Wichita, KS", 37.69, -97.34, "high"),
    ("Seattle, WA", 47.61, -122.33, "low"),
    ("Phoenix, AZ", 33.45, -112.07, "low"),
    ("Caribou, ME", 46.86, -68.01, "low"),
]


def report_reference_sites(ds: xr.Dataset, label: str):
    print(f"\n  Reference sites [{label}] -- expect high in the Plains/Dixie, low in the West/NE:")
    print(f"    {'site':18s} {'expect':6s} {'rate':>9s} {'q25':>8s} {'q75':>8s} {'pctile':>7s} {'n':>5s}")
    for name, la, lo, expect in REFERENCE_SITES:
        try:
            pt = ds.sel(lat=la, lon=lo, method="nearest")
        except Exception:  # noqa: BLE001
            print(f"    {name:18s} {expect:6s}  (outside grid)")
            continue
        med = float(pt["median"].values)
        if not np.isfinite(med):
            print(f"    {name:18s} {expect:6s}  (outside CONUS mask)")
            continue
        print(f"    {name:18s} {expect:6s} {med:9.3f} {float(pt['lower_ci']):8.3f} "
              f"{float(pt['upper_ci']):8.3f} {float(pt['percentile']):7.1f} "
              f"{int(pt['n_events']):5d}")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolution", type=float, default=0.25)
    ap.add_argument("--start-year", type=int, default=1950)
    ap.add_argument("--end-year", type=int, default=2025)
    ap.add_argument("--geometry", choices=["track", "touchdown"], default="track")
    ap.add_argument("--estimator", choices=["posterior", "mle"], default="posterior",
                    help="posterior (default): coherent Gamma median+quartiles. "
                         "mle: k/T with a clipped lower bound (the superseded hybrid).")
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    ap.add_argument("--label", default=None,
                    help="directory suffix; default derived from the window")
    args = ap.parse_args()

    if not SPC_CSV.exists():
        print(f"FAILED: {SPC_CSV} not found. Run scripts/download_tornado_spc.py first.")
        return 1

    res = args.resolution
    label = args.label or ("full" if args.start_year <= 1950 else f"from{args.start_year}")

    print(f"Tornado hazard layer -- SPC occurrence, CONUS, {args.start_year}-{args.end_year}")
    print(f"  resolution {res} deg, geometry={args.geometry}, label={label}")

    df = pd.read_csv(SPC_CSV, low_memory=False)
    n_all = len(df)
    df = df[(df.yr >= args.start_year) & (df.yr <= args.end_year)]
    df = df[(df.slat > LAT_MIN) & (df.slat < LAT_MAX)
            & (df.slon > LON_MIN) & (df.slon < LON_MAX)]
    n_years = args.end_year - args.start_year + 1
    print(f"  {n_all:,} records in file -> {len(df):,} in window and box ({n_years} years)")

    lat_edges, lon_edges, lat_c, lon_c = build_grid(res)
    nlat, nlon = len(lat_c), len(lon_c)
    area = cell_area_km2(lat_edges, lon_edges)
    print(f"  grid {nlat} x {nlon} (global {res} deg)")

    print("  building CONUS mask from Natural Earth ...")
    mask = build_conus_mask(lat_c, lon_c)
    a_in = area[mask]
    print(f"  CONUS land cells: {mask.sum():,} of {mask.size:,} global cells "
          f"({100*mask.sum()/mask.size:.1f}%); everything else is NaN, never 0")
    print(f"  in-mask cell area {a_in.min():,.0f}-{a_in.max():,.0f} km^2 "
          f"({100*(a_in.max()/a_in.min()-1):.0f}% variation -> area normalisation is required)")

    n_unrated = int((df.mag == -9).sum())
    written = []

    for rung, (threshold, description) in RUNGS.items():
        sub = df if threshold is None else df[df.mag >= threshold]
        counts, fallback = count_grid(sub, res, nlat, nlon, args.geometry)

        rate, lo, hi = rate_and_interval(counts, area, n_years, args.estimator)
        for arr in (rate, lo, hi):
            arr[~mask] = np.nan
        counts_f = counts.astype(np.float64)
        counts_f[~mask] = np.nan
        pct = two_tier_percentile(rate, counts, mask)

        occupied = int(np.nansum(counts_f > 0))
        in_mask = int(mask.sum())

        ds = xr.Dataset(
            {
                "median": (("lat", "lon"), rate.astype(np.float32)),
                "lower_ci": (("lat", "lon"), lo.astype(np.float32)),
                "upper_ci": (("lat", "lon"), hi.astype(np.float32)),
                "percentile": (("lat", "lon"), pct.astype(np.float32)),
                "n_events": (("lat", "lon"), counts_f.astype(np.float32)),
                "conus_mask": (("lat", "lon"), mask.astype(np.int8)),
            },
            coords={"lat": lat_c.astype(np.float64), "lon": lon_c.astype(np.float64)},
        )

        est_note = ("All three are percentiles of ONE Gamma(k+1/2, T) posterior."
                    if args.estimator == "posterior" else
                    "HYBRID (superseded): median is the MLE k/T and lower_ci is the Gamma "
                    "quartile CLIPPED to it, so the pair is not a posterior interval.")
        ds["median"].attrs = {
            "long_name": f"Tornado crossing rate ({description})",
            "units": "crossings per 1e4 km2 per year",
            "note": ("Posterior median of the rate." if args.estimator == "posterior"
                     else "MLE k/T; exactly 0 where nothing was reported.") + " " + est_note,
        }
        ds["lower_ci"].attrs = {
            "long_name": "25th percentile of the Gamma(k+1/2, T) rate posterior",
            "units": "crossings per 1e4 km2 per year",
            "note": ("NOT an empirical quantile of the annual count series -- that series "
                     "is 91.6% zeros at this resolution and its IQR is identically (0,0). "
                     + est_note),
        }
        ds["upper_ci"].attrs = {
            "long_name": "75th percentile of the Gamma(k+1/2, T) rate posterior",
            "units": "crossings per 1e4 km2 per year",
            "note": ("Where nothing was recorded this is the informative number: the rate "
                     "bound implied by seeing no report in the window. " + est_note),
        }
        ds["percentile"].attrs = {
            "long_name": "Two-tier percentile of tornado crossing rate",
            "units": "1",
            "note": ("Zeros -> 1; non-zeros ranked 2-100 against CONUS cells with at least "
                     "one crossing. Percentile 50 means the median AMONG CELLS THAT HAVE "
                     "EVER RECORDED ONE, not median risk."),
        }
        ds["n_events"].attrs = {
            "long_name": "Tornado damage paths crossing this cell in the window",
            "units": "count",
            "note": ("Sum over cells EXCEEDS the tornado count by design -- a multi-cell "
                     "track increments every cell it crosses."),
        }
        ds["conus_mask"].attrs = {
            "long_name": "1 where the cell centre lies inside the United States",
            "note": "Outside the mask every field is NaN, never 0 -- unobserved is not zero-risk.",
        }
        ds["lat"].attrs = {"units": "degrees_north", "axis": "Y"}
        ds["lon"].attrs = {"units": "degrees_east", "axis": "X"}

        ds.attrs = {
            "title": f"CONUS historical tornado hazard -- {description}",
            "hazard": "Tornado (severe convective storm)",
            # GLOBAL `units` as well as the per-variable one: the delivery plan and
            # layers.csv read ds.attrs["units"], so a layer that carries units only on its
            # variables reports "?" to the customer.
            "units": "crossings per 1e4 km2 per year",
            "source_dataset": "NOAA SPC severe weather database, actual tornadoes 1950-2025",
            "source_url": "https://www.spc.noaa.gov/wcm/",
            "source_licence": "US Government work -- public domain",
            "ingest_script": "scripts/download_tornado_spc.py",
            "processing_script": "scripts/process_tornado_spc.py",
            "processed_on": datetime.now(timezone.utc).strftime("%Y-%m-%d"),

            "output_contract": ("observational-historical-v1 -- NOT the OUTPUT-SPEC.md decadal "
                                "contract. No decade dimension, no ols_slope, no sen_slope, no "
                                "n_members, no n_models. Those are absent because they have no "
                                "meaning for a single observational record, NOT set to 1."),
            "field_nature": "point-process occurrence rate, zero-inflated, continuous where non-zero",
            "statistic": "rate = crossings / (area * years); Gamma(k+1/2, T) posterior quartiles",
            "statistic_rationale": (
                "MEASURED, not chosen. The pooled-median branch was computed first and is "
                "degenerate: within-cell quantiles of the annual count over OCCUPIED cells are "
                "q25=0, median=0, q75=0 at 0.25 deg (and at 0.5 deg annual); 91.6% of "
                "(cell,year) pairs are zero at 0.25 deg, 97.3% at 0.1 deg. Publishing that "
                "branch would have written (0,0,0) across the whole map."),

            "temporal_window": f"{args.start_year}-{args.end_year}",
            "n_years": n_years,
            "magnitude_rung": rung,
            "magnitude_threshold": "none" if threshold is None else f"mag >= {threshold}",
            "rung_family": "all | f1plus | f2plus | f3plus -- built in one pass; the "
                           "threshold is a judgement call and the comparison is the deliverable",

            "percentile_direction": "higher_is_worse",
            "estimator": args.estimator,
            "rate_is_resolution_dependent": (
                "THE VALUE IS A TRACK-INTERSECTION FREQUENCY PER UNIT AREA, NOT AN AREAL "
                "TORNADO INTENSITY, AND IT DOES NOT SURVIVE RESAMPLING. One track crossing "
                "several cells is counted in each, so the per-area figure rises as the grid "
                "is refined. Measured 2026-08-18: area-weighted aggregation of the four "
                "0.25 deg children onto their 0.5 deg parent reads +7.9% at the median cell "
                "(p90 +23%), and total crossings inflate 1.123x from 0.5 to 0.25 deg. Never "
                "resample this layer to another resolution -- recompute it there."),
            "overdispersion_measured": (
                "Variance-to-mean ratio of annual per-cell counts, 2026-08-18: median 0.99, "
                "p90 1.70; 16.4% of cells above 1.5 for all reports and 6.2% for F2+. The "
                "marginal Poisson assumption behind the interval therefore holds well at the "
                "typical cell and UNDERSTATES the interval in the clustered tail, where "
                "outbreak days put many tornadoes in one cell in one year."),
            "spatial_resolution": f"{res} degrees",
            "geometry_assignment": args.geometry,
            "records_without_endpoint": int(fallback),

            "resolution_caveat": (
                "SCREENING LAYER -- NOT SITE-SPECIFIC. A tornado damage path is ~100 m wide "
                f"against a ~{res*111:.0f} km cell. The value is the rate over the cell; the "
                "probability that a specific building is struck is far smaller, and is NOT the "
                "footprint-to-cell-area ratio -- it depends on swath width, path length and "
                "orientation, which this layer does not model. This layer ranks which sites merit "
                "investigation and cannot support a design load or an asset-level estimate."),
            "reporting_bias_caveat": (
                "TORNADO REPORTS MEASURE OBSERVING CAPABILITY AS WELL AS HAZARD. In this "
                "database, reports rose 2.8x from the 1950s to the 2020s, weak (F0) reports "
                "rose 7.7x, while strong (F2+) reports FELL 19% and show no modern trend. "
                "Report density tracks population, roads and radar coverage, so a rate built "
                "on all reports is biased toward populated areas -- which correlates with the "
                "exposure it is meant to be scored against. The f2plus rung is the "
                "reporting-stable subset and exists for this reason."),
            "coverage_caveat": (
                "CONUS ONLY. SPC observes the United States; cells outside the mask are NaN "
                "because they are UNOBSERVED, not because they are safe. No global tornado "
                "occurrence dataset with comparable detection exists -- Europe's pan-continental "
                "ESWD holds 9,563 tornadoes for 1800-2014 against 73,458 US records since 1950, "
                "and forbids commercial use. Do not deliver this layer for a non-US site."),
            "no_trend_rationale": (
                "No slope is emitted. NOT because 76 years cannot support one -- it can -- but "
                "because the fitted trend would measure the observing system (see "
                "reporting_bias_caveat). Do not add a slope to this layer without first "
                "resolving the inhomogeneity."),
            "unrated_reports_in_window": n_unrated,
            "unrated_handling": ("mag = -9 counted in the 'all' rung, excluded from every "
                                 "F-rung. ~204/yr in the 2020s, so F-rungs understate recent "
                                 "years; n_unrated is recorded so the size is auditable."),

            "cells_in_mask": in_mask,
            "cells_with_at_least_one_crossing": occupied,
            "zero_tier_share_of_mask": round(100.0 * (in_mask - occupied) / max(in_mask, 1), 2),
            "qa_reviewed_on": "null -- NOT YET REVIEWED BY A HUMAN",
        }

        out_dir = args.out_root / f"tornado-spc_{rung}_{label}"
        out_dir.mkdir(parents=True, exist_ok=True)
        # Scenario token is `observed`, NOT `historical`: delivery.NON_DELIVERABLE_SCENARIOS
        # excludes "historical" (it means the ISIMIP historical FORCING run there), so a
        # file named that way is silently dropped by discover_scenarios and the layer
        # would look empty rather than broken.
        out_path = out_dir / f"tornado-{rung}_observed_processed.nc"
        encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
        ds.to_netcdf(out_path, encoding=encoding)
        written.append(out_path)

        print(f"\n  [{rung}] {len(sub):,} records -> {occupied:,}/{in_mask:,} cells occupied "
              f"({100*occupied/in_mask:.1f}%), zero tier {ds.attrs['zero_tier_share_of_mask']}%")
        finite = rate[mask & np.isfinite(rate)]
        nz = finite[finite > 0]
        if nz.size:
            print(f"       rate over occupied cells: p50={np.median(nz):.3f} "
                  f"p90={np.quantile(nz,0.9):.3f} max={nz.max():.3f} per 1e4 km2 per yr")
        report_reference_sites(ds, rung)
        print(f"       -> {out_path}")

    print(f"\nWrote {len(written)} layers. NOT delivery-ready:")
    print("  * test_shared_baseline.py will reject these (different contract, by design)")
    print("  * qa_reviewed_on is null in every file -- a human must read the maps")
    print("  * config/layer_registry.yaml and config/hazard_taxonomy.yaml are untouched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
