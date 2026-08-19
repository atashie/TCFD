#!/usr/bin/env python3
"""Historical baseline frequency of very large hail (>= 5 cm), 0.25 deg global.

PRODUCT: TCFD/CDP, `observational-historical-v1` (see OUTPUT-SPEC.md), NOT the decadal
contract. One panel, no decade axis, no slopes, no ensemble. Verified by
`scripts/test_observational_baseline.py`; `test_shared_baseline.py` rejects it by design.

WHY A NON-ISIMIP SOURCE. ISIMIP cannot express this hazard, not merely lack it: every hail
environment index needs CAPE, deep-layer shear and a freezing level, and the ISIMIP
bias-adjusted forcing publishes 11 SURFACE variables and nothing aloft, so the shear term
cannot be formed at all. Enumerated 2026-08-18 across all four rounds -- see
`config/isimip_search_catalog.yaml` and the WORKFLOW-ISSUES entry of that date.

SOURCE. Battaglioli, Taszarek, Groenemeijer, Pucik & Radler, "Contrasting trends in very
large hail events and related economic losses across the globe", Nature Geoscience 19,
52-58 (2026), doi 10.1038/s41561-025-01868-0. AR-CHaMo -- an additive logistic hazard model
that predicts thunderstorm occurrence and then hail conditional on a thunderstorm -- applied
to ERA5 every 3 hours on the native 0.25 deg grid for all 74 years, 1950-2023. Over 80
billion profiles.

WHICH COPY, AND WHY IT MATTERS. The same fields exist in two places under DIFFERENT
LICENCES, and only one of them is usable here:

  * Zenodo 10.5281/zenodo.17064885 -- 74 annual fields plus climatology and trends, and
    it is CC BY-NC-ND 4.0 on all three versions. NonCommercial rules it out for this
    product; NoDerivatives separately rules out regridding or aggregating it.
  * The article's Source Data -- CC BY 4.0, and NOT figure summaries: each file carries
    the full 813,600-cell grid. `MOESM7` is the climatology, cell-for-cell identical to the
    Zenodo file; `MOESM2` is the trend and its p-value.

This processor reads the CC BY copy. That choice is the reason the baseline is a 74-year
mean rather than a recent-decade window: the annual fields exist only under NC-ND. The
licence question is documented at `reports/maps/hail-vlh/essl_licence_query.md` and
is deliberately parked.

A UNITS TRAP IN THE SOURCE, MEASURED. The climatology field is described by its own
metadata as "Mean annual number of hail >= 5 cm events per ERA5 grid box (1950-2023)". It is
not. It is the 74-year TOTAL: summing all 74 annual fields reproduces it exactly, ratio
1.0000 with p1 = p99 = 1.0000 across 332,164 cells with a non-trivial baseline. Taken at
face value it overstates the annual rate 74-fold. The paper's own text gives the global
maximum as "around 0.50 events per year ... a VLH event every 2 years per grid box", which
is the divided figure; this file's maximum is 48.83 total = 0.66 per year at 31.75S 65.0W,
north-west of Cordoba, Argentina. Every published value here is divided by 74.

WHAT THE VALUE IS. `median` is the expected number of >= 5 cm hail events per 0.25 deg grid
box per year. It is a MODELLED rate, not a count of observations, and it is per GRID BOX --
so it does not survive resampling and it is not comparable between latitudes without the
area normalisation carried alongside as `events_per_1e4km2_per_year` (a 0.25 deg box is
~770 km2 at the equator and ~310 km2 at 60 deg).

THE INTERVAL IS SAMPLING, NOT MODEL UNCERTAINTY. `lower_ci`/`upper_ci` are the 25th and
75th percentiles of a Gamma(k + 1/2, T) posterior with k the modelled expected count over
the window and T = 74 years -- the same estimator, and the same Jeffreys prior, as the
sibling `tornado-spc` layers, so the two observational layers are read the same way. It
answers "given this modelled rate, how much would the realised count vary over 74 years",
NOT "how uncertain is AR-CHaMo". The latter is not quantified anywhere in the source and is
larger: the model is calibrated on reports from Europe, the United States and Australia and
extrapolated everywhere else.

REFERENCE SITES (GUARDRAILS 12), from the published field, events per box per year:
    Cordoba / Mendoza, Argentina  0.66   global maximum, matches the paper's hotspot
    NE South Africa               > 0.30  the paper's second hotspot, matches
    US Great Plains, S Saskatchewan/Alberta, tri-border UY/PY/BR  all non-trivial
A field that is zero across northern Argentina would be wrong regardless of contract checks.

USAGE
    python3 scripts/process_hail_vlh_battaglioli.py
    python3 scripts/test_observational_baseline.py data/processed/hail-vlh-archamo_ge5cm
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.stats import gamma as gamma_dist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.land_mask import get_isimip_landmask  # noqa: E402

N_YEARS = 74
WINDOW = "1950-2023"
CELL = 0.25
EARTH_R_KM = 6371.0088

REFERENCE_SITES = {
    "Cordoba, Argentina": (-31.4, -64.2),
    "Mendoza, Argentina": (-32.9, -68.8),
    "Johannesburg, South Africa": (-26.2, 28.0),
    "Oklahoma City, USA": (35.5, -97.5),
    "Calgary, Canada": (51.0, -114.1),
    "Milan, Italy": (45.5, 9.2),
    "Brisbane, Australia": (-27.5, 153.0),
    "Beijing, China": (39.9, 116.4),
}


def read_source_grid(csv_path: Path, columns: list[str]):
    """Read a CC BY Source Data CSV into (lat, lon, {column: grid}).

    The files are one row per grid cell in an unspecified order, so the values are placed
    by coordinate rather than by reshape -- a reshape would be right only if the row order
    happened to be raster order, and nothing in the file promises that.
    """
    lats: list[float] = []
    lons: list[float] = []
    cols: dict[str, list[float]] = {c: [] for c in columns}
    with csv_path.open(newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        idx = {}
        for c in columns:
            if c not in header:
                raise SystemExit(f"{csv_path} has no column {c!r}; header is {header}")
            idx[c] = header.index(c)
        for row in rdr:
            lats.append(float(row[0]))
            lons.append(float(row[1]))
            for c in columns:
                v = row[idx[c]]
                cols[c].append(np.nan if v == "" else float(v))
    lat_u = np.array(sorted(set(lats)))
    lon_u = np.array(sorted(set(lons)))
    li = {v: i for i, v in enumerate(lat_u)}
    oi = {v: i for i, v in enumerate(lon_u)}
    ii = np.array([li[v] for v in lats])
    jj = np.array([oi[v] for v in lons])
    out = {}
    for c in columns:
        g = np.full((lat_u.size, lon_u.size), np.nan)
        g[ii, jj] = np.array(cols[c])
        out[c] = g
    return lat_u, lon_u, out


def to_global(lat_src: np.ndarray, lon_src: np.ndarray, grid: np.ndarray):
    """Pad the source domain onto a full global grid on the SAME lattice.

    The source covers 57S-84N only. Padding rather than reprojecting keeps every value
    exactly as published; the pad is NaN, which is what "outside the analysed domain"
    means and is never 0.
    """
    lat_g = np.round(np.arange(-90.0, 90.0, CELL), 4)          # 720, ascending
    lon_g = np.round(np.arange(-180.0, 180.0, CELL), 4)        # 1440
    out = np.full((lat_g.size, lon_g.size), np.nan)
    lat_idx = np.searchsorted(lat_g, np.round(lat_src, 4))
    lon_idx = np.searchsorted(lon_g, np.round(lon_src, 4))
    if not (np.all(np.isclose(lat_g[lat_idx], lat_src)) and np.all(np.isclose(lon_g[lon_idx], lon_src))):
        raise SystemExit("source coordinates are not on the global 0.25 deg lattice")
    out[np.ix_(lat_idx, lon_idx)] = grid
    return lat_g, lon_g, out, (lat_src.min(), lat_src.max())


def land_mask_on(lat_g: np.ndarray, lon_g: np.ndarray) -> np.ndarray:
    """ISIMIP3b land-sea mask, nearest-neighbour onto the 0.25 deg lattice.

    `.filled(np.nan)` is not optional: a masked array's `.data` still holds `_FillValue`,
    and comparing that directly once marked every cell on Earth as land.
    """
    path = get_isimip_landmask(version="3b")
    with xr.open_dataset(path, mask_and_scale=True) as ds:
        name = [v for v in ds.data_vars if "mask" in v.lower()] or list(ds.data_vars)
        da = ds[name[0]].squeeze()
        arr = np.ma.filled(np.ma.masked_invalid(da.values), np.nan)
        mlat = da["lat"].values.astype(float)
        mlon = da["lon"].values.astype(float)
    if mlat[0] > mlat[-1]:                      # the ISIMIP mask is published N->S
        mlat = mlat[::-1]
        arr = arr[::-1, :]
    i = np.clip(np.searchsorted(mlat, lat_g) - 0, 0, mlat.size - 1)
    i = np.array([int(np.argmin(np.abs(mlat - v))) for v in lat_g])
    j = np.array([int(np.argmin(np.abs(mlon - v))) for v in lon_g])
    up = arr[np.ix_(i, j)]
    return np.isfinite(up) & (up > 0.5)


def cell_area_km2(lat_g: np.ndarray, nlon: int) -> np.ndarray:
    """Area of each 0.25 deg cell, which varies by a factor of ~2.5 across the domain."""
    lat_edges = np.deg2rad(np.stack([lat_g - CELL / 2, lat_g + CELL / 2]))
    band = (2 * np.pi * EARTH_R_KM ** 2 * np.abs(np.sin(lat_edges[1]) - np.sin(lat_edges[0]))) / (360.0 / CELL)
    return np.repeat(band[:, None], nlon, axis=1)


def two_tier_percentile(rate: np.ndarray, expected: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Cells the model gives NO very large hail -> 1. The rest ranked 2..100 among themselves.

    Keyed on the modelled expected count, not on the posterior rate being zero: under a
    Gamma posterior the rate is small-but-positive everywhere, so `rate == 0` would select
    nothing and the tier would be a label rather than a partition. Identical in form to
    `process_tornado_spc.two_tier_percentile`, so the two observational layers tier the
    same way and can be read side by side.
    """
    pct = np.full(rate.shape, np.nan)
    inmask = mask & np.isfinite(rate)
    zero = inmask & (expected <= 0)
    nonzero = inmask & (expected > 0)
    pct[zero] = 1.0
    vals = rate[nonzero]
    if vals.size:
        order = vals.argsort()
        ranks = np.empty(vals.size, dtype=float)
        ranks[order] = np.arange(1, vals.size + 1, dtype=float)
        # average ranks for ties, so heavy ties cannot manufacture a spurious ordering
        uniq, inv, counts = np.unique(vals, return_inverse=True, return_counts=True)
        sums = np.zeros(uniq.size)
        np.add.at(sums, inv, ranks)
        ranks = (sums / counts)[inv]
        pct[nonzero] = 2.0 + 98.0 * (ranks - 1) / max(vals.size - 1, 1)
    return pct


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=Path("data/raw/hail-vlh-battaglioli"))
    ap.add_argument("--out-dir", type=Path, default=Path("data/processed/hail-vlh-archamo_ge5cm"))
    ap.add_argument("--qa-reviewed-on", default="null -- NOT YET REVIEWED BY A HUMAN",
                    help="date a human read the QA maps. Defaults to null: the review is of the "
                         "DATA, so a rebuild must not inherit a sign-off it has not earned.")
    args = ap.parse_args()

    src = args.raw_dir / "MOESM7.csv"
    if not src.exists():
        print(f"missing source: {src}", file=sys.stderr)
        return 2

    lat_s, lon_s, cols = read_source_grid(src, ["hail_ge5cm"])
    total = cols["hail_ge5cm"]
    print(f"source grid {total.shape}, lat {lat_s.min()}..{lat_s.max()}, lon {lon_s.min()}..{lon_s.max()}")
    print(f"  window total: max {np.nanmax(total):.4f} events over {N_YEARS} years "
          f"= {np.nanmax(total)/N_YEARS:.4f} per year")

    lat_g, lon_g, total_g, domain = to_global(lat_s, lon_s, total)
    land = land_mask_on(lat_g, lon_g)
    in_domain = np.isfinite(total_g)
    coverage = land & in_domain
    print(f"  coverage: {coverage.sum():,} cells (land AND inside {domain[0]:.0f}..{domain[1]:.0f} lat)")

    expected = np.where(coverage, total_g, np.nan)            # events over the whole window
    alpha = expected + 0.5                                    # Jeffreys prior
    with np.errstate(invalid="ignore"):
        med = gamma_dist.ppf(0.50, a=alpha, scale=1.0 / N_YEARS)
        lo = gamma_dist.ppf(0.25, a=alpha, scale=1.0 / N_YEARS)
        hi = gamma_dist.ppf(0.75, a=alpha, scale=1.0 / N_YEARS)
    for a in (med, lo, hi):
        a[~coverage] = np.nan

    pct = two_tier_percentile(med, expected, coverage)
    area = cell_area_km2(lat_g, lon_g.size)
    per_area = np.where(coverage, med / area * 1e4, np.nan)

    zero_share = 100.0 * np.mean(expected[coverage] <= 0)
    print(f"  zero tier: {zero_share:.2f}% of covered cells have no modelled VLH at all")
    print(f"  rate: median {np.nanmedian(med):.5f}, p99 {np.nanpercentile(med, 99):.4f}, "
          f"max {np.nanmax(med):.4f} events per box per year")

    # Trend field is INGESTED FOR CONTEXT ONLY and is not published -- see no_trend_rationale.
    trend_note = "not read"
    tp = args.raw_dir / "MOESM2.csv"
    if tp.exists():
        tlat, tlon, tcols = read_source_grid(tp, ["trend_hail_ge5cm", "pvalue_hail_ge5cm"])
        if not (np.array_equal(tlat, lat_s) and np.array_equal(tlon, lon_s)):
            raise SystemExit("the trend file is on a different grid from the climatology")
        _, _, tr, _ = to_global(tlat, tlon, tcols["trend_hail_ge5cm"])
        _, _, pv, _ = to_global(tlat, tlon, tcols["pvalue_hail_ge5cm"])
        # THE DENOMINATOR IS THE POINT. Over every finite cell the answer is dominated by the
        # vast low-rate interior where nothing can be detected; over cells that actually get
        # very large hail it is a different number with a different sign balance. Both are
        # stated, each with its mask named, because a share whose denominator is unnamed is
        # not a measurement.
        occupied = coverage & (expected > 0)
        strong = coverage & (expected / N_YEARS > 0.1)
        parts = []
        for label, m in (("covered land cells with any modelled VLH", occupied),
                         ("covered land cells above 0.1 events/box/yr", strong)):
            ok = m & np.isfinite(tr) & np.isfinite(pv)
            sig = ok & (pv < 0.05)
            parts.append(f"over {label} (n={int(ok.sum()):,}): {100*np.mean(sig[ok]):.1f}% significant "
                         f"at p<0.05, of which {100*np.mean(tr[sig] > 0):.1f}% positive")
        trend_note = ("MEASURED from the published trend field (MOESM2; NOT republished here) -- "
                      + "; ".join(parts) + ". Where a trend is significant this 74-year mean is "
                      "displaced from present-day conditions, and this layer cannot say by how much.")
    print(f"  {trend_note}")

    ds = xr.Dataset(
        {
            "median": (("lat", "lon"), med.astype("float32")),
            "lower_ci": (("lat", "lon"), lo.astype("float32")),
            "upper_ci": (("lat", "lon"), hi.astype("float32")),
            "percentile": (("lat", "lon"), pct.astype("float32")),
            "n_events": (("lat", "lon"), expected.astype("float32")),
            "events_per_1e4km2_per_year": (("lat", "lon"), per_area.astype("float32")),
            "coverage_mask": (("lat", "lon"), coverage.astype("int8")),
        },
        coords={"lat": lat_g, "lon": lon_g},
    )
    ds["median"].attrs = {"units": "events per 0.25 deg grid box per year",
                          "long_name": "Expected annual number of hail >= 5 cm events"}
    for _v in ("lower_ci", "upper_ci"):   # separate dicts, not one shared object
        ds[_v].attrs = {"units": "events per 0.25 deg grid box per year",
                        "long_name": f"{_v} of the Gamma(k+1/2, T) posterior on the rate"}
    ds["n_events"].attrs = {"units": "events per 0.25 deg grid box per 74 years",
                            "long_name": "Modelled expected count over the window -- NOT observed"}
    ds["events_per_1e4km2_per_year"].attrs = {"units": "events per 1e4 km2 per year"}
    ds["coverage_mask"].attrs = {"long_name": "1 where land AND inside the analysed 57S-84N domain"}

    ds.attrs = {
        "title": "Very large hail (>= 5 cm) historical baseline frequency, global 0.25 deg",
        "hazard": "Hail (severe convective storm)",
        "units": "events per 0.25 deg grid box per year",
        "source_dataset": ("Battaglioli et al., Nature Geoscience 19, 52-58 (2026), "
                           "AR-CHaMo applied to ERA5 3-hourly, 1950-2023 -- Source Data MOESM7 (CC BY 4.0)"),
        "source_url": "https://doi.org/10.1038/s41561-025-01868-0",
        "source_licence": ("CC BY 4.0 via the article's Source Data. The Zenodo deposit of the same "
                           "fields (10.5281/zenodo.17064885) is CC BY-NC-ND 4.0 and is NOT used."),
        "ingest_script": "curl of the article Source Data -- see data/raw/hail-vlh-battaglioli/",
        "processing_script": "scripts/process_hail_vlh_battaglioli.py",
        "output_contract": ("observational-historical-v1 -- NOT the OUTPUT-SPEC.md decadal contract. "
                            "No decade dimension, no ols_slope, no sen_slope, no n_members, no n_models."),
        "field_nature": "modelled occurrence rate, zero-inflated, continuous where non-zero",
        "statistic": "Gamma(k + 1/2, T) posterior median, k = modelled expected count over the window, T = 74 yr",
        "statistic_rationale": (
            "The same estimator as the sibling tornado-spc observational layers, so the two are read "
            "the same way. It is a SAMPLING interval on the realised count given the modelled rate, "
            "not an uncertainty on AR-CHaMo itself, which the source does not quantify and which is "
            "larger. The source field is the 74-YEAR TOTAL despite its own metadata calling it a mean "
            "annual number -- verified by summing all 74 annual fields, ratio 1.0000 across 332,164 "
            "cells -- so every value here is that total divided by 74."),
        "percentile_direction": "higher_is_worse",
        "percentile_basis": ("two-tier over COVERED cells only: cells the model gives no VLH score 1, "
                             "the rest are rank-scored 2..100 among themselves. Not a percentile of "
                             "all land, and 50 means median among cells that get very large hail."),
        "estimator": "posterior",
        "n_events_semantics": (
            "n_events is the MODELLED expected count over the window, not a count of observed "
            "events -- there is no observation in this layer at all. It carries the name the "
            "observational contract's verifier keys the percentile zero tier on, and the tier "
            "here means 'AR-CHaMo produces no very large hail in this cell', not 'nobody saw any'."),
        "temporal_window": WINDOW,
        "n_years": N_YEARS,
        "hail_size_threshold": ">= 5 cm diameter (very large hail)",
        "grid_lattice_note": (
            "0.25 deg on the ERA5 lattice (points at exact 0.25 multiples), latitude ASCENDING to match "
            "the tornado-spc observational layers. This is offset by 0.125 deg from those layers' "
            "cell-centre lattice, so the two grids share no coordinate values -- never union their masks."),
        "spatial_resolution_degrees": CELL,
        "trend_context_measured": trend_note,
        "resolution_caveat": (
            "SCREENING LAYER -- NOT SITE-SPECIFIC. A hail swath is a few km across against a ~28 km "
            "cell, so the value is the rate over the cell and not the probability that a given address "
            "is struck. The value is also PER GRID BOX and therefore resolution-dependent: it does not "
            "survive resampling, and cell area falls from ~770 km2 at the equator to ~310 km2 at 60 deg, "
            "which is why events_per_1e4km2_per_year is carried alongside for cross-latitude comparison."),
        "reporting_bias_caveat": (
            "AR-CHaMo IS CALIBRATED WHERE THE REPORTS ARE. Its lightning and conditional-hail components "
            "were trained on Europe, the United States and Australia; South America, Africa and most of "
            "Asia are extrapolation, and the paper itself flags equatorial Africa and the Sahel as places "
            "where sparse ground reports make the estimate hard to verify. A high value in an untrained "
            "region is a model statement, not an observation."),
        "coverage_caveat": (
            "DOMAIN IS 57S-84N AND LAND ONLY HERE. The source excludes maritime locations more than "
            "100 km from a coast, and this layer further restricts to land via the ISIMIP3b land-sea "
            "mask. Cells outside are NaN because they are UNANALYSED, never because they are safe."),
        "threshold_caveat": (
            "5 cm IS THE EXTREME TAIL, NOT THE DAMAGE THRESHOLD. Photovoltaic glass, roofing and vehicle "
            "panels are damaged from roughly 2-3 cm, so this layer counts a small subset of damaging "
            "hail. It ranks where the most violent hail occurs; it does not measure how often an asset "
            "is hit by hail that can hurt it."),
        "no_trend_rationale": (
            "No slope is emitted, for two reasons that are both external to the statistics. First, the "
            "user's decision of 2026-08-19: this layer is a current-conditions baseline and the trend is "
            "a separate product question. Second, the annual fields needed to fit a trend ourselves are "
            "CC BY-NC-ND and therefore unusable here, so any slope would have to be republished from the "
            "source's own trend field, which is a different thing from a slope we computed. What the "
            "source's trend says about this baseline is recorded in trend_context_measured."),
        "baseline_window_caveat": (
            "THE WINDOW IS THE WHOLE RECORD, 1950-2023, NOT A RECENT DECADE. It stands in for current "
            "conditions because the annual fields that would allow a 2014-2023 window are under a "
            "licence this product cannot use. Where the source's trend is significant this mean is "
            "displaced from present-day conditions -- see trend_context_measured for how much of the "
            "domain that affects."),
        "processed_on": "2026-08-19",
        "qa_reviewed_on": args.qa_reviewed_on,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "hail-vlh_observed_processed.nc"
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(out, encoding=enc)
    print(f"\nwrote {out} ({out.stat().st_size/1e6:.1f} MB)")

    print("\nREFERENCE SITES (events per box per year) -- GUARDRAILS 12:")
    for name, (la, lo_) in REFERENCE_SITES.items():
        i = int(np.argmin(np.abs(lat_g - la)))
        j = int(np.argmin(np.abs(lon_g - lo_)))
        v = med[i, j]
        p = pct[i, j]
        print(f"  {name:<28} {v:8.4f}   percentile {p:6.1f}"
              + ("   <-- ZERO/UNCOVERED" if not np.isfinite(v) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
