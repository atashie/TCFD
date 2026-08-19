#!/usr/bin/env python3
"""Aggregate the Arup global landslide hazard map from 30 arcsec to a 0.25 deg layer.

PRODUCT: TCFD/CDP. This layer answers `observational-historical-v1`, NOT the
OUTPUT-SPEC.md decadal contract -- there is no decade axis, no slopes and no ensemble,
because the source is a single deterministic historical map. Verifier is
`scripts/test_observational_baseline.py`. Options review and the ISIMIP-absence
receipt: docs/landslide-data-options-2026-08-19.md.

===========================================================================
WHAT THE SOURCE IS
===========================================================================
World Bank / GFDRR "Global landslide hazard map", produced by Arup (2021), rainfall
trigger, median over 1980-2018. EPSG:4326, 43201 x 15779 at 30 arcsec (~1 km),
float32, covering land between 59.49 S and 72.00 N. The value is an annual frequency
of *significant* landslides per km^2 -- "significant" meaning large enough to have
been reported had it occurred in a populated place, broadly >100 m^2.

It is the only global landslide product that publishes a RATE with physical units.
That is the whole reason it can be aggregated at all: the areal mean of a per-km^2
frequency is arithmetic, whereas the mean of an ordinal susceptibility class is not a
quantity.

===========================================================================
THE STATISTIC WAS MEASURED BEFORE IT WAS WRITTEN, AND THE OBVIOUS DESIGN FAILS
===========================================================================
Each 0.25 deg cell contains 900 native 1 km pixels, so the "obvious" design is the
within-cell spatial quartiles of those 900 values. That was computed first. Over
OCCUPIED cells only (>=1 hazard-bearing pixel -- the most favourable denominator
available, 89,871 cells):

    branch                                      median==0   q25==q75   lower_ci<0
    ------------------------------------------  ---------   --------   ----------
    A  quantiles over all 900 pixels               78.48%     64.73%          --
    B  areal mean +/- 1 SD  (the repo's
       pooled_mean_zero_inflated branch)            0.00%         --      88.55%
    C  quartiles of 5x5 sub-block (0.05 deg)
       means                                       49.38%     33.38%          --
    D  quantiles over HAZARD-BEARING pixels         0.00%      5.70%       0.00%

Branch A is the documented `pooled_mean_zero_inflated` failure in a new setting:
three of four published variables would be identically 0 across two-thirds of the
occupied map, in a file that still passes a structural contract check. Rejected.

Branch B is non-degenerate in its central value and is what the shared statistics
module does for `let`, but a landslide FREQUENCY cannot be negative and `mean - SD`
goes below zero in 88.55% of occupied cells. Publishing an unphysical lower bound in
nine of ten cells is worse than publishing a different statistic. Rejected.

Branch C was built specifically to rescue B by finding a non-negative interval that
brackets the areal mean. It does not: the within-cell field is strongly right-skewed,
so the areal mean exceeds its own sub-block upper quartile in 42.48% of occupied
cells. No spatial partition can fix that -- it is a property of the distribution, not
of the block size. Rejected, and recorded so it is not re-attempted.

Branch D is adopted. All three of `median`, `lower_ci` and `upper_ci` are quantiles
of ONE distribution -- the rate over the cell's hazard-bearing 1 km pixels -- so the
triple is coherent by construction and needs no clipping to make it order. (The
tornado build's first attempt mixed an MLE central value with posterior quartiles and
had to clip; that incoherence is the specific trap this design avoids.)

    median    = 50th percentile of the rate over hazard-bearing pixels in the cell
    lower_ci  = 25th percentile of the same
    upper_ci  = 75th percentile of the same

READ IT AS: "on the ground in this cell that carries landslide hazard, this is the
rate and its spread." It is NOT a cell average, and it must be read together with
`hazard_area_fraction` -- a cell that is 1.4% hazard-bearing and one that is 42%
hazard-bearing (the inter-quartile range of that fraction over occupied cells) can
carry the same median.

===========================================================================
WHY THE PERCENTILE RANKS ON A DIFFERENT VARIABLE, DELIBERATELY
===========================================================================
`percentile` ranks on `areal_mean_rate` (the cos-latitude-weighted mean over ALL 900
pixels), not on `median`. This is a declared departure and it is the one place this
layer's score and its central value do not track each other.

The reason is measured. The two orderings are nearly independent -- Spearman rank
correlation 0.34 over occupied cells -- so the choice changes essentially every
customer score, and the reference sites say the areal mean is the defensible one:

    site                 hazard-  cond.    areal      pct on     pct on
                         bearing  median   mean       cond.med   areal mean
    -------------------  -------  -------  --------   --------   ----------
    Baguio PH              0.999  0.22788  0.229637      100.0        100.0
    Medellin CO            0.869  0.03198  0.024498       98.8         99.2
    Shimla IN              0.994  0.02528  0.019195       97.6         98.3
    Kathmandu NP           0.648  0.02139  0.009102       95.9         91.5
    Cusco PE               0.878  0.02000  0.010989       95.0         93.5
    Bergen NO              0.576  0.00957  0.011918       87.5         94.3
    Apennines IT           0.257  0.00418  0.002763       58.9         74.2
    Cairo EG               0.001  0.00312  0.000003        5.3          1.4
    Amsterdam NL           0.000      ---  0.000000        ---          ---

The conditional median discards EXTENT, which is a first-order determinant of whether
an asset is exposed. It puts the Apennines -- the most landslide-mapped terrain in
Europe -- at the 59th percentile, and it gives Cairo the 5th percentile off a single
Nile-escarpment pixel covering 0.1% of the cell. The areal mean, which integrates
extent and intensity, ranks both correctly. Since `percentile` is what drives the
customer-facing score, it ranks on the areal mean and `areal_mean_rate` is published
so the ranking is auditable rather than asserted.

Consequence a reader must be told: a cell can show `median` = 0 with a high
percentile. That is not an inconsistency -- it means most of the cell is flat while
the cell as a whole still carries substantial landslide activity.

NOTE ON THE SHARED VERIFIER: its two-tier consistency check is keyed on a variable
named `n_events`. This layer has no event count -- it is a modelled rate -- so the
count is published as `n_hazard_pixels` and that check does not fire. The equivalent
assertions are made here instead, and printed, so the tier is still checked.

===========================================================================
ZERO IS AMBIGUOUS IN THE SOURCE AND MUST BE DISAMBIGUATED
===========================================================================
The source COG declares NO nodata value and writes exact 0.0 for ocean AND for flat
land. Measured 2026-08-19: the Pacific, the Atlantic, the Sahara, the Amazon
floodplain, the Netherlands and Greenland all read exactly 0.000000.

For a landslide field, 0 on land is a legitimate result -- flat ground has no slope
failure. 0 on water is not a result at all. Publishing both as 0 would rank ocean as
"least hazardous" and would break the repo rule that unobserved cells are NaN and
never 0. So a cell is published only if it is inside the source extent AND
(land by the ISIMIP3b land-sea mask, upsampled 0.5 -> 0.25 deg, OR carries at least
one hazard-bearing pixel). The second clause exists so that no cell the source itself
modelled as hazard-bearing can be masked away by a coarser coastline.

STATE THE MASK WITH EVERY SHARE (CLAUDE.md). Shares in this docstring and in the
file's attributes are over the published mask unless they name another denominator;
`occupied` always means cells with >=1 hazard-bearing pixel.

===========================================================================
WHY THERE IS NO SLOPE, AND NO SCENARIO
===========================================================================
The source is a single map summarising 1980-2018. There is no annual series to fit,
so no slope is emitted -- correctly ABSENT rather than faked, per the contract. There
is also no scenario axis and no forward-looking information of any kind in this
layer. It answers "where has significant landsliding been frequent", not "where will
it become more frequent". The projection options are in the options doc.

USAGE
    python3 scripts/process_landslide_arup.py
    python3 scripts/test_observational_baseline.py data/processed/landslide-arup_rf-median_hist
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

try:
    import rasterio
    from rasterio.windows import Window
except ImportError:  # pragma: no cover
    print("rasterio is required: .venv/bin/python3 -m pip install rasterio", file=sys.stderr)
    raise

from scipy.stats import rankdata

SRC = Path("data/raw/landslide-arup/LS_RF_Median_1980-2018_COG.tif")
LANDMASK = Path("data/masks/ISIMIP3b_landseamask.nc4")
OUT_DIR = Path("data/processed/landslide-arup_rf-median_hist")
OUT_FILE = OUT_DIR / "landslide-arup_observed_processed.nc"

NLAT, NLON = 720, 1440          # standard global 0.25 deg grid, matching the tornado layer
F = 30                          # native pixels per 0.25 deg edge (0.25 / (1/120))
SRC_TOP_LAT = 72.0              # source top edge; the file reports 72.00006, a 6.7 m offset
ROW0 = int((90.0 - SRC_TOP_LAT) / 0.25)   # target row index (north-down) of the source top

SOURCE_URL = ("https://datacatalogfiles.worldbank.org/ddh-published/0037584/"
              "DR0045419/LS_RF_Median_1980-2018_COG.tif")


def aggregate(src_path: Path) -> dict[str, np.ndarray]:
    """One streaming pass over the COG, north-down, emitting the per-cell fields."""
    s = rasterio.open(src_path)
    if abs(s.res[0] - 1.0 / 120.0) > 1e-8:
        raise ValueError(f"expected 30 arcsec pixels, got {s.res}")
    ncol = (s.width // F) * F           # drop the 43201st wrap column -> 43200 = 1440 * 30
    if ncol // F != NLON:
        raise ValueError(f"column count {ncol} does not divide into {NLON} cells")

    keys = ("median", "lower_ci", "upper_ci", "areal_mean", "max", "n_hazard", "n_pixels")
    out = {k: np.full((NLAT, NLON), np.nan, np.float32) for k in keys}

    nblocks = int(np.ceil(s.height / F))
    for b in range(nblocks):
        r0 = b * F
        nrows = min(F, s.height - r0)
        j = ROW0 + b
        if j >= NLAT:
            break

        a = s.read(1, window=Window(0, r0, ncol, nrows))          # (nrows, ncol)
        a = a.reshape(nrows, NLON, F).transpose(1, 0, 2).reshape(NLON, nrows * F)

        # Area weight: cos(lat) varies only down rows, so build it per native row and
        # tile across the F columns of each cell.
        lat = SRC_TOP_LAT - (r0 + np.arange(nrows) + 0.5) / 120.0
        w = np.repeat(np.cos(np.deg2rad(lat)), F)[None, :]

        out["areal_mean"][j] = (a * w).sum(1) / w.sum()
        out["max"][j] = a.max(1)
        out["n_hazard"][j] = (a > 0).sum(1)
        out["n_pixels"][j] = a.shape[1]

        # Branch D: quantiles over hazard-bearing pixels only. Cells with none give
        # NaN here and are set to 0 below -- their rate really is 0 everywhere.
        hazard_only = np.where(a > 0, a, np.nan)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            out["lower_ci"][j] = np.nanquantile(hazard_only, 0.25, axis=1)
            out["median"][j] = np.nanquantile(hazard_only, 0.50, axis=1)
            out["upper_ci"][j] = np.nanquantile(hazard_only, 0.75, axis=1)

        if b % 100 == 0:
            print(f"    block {b:4d}/{nblocks}  target row {j}")

    s.close()
    return out


def load_land_mask() -> np.ndarray:
    """ISIMIP3b 0.5 deg land-sea mask, upsampled x2 onto the 0.25 deg grid, north-down."""
    ds = xr.open_dataset(LANDMASK)
    var = "mask" if "mask" in ds.data_vars else list(ds.data_vars)[0]
    m = ds[var]
    if m.ndim == 3:
        m = m.isel({m.dims[0]: 0})
    # CLAUDE.md: a masked array's .data still holds _FillValue. Fill to NaN BEFORE
    # comparing, or every cell reads as land.
    arr = np.ma.filled(np.ma.masked_invalid(m.values), np.nan).astype(float)
    land = np.isfinite(arr) & (arr > 0.5)
    if float(ds["lat"][0]) < float(ds["lat"][-1]):      # store north-down
        land = land[::-1]
    ds.close()
    return np.repeat(np.repeat(land, 2, axis=0), 2, axis=1)


def two_tier_percentile(rate: np.ndarray, valid: np.ndarray, occupied: np.ndarray) -> np.ndarray:
    """Zeros -> tier 1; hazard-bearing cells ranked 2-100 on the areal mean rate.

    Ties share a rank ('average'), which matters here: the source field is strongly
    quantised, so an argsort-based rank would split identical rates arbitrarily.
    """
    pct = np.full(rate.shape, np.nan, np.float32)
    pct[valid & ~occupied] = 1.0
    v = rate[occupied]
    if v.size:
        r = rankdata(v, method="average")            # 1..n, ties averaged
        pct[occupied] = (2.0 + 98.0 * (r - 0.5) / v.size).astype(np.float32)
    return pct


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT_FILE)
    args = ap.parse_args()

    if not SRC.exists():
        print(f"Missing {SRC} -- run scripts/download_landslide_arup.py first", file=sys.stderr)
        return 2

    print(f"Aggregating {SRC} -> 0.25 deg")
    f = aggregate(SRC)

    in_extent = np.isfinite(f["n_pixels"])
    occupied = in_extent & (f["n_hazard"] > 0)
    land = load_land_mask()
    valid = in_extent & (land | occupied)

    median, lower, upper = f["median"].copy(), f["lower_ci"].copy(), f["upper_ci"].copy()
    # A valid cell with no hazard-bearing pixel has a rate of 0 everywhere in it. That is
    # a genuine zero, not a gap, so it publishes 0 rather than NaN -- and it is exactly
    # the tier-1 population of the percentile.
    zero_cells = valid & ~occupied
    for arr in (median, lower, upper):
        arr[zero_cells] = 0.0

    areal = np.where(valid, f["areal_mean"], np.nan).astype(np.float32)
    frac = np.where(valid, f["n_hazard"] / f["n_pixels"], np.nan).astype(np.float32)
    mx = np.where(valid, f["max"], np.nan).astype(np.float32)
    nhz = np.where(valid, f["n_hazard"], np.nan).astype(np.float32)

    pct = two_tier_percentile(f["areal_mean"], valid, occupied)

    for arr in (median, lower, upper):
        arr[~valid] = np.nan

    # --- the tier assertions the shared verifier cannot make for this layer ---
    print("\n  two-tier checks (the shared verifier's version is keyed on `n_events`):")
    t1 = np.isfinite(pct) & (pct == 1.0)
    print(f"    tier-1 cells == cells with no hazard-bearing pixel : "
          f"{bool(np.array_equal(t1, zero_cells))}")
    occ_pct = pct[occupied]
    print(f"    occupied cells score within [2, 100]               : "
          f"{float(occ_pct.min()):.2f} - {float(occ_pct.max()):.2f}")
    rho = float(np.corrcoef(rankdata(f["areal_mean"][occupied]), rankdata(occ_pct))[0, 1])
    print(f"    percentile monotonic in areal_mean_rate            : rho = {rho:.6f}")

    # --- flip to ascending latitude to match every other layer ---
    flip = lambda a: a[::-1, :].copy()  # noqa: E731
    median, lower, upper, pct = map(flip, (median, lower, upper, pct))
    areal, frac, mx, nhz = map(flip, (areal, frac, mx, nhz))

    lat = -90.0 + 0.25 * (np.arange(NLAT) + 0.5)
    lon = -180.0 + 0.25 * (np.arange(NLON) + 0.5)

    n_valid = int(valid.sum())
    n_occ = int(occupied.sum())
    units = "landslides per km2 per year"

    ds = xr.Dataset(
        {
            "median": (("lat", "lon"), median, {
                "long_name": "Rainfall-triggered landslide frequency on hazard-bearing ground",
                "units": units,
                "note": ("50th percentile over the cell's hazard-bearing 1 km pixels. All three "
                         "of median/lower_ci/upper_ci are quantiles of THIS ONE distribution. "
                         "Read with hazard_area_fraction: this is not a cell average."),
            }),
            "lower_ci": (("lat", "lon"), lower, {
                "long_name": "25th percentile of the rate over hazard-bearing ground",
                "units": units,
                "note": "Spatial spread WITHIN the cell, not ensemble or sampling uncertainty.",
            }),
            "upper_ci": (("lat", "lon"), upper, {
                "long_name": "75th percentile of the rate over hazard-bearing ground",
                "units": units,
                "note": "Spatial spread WITHIN the cell, not ensemble or sampling uncertainty.",
            }),
            "percentile": (("lat", "lon"), pct, {
                "long_name": "Two-tier percentile of landslide hazard",
                "units": "1",
                "note": ("RANKS ON areal_mean_rate, NOT on median -- declared departure, see the "
                         "processor docstring. Cells with no hazard-bearing pixel -> 1; the rest "
                         "ranked 2-100 among cells that carry hazard. A cell can show median = 0 "
                         "with a high percentile: most of it is flat, the cell still carries "
                         "substantial landslide activity."),
            }),
            "areal_mean_rate": (("lat", "lon"), areal, {
                "long_name": "Cos-latitude-weighted mean rate over the whole cell",
                "units": units,
                "note": "The quantity `percentile` ranks on. Published so the ranking is auditable.",
            }),
            "hazard_area_fraction": (("lat", "lon"), frac, {
                "long_name": "Fraction of the cell's area carrying non-zero landslide hazard",
                "units": "1",
                "note": ("Inter-quartile range over occupied cells is 0.014-0.424, so this varies "
                         "by a factor of ~30 between cells with comparable medians."),
            }),
            "max_rate": (("lat", "lon"), mx, {
                "long_name": "Maximum 1 km pixel rate within the cell",
                "units": units,
            }),
            "n_hazard_pixels": (("lat", "lon"), nhz, {
                "long_name": "Count of hazard-bearing 1 km pixels in the cell",
                "units": "1",
                "note": ("NOT an event count -- deliberately not named n_events, which the shared "
                         "verifier reserves for observed occurrences."),
            }),
        },
        coords={
            "lat": ("lat", lat, {"units": "degrees_north", "axis": "Y"}),
            "lon": ("lon", lon, {"units": "degrees_east", "axis": "X"}),
        },
        attrs={
            "title": "Global rainfall-triggered landslide hazard, 1980-2018, aggregated to 0.25 deg",
            "hazard": "Landslide (rainfall-triggered)",
            "units": units,
            "source_dataset": ("World Bank / GFDRR Global landslide hazard map (Arup, 2021), "
                               "rainfall trigger, median 1980-2018, 30 arcsec"),
            "source_url": SOURCE_URL,
            "source_licence": ("World Bank / GFDRR, Global landslide hazard map (Arup, 2021). "
                               "The two publisher-side catalogue records disagree -- the WB DDH "
                               "record states CC BY-NC 4.0 and the energydata.info mirror states "
                               "CC-BY-4.0, and the project report states neither (checked "
                               "2026-08-19). Retained because it is a real ambiguity in the "
                               "source, not because it is unresolved for us."),
            "source_licence_status": ("CLEARED FOR OUR LIMITED COMMERCIAL USE -- user "
                                      "determination 2026-08-19. Attribution to World Bank / "
                                      "GFDRR and Arup is required wherever this layer is "
                                      "published or quoted."),
            "attribution_required": ("World Bank / GFDRR Global landslide hazard map, produced "
                                     "by Arup (2021). Must appear wherever a value from this "
                                     "layer is published."),
            "ingest_script": "scripts/download_landslide_arup.py",
            "processing_script": "scripts/process_landslide_arup.py",
            "processed_on": dt.date.today().isoformat(),
            "output_contract": (
                "observational-historical-v1 -- NOT the OUTPUT-SPEC.md decadal contract. No "
                "decade dimension, no ols_slope, no sen_slope, no n_members, no n_models. Those "
                "are ABSENT rather than faked; the source is one deterministic historical map."),
            "field_nature": ("modelled areal occurrence rate, continuous where non-zero, "
                             "extremely zero-inflated and strongly right-skewed"),
            "statistic": ("quantiles (25/50/75) of the rate over the cell's HAZARD-BEARING 1 km "
                          "pixels -- one distribution, all three slots"),
            "statistic_rationale": (
                "MEASURED, not chosen; four branches were computed before one was written. Over "
                "occupied cells (n=89,871): quantiles over all 900 pixels give median==0 in "
                "78.48% and q25==q75 in 64.73% (the documented pooled_mean_zero_inflated "
                "failure); areal mean +/- 1 SD is non-degenerate but puts lower_ci below zero in "
                "88.55%, which is unphysical for a frequency; quartiles of 5x5 sub-block means "
                "fail to bracket the areal mean in 42.48% because the within-cell field is "
                "right-skewed, so no spatial partition rescues that branch. Quantiles over "
                "hazard-bearing pixels are degenerate in only 5.70% and never negative."),
            "percentile_basis": (
                "areal_mean_rate, NOT median. Spearman rank correlation between the two "
                "orderings is 0.34 over occupied cells, so this changes essentially every score. "
                "Reference sites decide it: the conditional median puts the Apennines at the "
                "58.9th percentile against 74.2nd on the areal mean, and gives Cairo 5.3 off a "
                "single pixel covering 0.1% of the cell against 1.4."),
            "percentile_direction": "higher_is_worse",
            "temporal_window": "1980-2018",
            "n_years": 39,
            "spatial_resolution": "0.25 degrees",
            "native_resolution": "30 arcsec (~1 km); 900 native pixels per published cell",
            "aggregation": ("cos-latitude-weighted; one streaming pass, no resampling of the "
                            "source values"),
            "land_mask": ("ISIMIP3b land-sea mask upsampled 0.5 -> 0.25 deg, UNION with any cell "
                          "carrying a hazard-bearing pixel so a coarser coastline cannot mask "
                          "away ground the source modelled as hazardous"),
            "cells_published": n_valid,
            "cells_occupied": n_occ,
            "zero_tier_share_of_mask": round(100.0 * (n_valid - n_occ) / n_valid, 2),
            "source_zero_is_ambiguous": (
                "The source COG declares no nodata and writes exact 0.0 for ocean AND for flat "
                "land -- Pacific, Atlantic, Sahara, Amazon, Netherlands and Greenland all read "
                "0.000000 (measured 2026-08-19). Disambiguated by the land mask above."),
            "median_branch_measured": (
                "median==0 in 78.48% and q25==q75 in 64.73% of occupied cells. Recorded as a "
                "measurement, not a caveat, so it cannot be promoted into a customer report."),
            "resolution_caveat": (
                "SCREENING LAYER -- NOT SITE-SPECIFIC. Slope stability turns on metres of "
                "terrain; a 0.25 deg cell is ~28 km and holds 900 km^2 of ground. The published "
                "value describes the hazard-bearing part of the cell, and the hazard-bearing "
                "part is a median 10.1% of it (inter-quartile range 1.4%-42.4%). This ranks "
                "which sites merit a geotechnical investigation and cannot support a design "
                "decision, a slope-stability judgement or an asset-level loss estimate."),
            "reporting_bias_caveat": (
                "THE SOURCE MODEL IS TRAINED ON REPORTED LANDSLIDES, so it inherits their "
                "geography. Global landslide inventories over-represent English-language and "
                "higher-GDP regions, and 'significant' is itself defined by reportability. "
                "Absence of hazard in a sparsely reported region is weaker evidence than absence "
                "in a densely reported one."),
            "coverage_caveat": (
                "RAINFALL-TRIGGERED ONLY, AND LAND BETWEEN 59.5 S AND 72.0 N. Earthquake-"
                "triggered landslides are published by the same project as a separate layer and "
                "are deliberately excluded here, because seismic triggering is stationary under "
                "every emissions pathway and mixing it in would attribute a tectonic hazard to a "
                "climate driver. Cells outside the mask are NaN because they are outside the "
                "modelled domain, never because they are safe."),
            "no_trend_rationale": (
                "Slopes are absent by contract. The source is a single map summarising 1980-2018 "
                "with no annual series behind it, so there is nothing to fit; a slope here would "
                "have to be invented rather than estimated."),
            "no_scenario_rationale": (
                "Historical only. This layer answers where significant landsliding has been "
                "frequent, and carries no forward-looking information. Projection options are in "
                "docs/landslide-data-options-2026-08-19.md."),
            "qa_reviewed_on": "null -- NOT YET REVIEWED BY A HUMAN",
        },
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in ds.data_vars}
    ds.to_netcdf(args.out, encoding=enc)

    print(f"\n  published cells : {n_valid:,}")
    print(f"  occupied cells  : {n_occ:,}  ({100.0 * n_occ / n_valid:.2f}% of the published mask)")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
