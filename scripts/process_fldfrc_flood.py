"""Process fldfrc (CaMa-Flood annual flooded area fraction) into the TCFD 8-value-class
format -- THREE PARALLEL LAYERS, one per flood-protection level.

fldfrc = the fraction of each grid cell inundated per year, from the CaMa-Flood v3.6.2
hydrodynamic model driven by ISIMIP2b GHM discharge (Willner's `flood-processing`,
doi 10.5281/zenodo.1241051). Unlike every other flood representation in ISIMIP this is a
genuine AREA share, which is why it was chosen -- see config/isimip_search_catalog.yaml
search_results.flooding for the full options review.

SOURCE: ISIMIP2b/DerivedOutputData/Zimmer2023, ingested and coarsened 150 arcsec -> 0.5 deg
by scripts/download_fldfrc_flood.py (area-preserving 12x12 block aggregation; see that
module for why the coarsening happens at ingest and how it is made auditable).

Ensemble: 24 members per scenario per protection level = 6 GHMs x 4 GCMs, complete.

    GHMs   clm45, cwatm, h08, lpjml, matsiro, watergap2   (clm50/mpi-hm/pcr-globwb
                                                           are historical-only)
    GCMs   gfdl-esm2m, hadgem2-es, ipsl-cm5a-lr, miroc5
    scen   rcp26, rcp60, rcp85     <- all three; the Lange2020 `ler` alternative has no rcp85
    soc    2005soc / co2, UNIFORM across all 216 source files (no soc/sens divergence
           to reconcile -- contrast csoil's mixed-CO2 ensemble, GUARDRAILS S9)

WHY THREE LAYERS AND NOT ONE
----------------------------
The protection level is not a parameter detail; it is the single biggest choice in this
layer, and it changes the direction of the story. Measured on one member
(h08/miroc5/rcp60), global flooded area, 2020s -> 2090s decadal means:

    none      6,104,969 -> 6,282,904 km2/yr    +2.9%   <- saturates: no risk signal
    100yr       260,894 ->   501,993 km2/yr   +92.4%
    flopros     743,970 -> 1,130,053 km2/yr   +51.9%

They answer three different questions and none is a substitute for another:

  none     every inundation event, no defenses credited. Largest absolute area, but it
           counts routine seasonal floodplain wetting that recurs every year in both
           decades, so the metric SATURATES and its trend is muted.
  100yr    only events rarer than 1-in-100 inundate. Best read as a SEVERITY THRESHOLD
           rather than "protection": it isolates severe floods WITHOUT relying on any
           real-world defense database.
  flopros  events below the local FLOod PROtection Standards threshold are contained --
           residual risk given actual defenses.

Protectiveness runs none < flopros < 100yr globally. "flopros ~= none where undefended"
is FALSE and was measured, not assumed: flopros retains only 19-36% of the undefended
signal even in the least-defended regions (Bangladesh 0.189, Niger/Nigeria 0.364,
Mekong 0.348), and 0.8-4% in well-defended ones (Netherlands 0.008, Mississippi 0.015,
Japan 0.043). FLOPROS is spatially real, though -- it is LESS protective than a uniform
100yr in Niger and Bangladesh and MORE protective in the Mississippi and Netherlands,
where defenses exceed 1-in-100. See WORKFLOW-ISSUES.md 2026-07-28.

PROCESSING DECISIONS
--------------------
statistic          decadal MEAN of the annual flooded fraction = expected annual flooded
                   area share. Never the median: this is an exposed-AREA frequency
                   measure, the same rule as the Lange 2020 exposure family. At 0.5 deg
                   a single year is ~93% exact zeros, so a decadal median would be 0
                   nearly everywhere and destroy the field.
trend              BASELINE-ANCHORED rate, not a within-decade slope. Annual flooded area
                   swings 17x between adjacent decades in one member (flopros, h08/miroc5:
                   2080 = 1,207,332 km2 vs 2100 = 69,751 km2), so an annual regression
                   inside a decade is noise. GUARDRAILS S10.
percentile         two-tier (zero-inflated), zeros -> 1, positives ranked against the
                   non-zero 2020s baseline -> [2,100]. higher_is_worse: flooding is a
                   hazard, so no inversion (contrast csoil, where stored carbon is an
                   asset).
normalization      none. All 24 members are the same physical quantity in the same unit,
                   produced by ONE hydrodynamic model (CaMa-Flood) fed by 6 GHMs, so the
                   inter-member spread is genuine GHM+GCM uncertainty and becomes the CI
                   ("model democracy").
spatial_smoothing  none. 24 members per scenario is thick (contrast `let`, 1 model x 4
                   GCMs, which needed 5x5 smoothing).
ci                 mean +/- 1 inter-member SD, clamped to [0, 100] %. The clamp is SAFE
                   here and is NOT the burntarea mistake: flooded fraction is a true
                   bounded share of cell area, so median <= 100 always and
                   min(median+sd, 100) >= median. Burnt area is CUMULATIVE and legitimately
                   exceeds 100%, which is why that layer leaves upper_ci unbounded.

COVERAGE
--------
62,066 cells carry data = 128.8 million km2 = 95.7% of global land excluding Antarctica
(Lange2020 `led` has 67,420 land cells; 61,546 shared). This is a NORMAL ISIMIP land
mask, not a sparse floodplain subset -- 79.6% of domain cells are >=99% inside the
CaMa-Flood domain and the median is 100%. Cells outside it stay NaN; nothing is
zero-filled. `floodplain_fraction` is emitted so a partly-covered cell is auditable
rather than merely looking low.

Because the coarsening divides by the FULL cell area, sum(median * cell_area) over any
region is exactly the flooded area in km2 -- the field is directly aggregatable, which
is what "area impacted per year" has to mean.

Output: {layer}/data/fldfrc_{scenario}_processed.nc with
{median, percentile, trend, lower_ci, upper_ci, n_members, n_models,
 floodplain_fraction} on (decade, lat, lon), units % of cell area.

Usage:
    python scripts/process_fldfrc_flood.py [--protection none,100yr,flopros] [--no-publish]
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.contact_sheet import render_contact_sheet  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402
from utils.trend_significance import (  # noqa: E402
    METHOD as SIGNIFICANCE_METHOD, TREND_METHOD, AnnualEnsembleMean, mk_expanding,
    significance_definition, theilsen_decadal, trend_definition_decadal,
)

VAR = "fldfrc"
PROTECTIONS = ["none", "100yr", "flopros"]
RAW_PATTERN = "*_fldfrc_*_halfdeg_global_annual_*.nc"
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
VMAX = 100.0                      # a fraction of cell area cannot exceed 100% -- a real
                                  # ceiling, unlike burntarea's cumulative percentage
TWO_TIER_ZERO_THRESHOLD = 0.02
EXPECTED_MEMBERS = 24

PROTECTION_MEANING = {
    "none": ("no flood protection assumed -- EVERY inundation event contributes, "
             "including routine seasonal floodplain wetting. Largest absolute area but "
             "the metric SATURATES, so its climate trend is muted (+2.9% globally by "
             "the 2090s on a test member vs +52-92% for the thresholded levels)."),
    "100yr": ("only events rarer than 1-in-100-years inundate. Best read as a SEVERITY "
              "THRESHOLD rather than protection: it isolates severe floods without "
              "relying on any real-world defense database. Most protective of the three "
              "in most regions."),
    "flopros": ("events below the local FLOod PROtection Standards (FLOPROS) threshold "
                "are contained -- residual risk given ACTUAL defenses. Spatially real: "
                "less protective than a uniform 100yr in Niger/Bangladesh, more "
                "protective in the Mississippi/Netherlands. NOT a proxy for "
                "'unprotected' -- it retains only 19-36% of the undefended signal even "
                "in the least-defended regions."),
}


def log(msg):
    print(msg, flush=True)


def layer_id(protection):
    return f"riverflood_fldfrc-{protection}_annual"


def parse_name(fpath):
    """Parse a coarsened member filename FROM THE END.

    cama-flood_{ghm}_{gcm}_ewembi_{rcp}_2005soc_co2_fldfrc_{protection}_halfdeg_global_annual_{y0}_{y1}.nc
              [-13]  [-12]          [-10]                    [-6]

    Parsed from the end because 'cama-flood' already contains a hyphen and any future
    prefix change would shift positive indices. Model/GCM names never contain '_'.
    """
    p = os.path.basename(fpath)[:-3].split("_")
    return dict(model=p[-13], gcm=p[-12], scenario=p[-10], protection=p[-6],
                member=f"{p[-13]}_{p[-12]}")


def load_member(fpath):
    """Load one coarsened member as (year, lat, lon) fraction, plus its domain mask."""
    ds = xr.open_dataset(fpath)
    da = ds[VAR].load()
    fpf = ds["floodplain_fraction"].load().values.astype(np.float32)
    ds.close()
    return da, fpf


def decade_mean_map(da, decade):
    """Mean over the decade window -> per-member expected annual flooded fraction."""
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return da.isel(year=sel).mean("year").values.astype(np.float32)


def anchored_trend(med_stack, i, b_idx):
    """trend[i] = (median[i] - median[baseline]) / (i - baseline), % per decade.

    Anchored at the 2020s origin, so each panel is the trend FROM the baseline TO that
    decade and is exactly the change map / elapsed decades -- spatially coherent by
    construction. A within-decade slope of the annual series is NOT used: annual flooded
    area swings ~17x between adjacent decades, which yields a spotty sign-flipping field.
    The baseline decade has no elapsed change -> 0, identical across scenarios so the
    shared baseline stays bit-identical.
    """
    span = i - b_idx
    if span == 0:
        return np.zeros(med_stack.shape[1:], np.float32)
    return ((med_stack[i] - med_stack[b_idx]) / span).astype(np.float32)


def make_pct_fn(baseline_flat):
    """Percentile-of-score against the 2020s baseline; higher = worse.

    Flooded area is strongly zero-inflated even as a decadal mean (large land areas
    never flood), so the two-tier branch is expected: exact zero -> 1, positives ranked
    against the NON-ZERO baseline -> [2,100]. A single-tier rank would put most of the
    world in the bottom percentile and crush the gradient where it matters.
    """
    frac_zero = float(np.mean(baseline_flat == 0.0))
    two_tier = frac_zero >= TWO_TIER_ZERO_THRESHOLD

    if two_tier:
        nz_sort = np.sort(baseline_flat[baseline_flat > 0])
        n_nz = len(nz_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            res = np.ones(vals.shape, np.float32)
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="left") / n_nz
                res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
            out[fin] = res
            return out.reshape(arr.shape)
    else:
        all_sort = np.sort(baseline_flat)
        n = len(all_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            if n > 0:
                frac = np.searchsorted(all_sort, vals, side="right") / n
                out[fin] = np.clip(100.0 * frac, 1.0, 100.0).astype(np.float32)
            return out.reshape(arr.shape)

    return pct, ("two_tier" if two_tier else "single_tier"), frac_zero


def cell_area_km2(lats, dlat=0.5, dlon=0.5, R=6371.0):
    half = np.deg2rad(dlat) / 2
    la = np.deg2rad(lats)
    return (R ** 2) * np.deg2rad(dlon) * (np.sin(la + half) - np.sin(la - half))


def process_protection(protection, publish=True):
    lid = layer_id(protection)
    log("=" * 78)
    log(f"fldfrc protection={protection}  ->  layer {lid}")
    log("=" * 78)

    files = [str(p) for p in storage.stage_raw(lid, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no members in s3://{storage.BUCKET}/{storage.raw_prefix(lid)}")
        log("Run scripts/download_fldfrc_flood.py first.")
        return None

    meta = {f: parse_name(f) for f in files}
    wrong = {f: m["protection"] for f, m in meta.items() if m["protection"] != protection}
    if wrong:
        raise ValueError(f"{len(wrong)} staged file(s) have the wrong protection level: "
                         f"{sorted(set(wrong.values()))} in prefix for {protection}")
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models}")
    log(f"GCMs: {gcms}")

    stage = storage.staging_dir(lid)
    out_dir = stage / "data"

    da0, fpf0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    A = cell_area_km2(lats)
    del da0

    # ---- Pass 1: per-member decadal-mean maps (as % of cell area) ---------------
    dec = {s: {} for s in scenarios}
    fpf_stack = []
    # trend_pvalue is tested on the ensemble-mean ANNUAL series (GUARDRAILS S15). Held
    # in % of cell area, matching the decadal maps below. `trend` uses only the decadal
    # medians (S10) -- vital here: a single year is ~93% exact zeros, so a median of
    # pairwise ANNUAL slopes would be exactly 0 almost everywhere.
    annual_acc = {s: AnnualEnsembleMean(DECADES[0],
                                        DECADES[-1] + WINDOW_YEARS - 1, (LAT, LON))
                  for s in scenarios}
    for f in sorted(files):
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da, fpf = load_member(f)
        annual_acc[s].add(da.year.values, da.values * 100.0)
        fpf_stack.append(fpf)
        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        for i, d in enumerate(DECADES):
            m = decade_mean_map(da, d)
            if m is not None:
                maps[i] = m * 100.0                     # fraction -> % of cell area
        dec[s][member] = maps
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            glob_km2 = float(np.nansum(maps[0] / 100.0 * A[:, None]))
        log(f"  {info['model']:<10} {info['gcm']:<13} {s}  years={da.year.size:<3d} "
            f"2020s global flooded area={glob_km2:>12,.0f} km2/yr")
        del da

    for s in scenarios:
        if len(dec[s]) != EXPECTED_MEMBERS:
            log(f"  WARNING: scenario {s} has {len(dec[s])} members, expected "
                f"{EXPECTED_MEMBERS} -- the ensemble is incomplete")

    # floodplain_fraction is static per member; they agree, so take the max as the
    # domain (a cell in ANY member's domain is in the layer's domain).
    # The CaMa-Flood river network is model-independent, so every member SHOULD carry
    # the same domain -- but that is an assumption, and an assumption about a land mask
    # is exactly the kind that has bitten this project before (csoil: 58,714-67,647
    # cells across 5 models). Measure it and say so; take the union as the layer domain
    # and let n_members carry per-cell depth.
    fpf_arr = np.stack(fpf_stack, 0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # cells NaN in every member -> all-NaN slice
        fpf_layer = np.nanmax(fpf_arr, axis=0).astype(np.float32)
    per_member_cells = np.isfinite(fpf_arr).sum(axis=(1, 2))
    if per_member_cells.min() != per_member_cells.max():
        log(f"  NOTE: member domains DIFFER: {per_member_cells.min():,}-"
            f"{per_member_cells.max():,} cells; union = {int(np.isfinite(fpf_layer).sum()):,}")
    else:
        log(f"  member domains identical across all {len(fpf_stack)} files: "
            f"{per_member_cells[0]:,} cells")

    # ---- Shared 2020s baseline, identical across scenarios ---------------------
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec[s]})
    shared_2020 = []
    for member in all_members:
        per_scen = [dec[s][member][b_idx] for s in scenarios if member in dec[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_2020.append(np.nanmean(np.stack(per_scen, 0), axis=0))
    shared_2020 = np.stack(shared_2020, 0).astype(np.float32)

    # ---- Per-member contact sheet (GUARDRAILS S11) -----------------------------
    sheet_path = None
    try:
        sheet_path = render_contact_sheet(
            {m: shared_2020[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", lid, BASELINE_DECADE, units="% of cell area",
            note=(f"{len(all_members)} members (6 GHMs x 4 GCMs), protection="
                  f"{protection}. All share one hydrodynamic model (CaMa-Flood v3.6.2), "
                  f"so panels should agree on WHERE rivers are and differ mainly in "
                  f"magnitude; a panel disagreeing on river GEOMETRY is a defect. "
                  f"Coarsened 150arcsec -> 0.5 deg, area-preserving."))
        log(f"  contact sheet: {sheet_path}  <-- REVIEW before trusting the layer")
    except Exception as e:                                          # noqa: BLE001
        log(f"  WARNING: contact sheet failed ({type(e).__name__}: {e}); the "
            f"per-member visual check of GUARDRAILS S11 has NOT been produced")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared_median = np.nanmean(shared_2020, axis=0)
        shared_sd = np.nanstd(shared_2020, axis=0)
    shared_lo = np.clip(shared_median - shared_sd, 0, None)
    shared_hi = np.clip(shared_median + shared_sd, None, VMAX)
    shared_nmem = np.sum(np.isfinite(shared_2020), axis=0).astype(np.float32)
    member_models = np.array([m.split("_")[0] for m in all_members])
    shared_nmodel = np.zeros_like(shared_nmem)
    for mdl in np.unique(member_models):
        shared_nmodel += np.any(np.isfinite(shared_2020[member_models == mdl]),
                                axis=0).astype(np.float32)
    land = shared_nmem > 0
    n_land = int(land.sum())
    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)

    base_area = float(np.nansum(shared_median / 100.0 * A[:, None]))
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, land cells "
        f"n={n_land:,}, exact-zero fraction={frac_zero:.2%}, percentile={pct_mode}")
    log(f"  global flooded area (ensemble mean) = {base_area:,.0f} km2/yr")
    log(f"  coverage: all-member {100*np.mean(shared_nmem[land]==len(all_members)):.2f}%, "
        f">=2 models {100*np.mean(shared_nmodel[land]>=2):.2f}%, "
        f"single-model {100*np.mean(shared_nmodel[land]==1):.2f}%")

    # ---- Phase A: per-scenario median / CI / percentile ------------------------
    med_by, lo_by, hi_by, pct_by, nmem_by, nmodel_by, mem_by = {}, {}, {}, {}, {}, {}, {}
    for s in scenarios:
        members = sorted(dec[s])
        mem_by[s] = members
        mdl_of = np.array([m.split("_")[0] for m in members])
        stack = np.stack([dec[s][m] for m in members], 0)
        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        lower, upper, percentile = (np.full_like(median, np.nan) for _ in range(3))
        nmem, nmodel = np.zeros_like(median), np.zeros_like(median)
        for i, d in enumerate(DECADES):
            layer = stack[:, i, :, :]
            nmem[i] = np.sum(np.isfinite(layer), axis=0)
            for mdl in np.unique(mdl_of):
                nmodel[i] += np.any(np.isfinite(layer[mdl_of == mdl]), axis=0)
            if d == BASELINE_DECADE:
                median[i], lower[i], upper[i], percentile[i] = (
                    shared_median, shared_lo, shared_hi, shared_pct)
                nmem[i], nmodel[i] = shared_nmem, shared_nmodel
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(layer, axis=0)
                sd = np.nanstd(layer, axis=0)
            lower[i] = np.clip(median[i] - sd, 0, None)
            # Safe clamp: flooded fraction is a bounded share, so median <= VMAX and
            # min(median+sd, VMAX) >= median. Not the burntarea trap.
            upper[i] = np.clip(median[i] + sd, None, VMAX)
            percentile[i] = pct(median[i])
        nmem[nmem == 0] = np.nan
        nmodel[~np.isfinite(nmem)] = np.nan
        med_by[s], lo_by[s], hi_by[s], pct_by[s] = median, lower, upper, percentile
        nmem_by[s], nmodel_by[s] = nmem, nmodel

    # ---- Phase B: Theil-Sen trend + significance + write -----------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    for s in scenarios:
        log(f"\n--- {s} ---")
        median, lower, upper, percentile = med_by[s], lo_by[s], hi_by[s], pct_by[s]
        # Theil-Sen on the DECADAL medians (S10); p-value on the ANNUAL series (S15).
        trend = theilsen_decadal(median, DECADES, window_years=WINDOW_YEARS,
                                 baseline_decade=BASELINE_DECADE).astype(np.float32)
        years_a, mean_annual = annual_acc[s].result()
        tpval, ttau, tnobs = mk_expanding(years_a, mean_annual, DECADES,
                                          window_years=WINDOW_YEARS,
                                          baseline_decade=BASELINE_DECADE)
        for arr in (trend, tpval, ttau):
            arr[~np.isfinite(median)] = np.nan
        tnobs[~np.isfinite(median)] = 0
        for i, d in enumerate(DECADES):
            area = float(np.nansum(median[i] / 100.0 * A[:, None]))
            log(f"  {d}s: global flooded area={area:>12,.0f} km2/yr "
                f"({100*(area/base_area-1):+6.1f}% vs 2020s)  "
                f"mean trend={np.nanmean(trend[i]):+.5f} %/dec")

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "trend_pvalue": (["decade", "lat", "lon"], tpval.astype(np.float32)),
                "trend_tau": (["decade", "lat", "lon"], ttau.astype(np.float32)),
                "trend_n_obs": (["decade", "lat", "lon"], tnobs.astype(np.float32)),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
                "n_members": (["decade", "lat", "lon"], nmem_by[s]),
                "n_models": (["decade", "lat", "lon"], nmodel_by[s]),
                "floodplain_fraction": (["lat", "lon"], fpf_layer),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "protection_level": protection,
                "protection_meaning": PROTECTION_MEANING[protection],
                "long_name": ("Annual river-flooded area (percent of grid cell), "
                              f"protection={protection}"),
                "units": "%",
                "statistic": "decadal_mean_annual_flooded_area_percent",
                "value_note": (
                    "median = equal-weighted ensemble MEAN over 6 GHMs x 4 GCMs of the "
                    "annual flooded area share of each grid cell, in % of cell area. "
                    "MEAN not median: this is an exposed-AREA frequency measure (same "
                    "rule as the Lange 2020 exposure family), and at 0.5 deg a single "
                    "year is ~93% exact zeros so a decadal median would be 0 nearly "
                    "everywhere. Because the source coarsening divides by the FULL cell "
                    "area, sum(median/100 * cell_area) over any region is exactly the "
                    "flooded area in km2 -- the field is directly aggregatable."),
                "source_grid": (
                    "native 150 arcsec (4320 x 8640) CaMa-Flood output, coarsened to the "
                    "0.5 deg ISIMIP grid by area-weighted 12x12 block aggregation "
                    "(exact alignment: block centres match ISIMIP centres to 1.4e-14 deg; "
                    "area-conserving to <1e-9 relative error). Each 0.5 deg value "
                    "therefore carries real sub-grid information rather than a flag."),
                "normalization": (
                    "none -- all 24 members are the same physical quantity in the same "
                    "unit from ONE hydrodynamic model (CaMa-Flood v3.6.2) fed by 6 GHMs, "
                    "so the spread is genuine GHM+GCM uncertainty and is retained as the "
                    "CI (model democracy)."),
                "spatial_smoothing": "none (24 members per scenario is thick)",
                "trend_definition": trend_definition_decadal(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "trend_method": TREND_METHOD,
                "significance_method": SIGNIFICANCE_METHOD,
                "significance_definition": significance_definition(
                    DECADES, window_years=WINDOW_YEARS,
                    baseline_decade=BASELINE_DECADE),
                "significance_pooling": (
                    "flat mean across members within each year; the p-value is "
                    "tested on that ANNUAL series while `trend` is fitted on the "
                    "DECADAL medians -- a single year is ~93%% exact zeros, so an "
                    "annual Theil-Sen would be exactly 0 almost everywhere (S10)"),
                "trend_units": "% of cell area decade-1",
                "ci_definition": (
                    "lower/upper_ci = ensemble mean -/+ 1 inter-member standard deviation "
                    "across the 24 GHM x GCM members, clamped to [0, 100] %. The upper "
                    "clamp is SAFE here -- flooded fraction is a true bounded share of "
                    "cell area, so median <= 100 and min(median+sd, 100) >= median. "
                    "(Contrast the burntarea layer, where annual burnt area is CUMULATIVE, "
                    "legitimately exceeds 100%, and clamping the upper CI would push it "
                    "below the median.)"),
                "coverage_note": (
                    "62,066 cells carry data = 128.8 million km2 = 95.7% of global land "
                    "excluding Antarctica, against 67,420 land cells for Lange2020 `led` "
                    "(61,546 shared). This is a NORMAL ISIMIP land mask, not a sparse "
                    "floodplain subset: 79.6% of domain cells lie >=99% inside the "
                    "CaMa-Flood domain and the median is 100%. Cells outside it are NaN "
                    "(no model data) and are NOT zero-filled. floodplain_fraction gives "
                    "the area share of each cell inside the domain, so a partly-covered "
                    "cell is auditable rather than merely looking low. n_members / "
                    "n_models give per-cell ensemble depth."),
                "percentile_baseline": (
                    f"{pct_mode}: ranked against the 2020s ensemble-mean land "
                    "distribution (two-tier reserves 1 for exact-zero cells and maps "
                    "positives to [2,100])."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "known_latitude_seams": "60.0",
                "known_issues": (
                    "SHARP LEVEL DISCONTINUITY AT 60.0 N, INHERITED FROM THE SOURCE. The "
                    "zonal-mean flooded fraction roughly HALVES across a single 0.5 deg row "
                    "boundary at exactly 60.0 N -- 11.06% just north vs 5.93% just south "
                    "(1.87x) in the NATIVE 150 arcsec data, so it is not an artefact of the "
                    "coarsening; it is the largest single-row jump anywhere between 45 and "
                    "80 N and it appears identically in every GHM. The cause is CaMa-Flood's "
                    "floodplain topography changing DEM at the SRTM/HydroSHEDS coverage "
                    "limit (SRTM spans 60 N-56 S), so absolute flooded area north of 60 N is "
                    "inflated by roughly 1.9x relative to south of it and is NOT comparable "
                    "across the boundary. Any zonal or regional aggregate straddling 60 N "
                    "(Scandinavia, Russia, Canada, Alaska) mixes two topographic datasets. "
                    "IMPORTANT -- this is a LEVEL bias, not a trend bias: because the trend "
                    "and change maps difference the same cells against their own 2020s "
                    "baseline, the static offset largely cancels, and the change field is "
                    "spatially coherent across the seam (rcp85 none, 2020s->2090s: 60-70 N "
                    "+2.2%, 50-60 N -0.8%, 30-50 N +0.4%, tropics +8.3%, S subtropics "
                    "+6.4%). Use trend/change across 60 N with confidence; compare absolute "
                    "levels across it only with this caveat in hand."),
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "window_years": WINDOW_YEARS,
                "n_members": len(mem_by[s]),
                "n_models": len(models),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_sens": ("2005soc / co2, uniform across all 216 source files -- no "
                             "soc/sens divergence to reconcile."),
                "source_dataset": (
                    "ISIMIP2b DerivedOutputData/Zimmer2023 (fldfrc, 150arcsec). "
                    "CaMa-Flood v3.6.2 driven by ISIMIP2b GHM discharge; GEV fit to "
                    "picontrol (fit_soc=1860soc), discharge_threshold 0.1mm/d, produced by "
                    "github.com/swillner/flood-processing, doi 10.5281/zenodo.1241051."),
                "round_note": (
                    "ISIMIP2b / CMIP5 / RCP. Deliberate: no SSP-era fractional "
                    "flood-AREA product exists. The only ISIMIP3b flood product "
                    "(Heinicke2026 floodedarea) is a BINARY occurrence flag despite "
                    "declaring long_name='Exposed Area Share' -- verified on 6 members "
                    "across 3 models, 4 GCMs and 3 SSPs -- and binarising at 0.5 deg "
                    "overstates flooded area 17-28x. Its mask also covers 94.7% of the "
                    "globe including ocean. See config/isimip_search_catalog.yaml "
                    "search_results.flooding."),
                "sibling_layers": (
                    "riverflood_fldfrc-none_annual / -100yr_annual / -flopros_annual are "
                    "three PARALLEL layers over identical members, differing only in the "
                    "flood-protection assumption. They answer different questions and are "
                    "not substitutes; compare them rather than picking one blindly."),
                "description": (
                    f"River-flood area processed to the TCFD 8-value-class format with a "
                    f"shared 2020s baseline; 24-member (6 GHM x 4 GCM) ensemble, "
                    f"protection={protection}, rcp26/60/85, no normalization, no spatial "
                    f"smoothing."),
            },
        )
        path = out_dir / f"{VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    if not publish:
        log(f"\n--no-publish: staged only at {stage}")
        return dict(layer_id=lid, version=None, stage=str(stage))

    # ---- Provenance chain back to the 150 arcsec originals ---------------------
    # STORAGE.md: inputs.files[].source_url + checksum in layer.json is what makes
    # storage.cleanup_raw safe. It matters more than usual here, because the raw prefix
    # holds the COARSENED 0.5 deg fields rather than what ISIMIP served -- so the
    # manifest is the only place the chain back to the original is written down. Each
    # ingested file carries its own source_url and the sha512 OF THE ORIGINAL in its
    # global attrs; lift them into the manifest.
    raw_entries = []
    for f in files:
        p = Path(f)
        with xr.open_dataset(p) as rds:
            ra = rds.attrs
        raw_entries.append({
            "name": p.name,
            "bytes": p.stat().st_size,
            "sha256": storage.sha256_file(p),
            "source_url": ra.get("source_url", ""),
            "source_sha512": ra.get("source_sha512", ""),
            "source_bytes": ra.get("source_size_bytes", ""),
            "transform": ra.get("transform", ""),
        })
    missing = [e["name"] for e in raw_entries if not e["source_url"]]
    if missing:
        log(f"  WARNING: {len(missing)} raw file(s) carry no source_url -- "
            f"cleanup_raw will refuse to delete this prefix, which is the correct "
            f"behaviour but means re-ingest cannot be verified: {missing[:3]}")
    else:
        log(f"  provenance: {len(raw_entries)} raw inputs, all with source_url + "
            f"original sha512 recorded")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        lid, stage,
        raw_entries=raw_entries,
        # Reprocessing on an unchanged tree resolves to the same version id, so let it
        # BUMP rather than error: the superseded version stays as history and
        # _VERSION.json records the chain (STORAGE.md; isimip-process-visualize skill).
        on_exists="bump",
        created_by="scripts/process_fldfrc_flood.py",
        notes=(f"ISIMIP2b Zimmer2023 fldfrc, protection={protection}, rcp26/60/85; "
               f"24 members/scenario (6 GHMs x 4 GCMs); CaMa-Flood 150arcsec coarsened "
               f"area-preserving to 0.5 deg. One of three parallel protection-level "
               f"layers. See WORKFLOW-ISSUES.md 2026-07-28."),
    )
    log("\nGenerating QA/QC report and maps...")
    res = finalize_layer(lid, version=version,
                         extra_maps=[sheet_path] if sheet_path else None)
    return dict(layer_id=lid, version=version, fit=res.get("fit_for_use"),
                qa=(res.get("qa") or {}).get("summary"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protection", default=",".join(PROTECTIONS))
    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args()
    want = [p.strip() for p in args.protection.split(",") if p.strip()]
    bad = set(want) - set(PROTECTIONS)
    if bad:
        raise SystemExit(f"unknown protection level(s): {sorted(bad)}")

    results = []
    for p in want:
        results.append(process_protection(p, publish=not args.no_publish))

    log("\n" + "=" * 78)
    log("SUMMARY")
    log("=" * 78)
    for r in results:
        if not r:
            continue
        qa = r.get("qa") or {}
        log(f"  {r['layer_id']:<38} {r.get('version') or 'staged'}  "
            f"QA={qa.get('verdict','-')} fit_for_use={r.get('fit')}")


if __name__ == "__main__":
    main()
