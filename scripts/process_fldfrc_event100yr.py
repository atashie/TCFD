"""Process the 1-in-100-year river-flood EVENT into the TCFD 6-value-class format.

THE FOURTH FLOOD LAYER, and the only one that is not an expected-annual-area. The three
`riverflood_fldfrc-{none,100yr,flopros}_annual` layers all answer "how much area is
flooded in an average year, given protection level P". This one answers the two questions
a 1-in-100-year flood product actually needs:

    HOW OFTEN is the (preindustrial) 1-in-100-year flood exceeded?   <- primary
    HOW MUCH AREA does it cover when it is?                          <- companion

WHY IT CANNOT BE `none` + `100yr`
---------------------------------
The obvious construction -- add the undefended annual area to the 1-in-100 residual --
double-counts, because `100yr` is a SUBSET of `none`, not its complement. Measured per
cell per year at native 150 arcsec (h08/miroc5/rcp60):

    100yr > none                          0.000% of cells, every year tested
    where 100yr > 0, 100yr == none        84-88% of cells (mean ratio 0.90-0.93)
    none + 100yr > 100% of a cell         6,867 -> 28,139 cells   <- impossible for an
                                                                     area fraction

`100yr` is a FILTERED COPY of `none`: the same inundation field, retained only in years
when the event overtopped a 1-in-100 defence (with a partial-reduction component in the
other ~13% of active cells). `none` is therefore already the total undefended extent --
the MAXIMUM of the three layers, not a floor to add to. The sum is also numerically
unstable: in the 2020s it lands 17% BELOW the correctly measured footprint, and the error
drifts with scenario because the two terms scale completely differently (`none` +5%,
`100yr` +226% at rcp85).

THE CORRECT CONSTRUCTION
------------------------
`fldfrc_100yr[cell, year] > 0` flags exactly the years in which the threshold was
exceeded, and `fldfrc_none[cell, year]` records that event's full undefended extent. So,
per member and decade:

    exceedance_frequency[cell] = years with 100yr > 0  /  valid years        -> % of years
    event_footprint[cell]      = mean of `none` over those years             -> % of cell

Both are then ensemble-averaged over the 24 (GHM x GCM) members.

WHICH IS THE PRIMARY VALUE, AND WHY
-----------------------------------
`median` carries the EXCEEDANCE FREQUENCY (% of years), not the footprint. The footprint
is nearly static -- it grows only +1.8% from the 2020s to the 2090s at rcp85, because the
extent of a given-magnitude flood is set by topography -- while the frequency more than
doubles (+122%): the preindustrial 1-in-100-year flood becomes roughly a 1-in-6-year
flood. Making the footprint primary would produce a confidently flat layer, the same trap
the `none` protection level falls into. The risk signal is the frequency.

This also explains the other three layers: `none` saturates (+5.1% at rcp85) because it
counts every flood and topography fixes the extent, while `100yr` explodes (+225.5%)
because it counts only the rare tail and the rare tail is becoming ordinary. All four
layers are measuring a frequency change; only this one states it directly.

CAVEAT CARRIED IN THE OUTPUT ATTRS
----------------------------------
The threshold is a GEV fit to **picontrol** (`fit_scenario='picontrol'`, `fit_soc=
'1860soc'`), so "1-in-100" is PREINDUSTRIAL, not present-day. Already by the 2020s the
mean exceedance frequency over cells that exceed at all is ~1-in-12 years. That figure is
also biased upward by a selection effect (cells never exceeding in a decade are excluded
from the mean), so the 2020s->2090s RATIO is the robust part; the absolute level is
indicative. Both caveats are recorded in `known_issues`.

INPUTS COME FROM TWO SIBLING LAYERS' RAW PREFIXES
-------------------------------------------------
This layer derives from `riverflood_fldfrc-none_annual` and
`riverflood_fldfrc-100yr_annual` raw staging, matched on (model, gcm, scenario). Its own
raw prefix stays empty by design; `raw_entries` records the inputs from both siblings so
layer.json still carries the full provenance chain back to the 150 arcsec originals.

Output: fldfrc_{scenario}_processed.nc with
{median, percentile, trend, lower_ci, upper_ci, n_members, n_models,
 event_footprint, event_footprint_lower_ci, event_footprint_upper_ci,
 n_members_footprint, floodplain_fraction} on (decade, lat, lon).

Usage:
    python scripts/process_fldfrc_event100yr.py [--no-publish]
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
# Reuse the sibling processor's primitives rather than re-deriving them: the percentile
# tiering, the anchored-trend definition and the cell-area weighting must match the other
# three flood layers exactly, or the four are not comparable.
from process_fldfrc_flood import (  # noqa: E402
    make_pct_fn, anchored_trend, cell_area_km2, parse_name,
    DECADES, BASELINE_DECADE, WINDOW_YEARS, VMAX, EXPECTED_MEMBERS,
)

VAR = "fldfrc"
LAYER_ID = "riverflood_fldfrc-event100yr_annual"
SRC_UNPROTECTED = "riverflood_fldfrc-none_annual"
SRC_THRESHOLD = "riverflood_fldfrc-100yr_annual"
RAW_PATTERN = "*_fldfrc_*_halfdeg_global_annual_*.nc"


def log(msg):
    print(msg, flush=True)


def member_key(path):
    m = parse_name(path)
    return (m["model"], m["gcm"], m["scenario"])


def load_pair(f_none, f_100):
    """Return (none_annual, exceeded_mask, years, floodplain_fraction) for one member."""
    dn = xr.open_dataset(f_none)
    dh = xr.open_dataset(f_100)
    n = dn[VAR].load()
    h = dh[VAR].load()
    if not np.array_equal(n.year.values, h.year.values):
        raise ValueError(f"year axes differ between {os.path.basename(f_none)} and "
                         f"{os.path.basename(f_100)}")
    fpf = dn["floodplain_fraction"].load().values.astype(np.float32)
    yrs = n.year.values
    na = n.values.astype(np.float32)
    ha = h.values.astype(np.float32)
    dn.close()
    dh.close()
    # A year counts as valid only where BOTH fields are present, so the frequency
    # denominator and the footprint numerator rest on the same year set.
    valid = np.isfinite(na) & np.isfinite(ha)
    exceeded = valid & (ha > 0)
    return na, exceeded, valid, yrs, fpf


def decade_stats(na, exceeded, valid, yrs, decade):
    """(frequency % of years, footprint % of cell area) for one member and decade."""
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None, None
    ex = exceeded[sel]
    va = valid[sel]
    n_valid = va.sum(axis=0)
    n_exc = ex.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        freq = np.where(n_valid > 0, n_exc / np.maximum(n_valid, 1), np.nan) * 100.0
        # footprint: mean of the UNDEFENDED extent over exceeding years only
        fp = np.where(n_exc > 0,
                      np.where(ex, na[sel], 0.0).sum(axis=0) / np.maximum(n_exc, 1),
                      np.nan) * 100.0
    return freq.astype(np.float32), fp.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args()

    log("=" * 78)
    log(f"1-in-100-year flood EVENT -> layer {LAYER_ID}")
    log("=" * 78)

    files_n = {member_key(str(p)): str(p)
               for p in storage.stage_raw(SRC_UNPROTECTED, RAW_PATTERN)}
    files_h = {member_key(str(p)): str(p)
               for p in storage.stage_raw(SRC_THRESHOLD, RAW_PATTERN)}
    if not files_n or not files_h:
        log(f"ERROR: need raw for BOTH {SRC_UNPROTECTED} and {SRC_THRESHOLD}. "
            f"Found {len(files_n)} / {len(files_h)}. Run download_fldfrc_flood.py first.")
        return None

    keys = sorted(set(files_n) & set(files_h))
    only_n, only_h = sorted(set(files_n) - set(files_h)), sorted(set(files_h) - set(files_n))
    if only_n or only_h:
        log(f"  NOTE: {len(only_n)} member(s) only in {SRC_UNPROTECTED}, "
            f"{len(only_h)} only in {SRC_THRESHOLD}; using the {len(keys)} matched pairs")
    scenarios = sorted({k[2] for k in keys})
    models = sorted({k[0] for k in keys})
    gcms = sorted({k[1] for k in keys})
    log(f"Matched pairs: {len(keys)} | scenarios: {scenarios}")
    log(f"Models: {models}")
    log(f"GCMs: {gcms}")

    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Pass 1: per-member decadal frequency + footprint ----------------------
    freq = {s: {} for s in scenarios}
    foot = {s: {} for s in scenarios}
    fpf_stack = []
    lats = lons = None
    for k in keys:
        model, gcm, scen = k
        na, exceeded, valid, yrs, fpf = load_pair(files_n[k], files_h[k])
        if lats is None:
            with xr.open_dataset(files_n[k]) as d0:
                lats, lons = d0.lat.values, d0.lon.values
        fpf_stack.append(fpf)
        LAT, LON = len(lats), len(lons)
        fq = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        fp = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        for i, d in enumerate(DECADES):
            a, b = decade_stats(na, exceeded, valid, yrs, d)
            if a is not None:
                fq[i], fp[i] = a, b
        member = f"{model}_{gcm}"
        freq[scen][member] = fq
        foot[scen][member] = fp
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            log(f"  {model:<10} {gcm:<13} {scen}  2020s: exceed={np.nanmean(fq[0]):6.2f}% of yrs"
                f"  footprint={np.nanmean(fp[0]):6.2f}% of cell")
        del na, exceeded, valid

    for s in scenarios:
        if len(freq[s]) != EXPECTED_MEMBERS:
            log(f"  WARNING: scenario {s} has {len(freq[s])} members, expected "
                f"{EXPECTED_MEMBERS}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fpf_layer = np.nanmax(np.stack(fpf_stack, 0), axis=0).astype(np.float32)

    A = cell_area_km2(lats)
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in freq[s]})

    # ---- Shared 2020s baseline (identical across scenarios) --------------------
    def shared(store):
        out = []
        for member in all_members:
            per = [store[s][member][b_idx] for s in scenarios if member in store[s]]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out.append(np.nanmean(np.stack(per, 0), axis=0))
        return np.stack(out, 0).astype(np.float32)

    sh_freq, sh_foot = shared(freq), shared(foot)

    sheet_path = None
    try:
        sheet_path = render_contact_sheet(
            {m: sh_freq[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", LAYER_ID, BASELINE_DECADE,
            units="% of years exceeding the preindustrial 1-in-100 flood",
            note=(f"{len(all_members)} members (6 GHMs x 4 GCMs). Panels show the PRIMARY "
                  f"value (exceedance frequency), not the footprint. All share one "
                  f"hydrodynamic model (CaMa-Flood v3.6.2), so they should agree on river "
                  f"GEOMETRY and differ in magnitude. The ~1.9x level step at 60N is a "
                  f"known source seam (SRTM/HydroSHEDS limit), not a member defect."))
        log(f"  contact sheet: {sheet_path}  <-- REVIEW before trusting the layer")
    except Exception as e:                                          # noqa: BLE001
        log(f"  WARNING: contact sheet failed ({type(e).__name__}: {e})")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sh_f_med, sh_f_sd = np.nanmean(sh_freq, 0), np.nanstd(sh_freq, 0)
        sh_p_med, sh_p_sd = np.nanmean(sh_foot, 0), np.nanstd(sh_foot, 0)
    sh_f_lo = np.clip(sh_f_med - sh_f_sd, 0, None)
    sh_f_hi = np.clip(sh_f_med + sh_f_sd, None, VMAX)
    sh_p_lo = np.clip(sh_p_med - sh_p_sd, 0, None)
    sh_p_hi = np.clip(sh_p_med + sh_p_sd, None, VMAX)
    sh_nmem = np.sum(np.isfinite(sh_freq), 0).astype(np.float32)
    sh_nmem_fp = np.sum(np.isfinite(sh_foot), 0).astype(np.float32)
    mdl_of_all = np.array([m.split("_")[0] for m in all_members])
    sh_nmodel = np.zeros_like(sh_nmem)
    for mdl in np.unique(mdl_of_all):
        sh_nmodel += np.any(np.isfinite(sh_freq[mdl_of_all == mdl]), 0).astype(np.float32)

    land = sh_nmem > 0
    pct, pct_mode, frac_zero = make_pct_fn(sh_f_med[np.isfinite(sh_f_med)].ravel())
    sh_pct = pct(sh_f_med)
    log(f"\nShared 2020s baseline: {sh_freq.shape[0]} members, land cells "
        f"{int(land.sum()):,}, exact-zero fraction={frac_zero:.2%}, percentile={pct_mode}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ever = np.isfinite(sh_p_med)
        log(f"  mean exceedance frequency over cells that exceed at all: "
            f"{np.nanmean(sh_f_med[ever]):.3f}% of years "
            f"(= 1-in-{100/max(np.nanmean(sh_f_med[ever]),1e-9):.1f} years)")
        log(f"  event footprint defined on {int(ever.sum()):,} cells; "
            f"global footprint area = "
            f"{np.nansum(sh_p_med/100*A[:, None]):,.0f} km2")

    # ---- Phase A/B: per scenario ---------------------------------------------
    for s in scenarios:
        members = sorted(freq[s])
        mdl_of = np.array([m.split("_")[0] for m in members])
        fq = np.stack([freq[s][m] for m in members], 0)
        fp = np.stack([foot[s][m] for m in members], 0)
        shp = (len(DECADES), len(lats), len(lons))
        median, lower, upper, percentile = (np.full(shp, np.nan, np.float32) for _ in range(4))
        f_med, f_lo, f_hi = (np.full(shp, np.nan, np.float32) for _ in range(3))
        nmem, nmodel, nmem_fp = (np.zeros(shp, np.float32) for _ in range(3))
        for i, d in enumerate(DECADES):
            lf, lp = fq[:, i], fp[:, i]
            nmem[i] = np.sum(np.isfinite(lf), 0)
            nmem_fp[i] = np.sum(np.isfinite(lp), 0)
            for mdl in np.unique(mdl_of):
                nmodel[i] += np.any(np.isfinite(lf[mdl_of == mdl]), 0)
            if d == BASELINE_DECADE:
                median[i], lower[i], upper[i], percentile[i] = (
                    sh_f_med, sh_f_lo, sh_f_hi, sh_pct)
                f_med[i], f_lo[i], f_hi[i] = sh_p_med, sh_p_lo, sh_p_hi
                nmem[i], nmodel[i], nmem_fp[i] = sh_nmem, sh_nmodel, sh_nmem_fp
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(lf, 0); sdf = np.nanstd(lf, 0)
                f_med[i] = np.nanmean(lp, 0);  sdp = np.nanstd(lp, 0)
            lower[i] = np.clip(median[i] - sdf, 0, None)
            upper[i] = np.clip(median[i] + sdf, None, VMAX)
            f_lo[i] = np.clip(f_med[i] - sdp, 0, None)
            f_hi[i] = np.clip(f_med[i] + sdp, None, VMAX)
            percentile[i] = pct(median[i])
        nmem[nmem == 0] = np.nan
        nmodel[~np.isfinite(nmem)] = np.nan
        nmem_fp[nmem_fp == 0] = np.nan

        trend = np.full(shp, np.nan, np.float32)
        log(f"\n--- {s} ---")
        for i, d in enumerate(DECADES):
            trend[i] = anchored_trend(median, i, b_idx)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ev = np.isfinite(f_med[i])
                mf = np.nanmean(median[i][ev]) if ev.any() else np.nan
                log(f"  {d}s: exceedance={mf:6.3f}% of yrs "
                    f"(1-in-{100/max(mf,1e-9):5.1f} yr)  "
                    f"footprint={np.nansum(f_med[i]/100*A[:, None]):12,.0f} km2  "
                    f"trend={np.nanmean(trend[i]):+.4f} %yrs/dec")

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
                "n_members": (["decade", "lat", "lon"], nmem),
                "n_models": (["decade", "lat", "lon"], nmodel),
                "event_footprint": (["decade", "lat", "lon"], f_med),
                "event_footprint_lower_ci": (["decade", "lat", "lon"], f_lo),
                "event_footprint_upper_ci": (["decade", "lat", "lon"], f_hi),
                "n_members_footprint": (["decade", "lat", "lon"], nmem_fp),
                "floodplain_fraction": (["lat", "lon"], fpf_layer),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": ("Exceedance frequency of the preindustrial 1-in-100-year "
                              "river flood (percent of years)"),
                "units": "%",
                "statistic": "decadal_mean_exceedance_frequency_percent_of_years",
                "primary_value": (
                    "median = EXCEEDANCE FREQUENCY, in % of years in which the "
                    "preindustrial 1-in-100-year flood is exceeded. The event FOOTPRINT "
                    "(% of cell area flooded when it is exceeded) is carried separately as "
                    "`event_footprint` with its own CI. The frequency is primary because "
                    "the footprint is nearly static -- +1.8% from the 2020s to the 2090s "
                    "at rcp85, since topography fixes the extent of a given-magnitude "
                    "flood -- while the frequency more than doubles (+122%). Making the "
                    "footprint primary would produce a confidently flat layer."),
                "derivation": (
                    "Per member and decade, from the paired annual series of the two "
                    "sibling layers: exceedance_frequency = (years with fldfrc_100yr > 0) "
                    "/ (valid years); event_footprint = mean of fldfrc_none over exactly "
                    "those years. A year counts as valid only where BOTH fields are "
                    "present, so numerator and denominator rest on the same year set. "
                    "Then equal-weighted across the 24 GHM x GCM members."),
                "why_not_a_sum": (
                    "This layer exists because `none` + `100yr` does NOT give a "
                    "1-in-100-year flood extent: `100yr` is a SUBSET of `none`, not its "
                    "complement. Measured per cell per year at native 150 arcsec -- "
                    "100yr > none in 0.000% of cells; where 100yr > 0 it EQUALS none in "
                    "84-88% of cells (mean ratio 0.90-0.93); and none + 100yr would exceed "
                    "100% of a cell in 6,867-28,139 cells, which is impossible for an area "
                    "fraction. `100yr` is a filtered copy of `none`, retained only in years "
                    "that overtopped a 1-in-100 defence. The sum is also numerically "
                    "unstable: 17% BELOW the measured footprint in the 2020s, drifting with "
                    "scenario because the terms scale differently (none +5%, 100yr +226%)."),
                "known_issues": (
                    "(1) '1-in-100' IS PREINDUSTRIAL, NOT PRESENT-DAY. The threshold is a "
                    "GEV fit to picontrol (fit_scenario='picontrol', fit_soc='1860soc'), so "
                    "this measures exceedance of the PREINDUSTRIAL 1-in-100 discharge. "
                    "Already in the 2020s the mean frequency over cells that exceed at all "
                    "is ~1-in-12 years, not 1-in-100. (2) THAT ABSOLUTE LEVEL IS BIASED "
                    "UPWARD by a selection effect -- cells that never exceed within a "
                    "decade are excluded from the conditional mean -- so read the "
                    "2020s->2090s RATIO as the robust quantity and the level as "
                    "indicative. (3) `event_footprint` is undefined where a member never "
                    "exceeds in a decade, so its ensemble depth is thinner and more "
                    "variable than the frequency's; use `n_members_footprint`, not "
                    "`n_members`, to qualify it. (4) SHARP ~1.9x LEVEL DISCONTINUITY AT "
                    "60.0 N inherited from the source (CaMa-Flood changes DEM at the "
                    "SRTM/HydroSHEDS coverage limit), present in the native 150 arcsec "
                    "data and identical across all 6 GHMs. It is a LEVEL bias, not a trend "
                    "bias -- trend and change difference the same cells against their own "
                    "baseline, so it largely cancels. See GUARDRAILS §14."),
                "known_latitude_seams": "60.0",
                "normalization": (
                    "none -- one hydrodynamic model (CaMa-Flood v3.6.2) fed by 6 GHMs, so "
                    "inter-member spread is genuine GHM+GCM uncertainty and becomes the CI."),
                "spatial_smoothing": "none (24 members per scenario is thick)",
                "trend_definition": (
                    "baseline-anchored rate on the PRIMARY value: trend[decade] = "
                    "(median[decade] - median[2020s]) / elapsed decades, in percentage "
                    "points of years per decade. Matches the three sibling flood layers so "
                    "all four are comparable. A within-decade slope is not used -- annual "
                    "flooded area swings ~17x between adjacent decades. The 2020s baseline "
                    "has no elapsed change -> 0, identical across scenarios."),
                "trend_units": "% of years decade-1",
                "ci_definition": (
                    "lower/upper_ci = ensemble mean +/- 1 inter-member SD over the 24 "
                    "GHM x GCM members, clamped to [0, 100]. Safe: both the frequency "
                    "(% of years) and the footprint (% of cell area) are bounded shares, "
                    "so min(mean+sd, 100) >= mean. event_footprint_lower_ci / "
                    "event_footprint_upper_ci are the same construction for the companion."),
                "percentile_baseline": (
                    f"{pct_mode}: exceedance frequency ranked against the 2020s "
                    "ensemble-mean land distribution."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "n_models": len(models),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_sens": "2005soc / co2, uniform across all source files",
                "source_dataset": (
                    "ISIMIP2b DerivedOutputData/Zimmer2023 fldfrc at protection levels "
                    "'0' (none) and '100' (100yr), 150 arcsec, coarsened area-preservingly "
                    "to 0.5 deg. CaMa-Flood v3.6.2, GEV fit to picontrol, "
                    "discharge_threshold 0.1mm/d, github.com/swillner/flood-processing, "
                    "doi 10.5281/zenodo.1241051."),
                "derived_from_layers": f"{SRC_UNPROTECTED} + {SRC_THRESHOLD}",
                "sibling_layers": (
                    "riverflood_fldfrc-{none,100yr,flopros}_annual are expected-annual-area "
                    "layers at three protection standards; this fourth layer is the only "
                    "one expressed as event frequency + event footprint. Together they "
                    "separate the two things climate change does to flooding: it barely "
                    "changes the extent of a given flood (+1.8%) but roughly doubles how "
                    "often it happens (+122% at rcp85), which is why `none` saturates at "
                    "+5.1% and `100yr` explodes at +225.5%."),
                "description": (
                    "Exceedance frequency of the preindustrial 1-in-100-year river flood, "
                    "with the undefended event footprint as a companion; TCFD "
                    "6-value-class format, shared 2020s baseline, 24-member (6 GHM x 4 GCM) "
                    "ensemble, rcp26/60/85, no normalization, no spatial smoothing."),
            },
        )
        path = out_dir / f"{VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    if args.no_publish:
        log(f"\n--no-publish: staged only at {stage}")
        return dict(layer_id=LAYER_ID, version=None)

    # Provenance from BOTH sibling raw prefixes.
    raw_entries = []
    for k in keys:
        for src in (files_n[k], files_h[k]):
            p = Path(src)
            with xr.open_dataset(p) as rds:
                ra = rds.attrs
            raw_entries.append({
                "name": p.name,
                "bytes": p.stat().st_size,
                "sha256": storage.sha256_file(p),
                "source_url": ra.get("source_url", ""),
                "source_sha512": ra.get("source_sha512", ""),
                "transform": ra.get("transform", ""),
            })
    missing = [e["name"] for e in raw_entries if not e["source_url"]]
    log(f"  provenance: {len(raw_entries)} raw inputs from 2 sibling prefixes, "
        f"{len(raw_entries)-len(missing)} with source_url + original sha512")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID, stage, on_exists="bump", raw_entries=raw_entries,
        created_by="scripts/process_fldfrc_event100yr.py",
        notes=("1-in-100-year flood EVENT: exceedance frequency (primary) + undefended "
               "event footprint, derived from the paired none/100yr annual series. "
               "24 members/scenario, rcp26/60/85. Built because none+100yr double-counts "
               "(100yr is a subset of none). See WORKFLOW-ISSUES.md 2026-07-29."),
    )
    log("\nGenerating QA/QC report and maps...")
    res = finalize_layer(LAYER_ID, version=version,
                         extra_maps=[sheet_path] if sheet_path else None)
    log(f"\n{LAYER_ID}  {version}  QA="
        f"{(res.get('qa') or {}).get('summary',{}).get('verdict','-')}  "
        f"fit_for_use={res.get('fit_for_use')}")
    return dict(layer_id=LAYER_ID, version=version)


if __name__ == "__main__":
    main()
