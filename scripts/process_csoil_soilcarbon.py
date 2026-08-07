"""Process csoil-total (soil organic carbon stock) into the TCFD output contract.

REFERENCE IMPLEMENTATION for OUTPUT-SPEC.md. Other processors should be ported onto
`utils/decadal_stats.py` following the structure here.

csoil-total = the total soil organic carbon pool of each land grid cell, reported
by ISIMIP3b biomes (vegetation) models in kg C m-2. This is the direct subsurface
carbon-STORAGE signal (distinct from the vegetation pools cveg/croot/cvegbg and
from the net-sink FLUX nbp; and distinct from the Lange 2020 exposure family).

Ensemble: 3 ISIMIP3b models x their CMIP6 GCMs x {ssp126, ssp370, ssp585} =
12 members per scenario. The three models were value-checked (2026-07-25), all in
the SAME unit (kg C m-2) with COMPARABLE magnitudes (medians within ~2x):

    model             GCMs  CO2 run     median  mean   max   %zero  long_name
    classic            2    default*    5.76    6.67   70.3  4.48   "Carbon Mass in Soil Pool"
    jules-es-vn6p3     5    2015co2**   10.28   10.69  67.9  0.00   "Carbon Mass in Soil Pool"
    mc2-usfs-r87g5c1   5    default*    7.67    8.42   37.8  0.00   "Carbon Mass in Soil Pool"

   * default        = transient CO2 (SSP-consistent, includes CO2 fertilization).
  ** 2015co2 (FIXED CO2) is the ONLY run JULES-ES-VN6P3 publishes for csoil-total,
     so the ensemble MIXES CO2 treatments -- JULES's soil-carbon TREND is muted
     (no fertilization signal) relative to the two transient models. Land use is
     held fixed for all (2015soc-from-histsoc for classic/jules; nat for mc2-usfs,
     a natural-vegetation biome model with no land-use forcing). User decision
     2026-07-25: retain all 12 members (accept the mixed CO2) rather than drop
     JULES and thin the ensemble to 2 structurally-diverse models.

All three share the SAME unit and comparable magnitudes, so -- like burntarea and
unlike the water-index TWS -- NO normalization is applied: members are equal-
weighted in raw kg C m-2 and the inter-member spread becomes the CI ("model
democracy"; note CLASSIC contributes only 2 GCMs vs 5 each for the others).

NOTE the measured consequence for the slopes: this ensemble's between-member level
offset is ~68.7x its interannual SD, which is the regime where `ols_slope` picks up
level offsets as trend if member coverage is uneven. `sen_slope` is the robust read
here; see OUTPUT-SPEC.md.

Direction: higher = BETTER (more stored carbon), so the risk is LOSS -- the
percentile is INVERTED (low stock / decline -> high risk percentile) and change maps
read red = soil-carbon decline. Thick 12-member ensemble => NO spatial smoothing.
Shared 2020s baseline across ssp126/370/585.

Time axis: all members use "days since 1601" but different calendars (365_day for
classic/mc2-usfs; proleptic_gregorian for jules) -> calendar-aware year parse.
ISIMIP3b csoil starts in 2015, so there is no full 2010s decade; the layer begins
at the 2020s baseline (decades 2020s-2090s).

Output: csoil_{scenario}_processed.nc on (decade, lat, lon), units kg m-2, with
{median, lower_ci, upper_ci, percentile, ols_slope, sen_slope, n_members, n_models}.
See OUTPUT-SPEC.md for the definition of each.
"""

import os
import re
import sys
import glob
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.utils.decadal_stats import (  # noqa: E402
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

VAR = "csoil-total"
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2015, 2099   # ISIMIP3b csoil covers 2015-2100 (2100 unused)
CI_FLOOR = 0.0                    # soil carbon is non-negative; no upper cap
TWO_TIER_ZERO_THRESHOLD = 0.02    # use two-tier percentile if >2% of baseline is exact 0
HIGHER_IS_BETTER = True           # more stored carbon = better; risk = loss -> invert pct

#: Slopes are fitted per YEAR and rescaled once, here, to the declared decade unit.
#: Do NOT also fit against a decade index -- that double-counts and inflates 10x.
SLOPE_PER_DECADE = 10.0

#: Model families: members that are two configurations of one model must not each get
#: a vote (GUARDRAILS). csoil's three models are structurally distinct, so each is its
#: own family; the map is explicit so a future orchidee/orchidee-dgvm pair is handled.
MODEL_FAMILY = {}


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """Extract (model, gcm, scenario, member) from a standard ISIMIP3b filename.

    e.g. classic_gfdl-esm4_w5e5_ssp585_2015soc-from-histsoc_default_csoil-total_global_annual_2015_2100.nc
         [0]=model  [1]=gcm  [2]=bias  [3]=scenario  [4]=soc  [5]=sens  [6]=variable ...
    (model/gcm names contain hyphens but no underscores, so the split is safe.)
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], member=f"{p[0]}_{p[1]}")


def years_from_time(ds):
    """Parse integer calendar years from a CF time axis (xarray cannot decode
    these non-standard-calendar files). All members use 'days since 1601' but
    disagree on the calendar:
      classic/mc2-usfs -> 365_day             (divide by 365)
      jules-es-vn6p3   -> proleptic_gregorian (divide by 365.25)
    so we handle years / days / months and use the calendar's days-per-year."""
    t = ds.time
    units = t.attrs.get("units", "")
    cal = t.attrs.get("calendar", "365_day")
    m = re.search(r"(years|days|months)\s+since\s+(\d+)", units)
    if not m:
        raise ValueError(f"cannot parse time units {units!r}")
    step, base = m.group(1), int(m.group(2))
    vals = t.values.astype("float64")
    if step == "years":
        yrs = base + vals
    elif step == "months":
        yrs = base + vals / 12.0
    else:  # days
        dpy = 360.0 if "360" in cal else 365.0 if "365" in cal else 365.25
        yrs = base + vals / dpy
    return np.round(yrs).astype(int)


def load_member(fpath):
    """Load one member file as (year, lat, lon), restricted to MIN..MAX_YEAR."""
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    return da


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the 2020s baseline, returned as a RISK score.

    The raw score ranks each value against the 2020s land distribution (higher raw
    value -> higher raw score). Because higher soil carbon is BETTER, the score is
    INVERTED to a risk percentile via (101 - raw) -> a cell with a LOW stock (or
    that declines below the baseline distribution) gets a HIGH risk percentile.

    Tier is chosen empirically from the baseline's exact-zero fraction:

    * SINGLE-tier (expected here): rank against the FULL 2020s land distribution
      -> [1, 100]. Correct when zeros are negligible.

    * TWO-tier (zero-inflated): value 0 -> lowest stock; value > 0 -> ranked
      against the NON-ZERO baseline. Used only if the baseline is materially
      zero-inflated (>2% exact zeros). After inversion, exact-zero (barren)
      cells map to the highest risk.
    """
    frac_zero = float(np.mean(baseline_flat == 0.0))
    two_tier = frac_zero >= TWO_TIER_ZERO_THRESHOLD

    def _invert(sc):
        return np.clip(101.0 - sc, 1.0, 100.0) if higher_is_better else sc

    if two_tier:
        nz_sort = np.sort(baseline_flat[baseline_flat > 0])
        n_nz = len(nz_sort)

        def pct(arr):
            flat = arr.ravel()
            out = np.full(flat.shape, np.nan, np.float32)
            fin = np.isfinite(flat)
            vals = flat[fin]
            res = np.ones(vals.shape, np.float32)  # zero-stock cells -> raw 1
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="left") / n_nz
                res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
            out[fin] = _invert(res)
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
                out[fin] = _invert(np.clip(100.0 * frac, 1.0, 100.0)).astype(np.float32)
            return out.reshape(arr.shape)

    return pct, ("two_tier" if two_tier else "single_tier"), frac_zero


def scatter(flat_land, land_idx, shape):
    """Put a land-cell vector back onto the full (lat, lon) grid; ocean stays NaN."""
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


def main():
    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / "soilcarbon_csoil_annual"
    out_dir = root / "data" / "processed" / "soilcarbon_csoil_annual"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / "*_csoil-total_global_annual_2015_2100.nc")))
    if not files:
        log(f"ERROR: no csoil-total member files found in {raw_dir}")
        return
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    log("=" * 66)
    log("Processing csoil-total (soil organic carbon) -> TCFD output contract")
    log("=" * 66)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models} | GCMs: {gcms}")
    log("No normalization (raw kg C m-2, equal-weight model democracy); no spatial smoothing.")
    log("Direction: higher_is_better (risk = soil-carbon LOSS) -> percentile inverted.")

    da0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {y: i for i, y in enumerate(years)}

    # ---- Pass 1: load every member's ANNUAL series ---------------------------
    # The expanding-window slopes need annual values, so we cannot collapse to
    # decadal means on load the way the pre-spec version did.
    log("\nLoading annual member series...")
    raw = {s: {} for s in scenarios}        # raw[scen][member] = (year, lat, lon)
    for f in files:
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
        if da.year.size == 0:
            log(f"  WARNING: {os.path.basename(f)} has no years in "
                f"{MIN_YEAR}-{MAX_YEAR}; skipping (check time-axis parse)")
            continue
        cube = np.full((n_years, LAT, LON), np.nan, np.float32)
        for k, y in enumerate(da.year.values):
            if y in y_index:
                cube[y_index[int(y)]] = da.values[k]
        raw[s][member] = cube
        log(f"  loaded {info['model']:<18} {info['gcm']:<13} {s}")
        del da

    # ---- Field nature: measured, never assumed (GUARDRAILS §9) ---------------
    probe = next(iter(raw[scenarios[0]].values()))
    boolean = is_boolean_field(probe)
    stat_name = "pooled_mean_boolean" if boolean else "pooled_median"
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'} "
        f"-> decadal statistic = {stat_name}")
    if boolean:
        log("  WARNING: csoil is a continuous stock; a boolean classification here "
            "means the input is not what this processor expects. Check the members.")

    # ---- Land mask: union of finite cells across every member ----------------
    finite_any = np.zeros((LAT, LON), bool)
    for s in scenarios:
        for cube in raw[s].values():
            finite_any |= np.isfinite(cube).any(axis=0)
    land_idx = np.flatnonzero(finite_any.ravel())
    n_land = land_idx.size
    log(f"Land cells (union over members): {n_land:,} of {LAT * LON:,}")

    # Repack to (member, year, n_land) -- ~273 MB/scenario at 12 members, vs 1.06 GB
    # on the full grid. Slope chunking already operates on a flat cell axis.
    annual = {}
    members_by_scen = {}
    for s in scenarios:
        mem = sorted(raw[s])
        members_by_scen[s] = mem
        arr = np.full((len(mem), n_years, n_land), np.nan, np.float32)
        for i, m in enumerate(mem):
            arr[i] = raw[s][m].reshape(n_years, -1)[:, land_idx]
        annual[s] = arr
        del raw[s]
    del raw

    # ---- Shared 2020s baseline ----------------------------------------------
    # Pooled over every (year, member, SCENARIO) observation in the 2020s window, so
    # the baseline panel is bit-identical in all three files by construction.
    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios; the shared "
            "baseline is only valid for a uniform ensemble. Declaring "
            "members_by_scenario so QA groups by composition.")
    base_pool = np.concatenate([annual[s] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, boolean=boolean, window_years=WINDOW_YEARS)
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, None)
    b_hi = np.clip(b_hi, CI_FLOOR, None)

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"\nShared 2020s baseline: land n={baseline_flat.size:,}, "
        f"exact-zero fraction={frac_zero:.2%}, percentile mode={pct_mode}, "
        f"global-mean csoil={np.nanmean(b_med):.4f} kg m-2")

    # ---- Per-scenario assembly ----------------------------------------------
    b_idx = DECADES.index(BASELINE_DECADE)
    for s in scenarios:
        log(f"\n{'=' * 66}\nAssembling scenario {s}\n{'=' * 66}")
        mem = members_by_scen[s]
        cube = annual[s]                                   # (member, year, n_land)
        fams = sorted({family_of(m.split("_")[0]) for m in mem})

        nd = len(DECADES)
        out = {k: np.full((nd, LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, boolean=boolean, window_years=WINDOW_YEARS)
                lo = np.clip(lo, CI_FLOOR, None)
                hi = np.clip(hi, CI_FLOOR, None)
                pc = pct(med)

            sl = expanding_slopes(cube, years, d, BASELINE_DECADE,
                                  window_years=WINDOW_YEARS)

            # Per-cell ensemble depth, from the decade window's own coverage.
            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)      # (member, cell)
            n_mem = present.sum(axis=0).astype(np.float32)
            fam_idx = {f: k for k, f in enumerate(fams)}
            fam_present = np.zeros((len(fams), n_land), bool)
            for mi, m in enumerate(mem):
                fam_present[fam_idx[family_of(m.split("_")[0])]] |= present[mi]
            n_mod = fam_present.sum(axis=0).astype(np.float32)
            n_mem[n_mem == 0] = np.nan
            n_mod[np.isnan(n_mem)] = np.nan

            for key, vec in (("median", med), ("lower_ci", lo), ("upper_ci", hi),
                             ("percentile", pc),
                             ("ols_slope", sl.ols_slope * SLOPE_PER_DECADE),
                             ("sen_slope", sl.sen_slope * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
            with warnings.catch_warnings():
                # The baseline panel's slopes are all-NaN by design; nanmean says so.
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                slope_txt = (f"ols={np.nanmean(out['ols_slope'][i]):+.4f}  "
                             f"sen={np.nanmean(out['sen_slope'][i]):+.4f} kg m-2/dec"
                             if d != BASELINE_DECADE else "slopes=NaN (baseline)")
            log(f"  {d}s: {tag:<15} global-mean={np.nanmean(out['median'][i]):.4f}  "
                f"{slope_txt}")

        # ---- GUARDRAIL: slope and median masks must agree ---------------------
        # A bare np.zeros() for the baseline panel would make the whole OCEAN a finite
        # zero, and the QA report does not catch it (it only checks that FINITE
        # baseline trends equal zero, never that the masks agree).
        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), "baseline ols must be NaN"
                assert np.all(np.isnan(out["sen_slope"][i])), "baseline sen must be NaN"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), (
                    f"{k} finite where median is NaN at {d}s ({extra.sum()} cells) "
                    "-- ocean leak")

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": "Soil organic carbon stock (total soil carbon pool)",
                "units": "kg m-2",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "field_nature": "boolean_01" if boolean else "continuous",
                "value_note": (
                    "median = MEDIAN over the pooled (year x member) sample inside the "
                    "decade window, across the 12 ISIMIP3b biomes members "
                    "(classic, jules-es-vn6p3, mc2-usfs x their CMIP6 GCMs); raw "
                    "ISIMIP3b values in kg C m-2."),
                "ci_definition": (
                    "lower_ci/upper_ci = 25th/75th percentile (interquartile range) of "
                    "the same pooled (year x member) sample, floored at 0. The IQR "
                    "therefore carries BOTH interannual variability and inter-model "
                    "spread; it is not a pure model-spread band."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope. Both "
                    "are fitted over an EXPANDING window from the start of the 2020s "
                    "baseline through the end of the target decade, stacking every "
                    "(year, member) observation as an independent point. The baseline "
                    "panel is NaN (no elapsed period). The two estimators fail in "
                    "OPPOSITE regimes -- sen collapses to exactly 0 on zero-inflated "
                    "fields, ols absorbs member level offsets as trend when member "
                    "coverage is uneven -- so disagreement between them is the signal "
                    "that a cell's trend is not robust. This ensemble's between-member "
                    "level offset is ~68.7x its interannual SD, so sen_slope is the "
                    "robust read here."),
                "slope_units": "kg m-2 decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: raw score ranks each cell against the 2020s "
                    "ensemble land distribution, then INVERTED to a risk percentile "
                    "(101 - raw) because higher soil carbon is better -> low stock / "
                    "decline = high risk percentile."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_better",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "decade_note": (
                    "ISIMIP3b csoil covers 2015-2100, so there is no full 2010s "
                    "decade; the layer begins at the 2020s baseline (2020s-2090s)."),
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "normalization": (
                    "none -- all 3 models report the SAME unit (kg C m-2) with "
                    "COMPARABLE magnitudes, so they are equal-weighted in raw "
                    "kg C m-2 (model-democracy decision). Note CLASSIC contributes "
                    "only 2 GCMs vs 5 each for jules/mc2-usfs, so it is underweighted "
                    "in the flat member pool."),
                "co2_treatment": (
                    "MIXED: classic & mc2-usfs use transient CO2 (default run); "
                    "jules-es-vn6p3 publishes ONLY the fixed-2015-CO2 run for "
                    "csoil-total, so its soil-carbon trend is muted. Retained to keep "
                    "12 members (user decision 2026-07-25)."),
                "spatial_smoothing": "none (12-member ensemble is thick enough)",
                "source_dataset": "ISIMIP3b OutputData/biomes (csoil-total, annual)",
                "description": (
                    "Soil organic carbon stock processed to the TCFD output contract "
                    "(OUTPUT-SPEC.md) with a shared 2020s baseline; 3-model x "
                    "CMIP6-GCM ensemble in raw kg C m-2, no normalization, no spatial "
                    "smoothing, higher_is_better (risk = loss)."),
            },
        )

        # Compression: these grids are ~7x smaller zlib-compressed and nothing
        # downstream changes. The pre-spec processors wrote them uncompressed.
        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"csoil_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        size_mb = path.stat().st_size / (1024 * 1024)
        log(f"  saved {path}  ({size_mb:.1f} MB)")

    log("\nDone.")


if __name__ == "__main__":
    main()
