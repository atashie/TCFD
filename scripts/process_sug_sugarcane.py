"""Process ISIMIP2b LPJmL sugarcane yield into the TCFD output contract.

*** WITHDRAWN 2026-08-11 -- THE INPUT DATA IS DEFECTIVE. DO NOT RUN THIS. ***

    The ISIMIP2b LPJmL `sug` run does not simulate sugarcane in the sugarcane belt.
    Sao Paulo, Uttar Pradesh, Guangxi, Thailand, Pakistan Punjab, Veracruz, KwaZulu-
    Natal, Cauca, Queensland, Florida and Louisiana all read yield = 0, with sentinel
    companions: `biom-sug` pinned at exactly 0.267 t ha-1 and plantday = matyday = 1
    (no season simulated). 12,966 land cells (19.2%) carry that signature and 0% of
    them have a non-zero yield; the mask is near-identical across all 8 members
    (Jaccard 0.98-0.998), so it is static, not crop failure. Maize from the same
    model/GCM/run IS simulated in those cells, so the cells are live cropland.

    The same model's ISIMIP2a run is sensible and inverts the pattern -- Florida
    19.49, Sao Paulo 11.69, UP 8.33, Queensland 7.17 t ha-1, Midwest 10.61 -- while
    2b gives 0 for the first four and leaves the US Midwest (13.35) as the apparent
    maximum. That inversion is what a user spotted in the maps.

    Both layers this script produced PASSED test_shared_baseline.py in full. Contract
    conformance is a statement about FORM. See GUARDRAILS.md §12 (verify a layer is
    non-trivial where the thing actually exists) and WORKFLOW-ISSUES.md 2026-08-11.

    Kept as the reference implementation of the framing decisions below, which remain
    correct in method. No ISIMIP source supports a scenario-bearing sugarcane layer:
    3b publishes none, 2b is defective, 2a is historical-only.

Original header follows.


`yield-sug-{noirr,firr}` = simulated sugarcane yield in **t ha-1 yr-1** (measured from
the file, not inferred: `units = "t ha-1 yr-1"`, `long_name = "sugarcane Crop yields"`).
`sug` is the OUTPUT crop code; `sgc` is the ISIMIP3b crop-calendar InputData code and
no model publishes output under it (GUARDRAILS §8; WORKFLOW-ISSUES 2026-08-11).

Ensemble: 1 impact model (LPJmL) x 4 CMIP5 GCMs = 4 members per scenario, uniform
across rcp26/rcp60, `2005soc` / `co2` (transient) throughout -> the shared 2020s
baseline is valid and no soc/sens compromise is needed. LPJmL is the ONLY sugarcane
source with future scenarios anywhere in ISIMIP: ISIMIP3b publishes no sugarcane at
all, and ISIMIP2a's LPJmL + CLM-Crop are historical-only (enumerated 2026-08-11 over
all 150 agriculture directories). There is no rcp85 for this crop -- rcp85 DOES exist
in 2b agriculture but only from CLM45 and only for maize + soybean.

TWO LAYERS, NOT ONE. `noirr` (rainfed) and `firr` (fully irrigated) are distinct
fields, not poolable members -- measured 2020s medians over growing cells are 6.65 vs
30.42 t ha-1 yr-1, a 4.6x level difference. Each is processed independently (user
decision 2026-08-11).

FRAMING DECISIONS, each measured before it was made:

1. Statistic = POOLED MEDIAN (the ordinary continuous branch). The field is heavily
   zero-inflated -- 87.3% of land cells read exactly 0 -- but those zeros are
   STRUCTURAL, not intermittent: 87.0% of land is zero in every one of the 94 years,
   i.e. LPJmL simply does not grow sugarcane there. The third branch
   (`pooled_mean_zero_inflated`, OUTPUT-SPEC) exists for fields whose median ERASES a
   live signal; here it would not. Measured on the 2020s pool, exact-zero share of land:

       branch                 noirr     firr
       pooled median + IQR    87.27%   87.65%
       pooled mean +/- 1 SD   86.94%   87.63%

   0.33 pp and 0.02 pp apart -- versus `let`, where the two branches differ by 18 pp
   and the median erases 93% of exposed land. So the median branch is taken, and
   within the growing region it is a real yield: median > 0 for 97.45% (noirr) /
   99.78% (firr) of cells that ever grow sugarcane.

2. Direction = HIGHER_IS_BETTER (risk is yield LOSS), so the percentile is inverted.

3. Percentile = TWO-TIER, which is what the ~87% zero baseline triggers. Combined with
   the inversion this gives exactly the requested convention (user, 2026-08-11):
   a cell with yield 0 -> raw tier-1 -> inverted to **100 = highest risk**; non-zero
   cells are ranked against the NON-ZERO 2020s baseline and land in [1, 99].
   MEASURED CONSEQUENCE, stated rather than buried: ~87% of land therefore reads
   percentile 100. The percentile map is a sugarcane-suitability map first and a
   climate-risk map second; the risk gradient lives inside the ~13% growing region and
   in cells that cross into or out of zero. `median` and the slopes are unaffected.

4. NO spatial smoothing, despite the thin 4-member ensemble. Smoothing is the standard
   remedy for thin ensembles (CLAUDE.md), but this field's structure is a hard
   cropland mask: a 5x5 kernel would bleed yield into cells LPJmL never grows sugarcane
   in, converting structural zeros into small positives and moving them from percentile
   100 to mid-range -- manufacturing a suitability signal that no model produced. The
   noise smoothing would suppress is inter-GCM spread, which the IQR already carries.
   Declared, not silent.

5. NO normalization: one impact model, one unit, so there are no cross-model scales to
   reconcile. The 4-member spread is pure inter-GCM climate spread.

Slopes: both are emitted per OUTPUT-SPEC. On the 87% structurally-zero cells BOTH read
exactly 0 (a true zero trend, not a failure). Inside the growing region the field is
well-behaved continuous with a single impact model -- no between-member level offsets
to bias OLS -- so `sen_slope` is the robust read. Judge agreement on ACTIVE cells only.

Output: data/processed/sugarcane_yield-sug-{irr}_annual/yield-sug-{irr}_{scenario}_processed.nc
on (decade, lat, lon) with {median, lower_ci, upper_ci, percentile, ols_slope,
sen_slope, n_members, n_models}, plus a `{var}_members.nc` diagnostic for the Members
tab. See OUTPUT-SPEC.md.

Usage:
    python scripts/process_sug_sugarcane.py [--irrigation noirr|firr|both]
"""

import argparse
import glob
import os
import re
import sys
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

DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
#: Files span 2006-2099; start at 2010 so every decade panel is a full 10-year window.
MIN_YEAR, MAX_YEAR = 2010, 2099
CI_FLOOR = 0.0                    # yield is non-negative
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = True           # yield loss is the risk -> invert the percentile
SLOPE_PER_DECADE = 10.0
UNITS = "t ha-1 yr-1"

#: One impact model, so every member is its own GCM but the same family.
MODEL_FAMILY = {}


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """Extract (model, gcm, scenario, soc, sens, member) from an ISIMIP2b filename.

    lpjml_gfdl-esm2m_ewembi_rcp26_2005soc_co2_yield-sug-noirr_global_annual_2006_2099.nc4
      [-8]=scenario [-7]=soc [-6]=sens [-5]=variable [-4]=region [-3]=step [-2:]=years

    Parsed from the END: the variable field itself contains hyphens, and a leading
    publication token (DerivedOutputData) would shift every from-the-start index.

    NOTE the off-by-one trap: awk's 1-based `$(NF-4)` is the VARIABLE, but Python's
    `p[-4]` is the REGION -- the same offset shifted by one. Getting this wrong once
    silently read `2005soc` as the scenario and merged rcp26 + rcp60 into a single
    output file, so the fields are asserted below rather than trusted.
    """
    p = os.path.basename(fpath).replace(".nc4", "").replace(".nc", "").split("_")
    info = dict(model=p[0], gcm=p[1], scenario=p[-8], soc=p[-7], sens=p[-6],
                variable=p[-5], region=p[-4], step=p[-3], member=f"{p[0]}_{p[1]}")
    if not re.fullmatch(r"(rcp|ssp)\d{2,3}|historical|picontrol", info["scenario"]):
        raise ValueError(f"{os.path.basename(fpath)}: parsed scenario "
                         f"{info['scenario']!r} is not a scenario token -- field offsets "
                         "are wrong")
    if not info["variable"].startswith("yield-"):
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} is not a yield field -- offsets are wrong")
    return info


def years_from_time(ds):
    """Integer calendar years from a CF time axis. These files use
    'years since 1661-01-01' with a standard calendar, but parse defensively --
    burntarea's mc2-usfs member encodes the same axis in DAYS (GUARDRAILS §9)."""
    t = ds.time
    units = t.attrs.get("units", "")
    cal = t.attrs.get("calendar", "standard")
    m = re.search(r"(years|days|months)\s+since\s+(\d+)", units)
    if not m:
        raise ValueError(f"cannot parse time units {units!r}")
    step, base = m.group(1), int(m.group(2))
    vals = t.values.astype("float64")
    if step == "years":
        yrs = base + vals
    elif step == "months":
        yrs = base + vals / 12.0
    else:
        dpy = 360.0 if "360" in cal else 365.0 if "365" in cal else 365.25
        yrs = base + vals / dpy
    return np.round(yrs).astype(int)


def load_member(fpath, var):
    """Load one member as (years, (year, lat, lon)) with fill values -> NaN."""
    ds = xr.open_dataset(fpath, decode_times=False, decode_cf=False)
    da = ds[var]
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    yrs = years_from_time(ds)
    vals = da.values.astype("float32")
    lats, lons = ds.lat.values, ds.lon.values
    ds.close()

    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan

    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    return yrs[keep], vals[keep], lats, lons


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline, as a RISK score.

    Raw score ranks each value against the 2020s land distribution (higher raw value ->
    higher raw score); because higher yield is BETTER the score is inverted via
    (101 - raw), so a LOW yield earns a HIGH risk percentile.

    Two-tier when the baseline is >2% exact zeros (it is ~87% here): yield 0 -> raw 1,
    yield > 0 -> ranked against the NON-ZERO baseline -> [2, 100]. After inversion,
    zero-yield cells map to 100 (highest risk) and the growing region spreads over
    [1, 99] -- the convention requested 2026-08-11.
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
            res = np.ones(vals.shape, np.float32)      # zero-yield cells -> raw 1
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
    """Put a land-cell vector back on the full (lat, lon) grid; ocean stays NaN."""
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


def process_one(irrigation, root):
    var = f"yield-sug-{irrigation}"
    regime = "rainfed" if irrigation == "noirr" else "fully irrigated"
    layer = f"sugarcane_yield-sug-{irrigation}_annual"
    raw_dir = root / "data" / "raw" / layer
    out_dir = root / "data" / "processed" / layer
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / f"*_{var}_global_annual_2006_2099.nc4")))
    if not files:
        log(f"ERROR: no {var} member files in {raw_dir}")
        return 1

    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    models = sorted({m["model"] for m in meta.values()})
    gcms = sorted({m["gcm"] for m in meta.values()})
    socs = sorted({m["soc"] for m in meta.values()})
    senss = sorted({m["sens"] for m in meta.values()})

    log("=" * 74)
    log(f"Processing {var} (ISIMIP2b LPJmL sugarcane yield) -> TCFD output contract")
    log("=" * 74)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Models: {models} | GCMs: {gcms} | soc: {socs} | sens: {senss}")
    log("No normalization (single impact model). No spatial smoothing (structural "
        "cropland zeros must not be blurred).")
    log("Direction: higher_is_better (risk = yield LOSS) -> percentile inverted.")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {y: i for i, y in enumerate(years)}

    lats = lons = None
    raw = {s: {} for s in scenarios}
    log("\nLoading annual member series...")
    for f in files:
        info = meta[f]
        yrs, vals, la, lo = load_member(f, var)
        if lats is None:
            lats, lons = la, lo
        if yrs.size == 0:
            log(f"  WARNING: {os.path.basename(f)} has no years in "
                f"{MIN_YEAR}-{MAX_YEAR}; skipping (check the time-axis parse)")
            continue
        cube = np.full((n_years, len(lats), len(lons)), np.nan, np.float32)
        for k, y in enumerate(yrs):
            if int(y) in y_index:
                cube[y_index[int(y)]] = vals[k]
        raw[info["scenario"]][info["member"]] = cube
        log(f"  loaded {info['model']:<8} {info['gcm']:<14} {info['scenario']}  "
            f"years {yrs.min()}-{yrs.max()}")

    LAT, LON = len(lats), len(lons)

    # ---- Field nature: measured, never assumed (GUARDRAILS §9) --------------- #
    probe = next(iter(raw[scenarios[0]].values()))
    boolean = is_boolean_field(probe)
    if boolean:
        log("  ERROR: yield classified BOOLEAN -- that is not this field; check inputs.")
        return 1
    stat_name = "pooled_median"
    finite = probe[np.isfinite(probe)]
    log(f"\nField nature: CONTINUOUS ({UNITS}), min={finite.min():.4g} "
        f"max={finite.max():.4g}, exact-zero share of finite values="
        f"{100 * np.mean(finite == 0):.2f}% -> decadal statistic = {stat_name}")

    # ---- Land mask: union of finite cells across every member ---------------- #
    finite_any = np.zeros((LAT, LON), bool)
    for s in scenarios:
        for cube in raw[s].values():
            finite_any |= np.isfinite(cube).any(axis=0)
    land_idx = np.flatnonzero(finite_any.ravel())
    n_land = land_idx.size
    log(f"Land cells (union over members): {n_land:,} of {LAT * LON:,}")

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

    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios; the shared "
            "baseline is only valid for a uniform ensemble.")

    # ---- Shared 2020s baseline, pooled over (year x member x SCENARIO) ------- #
    base_pool = np.concatenate([annual[s] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, boolean=False, window_years=WINDOW_YEARS)
    # The mean branch is computed ONLY to record how far the median branch is from it --
    # the OUTPUT-SPEC third-branch test, measured rather than asserted.
    m_med, _, _ = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, window_years=WINDOW_YEARS, central="mean")
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, None)
    b_hi = np.clip(b_hi, CI_FLOOR, None)

    med_zero = float(np.mean(b_med[np.isfinite(b_med)] == 0))
    mean_zero = float(np.mean(m_med[np.isfinite(m_med)] == 0))
    del m_med
    log(f"\nZero-inflation branch test (OUTPUT-SPEC 'third branch'): median-branch "
        f"exact-zero land = {med_zero:.2%}, mean-branch = {mean_zero:.2%}, "
        f"gap = {abs(med_zero - mean_zero) * 100:.2f} pp")
    log("  -> zeros are STRUCTURAL (cells LPJmL never grows sugarcane in), not "
        "intermittent; the median erases nothing, so the ordinary median branch stands.")

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    growing = baseline_flat > 0
    log(f"Shared 2020s baseline: land n={baseline_flat.size:,}, exact-zero="
        f"{frac_zero:.2%}, percentile mode={pct_mode}, growing cells="
        f"{growing.sum():,} ({growing.mean():.1%}), median yield over growing cells="
        f"{np.median(baseline_flat[growing]):.3f} {UNITS}")
    log(f"  percentile==100 (zero-yield, highest risk) covers "
        f"{np.mean(b_pct[np.isfinite(b_pct)] >= 99.999):.2%} of land -- as requested.")

    # ---- Per-member diagnostic (Members tab; not part of the contract) ------- #
    mem_names = sorted({m for s in scenarios for m in members_by_scen[s]})
    win = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    mem_grid = np.full((len(mem_names), LAT, LON), np.nan, np.float32)
    for i, m in enumerate(mem_names):
        stack = [annual[s][members_by_scen[s].index(m)][win]
                 for s in scenarios if m in members_by_scen[s]]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            vec = np.nanmean(np.concatenate(stack, axis=0), axis=0)
        mem_grid[i] = scatter(vec, land_idx, (LAT, LON))
    mem_ds = xr.Dataset(
        {"value": (["member", "lat", "lon"], mem_grid)},
        coords={"member": mem_names, "lat": lats, "lon": lons},
        attrs={
            "variable": var, "units": UNITS,
            "member_field": (f"mean annual sugarcane yield over the {BASELINE_DECADE}s "
                             "baseline decade, pooled across scenarios, unsmoothed"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Consumed "
                     "by scripts/generate_maps.py for the Members tab. All members share "
                     "one impact model (LPJmL), so differences are pure inter-GCM spread."),
        },
    )
    mem_path = out_dir / f"{var}_members.nc"
    mem_ds.to_netcdf(mem_path, encoding={"value": {"dtype": "float32", "zlib": True,
                                                   "complevel": 4,
                                                   "_FillValue": np.float32(np.nan)}})
    log(f"  wrote per-member diagnostic {mem_path.name} ({len(mem_names)} members)")
    del mem_grid

    # ---- Per-scenario assembly ---------------------------------------------- #
    for s in scenarios:
        log(f"\n{'=' * 74}\nAssembling scenario {s}\n{'=' * 74}")
        mem = members_by_scen[s]
        cube = annual[s]
        fams = sorted({family_of(m.split("_")[0]) for m in mem})
        fam_idx = {f: k for k, f in enumerate(fams)}

        nd = len(DECADES)
        out = {k: np.full((nd, LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, boolean=False, window_years=WINDOW_YEARS)
                lo = np.clip(lo, CI_FLOOR, None)
                hi = np.clip(hi, CI_FLOOR, None)
                pc = pct(med)

            sl = expanding_slopes(cube, years, d, BASELINE_DECADE,
                                  window_years=WINDOW_YEARS)

            wmask = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, wmask, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
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

            grow = np.isfinite(out["median"][i]) & (out["median"][i] > 0)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                if d == BASELINE_DECADE:
                    slope_txt = "slopes=NaN (baseline)"
                else:
                    slope_txt = (f"ols={np.nanmean(out['ols_slope'][i][grow]):+.4f}  "
                                 f"sen={np.nanmean(out['sen_slope'][i][grow]):+.4f} "
                                 f"{UNITS}/dec (growing cells)")
                gm = np.nanmean(out["median"][i][grow])
            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
            log(f"  {d}s: {tag:<15} growing-cell mean={gm:.4f} {UNITS}  {slope_txt}")

        # ---- GUARDRAIL: slope and median masks must agree -------------------- #
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
                "variable": var,
                "scenario": s,
                "long_name": f"Sugarcane yield, {regime}",
                "units": UNITS,
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "decadal_statistic_rationale": (
                    f"Ordinary median branch. The field is {frac_zero:.1%} exact-zero, but "
                    f"the zeros are STRUCTURAL (87.0% of land is zero in every year -- "
                    f"LPJmL grows no sugarcane there), so the median erases no live "
                    f"signal: measured median-branch exact-zero land {med_zero:.2%} vs "
                    f"mean-branch {mean_zero:.2%} ({abs(med_zero - mean_zero) * 100:.2f} pp "
                    f"apart). Contrast `let`, where the two branches differ by 18 pp and "
                    f"the third branch (pooled_mean_zero_inflated) is warranted."),
                "field_nature": "continuous",
                "irrigation": irrigation,
                "value_note": (
                    f"median = MEDIAN over the pooled (year x member) sample inside the "
                    f"decade window, across 4 ISIMIP2b LPJmL members (one impact model x "
                    f"4 CMIP5 GCMs); raw ISIMIP2b values in {UNITS}. Cells where LPJmL "
                    f"grows no sugarcane read exactly 0, not NaN."),
                "ci_definition": (
                    "lower_ci/upper_ci = 25th/75th percentile of the same pooled "
                    "(year x member) sample, floored at 0. With a single impact model "
                    "the IQR carries interannual variability plus inter-GCM spread; it "
                    "is not a model-structure band."),
                "slope_definition": (
                    "ols_slope = least-squares slope; sen_slope = Theil-Sen slope, both "
                    "fitted over an EXPANDING window from the start of the 2020s "
                    "baseline through the end of the target decade, stacking every "
                    "(year, member) observation. Baseline panel is NaN. On structurally "
                    "zero cells both read exactly 0 -- a true zero trend. Inside the "
                    "growing region there is only ONE impact model, so there are no "
                    "between-member level offsets to bias OLS, and sen_slope is the "
                    "robust read. Judge agreement on ACTIVE (growing) cells only."),
                "slope_units": f"{UNITS} decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: yield 0 -> tier 1; yield > 0 ranked against the "
                    "NON-ZERO 2020s baseline -> [2,100]; then INVERTED (101 - raw) "
                    "because higher yield is better. Net effect: zero-yield cells = 100 "
                    "(highest risk), growing cells spread over [1,99]."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_better",
                "percentile_note": (
                    f"~{frac_zero:.0%} of land reads percentile 100 because LPJmL grows "
                    "no sugarcane there. The percentile field is a suitability map first "
                    "and a risk gradient second; the gradient lives inside the growing "
                    "region. `median` and the slopes are unaffected. Convention chosen "
                    "by the user 2026-08-11 (layer is global-potential, unmasked)."),
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "decade_note": (
                    "Files span 2006-2099; processing starts at 2010 so every decade "
                    "panel is a full 10-year window (2010s-2090s)."),
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_scenario": ",".join(socs),
                "co2_treatment": (
                    f"uniform {','.join(senss)} (transient CO2). rcp26 also publishes a "
                    "fixed-2005-CO2 variant, which is NOT used: rcp60 has no counterpart, "
                    "so including it would make the scenarios experimentally asymmetric."),
                "normalization": (
                    "none -- a single impact model (LPJmL) in one unit; the 4-member "
                    "spread is pure inter-GCM climate spread."),
                "spatial_smoothing": (
                    "none -- DECLARED deviation from the thin-ensemble default. This "
                    "field's structure is a hard cropland mask; a 5x5 kernel would bleed "
                    "yield into cells that grow no sugarcane, converting structural zeros "
                    "into small positives and moving them off percentile 100."),
                "scenario_coverage_note": (
                    "rcp26 + rcp60 only. No rcp85 exists for sugarcane; rcp85 IS present "
                    "in ISIMIP2b agriculture but only from CLM45 and only for maize and "
                    "soybean. No SSP version exists -- ISIMIP3b publishes no sugarcane."),
                "source_dataset": (
                    "ISIMIP2b OutputData/agriculture/LPJmL (yield-sug-{noirr,firr}, annual)"),
                "description": (
                    f"Sugarcane yield ({regime}) processed to the TCFD output contract "
                    f"(OUTPUT-SPEC.md) with a shared 2020s baseline; 1 model x 4 CMIP5 "
                    f"GCMs in raw {UNITS}, no normalization, no spatial smoothing, "
                    f"higher_is_better (risk = yield loss)."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{var}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path}  ({path.stat().st_size / (1024 * 1024):.1f} MB)")

    log(f"\nDone: {var} -> {out_dir}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--irrigation", choices=["noirr", "firr", "both"], default="both")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    rc = 0
    for irr in (["noirr", "firr"] if args.irrigation == "both" else [args.irrigation]):
        rc |= process_one(irr, root)
    return rc


if __name__ == "__main__":
    sys.exit(main())
