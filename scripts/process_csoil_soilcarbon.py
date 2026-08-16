"""Process csoil-total (soil organic carbon stock) into the TCFD output contract.

REFERENCE IMPLEMENTATION for OUTPUT-SPEC.md. Other processors should be ported onto
`utils/decadal_stats.py` following the structure here.

csoil-total = the total soil organic carbon pool of each land grid cell, reported
by ISIMIP3b biomes (vegetation) models in kg C m-2. This is the direct subsurface
carbon-STORAGE signal (distinct from the vegetation pools cveg/croot/cvegbg and
from the net-sink FLUX nbp; and distinct from the Lange 2020 exposure family).

REBUILT 2026-08-14. The 2026-07-25 build is superseded and its output is gone.

Ensemble: 4 ISIMIP3b models x their CMIP6 GCMs x {ssp126, ssp370, ssp585} =
17 members per scenario, 51 files, 888.7 MiB, every one verified against the
publisher's sha512 sidecar. Ingest: scripts/download_csoil_isimip3b.py.

THE 12-MEMBER ENSEMBLE WAS WRONG BY OMISSION. LPJmL5-7-10-fire publishes
csoil-total annually for all three SSPs and was missing from the catalog's model
list, which read as a measured fact for three weeks (+5 members, +1 structurally
independent model). The sector trap was NOT the cause -- no sector holds a csoil
model `biomes` lacks -- the omission was inside the sector that had been walked,
which no completeness argument about sectors would have caught. Only projecting
the variable field over a full listing did. GUARDRAILS 11.

Value-checked 2026-08-14 by scripts/check_csoil_nature.py over all 51 files. All
four models report the SAME unit and long_name; the 2020s distributions are:

    model             GCMs  CO2 run     p50     p95    max    %zero  finite%
    classic            2    default*     5.70   17.90   69.5   4.50   25.73
    mc2-usfs-r87g5c1   5    default*     7.55   16.51   35.2   0.00   22.73
    jules-es-vn6p3     5    2015co2**   10.45   27.38   67.5   0.00   26.10
    lpjml5-7-10-fire   5    default*    16.64   91.09  179.1   6.19   26.01

   * default  = transient CO2 (SSP-consistent, includes CO2 fertilization).
  ** 2015co2 (FIXED CO2) is the ONLY run JULES-ES-VN6P3 publishes for csoil-total,
     so the ensemble MIXES CO2 treatments. Adding LPJmL improved the balance from
     7/12 transient to 12/17.

DO NOT SAY THE FIXED-CO2 MEMBER'S TREND IS "MUTED" -- IT IS THE STRONGEST LOSER.
The 2026-07-25 build asserted that JULES's soil-carbon trend was muted "because it
carries no fertilization signal". That is the expectation, and it is measurably
backwards. Measured 2026-08-15 on ssp585, 2090s minus 2020s over each member's own
footprint:

    model             sens      2020s mean    change    % change
    lpjml5-7-10-fire  default        28.07    -0.772      -2.75%
    jules-es-vn6p3    2015co2        10.74    -0.470      -4.37%   <-- fixed CO2
    mc2-usfs-r87g5c1  default         8.33    -0.004      -0.05%
    classic           default         6.70    +0.053      +0.79%

The fixed-CO2 member shows the LARGEST RELATIVE LOSS of the four. The mechanism is
real -- no fertilization means less litter input, so less soil carbon -- but it
makes the decline STRONGER, not weaker. The direction a CO2 treatment biases a
trend is a measurement, never a deduction from the mechanism's name (GUARDRAILS 9,
and the `/isimip-process-visualize` skill records this exact layer as the case that
established the rule).

Land use is held fixed in every member (2015soc-from-histsoc for classic/jules/
lpjml; nat for mc2-usfs, a natural-vegetation biome model with no land-use
forcing). No soc/sens combination is shared by all four models -- demanding
uniform sens drops JULES, demanding uniform soc drops MC2 as well, either way
back to 12 members and 3 models. The heterogeneity is RETAINED and DECLARED.

NORMALIZATION IS STILL DECLINED, BUT THE 3-MODEL EVIDENCE FOR IT DID NOT SURVIVE.
The old rationale was "comparable magnitudes, medians within ~2x". With LPJmL the
median spread is 2.92x and LPJmL carries the upper tail alone (p95 91.1 against
16.5-27.4). It is declined on different grounds: the four models form a GRADIENT
rather than clusters, their disagreement about the size of the stock IS the
structural uncertainty the CI exists to carry, and rescaling would destroy units
that can be checked against measured soil-carbon inventories.

WHICH DECADAL BRANCH, MEASURED RATHER THAN DEFAULTED. `pooled_median` is retained
because the pooled sample is unimodal: the largest gap between adjacent model
medians is well under the lower model's own IQR, so the median summarises the pool
instead of selecting a cluster. Contrast `permafrost-3b`, which needed
OUTPUT-SPEC's fourth branch with seven members at ~0.04 and five at ~0.95. The
measurement is re-run at build time and written into
`decadal_statistic_rationale`, so it cannot silently go stale.

GUARDRAILS 12 -- REFERENCE SITES, 2020s mean kg C m-2, measured on the raw files.
A contract PASS is not a sanity check; the sugarcane layers passed every check and
were meaningless. Soil carbon has a strong published prior -- organic soils and
chernozems are the global maxima, hot deserts the minima -- which is what makes
this able to fail:

    site                          classic  jules  lpjml   mc2
    Flow Country blanket bog, GB    51.62  11.15  31.10  20.68   organic, expect high
    Hudson Bay Lowlands, CA         17.16  19.04  56.71  13.84   organic, expect high
    W Siberian Lowland peat, RU     18.12   9.76  52.43   6.32   organic, expect high
    Ukraine chernozem               13.99  21.47  12.18   5.48   chernozem, expect high
    Iowa mollisol, US                8.52  21.47  11.48   8.24   mollisol, expect high
    Amazon terra firme, BR           7.20   5.19  18.12  14.23   mid
    Tamanrasset, Sahara, DZ          0.00   0.00   0.00   1.79   desert, expect ~0
    Atacama, CL                      0.00   0.00   0.03   2.12   desert, expect ~0
    Taklamakan, CN                   0.00   0.03   0.05   0.91   desert, expect ~0
    Pacific ocean (5N,150W)           NaN    NaN    NaN    NaN   control, expect NaN

No model puts a desert at or above an organic soil, and ocean is NaN in all four,
so `isfinite` IS a usable footprint here -- unlike `cropfailure`, whose publisher
zero-fills the globe. Finite share is 22.7-26.1% of the grid against ~27-29% land.

THE MODELS AGREE ON THE ORDERING AND DISAGREE ON THE MAGNITUDE, most sharply in
high-latitude organic soils, and the disagreement is directional rather than
noisy. At the North Slope of Alaska lpjml reads 65.60 against 3.67 / 4.94 / 5.69;
by the 2090s lpjml alone projects large northern peat losses (W Siberia -5.88,
North Slope -4.31, Hudson Bay -4.27, Sodankyla -4.20 kg m-2) where the other three
are flat or slightly positive. The ensemble central value therefore UNDERSTATES a
permafrost-carbon risk that one model of four considers large. That is in
`interpretation_caveat`, which is on the delivery allowlist, so it reaches the
customer instead of dying in this file.

WHICH SLOPE TO READ IS A MEASUREMENT, NOT AN INFERENCE. The 3-model build argued
`sen_slope` from a between-member level offset of ~68.7x the interannual SD -- the
regime where `ols_slope` absorbs member level offsets as trend when coverage is
uneven -- and that reasoning is still sound and still cited in OUTPUT-SPEC.md. It
is nevertheless only a prediction: the standing rule is to measure on ACTIVE cells
with `generate_customer_delivery.py --measure-slopes` and record the result in
DATASET-ATTRIBUTES.md, because `wildfire` looked safe for Sen on its zero fraction
and was not, and `permafrost-3b` had the second-lowest Sen-zero share of any layer
and Sen was still wrong. Do not register this layer's `recommended_slope` from the
paragraph above; register it from the measurement.

THE MASK IS TIME-VARYING, BY ONE CELL. Found on the 2026-08-15 rebuild, when the
mask-agreement guardrail refused to write ssp585: at the 2060s panel one cell has
observations early in the EXPANDING slope window and none inside the decade
itself, so it carried a finite slope against a NaN median -- a trend over a decade
the subject was absent from. The slopes are now masked to each decade's own median
mask, as OUTPUT-SPEC requires of any layer in this regime. One cell is trivial
next to `npp-tempnle`'s 53 rising to 374, and it is fixed rather than tolerated
because they are the same defect: on that layer dropping 53 of 25,821 cells moved
the mean slope from -1.89 to +0.64.

SMOOTHING IS DECLINED ON A SPLIT-HALF TEST, not on ensemble depth. Depth alone does
not settle it -- `cropfailure-3b` had 400 draws per cell-decade and roughness close
to raw `let`. The test runs at build time (halves stratified by model, alternating
GCMs within each model, so the halves have the same model composition) and its
measured numbers are written into the `spatial_smoothing` attribute.

Direction: higher = BETTER (more stored carbon), so the risk is LOSS -- the
percentile is INVERTED (low stock / decline -> high risk percentile) and change maps
read red = soil-carbon decline. Shared 2020s baseline across ssp126/370/585.

Time axis: all members use "days since 1601" but different calendars (365_day for
lpjml/mc2-usfs, proleptic_gregorian for jules, and CLASSIC publishes BOTH across its
two GCMs) -> calendar-aware year parse. ISIMIP3b csoil starts in 2015, so there is no
full 2010s decade; the layer begins at the 2020s baseline (decades 2020s-2090s).

Output: csoil_{scenario}_processed.nc on (decade, lat, lon), units kg m-2, with
{median, lower_ci, upper_ci, percentile, ols_slope, sen_slope, n_members, n_models}.
See OUTPUT-SPEC.md for the definition of each.
"""

import os
import re
import sys
import glob
import argparse
import warnings
import multiprocessing as mp
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

#: Theil-Sen is quadratic in stacked (year x member) points and is the only expensive term.
#: At 17 members x 85 years the final decade stacks 1,445 observations = ~1.04M pairs per
#: cell, roughly double the 12-member build this processor was written for, which put a
#: single-core run near seven hours. The work is embarrassingly parallel over cells, so it
#: is split across processes exactly as `process_thawdepth_permafrost.py` does.
#: Pair SUBSAMPLING is deliberately NOT used: at csoil's ~68x level-offset ratio a 100k cap
#: costs ~15% of the slope (OUTPUT-SPEC Performance).
SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

#: Set per worker by the fork, never pickled -- the cube is hundreds of MB.
_CUBE = None

#: Model families: members that are two configurations of one model must not each get
#: a vote (GUARDRAILS). csoil's four models are structurally distinct, so each is its
#: own family; the map is explicit so a future orchidee/orchidee-dgvm pair is handled.
MODEL_FAMILY = {}


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """Extract (model, gcm, scenario, soc, sens, member) from an ISIMIP3b filename.

    e.g. lpjml5-7-10-fire_gfdl-esm4_w5e5_ssp126_2015soc-from-histsoc_default_csoil-total_global_annual_2015_2100.nc

    Parsed from the END. A forward index happens to work for OutputData -- model and GCM
    tokens carry hyphens, never underscores -- but it breaks silently the moment a file
    gains a leading token, which is exactly what DerivedOutputData filenames do: there
    `ISIMIP3b` lands on the forcing field and gets reported as the scenario. The variable
    is asserted rather than assumed, so a grammar change fails loudly instead of
    mislabelling every member.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def roughness(field2d, mask2d):
    """mean |cell - its 4-neighbour mean| over the mask, normalized by the mean.

    Same definition as `let`'s smoothing measurement, so the numbers are comparable across
    layers: raw `let` 0.389, `cropfailure` 0.347, `permafrost-3b` 0.114.
    """
    a = np.where(mask2d, field2d, np.nan)
    with warnings.catch_warnings():
        # A cell whose four neighbours are all off-mask has no neighbour mean; nanmean
        # says so once per such cell and floods the log. NaN is the right answer and it
        # drops out of the nanmean below.
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        nb = np.nanmean(np.stack([
            np.roll(a, 1, 0), np.roll(a, -1, 0),
            np.roll(a, 1, 1), np.roll(a, -1, 1)]), axis=0)
    d = np.abs(a - nb)
    m = np.nanmean(a)
    return float(np.nanmean(d) / m) if np.isfinite(m) and m > 0 else float("nan")


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


def slope_chunk_cells(n_members, n_years, jobs=1):
    """Cells per block, sized so the budget is a TOTAL across workers, not per worker.

    Getting this wrong is how the rebuild died twice: every forked worker sizes its own
    scratch from the same constant, so 8 workers meant 8 x 400 MB = 3.2 GB of slope scratch
    on top of a 1.24 GB cube, on a 16 GB machine with ~5 GB free. Dividing by `jobs` keeps
    peak scratch at the budget no matter how wide the pool is.
    """
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    per_cell = 4 * pairs * 4
    budget = SLOPE_MEM_BUDGET_BYTES // max(jobs, 1)
    return max(4, min(512, int(budget // max(per_cell, 1))))


def _slope_block(task):
    s, e, years, decade, baseline, window, chunk = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk)
    # Plain arrays, NOT the SlopeResult: its `__getattr__ = dict.__getitem__` turns
    # pickle's __getstate__ probe into a KeyError and kills the pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, jobs, n_land):
    """Expanding-window slopes, split over cells across processes.

    Identical output to the serial call -- a slope is fitted from each cell's own
    observations and cells are independent, so blocking over the flat cell axis cannot
    change a result, only the wall clock.
    """
    if decade == BASELINE_DECADE or jobs <= 1:
        return expanding_slopes(cube, years, decade, BASELINE_DECADE,
                                window_years=WINDOW_YEARS, chunk_cells=chunk)
    n_blocks = max(jobs * 8, 1)
    edges = np.linspace(0, n_land, n_blocks + 1).astype(int)
    tasks = [(int(a), int(b), years, decade, BASELINE_DECADE, WINDOW_YEARS, chunk)
             for a, b in zip(edges[:-1], edges[1:]) if b > a]
    ols = np.full(n_land, np.nan, np.float32)
    sen = np.full(n_land, np.nan, np.float32)
    ctx = mp.get_context("fork")
    with ctx.Pool(jobs) as pool:
        for s, e, o, sn in pool.imap_unordered(_slope_block, tasks):
            ols[s:e] = o
            sen[s:e] = sn
    return {"ols_slope": ols, "sen_slope": sen}


def main():
    global _CUBE
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2),
                    help="processes for the Theil-Sen stage (cells are independent, so "
                         "this changes wall clock only, never a value)")
    ap.add_argument("--members-only", action="store_true",
                    help="write only the per-member diagnostic and exit, skipping the "
                         "expensive Theil-Sen stage")
    ap.add_argument("--scenarios", nargs="*", default=None,
                    help="Assemble and WRITE only these scenarios. Every discovered "
                         "scenario is still loaded and pooled into the shared 2020s "
                         "baseline regardless, so a subset run is bit-identical to the "
                         "full run -- this RESUMES an interrupted run, it does not narrow "
                         "the baseline. Restricting the baseline instead would silently "
                         "break the cross-scenario comparability the contract requires.")
    args = ap.parse_args()

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

    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    all_members = members_by_scen[scenarios[0]]

    # ---- Branch check: is the pooled sample UNIMODAL? ------------------------
    # OUTPUT-SPEC's fourth branch exists because the median SELECTS the larger cluster when
    # members separate, and jumps when the balance tips. Four models at four levels is
    # exactly the shape that needs testing, so `pooled_median` is declined-into rather than
    # defaulted-to: each model's own 2020s field is built separately and the gap between
    # adjacent model medians is measured against the within-model IQR. `permafrost-3b`
    # needed branch 4 at a separation of 0.035 / 0.046 / 0.951.
    per_model = {}
    for m in models:
        idx = [i for i, mem in enumerate(all_members) if mem.rsplit("_", 1)[0] == m]
        pool = np.concatenate([annual[s][idx][:, bwin, :] for s in scenarios], axis=0)
        h, _, _ = pooled_decadal_stat(pool, years[bwin], BASELINE_DECADE,
                                      boolean=boolean, window_years=WINDOW_YEARS)
        fin = h[np.isfinite(h)]
        q1, q3 = np.percentile(fin, [25, 75])
        per_model[m] = dict(p50=float(np.median(fin)), q1=float(q1), q3=float(q3),
                            iqr=float(q3 - q1), n=int(fin.size))
        del pool, h
    order = sorted(per_model, key=lambda m: per_model[m]["p50"])
    gaps = []
    for k in range(len(order) - 1):
        a, b = order[k], order[k + 1]
        iqr = per_model[a]["iqr"]
        gaps.append((a, b, per_model[b]["p50"] - per_model[a]["p50"],
                     (per_model[b]["p50"] - per_model[a]["p50"]) / iqr if iqr > 0 else np.inf))
    max_rel = max((g[3] for g in gaps), default=0.0)
    multimodal = max_rel >= 1.0
    log("\nBranch check -- per-model 2020s medians (is the pool unimodal?):")
    for m in order:
        d = per_model[m]
        log(f"  {m:<20} p50={d['p50']:7.3f}  IQR {d['q1']:7.3f}-{d['q3']:<8.3f}")
    for a, b, gap, rel in gaps:
        log(f"    gap {a} -> {b}: {gap:+.3f} = {rel:.2f} x IQR"
            + ("   <-- SEPARATED" if rel >= 1.0 else ""))
    stat_rationale = (
        f"pooled_median retained ON MEASUREMENT, not by default. The four models' 2020s "
        f"medians are "
        + ", ".join(f"{m} {per_model[m]['p50']:.2f}" for m in order)
        + f" kg m-2 -- a {per_model[order[-1]]['p50'] / max(per_model[order[0]]['p50'], 1e-9):.2f}x "
          f"spread, but a GRADIENT, not clusters. The largest gap between adjacent model "
          f"medians is {max_rel:.2f}x the lower model's own IQR, below the ~1.0 at which the "
          f"median stops summarising the pool and starts selecting a cluster. Contrast "
          f"`permafrost-3b`, which took pooled_mean_multimodel with seven members at ~0.04 "
          f"and five at ~0.95. Zero-inflation is also below the degenerate regime "
          f"(baseline exact-zero {frac_zero:.2%}), so branch 3 does not apply either.")
    if multimodal:
        stat_rationale = (
            f"WARNING: the pooled sample IS separated (largest adjacent-model gap "
            f"{max_rel:.2f}x IQR) and pooled_median may be selecting a cluster rather than "
            f"summarising. OUTPUT-SPEC branch 4 (pooled_mean_multimodel) should be "
            f"considered with the user before this layer ships. " + stat_rationale)
        log(f"\n  WARNING: max adjacent-model gap {max_rel:.2f}x IQR -- branch 4 may apply.")

    # ---- Smoothing decision, measured on THIS layer --------------------------
    # Ensemble depth alone does not settle it -- `cropfailure-3b` had 400 draws per
    # cell-decade and roughness close to raw `let`. The split-half test does: if two
    # disjoint halves give the same roughness AND reproduce each other, the roughness is
    # real spatial structure rather than sampling noise. Halves are STRATIFIED BY MODEL
    # (alternating GCMs within each model); an alphabetical split gives the halves
    # different model compositions and measures that instead of noise.
    halves_mem = [[], []]
    for m in models:
        gs = sorted(mem.rsplit("_", 1)[1] for mem in all_members
                    if mem.rsplit("_", 1)[0] == m)
        for j, g in enumerate(gs):
            halves_mem[j % 2].append(f"{m}_{g}")
    halves = []
    for half in halves_mem:
        idx = [all_members.index(mm) for mm in half]
        pool = np.concatenate([annual[s][idx][:, bwin, :] for s in scenarios], axis=0)
        h, _, _ = pooled_decadal_stat(pool, years[bwin], BASELINE_DECADE,
                                      boolean=boolean, window_years=WINDOW_YEARS)
        halves.append(h)
        del pool
    both = np.isfinite(halves[0]) & np.isfinite(halves[1]) & np.isfinite(b_med)
    common2d = np.zeros((LAT, LON), bool)
    common2d.ravel()[land_idx[both]] = True
    rough_all = roughness(scatter(b_med, land_idx, (LAT, LON)), common2d)
    rough_h = [roughness(scatter(h, land_idx, (LAT, LON)), common2d) for h in halves]
    corr = (float(np.corrcoef(halves[0][both], halves[1][both])[0, 1])
            if both.sum() > 2 else float("nan"))
    rough_stable = abs(rough_h[0] - rough_h[1]) < 0.25 * max(rough_all, 1e-9)
    reproduces = np.isfinite(corr) and corr >= 0.9 and rough_stable
    smoothing_note = (
        f"none. Roughness on the {BASELINE_DECADE}s panel is {rough_all:.3f} (raw `let` "
        f"0.389, `cropfailure` 0.347, `permafrost-3b` 0.114). SPLIT-HALF TEST, halves "
        f"stratified by model (alternating GCMs within each model, {len(halves_mem[0])} vs "
        f"{len(halves_mem[1])} members) with all three roughnesses measured on the "
        f"{int(both.sum()):,} cells both halves cover: halves read {rough_h[0]:.3f} and "
        f"{rough_h[1]:.3f}, Pearson r={corr:.3f}. "
        + ("The structure reproduces across independent halves, so it is real spatial "
           "structure -- soil-carbon gradients follow biome, permafrost and wetland "
           "boundaries -- and a kernel would blur exactly what the layer reports."
           if reproduces else
           "The halves do NOT reproduce each other well, so part of this roughness is "
           "sampling noise and a case for smoothing exists. It is still NOT applied here: "
           "at 17 members the pool is thick, and the disagreement is between MODELS about "
           "where carbon is (see interpretation_caveat), which spatial averaging would "
           "hide rather than fix. Revisit with the user if the maps look noisy."))
    log(f"\nSmoothing: roughness {rough_all:.3f}; split-half {rough_h[0]:.3f}/"
        f"{rough_h[1]:.3f} on {int(both.sum()):,} common cells, r={corr:.3f}, "
        f"reproduces={reproduces} -> spatial_smoothing=none")
    del halves

    # ---- Per-member diagnostic ------------------------------------------------
    # NOT part of the OUTPUT-SPEC contract. It exists because every statistic in a
    # value-check table is invariant under spatial rearrangement, so no table can see a
    # spatial defect -- and this layer has a documented one: CLASSIC declares 0.5 deg but
    # is natively 1 deg, replicated 2x2 with a one-cell longitude offset. LPJmL has never
    # been looked at at all. Consumed by generate_maps.py's Members tab and by
    # render_contact_sheet.py, which renders a PNG that can actually be opened here.
    mem_names = members_by_scen[scenarios[0]]
    mem_grid = np.full((len(mem_names), LAT, LON), np.nan, np.float32)
    for mi, mname in enumerate(mem_names):
        stack = [annual[s][members_by_scen[s].index(mname)][bwin]
                 for s in scenarios if mname in members_by_scen[s]]
        if not stack:
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            flat = np.nanmean(np.concatenate(stack, axis=0), axis=0)
        mem_grid[mi] = scatter(flat, land_idx, (LAT, LON))
    mem_ds = xr.Dataset(
        {"value": (["member", "lat", "lon"], mem_grid)},
        coords={"member": mem_names, "lat": lats, "lon": lons},
        attrs={
            "variable": VAR,
            "units": "kg m-2",
            "member_field": (f"mean soil organic carbon stock over the {BASELINE_DECADE}s "
                             "baseline decade, pooled across all scenarios"),
            "note": ("Diagnostic only -- not part of the OUTPUT-SPEC contract. Consumed "
                     "by scripts/generate_maps.py to render the Members tab and by "
                     "scripts/render_contact_sheet.py."),
        },
    )
    mem_path = out_dir / f"{VAR.split('-')[0]}_members.nc"
    mem_ds.to_netcdf(mem_path, encoding={"value": {"dtype": "float32", "zlib": True,
                                                   "complevel": 4,
                                                   "_FillValue": np.float32(np.nan)}})
    log(f"\nwrote per-member diagnostic {mem_path.name} ({len(mem_names)} members)")
    del mem_grid
    if args.members_only:
        log("--members-only: done.")
        return

    # ---- Per-scenario assembly ----------------------------------------------
    # NOTE the asymmetry, and keep it: loading, the shared baseline, the branch check and
    # the split-half above all run over EVERY discovered scenario. Only the write loop is
    # narrowed by --scenarios, so a resumed run reproduces a full run exactly.
    b_idx = DECADES.index(BASELINE_DECADE)
    write_scenarios = ([s for s in scenarios if s in args.scenarios]
                       if args.scenarios else scenarios)
    if args.scenarios:
        missing = sorted(set(args.scenarios) - set(scenarios))
        if missing:
            log(f"ERROR: --scenarios names {missing}, not present in {raw_dir}")
            return
        log(f"\nWriting only {write_scenarios} (baseline still pooled over {scenarios})")
        # The baseline, branch check and split-half are all done by here, so the scenarios
        # we are not writing are dead weight -- ~412 MB each at 17 members. Freeing them is
        # the difference between fitting in memory and not.
        for s_drop in [x for x in scenarios if x not in write_scenarios]:
            del annual[s_drop]
    for s in write_scenarios:
        log(f"\n{'=' * 66}\nAssembling scenario {s}\n{'=' * 66}")
        mem = members_by_scen[s]
        cube = annual[s]                                   # (member, year, n_land)
        # Published to the fork so workers inherit the cube instead of pickling it.
        _CUBE = cube
        chunk = slope_chunk_cells(len(mem), n_years, args.jobs)
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

            sl = compute_slopes(cube, years, d, chunk, args.jobs, n_land)

            # THE MASK IS TIME-VARYING, so the slopes must be masked to each decade's own
            # median mask (OUTPUT-SPEC). The slope window EXPANDS from the baseline while
            # the median window is the decade alone, so a cell with observations early in
            # the window and none inside the decade gets a finite slope against a NaN
            # median -- a trend over a decade the subject was absent from, which the
            # mask-agreement guardrail below rejects, correctly.
            # csoil's mask is very nearly static: this fires on ONE cell, at ssp585 2060s,
            # against `npp-tempnle`'s 53 rising to 374. It is fixed rather than
            # special-cased because one leaking cell and 374 are the same defect, and on
            # npp-tempnle dropping 53 of 25,821 moved the mean slope from -1.89 to +0.64.
            gone = ~np.isfinite(med)
            leak = int(np.count_nonzero(
                gone & (np.isfinite(sl["ols_slope"]) | np.isfinite(sl["sen_slope"]))))
            if leak:
                sl = {"ols_slope": np.where(gone, np.nan, sl["ols_slope"]),
                      "sen_slope": np.where(gone, np.nan, sl["sen_slope"])}
                log(f"    masked {leak} slope cell(s) with no observation inside the "
                    f"{d}s median window")

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
                             # Subscript, not attribute: the parallel path returns a plain
                             # dict while the serial path returns a SlopeResult. Both are
                             # dict-like; only the dict lacks attribute access.
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
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
                # Derived from the coordinates, never hardcoded: resolution became a
                # per-layer property on 2026-08-14 when the flood layers shipped at 0.25 deg,
                # and a constant would be right until the day it silently was not.
                "spatial_resolution_degrees": round(float(abs(lats[1] - lats[0])), 6),
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "field_nature": "boolean_01" if boolean else "continuous",
                "value_note": (
                    f"median = MEDIAN over the pooled (year x member) sample inside the "
                    f"decade window, across the {len(mem)} ISIMIP3b biomes members "
                    f"({', '.join(models)} x their CMIP6 GCMs); raw ISIMIP3b values in "
                    f"kg C m-2."),
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
                    "none -- all 4 models report the SAME unit (kg C m-2) and the values "
                    "are physically meaningful and comparable to field observations, so "
                    "they are equal-weighted in raw kg C m-2 (model-democracy decision). "
                    "RE-MEASURED 2026-08-14 with lpjml5-7-10-fire added and the spread is "
                    "WIDER than the 3-model build recorded: 2020s medians 5.70 (classic) / "
                    "7.55 (mc2-usfs) / 10.45 (jules) / 16.64 (lpjml), a 2.92x range against "
                    "the 1.84x measured on three models, and lpjml carries the upper tail "
                    "alone (p95 91.1 vs 16.5-27.4). Normalization is still declined: the "
                    "models form a GRADIENT rather than clusters (see "
                    "decadal_statistic_rationale), the disagreement about the size of the "
                    "stock IS the structural uncertainty the CI is meant to carry, and "
                    "rescaling would destroy units that can be checked against measured "
                    "soil-carbon inventories. CLASSIC contributes only 2 GCMs against 5 for "
                    "each of the others, so it is underweighted in the flat member pool."),
                "decadal_statistic_rationale": stat_rationale,
                # THE ONE ATTRIBUTE THAT MUST REACH THE CUSTOMER. `LAYER_ATTRS_EXPORTED` in
                # scripts/utils/delivery.py is a CLOSED allowlist: a caveat written under any
                # other name is silently dropped before layers.csv, the caveat generator and
                # both reports. The 2026-07-25 build put this text in `co2_treatment`, which
                # is not on that list, so this layer's defining limitation would have
                # vanished between the file and the filing.
                "interpretation_caveat": (
                    "SOIL ORGANIC CARBON ONLY, AND THE CLIMATE SIGNAL ONLY. (1) This layer "
                    "measures the soil organic carbon stock, which is ONE of the ten soil "
                    "degradation processes enumerated in recital 4 of Directive (EU) "
                    "2025/2360; it says nothing about sealing, compaction, contamination, "
                    "nutrient balance, acidification, salinisation or erosion. Do not report "
                    "it as 'soil degradation'. (2) LAND USE IS HELD FIXED in every member "
                    "(2015soc-from-histsoc for classic/jules/lpjml, nat for mc2-usfs), so "
                    "the layer CANNOT see management-driven carbon loss -- tillage, "
                    "conversion, residue removal -- which is most of observed soil-carbon "
                    "loss. That is deliberate: it isolates the climate signal, and it is a "
                    "hard limit on any claim about soil condition in general. (3) MIXED CO2: "
                    "jules-es-vn6p3 publishes ONLY the fixed-2015-CO2 run for csoil-total, "
                    "so 5 of 17 members carry no CO2-fertilisation signal; the other 12 are "
                    "transient. MEASURED on ssp585, those 5 show the LARGEST relative loss "
                    "of the four models (-4.37%, against -2.75% lpjml, -0.05% mc2-usfs and "
                    "+0.79% classic), because removing fertilisation removes litter input. "
                    "The mixed treatment therefore makes the ensemble decline somewhat "
                    "STRONGER than a uniformly transient one would, not weaker -- the "
                    "opposite of what this build asserted before it was measured. Retained "
                    "to keep the ensemble (user decision 2026-07-25, re-affirmed on rebuild "
                    "2026-08-15). (4) THE MODELS DISAGREE ABOUT WHERE CARBON IS, most "
                    "sharply in high-latitude organic soils: at the 2020s North Slope of "
                    "Alaska lpjml reads 65.6 kg m-2 against 3.7-5.7 for the other three, and "
                    "lpjml alone projects large northern peat losses by the 2090s (W "
                    "Siberia -5.9, Hudson Bay -4.3, Sodankyla -4.2 kg m-2) where the others "
                    "are flat. The ensemble central value therefore UNDERSTATES a "
                    "permafrost-carbon risk that one of four models considers large -- read "
                    "the CI, not the median alone, at high latitudes."),
                "co2_treatment": (
                    "MIXED: classic, lpjml5-7-10-fire & mc2-usfs use transient CO2 (default "
                    "run); jules-es-vn6p3 publishes ONLY the fixed-2015-CO2 run for "
                    "csoil-total. MEASURED 2026-08-15: those members show the LARGEST "
                    "relative loss of the four models, not a muted trend. Retained to keep 17 "
                    "members (user decision 2026-07-25, re-affirmed 2026-08-14). NOTE: this "
                    "attribute is NOT on LAYER_ATTRS_EXPORTED and does not reach delivery -- "
                    "the customer-facing statement is in `interpretation_caveat`."),
                "spatial_smoothing": smoothing_note,
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
