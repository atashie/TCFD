"""ISIMIP3b `thawdepth` -> permafrost layer, normalized per model by its own soil column.

WHAT THIS LAYER PUBLISHES

    thawfrac = annual maximum thaw depth / that model's soil column depth,  in [0, 1]

i.e. THE FRACTION OF THE MODELLED SOIL COLUMN THAT IS THAWED at the end of the thaw season.
0 means the column never thaws; 1 means nothing frozen is left. The customer-facing quantity
is the CHANGE in that fraction against the 2020s baseline -- the share of the column that
transitions from permafrost to non-permafrost -- which is the dashboard's Anomaly panel and
the `transition_summary.json` written beside the processed files.

WHY NORMALIZED, AND WHY THAT IS NOT COSMETIC (measured 2026-08-14, all 36 files, see
scripts/check_thawdepth_nature.py and data/raw/.../value_check.json)

The raw metres CANNOT be pooled across these three models. Each publishes a different soil
column and pins its permafrost-free cells at that column's depth:

    model              column   mass at column   exact 0   Yakutsk 2020s
    CLASSIC             61.4 m      85.2%          0.3%       8.03 m
    LPJmL5-7-10-fire    13.0 m      69.2%          2.6%       0.94 m
    JULES-ES-VN6P3       3.0 m      80.2%          4.5%       2.995 m

One site, one decade, three answers spanning an order of magnitude -- and JULES's 2.995 m is
99.8% of its own column, so it is nearly censored where the hazard lives. A pooled median
over those numbers is an average of three soil discretisations, not a depth. Dividing by each
model's column puts every member on the same [0,1] axis, where 1 has ONE meaning in all three
("no frozen ground left in the modelled column") and the change is interpretable.

THE HONEST LIMIT OF THAT CHOICE, which must be disclosed rather than solved here: 0.5 of a
3 m column is not physically 0.5 of a 61.4 m column. Normalisation makes the ENDPOINTS
commensurable and the interior only ordinally so. That is why the reported quantity is the
transition -- a move from "some permafrost" to "none" means the same thing in every model,
while an intermediate depth does not (user decision 2026-08-14).

NO MEMBER IS DROPPED, AND THE MEASUREMENT THAT SETTLED IT. JULES's normalized values sit near
1 and look, in normalized space, like a censored member worth removing. Checked in RAW metres
against observed active-layer thickness (~0.3-1.5 m in continuous permafrost, up to ~2-3 m at
the warm margin), the ranking reverses -- 2020s thaw depth over each model's own domain:

    model              column     p25     p50     p75     p95   Fairbanks
    LPJmL5-7-10-fire    13.0 m   0.52    0.83    1.31    4.57      1.49 m
    JULES-ES-VN6P3       3.0 m   0.93    2.85    2.99    3.00      3.00 m (at its column)
    CLASSIC             61.4 m   1.20    2.15    4.90   28.00     61.40 m (no permafrost)

LPJmL is closest to observation; JULES is deep but not absurd; CLASSIC's p95 of 28 m is not
an active layer at all, and it reports NO permafrost at Fairbanks, which sits in the
discontinuous zone. The model with the most headroom is the least physical, so "has room to
move" is not evidence of quality. JULES's baseline pattern also agrees with the other two
about as well as they agree with each other (Spearman rho 0.785 with CLASSIC, 0.657 with
LPJmL, against 0.631 between CLASSIC and LPJmL), so it carries pattern information its
compressed level hides.

WHAT THE ENSEMBLE DOES NOT AGREE ON, which the CI must carry rather than hide: WHERE
permafrost is lost. By the 2090s under ssp585, over each model's own domain, all three lose a
comparable SHARE (CLASSIC 31.0% of 12.50 M km2, JULES 58.0% of 16.32, LPJmL 32.3% of 29.95)
-- but on the cells CLASSIC and LPJmL both call permafrost, CLASSIC loses 28.6% and LPJmL
4.0%, a Jaccard overlap of 2.8%. LPJmL's losses sit in the fringe it alone claims. Comparable
totals, different maps.

DECISION 1 -- MEAN +/- 1 SD, NOT MEDIAN/IQR (user decision 2026-08-14). Because the members
form two separated clusters in normalized space, the median reports whichever cluster holds
more members and JUMPS when the balance tips: under the median branch the ssp585 spatial
median went 0.40 (2080s) -> 0.93 (2090s), which reads as accelerating thaw and is a
cluster-crossing artifact. The mean moves smoothly between them and the SD carries the
disagreement measured above. This is a DECLARED deviation from the spec's median default for
continuous fields, for a reason the spec does not yet name (multimodality, not
zero-inflation); the run-time counterfactual is recorded in decadal_statistic_rationale.

THE FOOTPRINT IS THE 2020s PERMAFROST DOMAIN, NOT THE FINITE MASK

`isfinite` covers 27.5% of the grid including Nairobi and Paris, which read exactly the
column depth -- a permafrost-free cell is published AT the column, never as NaN. Publishing
on the finite mask would put a constant 1.0 across the tropics into the percentile baseline
and report a hazard where there is no ground ice at all.

So a member contributes a cell only where IT has permafrost in the 2020s baseline decade
(2020s maximum thawfrac strictly below its column). Cells are published where at least
`--min-models` model families satisfy that. Defined on the BASELINE ALONE and then held
fixed: a member that loses its permafrost by 2090 keeps contributing, rising to 1.0, because
that transition is the entire signal. Defining the footprint over the whole record instead
would delete exactly the cells the layer exists to report.

Domain areas measured on the 2020s baseline (area-weighted; cell counts are not comparable
at these latitudes). Against Obu et al. 2019: ~14 M km2 permafrost AREA, ~21 M km2 permafrost
REGION in the Northern Hemisphere.

    CLASSIC          12.27 M km2   plausible
    JULES-ES-VN6P3   14.72 M km2   plausible
    LPJmL5-7-10-fire 29.08 M km2   2x the observed area, above the whole region

The models also disagree about WHERE: 8.44 M km2 is called permafrost by all 12 members,
against a >=1-member union of 30.95 M km2. That disagreement, not the ensemble mean, is why
`--min-models` defaults to 2 rather than the full union -- a 1-model cell here is usually
LPJmL alone, in its measured over-extension. The tier table is re-measured at run time and
logged; change the flag if the user prefers the inclusive reading.

SATURATION IS THIS LAYER'S DEFINING DEFECT, and it is the `heatwave` problem from the other
end. Once a cell's column is fully thawed the value pins at 1.0 and BOTH slopes go to ~0 --
reading as "no trend" while meaning "no permafrost left". Trend on this layer must be read
with the Anomaly panel, never alone. The saturated share is measured per panel at run time
and recorded in `slope_definition`.

WHAT THIS LAYER IS NOT

Not ground stability, and not solifluction. Thaw depth says nothing about ice content, excess
ground ice, thaw settlement or slope, and ISIMIP publishes no variable for any of them -- so
a cell losing its permafrost here is not thereby a subsidence forecast. Not a hazard score
outside its own footprint either: a site with no permafrost has no thaw hazard, and is NaN
here rather than zero.

GUARDRAILS 12 reference sites (2020s -> 2090s ssp585, ensemble median thawfrac; measured at
ingest in raw metres, see check_thawdepth_nature.py): continuous permafrost at Yakutsk,
Tiksi, Norilsk, Deadhorse, Inuvik, Resolute and Longyearbyen; discontinuous at Fairbanks,
Yellowknife, Churchill and Salekhard; alpine at Nagqu. CLASSIC carries no Tibetan-plateau
permafrost at all (Nagqu pinned at its column throughout), so the alpine domain rests on the
other two models -- visible in `n_models`, and the reason that field must be read.

Usage:
    python scripts/process_thawdepth_permafrost.py --footprint-only   # tiers + areas, exit
    python scripts/process_thawdepth_permafrost.py [--min-models 2] [--jobs N]
"""

import argparse
import glob
import json
import multiprocessing as mp
import os
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.decadal_stats import (  # noqa: E402
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

VAR = "thawdepth"
OUT_VAR = "thawfrac"
LAYER_ID = "permafrost-isimip3b_thawdepth_annual"

#: ISIMIP3b future runs start in 2015, so 2020 is decade index 0 -- no 2010s panel.
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

#: The value is a fraction of the soil column.
CI_FLOOR, CI_CAP = 0.0, 1.0
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = False          # more of the column thawed = worse
SLOPE_PER_DECADE = 10.0

#: Relative slack for "is this value at the column depth". An ABSOLUTE tolerance reported
#: LPJmL's at-column share as 0.0% while its median sat on the column: that model's maximum
#: is 13.001 and the pinned mass is at 13.000. Any ceiling test on a float field needs slack
#: proportional to the ceiling.
CEIL_RTOL = 1e-4
#: In normalized units, "the column is fully thawed".
FULL_TOL = 1e-4

#: See the docstring -- measured, not inherited from `led` (>=2) or `driedarea` (union).
MIN_MODELS = 2

#: DECLARED DEVIATION -- see decision 1 in the docstring. Mean +/- 1 SD on a continuous,
#: NOT zero-inflated field, because the ensemble is MULTIMODAL across models and the median
#: therefore jumps between clusters instead of tracking the field. User decision 2026-08-14.
STAT_NAME = "pooled_mean_multimodel"
CENTRAL = "mean"

#: Three distinct model codes; no two configurations of one code, so family == model.
MODEL_FAMILY = {}

SLOPE_MEM_BUDGET_BYTES = 400 * 1024**2

_CUBE = None


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, member) read from the END of the filename.

    e.g. lpjml5-7-10-fire_gfdl-esm4_w5e5_ssp126_2015soc-from-histsoc_default_thawdepth_global_annual_2015_2100.nc

    OutputData carries no leading publication token, but parsing from the end costs nothing
    and is the only form that survives a publication that adds one.
    """
    p = os.path.basename(fpath).split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def decode_years(ds, fpath):
    """Years per record, from the filename span and CHECKED against the axis.

    The `time` units differ across members ('days since 1601-01-01' with and without a
    time-of-day, and '1601-1-1'), so the epoch string is not a reliable key. What IS uniform
    is a 365-day step, verified here; the span comes from the filename, which is part of the
    ISIMIP grammar.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    y0, y1 = int(p[-2]), int(p[-1])
    t = np.asarray(ds["time"].values, dtype="float64")
    n = y1 - y0 + 1
    if t.size != n:
        raise ValueError(f"{os.path.basename(fpath)}: {t.size} records but the filename "
                         f"declares {y0}-{y1} ({n} years)")
    d = np.diff(t)
    if d.size and ((d <= 0).any() or not (359.0 <= float(np.median(d)) <= 367.0)):
        raise ValueError(f"{os.path.basename(fpath)}: time axis is not an annual sequence "
                         f"(median step {float(np.median(d)) if d.size else float('nan')})")
    return np.arange(y0, y1 + 1, dtype=int)


def read_values(fpath):
    """Raw (year, lat, lon) metres with the declared fill replaced by NaN."""
    with xr.open_dataset(fpath, decode_times=False) as ds:
        da = ds[VAR]
        yrs = decode_years(ds, fpath)
        vals = da.values.astype("float32")
        fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    if fill is not None and np.isfinite(fill):
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan
    return yrs, vals


def column_depth(model_files):
    """This model's soil column depth, taken as the MODE of its upper tail.

    Not the maximum: LPJmL's maximum is 13.001 m while the pinned mass sits at exactly
    13.000, so normalising by the max would leave every permafrost-free cell at 0.99992 and
    no "fully thawed" test would fire. The mode of the top percentile is the depth the model
    actually reports when nothing is frozen, and it is asserted to agree with the maximum.
    """
    tops = []
    for f in model_files:
        _, v = read_values(f)
        fin = v[np.isfinite(v)]
        if fin.size == 0:
            continue
        cut = np.percentile(fin, 99.0)
        tail = np.round(fin[fin >= cut], 4)
        tops.append((Counter(tail.tolist()).most_common(1)[0][0], float(fin.max())))
        del v, fin, tail
    modes = sorted({m for m, _ in tops})
    vmax = max(x for _, x in tops)
    if len(modes) != 1:
        raise ValueError(f"soil column depth is not consistent across this model's files: "
                         f"{modes} -- members are not on one vertical grid")
    col = float(modes[0])
    if not (col <= vmax <= col * (1.0 + 1e-3)):
        raise ValueError(f"column mode {col} is not the top of the range (max {vmax}) -- "
                         f"the pinned mass is not the column depth")
    return col


def to_frac(vals, col):
    """Metres -> fraction of the soil column, clipped to [0, 1]."""
    return np.clip(vals / np.float32(col), 0.0, 1.0).astype(np.float32)


def baseline_permafrost(fpaths, col, years_of):
    """Cells where this member still has frozen ground in the 2020s baseline decade.

    The test is on the member's own column: the 2020s MAXIMUM thaw stays strictly below it.
    Union across the member's scenario files -- the 2020s is the shared baseline decade, so
    the three scenarios describe the same period and should agree; a union is the inclusive
    reading if they do not.
    """
    acc = None
    for f in fpaths:
        yrs, v = read_values(f)
        win = (yrs >= BASELINE_DECADE) & (yrs <= BASELINE_DECADE + WINDOW_YEARS - 1)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            mx = np.nanmax(v[win], axis=0)
        has = np.isfinite(mx) & (mx < col * (1.0 - CEIL_RTOL))
        acc = has if acc is None else (acc | has)
        del v, mx
    years_of.append(None)
    return acc


def cell_area_km2(lats, lons):
    """Area of each 0.5 deg cell. Cell counts are meaningless at 60-80 deg N."""
    R = 6371.0088
    dlat = float(abs(lats[1] - lats[0]))
    dlon = float(abs(lons[1] - lons[0]))
    col = (np.deg2rad(dlat) * R) * (np.deg2rad(dlon) * R * np.cos(np.deg2rad(lats)))
    return np.repeat(col[:, None], len(lons), axis=1).astype(np.float64)


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score against the shared 2020s baseline PERMAFROST distribution.

    The reference population is the published permafrost footprint, not all land, so a
    percentile here is a rank among permafrost cells. Two-tier when the baseline carries
    >2% exact zeros -- a cell that never thaws at all is the least exposed there is, and
    tier 1 says so. The mode is chosen from the measured share at run time.
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
            res = np.ones(vals.shape, np.float32)
            pos = vals > 0
            if n_nz > 0:
                frac = np.searchsorted(nz_sort, vals[pos], side="right") / n_nz
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
    """Put a footprint-cell vector back on the (lat, lon) grid; everything else NaN."""
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


def roughness(field2d, mask2d):
    """mean |cell - its 4-neighbour mean| over the footprint, normalized by the mean.

    Same definition as `let`'s smoothing measurement, so the numbers are comparable:
    raw `let` 0.389, `cropfailure` 0.347.
    """
    a = np.where(mask2d, field2d, np.nan)
    nb = np.nanmean(np.stack([
        np.roll(a, 1, 0), np.roll(a, -1, 0),
        np.roll(a, 1, 1), np.roll(a, -1, 1)]), axis=0)
    d = np.abs(a - nb)
    m = np.nanmean(a)
    return float(np.nanmean(d) / m) if np.isfinite(m) and m > 0 else float("nan")


def slope_chunk_cells(n_members, n_years, max_pairs):
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    if max_pairs is not None:
        pairs = min(pairs, max_pairs)
    per_cell = 4 * pairs * 4
    return max(4, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def _slope_block(task):
    s, e, years, decade, baseline, window, chunk, max_pairs = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk, max_pairs=max_pairs)
    # Plain arrays, NOT the SlopeResult: its `__getattr__ = dict.__getitem__` turns
    # pickle's __getstate__ probe into a KeyError and kills the pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, max_pairs, jobs, n_land):
    if decade == BASELINE_DECADE or jobs <= 1:
        return expanding_slopes(cube, years, decade, BASELINE_DECADE,
                                window_years=WINDOW_YEARS, chunk_cells=chunk,
                                max_pairs=max_pairs)
    n_blocks = max(jobs * 8, 1)
    edges = np.linspace(0, n_land, n_blocks + 1).astype(int)
    tasks = [(int(a), int(b), years, decade, BASELINE_DECADE, WINDOW_YEARS, chunk, max_pairs)
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
    ap.add_argument("--scenarios", nargs="*", default=None)
    ap.add_argument("--min-models", type=int, default=MIN_MODELS)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--footprint-only", action="store_true",
                    help="measure the permafrost domain and its tiers, then exit")
    ap.add_argument("--skip-slopes", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    out_dir = root / "data" / "processed" / LAYER_ID
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc")))
    if not files:
        log(f"ERROR: no {VAR} files in {raw_dir}")
        return 2

    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    write_scenarios = args.scenarios or scenarios

    log("=" * 74)
    log(f"thawdepth (ISIMIP3b/SSP permafrost) -> TCFD contract [{LAYER_ID}]")
    log("=" * 74)
    log(f"{len(files)} files | scenarios {scenarios} | writing {write_scenarios}")

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    LAT, LON = len(lats), len(lons)
    area = cell_area_km2(lats, lons)

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {int(y): i for i, y in enumerate(years)}

    members_by_scen = {s: [] for s in scenarios}
    for f in files:
        members_by_scen[meta[f]["scenario"]].append(meta[f]["member"])
    for s in scenarios:
        members_by_scen[s] = sorted(members_by_scen[s])

    models = sorted({meta[f]["model"] for f in files})
    gcms = sorted({meta[f]["gcm"] for f in files})
    socs = sorted({meta[f]["soc"] for f in files})
    senss = sorted({meta[f]["sens"] for f in files})

    # ---- Pass 1: per-model soil column depth --------------------------------- #
    t0 = time.time()
    log("\nSoil column depth per model (mode of the upper tail, asserted to top the range):")
    column = {}
    for m in models:
        column[m] = column_depth([f for f in files if meta[f]["model"] == m])
        log(f"  {m:<18} {column[m]:>8.3f} m")

    # ---- Pass 2: 2020s permafrost footprint ---------------------------------- #
    all_members = sorted({meta[f]["member"] for f in files})
    files_by_member = {m: [f for f in files if meta[f]["member"] == m]
                       for m in all_members}
    mem_mask = {m: baseline_permafrost(files_by_member[m],
                                       column[m.rsplit("_", 1)[0]], [])
                for m in all_members}

    model_mask = {m: np.zeros((LAT, LON), bool) for m in models}
    for mem, mk in mem_mask.items():
        model_mask[mem.rsplit("_", 1)[0]] |= mk
    nmod_static = np.sum([model_mask[m] for m in models], axis=0).astype(np.int16)
    union = nmod_static > 0
    keep2d = nmod_static >= args.min_models

    model_area = {m: float((model_mask[m] * area).sum()) / 1e6 for m in models}
    log(f"\n{BASELINE_DECADE}s permafrost domain ({time.time() - t0:.0f}s) -- "
        f"area-weighted, Obu 2019 reference: ~14 M km2 permafrost area, ~21 M km2 region:")
    for m in models:
        log(f"  {m:<18} {int(model_mask[m].sum()):>7,} cells  {model_area[m]:>6.2f} M km2")
    for k in range(1, len(models) + 1):
        sel = nmod_static == k
        if sel.any():
            log(f"  exactly {k} model(s): {int(sel.sum()):>7,} cells  "
                f"{float((sel * area).sum()) / 1e6:>6.2f} M km2")
    log(f"  union {int(union.sum()):,} cells "
        f"({float((union * area).sum()) / 1e6:.2f} M km2)  ->  publishing >= "
        f"{args.min_models} model(s): {int(keep2d.sum()):,} cells "
        f"({float((keep2d * area).sum()) / 1e6:.2f} M km2); dropped "
        f"{int((union & ~keep2d).sum()):,}")
    if args.footprint_only:
        log("\n--footprint-only: done.")
        return 0

    land_idx = np.flatnonzero(keep2d.ravel())
    n_land = land_idx.size

    # ---- Pass 3: pack (member, year, footprint cell) as NORMALIZED fractions -- #
    annual = {s: np.full((len(members_by_scen[s]), n_years, n_land), np.nan, np.float32)
              for s in scenarios}
    slot = {s: {m: i for i, m in enumerate(members_by_scen[s])} for s in scenarios}
    # Each member is masked to ITS OWN baseline permafrost domain: a model with no permafrost
    # in a cell contributes NO observation there rather than a structural 1.0, which would
    # otherwise drag the pooled median to the ceiling and make n_models read 3 everywhere.
    mem_keep = {m: mem_mask[m].ravel()[land_idx] for m in all_members}
    for f in files:
        info = meta[f]
        s, m = info["scenario"], info["member"]
        col = column[info["model"]]
        yrs, cube = read_values(f)
        flat = to_frac(cube, col).reshape(cube.shape[0], -1)
        drop = ~mem_keep[m]
        for k, y in enumerate(yrs):
            yi = y_index.get(int(y))
            if yi is not None:
                row = flat[k, land_idx]
                row[drop] = np.nan
                annual[s][slot[s][m], yi] = row
        del cube, flat
    nmem_static = np.sum([mem_keep[m] for m in all_members], axis=0)
    log(f"\nPacked {len(files)} members over {n_land:,} footprint cells "
        f"({sum(a.nbytes for a in annual.values()) / 1024**3:.2f} GB resident) "
        f"[{time.time() - t0:.0f}s]")
    for lo, hi in ((1, 3), (4, 6), (7, 9), (10, 11), (12, 12)):
        c = int(((nmem_static >= lo) & (nmem_static <= hi)).sum())
        if c:
            log(f"    {lo:>2}-{hi:<2} members: {c:>7,}")

    # ---- Field nature: measured, never assumed (GUARDRAILS 9) ---------------- #
    boolean = is_boolean_field(annual[scenarios[0]])
    log(f"\nField nature: {'BOOLEAN {0,1}' if boolean else 'CONTINUOUS'}")
    if boolean:
        log("  ERROR: thawfrac measured BOOLEAN -- the inputs changed; re-run "
            "scripts/check_thawdepth_nature.py before processing.")
        return 3
    fin0 = annual[scenarios[0]][np.isfinite(annual[scenarios[0]])]
    frac_full_annual = float((fin0 >= 1.0 - FULL_TOL).mean())
    log(f"  annual cell-year share at a FULLY THAWED column, inside the footprint: "
        f"{frac_full_annual:.2%}  (0 by construction in the 2020s; this is 2020-2099)")
    log(f"  annual exact-zero share (never thaws): {float((fin0 == 0).mean()):.2%}")
    del fin0

    # ---- Shared 2020s baseline ------------------------------------------------ #
    uniform = len({tuple(members_by_scen[s]) for s in scenarios}) == 1
    if not uniform:
        log("\nWARNING: ensemble composition differs across scenarios.")
    bwin = (years >= BASELINE_DECADE) & (years <= BASELINE_DECADE + WINDOW_YEARS - 1)
    base_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years[bwin], BASELINE_DECADE, window_years=WINDOW_YEARS, central=CENTRAL)
    del base_pool
    b_lo = np.clip(b_lo, CI_FLOOR, CI_CAP)
    b_hi = np.clip(b_hi, CI_FLOOR, CI_CAP)

    # The declared deviation must be justified by NUMBERS recomputed at run time, never
    # quoted from a docstring: what would the median branch have published here?
    med_pool = np.concatenate([annual[s][:, bwin, :] for s in scenarios], axis=0)
    alt_med, _, _ = pooled_decadal_stat(
        med_pool, years[bwin], BASELINE_DECADE, window_years=WINDOW_YEARS, central="median")
    del med_pool
    alt_fin = alt_med[np.isfinite(alt_med)]
    this_fin = b_med[np.isfinite(b_med)]
    # A median of a two-cluster sample sits IN one cluster; the mean sits between them. The
    # tell is the gap: how far the published central value is from the nearer cluster.
    rationale = (
        f"MEAN +/- 1 SD on a continuous, NOT zero-inflated field -- a declared deviation "
        f"from the OUTPUT-SPEC median/IQR default for continuous layers (user decision "
        f"2026-08-14). The reason is MULTIMODALITY ACROSS MODELS, not zero-inflation: after "
        f"per-model column normalisation JULES-ES-VN6P3 sits at 0.95 of its 3 m column in "
        f"the 2020s while CLASSIC sits at 0.035 and LPJmL at 0.046, so the pooled sample is "
        f"two separated clusters and the MEDIAN reports whichever cluster holds more "
        f"members. That is not a rounding difference: under the median branch the ssp585 "
        f"spatial median jumped 0.40 (2080s) -> 0.93 (2090s) purely because the median "
        f"CROSSED between clusters, which reads as accelerating thaw and is not. Measured on "
        f"this run's {BASELINE_DECADE}s panel: the median branch would publish a spatial "
        f"median of {np.median(alt_fin):.4f} against the mean branch's "
        f"{np.median(this_fin):.4f}. The mean is a compromise between clusters -- a value no "
        f"single model produces -- and the SD is doing real work here: it carries genuine "
        f"inter-model disagreement, not noise. Read the CI on this layer.")
    del alt_med, alt_fin, this_fin

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"\nShared {BASELINE_DECADE}s baseline: n={baseline_flat.size:,}, "
        f"median thawed fraction {np.nanmedian(b_med):.4f}, "
        f"exact-zero {frac_zero:.2%}, percentile mode={pct_mode}")

    # Smoothing decision, measured on this layer rather than inherited. The split-half test
    # settles what the draw count alone cannot: if two disjoint halves give the same
    # roughness and correlate, the roughness is real spatial structure, not sampling noise.
    # The halves are STRATIFIED BY MODEL -- alternate GCMs within each model -- not split
    # down the sorted member list. An alphabetical split gave one half 3 JULES + 2 LPJmL and
    # the other 2 + 3, so the two halves had different permafrost DOMAINS and the test
    # measured model composition rather than sampling noise (r=0.376, and each half read
    # SMOOTHER than the full ensemble, which is the tell).
    halves_mem = [[], []]
    for m in models:
        gs = sorted(g for mem in all_members if mem.rsplit("_", 1)[0] == m
                    for g in [mem.rsplit("_", 1)[1]])
        for j, g in enumerate(gs):
            halves_mem[j % 2].append(f"{m}_{g}")
    halves = []
    for half in halves_mem:
        pool = np.concatenate(
            [annual[s][[members_by_scen[s].index(m) for m in half]][:, bwin, :]
             for s in scenarios], axis=0)
        h, _, _ = pooled_decadal_stat(pool, years[bwin], BASELINE_DECADE,
                                      window_years=WINDOW_YEARS, central=CENTRAL)
        halves.append(h)
        del pool
    # All three roughnesses on the SAME cells, or the comparison is between different maps.
    both = np.isfinite(halves[0]) & np.isfinite(halves[1]) & np.isfinite(b_med)
    common2d = np.zeros((LAT, LON), bool)
    common2d.ravel()[land_idx[both]] = True
    rough_all = roughness(scatter(b_med, land_idx, (LAT, LON)), common2d)
    rough_h = [roughness(scatter(h, land_idx, (LAT, LON)), common2d) for h in halves]
    corr = float(np.corrcoef(halves[0][both], halves[1][both])[0, 1]) if both.sum() > 2 \
        else float("nan")
    rough_stable = abs(rough_h[0] - rough_h[1]) < 0.25 * max(rough_all, 1e-9)
    reproduces = np.isfinite(corr) and corr >= 0.9 and rough_stable
    partial = np.isfinite(corr) and 0.7 <= corr < 0.9 and rough_stable
    smoothing_note = (
        f"none. Roughness on the {BASELINE_DECADE}s panel is {rough_all:.3f} "
        f"(raw `let` 0.389, `cropfailure` 0.347). SPLIT-HALF TEST, halves stratified by "
        f"model (alternating GCMs within each model, {len(halves_mem[0])} vs "
        f"{len(halves_mem[1])} members) and all three roughnesses measured on the "
        f"{int(both.sum()):,} cells both halves cover: halves read {rough_h[0]:.3f} and "
        f"{rough_h[1]:.3f}, Pearson r={corr:.3f}. "
        + ("The structure reproduces across independent halves, so it is real spatial "
           "structure rather than sampling noise."
           if reproduces else
           "Roughness is stable across the halves -- so it is NOT a thin-ensemble artifact "
           "that averaging would remove -- but they correlate less than `cropfailure`'s "
           "0.977, so a minority of it is sampling noise. Smoothing is still not applied, "
           "on a second and independent ground: the roughness concentrates at the permafrost "
           "boundary, which is precisely the line this layer reports moving, and a 5x5 "
           "kernel would smear permafrost into cells the models agree have none."
           if partial else
           "THE HALVES DO NOT REPRODUCE EACH OTHER, so part of this roughness IS sampling "
           "noise and a case for smoothing exists on those grounds. It is still not applied: "
           "the roughness concentrates at the permafrost boundary, which is precisely the "
           "line this layer reports moving, and a 5x5 kernel would smear permafrost into "
           "cells the models agree have none. Revisit with the user if the maps look noisy.")
    )
    log(f"  roughness {rough_all:.3f}; split-half {rough_h[0]:.3f}/{rough_h[1]:.3f}, "
        f"r={corr:.3f}, reproduces={reproduces} -> spatial_smoothing=none")
    del halves

    # ---- Per-member diagnostic for the dashboard Members tab ------------------ #
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
    xr.Dataset(
        {"value": (["member", "lat", "lon"], mem_grid)},
        coords={"member": mem_names, "lat": lats, "lon": lons},
        attrs={
            "variable": OUT_VAR, "units": "1",
            "member_field": (f"mean thawed fraction of the soil column over the "
                             f"{BASELINE_DECADE}s baseline, pooled across scenarios"),
            "note": ("Diagnostic only. Each member is normalized by ITS OWN soil column "
                     f"({'; '.join(f'{m} {column[m]:.1f} m' for m in models)}), so levels "
                     "are comparable but the underlying depths are not. Look for footprint "
                     "differences -- the models disagree about where permafrost IS."),
        },
    ).to_netcdf(out_dir / f"{OUT_VAR}_members.nc",
                encoding={"value": {"dtype": "float32", "zlib": True, "complevel": 4,
                                    "_FillValue": np.float32(np.nan)}})
    log(f"  wrote {OUT_VAR}_members.nc ({len(mem_names)} members)")
    del mem_grid

    for s_drop in [s for s in scenarios if s not in write_scenarios]:
        del annual[s_drop]

    chunk = slope_chunk_cells(max(len(m) for m in members_by_scen.values()),
                              n_years, args.max_pairs)
    jobs = max(1, args.jobs)
    log(f"Theil-Sen chunk_cells={chunk}, jobs={jobs}")

    area_flat = area.ravel()[land_idx]
    transition = {}

    # ---- Per-scenario assembly ------------------------------------------------ #
    for s in write_scenarios:
        log(f"\n{'=' * 74}\nAssembling scenario {s}\n{'=' * 74}")
        mem = members_by_scen[s]
        cube = annual[s]
        _CUBE = cube
        fams = sorted({family_of(m.rsplit("_", 1)[0]) for m in mem})
        fam_idx = {f: k for k, f in enumerate(fams)}

        nd = len(DECADES)
        out = {k: np.full((nd, LAT, LON), np.nan, np.float32)
               for k in ("median", "lower_ci", "upper_ci", "percentile",
                         "ols_slope", "sen_slope", "n_members", "n_models")}
        med_flat = {}

        for i, d in enumerate(DECADES):
            td = time.time()
            if d == BASELINE_DECADE:
                med, lo, hi, pc = b_med, b_lo, b_hi, b_pct
            else:
                med, lo, hi = pooled_decadal_stat(
                    cube, years, d, window_years=WINDOW_YEARS, central=CENTRAL)
                lo = np.clip(lo, CI_FLOOR, CI_CAP)
                hi = np.clip(hi, CI_FLOOR, CI_CAP)
                pc = pct(med)
            med_flat[d] = med

            if args.skip_slopes:
                sl = dict(ols_slope=np.full(n_land, np.nan, np.float32),
                          sen_slope=np.full(n_land, np.nan, np.float32))
            else:
                sl = compute_slopes(cube, years, d, chunk, args.max_pairs, jobs, n_land)

            win = (years >= d) & (years <= d + WINDOW_YEARS - 1)
            present = np.isfinite(cube[:, win, :]).any(axis=1)
            n_mem = present.sum(axis=0).astype(np.float32)
            fam_present = np.zeros((len(fams), n_land), bool)
            for mi, m in enumerate(mem):
                fam_present[fam_idx[family_of(m.rsplit("_", 1)[0])]] |= present[mi]
            n_mod = fam_present.sum(axis=0).astype(np.float32)
            n_mem[n_mem == 0] = np.nan
            n_mod[np.isnan(n_mem)] = np.nan

            for key, vec in (("median", med), ("lower_ci", lo), ("upper_ci", hi),
                             ("percentile", pc),
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                warnings.filterwarnings("ignore", message="All-NaN slice encountered")
                sat = float(np.nanmean(med >= 1.0 - FULL_TOL))
                slope_txt = ("slopes=NaN (baseline)" if d <= BASELINE_DECADE else
                             f"ols={np.nanmean(out['ols_slope'][i]):+.4f} "
                             f"sen={np.nanmean(out['sen_slope'][i]):+.4f} /dec")
                log(f"  {d}s: median thawed frac {np.nanmedian(med):.4f}  "
                    f"saturated {sat:.1%}  {slope_txt}  [{time.time() - td:.0f}s]")

        # ---- GUARDRAIL: slope and median masks must agree -------------------- #
        for i, d in enumerate(DECADES):
            if d <= BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), f"ols must be NaN at {d}s"
                assert np.all(np.isnan(out["sen_slope"][i])), f"sen must be NaN at {d}s"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), f"{k} finite where median is NaN at {d}s"
            assert np.all(out["lower_ci"][i][med_finite]
                          <= out["median"][i][med_finite] + 1e-5), f"lower_ci > median {d}s"
            assert np.all(out["upper_ci"][i][med_finite]
                          >= out["median"][i][med_finite] - 1e-5), f"upper_ci < median {d}s"

        # ---- The reported quantity: permafrost -> no permafrost -------------- #
        base = med_flat[BASELINE_DECADE]
        had = np.isfinite(base) & (base < 1.0 - FULL_TOL)
        rows = []
        for d in DECADES:
            m = med_flat[d]
            delta = m - base
            ok = np.isfinite(delta) & had
            lost = ok & (m >= 1.0 - FULL_TOL)
            rows.append(dict(
                decade=d,
                mean_transition=float(np.average(delta[ok], weights=area_flat[ok])),
                p90_transition=float(np.percentile(delta[ok], 90)) if ok.any() else np.nan,
                area_lost_Mkm2=float(area_flat[lost].sum()) / 1e6,
                area_lost_frac=float(area_flat[lost].sum() / area_flat[had].sum()),
                cells_lost=int(lost.sum()),
            ))
        transition[s] = dict(
            baseline_area_Mkm2=float(area_flat[had].sum()) / 1e6,
            baseline_cells=int(had.sum()), decades=rows)
        log(f"\n  Column fraction transitioning permafrost -> none, vs the "
            f"{BASELINE_DECADE}s (area-weighted over {float(area_flat[had].sum())/1e6:.2f} "
            f"M km2 of baseline permafrost):")
        log(f"    {'decade':>7} {'mean dFrac':>11} {'90th pct':>9} "
            f"{'fully lost':>12} {'of baseline':>12}")
        for r in rows:
            log(f"    {r['decade']:>7} {r['mean_transition']:>11.4f} "
                f"{r['p90_transition']:>9.4f} {r['area_lost_Mkm2']:>9.2f} M km2 "
                f"{r['area_lost_frac']:>11.1%}")

        fin_sen = out["sen_slope"][-1][np.isfinite(out["sen_slope"][-1])]
        fin_ols = out["ols_slope"][-1][np.isfinite(out["ols_slope"][-1])]
        sen_zero = float((fin_sen == 0).mean()) if fin_sen.size else float("nan")
        active = (fin_sen != 0) | (fin_ols != 0)
        sen_zero_active = (float((fin_sen[active] == 0).mean())
                           if active.any() else float("nan"))
        sat_final = float(np.nanmean(med_flat[DECADES[-1]] >= 1.0 - FULL_TOL))

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "source_variable": VAR,
                "scenario": s,
                "long_name": "Thawed fraction of the modelled soil column (annual maximum)",
                "units": "1",
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": STAT_NAME,
                "decadal_statistic_rationale": rationale,
                "field_nature": "continuous",
                "value_note": (
                    "median = MEAN of the pooled (year x member) sample inside the decade "
                    f"window over {len(mem)} members (3 land models x {{5,5,2}} CMIP6 GCMs). "
                    "The value is the ANNUAL MAXIMUM THAW DEPTH DIVIDED BY THAT MODEL'S OWN "
                    "SOIL COLUMN DEPTH, in [0,1]: 0 = the column never thaws, 1 = nothing "
                    "frozen is left. The raw ISIMIP field is metres (long_name 'Annual "
                    "Maximum Thaw Depth' / 'Maximum Thaw Depth'; cell_methods is ABSENT in "
                    "all 36 files, so the 'maximum' reading rests on long_name alone)."),
                "normalization": (
                    "PER-MODEL SOIL COLUMN DEPTH -- the load-bearing choice on this layer. "
                    f"Measured columns: {'; '.join(f'{m} {column[m]:.3f} m' for m in models)}. "
                    "The raw metres are NOT comparable across these models: each pins its "
                    "permafrost-free cells at its own column (mass at the column 85.2% / "
                    "69.2% / 80.2% for CLASSIC / LPJmL / JULES), and Yakutsk in the 2020s "
                    "reads 8.03 m, 0.94 m and 2.995 m in the three. Normalising makes the "
                    "ENDPOINTS commensurable -- 1 means 'no frozen ground left' in every "
                    "model -- and the interior only ORDINALLY so: 0.5 of a 3 m column is not "
                    "physically 0.5 of a 61.4 m column. That is why the reported quantity is "
                    "the CHANGE (see transition_metric), which is a transition between two "
                    "states that mean the same thing everywhere."),
                "transition_metric": (
                    "THE CUSTOMER-FACING QUANTITY IS THE CHANGE AGAINST THE 2020s: "
                    "median(decade) - median(2020s) is the fraction of the soil column that "
                    "transitions from permafrost to non-permafrost. It is the dashboard's "
                    "Anomaly panel, and the per-decade area-weighted figures are in "
                    "transition_summary.json beside this file. A cell reaching 1.0 has lost "
                    "its modelled permafrost entirely."),
                "ci_definition": (
                    "lower_ci/upper_ci = MEAN -/+ 1 standard deviation of the same pooled "
                    f"(year x member) sample, clipped to [{CI_FLOOR}, {CI_CAP}]. THE CI IS "
                    "NOT DECORATION ON THIS LAYER AND MUST BE READ. It is wide because the "
                    "three models genuinely disagree, in every way a model can: on the "
                    "normalized level (2020s medians 0.035 / 0.046 / 0.95 for CLASSIC / "
                    "LPJmL / JULES), on the domain (12.5 / 30.0 / 16.3 M km2 of 2020s "
                    "permafrost), and on WHERE permafrost is lost -- on the cells CLASSIC "
                    "and LPJmL share, CLASSIC loses 28.6% of them by the 2090s under ssp585 "
                    "and LPJmL loses 4.0%, a Jaccard overlap of 2.8%. A narrow CI on this "
                    "layer would mean the ensemble had been trimmed, not that the answer is "
                    "certain."),
                "slope_definition": (
                    "ols_slope = least-squares; sen_slope = Theil-Sen. Both fitted over an "
                    "EXPANDING window from the start of the 2020s baseline through the end "
                    "of the target decade, stacking every (year, member) observation. The "
                    "2020s panel is NaN. THIS LAYER SATURATES: a cell whose column is fully "
                    f"thawed pins at 1.0, and {sat_final:.1%} of published cells are there "
                    f"in the final panel, where BOTH slopes go to ~0 and mean 'no permafrost "
                    f"left', NOT 'no trend'. sen_slope is exactly 0 on {sen_zero:.1%} of "
                    f"finite cells ({sen_zero_active:.1%} of ACTIVE cells). Read the trend "
                    "WITH the Anomaly panel on this layer, never alone."),
                "slope_units": "1 decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: each cell's decadal thawed fraction ranked against the "
                    f"shared {BASELINE_DECADE}s ensemble distribution over the PERMAFROST "
                    "FOOTPRINT, not all land. A percentile here is a rank among permafrost "
                    "cells; a site outside the footprint is NaN, not zero, because a place "
                    "with no permafrost has no thaw hazard."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models),
                "gcms": ",".join(gcms),
                "soc_treatment": (
                    f"UNIFORM {','.join(socs)}. This is the ONLY soc all three models share; "
                    "the repository's usual 2015soc preference would drop JULES-ES-VN6P3 "
                    "entirely (7 members, 2 models instead of 12 and 3)."),
                "co2_treatment": f"UNIFORM {','.join(senss)}",
                "spatial_smoothing": smoothing_note,
                "minimum_models": args.min_models,
                "mask_rule": (
                    "THE FOOTPRINT IS THE 2020s PERMAFROST DOMAIN, NOT THE FINITE MASK. A "
                    "permafrost-free cell is published AT the column depth, never as NaN "
                    "(Nairobi and Paris read exactly the column in all three models), so "
                    "`isfinite` covers 27.5% of the grid and carries no domain information. "
                    "A member contributes a cell only where IT has frozen ground in the "
                    "2020s -- 2020s maximum thaw strictly below its own column -- and that "
                    "footprint is then HELD FIXED, so a member losing its permafrost later "
                    "keeps contributing and rises to 1.0. Measured 2020s domains (union "
                    "over each model's GCMs and scenarios): "
                    + "; ".join(f"{m} {model_area[m]:.2f}" for m in models)
                    + " M km2, against Obu et al. 2019's ~14 M km2 permafrost area and ~21 "
                    "M km2 permafrost region -- LPJmL alone exceeds the entire observed "
                    "region, which is why cells are published only where >= "
                    f"{args.min_models} of 3 models agree. Union is "
                    f"{float((union * area).sum()) / 1e6:.2f} M km2; published "
                    f"{float((keep2d * area).sum()) / 1e6:.2f} M km2, which sits between the "
                    "observed permafrost area and the observed permafrost region. Read "
                    "n_models: the "
                    "alpine domain rests on 2 models because CLASSIC carries no "
                    "Tibetan-plateau permafrost at all."),
                "interpretation_caveat": (
                    "THIS IS NOT A GROUND-STABILITY OR SUBSIDENCE LAYER. Thaw depth says "
                    "nothing about ice content, excess ground ice, thaw settlement or slope, "
                    "and ISIMIP publishes no variable for any of them, so a cell losing its "
                    "permafrost here is not thereby a foundation-damage forecast -- the "
                    "damage depends on ground ice this product does not carry. Solifluction "
                    "is a separate hazard with NO representation anywhere in ISIMIP. The "
                    "layer is also silent outside its footprint by design: a site with no "
                    "permafrost is NaN, which is an absence of this hazard, not a low score."),
                "source_dataset": (
                    "ISIMIP3b OutputData, annual `thawdepth`, w5e5, 2015-2100. Taken from "
                    "the `permafrost` sector for LPJmL5-7-10-fire and from `biomes` for "
                    "JULES-ES-VN6P3 and CLASSIC, which publish thawdepth but are NOT in the "
                    "permafrost sector at all -- a sector-only walk answers '1 model, 5 "
                    "members' to a question whose answer is 3 models and 12 members. The "
                    "files are byte-identical across the permafrost, biomes and water_global "
                    "sectors (Content-Length AND ETag verified 2026-08-14), so each member "
                    "is ingested exactly once. This product publishes NO sidecars, so inputs "
                    "are verified by Content-Length and the recorded sha512 is computed "
                    "locally -- it is not a publisher checksum."),
                "description": (
                    "ISIMIP3b/SSP permafrost thaw processed to the TCFD output contract "
                    f"(OUTPUT-SPEC.md) with a shared 2020s baseline; {len(mem)}-member "
                    f"3-model x CMIP6-GCM ensemble, {STAT_NAME} on the per-model normalized "
                    f"thawed column fraction, no smoothing, >= {args.min_models}-model "
                    f"permafrost footprint, {pct_mode} percentile, higher_is_worse."),
            },
        )

        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path.name} ({path.stat().st_size / 1024**2:.1f} MB)")
        _CUBE = None

    with open(out_dir / "transition_summary.json", "w") as fh:
        json.dump(dict(
            layer=LAYER_ID, variable=OUT_VAR, baseline_decade=BASELINE_DECADE,
            soil_column_m=column, min_models=args.min_models,
            metric=("median(decade) - median(2020s) of the thawed fraction of the soil "
                    "column; area-weighted over the 2020s permafrost footprint"),
            scenarios=transition), fh, indent=2)
    log(f"\nwrote transition_summary.json")
    log("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
