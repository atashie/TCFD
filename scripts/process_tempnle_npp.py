"""Process temperate needleleaf evergreen NPP (ISIMIP2b) into the TCFD output contract.

Temperate evergreen conifer productivity -- Douglas-fir, loblolly/ponderosa pine, Norway
spruce, radiata pine. Each model names the PFT differently; there is no shared code:

    CLM45     npp-needleleaf-evergreen-tree-temperate
    ORCHIDEE  npp-tendev
    LPJmL     npp-temperate-needleleaved-evergreen-tree

These are the three genuinely TEMPERATE classes with ISIMIP2b scenarios. ISIMIP3b has
exactly one temperate needleleaf class (LPJmL5-7-10-fire `tene`) and it publishes no
usable NPP -- npp-tene exists for ukesm1-0-ll ONLY, a 1-of-5 fragment (verified
2026-08-12 over every 3b model x GCM in biomes+fire+permafrost). CLASSIC `evgndltr` and
JULES `ndlevg` do have SSP NPP but merge boreal with temperate. Hence ISIMIP2b/RCP.

ENSEMBLE (2005soc / co2 transient, annual, 2006-2099):

    scenario   ORCHIDEE   LPJmL   CLM45                       total
    rcp26      4 GCMs     4 GCMs  hadgem2-es, miroc5          10
    rcp60      4 GCMs     4 GCMs  ipsl-cm5a-lr                 9
    rcp85      4 GCMs     4 GCMs  gfdl-esm2m, hadgem2-es      10

CLM45 publishes only 7 files for this PFT across all future scenarios and NO GCM of its
own appears in all three RCPs, so composition is scenario-dependent. Per user decision
2026-08-12 the shared 2020s baseline is pooled over EVERY (year x member x scenario)
observation in the window -- one panel, written bit-identically into all three files --
and `members_by_scenario` records the per-scenario identity.

*** DENOMINATOR HARMONIZATION -- the load-bearing choice, measured 2026-08-12 ***

The three models do NOT report NPP on the same denominator. Measured, not assumed:

    model     cover published        corr(cover, npp)   verdict
    CLM45     none (no pft- at all)  --                 per-tile (reports only on its tile)
    ORCHIDEE  yes, as a FRACTION     +0.279             per-tile
    LPJmL     yes, as true PERCENT   +0.898             COVER-SCALED (per grid cell)

Across cover quintiles LPJmL's NPP rises ~2300x while ORCHIDEE's rises 2.7x. Median NPP
over 2,488 common cells, g C m-2 yr-1:

    basis                                  CLM45   ORCHIDEE   LPJmL   spread
    raw as published                         586        471     162   3.62x
    PER-TILE (LPJmL / cover)   <-- CHOSEN    586        471     696   1.48x
    per-gridcell (ORCHIDEE x cover)      impossible      36     162   4.49x

Per-tile is both the tightest and the ONLY basis that can include CLM45, which publishes
no cover fraction for any PFT and therefore can never be converted to per-gridcell. The
layer is therefore **NPP per unit conifer-stand area** -- "how productive is a stand
here", not "how much stand is here".

MINIMUM COVER 2%. Dividing by a near-zero cover explodes: unthresholded, LPJmL's per-tile
p99 is 160,203 g C m-2 yr-1 (nonsense); at 2% it is 5,200 and 6,212 cells survive. The
same 2% presence rule is applied to ORCHIDEE (it publishes a fraction, so presence is
knowable) and CLM45 needs none -- its own reporting mask IS its presence mask.

NOTE the consequence for depth: a cell is covered by the models that place a stand there,
so `n_members` and `n_models` vary spatially. CLM45 is NaN outside its tile entirely --
including New Zealand and Chile, where only ORCHIDEE and LPJmL contribute.

UNITS: all three publish `kg m-2 s-1`; emitted as **g C m-2 yr-1** using the 365_day
calendar all three declare (x 1000 g/kg x 31,536,000 s/yr).

REFERENCE SITES (GUARDRAILS §12), 2020s, g C m-2 yr-1, raw as published:

    PNW Oregon (Douglas-fir)   CLM45 366  ORCHIDEE 294  LPJmL 305
    Georgia (loblolly)               815           442        177
    Bavaria (Norway spruce)          600           545         84
    Japan Honshu                     745           538          0 (LPJmL places no stand)
    NZ / Chile (radiata)             NaN       618/285    321/160
    Sahara / Siberia 65N / Greenland NaN             0          0   <- controls behave

    S Sweden reads 0 in all three, correctly: those stands are BOREAL in every one of
    these schemes (Siberia 65N is 0 too). This layer is the temperate class only.

Direction: higher = BETTER (productivity), so the percentile is INVERTED -- low or
declining NPP earns a high risk percentile. No spatial smoothing: the field is a stand
presence mask and a kernel would bleed productivity into cells where no model places a
stand (the same reasoning applied to sugarcane, and declared here too).

Output: data/processed/conifer-temperate_npp-tempnle_annual/npp-tempnle_{scenario}_processed.nc
"""

import argparse
import glob
import multiprocessing as mp
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

LAYER_ID = "conifer-temperate_npp-tempnle_annual"
OUT_VAR = "npp-tempnle"
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2010, 2099
TWO_TIER_ZERO_THRESHOLD = 0.02
HIGHER_IS_BETTER = True
SLOPE_PER_DECADE = 10.0
MIN_COVER = 0.02                  # presence threshold; see the docstring's p99 table
#: kg m-2 s-1 -> g C m-2 yr-1 on the 365_day calendar all three models declare.
UNIT_SCALE = 1000.0 * 365.0 * 86400.0
UNITS = "g C m-2 yr-1"

#: model -> (filename prefix, its own temperate-NLE class code, cover handling)
#: cover: "divide"  = published per grid cell, divide by cover to reach the per-tile basis
#:        "mask"    = already per-tile; cover used only as a presence mask
#:        "none"    = no cover published; the model's own reporting mask is presence
MODELS = {
    "clm45":    ("clm45", "needleleaf-evergreen-tree-temperate", "none"),
    "orchidee": ("orchidee", "tendev", "mask"),
    "lpjml":    ("lpjml", "temperate-needleleaved-evergreen-tree", "divide"),
}
MODEL_FAMILY = {}                 # three structurally distinct models; one vote each

REFERENCE_SITES = [
    ("PNW Oregon (Douglas-fir)", 45.0, -123.5), ("Sierra Nevada (ponderosa)", 38.5, -120.0),
    ("Georgia SE-US (loblolly)", 32.5, -82.5), ("Bavaria (Norway spruce)", 48.5, 11.5),
    ("Japan Honshu (sugi)", 36.0, 138.0), ("NZ North Is. (radiata)", -38.5, 176.0),
    ("Chile Bio-Bio (radiata)", -37.0, -72.5), ("Spain Iberian range", 40.5, -3.0),
    ("CONTROL Amazon", -3.0, -60.0), ("CONTROL Sahara", 23.0, 10.0),
    ("CONTROL Siberia 65N", 65.0, 100.0),
]


def log(msg):
    print(msg, flush=True)


def family_of(model):
    return MODEL_FAMILY.get(model, model)


def parse_name(fpath):
    """(model, gcm, scenario, soc, sens, variable) from an ISIMIP2b biomes filename.

    orchidee_gfdl-esm2m_ewembi_rcp26_2005soc_co2_npp-tendev_global_annual_2006_2099.nc4
      [-8]=scenario [-7]=soc [-6]=sens [-5]=variable [-4]=region [-3]=step [-2:]=years

    Parsed from the END. NOTE awk's `$(NF-4)` is Python's `p[-5]` -- transcribing awk
    offsets directly cost a silent scenario merge once already (WORKFLOW-ISSUES 2026-08-11).
    """
    p = os.path.basename(fpath).replace(".nc4", "").replace(".nc", "").split("_")
    info = dict(model=p[0], gcm=p[1], scenario=p[-8], soc=p[-7], sens=p[-6],
                variable=p[-5], member=f"{p[0]}_{p[1]}")
    if not re.fullmatch(r"(rcp|ssp)\d{2,3}|historical|picontrol", info["scenario"]):
        raise ValueError(f"{os.path.basename(fpath)}: parsed scenario "
                         f"{info['scenario']!r} is not a scenario token -- offsets wrong")
    return info


def years_from_time(ds):
    """Integer years from the CF axis. All three declare 'years since 1661' on a
    365_day calendar, but parse defensively -- a sibling layer encoded the same axis
    in DAYS (GUARDRAILS §9)."""
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
    else:
        dpy = 360.0 if "360" in cal else 365.0 if "365" in cal else 365.25
        yrs = base + vals / dpy
    return np.round(yrs).astype(int)


def resolve_var(ds, var):
    """Find `var` in `ds` tolerating the separator drift ISIMIP2b actually contains.

    Measured 2026-08-12: LPJmL names this variable with HYPHENS in its rcp26/rcp60 files
    (`npp-temperate-needleleaved-evergreen-tree`) and with UNDERSCORES in rcp85
    (`npp_temperate_needleleaved_evergreen_tree`) -- same model, same PFT, same run,
    different separator per scenario. Constructing the name from the filename therefore
    fails on a third of the LPJmL members.
    """
    if var in ds.variables:
        return var
    alt = var.replace("-", "_")
    if alt in ds.variables:
        return alt
    coords = {"lat", "lon", "time", "time_bnds", "lat_bnds", "lon_bnds"}
    cands = [v for v in ds.variables if v not in coords and ds[v].ndim == 3]
    if len(cands) == 1:
        return cands[0]
    raise KeyError(f"cannot resolve {var!r} in {list(ds.variables)}")


def _read(fpath, var):
    ds = xr.open_dataset(fpath, decode_times=False, decode_cf=False)
    da = ds[resolve_var(ds, var)]
    fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
    yrs = years_from_time(ds)
    vals = da.values.astype("float32")
    lats, lons = ds.lat.values, ds.lon.values
    attrs = dict(units=da.attrs.get("units"), long_name=da.attrs.get("long_name"))
    ds.close()
    if fill is not None:
        vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
    vals[~np.isfinite(vals)] = np.nan
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    return yrs[keep], vals[keep], lats, lons, attrs


def load_member(raw_dir, model_key, gcm, scen):
    """Load one member on the PER-TILE basis in g C m-2 yr-1, presence-masked.

    Returns (years, cube, diag) or None when the member does not exist.
    """
    prefix, pft, cover_mode = MODELS[model_key]
    stem = f"{prefix}_{gcm}_ewembi_{scen}_2005soc_co2"
    npp_f = raw_dir / f"{stem}_npp-{pft}_global_annual_2006_2099.nc4"
    if not npp_f.exists():
        return None
    yrs, npp, lats, lons, attrs = _read(npp_f, f"npp-{pft}")
    diag = dict(units_in=attrs["units"], raw_finite=int(np.isfinite(npp).any(axis=0).sum()))

    if cover_mode != "none":
        cov_f = raw_dir / f"{stem}_pft-{pft}_global_annual_2006_2099.nc4"
        if not cov_f.exists():
            raise FileNotFoundError(f"cover fraction missing for {model_key} {gcm} {scen}")
        cyrs, cov, _, _, cattrs = _read(cov_f, f"pft-{pft}")
        if not np.array_equal(cyrs, yrs):
            raise ValueError(f"{model_key} {gcm} {scen}: cover years != npp years")
        # `units` lies for ORCHIDEE ('%' but stored 0-1); decide from the VALUES.
        cmax = np.nanmax(cov)
        as_percent = bool(cmax > 1.5)
        cov = cov / 100.0 if as_percent else cov
        diag.update(cover_declared=cattrs["units"], cover_stored=("percent" if as_percent
                                                                  else "fraction"),
                    cover_max=float(np.nanmax(cov)))
        present = np.isfinite(cov) & (cov >= MIN_COVER)
        if cover_mode == "divide":
            npp = np.where(present, npp / np.clip(cov, MIN_COVER, None), np.nan)
        else:
            npp = np.where(present, npp, np.nan)
    npp = npp * np.float32(UNIT_SCALE)
    diag["present_cells"] = int(np.isfinite(npp).any(axis=0).sum())
    return yrs, npp.astype("float32"), lats, lons, diag


def make_pct_fn(baseline_flat, higher_is_better=HIGHER_IS_BETTER):
    """Percentile-of-score vs the shared 2020s baseline, as a RISK score.

    Higher NPP is better, so the raw score is inverted (101 - raw): low or declining
    productivity earns a HIGH risk percentile. Two-tier only if the baseline is
    materially zero-inflated, which a presence-masked productivity field should not be.
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
    out = np.full(shape, np.nan, np.float32)
    out.ravel()[land_idx] = flat_land
    return out


#: Per-worker Theil-Sen memory budget, as in process_let_cyclone / process_driedarea.
SLOPE_MEM_BUDGET_BYTES = 700 * 1024**2
_CUBE = None            # set in main(); forked workers inherit it copy-on-write


def slope_chunk_cells(n_members, n_years, max_pairs):
    """Chunk width that keeps Theil-Sen peak memory inside the per-worker budget."""
    obs = n_members * n_years
    pairs = obs * (obs - 1) // 2
    if max_pairs is not None:
        pairs = min(pairs, max_pairs)
    per_cell = 4 * pairs * 4
    return max(4, min(512, int(SLOPE_MEM_BUDGET_BYTES // max(per_cell, 1))))


def _slope_block(task):
    """Worker: expanding slopes for one contiguous block of land cells."""
    s, e, years, decade, baseline, window, chunk, max_pairs = task
    res = expanding_slopes(_CUBE[:, :, s:e], years, decade, baseline,
                           window_years=window, chunk_cells=chunk, max_pairs=max_pairs)
    # Return plain arrays, NOT the SlopeResult: it is a dict subclass whose
    # __getattr__ = dict.__getitem__ turns pickle's probe for __getstate__ into a
    # KeyError instead of an AttributeError, which kills the whole pool.
    return s, e, res["ols_slope"], res["sen_slope"]


def compute_slopes(cube, years, decade, chunk, max_pairs, jobs, n_land):
    """expanding_slopes over all cells, fanned across `jobs` forked workers.

    Single-threaded this layer needs ~1 h per scenario: the 2090s window stacks
    10 members x 80 years = 800 observations, i.e. 319,600 exact Theil-Sen pairs per
    cell over 27,377 cells.
    """
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


def site_table(grid, lats, lons, title):
    """GUARDRAILS §12 -- print the layer at named places the subject actually occupies."""
    log(f"\n  {title}")
    for name, la, lo in REFERENCE_SITES:
        i = int(np.abs(lats - la).argmin())
        j = int(np.abs(lons - lo).argmin())
        v = grid[i, j]
        log(f"      {name:28s} {'NaN' if not np.isfinite(v) else f'{v:9.1f}'}")


def main():
    global MIN_COVER, _CUBE
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cover", type=float, default=MIN_COVER)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()
    MIN_COVER = args.min_cover

    root = Path(__file__).resolve().parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    out_dir = root / "data" / "processed" / LAYER_ID
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(raw_dir / "*_npp-*_global_annual_2006_2099.nc4")))
    if not files:
        log(f"ERROR: no npp member files in {raw_dir}")
        return 1
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})

    log("=" * 78)
    log("Processing temperate needleleaf evergreen NPP -> TCFD output contract")
    log("=" * 78)
    log(f"Scenarios: {scenarios}")
    log(f"Basis: PER-TILE (NPP per unit conifer-stand area); min cover {MIN_COVER:.0%}; "
        f"units {UNITS}")
    log("Direction: higher_is_better (risk = productivity loss) -> percentile inverted.")
    log("No spatial smoothing (a stand-presence field must not be blurred).")

    years = np.arange(MIN_YEAR, MAX_YEAR + 1)
    n_years = years.size
    y_index = {y: i for i, y in enumerate(years)}

    lats = lons = None
    raw = {s: {} for s in scenarios}
    log("\nLoading members (harmonizing to the per-tile basis)...")
    for model_key in MODELS:
        for gcm in ("gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"):
            for scen in scenarios:
                got = load_member(raw_dir, model_key, gcm, scen)
                if got is None:
                    continue
                yrs, cube_in, la, lo, diag = got
                if lats is None:
                    lats, lons = la, lo
                cube = np.full((n_years, len(lats), len(lons)), np.nan, np.float32)
                for k, y in enumerate(yrs):
                    if int(y) in y_index:
                        cube[y_index[int(y)]] = cube_in[k]
                member = f"{MODELS[model_key][0]}_{gcm}"
                raw[scen][member] = cube
                cov_txt = ""
                if "cover_stored" in diag:
                    cov_txt = (f" cover declared={diag['cover_declared']!r} stored="
                               f"{diag['cover_stored']} max={diag['cover_max']:.3f}")
                log(f"  {model_key:9s} {gcm:13s} {scen}  raw_cells={diag['raw_finite']:6,d} "
                    f"-> present={diag['present_cells']:6,d}{cov_txt}")

    LAT, LON = len(lats), len(lons)
    for s in scenarios:
        log(f"  {s}: {len(raw[s])} members -> {sorted(raw[s])}")

    probe = next(iter(raw[scenarios[0]].values()))
    boolean = is_boolean_field(probe)
    if boolean:
        log("ERROR: NPP classified boolean -- wrong input.")
        return 1
    stat_name = "pooled_median"
    fin = probe[np.isfinite(probe)]
    log(f"\nField nature: CONTINUOUS ({UNITS}) min={fin.min():.4g} max={fin.max():.4g} "
        f"median={np.median(fin):.4g} negatives={100*np.mean(fin<0):.2f}% "
        f"exact-zero={100*np.mean(fin==0):.2f}% -> {stat_name}")

    finite_any = np.zeros((LAT, LON), bool)
    for s in scenarios:
        for cube in raw[s].values():
            finite_any |= np.isfinite(cube).any(axis=0)
    land_idx = np.flatnonzero(finite_any.ravel())
    n_land = land_idx.size
    log(f"Stand cells (union over members): {n_land:,} of {LAT*LON:,}")

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
    log(f"\nEnsemble uniform across scenarios: {uniform} "
        f"(CLM45 coverage is scenario-dependent by construction)")

    # ---- Shared 2020s baseline: pooled over EVERY (year x member x scenario) --------
    base_pool = np.concatenate([annual[s] for s in scenarios], axis=0)
    b_med, b_lo, b_hi = pooled_decadal_stat(
        base_pool, years, BASELINE_DECADE, boolean=False, window_years=WINDOW_YEARS)
    n_base_members = base_pool.shape[0]
    del base_pool

    baseline_flat = b_med[np.isfinite(b_med)]
    pct, pct_mode, frac_zero = make_pct_fn(baseline_flat)
    b_pct = pct(b_med)
    log(f"Shared 2020s baseline: pooled over {n_base_members} member-scenario series, "
        f"land n={baseline_flat.size:,}, exact-zero={frac_zero:.2%}, mode={pct_mode}")
    log(f"  median NPP over stand cells = {np.median(baseline_flat):.1f} {UNITS}  "
        f"(p05 {np.percentile(baseline_flat,5):.1f}, p95 {np.percentile(baseline_flat,95):.1f})")
    site_table(scatter(b_med, land_idx, (LAT, LON)), lats, lons,
               f"REFERENCE SITES (GUARDRAILS §12) -- 2020s baseline, {UNITS}:")

    # ---- Per-member diagnostic ------------------------------------------------------
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
    xr.Dataset(
        {"value": (["member", "lat", "lon"], mem_grid)},
        coords={"member": mem_names, "lat": lats, "lon": lons},
        attrs={"variable": OUT_VAR, "units": UNITS,
               "member_field": f"mean {BASELINE_DECADE}s NPP on the per-tile basis, "
                               "pooled across scenarios",
               "note": "Diagnostic only -- not part of the OUTPUT-SPEC contract."},
    ).to_netcdf(out_dir / f"{OUT_VAR}_members.nc",
                encoding={"value": {"dtype": "float32", "zlib": True, "complevel": 4,
                                    "_FillValue": np.float32(np.nan)}})
    log(f"\n  wrote per-member diagnostic {OUT_VAR}_members.nc ({len(mem_names)} members)")
    del mem_grid

    models_all = sorted({m.split("_")[0] for m in mem_names})
    gcms_all = sorted({m.split("_", 1)[1] for m in mem_names})

    for s in scenarios:
        log(f"\n{'='*78}\nAssembling scenario {s}\n{'='*78}")
        mem = members_by_scen[s]
        cube = annual[s]
        _CUBE = cube            # forked workers inherit this copy-on-write
        chunk = slope_chunk_cells(len(mem), n_years, args.max_pairs)
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
                pc = pct(med)
            sl = compute_slopes(cube, years, d, chunk, args.max_pairs, args.jobs, n_land)

            # The presence mask is TIME-VARYING (cover crosses the 2% threshold as
            # stands expand or retreat), so a cell can have observations inside the
            # EXPANDING slope window yet none inside the decade window itself -- finite
            # slope, NaN median. That is not an ocean leak but it violates the contract's
            # mask agreement, and a trend over a decade with no stand is meaningless.
            # Mask slopes to the decade's own median mask and record how many cells.
            slope_only = np.isfinite(sl["ols_slope"]) & ~np.isfinite(med)
            if slope_only.any():
                sl = {"ols_slope": np.where(np.isfinite(med), sl["ols_slope"], np.nan),
                      "sen_slope": np.where(np.isfinite(med), sl["sen_slope"], np.nan)}

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
                             ("ols_slope", sl["ols_slope"] * SLOPE_PER_DECADE),
                             ("sen_slope", sl["sen_slope"] * SLOPE_PER_DECADE),
                             ("n_members", n_mem), ("n_models", n_mod)):
                out[key][i] = scatter(vec, land_idx, (LAT, LON))

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                gm = np.nanmean(out["median"][i])
                txt = ("slopes=NaN (baseline)" if d == BASELINE_DECADE else
                       f"ols={np.nanmean(out['ols_slope'][i]):+.3f} "
                       f"sen={np.nanmean(out['sen_slope'][i]):+.3f} {UNITS}/dec")
            tag = "shared baseline" if d == BASELINE_DECADE else f"{len(mem)} members"
            extra_txt = (f"  [{int(slope_only.sum()):,} slope-only cells masked]"
                         if d != BASELINE_DECADE and slope_only.any() else "")
            log(f"  {d}s: {tag:<15} mean={gm:8.2f} {UNITS}  {txt}{extra_txt}")

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                assert np.all(np.isnan(out["ols_slope"][i])), "baseline ols must be NaN"
                assert np.all(np.isnan(out["sen_slope"][i])), "baseline sen must be NaN"
                continue
            med_finite = np.isfinite(out["median"][i])
            for k in ("ols_slope", "sen_slope"):
                extra = np.isfinite(out[k][i]) & ~med_finite
                assert not extra.any(), f"{k} finite where median is NaN at {d}s -- ocean leak"

        ds_out = xr.Dataset(
            {k: (["decade", "lat", "lon"], v) for k, v in out.items()},
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": OUT_VAR,
                "scenario": s,
                "long_name": "Net primary productivity, temperate needleleaf evergreen "
                             "(conifer) stands",
                "units": UNITS,
                "output_spec": "OUTPUT-SPEC.md",
                "decadal_statistic": stat_name,
                "field_nature": "continuous",
                "pft_classes": ("clm45:needleleaf-evergreen-tree-temperate; "
                                "orchidee:tendev; lpjml:temperate-needleleaved-evergreen-tree"),
                "value_note": (
                    f"median = MEDIAN over the pooled (year x member) sample in the decade "
                    f"window, on the PER-TILE basis (NPP per unit conifer-stand area) in "
                    f"{UNITS}. Converted from the published kg m-2 s-1 with the 365_day "
                    f"calendar all three models declare."),
                "denominator_harmonization": (
                    "MEASURED 2026-08-12, not assumed. corr(cover, npp): orchidee +0.279 "
                    "(per-tile), lpjml +0.898 (cover-scaled per grid cell); clm45 publishes "
                    "no pft- fraction for any PFT and reports only on its own tile. LPJmL is "
                    "therefore DIVIDED by its cover fraction to reach the per-tile basis; "
                    "orchidee and clm45 are used as published. Median over 2,488 common "
                    "cells: raw 586/471/162 (3.62x spread), per-tile 586/471/696 (1.48x), "
                    "per-gridcell impossible for clm45 (4.49x for the other two). NOTE "
                    "orchidee's pft- declares units='%' but stores a FRACTION; the scale is "
                    "decided from the values, never the attribute."),
                "presence_mask": (
                    f"cover >= {MIN_COVER:.0%} for models publishing a fraction (orchidee, "
                    f"lpjml); clm45's own reporting mask is its presence mask. Without a "
                    f"threshold the per-tile division explodes -- lpjml p99 reaches 160,203 "
                    f"{UNITS} unthresholded vs 5,200 at 2%. Cells outside every model's "
                    f"presence are NaN, never 0."),
                "ci_definition": (
                    "lower_ci/upper_ci = 25th/75th percentile of the same pooled "
                    "(year x member) sample. NOT floored at zero: negative NPP is physically "
                    "meaningful (respiration exceeding photosynthesis)."),
                "slope_definition": (
                    "ols_slope = least-squares, sen_slope = Theil-Sen, both over an "
                    "EXPANDING window from the start of the 2020s baseline through the end "
                    "of the target decade, stacking every (year, member) observation. "
                    "Baseline panel is NaN. Member coverage is uneven here (clm45 is absent "
                    "from most cells and from whole scenarios), which is the regime where "
                    "ols_slope absorbs level offsets as trend -- read sen_slope, and treat "
                    "disagreement as a robustness warning."),
                "slope_units": f"{UNITS} decade-1",
                "percentile_baseline": (
                    f"{pct_mode}: ranked against the shared 2020s stand distribution, then "
                    "INVERTED (101 - raw) because higher productivity is better -- low or "
                    "declining NPP = high risk percentile."),
                "percentile_zero_fraction": round(frac_zero, 5),
                "percentile_direction": "higher_is_better",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "baseline_pooling": (
                    f"pooled over EVERY (year x member x scenario) observation in the "
                    f"{BASELINE_DECADE}s window ({n_base_members} member-scenario series), so "
                    f"the panel is bit-identical across scenarios DESPITE scenario-dependent "
                    f"composition (user decision 2026-08-12)."),
                "members_by_scenario": ";".join(
                    f"{sc}:{','.join(members_by_scen[sc])}" for sc in scenarios),
                "ensemble_uniform_across_scenarios": str(uniform),
                "ensemble_note": (
                    "CLM45 publishes only 7 files for this PFT across all future scenarios "
                    "and no GCM of its own in all three RCPs: rcp26 hadgem2-es+miroc5, "
                    "rcp60 ipsl-cm5a-lr, rcp85 gfdl-esm2m+hadgem2-es. n_members and "
                    "n_models therefore vary BY SCENARIO and BY CELL."),
                "co2_treatment": "uniform co2 (transient) for every member",
                "soc_scenario": "2005soc",
                "window_years": WINDOW_YEARS,
                "n_members": len(mem),
                "impact_models": ",".join(models_all),
                "gcms": ",".join(gcms_all),
                "normalization": "none beyond the denominator harmonization described above",
                "spatial_smoothing": (
                    "none -- DECLARED. The field is a stand-presence mask; a 5x5 kernel "
                    "would bleed productivity into cells where no model places a stand."),
                "reference_sites_verified": (
                    "GUARDRAILS §12, 2020s raw: PNW Oregon 366/294/305, Georgia 815/442/177, "
                    "Bavaria 600/545/84, Japan 745/538/0, NZ -/618/321, Chile -/285/160 "
                    "(clm45/orchidee/lpjml). Controls Sahara, Siberia 65N, Greenland all 0. "
                    "S Sweden is 0 in all three CORRECTLY -- those stands are boreal in "
                    "every one of these schemes; this layer is the temperate class only."),
                "source_dataset": "ISIMIP2b OutputData/biomes (CLM45, ORCHIDEE, LPJmL)",
                "description": (
                    "Temperate needleleaf evergreen (conifer) stand productivity processed "
                    "to the TCFD output contract with a shared 2020s baseline; 3 models x "
                    "CMIP5 GCMs on a measured per-tile denominator, no smoothing, "
                    "higher_is_better (risk = productivity loss)."),
            },
        )
        encoding = {k: {"dtype": "float32", "zlib": True, "complevel": 4,
                        "_FillValue": np.float32(np.nan)} for k in out}
        path = out_dir / f"{OUT_VAR}_{s}_processed.nc"
        ds_out.to_netcdf(path, encoding=encoding)
        log(f"  saved {path.name}  ({path.stat().st_size/(1024*1024):.1f} MB)")

    log(f"\nDone -> {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
