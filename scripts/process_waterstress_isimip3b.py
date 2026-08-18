"""Water stress from ISIMIP3b WaterGAP2-2e — four metrics on the TCFD contract.

Demand over available water. Four metrics, each annual and worst-month:

    L  local   =  W(c)     / Q_runoff(c)      local self-sufficiency (a DIAGNOSTIC)
    R  routed  =  W_hat(c) / Q_river(c)       water stress (the reportable metric)

`W_hat` is demand accumulated over every cell draining through `c`, itself included, on
the DDM30 D8 network. Numerator and denominator then cover the same catchment, which is
the whole reason routing is here — a cell-local numerator over a routed denominator is
not withdrawal-to-availability in any sense the literature or the 40/80 bands support.

MODEL SCOPE: WaterGAP2-2e alone, 5 GCMs x 3 scenarios. Provisional (user decision
2026-08-17) and n_models = 1, so the interval carries NO structural uncertainty. See
docs/water-stress-status-2026-08-17.md §6 for the measurement of what that omits.

THE ARITHMETIC, all of it measured rather than assumed (§11):
  * flux -> volume uses cellarea [km2] AND contfrac/100. Reversed once; see §11.2-11.3.
  * `contfrac` is in PERCENT despite a `units` attribute of `1`.
  * calendar is `365_day` for every WaterGAP file (H08 differs -- read it, never assume).
  * `_FillValue` 1e20 survives `np.asarray()` of a masked array and is FINITE, so it
    passes `isfinite`. Mask |x| >= 1e19 explicitly before any arithmetic.

Supply is 1850soc (naturalised: no modern abstraction, no reservoirs), demand is 2015soc.
Pairing is member-wise. Dividing demand by post-abstraction discharge would put demand on
top and also subtract it from the bottom.

Usage:
    python scripts/process_waterstress_isimip3b.py --members 1     # one member, gates
    python scripts/process_waterstress_isimip3b.py                 # all 15
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import netCDF4

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.decadal_stats import pooled_decadal_stat, expanding_slopes, is_boolean_field  # noqa: E402

RAW = Path("data/raw/water_stress")
ROUTING = Path("data/raw/water_stress/routing")
FILL = 1e19
MODEL = "watergap2-2e"
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
YEARS = np.arange(2015, 2101)
#: 365_day calendar -- no leap years, ever. Verified on all 45 monthly files.
MONTH_SEC = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], float) * 86400.0
#: ARBITRARY, PRACTICAL, AND DECLARED (user decision 2026-08-18).
#:
#: The ratio is capped at 100. This is NOT a claim that demand cannot exceed 100x local
#: availability -- it demonstrably can, and 1000x is physically reasonable in a desert
#: settlement drawing on fossil groundwater or an aqueduct. The cap exists because the
#: ratio stops carrying information long before it stops being computable: SDG 6.4.2's
#: top class is ">100%" and Aqueduct's is ">80%", so a cell at 6.5e8 and a cell at 12 are
#: both "demand vastly exceeds renewable supply" and the extra orders of magnitude are an
#: artifact of dividing by a number that rounds to nothing.
#:
#: MEASURED before choosing it: the pathology is 0.016% of routed cell-years in 174 cells
#: (Sahara, hyperarid, genuine trickle baseline flow going to ~0 in a given year). Every
#: real >1 value survives -- 4.45% of routed cell-years. Three successive masks failed to
#: fix this because each time the cause moved; the cap is the honest treatment.
#:
#: CONSEQUENCE, which must reach the file attributes: this makes a CENSORED FIELD. Per
#: OUTPUT-SPEC, where a layer is pinned at a bound BOTH slope estimators go to ~0 and
#: AGREE, so the usual "disagreement means the trend is fragile" rule gives no warning.
#: The share of cells at the cap is measured per decade and declared, as heatwave-3b does.
RATIO_CAP = 100.0

#: DDM30 D8, for ASCENDING latitude. Unique cycle-free mapping of 16 candidates.
OFF = {1: (0, 1), 2: (-1, 1), 3: (-1, 0), 4: (-1, -1),
       5: (0, -1), 6: (1, -1), 7: (1, 0), 8: (1, 1)}


def read(path: Path, var: str, tslice=None):
    """Read a variable with fills normalised to NaN. See module docstring on why the
    explicit magnitude mask is required rather than trusting the mask machinery."""
    with netCDF4.Dataset(path) as nc:
        v = nc.variables[var]
        v.set_auto_mask(False)
        a = np.asarray(v[tslice] if tslice is not None else v[:], np.float64)
        lat = np.asarray(nc.variables["lat"][:], np.float64)
        lon = np.asarray(nc.variables["lon"][:], np.float64)
    return np.where(np.abs(a) >= FILL, np.nan, a), lat, lon


def to_ddm_grid(x, slat, slon, dlat, dlon):
    """Subset a model field onto the DDM30 grid by EXACT coordinate match.

    DDM30 latitude ascends (-55.75 -> 83.75) over 280 rows; model output descends
    (89.75 -> -89.75) over 360. Assuming either orientation silently flips the world.
    """
    ri = {round(float(v), 3): i for i, v in enumerate(slat)}
    ci = {round(float(v), 3): i for i, v in enumerate(slon)}
    rows = [ri[round(float(v), 3)] for v in dlat]
    cols = [ci[round(float(v), 3)] for v in dlon]
    return x[..., rows, :][..., :, cols]


def build_network(fd):
    """(next-cell index, topological order) for the DDM30 land cells."""
    nlat, nlon = fd.shape
    land = np.isfinite(fd)
    ii, jj = np.nonzero(land)
    codes = fd[ii, jj].astype(int)
    nxt = np.full(nlat * nlon, -1, np.int64)
    for c, (di, dj) in OFF.items():
        s = codes == c
        if not s.any():
            continue
        ni, nj = ii[s] + di, (jj[s] + dj) % nlon
        ok = (ni >= 0) & (ni < nlat)
        nic = ni.clip(0, nlat - 1)
        ok &= land[nic, nj]
        nxt[ii[s] * nlon + jj[s]] = np.where(ok, nic * nlon + nj, -1)
    nodes = ii * nlon + jj
    indeg = np.zeros(nlat * nlon, np.int32)
    t = nxt[nodes]
    np.add.at(indeg, t[t >= 0], 1)
    order, front = [], nodes[indeg[nodes] == 0]
    while front.size:
        order.append(front)
        tg = nxt[front]
        m = tg >= 0
        u, cn = np.unique(tg[m], return_counts=True)
        indeg[u] -= cn
        front = u[indeg[u] == 0]
    order = np.concatenate(order)
    if order.size != nodes.size:
        raise RuntimeError(f"DDM30 network not acyclic: {nodes.size - order.size} cells")
    return nxt, order, nodes, land


def accumulate(payload_flat, nxt, order):
    """Push every cell's payload downstream. payload_flat is (ncell_flat, ntime)."""
    acc = payload_flat.copy()
    for k in order:
        t = nxt[k]
        if t >= 0:
            acc[t] += acc[k]
    return acc


def member_river_mask(di, land, threshold):
    """River mask for ONE member, from that member's own 2020s baseline discharge.

    MEASURED 2026-08-18: a pooled (ensemble-mean) mask is WRONG. A cell at 70.75N
    154.25W carries a pooled baseline of 0.31 m3/s and passes a 0.1 m3/s cut, while
    THIS member's own discharge there is 3e-14 m3/s -- thirteen orders of magnitude
    smaller. The cell then divides real accumulated demand by essentially nothing and
    returns a ratio of 6e9 in every year. 686 cells pass the pooled mask while their own
    member baseline is below 0.1 m3/s; 188 are below 0.001.

    One member's water must not license another member's ratio. The mask stays fixed in
    time (2020s baseline, per the shared-baseline rule) but is now per member.
    """
    base = np.nanmean(di[60:180], axis=0)          # 2020-2029 monthly mean, m3/s
    return np.isfinite(base) & (base > threshold) & land


def member_annual(gcm, scenario, ddm, diag):
    """Four annual metric series for one member, with the precedence table applied.

    Returns dict of metric -> (year, ncell) on the DDM30 land axis.
    """
    dlat, dlon, land, nxt, order, nodes, area_eff, river_threshold = ddm
    nlat, nlon = land.shape
    ii, jj = np.nonzero(land)

    def load(var, soc):
        f = RAW / f"{MODEL}_{gcm}_w5e5_{scenario}_{soc}_default_{var}_global_monthly_2015_2100.nc"
        x, slat, slon = read(f, var)
        return to_ddm_grid(x, slat, slon, dlat, dlon)          # (1032, nlat, nlon)

    pt = load("ptotww", "2015soc")
    qt = load("qtot", "1850soc")
    di = load("dis", "1850soc")

    river = member_river_mask(di, land, river_threshold)
    diag.setdefault("river_cells_by_member", {})[f"{gcm}|{scenario}"] = int(river.sum())

    nt = pt.shape[0]
    dt = np.tile(MONTH_SEC, nt // 12)[:, None, None]

    # --- rule 3: negatives are physically impossible and must surface, not hide inside
    #     the zero or mask rules that follow.
    for nm, arr in (("ptotww", pt), ("qtot", qt), ("dis", di)):
        n_neg = int(np.nansum(arr < 0))
        diag[f"negative_{nm}"] = diag.get(f"negative_{nm}", 0) + n_neg
        if n_neg:
            arr[arr < 0] = np.nan

    W = pt * area_eff * 1000.0 * dt          # m3 per month, per cell
    Qr = qt * area_eff * 1000.0 * dt         # m3 per month, local runoff
    Qv = di * dt                             # m3 per month, routed discharge

    # --- routed numerator: accumulate demand downstream
    flat = np.zeros((nlat * nlon, nt))
    w2 = np.where(np.isfinite(W), W, 0.0)
    flat[nodes] = w2[:, ii, jj].T
    What = np.full((nlat * nlon, nt), np.nan)
    What[nodes] = accumulate(flat[nodes] if False else flat, nxt, order)[nodes]
    What = What.reshape(nlat, nlon, nt).transpose(2, 0, 1)

    def monthly_ratio(num, den, routed):
        r = np.full(num.shape, np.nan)
        valid = np.isfinite(num) & np.isfinite(den) & land[None, :, :]
        if routed:
            valid &= river[None, :, :]                 # rule 4, BEFORE the zero rules
        zero_num = valid & (num == 0)
        r[zero_num] = 0.0                              # rule 5: no demand is not stress
        pos = valid & (num > 0)
        bad = pos & (den == 0)                         # rule 6: undefined, counted
        diag["zero_denom_positive_num"] = diag.get("zero_denom_positive_num", 0) + int(bad.sum())
        good = pos & (den > 0)
        r[good] = num[good] / den[good]
        return r

    out = {}
    for key, num, den, routed in (("L", W, Qr, False), ("R", What, Qv, True)):
        mr = monthly_ratio(num, den, routed)
        mr3 = mr.reshape(-1, 12, nlat, nlon)
        complete = np.isfinite(mr3).all(axis=1)        # rule 8: 12 valid months or NaN
        with np.errstate(invalid="ignore"):
            num_a = np.where(np.isfinite(num), num, 0.0).reshape(-1, 12, nlat, nlon).sum(1)
            den_a = np.where(np.isfinite(den), den, 0.0).reshape(-1, 12, nlat, nlon).sum(1)
            ann = np.where(complete & (den_a > 0), num_a / np.where(den_a > 0, den_a, 1), np.nan)
            ann = np.where(complete & (num_a == 0), 0.0, ann)
            mx = np.where(complete, np.nanmax(np.where(np.isfinite(mr3), mr3, -np.inf), axis=1), np.nan)
        mx = np.where(np.isfinite(mx), mx, np.nan)
        for nm, arr in (("ann", ann), ("max", mx)):
            hit = np.isfinite(arr) & (arr > RATIO_CAP)
            diag[f"capped_{key}_{nm}"] = diag.get(f"capped_{key}_{nm}", 0) + int(hit.sum())
            diag[f"total_{key}_{nm}"] = diag.get(f"total_{key}_{nm}", 0) + int(np.isfinite(arr).sum())
            arr[hit] = RATIO_CAP
        out[f"{key}_ann"] = ann[:, ii, jj].astype(np.float32)
        out[f"{key}_max"] = mx[:, ii, jj].astype(np.float32)
        # Annual denominator, retained so "supply collapse" can be defined against each
        # cell's OWN climatology rather than an absolute cut that would penalise a small
        # river for being small.

    for k, v in out.items():                            # rule 10: postcondition
        if np.isinf(v).any():
            raise RuntimeError(f"{k}: {int(np.isinf(v).sum())} infinite values — "
                               f"precedence table failed")
    return out


def load_ddm(river_threshold_m3s):
    fd, dlat, dlon = read(ROUTING / "ddm30_flowdir_cru_neva.nc", "flowdirection")
    nxt, order, nodes, land = build_network(fd)
    ca, slat, slon = read(RAW / f"{MODEL}_gfdl-esm4_w5e5_ssp370_2015soc_default_cellarea_global_fixed_0000_0000.nc", "cellarea")
    cf, _, _ = read(RAW / f"{MODEL}_gfdl-esm4_w5e5_ssp370_2015soc_default_contfrac_global_fixed_0000_0000.nc", "contfrac")
    area_eff = to_ddm_grid(ca, slat, slon, dlat, dlon) * (to_ddm_grid(cf, slat, slon, dlat, dlon) / 100.0)

    # NOTE: the river mask is built PER MEMBER inside member_annual, not here. A pooled
    # mask lets one member's water justify another member's ratio -- see
    # member_river_mask() for the measurement that forced this.
    return (dlat, dlon, land, nxt, order, nodes, area_eff, river_threshold_m3s), None


DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE = 2020
METRICS = ["L_ann", "L_max", "R_ann", "R_max"]
LONG = {
    "L_ann": "Local water stress, annual (demand / locally generated runoff)",
    "L_max": "Local water stress, highest-stress month",
    "R_ann": "Water stress, annual (upstream-accumulated demand / river discharge)",
    "R_max": "Water stress, highest-stress month",
}


def write_contract(metric, scenario, ann_by_member, shared_base, land, dlat, dlon,
                   diag, out_dir):
    """Emit one OUTPUT-SPEC file: (decade, lat, lon) with the full variable set."""
    import xarray as xr

    A = np.stack(ann_by_member)                       # (member, year, cell)
    nlat, nlon = land.shape
    ii, jj = np.nonzero(land)
    shape = (len(DECADES), nlat, nlon)
    out = {k: np.full(shape, np.nan, np.float32)
           for k in ("median", "lower_ci", "upper_ci", "percentile",
                     "ols_slope", "sen_slope", "n_members", "n_models")}

    # Baseline reference for the percentile: this METRIC's own shared 2020s panel, never
    # another metric's -- the four are different estimands with different masks.
    ref = shared_base[np.isfinite(shared_base)]
    ref_nonzero = ref[ref > 0]
    zero_share = float((ref == 0).mean()) if ref.size else 0.0
    two_tier = zero_share > 0.02

    for di, dec in enumerate(DECADES):
        if dec == BASELINE:
            stat, lo, hi = shared_base, None, None
            statf = shared_base
        else:
            statf, lo, hi = pooled_decadal_stat(A, YEARS, dec, boolean=False)
        grid = np.full((nlat, nlon), np.nan, np.float32)
        grid[ii, jj] = statf
        out["median"][di] = grid
        if lo is not None:
            for nm, v in (("lower_ci", lo), ("upper_ci", hi)):
                g = np.full((nlat, nlon), np.nan, np.float32); g[ii, jj] = v
                out[nm][di] = g
        # percentile of score against the shared 2020s baseline distribution
        pct = np.full(statf.shape, np.nan)
        fin = np.isfinite(statf)
        if two_tier:
            pct[fin & (statf == 0)] = 1.0
            nz = fin & (statf > 0)
            if ref_nonzero.size:
                pct[nz] = 1.0 + 99.0 * np.searchsorted(np.sort(ref_nonzero), statf[nz]) / ref_nonzero.size
        else:
            if ref.size:
                pct[fin] = 100.0 * np.searchsorted(np.sort(ref), statf[fin]) / ref.size
        g = np.full((nlat, nlon), np.nan, np.float32); g[ii, jj] = pct
        out["percentile"][di] = g

        sl = expanding_slopes(A, YEARS, dec, BASELINE)
        for nm, v in (("ols_slope", sl["ols_slope"]), ("sen_slope", sl["sen_slope"])):
            g = np.full((nlat, nlon), np.nan, np.float32); g[ii, jj] = v
            out[nm][di] = g
        # per-cell, per-decade member coverage -- NOT a constant; the precedence table
        # reduces coverage unevenly and the contract requires the real count.
        w = (YEARS >= dec) & (YEARS < dec + 10)
        cnt = np.isfinite(A[:, w, :]).any(axis=1).sum(axis=0)
        g = np.full((nlat, nlon), np.nan, np.float32); g[ii, jj] = cnt
        out["n_members"][di] = g
        out["n_models"][di] = np.where(np.isfinite(g), 1.0, np.nan)

    cap_hits = diag.get(f"capped_{metric.split('_')[0]}_{metric.split('_')[1]}", 0)
    cap_tot = max(diag.get(f"total_{metric.split('_')[0]}_{metric.split('_')[1]}", 1), 1)
    ds = xr.Dataset(
        {k: (("decade", "lat", "lon"), v) for k, v in out.items()},
        coords={"decade": DECADES, "lat": dlat, "lon": dlon},
        attrs={
            "variable": metric,
            "long_name": LONG[metric],
            "units": "1 (dimensionless ratio)",
            "slope_units": "1 per year",
            "decadal_statistic": "pooled_median",
            "percentile_direction": "higher_is_worse",
            "relative_baseline": "false",
            "impact_models": MODEL,
            "n_members": len(ann_by_member),
            "n_models": 1,
            "scenario": scenario,
            "ratio_cap": RATIO_CAP,
            "censoring_note": (
                f"Ratio capped at {RATIO_CAP:g}. ARBITRARY and PRACTICAL, not a claim "
                f"that demand cannot exceed {RATIO_CAP:g}x supply -- it can, and 1000x is "
                f"physically reasonable. {100*cap_hits/cap_tot:.4f}% of member-years hit "
                f"the cap. THIS IS A CENSORED FIELD: where cells sit at the bound both "
                f"slope estimators go to ~0 and AGREE, so slope agreement is not evidence "
                f"of a robust trend there."),
            "single_model_note": (
                "One impact model. The interval reflects GCM and interannual spread only "
                "and carries NO structural uncertainty. A three-model comparison measured "
                "68.4% ordering reversal and 5-30x arid-basin spread; see "
                "docs/water-stress-status-2026-08-17.md section 6."),
            "supply_soc": "1850soc (naturalised: negligible abstraction, no modern dams)",
            "demand_soc": "2015soc",
            "calendar_note": "365_day; annual weighting uses no-leap month lengths",
            "area_conversion": "flux * cellarea[km2] * (contfrac/100) * 1000",
            "routing": "DDM30 D8; routed metrics use upstream-accumulated demand",
            "percentile_basis": f"this metric's own shared 2020s panel; two_tier={two_tier}",
        },
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{metric}_{scenario}_processed.nc"
    ds.to_netcdf(path)
    return path


def _write_all(results, land, dlat, dlon, diag) -> int:
    """Shared 2020s baseline then 12 contract files.

    Baseline order matters (OUTPUT-SPEC): average EACH MEMBER across scenarios first,
    then pool. Pooling first is a different number.
    """
    out_dir = Path("data/processed/waterstress-isimip3b")
    w = (YEARS >= BASELINE) & (YEARS < BASELINE + 10)
    for metric in METRICS:
        per_gcm = []
        for g in GCMS:
            stack = np.stack([results[(g, s)][metric][w] for s in SCENARIOS])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                per_gcm.append(np.nanmean(stack, axis=0))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            shared = np.nanmedian(np.concatenate(per_gcm, axis=0), axis=0)
        for s in SCENARIOS:
            ann = [results[(g, s)][metric] for g in GCMS]
            path = write_contract(metric, s, ann, shared, land, dlat, dlon, diag, out_dir)
            print(f"  wrote {path.name}")
    print(f"\nwrote 12 contract files to {out_dir}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-cache", action="store_true",
                    help="skip the member loop and rebuild contract files from the cache")
    ap.add_argument("--no-cache", action="store_true", help="do not write the cache")
    ap.add_argument("--members", type=int, default=0,
                    help="process only the first N members (0 = all 15)")
    ap.add_argument("--river-threshold", type=float, default=0.1,
                    help="baseline mean discharge [m3/s] above which a cell is a river. "
                         "0.1 m3/s is NOMINAL flow at 0.5 deg (user decision 2026-08-18) "
                         "and keeps 96.6%% of land. A large mask was considered and "
                         "REJECTED: 50 m3/s is the Thames at London and would delete the "
                         "Mediterranean, Australian, SW-US and African dryland rivers -- "
                         "precisely the basins this layer exists to find. Measurement "
                         "showed the extreme ratios do NOT come from tiny rivers (cells "
                         "exceeding 100 have a median baseline of 5.0 m3/s) but from "
                         "year-specific denominator collapse, handled separately.")
    args = ap.parse_args()

    if args.from_cache:
        z = np.load("data/processed/_waterstress_annual_cache.npz")
        land, dlat, dlon = z["land"], z["dlat"], z["dlon"]
        results, diag = {}, {}
        for k in z.files:
            if "|" not in k:
                continue
            g, s, m = k.split("|")
            results.setdefault((g, s), {})[m] = z[k]
        print(f"restored {len(results)} members from cache")
        return _write_all(results, land, dlat, dlon, diag)

    print("loading DDM30 network and building the fixed river mask...", flush=True)
    ddm, _ = load_ddm(args.river_threshold)
    dlat, dlon, land, nxt, order, nodes, area_eff, _thr = ddm
    print(f"  land {int(land.sum()):,}; river mask is PER MEMBER at >{args.river_threshold} m3/s")

    pairs = [(g, s) for s in SCENARIOS for g in GCMS]
    if args.members:
        pairs = pairs[: args.members]
    diag: dict = {}
    results = {}
    for i, (g, s) in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] {g} {s} ...", flush=True)
        results[(g, s)] = member_annual(g, s, ddm, diag)

    print("\n--- diagnostics ---")
    for k, v in sorted(diag.items()):
        if isinstance(v, dict):
            n = np.array(list(v.values()))
            print(f"  {k}: min {n.min():,} max {n.max():,} spread {n.max()-n.min():,}")
        else:
            print(f"  {k}: {v:,}")
    for m in METRICS:
        hits = diag.get(f"capped_{m.split('_')[0]}_{m.split('_')[1]}", 0)
        tot = max(diag.get(f"total_{m.split('_')[0]}_{m.split('_')[1]}", 1), 1)
        print(f"  CAP {m}: {hits:,}/{tot:,} = {100*hits/tot:.4f}% at {RATIO_CAP:g}")

    if args.members:
        print("\npartial run -- contract files not written")
        return 0

    # Cache before the writer runs. The member loop is ~10 minutes; a bug downstream of
    # it should cost seconds to retry, not the whole thing again.
    cache = Path("data/processed/_waterstress_annual_cache.npz")
    if not args.no_cache:
        np.savez(cache, **{f"{g}|{s}|{k}": v for (g, s), d in results.items()
                           for k, v in d.items()},
                 land=land, dlat=dlat, dlon=dlon)
        print(f"cached annual arrays -> {cache}")

    return _write_all(results, land, dlat, dlon, diag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
