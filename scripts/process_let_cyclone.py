"""Process let (tropical-cyclone exposure) into the TCFD 6-value-class format.

let (Lange 2020, ISIMIP2b DerivedOutputData, model ke-tg-meanfield) is the
tropical-cyclone member of the extreme-event exposure family. UNLIKE led
(drought), which is a binary {0,1} per-cell flag, let is a CONTINUOUS FRACTION
in [0,1): each cell holds the fraction of its land area exposed to tropical
cyclones that year (verified 2026-07-24 -- ~97-98% exact 0, remainder spread
smoothly over (0,1), never exactly 1). The decadal statistic is still the MEAN
(the expected annual exposed-area fraction); the median is useless here because
most cell-years are 0.

Thin ensemble: 1 impact model (ke-tg-meanfield) x 4 CMIP5 GCMs = 4 members,
vs led's 31. To make per-cell decadal estimates robust we apply SPATIAL
SMOOTHING to each per-member decadal map before ensemble aggregation: a 5x5
window, weight w = exp(-d / L) with L = 0.7 grid-cells (user-chosen "sharp"
decay), distance cos(lat)-scaled in longitude so the footprint stays ~isotropic
in km through the subtropics, normalized over the non-NaN (land) neighbours only,
longitude-wrapped at the antimeridian, latitude clamped at grid edges. Ocean
cells stay NaN (never filled).

Everything downstream derives from the smoothed inputs (median, lower/upper_ci,
percentile, shared 2020s baseline, within-decade trend) so the product is
internally consistent.

Percentile deviation from led: let's 2020s baseline is zero-inflated (~75% exact
0 after smoothing), so a single percentile-of-score would crush the real spatial
gradient. Instead we use a TWO-TIER percentile: cells with 0 exposure -> the 1st
percentile; cells with > 0 exposure -> ranked only against the NON-ZERO cells of
the 2020s baseline and mapped linearly to [2, 100], so the exposed coasts/tropics
use the full dynamic range. (led used a plain side="right" score, fine there
because its 31-member baseline is continuous.)

Output: let_{scenario}_processed.nc with variables
{median, percentile, trend, lower_ci, upper_ci} on (decade, lat, lon).
"""

import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402

VAR = "let"
LAYER_ID = "cyclone_let_annual"
RAW_PATTERN = "*_2006_2099.nc4"
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2010, 2099

# spatial-smoothing kernel
SMOOTH_L = 0.7      # exponential length scale (grid cells); "sharp" decay
SMOOTH_RADIUS = 2   # 5x5 window


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (impact_model, gcm, scenario, member) from a lange2020 filename."""
    p = os.path.basename(fpath).split("_")
    return dict(model=p[1], gcm=p[2], scenario=p[4], member=f"{p[1]}_{p[2]}")


def years_from_time(ds):
    """Parse integer calendar years from a 'years since YYYY-...' time axis."""
    units = ds.time.attrs.get("units", "")
    m = re.search(r"years since (\d+)", units)
    base = int(m.group(1)) if m else 0
    return (base + ds.time.values).astype(int)


def load_member(fpath):
    """Load one member file as (year, lat, lon), restricted to MIN..MAX_YEAR."""
    ds = xr.open_dataset(fpath, decode_times=False)
    yrs = years_from_time(ds)
    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).load()
    ds.close()
    return da


def decade_mean_map(da, decade):
    """Mean over a decade window = per-member exposed-fraction map (lat, lon)."""
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) == 0:
        return None
    return da.isel(year=sel).mean("year").values.astype(np.float32)


# ----------------------------- spatial smoothing ---------------------------- #
def make_kernel(lats, L=SMOOTH_L, radius=SMOOTH_RADIUS):
    """Precompute per-row weight vectors for each (di, dj) offset in the window.

    Distance is cos(lat)-scaled in longitude: d = sqrt(di^2 + (dj*cos(lat))^2),
    so at 35N a longitude step counts as ~0.8 of a latitude step (cells narrow
    toward the poles). Returns a list of (di, dj, w) with w shape (nlat,)."""
    coslat = np.cos(np.deg2rad(lats)).astype(np.float64)
    offs = []
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            d = np.sqrt(di * di + (dj * coslat) ** 2)
            offs.append((di, dj, np.exp(-d / L)))
    return offs


def _shift(a, di, dj):
    """Return a[i+di, j+dj] with longitude wrap and latitude-edge rows zeroed."""
    nlat, nlon = a.shape
    i_src = np.arange(nlat) + di
    lat_ok = (i_src >= 0) & (i_src < nlat)
    i_src = np.clip(i_src, 0, nlat - 1)
    j_src = (np.arange(nlon) + dj) % nlon
    out = a[i_src][:, j_src]
    if not lat_ok.all():
        out = out.copy()
        out[~lat_ok, :] = 0.0
    return out


def smooth_field(field, kernel):
    """5x5 exponential-decay smooth over non-NaN neighbours; ocean stays NaN."""
    valid = np.isfinite(field)
    f0 = np.where(valid, field, 0.0).astype(np.float64)
    vf = valid.astype(np.float64)
    wsum = np.zeros(field.shape, np.float64)
    wtot = np.zeros(field.shape, np.float64)
    for di, dj, w in kernel:
        wcol = w[:, None]
        wsum += wcol * _shift(f0, di, dj)
        wtot += wcol * _shift(vf, di, dj)
    out = np.full(field.shape, np.nan, np.float64)
    good = wtot > 0
    out[good] = wsum[good] / wtot[good]
    out[~valid] = np.nan  # never fill ocean
    return out.astype(np.float32)


def slope_over_years(annual_stack, years):
    """OLS slope per cell across `years` for a (nyear, lat, lon) stack."""
    if len(years) < 3:
        return np.full(annual_stack.shape[1:], np.nan, np.float32)
    t = years - years.mean()
    var_t = float(np.mean(t ** 2))
    if var_t <= 0:
        return np.full(annual_stack.shape[1:], np.nan, np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yc = annual_stack - np.nanmean(annual_stack, axis=0, keepdims=True)
        cov = np.nanmean(t[:, None, None] * yc, axis=0)
    return (cov / var_t).astype(np.float32)


def make_pct_fn(baseline_flat):
    """Two-tier percentile for a zero-inflated hazard (higher = worse):

    * value == 0  -> percentile 1  (never-exposed land)
    * value  > 0  -> ranked ONLY against the non-zero cells of the 2020s
                     baseline, linearly mapped to [2, 100] (exposed land uses
                     the full 2-100 range; the median non-zero baseline cell
                     sits near 51).

    Scoring the positives against the full (~75%-zero) baseline would crush the
    real gradient into the top quantiles; restricting to non-zero cells restores
    contrast among the coasts/tropics that actually see cyclones."""
    nz_sort = np.sort(baseline_flat[baseline_flat > 0])
    n_nz = len(nz_sort)

    def pct(arr):
        flat = arr.ravel()
        out = np.full(flat.shape, np.nan, np.float32)
        fin = np.isfinite(flat)
        vals = flat[fin]
        res = np.ones(vals.shape, np.float32)  # zero-exposure cells -> 1
        pos = vals > 0
        if n_nz > 0:
            frac = np.searchsorted(nz_sort, vals[pos], side="left") / n_nz
            res[pos] = np.clip(2.0 + 98.0 * frac, 2.0, 100.0)
        out[fin] = res
        return out.reshape(arr.shape)

    return pct


def main():
    files = [str(p) for p in storage.stage_raw(LAYER_ID, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no let member files found in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
        log("Upload members with storage.ingest_raw() first.")
        return
    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})
    log("=" * 64)
    log("Processing let (tropical-cyclone exposure) -> TCFD 6-value-class format")
    log("=" * 64)
    log(f"Members: {len(files)} | scenarios: {scenarios}")
    log(f"Smoothing: 5x5, w=exp(-d/L) L={SMOOTH_L}, cos(lat)-scaled, non-NaN-normalized")

    # grid geometry + kernel
    da0 = load_member(files[0])
    lats, lons = da0.lat.values, da0.lon.values
    LAT, LON = len(lats), len(lons)
    del da0
    kernel = make_kernel(lats)

    # ---- Pass 1: per-member SMOOTHED decadal maps + raw annual ensemble ------
    dec_sm = {s: {} for s in scenarios}  # dec_sm[scen][member] = (n_dec, lat, lon) smoothed
    years_all = np.arange(MIN_YEAR, MAX_YEAR + 1)
    yidx = {y: i for i, y in enumerate(years_all)}
    ens_sum = {s: np.zeros((len(years_all), LAT, LON), np.float64) for s in scenarios}
    ens_cnt = {s: np.zeros((len(years_all), LAT, LON), np.float64) for s in scenarios}

    raw_gm, sm_gm = [], []  # 2020s global-mean before/after smoothing (diagnostic)
    for f in files:
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da = load_member(f)
        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        for i, d in enumerate(DECADES):
            raw = decade_mean_map(da, d)
            if raw is not None:
                sm = smooth_field(raw, kernel)
                maps[i] = sm
                if d == BASELINE_DECADE:
                    raw_gm.append(np.nanmean(raw))
                    sm_gm.append(np.nanmean(sm))
        dec_sm[s][member] = maps
        # raw annual accumulation for the trend
        v = da.values.astype(np.float64)
        fin = np.isfinite(v)
        rows = np.array([yidx[y] for y in da.year.values])
        ens_sum[s][rows] += np.where(fin, v, 0.0)
        ens_cnt[s][rows] += fin
        log(f"  loaded+smoothed {info['gcm']:<13} {s}")
        del da, v, fin

    log(f"\nSmoothing check @2020s: raw global-mean={np.mean(raw_gm):.5f} "
        f"-> smoothed={np.mean(sm_gm):.5f} (should be ~equal; smoothing conserves mass)")

    # raw ensemble annual, then smooth each annual map (for the trend)
    ens_annual_sm = {}
    for s in scenarios:
        with np.errstate(invalid="ignore", divide="ignore"):
            ea = ens_sum[s] / ens_cnt[s]
        ea[ens_cnt[s] == 0] = np.nan
        sm = np.stack([smooth_field(ea[k], kernel) for k in range(ea.shape[0])], 0)
        ens_annual_sm[s] = sm.astype(np.float32)

    # ---- Shared 2020s baseline (from smoothed maps) ---------------------------
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec_sm[s]})
    shared_2020 = []
    for member in all_members:
        per_scen = [dec_sm[s][member][b_idx] for s in scenarios if member in dec_sm[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_2020.append(np.nanmean(np.stack(per_scen, 0), axis=0))
    shared_2020 = np.stack(shared_2020, 0).astype(np.float32)  # (member, lat, lon)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared_median = np.nanmean(shared_2020, axis=0)
        shared_sd = np.nanstd(shared_2020, axis=0)
    shared_lo = np.clip(shared_median - shared_sd, 0, 1)
    shared_hi = np.clip(shared_median + shared_sd, 0, 1)

    baseline_flat = shared_median[np.isfinite(shared_median)].ravel()
    pct = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)
    frac_zero = float(np.mean(baseline_flat == 0.0))
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"land cells n={len(baseline_flat):,}, exact-zero fraction={frac_zero:.2%}, "
        f"global-mean exposed fraction={np.nanmean(shared_median):.5f}")

    # shared 2020s within-decade trend (mean of both scenarios' smoothed annual)
    win = np.arange(BASELINE_DECADE, BASELINE_DECADE + WINDOW_YEARS)
    wrows = np.array([yidx[y] for y in win])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # all-ocean cells -> empty nanmean slice
        pooled = np.nanmean(np.stack([ens_annual_sm[s][wrows] for s in scenarios], 0), axis=0)
    shared_trend = slope_over_years(pooled.astype(np.float32), win)

    # ---- Per-scenario assembly + write ----------------------------------------
    for s in scenarios:
        log(f"\n{'='*64}\nAssembling scenario {s}\n{'='*64}")
        members = sorted(dec_sm[s])
        stack = np.stack([dec_sm[s][m] for m in members], 0)  # (member, dec, lat, lon)

        median = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        percentile = np.full_like(median, np.nan)
        trend = np.full_like(median, np.nan)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)

        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                median[i], percentile[i] = shared_median, shared_pct
                lower[i], upper[i], trend[i] = shared_lo, shared_hi, shared_trend
                log(f"  {d}s: shared baseline (identical across scenarios)")
                continue
            layer = stack[:, i, :, :]  # (member, lat, lon) smoothed
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median[i] = np.nanmean(layer, axis=0)
                sd = np.nanstd(layer, axis=0)
            lower[i] = np.clip(median[i] - sd, 0, 1)
            upper[i] = np.clip(median[i] + sd, 0, 1)
            percentile[i] = pct(median[i])
            win_d = np.arange(d, d + WINDOW_YEARS)
            rows = np.array([yidx[y] for y in win_d if y in yidx])
            trend[i] = slope_over_years(ens_annual_sm[s][rows], win_d[: len(rows)])
            log(f"  {d}s: {layer.shape[0]} members  "
                f"global-mean exposed fraction={np.nanmean(median[i]):.5f}")

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": "Land area fraction exposed to tropical cyclones",
                "units": "1",
                "statistic": "decadal_mean_exposed_area_fraction_spatially_smoothed",
                "value_note": ("median = ensemble mean over 4 GCMs of the annual "
                               "fraction of each cell's land area exposed to tropical "
                               "cyclones; raw ISIMIP values are fractional in [0,1) "
                               "(NOT a binary flag)."),
                "spatial_smoothing": (f"5x5 window, w=exp(-d/L) L={SMOOTH_L} grid-cells "
                                      "(sharp decay), cos(lat)-scaled longitude distance, "
                                      "normalized over non-NaN land neighbours, longitude-"
                                      "wrapped; applied to each per-member decadal map "
                                      "before ensemble aggregation."),
                "smoothing_rationale": ("thin ensemble (1 impact model ke-tg-meanfield x "
                                        "4 GCMs = 4 members); spatial pooling borrows "
                                        "strength from neighbouring cells."),
                "ci_definition": ("lower/upper_ci = ensemble mean -/+ 1 inter-GCM standard "
                                  "deviation of the smoothed per-member decadal maps, "
                                  "clamped to [0,1]."),
                "percentile_baseline": ("two-tier: value 0 -> percentile 1; value > 0 -> "
                                        "ranked against the NON-ZERO cells of the 2020s "
                                        "ensemble-mean baseline, linearly mapped to [2,100] "
                                        "(exposed land uses the full 2-100 range)."),
                "percentile_direction": "higher_is_worse",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "ensemble_note": "single impact model; ensemble spread is inter-GCM only",
                "impact_models": ",".join(sorted({meta[f]["model"] for f in files})),
                "gcms": ",".join(sorted({meta[f]["gcm"] for f in files})),
                "source_dataset": "lange2020 ISIMIP2b DerivedOutputData (ke-tg-meanfield)",
                "description": ("Tropical-cyclone exposure processed to TCFD 6-value-class "
                                "format with shared 2020s baseline and 5x5 spatial "
                                "smoothing to offset the thin 4-member ensemble."),
            },
        )
        path = out_dir / f"let_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID,
        stage,
        created_by="scripts/process_let_cyclone.py",
        notes="See WORKFLOW-ISSUES.md 2026-07-24; GUARDRAILS §8-§10.",
    )

    # Every ingest-and-process run leaves reviewable HTML evidence behind.
    # Never raises: the data is already published and gated, and these
    # artifacts are regenerable (see scripts/utils/finalize.py).
    log("\nGenerating QA/QC report and maps...")
    finalize_layer(LAYER_ID, version=version)

    log("\nDone.")


if __name__ == "__main__":
    main()
