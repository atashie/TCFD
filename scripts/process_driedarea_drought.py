"""Process ISIMIP3b `driedarea` (drought exposure) into the TCFD 6-value-class format.

`driedarea` (Heinicke2026, ISIMIP3b DerivedOutputData) is the SSP-era sibling of the
ISIMIP2b Lange2020 `led` member, and shares its data nature: a BINARY per-cell annual
flag, 1 = the cell was exposed to drought that year, 0 = not. Underneath, a cell is
exposed in a month when root-zone soil moisture falls below the 2.5th percentile of
the preindustrial baseline for at least seven consecutive months; the annual value is
the maximum of the monthly flag (`cell_methods: time: maximum`), i.e. "at any point
this year". The "fraction of area exposed" emerges only on aggregation -- across
years, members, or space. It is NOT a within-cell area share.

Verified across all 45 members before processing: n_unique == 2, range [0,1], no
intermediate values, calendar `proleptic_gregorian`, time units `days since 1601-01-01`.

Consequences handled here (differ from the continuous-variable reference, process_qg.py):
  * Decadal statistic is the MEAN over years and members (= drought frequency), NOT
    the median -- a median of binary flags collapses to 0/1. This is the standing rule
    for the whole exposure family.
  * Ensemble member = (impact_model, GCM), parsed from filename fields [0] and [1].
    NOTE the offset differs from lange2020, whose files are prefixed with the
    publication name and put the model at [1] / the GCM at [2]. Reusing the `led`
    parser here would mis-key every member.
  * Confidence bounds are the ensemble mean +/- 1 inter-member SD, clamped to [0,1] --
    the raw values are binary, so the mean is a probability-like quantity and bounds
    outside [0,1] are meaningless. This is an inter-member spread envelope, not a
    sampling confidence interval: members share GCMs and GHMs, and single-member cells
    collapse to zero width. Use n_members / n_models to judge a cell.

DELIBERATE CHOICES (do not "fix" without re-reading this block)
---------------------------------------------------------------
UNION MASK, NOTHING FILTERED. The three GHMs do not share a land mask -- h08 53,713
cells, jules-w2 57,523, watergap2-2e 56,410; union 63,455, intersection 46,024. Every
cell any member covers is retained, and per-cell `n_members` / `n_models` are emitted
so thin coverage is auditable downstream rather than silently dropped here.

SINGLE-TIER PERCENTILE, fixed rather than auto-selected. The burntarea processor picks
its tier from an exact-zero threshold; that rule is deliberately NOT applied here. On
the shared 2020s baseline the exact-zero mass is 3.59% over the union but only 0.18%
over fully-covered cells -- i.e. the zeros are an artefact of unequal model coverage
(a single-model cell rests on 30 samples, a fully covered one on 450), not evidence
that drought is rare. Auto-switching to two-tier on that basis would collapse every
never-dry cell onto percentile 1 for the wrong reason. The measured zero fraction is
recorded in the output attrs so the choice stays auditable.

BASELINE-ANCHORED TREND, not a within-decade slope. A decadal frequency built from
binary flags is far too noisy year-to-year for an OLS slope over ten annual points;
the anchored rate is spatially coherent by construction and consistent with the change
map. Matches burntarea and csoil. Unlike those two, the baseline decade's trend is NaN
off-land rather than a finite zero everywhere (they use `np.zeros`, which makes the
whole ocean a finite zero that QA does not catch).

NO SPATIAL SMOOTHING. 15 members per scenario; `let` needed 5x5 smoothing at 4.

Output files: driedarea_{scenario}_processed.nc with variables
{median, percentile, trend, lower_ci, upper_ci, n_members, n_models} on (decade, lat, lon).
"""

import hashlib
import os
import sys
import warnings
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402
from utils.layer_publish import publish_processed_layer  # noqa: E402
from utils.contact_sheet import render_contact_sheet  # noqa: E402
from utils.finalize import finalize_layer  # noqa: E402

VAR = "driedarea"
LAYER_ID = "drought_driedarea_annual"
RAW_PATTERN = "*_driedarea_global_annual_2015_2100.nc"
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
BASELINE_DECADE = 2020
WINDOW_YEARS = 10
MIN_YEAR, MAX_YEAR = 2020, 2099

EXPECTED_MODELS = {"h08", "jules-w2", "watergap2-2e"}
EXPECTED_GCMS = {"gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"}
EXPECTED_PER_SCENARIO = len(EXPECTED_MODELS) * len(EXPECTED_GCMS)  # 15

# Reconstructing the origin URL keeps `inputs.files` complete in layer.json, which is
# what later lets storage.cleanup_raw() delete raw safely.
SOURCE_BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/Heinicke2026"
MODEL_DIRS = {"h08": "H08", "jules-w2": "JULES-W2", "watergap2-2e": "WaterGAP2-2e"}


def log(msg):
    print(msg, flush=True)


def parse_name(fpath):
    """Extract (model, gcm, scenario, member) from a Heinicke2026 filename.

    {model}_{gcm}_w5e5_{ssp}_2015soc_default_driedarea_global_annual_2015_2100.nc
    Fields [0] and [1] -- NOT [1] and [2] as in lange2020.
    """
    p = os.path.basename(fpath).split("_")
    return dict(model=p[0], gcm=p[1], scenario=p[3], member=f"{p[0]}_{p[1]}")


def source_url(name):
    info = parse_name(name)
    return f"{SOURCE_BASE}/{MODEL_DIRS[info['model']]}/{info['gcm']}/future/{name}"


def load_member(fpath):
    """Load one member as (year, lat, lon) over MIN..MAX_YEAR, checking data nature.

    Decodes time with the file's own declared calendar rather than trusting xarray's
    default -- per-member calendar divergence has bitten this project before.
    """
    ds = xr.open_dataset(fpath, decode_times=False)
    units = ds.time.attrs["units"]
    calendar = ds.time.attrs.get("calendar", "standard")
    yrs = np.array([t.year for t in cftime.num2date(ds.time.values, units, calendar)])

    da = ds[VAR].assign_coords(year=("time", yrs)).swap_dims({"time": "year"})
    keep = np.where((yrs >= MIN_YEAR) & (yrs <= MAX_YEAR))[0]
    da = da.isel(year=keep).transpose("year", "lat", "lon").load()
    lats, lons = ds.lat.values.copy(), ds.lon.values.copy()
    ds.close()

    finite = da.values[np.isfinite(da.values)]
    n_unique = np.unique(finite).size
    return da, lats, lons, calendar, n_unique


def decade_freq_map(da, decade):
    """Mean over a full decade window = per-member drought frequency map (lat, lon).

    Requires the complete window: a partial decade would be averaged as if it were a
    shorter record, silently mixing sample sizes across the map.
    """
    yrs = da.year.values
    sel = np.where((yrs >= decade) & (yrs <= decade + WINDOW_YEARS - 1))[0]
    if len(sel) < WINDOW_YEARS:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return da.isel(year=sel).mean("year").values.astype(np.float32)


def anchored_trend(med_stack, i, b_idx, baseline_finite):
    """Baseline-anchored trend: (median[i] - median[baseline]) / elapsed decades.

    Units are value per decade. The baseline decade itself has no elapsed change, so
    it is 0 -- but only where the baseline median exists; off-land it stays NaN, so
    the trend's finite mask matches the median's.
    """
    span = i - b_idx
    if span == 0:
        out = np.zeros(med_stack.shape[1:], np.float32)
        out[~baseline_finite] = np.nan
        return out
    return ((med_stack[i] - med_stack[b_idx]) / span).astype(np.float32)


def make_pct_fn(baseline_flat):
    """Single-tier percentile-of-score against the 2020s baseline (higher = worse).

    Fixed single-tier by decision -- see the DELIBERATE CHOICES block. Returns the
    scorer plus the measured exact-zero fraction, which is recorded in the output
    attrs so the choice can be re-examined from the published file alone.
    """
    frac_zero = float(np.mean(baseline_flat == 0.0)) if baseline_flat.size else float("nan")
    bsort = np.sort(baseline_flat)
    n = len(bsort)

    def pct(arr):
        flat = arr.ravel()
        out = np.full(flat.shape, np.nan, np.float32)
        fin = np.isfinite(flat)
        if n > 0:
            frac = np.searchsorted(bsort, flat[fin], side="right") / n
            out[fin] = np.clip(100.0 * frac, 1.0, 100.0).astype(np.float32)
        return out.reshape(arr.shape)

    return pct, frac_zero


def main():
    files = [str(p) for p in storage.stage_raw(LAYER_ID, RAW_PATTERN)]
    if not files:
        log(f"ERROR: no driedarea member files found in "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
        log("Ingest members with scripts/download_driedarea_drought.py first.")
        return

    stage = storage.staging_dir(LAYER_ID)
    out_dir = stage / "data"
    meta = {f: parse_name(f) for f in files}
    scenarios = sorted({m["scenario"] for m in meta.values()})

    log("=" * 70)
    log("Processing driedarea (ISIMIP3b drought exposure) -> TCFD 6-value-class")
    log("=" * 70)
    log(f"Members: {len(files)} | scenarios discovered: {scenarios}")

    # ---- Membership assertions -------------------------------------------------
    # parse_name does no validation, and dec[s][member] = maps would silently
    # overwrite a duplicate, so verify the matrix is exactly what we expect.
    seen = {}
    for f in files:
        i = meta[f]
        key = (i["scenario"], i["member"])
        if key in seen:
            raise SystemExit(f"duplicate member file for {key}:\n  {seen[key]}\n  {f}")
        seen[key] = f
        if i["model"] not in EXPECTED_MODELS:
            raise SystemExit(f"unexpected impact model {i['model']!r} in {f}")
        if i["gcm"] not in EXPECTED_GCMS:
            raise SystemExit(f"unexpected GCM {i['gcm']!r} in {f}")
    for s in scenarios:
        n = sum(1 for (sc, _) in seen if sc == s)
        if n != EXPECTED_PER_SCENARIO:
            raise SystemExit(f"{s}: {n} members, expected {EXPECTED_PER_SCENARIO}")
    log(f"  membership OK: {len(scenarios)} scenarios x {EXPECTED_PER_SCENARIO} members")

    # ---- Pass 1: load members, build per-member decadal frequency maps ----------
    da0, lats, lons, _, _ = load_member(files[0])
    LAT, LON = len(lats), len(lons)
    del da0

    dec = {s: {} for s in scenarios}
    calendars, uniques = set(), set()

    for f in files:
        info = meta[f]
        s, member = info["scenario"], info["member"]
        da, mlats, mlons, calendar, n_unique = load_member(f)

        # Grid identity: values are written into a shared positional array, so a
        # flipped or shifted axis would corrupt a member invisibly.
        if not (np.array_equal(mlats, lats) and np.array_equal(mlons, lons)):
            raise SystemExit(f"grid mismatch in {os.path.basename(f)}")
        calendars.add(calendar)
        uniques.add(n_unique)

        maps = np.full((len(DECADES), LAT, LON), np.nan, np.float32)
        incomplete = []
        for i, d in enumerate(DECADES):
            fm = decade_freq_map(da, d)
            if fm is None:
                incomplete.append(d)
            else:
                maps[i] = fm
        if incomplete:
            raise SystemExit(f"{os.path.basename(f)}: incomplete decades {incomplete}")
        dec[s][member] = maps

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gm = float(np.nanmean(maps[DECADES.index(BASELINE_DECADE)]))
        log(f"  loaded {info['model']:<13} {info['gcm']:<14} {s}  "
            f"years={da.year.size:<3d} n_unique={n_unique} 2020s_mean={gm:.4f}")
        del da

    if uniques != {2}:
        log(f"  WARNING: not every member is binary -- n_unique values seen: {uniques}")
    else:
        log(f"\nAll {len(files)} members binary {{0,1}}; calendars: {sorted(calendars)}")

    # ---- Shared 2020s baseline -------------------------------------------------
    # Per member, average its 2020s map across scenarios, then ensemble-mean. By
    # 2020-2029 the SSPs have barely diverged in forcing: the per-member spread
    # across scenarios (mean |max-min| 0.104) sits BELOW the pure-sampling-noise
    # floor for a 10-year mean of a Bernoulli field (0.139), and the scenario global
    # means agree to 0.005. So this averages out weather noise rather than smearing
    # forced divergence, and it makes the baseline bit-identical across scenarios.
    b_idx = DECADES.index(BASELINE_DECADE)
    all_members = sorted({m for s in scenarios for m in dec[s]})
    shared_2020 = []
    for member in all_members:
        per_scen = [dec[s][member][b_idx] for s in scenarios if member in dec[s]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shared_2020.append(np.nanmean(np.stack(per_scen, 0), axis=0))
    shared_2020 = np.stack(shared_2020, 0).astype(np.float32)  # (member, lat, lon)

    # ---- Per-member contact sheet (GUARDRAILS S11) ------------------------------
    # Rendered here because this is the only point at which each member's own field
    # still exists separately; once pooled, a defective member is diluted and
    # invisible. The three GHM land masks differ visibly between panels.
    sheet_path = None
    try:
        sheet_path = render_contact_sheet(
            {m: shared_2020[i] for i, m in enumerate(all_members)},
            stage / "contact_sheet.html", LAYER_ID, BASELINE_DECADE, units="1",
            note=(f"{len(all_members)} members (3 GHMs x 5 GCMs), 2020s drought "
                  f"frequency. The three GHMs do NOT share a land mask (h08 53,713 / "
                  f"jules-w2 57,523 / watergap2-2e 56,410 cells; union 63,455, "
                  f"intersection 46,024) -- differing coastlines between panels are "
                  f"expected, not a defect. Nothing is masked or filtered."))
        log(f"  contact sheet: {sheet_path}  <-- REVIEW THIS before trusting the layer")
    except Exception as e:  # noqa: BLE001
        log(f"  WARNING: contact sheet failed ({type(e).__name__}: {e}); "
            f"the per-member visual check of GUARDRAILS S11 has NOT been produced")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared_median = np.nanmean(shared_2020, axis=0)
        shared_sd = np.nanstd(shared_2020, axis=0)
    shared_lo = np.clip(shared_median - shared_sd, 0, 1)
    shared_hi = np.clip(shared_median + shared_sd, 0, 1)
    baseline_finite = np.isfinite(shared_median)

    # Per-cell coverage. The three GHMs disagree on the land mask, so record how many
    # members and how many DISTINCT models back each cell. Union is retained; these
    # counts are what make a thin cell auditable instead of merely looking confident.
    member_models = np.array([m.split("_")[0] for m in all_members])
    shared_nmem = np.sum(np.isfinite(shared_2020), axis=0).astype(np.float32)
    shared_nmodel = np.zeros_like(shared_nmem)
    for model in sorted(set(member_models)):
        rows = np.where(member_models == model)[0]
        shared_nmodel += np.isfinite(shared_2020[rows]).any(axis=0)
    shared_nmem[~baseline_finite] = np.nan
    shared_nmodel[~baseline_finite] = np.nan

    baseline_flat = shared_median[baseline_finite].ravel()
    pct, frac_zero = make_pct_fn(baseline_flat)
    shared_pct = pct(shared_median)

    n_full = int((shared_nmodel == len(EXPECTED_MODELS)).sum())
    log(f"\nShared 2020s baseline: {shared_2020.shape[0]} members, "
        f"{len(baseline_flat):,} land cells (union), global-mean freq="
        f"{np.nanmean(shared_median):.4f}")
    log(f"  exact-zero fraction {frac_zero:.4f} ({100 * frac_zero:.2f}%); "
        f"single-tier percentile applied by decision (not auto-selected)")
    log(f"  coverage: {n_full:,} cells with all 3 models, "
        f"{len(baseline_flat) - n_full:,} partially covered (retained)")

    # ---- Per-scenario assembly + write -----------------------------------------
    for s in scenarios:
        log(f"\n{'=' * 70}\nAssembling scenario {s}\n{'=' * 70}")
        members = sorted(dec[s])
        stack = np.stack([dec[s][m] for m in members], 0)  # (member, dec, lat, lon)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            med_stack = np.nanmean(stack, axis=0).astype(np.float32)
            sd_stack = np.nanstd(stack, axis=0).astype(np.float32)
        med_stack[b_idx] = shared_median  # shared baseline, identical across scenarios

        median = med_stack
        percentile = np.full_like(median, np.nan)
        trend = np.full_like(median, np.nan)
        lower = np.full_like(median, np.nan)
        upper = np.full_like(median, np.nan)
        nmem = np.full_like(median, np.nan)
        nmodel = np.full_like(median, np.nan)

        scen_models = np.array([m.split("_")[0] for m in members])
        for i, d in enumerate(DECADES):
            if d == BASELINE_DECADE:
                percentile[i] = shared_pct
                lower[i], upper[i] = shared_lo, shared_hi
                nmem[i], nmodel[i] = shared_nmem, shared_nmodel
                trend[i] = anchored_trend(med_stack, i, b_idx, baseline_finite)
                log(f"  {d}s: shared baseline (identical across scenarios)")
                continue

            lower[i] = np.clip(median[i] - sd_stack[i], 0, 1)
            upper[i] = np.clip(median[i] + sd_stack[i], 0, 1)
            percentile[i] = pct(median[i])
            trend[i] = anchored_trend(med_stack, i, b_idx, baseline_finite)

            finite_i = np.isfinite(median[i])
            cnt = np.sum(np.isfinite(stack[:, i]), axis=0).astype(np.float32)
            mdl = np.zeros_like(cnt)
            for model in sorted(set(scen_models)):
                rows = np.where(scen_models == model)[0]
                mdl += np.isfinite(stack[rows, i]).any(axis=0)
            cnt[~finite_i] = np.nan
            mdl[~finite_i] = np.nan
            nmem[i], nmodel[i] = cnt, mdl

            log(f"  {d}s: {len(members)} members  "
                f"global-mean freq={np.nanmean(median[i]):.4f}  "
                f"change vs 2020s={np.nanmean(median[i]) - np.nanmean(shared_median):+.4f}")

        ds_out = xr.Dataset(
            {
                "median": (["decade", "lat", "lon"], median),
                "percentile": (["decade", "lat", "lon"], percentile),
                "trend": (["decade", "lat", "lon"], trend),
                "lower_ci": (["decade", "lat", "lon"], lower),
                "upper_ci": (["decade", "lat", "lon"], upper),
                "n_members": (["decade", "lat", "lon"], nmem),
                "n_models": (["decade", "lat", "lon"], nmodel),
            },
            coords={"decade": DECADES, "lat": lats, "lon": lons},
            attrs={
                "variable": VAR,
                "scenario": s,
                "long_name": "Annual frequency of exposure to drought",
                "units": "1",
                "statistic": "decadal_mean_exposure_frequency",
                "value_note": (
                    "median = ensemble-mean fraction of model-years in which the cell "
                    "was flagged as drought-exposed, built from BINARY per-cell annual "
                    "flags (1 = exposed at some point that year). This is an occurrence "
                    "FREQUENCY, not a within-cell area share: area affected is a "
                    "downstream derivation (cell area x frequency), and treating it as "
                    "such assumes the whole cell is affected when flagged."),
                "hazard_definition": (
                    "Root-zone soil moisture below the 2.5th percentile of the "
                    "preindustrial baseline for at least 7 consecutive months "
                    "(Lange et al. 2020, Earth's Future); annual value is the maximum "
                    "of the monthly flag."),
                "ci_definition": (
                    "lower/upper_ci = ensemble mean -/+ 1 inter-member standard "
                    "deviation, clamped to [0,1]. An inter-member SPREAD ENVELOPE, not "
                    "a sampling confidence interval: members share GCMs and GHMs, and "
                    "single-member cells collapse to zero width. Read with n_members."),
                "percentile_baseline": (
                    "percentile-of-score of the decadal ensemble-mean against the 2020s "
                    "ensemble-mean spatial distribution over the union land mask "
                    "(2020s centres near 50)"),
                "percentile_tier": "single_tier",
                "percentile_tier_note": (
                    f"Single-tier applied BY DECISION, not auto-selected from a "
                    f"threshold. Measured exact-zero fraction of the 2020s baseline: "
                    f"{frac_zero:.4f} over the union mask. The zero mass is "
                    f"concentrated in partially covered cells (a single-model cell "
                    f"rests on 30 samples vs 450 for a fully covered one), so it "
                    f"reflects uneven model coverage rather than drought being rare."),
                "percentile_direction": "higher_is_worse",
                "trend_definition": (
                    "baseline-anchored rate: (median[decade] - median[2020s]) / elapsed "
                    "decades, i.e. proportional to the change map. NOT a within-decade "
                    "slope -- a decadal frequency built from binary flags is too noisy "
                    "year-to-year. Baseline decade = 0 where the baseline exists, NaN "
                    "off-land."),
                "trend_units": "1 decade-1",
                "baseline_decade": BASELINE_DECADE,
                "baseline_source": "shared_across_all_scenarios",
                "baseline_note": (
                    "ISIMIP3b begins in 2015, so the layer starts at the 2020s "
                    "baseline; there is no full 2010s decade. Per member the 2020s map "
                    "is averaged across all three SSPs before pooling: by 2020-2029 the "
                    "cross-SSP spread is below the pure-sampling-noise floor of a "
                    "10-year Bernoulli mean, so this suppresses weather noise rather "
                    "than smearing forced divergence."),
                "spatial_mask": (
                    "UNION of the three GHM land masks -- h08 53,713 / jules-w2 57,523 "
                    "/ watergap2-2e 56,410 cells; union 63,455, intersection 46,024. "
                    "Partially covered cells are RETAINED, not masked; per-cell "
                    "n_members and n_models are emitted so they can be filtered "
                    "downstream. Nothing is zero-filled."),
                "smoothing": "none (15 members per scenario)",
                "normalization": "none (binary flags share a common definition)",
                "window_years": WINDOW_YEARS,
                "n_members": len(members),
                "n_models": len(set(scen_models)),
                "impact_models": ",".join(sorted(set(scen_models))),
                "gcms": ",".join(sorted({meta[f]["gcm"] for f in files})),
                "soc_sens": "2015soc/default for all 3 models (uniform)",
                "source_dataset": (
                    "ISIMIP3b DerivedOutputData/Heinicke2026 (driedarea); the SSP-era "
                    "sibling of ISIMIP2b Lange2020 `led`"),
                "description": (
                    "Drought exposure processed to TCFD 6-value-class format with a "
                    "shared 2020s baseline; 3 GHMs x 5 CMIP6 GCMs = 15 members per "
                    "scenario, ssp126/370/585. Union land mask, no filtering, no "
                    "smoothing, single-tier percentile, baseline-anchored trend."),
            },
        )
        path = out_dir / f"driedarea_{s}_processed.nc"
        ds_out.to_netcdf(path)
        log(f"  staged {path.name}")

    # ---- Publish ---------------------------------------------------------------
    # raw_entries must be recorded or storage.cleanup_raw() will later refuse to
    # delete raw staging for want of a source_url per input.
    raw_entries = []
    for f in sorted(files):
        name = os.path.basename(f)
        raw_entries.append({
            "name": name,
            "bytes": os.path.getsize(f),
            "sha256": hashlib.sha256(Path(f).read_bytes()).hexdigest(),
            "source_url": source_url(name),
        })

    log("\nPublishing to S3...")
    version = publish_processed_layer(
        LAYER_ID,
        stage,
        created_by="scripts/process_driedarea_drought.py",
        raw_entries=raw_entries,
        notes=("ISIMIP3b Heinicke2026 driedarea, ssp126/370/585; 15 members/scenario "
               "(h08, jules-w2, watergap2-2e x 5 CMIP6 GCMs). Binary exposure flags -> "
               "decadal MEAN frequency. Union land mask, nothing filtered; single-tier "
               "percentile by decision; baseline-anchored trend; no smoothing."),
    )

    # Every ingest-and-process run leaves reviewable HTML evidence behind.
    # Never raises: the data is already published and gated, and these
    # artifacts are regenerable (see scripts/utils/finalize.py).
    log("\nGenerating QA/QC report and maps...")
    finalize_layer(LAYER_ID, version=version,
                   extra_maps=[sheet_path] if sheet_path else None)

    log("\nDone.")


if __name__ == "__main__":
    main()
