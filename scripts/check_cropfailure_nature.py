"""Measure the data nature and cropland coverage of raw `cropfailure` before processing.

Two INDEPENDENT guardrail checks, run together because they need the same open files but
answer different questions:

  * GUARDRAILS 9 -- WHAT the values are. Data nature is taken from the values, never from
    the variable name, the CF `long_name`, a sibling variable, or another layer. This
    decides the decadal-statistic branch (pooled_median / pooled_mean_boolean /
    pooled_mean_zero_inflated), the CI definition and the percentile tiering. Within the
    Lange 2020 family `led` is binary and `let` is a continuous fraction despite identical
    naming; within Heinicke2026, `floodedarea` is binary while `driedarea` is not and the
    two also differ on masking. Nothing here is inheritable.

  * GUARDRAILS 12 -- WHERE the values are. A layer about crops must be non-trivial where
    crops demonstrably grow. The withdrawn ISIMIP2b sugarcane layers passed every
    contract check while reading exactly 0 across the entire cane belt, and the 87%
    zero-fraction was measured, classified as structural, and written up as fine. Counting
    zeros does not locate them.

Reads every member in data/raw/, so it is also the mask survey the processor needs: the
2.8x file-size spread across models (promet 3.13 MB -> isam 8.74 MB at the same GCM and
scenario) predicts heterogeneous coverage, and the minimum-model rule must be measured
here rather than inherited from `led` (>= 2) or `driedarea` (full union).

Usage:
    python scripts/check_cropfailure_nature.py [--limit N] [--sites-only]
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr

VAR = "cropfailure"
LAYER_ID = "cropfailure-isimip3b_cropfailure_annual"

#: GUARDRAILS 12 reference sites: major crop-producing regions where field crops
#: demonstrably grow, spread across continents, hemispheres and crop systems. A NaN here
#: means the layer does not cover cropland; an all-zero-across-every-decade-and-scenario
#: reading means it carries no signal there. Either is a STOP.
REFERENCE_SITES = [
    ("Iowa, US Corn Belt",          42.0,  -93.5),
    ("Mato Grosso, BR soy",        -13.0,  -56.0),
    ("Pampas, AR",                 -34.0,  -61.5),
    ("Punjab, IN wheat/rice",       30.5,   75.5),
    ("North China Plain",           35.0,  114.5),
    ("Ukraine steppe wheat",        49.0,   32.0),
    ("Beauce, FR wheat",            48.3,    1.8),
    ("WA wheatbelt, AU",           -31.5,  117.5),
    ("Guinea savanna, NG",          10.5,    8.0),
    ("Ethiopian highlands",          9.0,   38.5),
    ("Sahel millet, ML",            14.0,   -5.0),
    ("Java rice, ID",               -7.0,  110.0),
]


def parse_name(fpath):
    """(model, gcm, scenario, member) read from the END of the filename.

    Zantout2025 carries a LEADING publication token (`zantout2025_`) where its 3b sibling
    Heinicke2026 does not, so a forward index reads the wrong column. From the end it is
    invariant under either convention.
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
    """Years for each record, WITHOUT trusting the time axis to declare itself.

    Zantout2025's `time` coordinate carries `long_name`/`standard_name`/`axis` and NO
    `units` attribute, so xarray cannot decode it (`.dt` raises) and there is no epoch to
    read. The values are a bare contiguous integer sequence -- 165..250 for a file whose
    name declares 2015_2100 -- i.e. `years since 1850`, undeclared.

    Rather than hardcode 1850, the span is taken from the filename (which IS part of the
    ISIMIP grammar) and the axis is checked to be consistent with it: right length,
    contiguous unit steps. The implied epoch is returned so a member that disagrees with
    the rest is visible instead of silently mis-binned by a decade.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    y0, y1 = int(p[-2]), int(p[-1])
    t = np.asarray(ds["time"].values, dtype="float64")
    n = y1 - y0 + 1
    if t.size != n:
        raise ValueError(f"{os.path.basename(fpath)}: {t.size} records but the filename "
                         f"declares {y0}-{y1} ({n} years)")
    steps = np.unique(np.diff(t))
    if steps.size != 1 or not np.isclose(steps[0], 1.0):
        raise ValueError(f"{os.path.basename(fpath)}: time axis is not a contiguous "
                         f"annual sequence (steps {steps})")
    return np.arange(y0, y1 + 1, dtype=int), int(round(y0 - t[0]))


def site_indices(lats, lons, lat, lon):
    """Nearest grid cell, handling either a -180..180 or a 0..360 longitude axis."""
    lon_q = lon if lons.min() < 0 else (lon % 360.0)
    return int(np.abs(lats - lat).argmin()), int(np.abs(lons - lon_q).argmin())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="check only the first N members (iteration only -- the shipped "
                         "check must be exhaustive, a sampled scan can step over a rare "
                         "continuous tail and silently switch the published statistic)")
    ap.add_argument("--sites-only", action="store_true")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc")))
    if not files:
        print(f"ERROR: no {VAR} files in {raw_dir}")
        return 2
    if args.limit:
        files = files[:args.limit]

    print("=" * 78)
    print(f"cropfailure raw value check -- {len(files)} members")
    print("=" * 78)

    with xr.open_dataset(files[0]) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    print(f"grid: lat {lats.min():.2f}..{lats.max():.2f} ({len(lats)}), "
          f"lon {lons.min():.2f}..{lons.max():.2f} ({len(lons)})")
    sites = [(name, *site_indices(lats, lons, la, lo), la, lo)
             for name, la, lo in REFERENCE_SITES]

    # ---- per-member scan ------------------------------------------------------ #
    rows = []
    attr_sets = defaultdict(set)
    model_mask = {}
    nonzero_mask = {}
    site_series = defaultdict(dict)     # site -> member -> mean over all years
    global_min, global_max = np.inf, -np.inf
    any_non_binary = False

    for i, f in enumerate(files, 1):
        info = parse_name(f)
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[VAR]
            attr_sets["units"].add(str(da.attrs.get("units", "<absent>")))
            attr_sets["long_name"].add(str(da.attrs.get("long_name", "<absent>")))
            attr_sets["time_units"].add(str(ds.time.attrs.get("units", "<absent>")))
            attr_sets["calendar"].add(str(ds.time.attrs.get("calendar", "<absent>")))
            yrs, epoch = decode_years(ds, f)
            attr_sets["implied_time_epoch"].add(str(epoch))
            vals = da.values.astype("float32")
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))

        if fill is not None:
            vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
        vals[~np.isfinite(vals)] = np.nan

        fin = vals[np.isfinite(vals)]
        if fin.size == 0:
            print(f"  [{i}/{len(files)}] {info['member']:<28} {info['scenario']}  "
                  "ALL-NaN MEMBER -- STOP")
            return 3

        n_uniq = np.unique(fin).size
        exact0 = float((fin == 0.0).mean())
        exact1 = float((fin == 1.0).mean())
        binary = bool(n_uniq <= 2 and np.all((fin == 0.0) | (fin == 1.0)))
        if not binary:
            any_non_binary = True
        global_min = min(global_min, float(fin.min()))
        global_max = max(global_max, float(fin.max()))

        valid = np.isfinite(vals)
        mask2d = valid.any(axis=0)
        time_invariant = bool(np.array_equal(
            valid, np.broadcast_to(mask2d, valid.shape)))
        # The FINITE mask is useless on this product -- the publisher zero-fills the whole
        # globe, so `mask2d` is ~100% everywhere. What carries the footprint is where the
        # field is ever NON-ZERO, which is the cropland the model actually simulates.
        nz2d = (np.nan_to_num(vals) > 0).any(axis=0)
        model_mask.setdefault(info["model"], np.zeros(mask2d.shape, bool))
        model_mask[info["model"]] |= mask2d
        nonzero_mask.setdefault(info["model"], np.zeros(nz2d.shape, bool))
        nonzero_mask[info["model"]] |= nz2d

        for name, ilat, ilon, _, _ in sites:
            v = vals[:, ilat, ilon]
            site_series[name][f"{info['member']}|{info['scenario']}"] = (
                float(np.nanmean(v)) if np.isfinite(v).any() else float("nan"))

        rows.append(dict(
            member=info["member"], model=info["model"], scenario=info["scenario"],
            n_unique=int(n_uniq), vmin=float(fin.min()), vmax=float(fin.max()),
            exact0=exact0, exact1=exact1, binary=binary,
            land_cells=int(mask2d.sum()),
            land_frac=float(mask2d.mean()),
            nonzero_cells=int(nz2d.sum()),
            time_invariant_mask=time_invariant,
            n_years=int(yrs.size), y0=int(yrs.min()), y1=int(yrs.max()),
        ))
        if i % 10 == 0 or i == len(files):
            print(f"  scanned {i}/{len(files)}", flush=True)

    # ---- 1. metadata consistency --------------------------------------------- #
    print("\n" + "-" * 78)
    print("1. METADATA (per-member, not assumed uniform)")
    print("-" * 78)
    for k, v in attr_sets.items():
        flag = "" if len(v) == 1 else "   <-- DIVERGES ACROSS MEMBERS"
        print(f"  {k:<12} {sorted(v)}{flag}")

    # ---- 2. data nature ------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("2. DATA NATURE (GUARDRAILS 9 -- from the values)")
    print("-" * 78)
    n_bin = sum(r["binary"] for r in rows)
    uniq_counts = sorted({r["n_unique"] for r in rows})
    print(f"  binary members: {n_bin}/{len(rows)}")
    print(f"  distinct-value counts across members: {uniq_counts[:8]}"
          + (" ..." if len(uniq_counts) > 8 else ""))
    print(f"  global value range: [{global_min:.6g}, {global_max:.6g}]")
    z = np.array([r["exact0"] for r in rows])
    o = np.array([r["exact1"] for r in rows])
    print(f"  exact-0 share of finite cell-years: min {z.min():.2%}  "
          f"median {np.median(z):.2%}  max {z.max():.2%}")
    print(f"  exact-1 share: min {o.min():.2%}  median {np.median(o):.2%}  "
          f"max {o.max():.2%}")
    if n_bin == len(rows):
        verdict = ("BOOLEAN {0,1} in every member -> pooled_mean_boolean branch "
                   "(mean +/- 1 SD)")
    elif any_non_binary and n_bin > 0:
        verdict = ("MIXED -- some members binary, some continuous. This is a STOP: the "
                   "contract takes ONE branch per layer and a mixed ensemble means the "
                   "members are not the same quantity.")
    else:
        verdict = ("CONTINUOUS -> median/IQR branch UNLESS the decade pool is degenerate "
                   "at zero; see the zero-inflation section below")
    print(f"  VERDICT: {verdict}")

    # ---- 3. zero-inflation regime -------------------------------------------- #
    print("\n" + "-" * 78)
    print("3. ZERO-INFLATION (decides whether the third branch is even in scope)")
    print("-" * 78)
    print(f"  annual exact-zero share, ensemble median: {np.median(z):.2%}")
    print("  reference points: let 97.84% (took pooled_mean_zero_inflated),")
    print("                    burntarea 29.2% (did NOT qualify, median branch)")
    print("  NOTE: this is the ANNUAL share. The third branch is justified only by the "
          "exact-zero share of the MEDIAN of the decade POOL, which the processor must "
          "measure on the 2020s panel. Do not decide from this number alone.")

    # ---- 4. mask survey ------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("4. LAND / CROPLAND MASK")
    print("-" * 78)
    bad_mask = [r for r in rows if not r["time_invariant_mask"]]
    print(f"  time-invariant per member: {len(rows) - len(bad_mask)}/{len(rows)}")
    if bad_mask:
        print(f"  <-- {len(bad_mask)} member(s) have a TIME-VARYING mask; "
              "load_member() must not assume otherwise")
    by_model = defaultdict(list)
    for r in rows:
        by_model[r["model"]].append(r["land_cells"])
    print(f"  {'model':<14} {'finite min':>11} {'finite max':>11} {'globe %':>9}")
    for m in sorted(by_model):
        lo, hi = min(by_model[m]), max(by_model[m])
        print(f"  {m:<14} {lo:>11,} {hi:>11,} "
              f"{100.0 * hi / (len(lats) * len(lons)):>8.1f}%")
    print("  The FINITE mask is ~100% of the grid: this product ZERO-FILLS ocean, ice and "
          "non-cropland rather than masking them. So `isfinite` carries no footprint here "
          "and the publication mask must come from where the field is ever NON-ZERO.")

    # Is the global coverage a real ocean leak (the `floodedarea` failure mode) or just
    # zero-fill? Those look identical in a finite-mask survey and are completely different
    # defects. Decide it against the official ISIMIP3b land-sea mask.
    land = None
    try:
        sys.path.insert(0, str(root / "scripts"))
        from utils.land_mask import get_isimip_landmask  # noqa: E402
        lm = xr.open_dataset(get_isimip_landmask("3b"))["mask"].values
        land = np.nan_to_num(lm) > 0
    except Exception as e:  # noqa: BLE001
        print(f"  (land mask unavailable: {type(e).__name__}: {e})")

    if nonzero_mask:
        print(f"\n  NON-ZERO-EVER footprint (the cropland each model actually simulates):")
        hdr = f"  {'model':<14} {'cells':>9}"
        if land is not None:
            hdr += f" {'on land':>9} {'OFF LAND':>9}"
        print(hdr)
        for m in sorted(nonzero_mask):
            a = nonzero_mask[m]
            line = f"  {m:<14} {int(a.sum()):>9,}"
            if land is not None:
                line += f" {int((a & land).sum()):>9,} {int((a & ~land).sum()):>9,}"
            print(line)

        stack = np.array([nonzero_mask[m] for m in sorted(nonzero_mask)])
        nmod = stack.sum(axis=0)
        union = nmod > 0
        print(f"\n  union footprint: {int(union.sum()):,} cells")
        if land is not None:
            off = int((union & ~land).sum())
            print(f"    on land : {int((union & land).sum()):,} "
                  f"({100.0 * (union & land).sum() / land.sum():.1f}% of the "
                  f"{int(land.sum()):,} ISIMIP3b land cells)")
            print(f"    OFF land: {off:,}")
            if off == 0:
                print("    -> NOT an ocean leak. Every non-zero value sits on land; the "
                      "global extent is zero-fill only. This is the key difference from "
                      "`floodedarea`, which carries REAL values over open ocean.")
            else:
                print("    -> OCEAN LEAK: non-zero values off the land mask. STOP.")
            allzero_land = int((land & ~union).sum())
            print(f"    land cells all-zero in EVERY member: {allzero_land:,} "
                  f"({100.0 * allzero_land / land.sum():.1f}% of land) -- desert, ice, "
                  f"forest and any cropland that never fails")
        print("  coverage tiers (cells inside exactly K models' footprints):")
        for k in range(1, len(nonzero_mask) + 1):
            c = int((nmod == k).sum())
            if c:
                print(f"    {k:>2} model(s): {c:>8,}")

    # ---- 5. reference sites (GUARDRAILS 12) ----------------------------------- #
    print("\n" + "-" * 78)
    print("5. REFERENCE SITES (GUARDRAILS 12 -- is it non-trivial where crops grow?)")
    print("-" * 78)
    print(f"  {'site':<26} {'lat':>7} {'lon':>8} {'members':>8} {'NaN':>6} "
          f"{'mean':>9} {'max':>9}")
    stops = []
    for name, ilat, ilon, la, lo in sites:
        vals = np.array(list(site_series[name].values()), dtype=float)
        n_nan = int(np.isnan(vals).sum())
        finite = vals[np.isfinite(vals)]
        mean = float(finite.mean()) if finite.size else float("nan")
        mx = float(finite.max()) if finite.size else float("nan")
        flag = ""
        if finite.size == 0:
            flag = "  <-- STOP: NaN in EVERY member"
            stops.append(name)
        elif np.allclose(finite, 0.0):
            flag = "  <-- STOP: exactly 0 in every member"
            stops.append(name)
        elif n_nan:
            flag = f"  <-- masked in {n_nan} member(s)"
        print(f"  {name:<26} {la:>7.1f} {lo:>8.1f} {len(vals):>8} {n_nan:>6} "
              f"{mean:>9.4f} {mx:>9.4f}{flag}")

    print("\n  Reading: these are means over ALL years 2015-2100, so a small non-zero "
          "value is expected for an 'unprecedented event' index that is rare early and "
          "rises later. A NaN or an exact zero is the failure mode.")

    # ---- summary -------------------------------------------------------------- #
    out = root / "data" / "raw" / LAYER_ID / "value_check.json"
    with open(out, "w") as fh:
        json.dump(dict(
            n_members=len(rows), binary_members=n_bin,
            value_range=[global_min, global_max],
            exact0_median=float(np.median(z)), exact1_median=float(np.median(o)),
            attrs={k: sorted(v) for k, v in attr_sets.items()},
            rows=rows,
            sites={n: site_series[n] for n in site_series},
        ), fh, indent=2)
    print(f"\nwrote {out}")

    if stops:
        print(f"\nSTOP -- {len(stops)} reference site(s) are NaN or identically zero: "
              f"{stops}. Investigate before writing any output (GUARDRAILS 12).")
        return 4
    print("\nNo reference-site STOP condition. Data nature recorded above governs the "
          "processor's branch -- copy the measured numbers into its docstring.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
