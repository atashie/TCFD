#!/usr/bin/env python
"""Verify a processed layer against the TCFD output contract (OUTPUT-SPEC.md).

    python scripts/test_shared_baseline.py {processed_dir} [--var csoil]

Restored and extended. The original checked only the shared 2020s baseline; it was
dropped during the S3 era and no layer built since carried it, so the checks that had
to be re-derived by hand each time are folded in here.

Run this after every processing run -- this script is the contract;
`scripts/generate_layer_qa.py` writes the human-readable QA record beside it.

Exit code 0 = all checks passed, 1 = at least one failed.
"""

import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

CONTRACT_VARS = ["median", "lower_ci", "upper_ci", "percentile",
                 "ols_slope", "sen_slope"]

_fail = 0


def check(label, passed, detail=""):
    global _fail
    if not passed:
        _fail += 1
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {label}" + (f"  -- {detail}" if detail else ""))
    return passed


def load(processed_dir, var):
    pattern = os.path.join(processed_dir, f"{var}_*_processed.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"ERROR: no files matching {pattern}")
        sys.exit(2)
    out = {}
    for f in files:
        base = os.path.basename(f)
        scen = base[len(var) + 1:-len("_processed.nc")]
        out[scen] = xr.open_dataset(f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_dir")
    ap.add_argument("--var", default=None,
                    help="variable prefix; inferred from the directory if omitted")
    args = ap.parse_args()

    var = args.var
    if var is None:
        cands = sorted({os.path.basename(f).split("_")[0]
                        for f in glob.glob(os.path.join(args.processed_dir,
                                                        "*_processed.nc"))})
        if len(cands) != 1:
            print(f"ERROR: could not infer --var (candidates: {cands})")
            sys.exit(2)
        var = cands[0]

    ds = load(args.processed_dir, var)
    scen = sorted(ds)
    d0 = ds[scen[0]]
    print(f"Layer '{var}' -- scenarios {scen}, dims {dict(d0.sizes)}")

    print("\nSchema")
    missing = [v for v in CONTRACT_VARS if v not in d0.data_vars]
    check("all contract variables present", not missing,
          f"missing {missing}" if missing else "")
    check("legacy single 'trend' variable is gone", "trend" not in d0.data_vars,
          "found 'trend' -- layer predates OUTPUT-SPEC.md" if "trend" in d0.data_vars else "")
    check("decadal_statistic declared", "decadal_statistic" in d0.attrs,
          d0.attrs.get("decadal_statistic", ""))
    check("slope_units declared", "slope_units" in d0.attrs,
          d0.attrs.get("slope_units", ""))
    check("percentile_direction declared", "percentile_direction" in d0.attrs,
          d0.attrs.get("percentile_direction", ""))

    # ---- Grid ------------------------------------------------------------------------ #
    # 0.5 deg is the DEFAULT and preferred grid, not the only permitted one. What is
    # required is that the grid be REGULAR, GLOBAL and IDENTICAL ACROSS SCENARIOS, and that
    # the file say what it is. Until 2026-08-14 nothing here asserted anything about the
    # grid at all, so the first 0.25 deg layer passed every check while contradicting
    # OUTPUT-SPEC.md -- a silent mismatch is worse than a loud one.
    print("\nGrid")
    lat0, lon0 = d0["lat"].values, d0["lon"].values
    dlat, dlon = np.abs(np.diff(lat0)), np.abs(np.diff(lon0))
    regular = (dlat.size and dlon.size
               and np.allclose(dlat, dlat[0], rtol=0, atol=1e-6)
               and np.allclose(dlon, dlon[0], rtol=0, atol=1e-6))
    check("lat/lon spacing is regular", bool(regular),
          "" if regular else f"lat {dlat.min():.6f}-{dlat.max():.6f}, "
                             f"lon {dlon.min():.6f}-{dlon.max():.6f}")
    cell = float(dlat[0]) if regular else float("nan")
    square = regular and abs(cell - float(dlon[0])) <= 1e-6
    check("cells are square", bool(square),
          f"{cell:.4f} deg" if square else f"lat {cell:.6f} vs lon {float(dlon[0]):.6f}")
    check("grid covers the globe",
          bool(regular) and abs(lat0.size * cell - 180.0) < 1e-3
          and abs(lon0.size * float(dlon[0]) - 360.0) < 1e-3,
          f"{lat0.size} x {lon0.size}")
    check("latitude is descending (product convention)",
          lat0.size > 1 and lat0[0] > lat0[-1],
          f"{lat0[0]:.3f} -> {lat0[-1]:.3f}")
    # Every scenario must sit on the SAME grid: the shared-baseline and percentile
    # machinery compares panels cell-by-cell across scenarios.
    same_grid = all(np.array_equal(ds[s]["lat"].values, lat0)
                    and np.array_equal(ds[s]["lon"].values, lon0) for s in scen)
    check("all scenarios share one grid", same_grid)
    declared = d0.attrs.get("spatial_resolution_degrees")
    if declared is None:
        check("spatial_resolution_degrees declared", False,
              f"<absent> -- measured {cell:.4f} deg; add it to the processor's attrs")
    else:
        check("spatial_resolution_degrees matches the coordinates",
              abs(float(declared) - cell) <= 1e-6, f"declared {declared}, measured {cell:.4f}")
    if regular and abs(cell - 0.5) > 1e-6:
        print(f"  [INFO] non-default grid: {cell:.4f} deg ({lat0.size}x{lon0.size}). "
              "0.5 deg is the product default; this layer must declare why in its "
              "processor docstring, and delivery geometry follows its own grid.")

    # The baseline is NOT always decade index 0. Layers sourced from ISIMIP3b start at
    # 2020 (baseline == index 0), but ISIMIP2b layers such as `let`/`led` carry a full
    # 2010s panel first. Locate it from the declared attribute; assuming index 0 reads
    # the 2010s panel, which legitimately DIFFERS across scenarios, and reports a
    # spurious failure.
    decades = list(d0["decade"].values)
    b_dec = int(d0.attrs.get("baseline_decade", decades[0]))
    b_i = decades.index(b_dec) if b_dec in decades else 0
    post = range(b_i + 1, len(decades))          # panels with an elapsed period
    pre = range(0, b_i + 1)                      # baseline and anything before it

    print(f"\nShared {b_dec}s baseline (decade index {b_i})")
    base_ok = True
    for a in scen:
        for b in scen:
            for v in ("median", "lower_ci", "upper_ci", "percentile"):
                if v in d0.data_vars and not np.array_equal(
                        ds[a][v].values[b_i], ds[b][v].values[b_i], equal_nan=True):
                    base_ok = False
    check(f"{b_dec}s panel bit-identical across scenarios", base_ok)
    check("baseline_source recorded",
          d0.attrs.get("baseline_source") == "shared_across_all_scenarios",
          d0.attrs.get("baseline_source", "<absent>"))

    if len(scen) > 1:
        last = d0.sizes["decade"] - 1
        a, b = scen[0], scen[-1]
        ma, mb = ds[a]["median"].values[last], ds[b]["median"].values[last]
        f = np.isfinite(ma) & np.isfinite(mb)
        check("final decade DIFFERS across scenarios",
              f.any() and not np.allclose(ma[f], mb[f]),
              f"{a} vs {b}")
    else:
        print("  [SKIP] only one scenario -- cross-scenario checks NOT TESTED")

    print("\nValue classes")
    viol = 0
    for s in scen:
        m, lo, hi = (ds[s][k].values for k in ("median", "lower_ci", "upper_ci"))
        f = np.isfinite(m) & np.isfinite(lo) & np.isfinite(hi)
        viol += int(((lo[f] > m[f] + 1e-5) | (m[f] > hi[f] + 1e-5)).sum())
    check("lower_ci <= median <= upper_ci", viol == 0, f"{viol} violations")

    p = np.concatenate([ds[s]["percentile"].values.ravel() for s in scen])
    p = p[np.isfinite(p)]
    check("percentile within [1, 100]", p.size and p.min() >= 1 - 1e-4
          and p.max() <= 100 + 1e-4, f"min {p.min():.2f} max {p.max():.2f}")

    direction = d0.attrs.get("percentile_direction", "higher_is_worse")
    m0, p0 = d0["median"].values[b_i], d0["percentile"].values[b_i]
    f = np.isfinite(m0) & np.isfinite(p0)
    r = float(np.corrcoef(m0[f], p0[f])[0, 1]) if f.sum() > 2 else float("nan")
    want_neg = direction == "higher_is_better"
    check(f"percentile orientation matches '{direction}'",
          (r < 0) if want_neg else (r > 0), f"corr(median, percentile) = {r:+.3f}")

    print("\nSlopes")
    # Every panel AT OR BEFORE the baseline has no elapsed period, so all of them must
    # be NaN -- not just the baseline itself.
    b_nan = all(np.all(np.isnan(ds[s][k].values[i]))
                for s in scen for k in ("ols_slope", "sen_slope") for i in pre)
    check(f"panels through {b_dec}s have NaN slopes (not 0)", b_nan,
          "a finite 0 here makes the whole ocean a finite zero" if not b_nan else "")

    later = all(np.isfinite(ds[s][k].values[list(post)]).any()
                for s in scen for k in ("ols_slope", "sen_slope"))
    check("post-baseline slopes are populated", later)

    leak = 0
    for s in scen:
        for i in post:
            mf = np.isfinite(ds[s]["median"].values[i])
            for k in ("ols_slope", "sen_slope"):
                leak += int((np.isfinite(ds[s][k].values[i]) & ~mf).sum())
    check("no slope finite where median is NaN", leak == 0, f"{leak} cells (ocean leak)")

    for s in scen:
        o = ds[s]["ols_slope"].values[list(post)]
        se = ds[s]["sen_slope"].values[list(post)]
        f = np.isfinite(o) & np.isfinite(se)
        if not f.any():
            continue
        # Report on ACTIVE cells -- those where at least one estimator is non-zero.
        # A cell that is permanently 0 (never burns, never sees a cyclone) has a
        # genuinely zero slope under BOTH estimators, so including it inflates
        # "agreement" toward 100% and dilutes the sen zero-fraction. On `let` the
        # all-cell view reads 73% agreement / 99.2% sen-zero, while the active-cell
        # view reads 4.9% / 96.9% -- opposite conclusions from the same array.
        act = f & ((o != 0) | (se != 0))
        zero_sen_all = float((se[f] == 0).mean())
        if act.any():
            agree = float((np.sign(o[act]) == np.sign(se[act])).mean())
            zero_sen = float((se[act] == 0).mean())
            print(f"  [INFO] {s}: over ACTIVE cells ({act.sum():,} of {f.sum():,}) "
                  f"ols/sen sign agreement {agree:.1%}, sen exactly zero {zero_sen:.1%} "
                  f"({zero_sen_all:.1%} over all finite cells)")
            if zero_sen > 0.5:
                print("         -> zero-inflated field: read ols_slope, not sen_slope")

    if "n_members" in d0.data_vars:
        print("\nEnsemble depth")
        nm = d0["n_members"].values
        nmod = d0["n_models"].values
        check("n_members >= 1 where finite", np.nanmin(nm) >= 1,
              f"range {np.nanmin(nm):.0f}-{np.nanmax(nm):.0f}")
        check("n_models <= n_members", np.nanmax(nmod - nm) <= 0,
              f"models {np.nanmin(nmod):.0f}-{np.nanmax(nmod):.0f}")

    print(f"\n{'ALL CHECKS PASSED' if _fail == 0 else f'{_fail} CHECK(S) FAILED'}")
    return 1 if _fail else 0


if __name__ == "__main__":
    sys.exit(main())
