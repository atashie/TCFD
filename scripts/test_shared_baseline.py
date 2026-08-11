#!/usr/bin/env python
"""Verify a processed layer against the TCFD output contract (OUTPUT-SPEC.md).

    python scripts/test_shared_baseline.py {processed_dir} [--var csoil]

Restored and extended. The original checked only the shared 2020s baseline; it was
dropped during the S3 era and no layer built since carried it, so the checks that had
to be re-derived by hand each time are folded in here.

Run this after every processing run. `generate_qa_report.py` re-checks some of the
same invariants -- this script is the contract, that report is the safety net.

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
