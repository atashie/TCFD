#!/usr/bin/env python3
"""Contract check for `observational-historical-v1` layers.

WHY THIS EXISTS AS A SEPARATE VERIFIER
    `test_shared_baseline.py` asserts the OUTPUT-SPEC decadal contract: a `decade`
    dimension, two slopes, and an ensemble. An observational layer built from a
    single historical record has none of those and never will, so it fails that
    verifier BY DESIGN. The wrong fix is to relax the shared verifier -- it guards
    every projected layer in the registry, and its strictness is the product. The right fix is a second
    contract with its own check, which is this file.

    Neither verifier substitutes for the other, and neither is a sanity check on
    the DATA. A file can pass every check here and still be meaningless; both
    sugarcane layers passed everything and were empty. Reference-site inspection
    and a human reading the maps are separate obligations.

WHAT THIS CONTRACT REQUIRES
    dims        exactly (lat, lon). No decade, no scenario, no member axis.
    variables   median, lower_ci, upper_ci, percentile  (+ optional diagnostics)
    ordering    lower_ci <= median <= upper_ci everywhere finite
    percentile  within [1, 100], and consistent with the declared zero tier
    masking     unobserved cells are NaN, NEVER 0 -- "not observed" is not "no hazard"
    absence     ols_slope / sen_slope / n_members / n_models must be ABSENT, not
                present-and-empty and not faked to 1. A half-populated decadal
                contract is worse than an honestly different one, because a reader
                who sees n_models=1 concludes "thin ensemble" rather than "no
                ensemble exists".
    attributes  the declared caveat set, each non-empty AND non-negative

THE NON-NEGATIVE CAVEAT CHECK
    CLAUDE.md: a caveat attribute is promoted into a report on being NON-EMPTY, not
    on saying anything. Writing `resolution_caveat: "none -- the grid is fine
    enough"` therefore publishes a must-disclose caveat whose body reads "none".
    This verifier fails any caveat whose text begins with a negation, so that trap
    is caught at build time rather than in a customer's filing. Park a negative
    MEASUREMENT under a name that is not a caveat attribute.

USAGE
    python3 scripts/test_observational_baseline.py data/processed/tornado-spc_f2plus_full
    python3 scripts/test_observational_baseline.py data/processed          # all matching
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr

CONTRACT = "observational-historical-v1"

REQUIRED_VARS = ["median", "lower_ci", "upper_ci", "percentile"]

FORBIDDEN_VARS = ["ols_slope", "sen_slope", "n_members", "n_models", "trend", "significance"]

FORBIDDEN_DIMS = ["decade", "scenario", "member", "time"]

REQUIRED_ATTRS = [
    "output_contract",
    "hazard",
    "source_dataset",
    "source_url",
    "field_nature",
    "statistic",
    "statistic_rationale",
    "percentile_direction",
    "temporal_window",
    "n_years",
]

# Caveats that MUST be present and must actually say something.
REQUIRED_CAVEATS = [
    "resolution_caveat",
    "reporting_bias_caveat",
    "coverage_caveat",
    "no_trend_rationale",
]

NEGATION_OPENERS = ("none", "n/a", "na", "nan", "not applicable", "no ", "nil", "-", "")


class Result:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.checks = 0

    def check(self, label: str, ok: bool, detail: str = "") -> None:
        self.checks += 1
        if ok:
            print(f"    PASS  {label}")
        else:
            print(f"    FAIL  {label}" + (f" -- {detail}" if detail else ""))
            self.failures.append(f"{label}: {detail}" if detail else label)


def verify(path: Path, r: Result) -> None:
    print(f"\n  {path}")
    ds = xr.open_dataset(path)

    r.check("declares the observational contract",
            str(ds.attrs.get("output_contract", "")).startswith(CONTRACT),
            f"got {ds.attrs.get('output_contract', '<missing>')!r}")

    r.check("dims are exactly (lat, lon)",
            set(ds.sizes) == {"lat", "lon"}, f"got {dict(ds.sizes)}")

    present_forbidden_dims = [d for d in FORBIDDEN_DIMS if d in ds.sizes]
    r.check("no decadal/ensemble dimension", not present_forbidden_dims,
            f"present: {present_forbidden_dims}")

    missing = [v for v in REQUIRED_VARS if v not in ds.data_vars]
    r.check("required variables present", not missing, f"missing: {missing}")

    present_forbidden = [v for v in FORBIDDEN_VARS if v in ds.data_vars]
    r.check("decadal-contract variables are ABSENT, not faked",
            not present_forbidden, f"present: {present_forbidden}")

    if missing:
        return

    med = ds["median"].values
    lo = ds["lower_ci"].values
    hi = ds["upper_ci"].values
    pct = ds["percentile"].values
    finite = np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)

    viol = int(np.sum(lo[finite] > med[finite]) + np.sum(med[finite] > hi[finite]))
    r.check("lower_ci <= median <= upper_ci", viol == 0, f"{viol} violations")

    pf = np.isfinite(pct)
    r.check("percentile within [1, 100]",
            bool(pf.any()) and float(np.nanmin(pct)) >= 1.0 and float(np.nanmax(pct)) <= 100.0,
            f"range {np.nanmin(pct):.2f}-{np.nanmax(pct):.2f}" if pf.any() else "all NaN")

    # The three fields must be missing in the SAME cells. Without this, `finite` above
    # is an intersection that can silently excuse a NaN lower_ci under a finite median
    # from the ordering check -- the check would pass by not looking.
    masks = {v: np.isfinite(ds[v].values) for v in ("median", "lower_ci", "upper_ci")}
    r.check("median / lower_ci / upper_ci share one missingness pattern",
            bool(np.array_equal(masks["median"], masks["lower_ci"])
                 and np.array_equal(masks["median"], masks["upper_ci"])),
            "the three fields are finite in different cells")

    # Two-tier consistency, keyed on the OBSERVED COUNT rather than on the rate being
    # exactly zero. Under the posterior estimator the rate is positive everywhere, so a
    # `median <= 0` test would select nothing and pass vacuously -- a check that cannot
    # fail is worse than no check, because it reads as coverage.
    if "n_events" in ds.data_vars:
        cnt = ds["n_events"].values
        unobserved = finite & np.isfinite(cnt) & (cnt <= 0)
        tier_one = pf & (pct == 1.0)
        r.check("cells with no observed event sit in the percentile zero tier",
                bool(np.all(tier_one[unobserved])) if unobserved.any() else True,
                f"{int(np.sum(unobserved & ~tier_one))} unobserved cells not at percentile 1")
        r.check("percentile-1 cells are exactly the unobserved cells",
                bool(np.all(unobserved[tier_one])) if tier_one.any() else True,
                f"{int(np.sum(tier_one & ~unobserved))} tier-1 cells carry an observed event")

        # The positive tier must actually occupy [2, 100] and must be monotonic in rate,
        # otherwise "two-tier" describes the intent rather than the data.
        occupied = pf & np.isfinite(cnt) & (cnt > 0)
        if occupied.sum() > 1:
            r.check("occupied cells score in [2, 100]",
                    float(np.min(pct[occupied])) >= 2.0 and float(np.max(pct[occupied])) <= 100.0,
                    f"range {np.min(pct[occupied]):.3f}-{np.max(pct[occupied]):.3f}")
            rho = float(np.corrcoef(
                np.argsort(np.argsort(med[occupied])),
                np.argsort(np.argsort(pct[occupied])))[0, 1])
            r.check("percentile is monotonic in the rate over occupied cells",
                    rho > 0.999, f"rank correlation {rho:.6f}")

    # Unobserved must be NaN in EVERY field, never 0. Checking only `median` would let
    # an out-of-mask percentile or interval through under a label claiming otherwise.
    if "conus_mask" in ds.data_vars:
        outside = ds["conus_mask"].values == 0
        leaks = {v: int(np.sum(~np.isnan(ds[v].values[outside])))
                 for v in ("median", "lower_ci", "upper_ci", "percentile")}
        r.check("outside the coverage mask EVERY published field is NaN, never 0",
                all(n == 0 for n in leaks.values()),
                f"non-NaN out-of-mask cells: { {k: v for k, v in leaks.items() if v} }")

    # Attributes.
    missing_attrs = [a for a in REQUIRED_ATTRS if not str(ds.attrs.get(a, "")).strip()]
    r.check("required attributes present and non-empty", not missing_attrs,
            f"missing/blank: {missing_attrs}")

    for cav in REQUIRED_CAVEATS:
        text = str(ds.attrs.get(cav, "")).strip()
        r.check(f"{cav} present", bool(text), "missing or blank")
        if text:
            head = text.lower().lstrip("\"' ").split(".")[0].strip()
            negative = head in NEGATION_OPENERS or head.startswith(("none", "not applicable", "no caveat"))
            r.check(f"{cav} is not a NEGATIVE caveat", not negative,
                    f"reads {text[:60]!r} -- a non-empty caveat is published verbatim; "
                    "park negative measurements under a non-caveat attribute name")

    qa = str(ds.attrs.get("qa_reviewed_on", "")).strip().lower()
    if not qa or qa.startswith("null"):
        print("    NOTE  qa_reviewed_on is null -- no human has read the maps for this layer. "
              "Not a contract failure; a delivery blocker.")


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    root = Path(sys.argv[1])
    if root.is_file():
        targets = [root]
    else:
        targets = sorted(root.rglob("*_processed.nc"))
        keep = []
        for p in targets:
            try:
                with xr.open_dataset(p) as ds:
                    if str(ds.attrs.get("output_contract", "")).startswith(CONTRACT):
                        keep.append(p)
            except Exception:  # noqa: BLE001
                continue
        targets = keep

    if not targets:
        print(f"No {CONTRACT} files found under {root}. "
              "An empty result is not a pass -- check the path.")
        return 2

    print(f"Checking {len(targets)} file(s) against {CONTRACT}")
    r = Result()
    for p in targets:
        verify(p, r)

    print(f"\n{r.checks} checks run across {len(targets)} file(s)")
    if r.failures:
        print(f"FAILED -- {len(r.failures)} violation(s):")
        for f in r.failures:
            print(f"  * {f}")
        return 1
    print("PASS -- every file satisfies the observational contract.")
    print("This says the files are SHAPED right. It does not say the data means what its "
          "name says; that needs the reference sites and a human on the maps.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
