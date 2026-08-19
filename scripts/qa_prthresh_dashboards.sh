#!/usr/bin/env bash
# Generate the QA/QC artefacts for every precipitation layer, once stage 2 has finished.
#
#   1. test_shared_baseline.py per metric -- the OUTPUT-SPEC contract
#   2. generate_maps.py per metric        -> reports/maps/{metric}/    (6-tab dashboard)
#   3. render_contact_sheet.py per metric -> reports/contact_sheets/   (per-member PNG)
#   4. a cross-metric containment check that no single-layer verifier can do
#
# Contact sheets are not optional polish. Every statistic the contract checks is INVARIANT
# UNDER SPATIAL REARRANGEMENT, so none of them can see a spatial defect: a ~4x5 degree
# member once passed the full table twice and 37 algebraic checks, and a human looking at a
# picture caught it.
#
# Step 4 is the check no per-layer verifier can perform. A higher depth must be a strict
# subset of a lower one, so wetdays >= R10mm >= R20mm >= R50mm >= R100mm, R95pD >= R99pD,
# and Rx5day >= Rx1day (a 5-day window contains its own wettest day). A flipped comparison
# or a swapped threshold cannot survive it, and ten passing contract checks would not catch
# it. The same class of check found the flood layers' "additive decomposition that was
# actually containment".
set -u
cd "$(dirname "$0")/.."
PY=.venv/bin/python3
METRICS="R10mm R20mm R50mm R100mm R95pD R99pD wetdays Rx1day Rx5day prcptot"
mkdir -p reports/contact_sheets

echo "=== 1. OUTPUT-SPEC contract ==="
fails=0
for m in $METRICS; do
  d="data/processed/pluvial-isimip3b_${m}_annual"
  if [ ! -d "$d" ]; then
    printf "  %-8s MISSING (%s)\n" "$m" "$d"; fails=1; continue
  fi
  out=$($PY scripts/test_shared_baseline.py "$d" --var "$m" 2>&1)
  last=$(echo "$out" | tail -1)
  nf=$(echo "$out" | grep -c '\[FAIL\]')
  printf "  %-8s FAILs=%-3s %s\n" "$m" "$nf" "$last"
  if [ "$last" != "ALL CHECKS PASSED" ]; then fails=1; fi
done

echo
echo "=== 2. dashboards ==="
for m in $METRICS; do
  d="data/processed/pluvial-isimip3b_${m}_annual"
  [ -d "$d" ] || continue
  if $PY scripts/generate_maps.py "$m" "$d" "reports/maps/$m" >/dev/null 2>&1; then
    printf "  %-8s reports/maps/%s/%s/index.html\n" "$m" "$m" "$m"
  else
    printf "  %-8s DASHBOARD FAILED\n" "$m"; fails=1
  fi
done

echo
echo "=== 3. per-member contact sheets ==="
for m in $METRICS; do
  f="data/processed/pluvial-isimip3b_${m}_annual/${m}_members.nc"
  [ -f "$f" ] || continue
  if $PY scripts/render_contact_sheet.py "$f" "reports/contact_sheets/${m}_members.png" \
       --cols 5 >/dev/null 2>&1; then
    printf "  %-8s reports/contact_sheets/%s_members.png\n" "$m" "$m"
  else
    printf "  %-8s CONTACT SHEET FAILED\n" "$m"
  fi
done

echo
echo "=== 4. containment: a higher depth must be a SUBSET of a lower one ==="
$PY - <<'PYEOF'
import numpy as np, xarray as xr
from pathlib import Path
P = Path("data/processed")


def load(m, s):
    f = P / f"pluvial-isimip3b_{m}_annual" / f"{m}_{s}_processed.nc"
    return xr.open_dataset(f)["median"].values if f.exists() else None


NEST = [("wetdays", "R10mm"), ("R10mm", "R20mm"), ("R20mm", "R50mm"),
        ("R50mm", "R100mm"), ("R95pD", "R99pD"), ("Rx5day", "Rx1day")]
bad = 0
for s in ("ssp126", "ssp370", "ssp585"):
    for lo, hi in NEST:
        A, B = load(lo, s), load(hi, s)
        if A is None or B is None:
            print(f"  {s} {lo:>7s} >= {hi:<7s} SKIPPED (not built)")
            continue
        fin = np.isfinite(A) & np.isfinite(B)
        # The counts are exact integers pooled to a mean, so a hair of float slack is
        # enough. Rx5day/Rx1day are float sums of daily values, hence a larger tolerance.
        tol = 0.05 if lo == "Rx5day" else 1e-4
        viol = (B > A + tol) & fin
        n = int(viol.sum())
        bad += n
        worst = float((B - A)[viol].max()) if n else 0.0
        flag = "" if n == 0 else f"   <-- {n} cells, worst +{worst:.3f}"
        print(f"  {s} {lo:>7s} >= {hi:<7s} violations {n:6d}{flag}")
print(f"\nTOTAL CONTAINMENT VIOLATIONS: {bad} -> {'PASS' if bad == 0 else 'FAIL'}")
PYEOF

echo
if [ "$fails" -eq 0 ]; then
  echo "ALL QA STEPS PASSED"
else
  echo "SOME QA STEPS FAILED -- see above"
fi
