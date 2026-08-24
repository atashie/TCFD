#!/bin/bash
# A3 processing driver: hardened wrapper + fatal validator per variable, in sequence.
# potevap is EXCLUDED pending the Gate A3 potevap-H08 user decision.
set -u
PY=/Users/arik/github/TCFD/.venv/bin/python
S=/Users/arik/github/TCFD/scripts
M=/Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_water_download.csv
R=/Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_water_download_results.csv
OUT=/Users/arik/github/TCFD/data/processed/waterindex-ssp
LOG=/Users/arik/github/TCFD/data/raw/a3_processing.log
mkdir -p "$OUT"
run_var () {
  local var=$1; shift
  local expect=$1; shift
  echo "=== A3: $var start $(date -u +%H:%M:%SZ) ===" >> "$LOG"
  /usr/bin/caffeinate -i $PY $S/process_water_hardened.py "$var" \
     --manifest "$M" --results "$R" \
     --output "$OUT/waterIndexUnderlyingData_${var}_ssp.nc" "$@" >> "$LOG" 2>&1 \
  && /usr/bin/caffeinate -i $PY $S/validate_water_variable.py \
     "$OUT/waterIndexUnderlyingData_${var}_ssp.nc" \
     "$(if [ "$var" = precip ]; then echo pr; else echo $var; fi)" \
     --expect-units "$expect" >> "$LOG" 2>&1
  local rc=$?
  echo "=== A3: $var exit $rc $(date -u +%H:%M:%SZ) ===" >> "$LOG"
  return $rc
}
overall=0
run_var dis "m3 s-1" || overall=1
run_var qr "kg m-2 s-1" || overall=1
run_var rootmoist "% max capacity" || overall=1
run_var tws "1 (robust-z synthetic, target 1000/200)" --normalize --units-override "1 (robust-z synthetic, target 1000/200)" || overall=1
echo "=== A3 driver finished overall=$overall ===" >> "$LOG"
exit $overall
