#!/bin/bash
# Gate A3 decided 2026-08-25: potevap RAW, H08 included un-normalized (user decision:
# subwatershed HydroATLAS anchoring re-bases levels; raw spread = structural uncertainty).
PY=/Users/arik/github/TCFD/.venv/bin/python
S=/Users/arik/github/TCFD/scripts
M=/Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_water_download.csv
R=/Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_water_download_results.csv
OUT=/Users/arik/github/TCFD/data/processed/waterindex-ssp
LOG=/Users/arik/github/TCFD/data/raw/a3_processing.log
echo "=== A3: potevap start $(date -u +%H:%M:%SZ) (Gate A3: raw, H08 included) ===" >> "$LOG"
/usr/bin/caffeinate -i $PY $S/process_water_hardened.py potevap \
  --manifest "$M" --results "$R" \
  --output "$OUT/waterIndexUnderlyingData_potevap_ssp.nc" >> "$LOG" 2>&1 \
&& /usr/bin/caffeinate -i $PY $S/validate_water_variable.py \
  "$OUT/waterIndexUnderlyingData_potevap_ssp.nc" potevap \
  --expect-units "kg m-2 s-1" >> "$LOG" 2>&1
rc=$?
echo "=== A3: potevap exit $rc $(date -u +%H:%M:%SZ) ===" >> "$LOG"
exit $rc
