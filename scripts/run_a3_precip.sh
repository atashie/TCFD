#!/bin/bash
# Precip processing: the chunk receipts (interim/water_precip_monthly/chunk_receipts.jsonl)
# play the role the download-results sha256 receipts play for model-output variables,
# so the manifest raw-audit is replaced by receipt-completeness + the fatal validator.
PY=/Users/arik/github/TCFD/.venv/bin/python
OUT=/Users/arik/github/TCFD/data/processed/waterindex-ssp
LOG=/Users/arik/github/TCFD/data/raw/a3_processing.log
echo "=== A3: precip start $(date -u +%H:%M:%SZ) ===" >> "$LOG"
/usr/bin/caffeinate -i $PY /Users/arik/github/TCFD/scripts/process_water_variable.py precip \
  --data-dir /Users/arik/github/TCFD/data/interim/water_precip_monthly \
  --output "$OUT/waterIndexUnderlyingData_precip_ssp.nc" >> "$LOG" 2>&1 || { echo "=== A3: precip PROCESS FAILED ===" >> "$LOG"; exit 1; }
# provenance from the resolved chunk receipts
/usr/bin/caffeinate -i $PY - >> "$LOG" 2>&1 <<'PYEOF'
import json, hashlib, subprocess
from netCDF4 import Dataset
R="/Users/arik/github/TCFD/data/interim/water_precip_monthly/chunk_receipts.jsonl"
recs=[json.loads(l) for l in open(R)]
members=sorted({r["monthly_out"].split("/")[-2] for r in recs})
nc=Dataset("/Users/arik/github/TCFD/data/processed/waterindex-ssp/waterIndexUnderlyingData_precip_ssp.nc","a")
nc.setncattr("impact_models","w5e5 (bias-adjusted climate forcing; not a model ensemble)")
nc.setncattr("gcms", ", ".join(sorted({m.rsplit("_",1)[0] for m in members})))
nc.setncattr("members_resolved", json.dumps(members))
nc.setncattr("source_chunk_receipts", R)
nc.setncattr("source_chunk_receipts_sha256", hashlib.sha256(open(R,'rb').read()).hexdigest())
nc.setncattr("n_source_chunks", len(recs))
nc.setncattr("aggregation_note","daily InputData -> monthly-MEAN flux (calendar-aware resample), B0-verified")
nc.setncattr("processing_code_commit", subprocess.check_output(["git","-C","/Users/arik/github/TCFD","rev-parse","--short","HEAD"],text=True).strip())
nc.close(); print("precip provenance written:", len(recs), "chunks,", len(members), "members")
PYEOF
/usr/bin/caffeinate -i $PY /Users/arik/github/TCFD/scripts/validate_water_variable.py \
  "$OUT/waterIndexUnderlyingData_precip_ssp.nc" pr --expect-units "kg m-2 s-1" >> "$LOG" 2>&1
echo "=== A3: precip exit $? $(date -u +%H:%M:%SZ) ===" >> "$LOG"
