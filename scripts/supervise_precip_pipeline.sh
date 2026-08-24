#!/bin/bash
# Auto-restarting caffeinated supervisor for the B1 precip pipeline (receipt-resumable).
LOG=/Users/arik/github/TCFD/data/raw/precip_pipeline.log
for cycle in $(seq 1 60); do
  echo "=== precip supervisor cycle $cycle $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" >> "$LOG"
  /usr/bin/caffeinate -i /Users/arik/github/TCFD/.venv/bin/python \
    /Users/arik/github/TCFD/scripts/download_and_aggregate_precip.py \
    --manifest /Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_precip_daily.csv >> "$LOG" 2>&1
  if [ $? -eq 0 ]; then echo "=== precip supervisor: CLEAN, done ===" >> "$LOG"; exit 0; fi
  echo "=== precip supervisor: pass $cycle dirty, cooling 180s ===" >> "$LOG"
  sleep 180
done
echo "=== precip supervisor: cycle cap reached ===" >> "$LOG"
exit 1
