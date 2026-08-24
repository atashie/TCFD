#!/bin/bash
# Auto-restarting supervisor for the manifest download (idempotent per pass).
# caffeinate -i prevents idle sleep while a pass runs; passes are capped, and a
# clean reconciliation (exit 0) ends the loop.
LOG=/Users/arik/github/TCFD/data/raw/water_manifest_download.log
for cycle in $(seq 1 40); do
  echo "=== supervisor cycle $cycle $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" >> "$LOG"
  /usr/bin/caffeinate -i /Users/arik/github/TCFD/.venv/bin/python \
    /Users/arik/github/TCFD/scripts/download_water_manifest.py \
    /Users/arik/github/waterRiskIndex_beta/v2/enumeration/A1_20260821/manifest_water_download.csv \
    --base-dir /Users/arik/github/TCFD --burst 20 --cooldown 90 >> "$LOG" 2>&1
  if [ $? -eq 0 ]; then echo "=== supervisor: CLEAN RECONCILIATION, done ===" >> "$LOG"; exit 0; fi
  echo "=== supervisor: pass $cycle ended dirty, cooling 180s ===" >> "$LOG"
  sleep 180
done
echo "=== supervisor: cycle cap reached without clean pass ===" >> "$LOG"
exit 1
