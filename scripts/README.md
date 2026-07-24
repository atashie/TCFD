# TCFD Utility Scripts

Standalone Python scripts for data processing and visualization.

## Scripts

### generate_maps.py

Generates interactive Plotly maps from processed climate data.

```bash
python scripts/generate_maps.py
```

**Input**: `data/processed/*.nc` (processed NetCDF files)
**Output**: `reports/maps/{variable}/` (HTML map files)

Features:
- Per-scenario map files (SSP126, SSP370, SSP585)
- Multiple metrics: median, percentile, trend, confidence, change, anomaly
- Scientific notation on colorbars
- Index page with navigation grid

### process_qg.py

Processes raw groundwater runoff (qg) data from ISIMIP.

```bash
python scripts/process_qg.py
```

**Input**: `data/raw/*.nc` (raw NetCDF files from ISIMIP)
**Output**: `data/processed/qg_*.nc` (processed files by scenario)

### process_led_drought.py

Processes `led` drought exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. `led` is a **binary** per-cell annual flag, so the decadal statistic is the **mean** (drought frequency), not the median. Handles the `years since 1661` / `360_day` calendar, `(model, gcm)` ensemble members, shared 2020s baseline, and zero-inflation-aware percentile/CI. See WORKFLOW-ISSUES.md 2026-07-24.

```bash
python scripts/process_led_drought.py
```

**Input**: `data/raw/drought_led_annual/*_2006_2099.nc4`
**Output**: `data/processed/drought_led_annual/led_{rcp26,rcp60}_processed.nc`

### process_let_cyclone.py

Processes `let` tropical-cyclone exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. Unlike `led`, `let` is a **continuous fraction** [0,1) per cell (area of the cell exposed to TCs), so the decadal **mean** is the expected annual exposed-area fraction. Because the ensemble is thin (1 impact model `ke-tg-meanfield` × 4 GCMs = 4 members), it applies **5×5 exponential-decay spatial smoothing** (L=0.7, cos(lat)-scaled, non-NaN-normalized) to each per-member decadal map, and a **two-tier percentile** for the zero-inflated field (zeros→1; non-zeros ranked vs the non-zero 2020s baseline → [2,100]). See WORKFLOW-ISSUES.md 2026-07-24.

```bash
python scripts/process_let_cyclone.py
```

**Input**: `data/raw/drought_let_annual/*_2006_2099.nc4`
**Output**: `data/processed/cyclone_let_annual/let_{rcp26,rcp60}_processed.nc`

### generate_qa_report.py

Generates QA reports for processed data.

```bash
python scripts/generate_qa_report.py
```

**Output**: `reports/qg_qa_report.json`

## Dependencies

These scripts require the same dependencies as the `isimip-pipeline` package:
- xarray
- plotly
- numpy
- pandas

Install with: `pip install -e isimip-pipeline/`
