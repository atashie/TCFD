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
