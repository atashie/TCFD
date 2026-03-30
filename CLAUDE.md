# ISIMIP Climate Data Pipeline - Development Guide

This codebase pulls, processes, and validates ISIMIP climate model data for two distinct downstream products. Each product has its own processing logic, output format, and tooling. Workflows are consistent within each product class but differ across them. **Never confuse them.**

## Two Data Products

### 1. TCFD / CDP Reporting (6 value classes, annualized)

Processes ISIMIP data into annualized decadal statistics for physical climate risk assessment.

- **Purpose**: TCFD and CDP climate risk disclosures (timber, fisheries, health, agriculture, etc.)
- **Output**: `{variable}_{scenario}_processed.nc` — 6 value classes (median, percentile, trend, significance, lower_ci, upper_ci)
- **Tooling**: `isimip-pipeline` CLI + `scripts/process_*.py` (e.g., `process_qg.py`, `process_fish_b30cm.py`)
- **Skills**: `/isimip-search-download`, `/isimip-process-visualize`, `/isimip-extract-aggregate`
- **Visualization**: `scripts/generate_maps.py`
- **Key concepts**: Shared 2020s baseline, adaptive windowing, kernel smoothing, Theil-Sen trends, percentile-of-score ranking

### 2. Water Risk Index (20 value types, monthly)

Processes monthly ISIMIP data for 6 water variables into per-month ensemble means plus annual quantile breakpoints, feeding a dedicated Water Risk Index product.

- **Purpose**: Water risk scoring (total water storage, discharge, runoff, evapotranspiration, soil moisture, precipitation)
- **Output**: `C:\Cai_data\WaterIndex\waterIndexUnderlyingData_{var}_ssp.nc` — dimensions `(lat=360, lon=720, scenario=3, value_type=20, decade=9)`
- **Tooling**: Standalone scripts only (NOT the `isimip-pipeline` CLI)
- **Key concepts**: Per-month ensemble means (vt 0-11), annual mean (vt 12), annual quantile breakpoints Q05-Q95 (vt 13-19). No trends, no percentile scoring, no kernel smoothing.
- **Normalization**: Robust z-score per impact model (median/IQR from 2015-2024 reference period → target mean=1000, SD=200) applied before ensemble averaging to align models with different baseline assumptions.
- **QA/QC**: `validate_water_tws.py` (quantile ordering, annual mean consistency, seasonal sanity, cross-scenario checks); `compare_water_index.py` (trend-focused RCP vs SSP HTML comparison with Theil-Sen slope maps and spatial Spearman R²)

| Variable | Aggregation | Notes |
|----------|-------------|-------|
| tws, rootmoist | **mean** | Stock variables |
| dis, qr, potevap, precip | **sum** | Flux variables (kg/m2/s → mm/month) |

## Project Structure

```
TCFD/
├── isimip-pipeline/          # Python package (CLI) — TCFD/CDP workflow
│   ├── src/isimip_pipeline/  # cli, search, download, processing, visualization
│   └── tests/                # 299 tests
│
├── scripts/                  # Standalone scripts — both workflows
│   ├── utils/                # Shared utilities (land_mask, water_index_compare)
│   ├── process_qg.py         # TCFD: example annualized processor
│   ├── generate_maps.py      # TCFD: interactive Plotly maps
│   ├── config_water_*.py     # Water: variable configuration
│   ├── process_water_*.py    # Water: processing scripts
│   ├── validate_water_tws.py # Water: QA/QC validation
│   ├── compare_water_index.py # Water: RCP vs SSP comparison report
│   └── diagnose_tws_models.py # Water: model distribution diagnostics
│
├── config/                   # ISIMIP search catalog cache
├── data/                     # Raw + processed data (gitignored)
├── reports/                  # Generated reports (gitignored)
└── _deprecated/              # Archived legacy R code
```

## CLI Quick Reference (TCFD/CDP only)

```bash
isimip-pipeline search         # Search ISIMIP repository
isimip-pipeline download       # Batch download
isimip-pipeline process        # Process raw data
isimip-pipeline report         # Generate QA report
isimip-pipeline run            # Complete pipeline
isimip-pipeline find           # Search local processed datasets
isimip-pipeline catalog        # Manage ISIMIP catalog
isimip-pipeline cleanup        # Delete raw data after verification
```

## Skills

| Skill | When to Use | Product |
|-------|-------------|---------|
| `/isimip-search-download` | Searching ISIMIP repository, downloading data | TCFD/CDP |
| `/isimip-process-visualize` | Processing annualized NetCDF, generating QA reports | TCFD/CDP |
| `/isimip-extract-aggregate` | Extracting data by location/region for CSV export | TCFD/CDP |

## Key Documentation

| Document | Purpose |
|----------|---------|
| [GUARDRAILS.md](GUARDRAILS.md) | Critical rules that must never be violated |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log and resolutions |
