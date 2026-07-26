# ISIMIP Climate Data Pipeline - Development Guide

This codebase pulls, processes, and validates ISIMIP climate model data for two distinct downstream products. Each product has its own processing logic, output format, and tooling. Workflows are consistent within each product class but differ across them. **Never confuse them.**

## Two Data Products

### 1. TCFD / CDP Reporting (6 value classes, annualized)

Processes ISIMIP data into annualized decadal statistics for physical climate risk assessment.

- **Purpose**: TCFD and CDP climate risk disclosures (timber, fisheries, health, agriculture, drought, tropical cyclones, wildfire, soil/subsurface carbon, etc.)
- **Output**: `{variable}_{scenario}_processed.nc` — 6 value classes (median, percentile, trend, significance, lower_ci, upper_ci)
- **Tooling**: `isimip-pipeline` CLI + `scripts/process_*.py` (e.g., `process_qg.py`, `process_fish_b30cm.py`, `process_led_drought.py`, `process_let_cyclone.py`, `process_burntarea_fire.py`, `process_csoil_soilcarbon.py`)
- **Skills**: `/isimip-search-download`, `/isimip-process-visualize`, `/isimip-extract-aggregate`
- **Visualization**: `scripts/generate_maps.py` (diverging trend/change maps auto-reverse to red=worse when a processed file sets `percentile_direction: higher_is_worse`)
- **Key concepts**: Shared 2020s baseline, adaptive windowing, kernel smoothing, Theil-Sen trends, percentile-of-score ranking
- **Lange 2020 exposure family** (`led` drought, `leh` heatwave, `lew` wildfire, `ler` river-flood, `lec` crop-failure, `let` tropical-cyclone — ISIMIP2b annual, rcp26/rcp60; enumerated in `config/isimip_search_catalog.yaml` → `drought.exposure_lange2020.family`): decadal statistic is the **mean** (exposed-area frequency/fraction), never the median. **Data nature varies by member and MUST be value-checked, not assumed** — `led` is a **binary** {0,1} per-cell flag; `let` is a **continuous fraction** [0,1) (both verified). Dedicated processors (`process_led_drought.py`, `process_let_cyclone.py`); `process_qg.py` cannot be reused (it parses GCM from filename field [1], but lange2020 puts the impact model there). For **thin ensembles** (e.g. `let` = 1 model × 4 GCMs), apply **5×5 exponential-decay spatial smoothing** to per-member decadal maps to borrow strength from neighbors. For **zero-inflated** hazards, use a **two-tier percentile** (zeros→1; non-zeros ranked vs the non-zero 2020s baseline → [2,100]). CIs = mean ± 1 inter-member SD, clamped. See WORKFLOW-ISSUES.md 2026-07-24 and GUARDRAILS.md §8–§10.
- **Wildfire has two representations** — don't conflate them: the Lange 2020 `lew` *exposure* member (above), and the **direct `burntarea` burnt-area fraction** (ISIMIP2b `biomes`, in **%**), processed via `process_burntarea_fire.py` (3 fire models `lpj-guess/lpjml/mc2-usfs` × 4 GCMs × rcp26/60/85 = 12 members/scenario; raw %, **no normalization** since all models share the % unit, **no spatial smoothing** — the ensemble is thick). Its `trend` is a **baseline-anchored rate** `(median[decade] − median[2020s]) / elapsed decades` (% decade⁻¹, ∝ the change map), **not** a within-decade slope — fire is too noisy year-to-year. Per-member metadata **diverges** (mc2-usfs `days since` vs others' `years since`; lpj-guess mislabels `long_name` "Fire Return Interval" and floors at 0.1%) — value-check each. See WORKFLOW-ISSUES.md 2026-07-24 and GUARDRAILS.md §9–§10.
- **Soil / subsurface carbon storage** — `csoil-total` (soil organic carbon, ISIMIP3b `biomes`, **kg C m⁻²**), processed via `process_csoil_soilcarbon.py` (3 models `classic`/`jules-es-vn6p3`/`mc2-usfs` × CMIP6 GCMs (2+5+5) × ssp126/370/585 = 12 members/scenario). The direct **storage** pool — distinct from the vegetation pools `croot`/`cvegbg` (root-zone biomass, ISIMIP3b-only) and the net-sink **flux** `nbp` (sequestration rate). **`higher_is_better`**: stored carbon is an asset, so the risk is **loss** (percentile **inverted**; red = decline). Shared unit + comparable magnitudes (2020s medians ~5.8/7.7/10.3) → **no normalization** (model democracy); no smoothing (thick); baseline-anchored trend. **Mixed CO₂**: `jules` publishes only its fixed-2015-CO₂ csoil run (trend muted), `classic`/`mc2` transient — retained as 12 members (value-check experiment tokens per model; a uniform CO₂ treatment isn't guaranteed). Layer begins at the 2020s baseline (3b starts 2015). See WORKFLOW-ISSUES.md 2026-07-25 and GUARDRAILS.md §8–§9.

### 2. Water Risk Index (20 value types, monthly)

Processes monthly ISIMIP data for 6 water variables into per-month ensemble means plus annual quantile breakpoints, feeding a dedicated Water Risk Index product.

- **Purpose**: Water risk scoring (total water storage, discharge, runoff, evapotranspiration, soil moisture, precipitation)
- **Output**: `C:\Cai_data\WaterIndex\waterIndexUnderlyingData_{var}_ssp.nc` — dimensions `(lat=360, lon=720, scenario=3, value_type=20, decade=9)`
- **Tooling**: Standalone scripts only (NOT the `isimip-pipeline` CLI)
- **Key concepts**: Per-month ensemble means (vt 0-11), annual mean (vt 12), annual quantile breakpoints Q05-Q95 (vt 13-19). Quantile annual aggregation always uses mean (not sum) to keep vt 13-19 in same units as vt 0-12. No trends, no percentile scoring, no kernel smoothing.
- **Normalization**: Robust z-score per impact model (median/IQR from 2015-2024 reference period → target mean=1000, SD=200) applied before ensemble averaging **only when model scales diverge significantly** (e.g., TWS). Per-variable decision documented in memory.
- **QA/QC**: `validate_water_tws.py` (quantile ordering, annual mean consistency, seasonal sanity, cross-scenario checks); `compare_water_index.py` (trend-focused RCP vs SSP HTML comparison with Theil-Sen slope maps and spatial Spearman R²)
- **Units**: Output units match the original RCP files for each variable. See `config_water_variables.py` for per-variable `unit_conversion_factor`.

| Variable | Aggregation | Output Units | Notes |
|----------|-------------|--------------|-------|
| tws | **mean** | kg m-2 (normalized) | Stock; 4 models, normalized to synthetic units |
| rootmoist | **mean** | % max capacity | Stock; WEB-DHM-SG only (÷1187.29 × 100) |
| qr | **sum** | kg m-2 s-1 | Flux; 4 models, raw ISIMIP units, no normalization |
| dis | **mean** | m3 s-1 | Stock; 5 models, no normalization, raw ISIMIP units |
| potevap | **sum** | kg m-2 s-1 | Flux; 4 models, h08 selectively normalized to reference ensemble (cwatm/miroc/watergap2-2e) |
| precip | **sum** | TBD | TODO — climate forcing InputData, not model output |

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
│   ├── diagnose_tws_models.py # Water: model distribution diagnostics
│   └── download_water_*.py    # Water: per-variable ISIMIP download scripts
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
