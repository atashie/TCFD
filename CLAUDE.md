# ISIMIP Climate Data Pipeline - Development Guide

This codebase pulls, processes, and validates ISIMIP climate model data for two distinct downstream products. Each product has its own processing logic, output format, and tooling. Workflows are consistent within each product class but differ across them. **Never confuse them.**

## Two Data Products

### 1. TCFD / CDP Reporting (annualized decadal statistics)

Processes ISIMIP data into annualized decadal statistics for physical climate risk assessment.

- **Purpose**: TCFD and CDP climate risk disclosures (timber, fisheries, health, agriculture, drought, tropical cyclones, wildfire, soil/subsurface carbon, etc.)
- **Output contract**: [OUTPUT-SPEC.md](OUTPUT-SPEC.md) is authoritative — do not restate it elsewhere. `{variable}_{scenario}_processed.nc` on `(decade, lat, lon)`, carrying `median`, `lower_ci`, `upper_ci`, `percentile`, `ols_slope`, `sen_slope`, `n_members`, `n_models`.
- **Shared statistics**: `scripts/utils/decadal_stats.py` implements the contract; every processor calls it rather than rolling its own. `scripts/process_csoil_soilcarbon.py` is the **reference implementation**.
- **Verification**: `python scripts/test_shared_baseline.py {processed_dir}` after every run — exits non-zero on any contract violation.
- **Tooling**: `isimip-pipeline` CLI + `scripts/process_*.py`
- **Skills**: `/isimip-search-download`, `/isimip-process-visualize`, `/isimip-extract-aggregate`. **Invoke the skill before improvising.** Any question of the form "what data exists for X / what could we process / is X available" is `/isimip-search-download` — load it *first*. It carries the enumeration mechanics (serial harvest, `.nc4?` matching, `DerivedOutputData/`) and coverage facts that this file does not, and it is **authoritative over this file** where they disagree. Hand-rolled `curl` sweeps that skip it have re-derived documented findings and drawn absence conclusions from rate-limited empty listings.
- **Visualization**: `scripts/generate_maps.py` — six tabs (`Median | Percentile | Trend | Confidence | Anomaly | Members`); diverging slope/change panels are zero-centred on the 95th percentile of |value| and auto-reverse to red=worse when a file sets `percentile_direction: higher_is_worse`. Conventions and the browser-payload limits live in the `/isimip-process-visualize` skill — read it before changing a dashboard.
- **Key concepts**: Shared 2020s baseline, adaptive windowing, pooled (year × member) decadal statistics, dual OLS + Theil-Sen slopes on an expanding window, percentile-of-score ranking
- **Two slopes, not one** — they fail in **opposite** regimes, so neither is safe alone. `sen_slope` collapses to exactly 0 on zero-inflated hazards (91.3% of `driedarea` ssp126 cells); `ols_slope` absorbs between-member level offsets as trend when member coverage is uneven (measured +40% bias; `csoil-total`'s offset is 68.7× its interannual SD). Read `ols_slope` on zero-inflated fields, `sen_slope` otherwise, and treat **disagreement between them as the signal that a cell's trend is not robust**. This supersedes the retired single `trend` variable and the never-emitted `significance` class.
- **Per-hazard framing decisions live with their layer, not here.** Every shipped layer's specifics (units, ensemble composition, normalization/smoothing/percentile choices, per-member metadata quirks) are recorded in **its processor's module docstring** and its dated **WORKFLOW-ISSUES.md** entry. Read the one layer you are working on. Do **not** read the others for comparison unless the current task genuinely depends on them — see *Scope discipline* below.

#### Shipped TCFD layers

| Hazard | Variable | Round / scenarios | Processor |
|---|---|---|---|
| Drought (exposure) | `led` | 2b, rcp26/60 | `process_led_drought.py` |
| **Tropical cyclone (exposure)** | `let` | 2b, rcp26/60 | `process_let_cyclone.py` (rebuilt 2026-08-11) |
| **Wildfire (burnt area)** | `burntarea-total` | 3b `fire`+`biomes`, ssp126/370/585 | `process_burntarea_isimip3b.py` |
| ~~Wildfire, ISIMIP2b~~ — **SUPERSEDED 2026-08-10** | `burntarea` | 2b `biomes`, rcp26/60/85 | ~~`process_burntarea_fire.py`~~ |
| Soil / subsurface carbon | `csoil-total` | 3b `biomes`, ssp126/370/585 | `process_csoil_soilcarbon.py` |
| Timber, fisheries, health, … | various | — | `process_*.py` |

- **Availability questions stop at the inventory.** Report the model × GCM × scenario matrix with cadence, volume and the open decisions, then **end the turn**. No recommendation, no ranking, no `AskUserQuestion` to force a pick — dataset choice is a discussion the user opens, and we decide together from the table. Give your read on the trade-offs when asked. (Rule added 2026-08-08; see the skill's *Answering "what data exists for {hazard}?"*.)
- **Scope discipline.** When the user asks about **one** hazard — what data exists, what to process, how to frame it — answer within that hazard. Enumerate its options from the repository and judge each on **its own** properties: cadence, ensemble depth, units, data nature, scenario coverage, volume. Do **not** justify or rank a dataset by reference to an unrelated shipped layer ("aligns with `csoil`", "same shape as `let`"). Cross-layer comparison is warranted only when the user asks for it, or when it is a real constraint (a shared processor, the OUTPUT-SPEC contract, a guardrail). Precedent from another hazard is a *hypothesis to re-verify*, never a reason to prefer a dataset.
- **Cross-cutting framing rules** (they generalize; the per-layer instances do not): data nature is **measured, never inferred** — from the variable name, its `long_name`, a sibling member, or another layer (GUARDRAILS §9). The measured nature then determines the statistic, the CI, and the percentile treatment per [OUTPUT-SPEC.md](OUTPUT-SPEC.md). Thin ensembles get 5×5 exponential-decay spatial smoothing — but **the decay length is a per-layer measurement, not a constant**: `L=0.7` keeps 32% of the weight on the centre cell and preserves sparse structure, `L=2.5` keeps 8% and dissolves it. Zero-inflated hazards get a two-tier percentile (zeros→1, non-zeros ranked vs the non-zero 2020s baseline). Variables where a decline is the risk set `higher_is_better` and invert the percentile.
- **Three decadal-statistic branches, not two.** Continuous → pooled median/IQR; boolean {0,1} → pooled mean±SD; and **extreme zero-inflation** → `pooled_mean_zero_inflated`, mean±SD on a *continuous* field. The third exists because the boolean/continuous split is only a proxy for "is the decade pool degenerate at zero", and it fails hard: at `let`'s 97.84% annual zeros the median branch erases 93% of exposed land. It is a **declared** deviation — measure the median's exact-zero share and record it in `decadal_statistic_rationale`, never take this branch to improve contrast. `burntarea` at 29.2% zeros does not qualify; only `let` does today.
- **Extreme-event exposure families** — enumerated in `config/isimip_search_catalog.yaml` and in the `/isimip-search-download` skill, which is **authoritative over this file** for repository coverage. Two generations exist: **Lange 2020** (ISIMIP2b, rcp26/rcp60) — twelve members, `le{d,r,w,c,h,t}` land-area-exposed paired with `pe{d,r,w,c,h,t}` population-exposed twins; and its **ISIMIP3b/SSP re-issue** under `DerivedOutputData/`, split across `Heinicke2026` (`driedarea`, `floodedarea`) and `Zantout2025` (`heatwave`, `wildfire`, `cropfailure`) — named by hazard word, not `le*` code. Never conclude a family has no newer-round version without listing `DerivedOutputData/`. **And a family having no re-issue is not the hazard having no newer data** — ISIMIP3b's newest tropical-cyclone product sits in `InputData/climate/tropical_cyclones/`, not in `DerivedOutputData/` at all. List every directory level; scope each absence claim to exactly what was enumerated (2026-08-11).

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
