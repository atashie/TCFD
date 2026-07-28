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
- **Wildfire has several representations** — don't conflate them: the Lange 2020 `lew` *exposure* member (above; rcp26/rcp60 **only**, no rcp85), the `ffire` carbon-**emissions** flux, the ISIMIP3a-only `fire`-sector diagnostics (`firesize`/`firenr`/`fireints`/…, historical only — an evaluation reference, never a layer source), and the **direct `burntarea` burnt-area fraction**, which is the one processed. **Since 2026-07-28 the layer is ISIMIP3b `biomes` `burntarea-total`, ssp126/370/585** (it was ISIMIP2b/RCP before; newer round and scenario family win wherever the newer data is viable). `process_burntarea_fire.py`: **12 members/scenario** = `mc2-usfs` (annual, 5 GCMs) + `visit` (monthly, 5) + `classic` (monthly, 2); raw %, **no normalization** (same unit *and* within 1.6×), **no spatial smoothing** (thick), shared 2020s baseline (3b starts 2015 → no 2010s decade). `trend` is a **baseline-anchored rate** `(median[decade] − median[2020s]) / elapsed decades` (% decade⁻¹, ∝ the change map), **not** a within-decade slope — fire is too noisy year-to-year. Three traps: **monthly burntarea annualizes by SUM, not mean** (burnt area *accumulates* — verified against `classic`'s daily output, 1e-6 agreement; the csoil mean-annualization precedent is wrong here and would under-scale 12×); **annual values legitimately exceed 100%** where a cell reburns, so the CI is floored at 0 and **unbounded above** (clamping to 100 would push `upper_ci` below the median); and **`classic` runs at an effective 1.0°** (100% constant 2×2 blocks) — kept deliberately, while **`elm-eca` is EXCLUDED** at ~4°×5°. Also note the **declared unit cannot gate pooling**: 2b `clm45`/`orchidee` declare `%` while sitting ~1000× low on a 0–1 fraction scale. Emits **`n_members`/`n_models`** per cell (the 3 models don't share a land mask). See WORKFLOW-ISSUES.md 2026-07-28 and GUARDRAILS.md §9–§11.
- **Soil / subsurface carbon storage** — `csoil-total` (soil organic carbon, ISIMIP3b `biomes`, **kg C m⁻²**), processed via `process_csoil_soilcarbon.py`. **17 members/scenario**: 4 of the 5 models that publish it — `classic`(2 GCMs)/`jules-es-vn6p3`(5)/`mc2-usfs`(5)/`visit`(5) × ssp126/370/585. `visit` publishes csoil-total **only monthly** and is annualized by the **mean of each year's 12 months** (verified immaterial: within-year CV 0.11%). **`elm-eca` is EXCLUDED** — it declares 0.5° but is effectively **~4°×5°** (seams every 10 columns / 8 rows) and has the fattest tail, so it rendered as large bright rectangles; dropping it cost 5 land cells (0.01%). See `EXCLUDED_MODELS` and GUARDRAILS §11. Note the round-dependent name: ISIMIP**2b** uses bare `csoil`, ISIMIP**3b** uses `csoil-total`. The direct **storage** pool — distinct from the vegetation pools `croot`/`cvegbg` and the net-sink **flux** `nbp` (sequestration rate, 3b-only). **`higher_is_better`**: stored carbon is an asset, so the risk is **loss** (percentile **inverted**; red = decline). Shared unit, medians within 1.8× → **no normalization** (model democracy); no smoothing (thick); baseline-anchored trend. Also emits **`n_members`/`n_models`** per cell — the 4 models do **not** share a land mask, so ~81% of land carries all 17 members and 9.2% is single-model; those cells are **kept, not masked**. **Mixed CO₂**: `jules` publishes only its fixed-2015-CO₂ run — and contrary to the old note this does **not** mute its trend, it gives `jules` the **strongest loss** (−4.4% at ssp585) while the four transient models run flat-to-positive. Ensemble membership drives the headline: `jules`'s weight fell 42%→29%, and the global mean by 2090s is now **+2.81/+1.18/+0.72%** (ssp126/370/585) vs the 12-member run's +1.1/−1.4/−2.2%, which reproduces exactly from the same files — so this is membership, not processing. Read the sign as **contested across models**. Layer begins at the 2020s baseline (3b starts 2015). See WORKFLOW-ISSUES.md 2026-07-27/28 and GUARDRAILS.md §8–§9, §11.

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
└── _deprecated/              # Archived legacy R code
```

## Storage: S3 is canonical

**Data lives in S3, never on local disk.** See [STORAGE.md](STORAGE.md) for the full layout contract. Local `data/` and `reports/` directories are gone; the local cache is ephemeral scratch under `/tmp/tcfd-cache` (override with `TCFD_CACHE_ROOT`) and deleting it never loses data.

- **Bucket/prefix**: `s3://climate-ai-data-science-shiny-app-data/TCFD/` (`us-east-2`)
- **All keys are built by `isimip_pipeline/storage.py`** — never hardcode an S3 key.
- **A published layer version is immutable.** Reprocessing creates a new version `v{YYYY-MM-DD}_{git-sha}` (plus `-dirty` on an uncommitted tree). `publish_layer_version` refuses to overwrite an existing version unless told to.
- **Consumers read `…/layers/{layer_id}/current/`** and must call `storage.verify_complete()` first — a version prefix whose `_COMPLETE.json` is absent or mismatched is in-flight or corrupt, never data.
- `_COMPLETE.json` gates `data/` + `layer.json` only. `qa/` and `maps/` live in the version prefix (evidence stays pinned to its data) but are regenerable and ungated, so re-running `generate_maps.py` never invalidates a layer.
- **Credentials**: never pin static `AWS_*` env vars — `storage.s3_filesystem()` drops them so the SageMaker container provider auto-refreshes. Pinning a ~1h token kills long jobs mid-run.

Migrated to S3 so far: `process_led_drought.py`, `process_let_cyclone.py`, `process_burntarea_fire.py`, `process_csoil_soilcarbon.py`, `generate_maps.py`, `generate_qa_report.py`. The remaining `scripts/` still reference local paths and are pending the Python rebuild.

**Every ingest-and-process run must leave reviewable HTML evidence — and the maps are a QA gate, not a deliverable.** A layer is **unreviewed, not verified**, until someone has actually viewed an image of it: distribution statistics are invariant under spatial rearrangement and cannot see a spatial defect (GUARDRAILS §11 — a coarse `~4°×5°` member shipped past a full tabular value-check and 37 algebraic checks). Generating maps is not reviewing them. Do not declare a layer good or run `storage.cleanup_raw` before the maps have been looked at.

The four migrated processors call `finalize_layer(LAYER_ID, version=version)` (`scripts/utils/finalize.py`) straight after `publish_processed_layer`, which emits both derived artifacts: `qa/` (`generate_qa_report.py` — invariant checks + statistics, layer-generic, driven by each file's declared attrs) and `maps/` (`generate_maps.py`). Both are regenerable and ungated, so re-running either never invalidates a layer; `finalize_layer` therefore reports failures instead of raising, since the data is already published and gated by that point. Either can be re-run standalone:

```bash
python scripts/generate_qa_report.py {layer_id} [--version V] [--local-only]
python scripts/generate_maps.py     {layer_id} [--version V] [--local-only]
```

Review order: **`maps/contact_sheet.html`** first — one panel per ensemble member at full 0.5°, which is the only view that can show a spatial defect in an individual member (§11) — then the pooled maps. It is linked from the map index and included in the bundle.

For manual review, grab **`maps/maps_bundle.zip`** (~9 MB) — one object holding the whole interlinked collection, since the S3 console downloads a single file at a time. Unzip, open `index.html`. `qa/qa_report.html` is a standalone page and needs no bundle.

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

Skill definitions live in **`.claude/skills/{name}/SKILL.md` — inside this repo, so they are version-controlled.** (They previously lived in the untracked `~/.claude/skills/` and were lost when the environment changed; recreated 2026-07-28. Do not move them back out of the repo.)

| Skill | When to Use | Product |
|-------|-------------|---------|
| `/isimip-search-download` | Searching ISIMIP repository, enumerating availability, ingesting raw to S3 | TCFD/CDP |
| `/isimip-process-visualize` | Processing annualized NetCDF, publishing a layer + its QA report and maps | TCFD/CDP |
| `/isimip-extract-aggregate` | Extracting data by location/region for CSV export | TCFD/CDP |

## Key Documentation

| Document | Purpose |
|----------|---------|
| [GUARDRAILS.md](GUARDRAILS.md) | Critical rules that must never be violated |
| [STORAGE.md](STORAGE.md) | S3 layout contract, versioning, publish protocol, credentials |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log and resolutions |
