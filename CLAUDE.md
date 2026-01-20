# TCFD Climate Pipeline - Development Guide

Development guidance for Claude Code when working on this project.

## Project Structure

```
TCFD/
├── README.md                 # Project overview
├── CLAUDE.md                 # This file (development guide)
│
├── isimip-pipeline/          # Python package (CLI tool)
│   ├── src/isimip_pipeline/  # Source code
│   │   ├── cli.py            # 9 CLI commands
│   │   ├── search/           # LLM + ISIMIP API
│   │   ├── download/         # Batch downloads
│   │   ├── processing/       # Data alignment, validation, features
│   │   └── visualization/    # HTML reports
│   ├── tests/                # 299 tests (100% passing)
│   └── docs/                 # Documentation + dev history
│
├── scripts/                  # Standalone utility scripts
│   ├── generate_maps.py      # Interactive Plotly maps
│   ├── process_qg.py         # QG data processing
│   └── generate_qa_report.py # QA reports
│
├── config/                   # Configuration files
│   └── isimip_search_catalog.yaml  # ISIMIP search results cache
├── data/                     # Data (gitignored)
│   ├── raw/                  # Downloaded NetCDF
│   └── processed/            # Processed outputs
├── reports/                  # Generated reports (gitignored)
│
└── _deprecated/              # Archived legacy code (tracked in git)
```

## CLI Quick Reference

```bash
isimip-pipeline find           # Search local processed datasets
isimip-pipeline search         # Search ISIMIP repository
isimip-pipeline interactive    # Guided 4-step workflow
isimip-pipeline download       # Batch download with temp management
isimip-pipeline process        # Process raw data
isimip-pipeline report         # Generate QA report
isimip-pipeline run            # Complete pipeline
isimip-pipeline catalog        # Manage ISIMIP catalog
isimip-pipeline cleanup        # Delete raw data after verification
```

## Recommended Workflow

```bash
# 1. Search and download data
isimip-pipeline run "groundwater runoff" --name gw-runoff --keep-raw

# 2. Review QA report in browser
#    Open: ./data/processed/ outputs

# 3. Generate visualization maps
python scripts/generate_maps.py

# 4. After verifying processed data is correct, delete raw files
isimip-pipeline cleanup ./data/raw
```

## Development Skills

Always invoke relevant skills before starting work:

| Skill | When to Use |
|-------|-------------|
| `/superpowers:test-driven-development` | Before implementing features or bugfixes |
| `/superpowers:systematic-debugging` | When encountering bugs or test failures |
| `/superpowers:verification-before-completion` | Before claiming work is done |
| `/superpowers:writing-plans` | For multi-step implementations |
| `/process-metrics` | When processing ISIMIP data |

## Data Processing Parameters

When processing ISIMIP data, use these parameters:

- **Temporal binning**: 2010s-2090s (no data before 2010 or after 2099)
- **Adaptive windowing**: Minimum 100 data points per decade-bin
- **Percentile baseline**: Always use 2020s as reference distribution
- **Sub-annual data**: Always ask user for aggregation method. See GUARDRAILS.md §2.

### Shared 2020s Baseline

All climate projections share **identical 2020s baseline values**, computed as:
- Average of ALL scenarios (projection + historical/picontrol) that have overlapping data for the 2020s period
- For 2030s-2090s decades, values are computed per-scenario as before

This ensures consistent baseline comparisons across scenarios while allowing divergent projections for future decades. The `baseline_source` attribute in processed NetCDF files indicates whether shared baseline was used.

**Verification**: Run `python scripts/test_shared_baseline.py` after processing to verify:
1. 2020s values are identical across all scenarios
2. 2030s+ values differ across scenarios (as expected)
3. Files have `baseline_source: shared_across_all_scenarios` attribute

### Shared Baseline Implementation

Processing scripts must implement shared 2020s baseline following this pattern:

1. Collect 2020s data from ALL available scenarios (projection + historical/picontrol)
2. Average across scenarios per GCM to create shared baseline
3. Use shared baseline for ALL scenarios' 2020s decade
4. Use scenario-specific data for 2030s-2090s
5. Set `baseline_source: "shared_across_all_scenarios"` attribute
6. Output separate files per scenario

**Reference implementations**:
- `scripts/process_qg.py` - groundwater runoff processing
- `scripts/process_timber.py` - timber/wood carbon processing

**Verification**: Run `python scripts/test_shared_baseline.py {output_dir}` after processing.

### Output Organization

**Folder naming**: `{descriptive-name}_{variable}_{timestep}/`
- Example: `wildfire-burntarea_burntarea_monthly/`
- Example: `drought-severity_led_annual/`

**Multi-scenario rule**: All climate scenarios belong in ONE folder and ONE processed file.
- Processed NetCDF dimensions: `(lon, lat, decade, scenario, value_class)`
- Never create separate folders per scenario (e.g., avoid `burntarea-rcp26/`, `burntarea-rcp60/`)

**CLI workflow**: Always use `isimip-pipeline run` for downloading data.
- The `run` command automatically handles multi-scenario downloads into a single output directory
- Avoid manual `search` + `download` workflows that may fragment scenarios into separate folders

### QA Report Generation

**Always use `scripts/generate_maps.py`** for generating visualization reports.

Never write ad-hoc inline scripts for report generation. The established script provides:
- Per-scenario HTML files (one file per metric × scenario combination)
- 2020s vs 2090s comparison structure in each file
- Cross-navigation between metrics and scenarios
- Master index page with grid layout
- Browser-safe file sizes (~4MB each)
- Anomaly detection and JSON summaries
- Auto-detection of SSP vs RCP scenarios

**Report workflow**:
```bash
# After processing data, generate maps
python scripts/generate_maps.py {variable} {processed_dir} {output_dir}

# Example for heatwave data:
python scripts/generate_maps.py leh ./outputs/heatwave-exposure_leh-annual/processed ./reports/maps
```

**Do NOT**:
- Write custom inline Python to generate HTML visualizations
- Create per-decade files instead of per-scenario files
- Generate monolithic single-file reports (>10MB)

## ISIMIP Data Context

**Simulation Rounds:**
- ISIMIP3a/3b: SSP scenarios (ssp126, ssp370, ssp585)
- ISIMIP2a/2b: RCP scenarios (rcp26, rcp60, rcp85)

**Key Variables:** `led` (drought), `leh` (heatwave), `lew` (wildfire), `qg` (groundwater runoff), `burntarea`, `potevap`, `npp`, `gpp`, `cveg`, `b30cm` (large fish biomass)

**Documentation Updates:** When new simulation rounds, variables, or data sources are discovered in ISIMIP, work with the user to update this documentation. The ISIMIP repository evolves; keep this guide current.

### Scenario Handling

**Projection scenarios** (included in reports):
- SSP: ssp126, ssp245, ssp370, ssp585
- RCP: rcp26, rcp45, rcp60, rcp85

**Scenario discovery**: Scripts must discover scenarios dynamically from filesystem, not hardcoded lists. See GUARDRAILS.md §3.

**Non-projection scenarios** (excluded from reports):
- `picontrol` (Pre-Industrial Control)
- `historical`

Non-projection scenarios may be downloaded and processed alongside projection data to enhance baseline robustness (e.g., improving 2020s reference distributions), but they are **automatically excluded** from `generate_maps.py` report generation. Only actual climate projections appear in the final visualization reports.

**New Scenarios:** If you encounter scenarios not listed above (e.g., ssp245, rcp45), discuss with the user whether to add them to this documentation and update the processing/reporting code accordingly.

### Experiment & SOC Preferences

**Sensitivity experiments:** Prefer `default` (transient CO2). Ask user if `2015co2` or others should be included.

**SOC scenarios:** Prefer `2015soc` (isolates climate signal). Ask user if `histsoc` or `2015soc-from-histsoc` should be included.

**General principle:** Favor data abundance over restrictiveness—when in doubt, include more datasets.

## Search Workflow Guidelines

### Search Catalog Reference

**Before starting any ISIMIP search**, consult `config/isimip_search_catalog.yaml` to check if the variable, model, or scenario has already been investigated. This catalog contains empirical results from previous searches including: variable definitions and units, available PFT (Plant Functional Type) codes, which impact models and climate forcings have data, temporal and spatial coverage, and scenario availability. Using this catalog avoids redundant API queries and helps contextualize new search results against known data availability. When conducting new searches that yield useful information not already in the catalog, update the file with the new findings including the search date, dataset counts, and any relevant notes about data quality or limitations.

Summary tables should include:
- Variable code and description
- Simulation round (ISIMIP2a/2b/3a/3b)
- Time step (annual, monthly, daily)
- Land surface models available
- Climate scenarios
- Dataset counts

**Principles:**
1. Show all data sources (include all simulation rounds)
2. Prefer newer data (ISIMIP3b over ISIMIP2a)
3. Emphasize multi-model ensemble availability
4. **Maximize dataset count:** Always prioritize including the maximum number of datasets/models/GCMs available. More data points improve statistical robustness and uncertainty quantification. Never settle for single-model analysis when multi-model ensembles are available.

### Search Principles

**ISIMIP First:** Always search the ISIMIP repository before considering external data sources. Push back if the user requests non-ISIMIP data unless there's a compelling reason (e.g., ISIMIP has no relevant data after thorough search).

**Global Over Regional:** Prefer global-coverage datasets over regional/local data when available:
- Global biome models (LPJmL, ORCHIDEE, VISIT, CLM) over site-specific models
- This ensures consistent methodology across geographies

**No Direct Match Strategy:** When exact data isn't available in ISIMIP:
1. List close-enough proxy variables (e.g., NPP/GPP as proxy for growth, vegetation carbon for biomass)
2. Identify relevant plant functional types (e.g., "temperate broadleaf deciduous" for oak)
3. Suggest which ISIMIP sectors/models might contain usable approximations
4. Only after exhausting ISIMIP options, discuss external alternatives with user approval

## Archived Code

Legacy R code and deprecated files are in `_deprecated/`. This directory is tracked in git for historical reference. Contents include:
- Original R-based TCFD analysis workflows
- Early Python data converters
- Historical download lists

Do not modify files in `_deprecated/` unless specifically restoring functionality.

See [GUARDRAILS.md](GUARDRAILS.md) for critical rules that must never be violated.
