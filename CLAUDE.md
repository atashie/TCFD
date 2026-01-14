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
├── config/                   # Example configuration files
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

## ISIMIP Data Context

**Simulation Rounds:**
- ISIMIP3a/3b: SSP scenarios (ssp126, ssp370, ssp585)
- ISIMIP2a/2b: RCP scenarios (rcp26, rcp60, rcp85)

**Key Variables:** `led` (drought), `leh` (heatwave), `lew` (wildfire), `qg` (groundwater runoff), `burntarea`, `potevap`

## Search Workflow Guidelines

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

## Archived Code

Legacy R code and deprecated files are in `_deprecated/`. This directory is tracked in git for historical reference. Contents include:
- Original R-based TCFD analysis workflows
- Early Python data converters
- Historical download lists

Do not modify files in `_deprecated/` unless specifically restoring functionality.
