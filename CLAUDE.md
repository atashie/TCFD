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
│   ├── process_*.py          # Data processing scripts
│   └── test_shared_baseline.py # Baseline verification
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

## Development Skills

Always invoke relevant skills before starting work:

| Skill | When to Use |
|-------|-------------|
| `/isimip-search-download` | Searching ISIMIP repository, downloading data |
| `/isimip-process-visualize` | Processing NetCDF files, generating QA reports |
| `/superpowers:test-driven-development` | Before implementing features or bugfixes |
| `/superpowers:systematic-debugging` | When encountering bugs or test failures |
| `/superpowers:verification-before-completion` | Before claiming work is done |

## Key Documentation

| Document | Purpose |
|----------|---------|
| [GUARDRAILS.md](GUARDRAILS.md) | Critical rules that must never be violated |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log and resolutions |
| [config/isimip_search_catalog.yaml](config/isimip_search_catalog.yaml) | ISIMIP search results cache |

## Archived Code

Legacy R code and deprecated files are in `_deprecated/`. This directory is tracked in git for historical reference. Do not modify unless specifically restoring functionality.
