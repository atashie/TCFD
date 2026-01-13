# ISIMIP Climate Data Pipeline - Project Guidance

## Current Implementation Status

**All 14 development phases complete. 299/299 tests passing (100%).**

### Completed Modules

| Module | Purpose | Test Coverage |
|--------|---------|---------------|
| `cli.py` | 8 CLI commands | Full |
| `catalog.py` | ISIMIP metrics catalog | Full |
| `discovery.py` | Local dataset search | 15 tests |
| `duplicate_handler.py` | Duplicate detection | Full |
| `interactive.py` | Guided workflow | Full |
| `processing_log.py` | YAML processing history | Full |
| `alignment.py` | Multi-model data alignment | 21 tests |
| `validation.py` | Data quality checks | 31 tests |
| `features.py` | Statistical extraction | Full |
| `output.py` | CF-compliant NetCDF output | Full |
| `qa_report.py` | Interactive HTML reports | Full |

### CLI Commands Available

```bash
isimip-pipeline find           # Search local processed datasets
isimip-pipeline search         # Search ISIMIP repository
isimip-pipeline interactive    # Guided 4-step workflow
isimip-pipeline download       # Batch download with temp management
isimip-pipeline process        # Process raw data
isimip-pipeline report         # Generate QA report
isimip-pipeline run            # Complete pipeline
isimip-pipeline catalog        # Manage ISIMIP catalog
```

---

## Project Overview

This project automates the discovery, download, processing, and visualization of climate datasets from ISIMIP.org (Inter-Sectoral Impact Model Intercomparison Project). It translates an existing manual R workflow into a robust, automated Python CLI application.

## Key Goals

1. **Natural Language Search**: Accept queries like "drought exposure metrics" and find relevant ISIMIP datasets
2. **Automated Download**: Batch download NetCDF files with resume capability
3. **Data Processing**: Extract statistical features (decadal means, trends, significance, percentiles)
4. **Quality Assurance**: Validate data alignment and detect anomalies
5. **Visualization**: Generate interactive HTML reports for QA review
6. **Standardized Output**: Produce processed NetCDF files in a consistent format

## Technology Stack

- **Language**: Python 3.10+
- **CLI Framework**: Typer
- **NetCDF Handling**: xarray, netCDF4
- **Statistics**: scipy, pandas
- **Visualization**: Plotly
- **ISIMIP API**: isimip-client library
- **LLM Integration**: you.com Agent API

## Project Structure

```
isimip-pipeline/
├── pyproject.toml              # Package config with CLI entry points
├── config.example.yaml         # Template configuration
├── src/
│   └── isimip_pipeline/
│       ├── cli.py              # CLI entry points
│       ├── config.py           # Config loading
│       ├── search/             # LLM + ISIMIP API integration
│       ├── download/           # Batch download module
│       ├── processing/         # Data alignment, validation, features
│       └── visualization/      # Interactive HTML reports
└── outputs/
```

## Configuration

Configuration is stored in `~/.isimip-pipeline/config.yaml`:
- `paths.download_dir`: Raw NetCDF storage
- `paths.processed_dir`: Processed output location
- `api.you_api_key`: you.com API key
- `api.you_agent_id`: LLM agent ID for query parsing

## CLI Commands

```bash
isimip-pipeline find "<query>"                # Search local processed datasets
isimip-pipeline search "<query>"              # Search ISIMIP repository
isimip-pipeline interactive "<query>"         # Guided discovery workflow
isimip-pipeline download --selection <file>   # Download datasets
isimip-pipeline process <input_dir>           # Process raw data
isimip-pipeline report <processed_dir>        # Generate QA report
isimip-pipeline run "<query>"                 # Full pipeline
isimip-pipeline catalog --update              # Update ISIMIP catalog
```

## Reference Files

- **WORKFLOW_STATUS.md**: Detailed implementation plan and progress tracking
- **TROUBLESHOOTING.md**: Error log and solutions encountered during development

## ISIMIP Data Context

### Simulation Rounds
- ISIMIP2a/2b: RCP scenarios (rcp26, rcp60, rcp85)
- ISIMIP3a/3b: SSP scenarios (ssp126, ssp370, ssp585)

### Key Variables
- `led`: Land area exposed to drought
- `leh`: Land area exposed to heatwave
- `lew`: Land area exposed to wildfire
- `burntarea`: Fire burnt area fraction
- `potevap`: Potential evapotranspiration

### Data Structure
NetCDFs organized by: simulation_round / model / scenario / variable / timestep

## Processing Pipeline (from R workflow)

1. Load multi-model ensembles (GFDL, HadGEM, IPSL, MIROC, etc.)
2. Apply scalar conversions (e.g., fraction to percentage)
3. Aggregate monthly to yearly
4. Smooth with kernel smoother (bandwidth ~15 years)
5. Calculate decadal statistics:
   - Median of smoothed values across models
   - Theil-Sen trend slope
   - Spearman correlation p-value
   - IQR-based confidence bounds
6. Compute percentile ranks against baseline distribution
7. Output 5D array: (lon, lat, decade, scenario, value_class)

## Value Classes in Output

1. Absolute smoothed median value
2. Percentile rank (1-100)
3. Decadal trend (slope)
4. Trend significance (p-value)
5. Lower confidence bound
6. Upper confidence bound

## Development Workflow

### Test-Driven Development (TDD)

This project uses **Test-Driven Development**. Before implementing any feature or bugfix, invoke the TDD skill:

```
/superpowers:test-driven-development
```

The TDD workflow:
1. Write failing tests that define expected behavior
2. Implement minimal code to make tests pass
3. Refactor while keeping tests green

This ensures robust, well-tested code and catches regressions early.

### Other Development Skills

- `/superpowers:systematic-debugging` - Use when encountering bugs or test failures
- `/superpowers:verification-before-completion` - Use before claiming work is done
- `/superpowers:writing-plans` - Use for planning multi-step implementations
- `/superpowers:requesting-code-review` - Use after completing features

## Development Notes

- Prefer xarray for NetCDF operations (lazy loading, dask integration)
- Use httpx for async downloads with progress tracking
- Match existing R output format for validation against reference data
- Interactive maps should support decade/scenario selection
- **Always write tests first** - invoke TDD skill before implementing new modules
