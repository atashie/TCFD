# ISIMIP Climate Data Pipeline

Automated discovery, download, processing, and visualization of ISIMIP (Inter-Sectoral Impact Model Intercomparison Project) climate datasets.

## Features

- **Natural Language Search**: Query ISIMIP datasets using plain English (e.g., "drought exposure metrics")
- **Local Dataset Discovery**: Track and search previously processed datasets
- **Interactive Workflow**: Guided 6-step process from discovery to visualization
- **Duplicate Detection**: Prevent redundant processing based on variable+timestep
- **Multi-Model Processing**: Ensemble statistics across climate models
- **Statistical Features**: Decadal means, Theil-Sen trends, percentile ranks, confidence intervals
- **Quality Assurance**: Interactive HTML reports with geographic visualizations

## Installation

```bash
# Clone the repository
git clone https://github.com/your-repo/isimip-pipeline.git
cd isimip-pipeline

# Install in development mode
pip install -e ".[dev]"

# Verify installation
isimip-pipeline --version
```

### Requirements

- Python 3.10+
- Dependencies: typer, xarray, scipy, plotly, httpx, isimip-client

## Quick Start

### 1. Configure the Pipeline

Create `~/.isimip-pipeline/config.yaml`:

```yaml
paths:
  download_dir: "./data/raw"
  processed_dir: "./data/processed"
  reports_dir: "./reports"

api:
  you_api_key: "your-key-here"  # Optional: for LLM query parsing
  you_agent_id: "agent-id"       # Optional: you.com agent ID

processing:
  smoothing_bandwidth: 15        # Kernel smoothing bandwidth (years)
  percentile_bins: 100           # Percentile ranking bins
```

### 2. Search Local Datasets

```bash
# List all processed datasets
isimip-pipeline find

# Search by keyword
isimip-pipeline find drought

# Filter by variable and timestep
isimip-pipeline find -v led -t monthly

# Show detailed information
isimip-pipeline find drought --detailed
```

### 3. Interactive Workflow (Recommended)

```bash
# Start guided dataset discovery
isimip-pipeline interactive "drought exposure metrics"
```

This launches a 4-step interactive process:
1. **Local Search** - Check for existing datasets
2. **Remote Search** - Query ISIMIP API if not found locally
3. **User Selection** - Choose from grouped results (by variable+timestep)
4. **Duplicate Check** - Handle existing datasets

### 4. Download and Process

```bash
# Download selected datasets
isimip-pipeline download -s outputs/drought-severity_led-monthly/selection.json

# Process downloaded data
isimip-pipeline process outputs/drought-severity_led-monthly/raw
```

### 5. Full Pipeline (One Command)

```bash
# Run complete workflow end-to-end
isimip-pipeline run "drought exposure metrics" -o outputs/drought-data
```

## CLI Commands

| Command | Description | Example |
|---------|-------------|---------|
| `find` | Search local processed datasets | `isimip-pipeline find -v led -t monthly` |
| `interactive` | Guided discovery workflow | `isimip-pipeline interactive "drought"` |
| `search` | Search ISIMIP repository | `isimip-pipeline search "wildfire"` |
| `download` | Batch download datasets | `isimip-pipeline download -s selection.json` |
| `process` | Process NetCDF files | `isimip-pipeline process ./raw` |
| `report` | Generate QA report | `isimip-pipeline report ./processed` |
| `run` | Complete pipeline | `isimip-pipeline run "query" -o ./output` |
| `catalog` | Manage ISIMIP catalog | `isimip-pipeline catalog --update` |

## Output Structure

```
outputs/{name}_{variable}-{timestep}/
├── selection.json               # Dataset selection metadata
├── raw/                         # Downloaded NetCDF files
│   └── *.nc
├── processed/                   # Processed feature files
│   └── {variable}_processed.nc  # Standardized output
└── reports/
    └── qa_report.html           # Interactive visualization
```

### Folder Naming Convention

Folders are named `{descriptive-name}_{variable}-{timestep}`:
- `drought-severity_led-monthly`
- `fire-risk_burntarea-annual`
- `heatwave-exposure_leh-daily`

## Processed Data Format

Output NetCDF files have dimensions: `(lon, lat, decade, scenario, value_class)`

### Value Classes

| Index | Name | Description |
|-------|------|-------------|
| 1 | `smoothed_median` | Kernel-smoothed ensemble median |
| 2 | `percentile` | Rank against baseline (1-100) |
| 3 | `trend` | Theil-Sen slope estimate |
| 4 | `significance` | Spearman correlation p-value |
| 5 | `lower_bound` | 25th percentile bound |
| 6 | `upper_bound` | 75th percentile bound |

### Decades

Data is aggregated into 9 decades: 2010s, 2020s, ..., 2090s

## Processing Log

All processed datasets are tracked in `outputs/processed_data_log.yaml`:

```yaml
datasets:
  - descriptive_name: "drought-severity"
    variable: "led"
    timestep: "monthly"
    created_date: "2026-01-13T10:30:00"
    output_path: "./outputs/drought-severity_led-monthly"
    file_count: 15
    time_periods: ["2006-2100"]
    climate_scenarios: ["ssp126", "ssp370", "ssp585"]
    gcm_models: ["gfdl-esm4", "ukesm1-0-ll"]
    simulation_round: "ISIMIP3b"
```

## Architecture

```
src/isimip_pipeline/
├── cli.py                    # CLI entry points (8 commands)
├── config.py                 # Configuration loading
├── catalog.py                # ISIMIP metrics catalog
├── discovery.py              # Local dataset discovery
├── duplicate_handler.py      # Duplicate detection
├── interactive.py            # Interactive workflow
├── processing_log.py         # Processing log management
├── search/
│   ├── isimip_query.py       # ISIMIP API wrapper
│   ├── llm_parser.py         # LLM query parsing
│   └── result_table.py       # Results display & grouping
├── download/
│   └── downloader.py         # Async batch downloads
├── processing/
│   ├── processor.py          # Main processing logic
│   ├── features.py           # Statistical feature extraction
│   ├── alignment.py          # Multi-model data alignment
│   ├── validation.py         # Data quality validation
│   └── output.py             # NetCDF output generation
└── visualization/
    └── qa_report.py          # HTML report generation
```

## Data Validation

The pipeline validates input data for:
- **Fill Values**: Detects standard (1e20, -9999) and custom patterns
- **Spatial Gaps**: Identifies missing grid cells
- **Temporal Gaps**: Finds periods with data loss
- **Outliers**: IQR and z-score detection methods

Quality levels: `EXCELLENT` (<1% loss) | `GOOD` (1-5%) | `ACCEPTABLE` (5-10%) | `POOR` (10-25%) | `UNUSABLE` (>25%)

## Development

This project uses Test-Driven Development (TDD).

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_processor.py -v

# Run with coverage
pytest --cov=isimip_pipeline

# Run integration tests only
pytest tests/test_integration.py -v
```

### Test Coverage

- 299 total tests
- 299 passing (100% pass rate)
- Covers: search, download, processing, alignment, validation, output, visualization

## ISIMIP Data Context

### Simulation Rounds
- **ISIMIP2a/2b**: RCP scenarios (rcp26, rcp60, rcp85)
- **ISIMIP3a/3b**: SSP scenarios (ssp126, ssp370, ssp585)

### Common Variables
| Variable | Description |
|----------|-------------|
| `led` | Land area exposed to drought |
| `leh` | Land area exposed to heatwave |
| `lew` | Land area exposed to wildfire |
| `burntarea` | Fire burnt area fraction |
| `potevap` | Potential evapotranspiration |

## License

MIT License

## Contributing

1. Fork the repository
2. Create a feature branch
3. Write tests first (TDD)
4. Implement the feature
5. Submit a pull request

See `WORKFLOW_STATUS.md` for implementation details and roadmap.
