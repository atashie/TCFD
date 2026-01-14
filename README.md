# TCFD Climate Risk Analysis Pipeline

Automated discovery, download, processing, and visualization of climate datasets from [ISIMIP.org](https://www.isimip.org/) (Inter-Sectoral Impact Model Intercomparison Project) for TCFD (Task Force on Climate-related Financial Disclosures) risk assessments.

## Overview

This project provides two workflows for climate data analysis:

1. **Python Pipeline** (`isimip-pipeline/`): Automated CLI tool for searching, downloading, processing, and visualizing ISIMIP climate data
2. **Utility Scripts** (`scripts/`): Standalone scripts for map generation and data processing

## Quick Start

### Installation

```bash
# Install the pipeline package
cd isimip-pipeline
pip install -e .

# Configure (optional - create config file)
cp config.example.yaml ~/.isimip-pipeline/config.yaml
```

### Basic Usage

```bash
# Run complete pipeline
isimip-pipeline run "groundwater runoff" --name gw-runoff --keep-raw

# Generate visualization maps
python scripts/generate_maps.py

# Clean up raw data after verification
isimip-pipeline cleanup ./data/raw
```

## Project Structure

```
TCFD/
├── README.md                 # This file
├── CLAUDE.md                 # Development guide (for Claude Code)
│
├── isimip-pipeline/          # Python package
│   ├── src/isimip_pipeline/  # Source code
│   ├── tests/                # Test suite (299 tests, 100% passing)
│   └── docs/                 # Documentation
│
├── scripts/                  # Standalone utility scripts
│   ├── generate_maps.py      # Interactive map generation
│   ├── process_qg.py         # Data processing
│   └── generate_qa_report.py # QA report generation
│
├── config/                   # Configuration examples
│   └── _example_*.json       # Example selection files
│
├── data/                     # Data directory (gitignored)
│   ├── raw/                  # Downloaded NetCDF files
│   └── processed/            # Processed outputs
│
├── reports/                  # Generated reports (gitignored)
│   └── maps/                 # Interactive HTML maps
│
└── _deprecated/              # Archived legacy code
```

## CLI Commands

| Command | Description |
|---------|-------------|
| `isimip-pipeline search` | Search ISIMIP repository |
| `isimip-pipeline download` | Download datasets |
| `isimip-pipeline process` | Process raw data |
| `isimip-pipeline report` | Generate QA report |
| `isimip-pipeline run` | Complete pipeline |
| `isimip-pipeline cleanup` | Delete raw data after verification |
| `isimip-pipeline find` | Search local datasets |
| `isimip-pipeline catalog` | Manage ISIMIP catalog |

## Data Sources

This project uses data from ISIMIP (Inter-Sectoral Impact Model Intercomparison Project):

- **Simulation Rounds**: ISIMIP3b (SSP scenarios), ISIMIP2b (RCP scenarios)
- **Climate Scenarios**: SSP126, SSP370, SSP585
- **Variables**: Groundwater runoff (qg), burnt area, evapotranspiration, and more
- **Models**: Multi-model ensembles (GFDL, IPSL, MPI, MRI, UKESM)

## Documentation

- [Pipeline README](isimip-pipeline/README.md) - Detailed package documentation
- [Scripts README](scripts/README.md) - Utility scripts documentation
- [Development Guide](CLAUDE.md) - Claude Code development instructions
- [Development History](isimip-pipeline/docs/WORKFLOW_STATUS.md) - Implementation phases

## License

[Add license information]
