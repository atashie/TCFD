# TCFD Climate Risk Analysis Pipeline

Automated discovery, download, processing, and visualization of climate datasets from [ISIMIP.org](https://www.isimip.org/) (Inter-Sectoral Impact Model Intercomparison Project) for TCFD (Task Force on Climate-related Financial Disclosures) risk assessments.

## Overview

Two **data products** are built here, with different output contracts and different tooling.
They must not be mixed — see [CLAUDE.md](CLAUDE.md).

1. **TCFD / CDP reporting** — ISIMIP data processed into annualized decadal statistics per
   [OUTPUT-SPEC.md](OUTPUT-SPEC.md), then extracted per customer site into a CSV star schema,
   a QA dashboard, and two client-facing reports.
2. **Water Risk Index** — monthly processing of six water variables into 20 value types.
   Standalone scripts only; it does **not** use the `isimip-pipeline` CLI.

And two **layers of tooling**:

- **Python pipeline** (`isimip-pipeline/`) — CLI for searching, downloading and processing.
- **Standalone scripts** (`scripts/`) — per-layer processors, map generation, and the whole
  customer-delivery workflow.

## Customer delivery

Turning processed layers into something a client can file. The
**`/customer-delivery` skill is the entry point**; the commands below are what it drives.

```bash
# 1. Plan — resolves each asset to its hazard layers and STOPS. Touches no data.
python scripts/generate_customer_delivery.py --customer "Acme" --input sites.csv

# 2. Extract + dashboard + caveats + compliance report
python scripts/generate_customer_delivery.py --customer "Acme" --input sites.csv --run --reports

# 3. Bespoke report — needs facet profiles chosen and a narrative written
python scripts/generate_bespoke_report.py deliveries/acme/<date> --scaffold
python scripts/generate_bespoke_report.py deliveries/acme/<date>

# Always
python scripts/test_customer_delivery.py deliveries/acme/<date>
```

Stages run `inputs → extract → dashboard → caveats → compliance report → bespoke report`.
**Caveats runs before the reports** because the caveat set is an input to them: each report
must carry every must-disclose entry, and the verifier re-checks it. Reference documentation
is [ASSET-CATALOG.md](ASSET-CATALOG.md) for stages 1–2 and
[docs/reporting/](docs/reporting/README.md) for stages 3–4.

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
├── scripts/                  # Standalone scripts — both products
│   ├── process_*.py          # Per-layer processors (one per shipped layer)
│   ├── generate_maps.py      # Gridded interactive maps
│   ├── generate_customer_delivery.py   # Customer delivery driver
│   ├── generate_delivery_dashboard.py  # Per-delivery QA dashboard
│   ├── generate_delivery_caveats.py    # Stage 4 — caveat set
│   ├── generate_compliance_report.py   # Stage 3a — IFRS S2 spine
│   ├── generate_bespoke_report.py      # Stage 3b — facet-composed
│   ├── test_customer_delivery.py       # Delivery verifier
│   ├── test_shared_baseline.py         # Layer contract verifier
│   └── utils/                # delivery, report_common, report_figures, viz_common, …
│
├── docs/reporting/           # Stage 3–4 guidance + the facet profile library
├── config/                   # layer_registry, asset_catalog, hazard_taxonomy, search catalog
├── location-analyses/        # Example customer site inputs
├── deliveries/               # Customer deliveries (gitignored)
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
- **Variables**: Groundwater runoff (qg), drought exposure (led), burnt area, evapotranspiration, and more
- **Models**: Multi-model ensembles (GFDL, IPSL, MPI, MRI, UKESM)

## Documentation

| Document | Purpose |
|---|---|
| [CLAUDE.md](CLAUDE.md) | Development guide — start here |
| [OUTPUT-SPEC.md](OUTPUT-SPEC.md) | The processed-layer output contract |
| [GUARDRAILS.md](GUARDRAILS.md) | Rules that must never be violated |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log — what went wrong and what it cost |
| [ASSET-CATALOG.md](ASSET-CATALOG.md) | Customer delivery, stages 1–2 |
| [docs/reporting/](docs/reporting/README.md) | Customer delivery, stages 3–4 |
| [Pipeline README](isimip-pipeline/README.md) | CLI package documentation |
| [Scripts README](scripts/README.md) | Standalone scripts |

## License

[Add license information]
