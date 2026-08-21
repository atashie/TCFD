# TCFD Climate Risk Analysis Pipeline

Automated discovery, download, processing, and visualization of climate datasets from [ISIMIP.org](https://www.isimip.org/) (Inter-Sectoral Impact Model Intercomparison Project) for TCFD (Task Force on Climate-related Financial Disclosures) risk assessments.

## Overview

Two **data products** are built here, with different output contracts and different tooling.
They must not be mixed — see [CLAUDE.md](CLAUDE.md).

1. **TCFD / CDP reporting** — ISIMIP data processed into annualized decadal statistics per
   [OUTPUT-SPEC.md](OUTPUT-SPEC.md), then extracted per customer site into a CSV star schema,
   a QA dashboard, and two client-facing reports.
2. **Water Risk Index** — monthly processing of six water variables into 20 value types.
   Standalone scripts only.

All tooling lives in **`scripts/`** — per-layer processors, map generation, and the whole
customer-delivery workflow. (The former `isimip-pipeline` CLI package was archived from
HEAD on 2026-08-21; recover it via git history if ever needed.)

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

```bash
# Repo venv (Python 3.9) already carries the scientific stack
source .venv/bin/activate

# Process a layer (per-layer processors are indexed in scripts/README.md)
python scripts/process_csoil_soilcarbon.py

# Verify it against the output contract
python scripts/test_shared_baseline.py data/processed/soilcarbon_csoil_annual

# Generate visualization maps
python scripts/generate_maps.py {variable} {processed_dir} {output_dir}
```

## Project Structure

```
TCFD/
├── README.md                 # This file
├── CLAUDE.md                 # Development guide (for Claude Code)
│
├── scripts/                  # Standalone scripts — both products
│   ├── process_*.py          # Layer processors (the registry maps processor → layers)
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

## Data Sources

This project uses data from ISIMIP (Inter-Sectoral Impact Model Intercomparison Project):

- **Simulation Rounds**: ISIMIP3b (SSP scenarios), ISIMIP2b (RCP scenarios); plus non-ISIMIP
  observational sources only where ISIMIP is structurally unable to carry a hazard (NOAA SPC
  tornado, World Bank/GFDRR/Arup landslide, the Nature 2026 hail deposit)
- **Climate Scenarios**: ssp126/370/585 and rcp26/45/60/85, per layer
- **Layers**: the authoritative roster is `config/layer_registry.yaml`; what each layer is
  lives in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)
- **Models**: multi-model, multi-GCM ensembles, enumerated per layer

## Documentation

| Document | Purpose |
|---|---|
| [CLAUDE.md](CLAUDE.md) | Development guide — start here |
| [OUTPUT-SPEC.md](OUTPUT-SPEC.md) | The processed-layer output contract |
| [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md) | What each shipped layer is — per-dataset facts |
| [GUARDRAILS.md](GUARDRAILS.md) | Rules that must never be violated |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log — what went wrong and what it cost |
| [ASSET-CATALOG.md](ASSET-CATALOG.md) | Customer delivery, stages 1–2 |
| [docs/reporting/](docs/reporting/README.md) | Customer delivery, stages 3–4 |
| [docs/](docs/README.md) | Living reporting reference + frozen dated decision journals |
| [Scripts README](scripts/README.md) | Standalone scripts |

## License

[Add license information]
