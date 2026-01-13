# ISIMIP Pipeline

Automated discovery, download, processing, and visualization of ISIMIP climate datasets.

## Current Status

**Version**: 0.1.0 (Development)

### Working Components
- **Search**: Natural language queries via LLM or keyword fallback
- **Download**: Batch download with resume support and progress tracking
- **Processing**: Vectorized feature extraction (median, percentile, trend, significance, bounds)
- **Visualization**: Interactive HTML reports with Plotly geographic maps

### Known Issues (v0.1.0)
1. **Percentile ranking may show extreme values** - Under investigation; may be comparing smoothed values against raw historical distribution
2. **Raw files stored locally** - Should use temp folder and cleanup after processing
3. **Limited metadata tracking** - Need persistent ISIMIP metrics catalog
4. **No anomaly detection** - Missing QA validation for processed data

## Installation

```bash
pip install -e ".[dev]"
```

## Configuration

Create `~/.isimip-pipeline/config.yaml`:
```yaml
paths:
  download_dir: "./data/raw"
  processed_dir: "./data/processed"
  reports_dir: "./reports"

api:
  you_api_key: "your-key-here"  # Optional: for LLM query parsing
  you_agent_id: "52a88ef1-ff18-4485-bba8-cb5977baf136"

processing:
  smoothing_bandwidth: 15
  percentile_bins: 100
```

## Usage

### Individual Commands
```bash
# Search ISIMIP repository
isimip-pipeline search "drought exposure metrics"

# Download datasets
isimip-pipeline download --selection results.json

# Process NetCDF files
isimip-pipeline process ./data/raw

# Generate QA report
isimip-pipeline report ./data/processed
```

### Full Pipeline
```bash
# Run complete pipeline
isimip-pipeline run "drought exposure" --limit 20

# With specific output directory
isimip-pipeline run "wildfire burnt area" -o ./my_output

# Skip download (reuse existing files)
isimip-pipeline run "drought exposure" --skip-download
```

## Output Structure

```
output_dir/
├── all_available_datasets.json  # Full list of matching datasets
├── selection.json               # Selected datasets for processing
├── raw/                         # Downloaded NetCDF files
├── processed/                   # Processed feature files
│   └── {variable}_processed.nc  # Standardized output
└── reports/
    └── qa_report.html           # Interactive visualization
```

## Processed Data Format

Output NetCDF dimensions: `(lon, lat, decade, scenario, value_class)`

**Value Classes:**
1. `smoothed_median` - Kernel-smoothed median value
2. `percentile` - Rank against baseline distribution (1-100)
3. `trend` - Theil-Sen slope estimate
4. `significance` - Spearman correlation p-value
5. `lower_bound` - 25th percentile bound
6. `upper_bound` - 75th percentile bound

## Report Structure

The HTML QA report includes:
1. **Summary Statistics** - Data coverage, valid cell counts
2. **Baseline vs End-of-Century Comparison** - 2010s vs 2090s for all value classes
3. **Scenario Comparison** - 2090s median across climate scenarios (if multiple)

## Development

This project uses Test-Driven Development.

```bash
# Run tests
pytest

# Run with coverage
pytest --cov=isimip_pipeline

# Run specific test file
pytest tests/test_processor.py -v
```

## Architecture

```
src/isimip_pipeline/
├── cli.py                    # CLI entry points (Typer)
├── config.py                 # Configuration loading
├── search/
│   ├── llm_parser.py         # you.com API integration
│   ├── isimip_query.py       # ISIMIP API wrapper
│   └── result_table.py       # Results display
├── download/
│   └── downloader.py         # Async batch downloads
├── processing/
│   ├── processor.py          # Main processing logic
│   ├── features.py           # Statistical feature extraction
│   └── output.py             # NetCDF output generation
└── visualization/
    └── qa_report.py          # HTML report generation
```

## Roadmap

See `docs/ROADMAP.md` for planned improvements including:
- Percentile calculation fixes
- Temp file management
- ISIMIP metrics catalog
- Anomaly detection
- Conversational workflow support
