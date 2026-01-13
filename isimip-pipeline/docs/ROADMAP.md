# ISIMIP Pipeline Roadmap

## Current Status (v0.4.0 - Feature Complete)

**All 14 development phases complete. 299/299 tests passing (100%).**

The pipeline successfully:
- Searches ISIMIP repository via keyword fallback or LLM parsing
- Discovers and searches local processed datasets
- Downloads NetCDF files with temp file management and resume support
- Aligns multi-model data (spatial, temporal, calendar conversion)
- Validates data quality (fill values, gaps, outliers)
- Processes data into standardized 6-class output format
- Generates interactive HTML reports with Plotly maps
- Tracks processing history with duplicate detection

## Known Issues

### P1: Critical

#### 1. Percentile Ranking Shows Only Extreme Values

**Symptom**: Percentile maps display only values near 1 or 100, with no intermediate values.

**Root Cause Analysis**:
- Currently comparing smoothed median values against raw historical distribution
- The smoothing operation changes the value distribution significantly
- Ocean/NaN cells may be included in baseline calculation, skewing percentiles
- Using `percentileofscore` with `kind="weak"` may not be appropriate for smoothed data

**Investigation Steps**:
1. Add diagnostic logging to show distribution of percentile values
2. Verify historical baseline excludes NaN/fill values correctly
3. Compare smoothed vs raw value distributions
4. Test with land-only mask applied

**Proposed Fix**:
```python
# Option A: Compare smoothed against smoothed historical
smoothed_hist = apply_smoothing(historical_data)
percentile = percentileofscore(smoothed_hist.flatten(), smoothed_median, kind="rank")

# Option B: Use quantile-based mapping
hist_quantiles = np.nanpercentile(historical_data, np.arange(1, 101))
percentile = np.searchsorted(hist_quantiles, smoothed_median)
```

### P2: High Priority

#### 2. Raw Files Stored in Git Directory

**Current**: Downloads saved to `./data/raw` or output directory within repo.

**Issue**: Large NetCDF files should not be tracked in git.

**Proposed Solution**:
- Default download to system temp directory: `tempfile.mkdtemp(prefix="isimip_")`
- Auto-cleanup after processing completes
- Add `--keep-raw` flag to preserve downloads if needed
- Store only metadata about downloaded files

```python
# In downloader.py
import tempfile
import shutil

class Downloader:
    def __init__(self, keep_raw: bool = False):
        self.temp_dir = tempfile.mkdtemp(prefix="isimip_")
        self.keep_raw = keep_raw

    def cleanup(self):
        if not self.keep_raw and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
```

#### 3. Limited Metadata Tracking

**Current**: Only `selection.json` with dataset URLs.

**Proposed Enhancement**: Build persistent ISIMIP catalog.

**Schema** (`~/.isimip-pipeline/catalog.yaml`):
```yaml
last_updated: "2024-01-15T10:30:00Z"
variables:
  led:
    long_name: "Land Exposed to Drought"
    unit: "fraction"
    simulation_rounds: [ISIMIP2b, ISIMIP3b]
    scenarios: [historical, rcp26, rcp60, rcp85, ssp126, ssp370, ssp585]
    models: [gfdl-esm2m, hadgem2-es, ...]
    file_count: 42
    last_seen: "2024-01-15"
  burntarea:
    # ...
```

**Implementation**:
- Update catalog on every search/download
- Merge new findings with existing entries
- Track file counts, scenarios, and models per variable
- Export catalog statistics for user review

### P3: Medium Priority

#### 4. No Anomaly Detection in Processing

**Current**: Processing continues silently even with suspect data.

**Proposed QA Checks**:
- Flag cells with values outside expected physical bounds
- Detect spatial discontinuities (edge artifacts, model seams)
- Identify temporal jumps within decades
- Report percentage of valid/invalid cells per output

**Implementation in `validation.py`**:
```python
def detect_anomalies(data: xr.DataArray, variable: str) -> AnomalyReport:
    bounds = VARIABLE_BOUNDS.get(variable, (-inf, inf))

    out_of_bounds = (data < bounds[0]) | (data > bounds[1])
    spatial_gradient = compute_gradient_magnitude(data)
    high_gradient = spatial_gradient > threshold

    return AnomalyReport(
        out_of_bounds_pct=out_of_bounds.mean() * 100,
        high_gradient_cells=high_gradient.sum(),
        coverage_pct=(~np.isnan(data)).mean() * 100,
    )
```

#### 5. Missing Processed Data Documentation

**Current**: Processed NetCDF files have minimal attributes.

**Proposed**: Generate companion `.txt` metadata files.

**Content** (`led_processed.txt`):
```
ISIMIP Pipeline Processing Report
Generated: 2024-01-15 10:30:00 UTC
Variable: led (Land Exposed to Drought)

INPUT FILES (3)
  - lpjml_gfdl-esm4_historical_led_global_annual_1850_2014.nc
  - lpjml_gfdl-esm4_ssp126_led_global_annual_2015_2100.nc
  - lpjml_gfdl-esm4_ssp585_led_global_annual_2015_2100.nc

PROCESSING PARAMETERS
  Smoothing bandwidth: 15 years
  Trend method: Theil-Sen
  Significance method: Spearman correlation

OUTPUT DIMENSIONS
  Longitude: -180 to 180 (0.5° resolution, 720 cells)
  Latitude: -90 to 90 (0.5° resolution, 360 cells)
  Decades: [10, 20, 30, 40, 50, 60, 70, 80, 90]
  Scenarios: [historical, ssp126, ssp585]
  Value Classes: [smoothed_median, percentile, trend, significance,
                  lower_bound, upper_bound]

DATA QUALITY
  Valid cells: 67,234 / 259,200 (25.9%)
  Ocean/land mask: Inferred from data
  Fill value handling: 1e+20 → NaN

VALUE RANGES
  smoothed_median: 0.0 to 0.87 (mean: 0.12)
  percentile: 1 to 100 (mean: 48.2)
  trend: -0.05 to 0.08 per decade (mean: 0.01)
  significance: 0.001 to 0.95 (32% cells p<0.05)

ANOMALY FLAGS
  - 0 cells outside physical bounds (0-1 for fractions)
  - 142 cells with high spatial gradient (potential artifacts)
  - Coverage consistent across decades
```

## Feature Roadmap

### v0.2.0: Data Quality Focus ✅ COMPLETE

1. ✅ **Fix percentile calculation** - Implemented in features.py
2. ✅ **Temp file management** - Auto-cleanup with `--keep-raw` flag
3. ✅ **ISIMIP catalog** - Persistent metrics catalog in catalog.py
4. ✅ **Metadata generation** - CF-compliant NetCDF attributes in output.py

### v0.3.0: Validation & Robustness ✅ COMPLETE

1. ✅ **Anomaly detection** - Implemented in validation.py (outlier detection, fill values)
2. ✅ **Download verification** - File integrity checks in downloader.py
3. ✅ **Land/ocean masking** - Handled via NaN/fill value masking
4. ✅ **Fill value standardization** - Multi-format detection (1e20, -9999, NaN)

### v0.4.0: Conversational Workflow ✅ COMPLETE

The pipeline supports interactive refinement with human-in-the-loop processing.

**Implemented in interactive.py**:
```
[Search] → User reviews datasets → [Select]
    ↓
[Download] → User monitors progress → [Verify]
    ↓
[Process] → User reviews QA flags → [Accept/Reject]
    ↓
[Report] → User explores visualizations → [Export/Refine]
```

**Available Features**:
- ✅ `isimip-pipeline interactive` - Guided 4-step workflow
- ✅ `isimip-pipeline find` - Search local processed datasets
- ✅ Duplicate detection - Prevents redundant processing
- ✅ Processing log - Tracks all datasets in YAML format
- ✅ Natural language search - Via LLM query parsing

### v0.5.0: Advanced Analytics

1. **Multi-variable comparison** - Cross-correlation analysis
2. **Ensemble statistics** - Model agreement maps
3. **Regional aggregation** - Country/basin level summaries
4. **Time series export** - Location-based trajectories

## Architecture Improvements

### Current Pain Points

1. **Processor.py is too large** (~400 lines) - Split into focused modules
2. **No caching layer** - Reprocesses everything on each run
3. **Hard-coded parameters** - Should be configurable per variable
4. **Limited error recovery** - One failure stops entire pipeline

### Proposed Refactoring

```
processing/
├── processor.py       # Orchestration only
├── loader.py          # NetCDF loading, alignment
├── features.py        # Statistical calculations (existing)
├── validation.py      # Anomaly detection, QA
├── output.py          # NetCDF writing (existing)
└── cache.py           # Processed data caching
```

## Testing Strategy

### Current: 299 tests passing (100% pass rate)

Comprehensive test coverage including:
- **Unit tests** - All modules with isolated function testing
- **Integration tests** - Full pipeline workflows
- **Validation tests** - Fill values, gaps, outliers (31 tests)
- **Alignment tests** - Multi-model data harmonization (21 tests)
- **Discovery tests** - Local dataset search (15 tests)

### Future Enhancements:
- **Regression tests** - Compare outputs to reference R code
- **Performance tests** - Memory/time benchmarks for large datasets
- **Visual regression** - Snapshot testing for reports

## Contributing

See `CONTRIBUTING.md` for development setup and guidelines.

## Version History

- **v0.1.0** (Current): Initial working pipeline
  - Search, download, process, visualize
  - 6-class output format
  - Interactive HTML reports

## Next Steps

**v0.5.0 Goals (Advanced Analytics)**:
1. Multi-variable cross-correlation analysis
2. Model ensemble agreement maps
3. Regional aggregation (country/basin summaries)
4. Time series export for location-based analysis

**Documentation**:
- See `WORKFLOW_STATUS.md` for implementation details
- See `TROUBLESHOOTING.md` for common issues and solutions
- See `README.md` for CLI usage examples
