# Workflow Issues Log

This document tracks workflow mistakes and their resolutions. Each incident documents what went wrong, why it matters, and how to prevent recurrence.

See [GUARDRAILS.md](GUARDRAILS.md) for the rules derived from these incidents.

---

## Incident Log

### 2026-01-16: Fish TCB Downloaded Without Resolution Choice

**What happened**: When searching for fish catch abundance data, I found both monthly (~135-158 MB/file) and annual (~8-12 MB/file) data available. I proceeded to download 28 monthly files (3.62 GB total) without asking the user which resolution they preferred.

**Impact**: Downloaded 12x more data than potentially needed if annual resolution was sufficient.

**Available options that should have been presented**:

| Resolution | Files | Size per file | Total |
|------------|-------|---------------|-------|
| Monthly | 24 files | 135-158 MB | ~3.5 GB |
| Annual | 4 files (ZooMSS only) | 8-12 MB | ~40 MB |

**Correct action**: Should have used `AskUserQuestion` to present both options and let user decide based on their analysis needs.

**Rule created**: GUARDRAILS.md Section 1 - Temporal Resolution Must Be User-Selected

---

### 2026-01-20: Fish b30cm Processed Without Aggregation Choice

**What happened**: When processing b30cm (large fish biomass) monthly data, I automatically aggregated monthly values to annual using `groupby("time.year").mean()` without asking the user which aggregation method they preferred.

**Impact**: Processing proceeded with a reasonable method (mean for density data), but user was not consulted. In other cases (e.g., count data, extremes), this could have produced wrong results.

**What I did**:
- Monthly data to Annual: arithmetic mean of 12 months
- Annual to Decade: median across years in adaptive window
- Models to Ensemble: median across model outputs

**Correct action**: Should have:
1. Noted the units (g C m^-2) are a density metric
2. Presented aggregation options: mean, median, sum, min, max
3. Used `AskUserQuestion` to get explicit user choice
4. In this case, user confirmed mean was correct for density units

**Outcome**: User approved mean aggregation post-hoc, but workflow was incorrect.

**Rule created**: GUARDRAILS.md Section 2 - Sub-Annual Data Aggregation Requires User Choice

---

### 2026-01-20: Scenario rcp45 Excluded from Visualization

**What happened**: After processing b30cm data with all 4 RCP scenarios (rcp26, rcp45, rcp60, rcp85), the visualization script `generate_maps.py` only generated maps for 3 scenarios (rcp26, rcp60, rcp85), silently excluding rcp45.

**Root cause**: The script had hardcoded scenario lists:
```python
RCP_SCENARIO_LABELS = {
    "rcp26": "RCP2.6 (Low Emissions)",
    "rcp60": "RCP6.0 (Intermediate)",
    "rcp85": "RCP8.5 (High Emissions)",
    ...
}
```
The per-scenario file loading loop only tried scenarios in this dictionary, so rcp45 files were never discovered.

**Impact**: 25% of processed data was invisible in visualizations.

**Fix applied**: Changed `generate_maps.py` to dynamically discover scenarios from filesystem using glob pattern `{variable}_*_processed.nc` instead of iterating hardcoded lists.

**Rule created**: GUARDRAILS.md Section 3 - Scenario Discovery Must Be Dynamic

---

### 2026-01-21: Loblolly Pine Search Missed Climate-Specific PFT Datasets

**What happened**: When searching for loblolly pine (southern temperate conifer) timber data, I only found CLASSIC model's generic `evgndltr` PFT from ISIMIP3b. This PFT combines ALL evergreen needleleaf trees (temperate AND boreal) into one category, making it a poor proxy for a temperate-specific species.

**What was missed**:
- **MC2-USFS** (ISIMIP3a): `mesictemperateneedleleafforest`, `subtropicalevergreenneedleleafforest` - 47 detailed biome-specific PFTs
- **CLM45** (ISIMIP2b): `needleleaf-evergreen-tree-temperate` - climate-zone specific with RCP scenarios
- **LPJmL** (ISIMIP2b): `temperate-needleleaved-evergreen-tree` - climate-zone specific

**Impact**: Analysis using generic `evgndltr` included boreal conifer data (e.g., black spruce at 60N) mixed with target temperate species data. Climate projections for loblolly pine (southeastern US, 30-35N) were potentially confounded by boreal zone dynamics where climate change impacts differ substantially.

**Root causes**:
1. **Single-model tunnel vision**: Started with ISIMIP3b CLASSIC, assumed its PFT scheme was representative
2. **Incomplete catalog**: Only documented CLASSIC and CARAIB PFT schemes, missing CLM45, LPJmL, MC2-USFS
3. **API unreliability unrecognized**: ISIMIP API path-based searches returned unrelated flood data; didn't pivot to file server exploration
4. **No model documentation lookup**: Didn't check impactmodels.org pages for each model's PFT definitions
5. **Simulation round blindness**: Best climate-specific data was in ISIMIP2b and ISIMIP3a, not ISIMIP3b
6. **PFT not treated as searchable dimension**: Workflow asked for variable but not PFT subcategory

**Correct action**: Should have:
1. Enumerated ALL biomes sector models across simulation rounds (3b, 3a, 2b)
2. Checked each model's PFT documentation via impactmodels pages
3. Used file server (`files.isimip.org`) when API returned unreliable results
4. Presented user with climate-specificity options and explicit trade-offs
5. Asked user to select PFT subcategory before downloading

**Fix applied**:
- Added `biomes_models` registry to `config/isimip_search_catalog.yaml` documenting all models and their PFT schemes
- Added `pft_equivalences` table mapping equivalent PFT concepts across models
- Added Step 0 (Model Enumeration) to `/isimip-search-download` skill
- Added Step 2.5 (Subcategory Selection) to skill workflow
- Added API Fallback Strategy section for file server exploration
- Created GUARDRAILS.md Section 4 for multi-model vegetation searches

**Rule created**: GUARDRAILS.md Section 4 - Multi-Model Search for Vegetation Variables

---

### 2026-01-22: QA Report Map Alignment and Colorscale Issues

**What happened**: During tebrsu (temperate broadleaf summergreen) data visualization, two issues were identified:

1. **Geospatial alignment**: Map outputs don't perfectly align with basemap coastlines
2. **Trend colorscale not centered on zero**: Diverging colorscale doesn't have white=0

**Impact**: Visual interpretation issues. Trend maps may show positive trends as partially red (or negative as blue) due to asymmetric color scaling.

**Root causes**:

**Alignment issue**:
- ISIMIP data uses 0.5° × 0.5° grid with cell-center registration (coordinates at pixel centers, not edges)
- Visualization uses point markers (`go.Scattergeo` with marker size=2) at cell centers, not filled grid cells
- At 0.5° resolution (~55 km/cell), coastal cells straddle land/water but plot at single center point
- Basemap coastlines are vector data at much finer resolution than the raster grid

**Colorscale issue**:
- Trend maps use percentile-based scaling (`cmin = np.percentile(all_values, 2)`, `cmax = np.percentile(all_values, 98)`)
- This does NOT center on zero; if trends are mostly positive (e.g., 0.1 to 0.5), white falls at ~0.3 instead of 0
- Change maps correctly use symmetric scaling: `max_abs = np.percentile(np.abs(all_changes), 98); cmin, cmax = -max_abs, max_abs`

**Location**: `scripts/generate_maps.py`
- Lines 602-606: Trend colorscale (NOT symmetric)
- Lines 757-761: Change colorscale (correctly symmetric)

**Correct action for future**:

*Alignment*:
| Approach | Tradeoff |
|----------|----------|
| Use `go.Densitymapbox` with larger radius | Better gap filling, slower rendering |
| Increase marker size to ~4-5 | Fills gaps, may look blocky |
| Use `go.Heatmap` with image trace | True raster display, loses interactivity |
| Overlay as GeoTIFF on Leaflet | Pixel-perfect, requires different stack |

*Colorscale*:
- Trend maps should use symmetric scaling around 0: `max_abs = np.percentile(np.abs(all_values), 98); cmin, cmax = -max_abs, max_abs`
- Convention: RdBu with 0=white, positive (good)=blue, negative (bad)=red
- For "lower is better" variables, may need to invert or use RdBu_r

**Status**:
- Colorscale issue: **FIXED** (2026-01-22) - Added conditional symmetric scaling for trend metric in `generate_maps.py` lines 602-610
- Alignment issue: Documented for future improvement (low priority, reports are usable)

---

### 2026-03-31: Comparison Report Maps Blank for Small-Valued Variables

**What happened**: After reprocessing qr in raw kg/m²/s units (values ~10⁻⁶), the `compare_water_index.py` HTML report showed all maps as blank/zero.

**Root cause**: `_subsample()` used `np.round(arr, 2)` — fixed 2 decimal places. Values of 0.000005 kg/m²/s round to 0.00.

**Impact**: All heatmap panels in the comparison report appeared empty. Statistical tables and correlations were unaffected (computed before rounding).

**Fix applied**: Changed `_subsample()` to use adaptive significant-figure rounding (4 sig figs based on `max_abs` of the array) instead of fixed decimal places. `compare_water_index.py` line 170.

**Lesson**: Avoid fixed decimal rounding in generic visualization code — use significant figures when data scale varies across variables.

---

## Adding New Incidents

When documenting a new incident, include:

1. **Date and title**: Brief description of what went wrong
2. **What happened**: Detailed narrative of the error
3. **Impact**: What was the consequence (data loss, incorrect results, wasted time)
4. **Root cause**: Technical explanation of why it happened (if applicable)
5. **Correct action**: What should have been done instead
6. **Fix applied**: How was it resolved (if applicable)
7. **Rule created**: Reference to the guardrail added/updated
