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

## Adding New Incidents

When documenting a new incident, include:

1. **Date and title**: Brief description of what went wrong
2. **What happened**: Detailed narrative of the error
3. **Impact**: What was the consequence (data loss, incorrect results, wasted time)
4. **Root cause**: Technical explanation of why it happened (if applicable)
5. **Correct action**: What should have been done instead
6. **Fix applied**: How was it resolved (if applicable)
7. **Rule created**: Reference to the guardrail added/updated
