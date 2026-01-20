# Guardrails - Critical Rules

Rules in this file represent mistakes that were made and must never be repeated. These are non-negotiable requirements that override convenience or efficiency considerations.

See [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) for detailed incident documentation.

## 1. Temporal Resolution Must Be User-Selected

**Rule**: When ISIMIP data is available at multiple temporal resolutions (daily, monthly, annual), **ALWAYS ask the user which resolution they want before downloading**.

**Why this matters**:
- Monthly data is 12x larger than annual data
- Daily data is 365x larger than annual data
- A 10 MB annual file becomes 120 MB monthly or 3.6 GB daily
- Downloading the wrong resolution wastes significant time and storage

**Required workflow**:
1. Search ISIMIP and identify ALL available temporal resolutions
2. Present a summary table showing:
   - Resolution options (daily, monthly, annual)
   - Approximate file sizes for each
   - Number of datasets available at each resolution
3. Use `AskUserQuestion` to have the user explicitly choose
4. Only proceed with download after user selection

**Handling higher-frequency data (monthly/daily)**:
- See Rule §2 below for aggregation requirements
- Current processing scripts are designed for annual data aggregated to decadal metrics
- If monthly/daily data is selected, must ask user how to aggregate before processing

---

## 2. Sub-Annual Data Aggregation Requires User Choice

**Rule**: When processing monthly or daily data, **ALWAYS ask the user** how to aggregate before processing. Never auto-aggregate without explicit user approval.

**Why this matters**:
- Different variables require different aggregation methods based on their physical meaning
- **Density metrics** (g/m², kg/m², individuals/km²) → typically **mean**
- **Count/flux metrics** (events, kg/year, mm/day) → typically **sum**
- **Extreme metrics** (max temperature, min precipitation) → typically **max** or **min**
- **Categorical metrics** (dominant type, mode) → typically **mode**
- Wrong aggregation produces physically meaningless results

**Required workflow**:
1. Identify that data is sub-annual (monthly or daily)
2. **Tell user the variable units** (critical context for choosing aggregation)
3. Present aggregation options: mean, median, sum, min, max, mode
4. Use `AskUserQuestion` to get explicit user choice
5. Only proceed with processing after user selection
6. Document the aggregation method used in output metadata attributes

**Example aggregation guidance**:
| Variable | Units | Recommended Aggregation |
|----------|-------|------------------------|
| Biomass density | g C m⁻² | mean |
| Precipitation | mm/month | sum (for annual total) |
| Temperature | °C | mean (for annual avg), max/min (for extremes) |
| Fire events | count | sum |
| Species richness | count | mean or max |

---

## 3. Scenario Discovery Must Be Dynamic

**Rule**: Visualization and processing scripts must **discover scenarios from the filesystem**, not hardcoded lists. Never assume a fixed set of scenarios exists.

**Why this matters**:
- ISIMIP data includes varying scenarios: SSP (126, 245, 370, 585) and RCP (26, 45, 60, 85)
- Not all scenarios are always available for every variable/model combination
- Hardcoded lists silently exclude data that exists but isn't in the list
- Users may have custom scenario subsets

**Required behavior**:
- Scripts should glob for `{variable}_*_processed.nc` pattern
- Extract scenario names from filenames dynamically
- Load ALL matching files, not just predefined scenarios
- Use fallback labels/colors for unknown scenarios (e.g., "ssp245" → "SSP2-4.5")

