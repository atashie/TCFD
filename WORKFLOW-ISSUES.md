# Workflow Issues Log

This document tracks workflow mistakes and their resolutions. Each incident documents what went wrong, why it matters, and how to prevent recurrence.

See [GUARDRAILS.md](GUARDRAILS.md) for the rules derived from these incidents.

---

## Incident Log

**Index** (44 entries, chronological):

- 2026-01-16: Fish TCB Downloaded Without Resolution Choice
- 2026-01-20: Fish b30cm Processed Without Aggregation Choice
- 2026-01-20: Scenario rcp45 Excluded from Visualization
- 2026-01-21: Loblolly Pine Search Missed Climate-Specific PFT Datasets
- 2026-01-22: QA Report Map Alignment and Colorscale Issues
- 2026-03-31: Comparison Report Maps Blank for Small-Valued Variables
- 2026-04-07: Quantile Breakpoints 12x Too Large for Flux Variables
- 2026-07-24: led Drought Exposure — Binary Variable Needed Dedicated Processing
- 2026-07-24: Missed `let` Tropical-Cyclone Exposure by Guessing a Specifier Code
- 2026-07-24: let Tropical-Cyclone Exposure — Fractional Data, Thin Ensemble, Two-Tier Percentile
- 2026-07-24: burntarea Wildfire — Per-Member Metadata Divergence, No Normalization, Anchored Trend
- 2026-07-25: csoil Soil-Carbon — Catalog Incompleteness, API Anti-Bot, Mixed-CO₂ Ensemble
- 2026-08-07: Wildfire Availability Review Answered Against Other Layers, and Skipped Its Own Skill
- 2026-08-08: Wildfire Inventory — Recommendation Stapled to the Inventory, and 75 Minutes of Avoidable Enumeration
- 2026-08-10: ISIMIP3b burntarea Wildfire Layer — Memory Thrashing Cost 12 Hours, and Two Silent-Zero Bugs
- 2026-08-10: HTML Dashboard Rework — Tab Structure, Zero-Centred Scales, and a Browser-Payload Ceiling
- 2026-08-11: Tropical-Cyclone Inventory Missed the Newest Dataset by Path-Guessing; `let` Rebuilt With a Third Statistic Branch
- 2026-08-11: `led` Drought Rebuilt — Mixed soc, a Mask Rule Worth 10% of Land, and a Hazard That Does Not Mean What Its Name Says
- 2026-08-11: An Unverified Negative in Our Own Catalog Outlived Its Own Correction by Four Days
- 2026-08-11: A Complete Enumeration Filtered by an Unverified Token — Sugarcane Nearly Reported Absent
- 2026-08-11: `driedarea` — the ISIMIP3b/SSP Drought Layer the Catalog Said Did Not Exist
- 2026-08-11: Sugarcane Layers WITHDRAWN — The Model Does Not Grow Cane in the Cane Belt, and the Contract Test Passed Anyway
- 2026-08-11: Sugarcane Yield Layers Built (superseded by the withdrawal above) — Structural Zeros, and a From-the-END Offset That Merged Two Scenarios
- 2026-08-12: Temperate Conifer NPP Layer — A Measured Denominator, a Per-Scenario Variable-Name Drift, and a Time-Varying Presence Mask
- 2026-08-12: Customer Delivery Pipeline Built — External Review Found Eight Defects, One of Them in Shared Extraction Code
- 2026-08-13: Climate Score — A Composition Artifact That Reappeared at a Second Level of Rollup After Being "Fixed" at the First
- 2026-08-13: Delivery Dashboard — Customer Text That Could Execute, a Silent Plotly Failure, and Two Phantom Bugs From Browser Cache
- 2026-08-13: Stages 3 and 4 Built — The Dangerous Failure Is a Correct Report That Reads as a Complete One
- 2026-08-13: First Report Review — A Regional Narrative for the Wrong Region, a Metric That Should Not Have Been Published, and Two Figures Nobody Could Read
- 2026-08-13: External Review of the Report Tooling — Nine Confirmed Overstatements, Including Three Regulatory Claims That Were Simply Wrong
- 2026-08-13: Crop-Failure Ingest — A Recorded Negative That Blocked a Hazard for Weeks, a Product That Zero-Fills the Ocean, and a Contract Variable I Nearly Published as a Constant
- 2026-08-14: The CaMa-Flood Suite We Never Saw — Three Weeks of Treating a Known Path as a Search
- 2026-08-14: Heatwave Ingest — A Sidecar Trap That Downgraded a Shipped Layer's Provenance, a Layer Whose Flat Trend Means the Opposite, and an Index That Changed Between Rounds
- 2026-08-14: Permafrost Layer — Four Estimator Traps, Three of Them Mine, and a Model I Nearly Dropped for the Wrong Reason
- 2026-08-15: `csoil` rebuilt — a lost artifact, an ensemble understated by a third, and a parallel port that OOM-killed itself twice
- 2026-08-16: The threshold ladder — a censoring result that was a continent, an invented caveat name, and a `git checkout` I should not have run
- 2026-08-18: Nature 2026 Hail Deposit -- Evidence Gate Before Any Layer Decision
- 2026-08-18: A Spatial Average Has a WEIGHTING, and Choosing It Wrong Flipped the Sign in Six Regions
- 2026-08-18: The precipitation layers — a mask that published as a real zero, and a slope I chose for consistency instead of measurement
- 2026-08-18: Water-stress build — five defects, three external review rounds
- 2026-08-19: Landslide — the obvious aggregation was degenerate, and the licence is the real blocker
- 2026-08-20: Tornado QA pages were never human-reviewable
- 2026-08-21: Date rollover split a delivery across two dated folders
- 2026-08-21: Dashboard table rebuilt itself out from under its own dropdowns

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

### 2026-04-07: Quantile Breakpoints 12x Too Large for Flux Variables

**What happened**: The annual quantile breakpoints (vt13-19) for flux variables (potevap, qr) were computed by **summing** 12 monthly rate values instead of averaging them. Since the output is stored in raw rate units (kg/m²/s), this produced quantiles ~12× larger than the annual mean (vt12).

**Root cause**: In `process_water_variable.py` lines 436-439, the quantile annual aggregation branched on `self.var.aggregation`:
```python
if self.var.aggregation == "mean":
    annual = subset.groupby("time.year").mean(dim="time")
else:
    annual = subset.groupby("time.year").sum(dim="time")
```
For flux variables with `aggregation="sum"` (potevap, qr), this summed the 12 monthly rates. But vt12 (annual mean) was correctly computed as `nanmean(vt0-11)`, creating an internal units inconsistency.

**Evidence**: At (0°N, 30°E) for potevap ssp370 decade 2050:
- vt12 (Annual Mean) = 4.252e-05 kg/m²/s
- vt16 (Annual Q50) = 5.064e-04 kg/m²/s (before fix)
- Ratio vt16/vt12 = 11.91 ≈ 12 (number of months)
- Old RCP file had consistent Q50 ≈ Annual Mean

**Impact**: Broken quantile breakpoints (vt13-19) in both potevap and qr SSP output files. Monthly means (vt0-11) and annual mean (vt12) were unaffected.

**Fix applied**: Changed quantile computation to always use `.mean()`, keeping all 20 value types in consistent rate units. Both potevap and qr reprocessed and comparison reports regenerated. Confirmed post-fix: potevap vt16/vt12 = 0.992, qr vt16/vt12 = 0.761 (right-skewed, expected).

**Rule updated**: GUARDRAILS.md Section 7 — clarified that quantile annual aggregation always uses mean to match vt12 units.

---

### 2026-07-24: led Drought Exposure — Binary Variable Needed Dedicated Processing

**What happened**: Building the first TCFD drought hazard layer from `led` (Lange 2020 drought exposure, ISIMIP2b) surfaced several ways the standard continuous-variable pipeline would have produced wrong output. Each was caught before it reached the deliverable.

**Findings and resolutions**:

1. **`led` is a BINARY per-cell annual flag {0,1}**, not a fractional field. Verified across all 8 impact models (0.00% intermediate values). Its CF name "land area fraction exposed to drought" refers to the *spatial/ensemble aggregate*, not the grid cell. The standard decadal **median**-over-years would collapse a binary series to 0/1 — meaningless. **Resolution**: decadal statistic = **mean** (drought frequency). The ensemble mean across 31 members × 10 years is the fractional "land area impacted" field.

2. **`process_qg.py` cannot be reused.** It parses the GCM from filename field `[1]`, but lange2020 files are `lange2020_{model}_{gcm}_ewembi_{scenario}_...` — field `[1]` is the impact model, `[2]` the GCM. Reuse would have grouped 8 models as if they were GCMs and mis-concatenated overlapping years. **Resolution**: dedicated `scripts/process_led_drought.py`; ensemble member = `(model, gcm)`.

3. **Time axis is `years since 1661` on a `360_day` calendar** — xarray/cftime cannot decode it (`OutOfBoundsDatetime` / "years not recognized"). **Resolution**: open with `decode_times=False` and parse `year = 1661 + time_value`.

4. **Percentile-of-score compressed to 79–100** (should center ~50 at 2020s). Cause: scoring the ensemble-*mean* map against a per-*member* baseline that is ~79% zeros (zero-inflation). **Resolution**: score against the 2020s **ensemble-mean spatial** distribution (apples-to-apples) → proper 14–100 spread, 2020s median ≈51.

5. **CI ordering failed in 64% of cells** (upper_ci < median). Cause: Q25/Q75 across members don't bracket the *mean* of a right-skewed zero-inflated distribution. **Resolution**: CIs = ensemble **mean ± 1 inter-model SD, clamped [0,1]** → 100% ordered.

6. **Trend/change colorscale read backwards for a hazard, plus a latent caption bug.** `generate_maps.py` hardcoded `RdBu` (Plotly convention: low=red, high=blue), so *increasing* drought rendered blue; the change caption also claimed "Positive values (red)" while actually rendering blue. **Fix applied**: made `create_map_figure` accept `reversescale`; the generator now reverses trend/change to red=worse when the processed file sets `percentile_direction: higher_is_worse`, and computes the caption from the real direction. Variables without the attribute are unchanged. See GUARDRAILS.md §5.

**Also noted**: `led` is ISIMIP2b/RCP only (rcp26, rcp60 — no rcp85, no SSP), a mismatch to keep in mind against an SSP-based TCFD product. `scripts/test_shared_baseline.py` is hardcoded to `qg`/SSP filenames and skips `led`/RCP — shared-baseline identity was instead verified directly (max|diff|@2020s = 0 across all five fields).

**Rule updated**: GUARDRAILS.md §5 — colorscale direction keyed to `percentile_direction`.

---

### 2026-07-24: Missed `let` Tropical-Cyclone Exposure by Guessing a Specifier Code

**What happened**: Asked what TC data ISIMIP offers, I guessed the Lange 2020 exposure code for tropical cyclones was `letc` (le + "tc"). `variable=letc` returned `count=0`, and I concluded the Lange framework has no TC member and that TC data is "input-side only." A free-text fallback (`query="tropical cyclone"`) returned `count=1001` — the API maximum, i.e. truncated — dominated by GeoClaw storm-surge and 6-hourly wind/rain *forcing*, which reinforced the wrong conclusion. The user corrected me: the code is `let`, and TC exposure exists.

**Impact**: Would have missed a ready-to-process TCFD tropical-cyclone hazard layer (`let`, 28 datasets, ISIMIP2b annual, binary — structurally identical to `led`), and mischaracterized ISIMIP's TC holdings, had the user not known the correct code.

**Root cause**: Substituted recall for enumeration on a **controlled-vocabulary** field, then accepted a negative/truncated API result as authoritative. Three guardrails already in the `/isimip-search-download` skill were read past:
- "count = 0 for variables that definitely exist" is a listed red flag → should have triggered file-server/family enumeration, not a "doesn't exist" conclusion.
- "The API returns max 1001 results … results may be truncated" → `count=1001` means the result set was incomplete; generalizing from it was invalid.
- I had already written in the catalog that Lange 2020 has "**6 hazard categories**" and had just processed one (`led`) — the five sibling codes were never resolved.

**Correct action**: Enumerate the family instead of guessing. The Lange 2020 code is always `le` + ONE letter: `led` drought, `leh` heatwave, `lew` wildfire, `ler` river flood, `lec` crop failure, `let` tropical cyclone. All six verified present (28–219 datasets each); `letc`/`lecdd`/`letr` return 0.

**Fix applied**: (1) Authoritative enumerated family table written into `config/isimip_search_catalog.yaml` (`drought.exposure_lange2020.family` + `not_real`) so codes are never re-derived from memory. (2) GUARDRAILS.md §8 added.

**Rule created**: GUARDRAILS.md §8 — Never guess controlled-vocabulary specifier codes; enumerate. Treat `count=0` and `count=1001` as verification triggers, not answers.

---

### 2026-07-24: let Tropical-Cyclone Exposure — Fractional Data, Thin Ensemble, Two-Tier Percentile

**What happened**: Building the TCFD tropical-cyclone hazard layer from `let` (Lange 2020, ISIMIP2b, model `ke-tg-meanfield`) surfaced three design points, two prompted by user challenges.

**Findings and resolutions**:

1. **`let` is a CONTINUOUS FRACTION [0,1), not binary** — unlike `led`. Verified across GCMs: ~97–98% exact 0, the remainder spread smoothly over (0,1) with 100k+ unique values, never exactly 1. Each cell = the fraction of its land area exposed to tropical cyclones that year. I was about to assume it inherited `led`'s binary nature; a user challenge ("I'm pretty sure the raw values deliver a fractional area impacted") prompted the empirical check. **Resolution**: decadal statistic stays MEAN (now = expected annual exposed-area fraction, not a binary-flag frequency); all metadata reframed as fractional. **Rule created**: GUARDRAILS.md §9 — verify data nature empirically.

2. **Thin ensemble → spatial smoothing.** `let` has only 1 impact model × 4 GCMs = 4 members (vs `led`'s 31), so per-cell decadal estimates are noisy. **Resolution**: a 5×5 exponential-decay spatial kernel (w = exp(−d/L), L = 0.7 grid-cells, cos(lat)-scaled longitude distance, normalized over non-NaN land neighbours, longitude-wrapped) applied to each per-member decadal map before ensemble aggregation, so median/CIs/percentile/baseline/trend all derive from smoothed inputs. Smoothing conserves global mass (0.00427 → 0.00434). Ensemble spread is inter-GCM only — flagged in metadata. (`led` is left unsmoothed — its 31-member ensemble doesn't need it.)

3. **Two-tier percentile for a zero-inflated field.** `let`'s 2020s baseline is ~75% exact 0 after smoothing; a single percentile-of-score crushes the real gradient into the top quantiles. **Resolution**: zeros → percentile 1; non-zeros ranked only against the non-zero cells of the 2020s baseline, linearly mapped to [2,100] (2020s non-zero median = 51). Restores full dynamic range for the exposed coasts/tropics.

**Verified**: CI ordering 100% (614k cells); shared 2020s baseline bit-identical across scenarios (max|diff| = 0, all 5 fields); rcp60 global-mean exposed fraction rises ~52% by the 2090s vs rcp26's ~19%.

**Processor**: `scripts/process_let_cyclone.py`. **Rule created**: GUARDRAILS.md §9; the smoothing and two-tier-percentile techniques are documented here.

---

### 2026-07-24: burntarea Wildfire — Per-Member Metadata Divergence, No Normalization, Anchored Trend

**What happened**: Building the TCFD wildfire layer from `burntarea` (ISIMIP2b `biomes`, the direct burnt-area-fraction fire output — distinct from the Lange 2020 `lew` *exposure* family) surfaced three issues: two caught by value-checking, one by a user review of the maps.

**Findings and resolutions**:

1. **Per-member metadata diverges — value-check each (not just the data, the metadata).** All 3 annual models report `burntarea` in **% [0,100]**, but the file metadata disagrees: `lpj-guess` mislabels its `long_name` as "Fire Return Interval" (the values are unambiguously burnt-area %, floored at 0.1% so it never emits a true zero); `mc2-usfs` uses a **`days since 1661`** time axis while `lpj-guess`/`lpjml` use `years since 1661`. A year parser matching only `years since` silently produced garbage years for mc2-usfs, filtered every year out of the 2010–2099 window, and crashed the indexing. **Resolution**: generalized the time parser to handle years/days/months per the file's calendar (365_day → ÷365); treated lpj-guess as burnt-% and flagged the mislabel/floor in output metadata. Reinforces GUARDRAILS §9 — verify data nature **and metadata** empirically.

2. **Same units → no normalization (unlike TWS).** The 3 models agree on the unit (% burnt area) but disagree in magnitude: `mc2-usfs` (a coarse biome model) runs ~5–7× hotter in the mean and is 45% zero-inflated, vs ~0.6% for lpj-guess/lpjml. Because the unit is shared, the magnitude spread is genuine model uncertainty, not a scale artifact. **Resolution** (user decision): **no normalization**, equal-weight "model democracy" in raw %, inter-member spread retained as the CI. Contrast the water-index TWS, which needed robust z-score normalization because its models were on genuinely different scales.

3. **Trend semantics: anchored baseline→decade, not within-decade.** The first cut reused the within-decade annual slope (as `process_qg`/`led`/`let` do). A user review flagged it as spatially spotty/sign-flipping while the change map was coherent — correct, because fire is extremely noisy year-to-year, so a 10-point within-decade slope is mostly interannual noise. **Resolution**: `trend[decade] = (median[decade] − median[2020s]) / elapsed decades` (% decade⁻¹) — the rate *from the 2020s baseline to that decade*, built on decadal means and anchored at the baseline, so it is exactly the (decade−2020s) change map ÷ elapsed decades (corr = 1.0000 with the change map at 2090s). The 2020s baseline has no elapsed change → trend 0 (identical across scenarios). **Rule created**: GUARDRAILS.md §10.

**Verified**: CI ordering 100% (607k cells); shared 2020s baseline bit-identical across scenarios (max|diff| = 0, all 5 fields); global-mean burnt% flat under rcp26, +6.6% rcp60, +14.8% rcp85 by the 2090s.

**Processor**: `scripts/process_burntarea_fire.py`. Ensemble = 3 models × 4 GCMs × {rcp26, rcp60, rcp85} = 12 members/scenario, raw %, no normalization, no spatial smoothing (the 12-member ensemble is thick).

---

### 2026-07-25: csoil Soil-Carbon — Catalog Incompleteness, API Anti-Bot, Mixed-CO₂ Ensemble

**Context**: Building a subsurface/root-zone carbon-storage TCFD layer. Chose the direct soil-carbon pool `csoil-total` (ISIMIP3b biomes) over the vegetation pools (`croot`/`cvegbg`) and the net-sink flux (`nbp`) — all enumerated and presented first.

**What surfaced (and how each was handled)**:

1. **The catalog under-listed the biomes models.** `config/isimip_search_catalog.yaml` documented 5 ISIMIP2b biomes models; the file server has **11** (adds DLEM, JULES-ES-55, LPJ-GUESS, ORCHIDEE, ORCHIDEE-DGVM, VEGAS) and **5** in ISIMIP3b. The per-model `variables:` lists had been compiled for earlier timber/fire searches and are not exhaustive. → Enumerated the file server directly; the catalog is not authoritative for coverage. (GUARDRAILS §8.)

2. **The repository API is behind Anubis anti-bot.** `https://data.isimip.org/api/...` returned an "Access Denied" challenge page to WebFetch. → Fell back to the file server `https://files.isimip.org/...`. WebFetch's page-summarizer also **truncated** long autoindex listings (repeatedly claimed "only one csoil-total file" when there were many). → Switched to `curl -s <dir> | grep -oE` for exact, complete enumeration. (GUARDRAILS §8.)

3. **CO₂-treatment heterogeneity broke the requested config.** The user asked for transient CO₂, but `jules-es-vn6p3` publishes ONLY its fixed-2015-CO₂ (`2015co2`) run for `csoil-total` — no transient `default` — so a strict transient ensemble collapses to 2 structurally-diverse annual models (CLASSIC + MC2-USFS, MC2 being a coarse outlier). Surfaced this rather than silently proceeding; the user chose to **retain all 12 members and accept mixed CO₂** (JULES fixed, CLASSIC/MC2 transient), recorded in `co2_treatment` (JULES's soil-carbon trend is muted — no fertilization). (GUARDRAILS §9.)

**Value-check (clean, unlike burntarea)**: all 3 models report **kg C m⁻²** with comparable magnitudes (2020s medians ~5.8/7.7/10.3) → **no normalization** (model democracy); consistent, correct `long_name`; time axis "days since 1601" for all, only the calendar differs (365_day vs proleptic_gregorian → calendar-aware parse). Direction is **`higher_is_better`** (stored carbon is an asset; risk = loss) → percentile **inverted**. ISIMIP3b csoil starts 2015 → the layer begins at the 2020s baseline (no full 2010s decade).

**Verified**: CI ordering 100% (566k cells); 2020s baseline bit-identical across scenarios; trend↔change corr = 1.0 (ratio = span = 7); global-mean soil carbon +1.1% ssp126 / −1.4% ssp370 / −2.2% ssp585 by 2090s (loss scales with forcing; 37 / 53 / 55 % of land loses carbon).

**Processor**: `scripts/process_csoil_soilcarbon.py`. Ensemble = 3 models × CMIP6 GCMs (2+5+5) × {ssp126, ssp370, ssp585} = 12 members/scenario, raw kg C m⁻², no normalization, no spatial smoothing, `higher_is_better`, baseline-anchored trend. **[Superseded 2026-08-07: slopes now follow OUTPUT-SPEC.md (`ols_slope` + `sen_slope`, expanding window).]** See memory `project_csoil_soilcarbon_layer`, `reference_isimip_enumeration_curl`, `feedback_verify_experiment_tokens`.

---

### 2026-08-07: Wildfire Availability Review Answered Against Other Layers, and Skipped Its Own Skill

**What happened**: Asked to "process wildfire data, review the isimip repository and suggest possible data", I never invoked `/isimip-search-download` — I improvised a `curl` sweep — and then framed the findings comparatively against previously shipped layers, recommending a dataset partly because it was "SSP-aligned with `csoil`" and "structurally identical to `let`/`led`". The user's ask was about wildfire.

**Impact**:

1. **Wrong basis for a recommendation.** `csoil` is a soil-carbon stock; its round and processor shape are not evidence about a fire dataset. The ranking should have rested on cadence, ensemble depth, units, measured data nature, coverage and volume — properties of the candidates.
2. **Re-derived documented findings.** `Zantout2025` and the ISIMIP3a `fire` sector (10 models) are already written into `isimip-search-download/SKILL.md`. I presented both as new discoveries after ~15 tool calls of hand enumeration.
3. **Missed live warnings the skill carries.** `isimip-process-visualize/SKILL.md` already records that `burntarea` **accumulates** (monthly→annual = **SUM**, and copying the `csoil` mean precedent under-scales fire 12×), that `clm45`/`orchidee` declare `%` on a 0–1 fraction scale, and that `classic`/`2015soc-from-histsoc` is **mis-scaled across GCMs within one model**. None reached the user's answer.
4. **Unsound absence claims.** I harvested with **parallel** curls and matched only `\.nc"`. The skill warns that parallel curls get rate-limited and return empty listings indistinguishable from "no data", and that a `.nc`-only filter silently drops the entire ISIMIP2b round (2b publishes `.nc4`). Every "not present" in that review is therefore unverified.

**Root cause**: two compounding faults.

- **No trigger discipline.** Nothing instructed that a repository-availability question *is* the search skill. CLAUDE.md listed the skills by name only.
- **Always-on context was a per-layer changelog.** CLAUDE.md's TCFD section carried three dense narrative paragraphs on `led`/`let`, `burntarea` and `csoil` — loaded every session. Asked about wildfire, the most salient material in context was a comparative essay about other layers, so the answer became one. The `csoil` reference was not a slip; it was what the context primed.

**Correct action**: invoke `/isimip-search-download` first; enumerate serially matching `\.nc4?$`; report the full model × GCM × scenario matrix as inventory; then recommend separately, justified on the candidates' own properties, naming the pending decisions (monthly→annual aggregation §1–§2, soc/sens harmonization, data-nature value check §9).

**Fix applied**:

- CLAUDE.md — replaced the three per-layer narrative bullets with a shipped-layer table plus pointers to each processor docstring and incident entry; added a **Scope discipline** rule (answer within the hazard asked about; never rank a dataset by resemblance to an unrelated layer); made the skill-invocation trigger explicit and declared the skill authoritative over CLAUDE.md for coverage.
- Corrected two stale CLAUDE.md facts that caused the false "new discovery" framing: the Lange 2020 family is **twelve** (`le*` + `pe*` twins), not six, and its ISIMIP3b/SSP re-issue **does** exist (`Heinicke2026`, `Zantout2025`).
- `isimip-search-download/SKILL.md` — new section **"Answering 'what data exists for {hazard}?'"**: stay inside the hazard, separate inventory from recommendation, treat precedent as a hypothesis to re-verify rather than a reason to prefer, and never pre-answer framing choices that depend on measurements not yet taken.

**Rule created**: CLAUDE.md *Scope discipline*; `isimip-search-download` §"Answering 'what data exists for {hazard}?'". Reinforces GUARDRAILS §8 (enumerate, don't conclude from empty results) and §9 (measure, never inherit).

---

### 2026-08-08: Wildfire Inventory — Recommendation Stapled to the Inventory, and 75 Minutes of Avoidable Enumeration

**What happened**: Re-running the wildfire availability review (the 2026-08-07 entry above fixed the *scope* fault: the skill was invoked, the harvest was serial, `\.nc4?$` was matched, and the inventory stayed inside the hazard). Two new faults surfaced.

**Fault 1 — did not stop at the inventory.** The matrix was correct, but the same turn attached a ranked recommendation ("primary" / "runner-up" / "too thin") and an `AskUserQuestion` asking the user to pick a dataset and decide how it related to the shipped layer. The user's correction: *"we should always pause at the inventory for a table for a discussion, and only thereafter should we work together to decide which dataset(s) we may want to process."* Delivering the inventory with a recommendation attached forecloses the discussion the inventory exists to open; the ranking is also the least reliable part of the answer, since it rests on framing choices whose measurements have not been taken yet.

**Fault 2 — ~75 minutes of wall clock on enumeration** (19:19 → 20:32, plus a ~5 min sizing tail), most of it avoidable. Measured from the harvest logs:

| Phase | Wall time | Verdict |
|---|---|---|
| 3b `fire` + `Zantout2025` | 8 min | necessary |
| 3b `biomes` full sector | **28 min** | waste — its `burntarea` is byte-identical to `fire`, already harvested |
| Lange2020 attempts 1 + 2 (empty, then wrong depth) | **~35 min** | waste — one probe request settles depth |
| Lange2020 attempt 3 (correct) | 14 min | necessary |

Root causes, all quantified:

1. **Harvested whole sectors instead of one variable.** 86,065 filenames stored across three sectors; **87–93% were irrelevant** variables (`qtot`, `snd`, `lai`, `tsl`, `trans-*`) for a `burntarea` question. Grep the target variable during the harvest.
2. **Listed directories the inventory never used.** 93 of 142 listings (**65%**) were `historical/` + `pre-industrial/`. List `future/` first.
3. **Harvested a second sector without testing for duplication.** ISIMIP3b `fire` and `biomes` publish the *same* `burntarea` files — 2001/2031 basenames shared, spot-checked pair identical on `Content-Length` **and** `ETag`. A 2-request check would have replaced a 28-minute walk. (Also a real correctness trap: ingesting both double-weights those models.)
4. **Guessed directory depth twice.** `DerivedOutputData/Lange2020` is `{MODEL}/{gcm}/{future,historical,pre-industrial}/`. A 3-level walk returned empty (read as rate-limiting) and a 2-level walk found 0 files at `http 200`. **An `http 200` with 0 matched files means wrong depth, not "no data."**
5. **Sized 36 files sequentially** that were all *exactly* 425 MB. Sample 2–3 per (model, variable, cadence, span) group.
6. **Two off-by-one filename parses.** `DerivedOutputData/` names carry a leading publication token (`zantout2025_`, `lange2020_`), shifting every forward index by one so `$4` reads the forcing and is reported as the scenario. Both produced plausible-looking matrices that had to be re-derived. Parse from the END.
7. **A grep that can never match** — `_lew_.*rcp26`; scenario precedes variable in the filename. A zero-match grep looks exactly like absent data.
8. **Dead-end header chase.** Tried to read CF attributes via 1 MB and 4 MB HTTP range reads; HDF5 puts attribute headers at unpredictable offsets and neither read exposed `units`/`long_name`. Then discovered `python` is not on PATH (only `python3`) and `.venv` lacks `fsspec`, `h5netcdf` and `yaml` — so no lazy remote open and no YAML validation for the catalog edit.

**Fix applied**:

- `isimip-search-download/SKILL.md` — "Separate inventory from recommendation" replaced by **"STOP AT THE INVENTORY"** as a hard stop with an explicit do-not list (no ranking, no `AskUserQuestion`, no pre-answering, no downloading); new section **"Harvest only what the question needs"** carrying all six efficiency rules with their measured costs; new subsections on the leading publication token and on grep token order; the lazy-open bullet now records that `fsspec`/`h5netcdf` are absent and that range-reading CF attributes does not work.
- CLAUDE.md — new first bullet under the TCFD section: availability questions stop at the inventory.
- `config/isimip_search_catalog.yaml` — wildfire section re-enumerated (see below). Edited **unvalidated**: no YAML parser is available in any interpreter on this machine, so structure was verified by reading it back rather than by parsing.

**Inventory recorded** (2026-08-08, for the record — five representations, not three): 3b `burntarea-total` 12–17 members/ssp monthly (4.4–5.8 GB); 3b `mc2-usfs` **the only annual `burntarea` publisher in the whole round**, 5 members (180 MB); 3b `Zantout2025` `wildfire` exposure 12 members annual (14 GB, units unverified); 2b `Lange2020` `lew`/`pew` 19 members annual but **rcp26/rcp60 only** (456 MB); 3b `ffire` emissions; and the 3a `fire` sector's 10-model pool which has **zero ssp/rcp tokens** and cannot yield projections. `jules-es-vn6p3` publishes no `burntarea` despite listing in the `fire` sector.

**Rule created**: CLAUDE.md *availability questions stop at the inventory*; `isimip-search-download` §"STOP AT THE INVENTORY" and §"Harvest only what the question needs".

---

### 2026-08-10: ISIMIP3b burntarea Wildfire Layer — Memory Thrashing Cost 12 Hours, and Two Silent-Zero Bugs

**Context**: Building the SSP-generation wildfire layer from ALL ISIMIP3b `burntarea-total` (user decision: maximize members across both GCMs and impact models, prefer ssp over rcp). Ensemble = 5 impact models × their CMIP6 GCMs = **22 members/scenario** × {ssp126, ssp370, ssp585}, 66 files, 6.03 GiB. Monthly members summed to annual totals; `mc2-usfs` (the only annual publisher in the round) pooled in directly.

**1. Parallelising three scenarios drove a 16 GB machine into swap and cost ~12 hours.**

The loading pass — structurally copied from `process_csoil_soilcarbon.py` — materialised every member as a full 360×720 grid before repacking to land cells: 83 MB × 66 members = **5.5 GB per process**. Running one process per scenario tripled that. Measured at the point of intervention: swap 14.4 GB of 15.4 GB used, 19.5M swapouts, 1.5 billion decompressions, load average 13.7, processes at 55% CPU (I/O-bound, not computing).

The diagnostic tell was **identical work diverging**: the 2070s panel took 2,005 s / 9,217 s / 33,044 s in the three scenarios. Equal workloads differing 16× is a resource problem, never an algorithmic one. I had extrapolated from an uncontended 300-cell benchmark and did not re-check memory when the first panel times came in high.

**Resolution**: two-pass loading — Pass A scans for the land mask holding **one member at a time**, Pass B packs directly into `(member, year, land_cell)`; the full-grid stack is never built. The baseline also now slices its 10-year window *before* concatenating (180 MB, not 1.44 GB). Peak fell 5.5 GB → **2.9 GB**, and the run was done serially in one process. Verified as a pure memory change: the refactor reproduced the pre-refactor baseline exactly (0.7336%, 28.33% zeros on the 300-cell probe).

**Rule**: a processor's memory scales with `n_members`, and `csoil`'s 12-member shape does not transfer to 22. Size the resident set *before* choosing concurrency; when equal-cost panels diverge, check `vm.swapusage` first.

**2. `np.nansum` turned the ocean into land.** The monthly→annual sum used `np.nansum`, which maps an all-NaN ocean cell to **0.0**. The value check then reported **259,200** "land" cells — the entire grid — and every land mean was diluted by ~70% ocean zeros, which made `mc2-usfs` look 4–7× hotter than its siblings. Fixed by requiring all 12 months finite before emitting a year. Same failure class as the `np.zeros()` baseline-panel defect already in GUARDRAILS: **a silent finite zero where NaN belongs.**

**3. `--scenarios` silently changed the "shared" baseline.** Subsetting scenarios pooled the 2020s baseline over only the *selected* ones, so per-scenario runs would each have emitted a different shared baseline (measured drift: 0.7293% vs 0.7336%). Fixed to always pool all scenarios found on disk and only *write* the requested subset.

**Value-check (GUARDRAILS §9), unusually clean for this variable**: all 5 models declare `units="%"`, `long_name="Burnt Area Fraction"`, `days since` on `365_day`, `_FillValue=1e20`. Annual-total land means 1.89–3.72% — same unit, comparable magnitude → **no normalization**. Within-model cross-GCM spread 1.0–1.2×, i.e. no mis-scaling: `classic` was deliberately taken from `2015soc` to avoid its `2015soc-from-histsoc` run, which is documented mis-scaled across GCMs within the one model. Annual totals legitimately exceed 100% where a cell reburns (`classic` ~151%, `elm-eca` ~575%) and are **not** clipped.

**soc is MIXED by design**: no single token spans the ensemble (`elm-eca` only `2015soc-from-histsoc`, `visit` only `2015soc`, `mc2-usfs` only `nat`), so a uniform filter would drop 5–15 of 22 members. Recorded in `soc_by_model` / `soc_treatment`. CO₂ treatment *is* uniform (`default` transient for every member) — unlike `csoil`.

**Verified**: `test_shared_baseline.py` → **ALL CHECKS PASSED**. 2020s panel bit-identical across all three scenarios (max|diff| = 0.0 on median/lower_ci/upper_ci/percentile); 0 CI-ordering violations; 0 ocean leak in every decade; baseline slopes NaN not 0; `n_members` 1–22, `n_models` 1–5 per cell. Global-mean burnt area 0.766% (2020s) → 0.825 / 1.056 / 1.181 % by the 2090s under ssp126 / ssp370 / ssp585; 2090s `ols_slope` +0.031 / +0.153 / +0.212 % decade⁻¹ — a 6.8× trend separation across the forcing range.

**The two-slope rule earned its keep**: `sen_slope` is **exactly 0 on 66–76% of finite cells** (75.8% ssp126, 69.7% ssp370, 66.2% ssp585) and its sign agrees with `ols_slope` only 31–40% of the time. This is the zero-inflation regime OUTPUT-SPEC predicts; the layer records `ols_slope` as the one to read.

**Cost**: exact Theil-Sen at 22 members is ~1.53M pairs/cell at the 2090s panel, 3.4× the `csoil` ensemble. Serial run ≈ 4.5 h for all three scenarios on 10 cores' worth of single-threaded work. `--max-pairs` exists for iteration but was not used in production.

**Processor**: `scripts/process_burntarea_isimip3b.py`; downloader `scripts/download_burntarea_isimip3b.py` (resumable, size-verified against `Content-Length`, writes `download_provenance.csv`). Output `data/processed/wildfire-isimip3b_burntarea-total_annual/burntarea_{ssp126,ssp370,ssp585}_processed.nc`.

**Supersession (user decision 2026-08-10)**: this layer **replaces** the ISIMIP2b/RCP wildfire layer (`process_burntarea_fire.py`) as the shipped wildfire hazard — newer experiment generation, deeper ensemble in both dimensions (22 vs 12 members, 5 vs 3 impact models), and a high-forcing scenario. The 2b processor is kept for provenance with a SUPERSEDED banner: it documents that generation's per-member metadata defects (`lpj-guess`'s "Fire Return Interval" mislabel, `mc2-usfs`'s `days since` axis), which are real findings about *that* data and do **not** recur in 3b. Do not extend it, and do not read its framing decisions as precedent (GUARDRAILS §9). Marked in CLAUDE.md, `scripts/README.md`, and the catalog.

---

### 2026-08-10: HTML Dashboard Rework — Tab Structure, Zero-Centred Scales, and a Browser-Payload Ceiling

**Context**: User review of the first wildfire dashboard. Seven changes, all now standing conventions in `isimip-process-visualize/SKILL.md` — they apply to **every** layer's dashboard, not just this one.

**What changed**

1. **Colorbar titles removed.** A title like "Annual burnt area (percent of grid cell burned per year) [%]" reserved more horizontal space than the map itself. The text became the figure title directly above the map, with the decade as a smaller second line.
2. **Per-metric hover formats** (`HOVER_FORMATS`). Percentiles are integers on [1,100]; `.3e` there was noise, not precision.
3. **Tabs collapsed** from ten to six: `Median | Percentile | Trend | Confidence | Anomaly | Members`. `Lower_Ci`/`Upper_Ci` duplicated Confidence; `Ols_Slope`/`Sen_Slope`/`Change` all belong on Trend, where the two estimators can be read against each other — their *disagreement* is the diagnostic.
4. **Trend tab uses the 2030s and 2090s, not the 2020s.** The baseline's slopes are NaN by contract and its change-from-itself is identically 0, so a 2020s panel is guaranteed blank.
5. **New Members tab** — all 22 LSM × GCM members on one shared colour scale plus a per-member stats table, fed by an optional `{variable}_members.nc` diagnostic `(member, lat, lon)`. Skipped gracefully when absent so older layers still render. `--members-only` rebuilds it from cache in seconds.

**Bug found while doing it**: only the legacy `trend` metric received `reversescale`, so on a `higher_is_worse` layer the `ols_slope`/`sen_slope`/`change` maps rendered **increases in blue**. Now applied to all diverging metrics.

**Zero-centred scales, and why true max failed.** Panels are symmetric about zero (`cmin = -L`, `cmax = +L`). Implemented first with `L = max|value|` as specified, which proved unusable: `max|ols_slope|` is 38.2 %/dec against a **median of 0.049**, so 99.67% of cells sat inside the middle 10% of the colour range and the maps read as blank. Surfaced with numbers rather than silently overridden; the user chose the 95th percentile. Result:

| | true max | 95th pct |
|---|---|---|
| `ols_slope` limit | 38.2 | 1.16 |
| cells using >10% of colour range | 0.33% | **37.5%** |
| clamped at ends | 0% | 5.0% |

Cells beyond the limit are clamped to the endpoint colour by Plotly, never blanked, and each panel now states its scale and clamped share. **Caveat to carry into any writeup**: `sen_slope`'s limit is 14× smaller than `ols_slope`'s (0.080 vs 1.16), so equal redness across those two rows does **not** mean an equal rate.

**Browser payload — the constraint is marker count, not file size.** Every map is an SVG `Scattergeo` with one DOM marker per land cell (~70,849 at 0.5°). The Members page was asking the browser for **1.56 million** markers at 37.6 MB — genuinely past safe. Reduced to 6.2 MB / ~390k, and the layer total from **127 MB to 69 MB**, via `COORD_DECIMALS=2` (exact on a 0.5° grid), `VALUE_SIGFIGS=4`, and `MEMBERS_GRID_STRIDE=2` (Members tab only).

**Subtlety**: downsampling uses a NaN-aware `block_mean()`, not `values[::2, ::2]`. Slicing *samples* every other cell and discards the rest, which on a sparse, zero-inflated hazard silently deletes burning cells — the same "silent data loss that looks like a rendering choice" family as the `np.nansum` ocean bug.

**Known remaining ceiling**: Trend is now the heaviest page at ~425k markers. The next real step, if pages still feel slow, is a raster `go.Heatmap` (equirectangular is exactly linear in lon/lat, so it maps 1:1) at the cost of the coastline overlay. Not done.

**Also confirmed (no change made)**: the wildfire percentile does rank against **non-zero baseline values only**, verified by reconstructing the stored values from the 50,144 non-zero baseline cells to within 4e-6. Zeros map to **1** and non-zeros to **[2,100]**, per OUTPUT-SPEC.md.

> **RESOLVED 2026-08-10.** The user initially described the encoding as 0 / [1,100], then confirmed after review that **1 / [2,100] is correct — keep as is**. No code change. Recorded because the question will recur: the tiers exist so a never-burning cell is distinguishable from the lowest-ranked burning cell, and [1,100] keeps every layer's percentile on one scale regardless of whether it is zero-inflated. Do not "fix" this to 0-based.

**Files**: `scripts/generate_maps.py` (tab structure, scaling, payload), `scripts/process_burntarea_isimip3b.py` (`--members-only`, members diagnostic), `.claude/skills/isimip-process-visualize/SKILL.md` (conventions).

---

### 2026-08-11: Tropical-Cyclone Inventory Missed the Newest Dataset by Path-Guessing; `let` Rebuilt With a Third Statistic Branch

**What happened (the miss)**: Asked to review ISIMIP for tropical-cyclone data, I delivered an inventory that concluded "there is no SSP tropical-cyclone product." That was **wrong**. ISIMIP3b ships `InputData/climate/tropical_cyclones/MIT/` — per-storm wind footprints (Frieler et al. 2025), the newest and most direct TC hazard in the repository. The user had to tell me it existed.

**Root cause**: I listed `ISIMIP3b/InputData/` (which returned `climate/`, `composition_atmosphere/`, `geo_conditions/`, `socioeconomic/`) and then **jumped straight to `climate/atmosphere/`** without ever listing `climate/`. `tropical_cyclones/` sits at exactly that skipped level. Every other branch I walked was enumerated properly; this one I path-guessed because "atmosphere" was the obvious next hop.

**Why the absence claim was doubly unsafe**: I *did* correctly establish that the Lange 2020 exposure family has no 3b re-issue, then let that true finding stand in for the broader claim that no SSP TC data exists. A verified negative about one product family is not a negative about a hazard.

**Rules created**: GUARDRAILS.md §8 — *list every intermediate directory level, never path-guess past one*; *scope every absence claim to exactly the family enumerated*; *open a directory before inferring from its name*. Also recorded in `config/isimip_search_catalog.yaml` → `tropical_cyclone.enumerated` and in the `/isimip-search-download` skill. GUARDRAILS §8 was additionally corrected while here: it described the Lange 2020 family as "exactly six", contradicting CLAUDE.md and the skill — it is **twelve**, six `le*` land-exposed codes each paired with a `pe*` population-exposed twin.

**Also found on the corrected pass**: `ISIMIP3b/DerivedOutputData/TipESM2025/MIT/` is a **name collision** — it holds water models (CWATM, H08, JULES-W2, …), not tropical cyclones. I had flagged it as promising on the first pass purely because MIT is Emanuel's institution.

**Two decisions on the rebuilt `let` layer** (both user calls, both measured first):

1. **The decadal statistic is the pooled MEAN — a declared deviation from OUTPUT-SPEC.** `let` is 97.84% exact-zero at annual resolution, far past anything shipped (`burntarea` is 29.2%). The spec's median/IQR branch left **2,684** exposed cells against **15,122** under the mean — 93% of exposed land erased — with `lower_ci` 98.7% zero and the two-tier percentile assigning tier 1 to 96% of land. Smoothing does not rescue it (2,684 → still 96% zero). Resolution: a **third branch**, `pooled_mean_zero_inflated`, added to `decadal_stats.py` via a new `central=` parameter and documented in OUTPUT-SPEC.md. The justification is not convenience — the boolean/continuous split is a *proxy* for "is the decade pool degenerate at zero", and `let` is the continuous analogue of `led`, which the spec already grants the mean.

2. **The smoothing kernel was the problem, not the window.** Raw `let` renders as one-cell-wide storm tracks, not an impact field. The existing 5×5 kernel is `L=0.7` — labelled "sharp" in its own code — and keeps **32.1%** of the weight on the centre cell, removing only 63% of the track structure. Measured roughness (mean |cell − 4-neighbour mean| over exposed land, normalized):

   | kernel | centre wt | roughness | exposed cells |
   |---|---|---|---|
   | none (raw) | 100.0% | 0.389 | 10,794 |
   | 5×5 L=0.7 (previous) | 32.1% | 0.142 | 15,122 |
   | 5×5 **L=2.5** (chosen) | 8.1% | **0.044** | 15,122 |
   | 7×7 L=2.0 | 6.8% | 0.035 | 16,664 |

   `L=2.5` was chosen over widening because it costs **no extra spatial reach** — it adds no newly exposed cell, whereas every radius increase asserts that land further from any track is exposed. Radius is physically anchored: 2 cells = 111 km ≈ the hurricane-force wind radius. Applied per member per **year** before pooling.

**Two bugs found and fixed in the verification tooling** (neither in the layer):

- **`test_shared_baseline.py` hardcoded decade index 0 as the baseline.** True for 3b-sourced layers that start at 2020, false for `let`/`led`, which carry a full 2010s panel first. It read the 2010s panel — which legitimately differs across scenarios — and reported a spurious "2020s panel bit-identical" FAILURE. Fixed to locate the baseline from the declared `baseline_decade` attribute. The wildfire layer re-verifies unchanged.
- **The slope-agreement INFO line was diluted by permanently-zero cells.** A cell that never sees a cyclone has a genuinely zero slope under *both* estimators, so counting it inflates agreement. All-cell view: 73% sign agreement, 99.2% Sen-zero. Active-cell view (either slope non-zero): **3.0% and 97.0%** — opposite conclusions from the same array. Same dilution family as the earlier ocean-diluted `sen_slope == 0` misreport. The test now prints the active-cell figure with the all-cell figure in parentheses.

**Delivered**: 4 members/scenario × {rcp26, rcp60}, 68,249 land cells. Global-mean exposed fraction 0.004274 → **0.005116** (rcp26, +19.7%) and 0.004001 → **0.006068** (rcp60, +51.7%) by the 2090s. `test_shared_baseline.py`: ALL CHECKS PASSED. Dashboard 38 MB (wildfire was 69 MB). **Read `ols_slope`** — on baseline-exposed land Sen is exactly 0 on 96.9% of cells vs OLS's 1.8%.

**Known limitation, recorded not hidden**: a single impact model means `n_models` is 1 everywhere and the CI carries **inter-GCM spread only** — this layer represents no structural impact-model uncertainty.

**Rules created for the statistics half**: GUARDRAILS.md §9 — the third decadal-statistic branch; *the smoothing decay length `L` is a per-layer measurement, not a constant*; *prefer re-weighting inside the existing footprint over widening it*; *judge slope agreement on ACTIVE cells only*. OUTPUT-SPEC.md gains "The third branch: extreme zero-inflation".

**Files**: `scripts/process_let_cyclone.py` (rewritten), `scripts/download_let_cyclone.py` (new), `scripts/utils/decadal_stats.py` (`central=`), `scripts/test_shared_baseline.py` (two fixes), `isimip-pipeline/tests/test_decadal_stats.py` (+4 tests for the branch), `OUTPUT-SPEC.md`, `GUARDRAILS.md`, `CLAUDE.md`, `config/isimip_search_catalog.yaml`, `scripts/README.md`, both skills.

---

### 2026-08-11: `led` Drought Rebuilt — Mixed soc, a Mask Rule Worth 10% of Land, and a Hazard That Does Not Mean What Its Name Says

**What happened**: The shipped drought layer was rebuilt onto the OUTPUT-SPEC contract. The 2026-07-24 build was family-B — it emitted the retired single `trend` (an OLS slope fitted *inside* each decade, i.e. mostly interannual noise), built the CI from the spread of per-member decadal means rather than a pooled sample, and shipped no `ols_slope` / `sen_slope` / `n_members` / `n_models` and no `members.nc`. Both raw and processed data were also gone from disk, so the layer was re-ingested from source.

**The S3 trap**: `data/` was empty and the obvious explanation was the S3 store — commit `1897aae` "Redirect data storage to S3" migrated exactly these processors. It is **not on `main`**. That commit and the csoil expansion live only on `backup/origin-main-pre-force`, orphaned by a force-push, so `storage.py`, `STORAGE.md` and `layer_publish.py` do not exist on the working branch and nothing in the tree references S3. Checking `git log --all` and assuming reachability would have sent the re-ingest looking in a bucket that this branch cannot address. **Check `git merge-base --is-ancestor` before treating a commit's changes as present.**

**Three things enumeration caught that assumption would have gotten wrong:**

1. **soc is MIXED across the ensemble, and per-model.** `led` is `2005soc` for clm45, h08, lpjml, mpi-hm, pcr-globwb, watergap2 and `nosoc` for jules-w1, orchidee. A single hardcoded token — the natural copy from `let`, which is uniformly `nosoc` — would have 404'd on two of eight models and silently shipped a 6-model ensemble. Each model publishes exactly ONE variant, so harmonizing is not switching a token, it is dropping 2 of 8 impact models. All 31 members kept, declared in `soc_treatment`.
2. **The gap is one pair, not two.** The catalog said "mpi-hm missing 2 GCM combos"; the real gap is `mpi-hm × hadgem2-es` alone (a genuine directory 404), absent from **both** scenarios — so composition is uniform at 31/31 and the shared 2020s baseline is valid.
3. **`led` is binary — re-verified, not inherited.** All 62 members: exactly 2 unique values, exact-0 + exact-1 = 100.00%. The processor now *asserts* this at runtime and exits rather than proceeding if the field ever reads continuous. Its CF `long_name` ("Land area fraction exposed to drought") describes the spatial/ensemble aggregate, and taking it at face value is the classic §9 error.

**The mask decision (user call, 10% of land):** the 8 models disagree on the land mask — clm45 60,695 cells, orchidee 75,438, union 75,673, all-8 intersection 60,501. On the 7,557 cells covered by exactly ONE model the 2020s frequency read **1.63×** the all-8 level (0.0594 vs 0.0365). That is an ensemble-**composition** artifact, not climate: inter-model spread on this hazard is 7.8× (h08 0.0085 → orchidee 0.0660 global-mean frequency), and orchidee reads 0.0607 on its solo cells against 0.0666 on shared ones — its own field is not elevated there, the dilution is simply missing. Resolution: publish only where **≥ 2 models** have data → 68,116 cells. Residual after the cut is 696 two-model cells at 0.0594, with **no systematic gradient** below that (5/6/7-model tiers read 0.0088/0.0396/0.0307, *lower* than the 8-model level), so it is 1% of the layer rather than a trend. Recorded in `mask_rule` and `residual_level_step`; `n_models` is emitted per cell for a stricter downstream cut.

**The finding that matters most to a consumer**: **this variable is departure from preindustrial, not aridity.** The flag fires when soil moisture drops below the 2.5th percentile of the cell's *own* preindustrial distribution, so a permanently arid cell with a stable regime scores LOW, and the strongest signal appears where the climate has moved furthest. Cells above 70°N read **0.0754** mean frequency against a 0.0362 global mean (**2.08×**) — and that is real, not thin coverage: mean `n_models` there is **7.42**. Reading a high Siberian or northern-Canadian value as "arid" inverts the variable's meaning. Now carried in the layer as `interpretation_caveat`.

**`sen_slope` is exactly 0 on 100.0% of cells** — the cleanest instance yet of the documented zero-inflation collapse, and a *structural* one rather than a matter of degree: on a binary field every pairwise difference is −1, 0 or +1, and same-valued pairs dominate wherever drought is uncommon, so the median pairwise slope is 0 by construction. Both estimators are still emitted per contract; **read `ols_slope` on this layer**.

**A pickling bug in the parallel path**: at 31 members the final panel is ~3.9M Theil-Sen pairs per cell (~2.2 h single-core), so the slope stage was fanned across forked workers. Returning `expanding_slopes`' result from a worker fails the whole pool with `KeyError('__getstate__')` — `SlopeResult` is a `dict` subclass whose `__getattr__ = dict.__getitem__` turns pickle's probe for `__getstate__` into a KeyError instead of an AttributeError. Workers now return plain arrays. **Any processor parallelizing this module must do the same**; the shared module was left untouched.

**Delivered**: 31 members/scenario × {rcp26, rcp60}, 68,116 cells, 35 min wall (4.5 CPU-hours, 8 workers). Global-mean drought frequency 0.0362 (shared 2020s) → **0.0427** rcp26 (+18%) and **0.0646** rcp60 (+78%) by the 2090s; rcp26's `ols_slope` *decays* (+3.3e-3 → +0.9e-3 dec⁻¹) while rcp60's *accelerates* (+3.5e-3 → +4.4e-3). `test_shared_baseline.py`: ALL CHECKS PASSED. Maps reviewed — the change field is the classic subtropical drying belt (Mediterranean/Middle East, Sahel, southern Africa, Australia, NE Brazil/Central America); the 31-member contact sheet shows correct coastlines, no block structure and no hemisphere flips.

**Files**: `scripts/process_led_drought.py` (rewritten), `scripts/download_led_drought.py` (new), `CLAUDE.md`, `config/isimip_search_catalog.yaml` (drought section — enumerated matrix, soc mix, land masks; plus a **corrected** `simulation_round` that had claimed the Lange 2020 family has no ISIMIP3b re-issue, with a pointer to `driedarea` under Heinicke2026).

---

### 2026-08-11: An Unverified Negative in Our Own Catalog Outlived Its Own Correction by Four Days

**What happened**: Asked to double-check whether ISIMIP3b/SSP really has no drought data, enumeration found **two** complete SSP representations. `driedarea` under `DerivedOutputData/Heinicke2026` is a **gap-free 3 impact models × 5 CMIP6 GCMs × 3 SSPs = 45-file** matrix, uniform `2015soc`/`default`, annual 2015–2100, ~168 MB, with sha512 sidecars. Separately, **9 of 10** `water_global` models publish `soilmoist` under all three SSPs. The catalog had asserted the opposite.

**Precise timeline** (from `git log -S`):

| date | commit | event |
|---|---|---|
| 2026-07-24 | `edeb174` | `led` drought layer ships on rcp26/rcp60 |
| 2026-07-24 | `8fac010` | catalog written: `simulation_round: "ISIMIP2b only (NO ISIMIP3b/SSP version of this family)"` and *"NOT found (0 hits)"* — **both false** |
| 2026-08-07 | `ac524f5` | the `/isimip-search-download` skill gains the CORRECT fact (Heinicke2026 `driedarea` exists) |
| 2026-08-10 | `cf335fb` | CLAUDE.md gains the correct fact too |
| 2026-08-11 | `21d8a5a` | the `let` rebuild edits the very string *"NOT found (0 hits)"* in the TC block — and leaves the drought false negative standing one block away |
| 2026-08-11 | today | enumeration proves the catalog wrong; 45 files, no gaps |

So the false negative survived **18 days**, and for the last **4 of them the repository openly contradicted itself** — skill and CLAUDE.md said the re-issue exists, the catalog said it does not — with no mechanism that treats a contradiction as a defect.

**Root cause — negatives are never tested the way positives are.** A positive catalog claim is self-correcting: you go to the path, and the files are either there or they are not. A negative *ends the search before it starts*, so nothing ever exercises it. Compounding it:

1. **The "0 hits" was real but misinterpreted.** Searching the *code* `led` in ISIMIP3b genuinely returns nothing — the token is 2b-only. That is a fact about a controlled-vocabulary string, not about the hazard, because re-issues are renamed by hazard word (`driedarea`, `floodedarea`, `heatwave`, `wildfire`, `cropfailure`). The true narrow finding was recorded as a broad one.
2. **No receipt.** The claim named no directory and no date, so no later session could tell whether it rested on a listing of `DerivedOutputData/` (it did not) or on a code search (it did).
3. **Precedence was treated as sufficient.** The skill outranks the catalog, so once the skill was right, the answer was reachable — and it was, which is why `driedarea` was flagged in today's readiness check before any of this. But "believe the higher authority" left the wrong text in place for the next session that consults the catalog first.

**Also on me, separately from the docs**: `driedarea` was raised as an option but **not enumerated** before the fork was put to the user — described only as "a different product with a different ensemble", against a fully quantified `led`. Enumerating it costs 15 listings and ~2 minutes, and would have put 3 × 5 × 3 = 45, uniform soc, 168 MB on the table at decision time. A fork with numbers on one side only is a recommendation wearing a question mark. **Rule added to the skill**: if you name an alternative dataset, enumerate it *before* offering the choice; offering to enumerate "if you'd rather" pushes the cost of an informed choice onto the user.

**Rules created**: **GUARDRAILS.md §11 — a recorded negative is a claim, not evidence.** Every negative in the catalog now must carry `verified_absent_on: "<date> — listed <URL>"` naming the directory actually enumerated, or be written as `UNVERIFIED`. A negative about a CODE is not a negative about a HAZARD. Two of our own documents disagreeing is a **stop-and-resolve trigger**, not a precedence question — enumerate, then write the result back to *both*. When editing a doc region, re-read the surrounding claims about the same fact (the `21d8a5a` near-miss). Never repeat a recorded negative to the user without saying which enumeration and date it rests on.

**Applied, not just written**: all four remaining negatives in the catalog were audited — the Lange-2020-no-3b-re-issue claim now carries today's receipt (listed `ISIMIP3b/DerivedOutputData/`, exactly four publications, Heinicke2026 walked in full), and the three I did **not** enumerate (3a TC-flooding "no 3b equivalent", BioScen1.5 biodiversity, protocol-only crop codes) are now explicitly marked `UNVERIFIED` rather than left to look authoritative. The catalog header carries the convention.

**Files**: `GUARDRAILS.md` (§11 new), `.claude/skills/isimip-search-download/SKILL.md` (catalog-negatives section + enumerate-before-forking rule), `config/isimip_search_catalog.yaml` (header convention, receipts, `driedarea` and `water_global` SSP matrices), `WORKFLOW-ISSUES.md`.

---

### 2026-08-11: A Complete Enumeration Filtered by an Unverified Token — Sugarcane Nearly Reported Absent

**What happened**: Asked whether ISIMIP has coffee, and failing that sugarcane, the harvest walked the agriculture sector correctly — 63 ISIMIP3b model×GCM directories, then 2b/3a/2a — and applied a `grep -E '(sgc|cof|coffee)'` filter to each listing as it streamed. It matched **zero files in all four rounds**. Sugarcane exists: `yield-sug-{firr,noirr}`, ISIMIP2b LPJmL, 4 GCMs, rcp26/rcp60, 2006–2099. It surfaced only because the same pass *also* projected the variable field `$(NF-4)` into a distinct vocabulary, and `sug` appeared there. Had the harvest followed the skill's then-current instruction ("**grep the target variable during the harvest**") without that projection, the answer would have been "sugarcane does not exist in ISIMIP" — delivered with a receipt of 150 directories and zero empty listings.

**Root cause — the filter token came from our own catalog, and was real but from the wrong product.** `crop_codes.isimip3b` listed `sgc: "Sugarcane"`. That code is genuine ISIMIP vocabulary — it is one of the 20 crops in `InputData/socioeconomic/crop_calendar/` — but **no ISIMIP3b model publishes output for it**, and the rounds that *do* publish sugarcane (2a/2b) call it `sug`. Three failure surfaces stacked:

1. **`InputData` vocabulary was recorded as if it were `OutputData` vocabulary.** The protocol defines 20 crops; 3b models publish 11. Nothing in the entry said which product it had been observed in.
2. **A pre-filter converts a good enumeration into a false negative that looks authoritative.** The walk was genuinely complete — that is what makes this failure worse than §8's original `letc` case, where the *search* was thin. Here the search was exhaustive and the **reduction** was wrong, so every quality signal (dirs listed, retries, non-empty assertions) read green.
3. **The skill's wall-time rule pointed straight at it.** "Grep the target variable *during* the harvest" (added 2026-08-08 after a 75-minute enumeration) is a real optimisation, but as written it required knowing the answer's name before finding it. It collided with GUARDRAILS §8 and nothing said which wins.

**Impact**: none delivered — the vocabulary projection caught it and the inventory that went to the user included sugarcane. But it was luck, not design, and the user asked the right question afterwards.

**Two factual errors found while re-enumerating properly** (both had been written into the catalog earlier the same day):
- *"NO rcp85 for any 2b crop"* — false. **CLM45 publishes rcp85 for `mai` and `soy`.** The observation was true of LPJmL, the one model whose listing had been read, and was generalised to the round. Sugarcane itself genuinely has no rcp85.
- *"cassava has NO ssp370"* — true of `yield-cas`, false of the crop: `biom-cas` and `nleachcum-cas` do carry ssp370. **Scenario coverage is per (crop, metric), not per crop**, and aggregating over metrics hides the gap in exactly the variable a layer needs.

**Rules created** — GUARDRAILS §8 gains two bullets, and the skill's "Harvest only what the question needs" bullet was rewritten:

- **Project the variable FIELD; never let a believed-in token narrow a harvest.** Reduce with `awk -F'_' '{print $(NF-7), $(NF-5), $(NF-4)}' | sort -u` — same single pass, output is ~100 tokens instead of 86,065 filenames, and the vocabulary is auditable. Match the target against that vocabulary **offline**. A pre-filter can confirm presence; it can never establish absence, and its empty result is `UNVERIFIED` (§11), not a negative.
- **Vocabulary is per (round, product) — record which one you observed.** Codes drift across rounds for the same quantity: `sug`/`sgc`, `ben`/`bea`, `whe`/`swh`+`wwh`, `ric`/`ri1`+`ri2`, `csoil`/`csoil-total`. An inherited code is a hypothesis, not a lookup key.
- **`$(NF-7)` is not the scenario in every round.** ISIMIP2a filenames carry a bias-adjustment token and no soc field, so the from-the-end scenario offset differs; 2a/3a are historical-only and were labelled as such rather than reported with a mis-parsed token.

**Applied**: re-enumerated all 150 agriculture directories with no target filter (0 empty listings) and wrote `crop_availability_matrix` into the catalog — 24 crop codes × round × model × GCM × scenario, yield-only with each crop's other metrics listed. The `crop_codes.isimip3b` protocol-only claim, which had been sitting as `UNVERIFIED` since the audit four days earlier, now carries its receipt.

**Files**: `GUARDRAILS.md` (§8, two bullets), `.claude/skills/isimip-search-download/SKILL.md` (harvest-reduction rule rewritten with worked example; per-round code drift; InputData≠OutputData), `config/isimip_search_catalog.yaml` (`crop_availability_matrix`, `sugarcane`, `cassava`, `coffee`, `perennial_crops`, corrected rcp85/ssp370 claims), `WORKFLOW-ISSUES.md`.

---

### 2026-08-11: `driedarea` — the ISIMIP3b/SSP Drought Layer the Catalog Said Did Not Exist

**What happened**: Built the ISIMIP3b/SSP drought layer from `driedarea` (Heinicke2026), the dataset a stale catalog negative had hidden for 18 days (see the entry above). It ships **alongside** `led`, not instead of it: 3 impact models × 5 CMIP6 GCMs × 3 SSPs versus `led`'s 8 models × 4 CMIP5 GCMs × 2 RCPs. Deeper ensemble against newer scenarios — neither supersedes the other.

**What the §9 value-check caught that assumption would not have:**

1. **The mask defect is per-VARIABLE, not per-product.** `driedarea`'s sibling `floodedarea` — same publication, same three models — is non-NaN over 94.7% of the globe *including ocean*. `driedarea` is 20.4–22.2%, i.e. genuine land. Had the known-bad sibling been used to infer the product, this layer would have been rejected or "fixed" for a defect it does not have.
2. **Time is `days since 1601-01-01` on `proleptic_gregorian`**, not the `years since` axis the ISIMIP2b Lange2020 files use. Decoded with cftime through xarray; days ÷ 365 drifts ~0.4 yr over a 400-year offset, enough to push a January-stamped record into the wrong decade bin.
3. **The `long_name` is the generic "Exposed Area Share"** — it names neither the hazard nor the per-cell nature. Binary {0,1} was established from values (exactly 2 unique values, exact-0 + exact-1 = 100% in all 45 members), independently of `led`. The two agree; that was verified, not assumed.
4. **Filename grammar is per-PUBLICATION, not per-product.** Heinicke2026 files carry **no** leading publication token (`h08_gfdl-esm4_w5e5_ssp126_…`) while Lange2020's do (`lange2020_lpjml_…`). Both live under `DerivedOutputData/`. Parsing from the end is what makes one parser work for both.

**The mask rule was re-measured, not inherited — and came out the other way.** `led` publishes only where ≥2 of 8 models have data because its solo cells read **1.63×** the all-model level across a **7.8×** inter-model spread, one high model undiluted. Measured here, that mechanism is absent: the step is **1.03×** (1 model 0.0745 / 2 models 0.0599 / 3 models 0.0722), inter-model spread is **2.69×**, and solo cells are split across all three models (1,130 / 2,968 / 1,298) rather than dominated by an outlier. So the full union is published — 63,455 cells — and masking would have cost 8.5% of coverage to remove an artifact that is not there. **Two sibling layers of the same hazard legitimately resolve the same knob differently; that is the rule working, not an inconsistency.**

**Two findings from looking at the maps** (neither is an algebraic failure; the contract test passed before and after):

- **Thin coverage concentrates in the arid belt.** 1- and 2-model cells sit at median |lat| 34°/33° against 46° for 3-model cells — hydrological models disagree about desert cells. Mean `n_models`: Sahara 1.57, Arabian 1.87, Gobi 2.10, Australian interior 2.28, versus Amazon 3.00, Europe/Med 2.77. Crucially the **largest changes land on the best-covered cells** (ssp585 2090s−2020s: Amazon **+0.436** at n_models 3.00, Europe/Med **+0.323** at 2.77) while thin arid cells barely move (Sahara +0.034, Arabian −0.017), so the headline signal does not rest on thin coverage. A desert site-level query should still read `n_models` first.
- **The high-latitude signal is largely one model.** In the 60–80°N band: jules-w2 **0.2397**, watergap2-2e 0.0453, h08 0.0272 — an **8.8×** spread against a 2.69× global one. Above 60°N jules-w2 is **2.14×** the ensemble. Visible as a bright Arctic band on the Members tab. Do not report a boreal drought trend from this layer without saying it is one model of three.

Both are recorded as `coverage_geography` and `single_model_dominance` global attributes, patched onto the written files header-only with the data verified bit-identical (netCDF4 append mode), so the 26-minute slope run was not repeated.

**Delivered**: 15 members/scenario × {ssp126, ssp370, ssp585}, 63,455 cells, 26 min wall (8 workers). Shared 2020s baseline **0.0707** → **0.0936** ssp126 (+32%), **0.1687** ssp370 (+139%), **0.1920** ssp585 (+172%) by the 2090s, with `ols_slope` decaying under ssp126 (+7.6e-03 → +3.2e-03 dec⁻¹) and holding steady under ssp585 (~+1.75e-02). `test_shared_baseline.py`: ALL CHECKS PASSED. Maps reviewed — 15-member contact sheet structurally clean; the ssp585 change field is the CMIP6 drying pattern with a markedly stronger Amazon signal than the CMIP5-based `led` layer shows. `sen_slope` is exactly 0 on 100.0% of final-panel cells, the same structural binary-field collapse as `led` — **read `ols_slope`**.

**Files**: `scripts/process_driedarea_isimip3b.py` (new), `scripts/download_driedarea_isimip3b.py` (new), `CLAUDE.md` (shipped-layers table + the sibling-layers note), `config/isimip_search_catalog.yaml`.

---

### 2026-08-11: Sugarcane Layers WITHDRAWN — The Model Does Not Grow Cane in the Cane Belt, and the Contract Test Passed Anyway

**Outcome first**: the two layers below were built, passed every contract check, were mapped — and are **invalid**. The user asked one question about the maps ("why are there higher yields in the American Midwest than in Florida?") and the answer turned out to be an upstream data defect that no check in this repository was looking for.

**What the data does**: ISIMIP2b LPJmL's `yield-sug-*` reads exactly 0 across the entire real sugarcane belt — São Paulo, Uttar Pradesh, Guangxi, Thailand, Pakistan Punjab, Veracruz, KwaZulu-Natal, Cauca, Queensland, Florida, Louisiana. Those cells carry a sentinel signature in the companions: `biom-sug` pinned at exactly **0.267 t ha-1** and **`plantday = matyday = 1`** — planted and matured on day 1, i.e. no season was simulated. 12,966 land cells (19.2%) carry it and **0%** of them have a non-zero yield. The mask is near-identical across all 8 members (Jaccard 0.98–0.998), so it is static, not climate-driven crop failure. Maize from the same model, GCM and run *is* simulated in those cells (Florida 5.65, São Paulo 1.64 t ha-1) — the cells are live cropland; only cane is missing.

**The same model gets it right in ISIMIP2a**, which is what makes this a defect rather than a model quirk:

| region | 2a cells / mean t ha-1 | 2b cells / mean t ha-1 |
|---|---|---|
| Florida | 84 / **19.49** | 0 / 0.00 |
| Gulf LA-TX | 200 / 12.77 | 3 / 6.04 |
| São Paulo | 140 / 11.69 | 0 / 0.00 |
| Uttar Pradesh | 157 / 8.33 | 0 / 0.00 |
| Queensland | 193 / 7.17 | 2 / 12.63 |
| US Midwest | 276 / 10.61 | 307 / **13.35** |

75.2% of the cells where the 2a run grows cane are zero in 2b; 90.3% of 2b's sentinel cells are 2a cane land (mean 10.57 there). Both rounds report marginal *potential* yield in temperate cells — LPJmL simulates a crop wherever it can, not only where it is farmed, and that part is a model property, not the bug. The bug is that deleting the real belt promotes those marginal temperate cells to the visible maximum, which is exactly the inversion the user spotted.

**Root cause of MISSING IT — every check we ran was about form, not meaning.** `test_shared_baseline.py` verified schema, shared baseline, CI ordering, percentile range and orientation, slope masks, ensemble depth — and passed both layers cleanly, because all of that was true. GUARDRAILS §9 made me measure the data's *nature* (continuous, t ha-1, zero-inflated) and I did, carefully, and even wrote a paragraph rationalising the 87% zeros as "structural — LPJmL grows no sugarcane there". That sentence was **literally true and completely wrong about what it implied**: I checked *how many* zeros and *whether they moved*, never *where they were*. A single question — "does the non-zero mask contain the places this crop is actually grown?" — would have caught it in one minute.

**Rule created — GUARDRAILS §12, plausibility of PLACE, not just of value.** For any layer describing a thing with a known real-world geography (a crop, a biome, a fishery, an industry), verify the field is non-trivial **in a named list of reference locations where that thing demonstrably exists**, before processing and again before shipping. Record the sites and values in the processor docstring. Where a second round or a second model publishes the same quantity, spot-check the same sites in both: a cross-round contradiction of this size (Florida 19.49 vs 0.00) is visible in two lookups. A layer can satisfy every line of OUTPUT-SPEC.md and still be about nothing.

**Disposition**: both processed layers and both map sets retained on disk as evidence, each carrying `INVALID-DO-NOT-USE.md`. CLAUDE.md's shipped-layers row marked WITHDRAWN. The catalog carries a `DATA_DEFECT` block with the measurements. **No usable future-scenario sugarcane yield exists in ISIMIP**: 2b LPJmL is the only scenario-bearing source and it is defective, ISIMIP3b publishes no sugarcane, and the correct ISIMIP2a run is historical-only.

The build record below is kept because the framing work and the parse bug are still worth having.

---

### 2026-08-11: Sugarcane Yield Layers Built (superseded by the withdrawal above) — Structural Zeros, and a From-the-END Offset That Merged Two Scenarios

**Delivered**: two layers from ISIMIP2b LPJmL, `yield-sug-noirr` (rainfed) and `yield-sug-firr` (fully irrigated), each 4 members (1 impact model × 4 CMIP5 GCMs) × {rcp26, rcp60}, annual 2010–2099, uniform `2005soc`/`co2`, units **t ha-1 yr-1** (read from the file, not inferred). 16 raw files / 61.6 MiB, all sha512-verified. `test_shared_baseline.py`: **ALL CHECKS PASSED** on both layers. LPJmL is the only sugarcane source with future scenarios anywhere in ISIMIP; there is no SSP version, because ISIMIP3b publishes no sugarcane at all.

**Bug caught before delivery — awk offsets used as Python indices.** `parse_name` documented itself as parsing from the END (correctly) and then used `p[-7]` for the scenario, transcribing awk's 1-based `$(NF-7)`. Python's `p[-7]` is one field further right, so `2005soc` was read as the scenario, `co2` as the soc token, and the variable name as the sens token. The run completed cleanly and wrote a single plausible-looking `yield-sug-noirr_2005soc_processed.nc`: because members are keyed by `{model}_{gcm}`, rcp26 and rcp60 collided on the same four keys and **last-write-won**, so the "8 member" file was silently rcp60 only. Nothing failed — the log said `Members: 8 | scenarios: ['2005soc']`, which is the only visible tell. Correct offsets are `[-8]=scenario [-7]=soc [-6]=sens [-5]=variable`; awk's `$(NF-4)` and Python's `p[-5]` are the same field. The processor now **asserts** the parsed scenario matches `(rcp|ssp)\d{2,3}|historical|picontrol` and the variable starts with `yield-`, so the offsets cannot silently shift again.

**Framing decisions, each measured first:**

- **Statistic = ordinary `pooled_median`, not the third branch.** The field is 87.3% exact-zero, which looks like `let` territory, but the zeros are STRUCTURAL: 87.0% of land is zero in *every* one of the 94 years — LPJmL grows no sugarcane there. The OUTPUT-SPEC test is whether the median *erases* a signal, and it does not: median-branch exact-zero land 87.27% vs mean-branch 86.94%, a **0.33 pp** gap (`let`'s is 18 pp). Inside the growing region the median is a real yield — positive for 97.45% (noirr) / 99.78% (firr) of cells that ever grow.
- **Percentile = two-tier + `higher_is_better`**, giving the convention the user asked for: zero-yield → **100 (highest risk)**, growing cells ranked against the non-zero baseline over [1, 99]. Measured consequence, stated up front: **87.3% of land reads percentile 100**, so that field is a suitability map first and a risk gradient second. `median` and both slopes are unaffected.
- **NO spatial smoothing — a declared deviation** from the thin-ensemble default. The structure here is a hard cropland mask; a 5×5 kernel would bleed yield into cells that grow none, converting structural zeros into small positives and moving them off percentile 100. The noise smoothing would suppress is inter-GCM spread, which the IQR already carries.
- **`2005co2` sensitivity excluded**: rcp26 publishes it, rcp60 does not, so including it would make the two scenarios experimentally asymmetric.

**Results** (2020s → 2090s, over growing cells): irrigated rcp60 **−15.0%** with 99.5% of cells declining — with water limitation removed the residual signal is heat, and it is nearly unanimous. Rainfed rcp60 is **−0.2% net** but strongly dipolar (Brazil −5.1%, SE Asia −5.8%, US-Gulf −10.6% against Africa +14.2%, Australia +13.4%), i.e. a rainfall-plus-CO₂ signal rather than a uniform one. Rainfed rcp26 (−2.2%) declines *more* than rcp60 in the global mean, consistent with less CO₂ fertilization in the low-forcing scenario. Slope agreement over active cells: 81.6–89.2% (noirr), 95.7–98.3% (firr).

**Files**: `scripts/download_sug_sugarcane.py` (new), `scripts/process_sug_sugarcane.py` (new), `CLAUDE.md` (shipped-layers row), `config/isimip_search_catalog.yaml`, `WORKFLOW-ISSUES.md`.

---

### 2026-08-12: Temperate Conifer NPP Layer — A Measured Denominator, a Per-Scenario Variable-Name Drift, and a Time-Varying Presence Mask

**Delivered**: `npp-tempnle`, temperate needleleaf evergreen (conifer) stand productivity from ISIMIP2b `biomes` — CLM45 + ORCHIDEE + LPJmL × 4 CMIP5 GCMs × {rcp26, rcp60, rcp85}, `2005soc`/`co2`, annual 2010–2099, emitted in **g C m-2 yr-1**. 53 raw files / 454 MiB, all sha512-verified. `test_shared_baseline.py`: **ALL CHECKS PASSED**. Ported from the pre-spec `process_loblolly_npp.py` onto `decadal_stats`.

**Three defects found and fixed, none of which the contract test would have caught:**

1. **The denominator is not shared across models — measured, not assumed.** `corr(cover, npp)`: ORCHIDEE **+0.279** (per-tile), LPJmL **+0.898** (cover-scaled per grid cell); across cover quintiles LPJmL's NPP rises ~2,300× against ORCHIDEE's 2.7×. CLM45 publishes **no `pft-` fraction for any PFT** and reports only on its own tile. Median over 2,488 common cells: raw **586 / 471 / 162** (3.62× spread), per-tile **586 / 471 / 696** (**1.48×**), per-gridcell impossible for CLM45 (4.49× for the other two). The layer therefore ships on the **per-tile basis** — NPP per unit conifer-stand area — with LPJmL divided by its cover and a **2% minimum-cover presence mask**: unthresholded, LPJmL's per-tile p99 reaches **160,203 g C m-2 yr-1**, at 2% it is 5,200. Also confirmed the documented lying-units trap: ORCHIDEE's `pft-tendev` declares `units='%'` and stores a **fraction** (max 0.618) while LPJmL stores true **percent** (max 95) — the scale is decided from the values, never the attribute.

2. **The variable name changes separator BY SCENARIO within one model.** LPJmL writes `npp-temperate-needleleaved-evergreen-tree` (hyphens) in rcp26/rcp60 and `npp_temperate_needleleaved_evergreen_tree` (**underscores**) in rcp85 — same model, same PFT, same run. Building the variable name from the filename crashed on exactly those 8 of 53 files. `resolve_var()` now tries the exact name, then the underscore form, then the single 3-D data variable. Checked across all 53 files: isolated to LPJmL rcp85.

3. **A time-varying presence mask breaks slope/median mask agreement.** The cover threshold is applied per year, so a cell can carry observations inside the *expanding* slope window (2020–2039) yet none inside the decade window itself — finite slope, NaN median. The contract's mask-agreement assertion fired on the first run (`ols_slope finite where median is NaN at 2030s`). That is not an ocean leak; it is stands appearing and retreating. Slopes are now masked to the decade's own median mask, with the count logged (53 cells at rcp26 2030s rising to 374 by 2090s). The masked cells were the artifacts: removing 53 of 25,821 moved the rcp26 2030s mean slope from **−1.89 to +0.64**, because a stand vanishing mid-window produces a wild per-tile trend.

**Also**: slopes are now computed across forked workers (`--jobs`, mirroring `process_driedarea_isimip3b`). Single-threaded, the 2090s panel stacks 10 members × 80 years = 800 observations → 319,600 exact Theil-Sen pairs per cell over 27,377 cells, and one scenario took ~1 h.

**Ensemble is deliberately ragged.** CLM45 publishes only **7 files** for this PFT across all future scenarios and **no GCM of its own in all three RCPs** (rcp26 hadgem2-es+miroc5, rcp60 ipsl-cm5a-lr, rcp85 gfdl-esm2m+hadgem2-es), so composition is scenario-dependent: 10 / 9 / 10 members. Per user decision the shared 2020s baseline pools **all 29 member-scenario series**, so the panel stays bit-identical across scenarios; `members_by_scenario` records identity. `n_members` 1–10 and `n_models` 1–3 vary by cell.

**Results** (2020s → 2090s over ~25,400 stand cells): rcp26 **−0.2%** (46% of cells increasing, median `sen_slope` −0.25/dec); rcp60 **+32.6%** (91% increasing, +14.2/dec); rcp85 **+40.3%** (90% increasing, +17.6/dec). Reference sites 2020s→2090s rcp85: PNW Oregon 355→449, Georgia 502→663, Bavaria 590→683, Japan 560→863, NZ 1697→1743.

**Interpretive caveat, recorded because it will mislead otherwise**: under `higher_is_better` a rising NPP reads as *falling* risk, so this layer's rcp85 percentile map is greener than rcp26's almost everywhere. That is the CO₂-fertilization plus growing-season response of DGVMs run with transient CO₂ (`co2`), a well-known high-uncertainty regime — nutrient limitation is incompletely represented. The layer measures **stand productivity only**: it contains no drought mortality, fire, pest, or windthrow risk, and productivity is not timber value. The `2005co2` fixed-CO₂ runs exist for all three models and are the natural sensitivity test if the fertilization response needs isolating.

**Files**: `scripts/download_tempnle_npp.py` (new), `scripts/process_tempnle_npp.py` (new), `CLAUDE.md`, `config/isimip_search_catalog.yaml`, `WORKFLOW-ISSUES.md`.

### 2026-08-12: Customer Delivery Pipeline Built — External Review Found Eight Defects, One of Them in Shared Extraction Code

**Delivered**: the customer-delivery pipeline — `config/asset_catalog.yaml` (asset type → hazard layers), `config/layer_registry.yaml` (layer → disk location, status, which slope to read), `scripts/utils/delivery.py`, `scripts/generate_customer_delivery.py`, and a normalized CSV star schema (`locations` / `assets` / `layers` / `values`). The 28-column Looker `Export-Key.csv` contract is **retired** (user decision): no column may be derivable from another, so the hazard bands, the long-term-trend column (it *is* the final decade's slope under an expanding window) and the two never-emitted significance columns are gone.

**An external review (Codex) found eight defects in the first cut. All eight were reproduced independently before acting, and all eight were real.** The two that would have shipped wrong numbers:

1. **`spatial_extract.extract_by_point` did not wrap longitude at the antimeridian.** The search window masked on a raw `|lons - lon|`, so it was one-sided at the seam and **180°E and 180°W — the same meridian — returned different values**. Measured on `burntarea` ssp585 at 17°S, 2090s: **0.775 vs 0.962**; at 67°N the two answers were **62× apart** (2.37e-7 vs 1.47e-5). After the wrap both return 0.868, the average of the two one-sided answers. **This was a pre-existing defect in shared code**, so `extract_timber_locations.py` had been carrying it too. Inland values are bit-identical before and after.
2. **`data_status` scoped its domain mask to the delivery's own layers**, so a conifer-only delivery told an Amazon site at (−3, −60) it was `OUTSIDE_DOMAIN` — "offshore" — when the truth was that no conifer stand is modelled on perfectly good land. The mask is now the union over **every registry layer**, so a status cannot depend on what else the customer ordered.

The other six: `layers.csv` published one scenario's `n_members` as the layer's (`conifer-npp` is **10/9/10** across rcp26/60/85, so "10" is false for rcp60 — replaced by `n_members_by_scenario` + `members_by_scenario`); `slopes_agree` returned `False` on inactive 0/0 cells, which would make a downstream "unreliable trend" filter flag every quiet site (now `None` — OUTPUT-SPEC requires the active-cell view); non-numeric coordinates raised a bare `TypeError` past the CLI's error handler; the manifest carried no content hashes; QA review status was documented but never recorded (`qa_reviewed_on`, `null` on all five layers, surfaced as `NOT CONFIRMED`); and the percentile "regression check" was an orientation spot-check, not proof of passthrough.

**`scripts/test_customer_delivery.py` is the fix for that last one and the general guard.** It recomputes **every** delivered metric from the source NetCDF with a Gaussian weighting **written from the spec rather than imported** — calling the same function the delivery called would prove only determinism. Confirmed to fail on injected corruption: a double-inverted `conifer-npp` percentile, a ×10 `wildfire` slope and one wrong value produced **161 violations across 3,106 checks**.

**`recommended_slope` is measured per layer, not inferred.** `wildfire` needs `ols_slope`: Sen collapses to exactly 0 on **74.0%** of its active cells at the 2090s, despite a grid-level zero fraction of only 29.2% — the collapse lives in the year-pair *differences*, not the values. Only `conifer-npp` (2.1%) is safe on Sen. `--measure-slopes` re-checks this and flags registry disagreement.

**Spatial averaging was measured and documented end to end** (ASSET-CATALOG.md, "Spatial averaging — the complete picture"), because a delivered number passes through **two** kernels on some layers and nothing connected them. Stage 1 is the per-layer processing kernel (only `cyclone`, L=2.5, 8.1% on the centre cell). Stage 2 is extraction, and it is **a 4-cell blend, not a 3×3 smoother** — the 3×3 stencil requires the site to sit exactly on a cell centre, and **100% of 20,000 random sites got 2×2**, a 1°×1° footprint ≈111 km N-S. The kernel is truncated at 2σ (13.5% of peak, chopped not tapered), works in degree space so it reaches **2.6× further E-W than N-S at 67°N**, and drops NaN neighbours with renormalization — so a coastal site returns a land-derived value with `data_status` reading `OK` throughout. That last is **accepted**, not a gap: customer locations are masked to land upstream (user decision). Coordinate precision matters accordingly: moving a site ±0.25° changed 2090s burnt area at Shasta from **1.248 to 3.979**, 166% of the centre value.

**Rule created**: `ASSET-CATALOG.md` (new, the reference doc); `/customer-delivery` skill (new). An unknown asset type is an **error, never a default** — a silently defaulted bundle delivers a hazard with no transmission channel to the asset and the customer cannot tell it was a guess.

**Files**: `ASSET-CATALOG.md`, `.claude/skills/customer-delivery/SKILL.md`, `config/asset_catalog.yaml`, `config/layer_registry.yaml`, `scripts/utils/delivery.py`, `scripts/generate_customer_delivery.py`, `scripts/test_customer_delivery.py`, `scripts/utils/spatial_extract.py`, `CLAUDE.md`.

---

### 2026-08-13: Climate Score — A Composition Artifact That Reappeared at a Second Level of Rollup After Being "Fixed" at the First

**Delivered**: `climate_score.csv` — the unweighted mean of `percentile` across an asset's hazards, per forcing tier per decade. Percentile is the only cross-hazard comparable axis (`value` is in native units that differ per layer) and is already oriented so 100 is worst on every layer, `higher_is_better` ones included, which is what makes the mean legitimate.

**It is keyed on a forcing TIER, not a scenario code, and that is forced by the data.** No native ISIMIP code spans both rounds. Measured on a timber asset carrying three hazards: `rcp26` sees **1 of 3** (conifer-npp), `ssp126` sees **2 of 3** (drought-3b, wildfire), tier `low` sees **3 of 3**. A score keyed on `ssp126` would average two hazards and be labelled "across all hazards". This is the one place harmonization enters the CSV; `values.csv` still carries native codes only, and the `scenarios` column records which contributed. Scenario codes are averaged **within** a hazard before hazards are averaged, so a hazard contributing two codes to one tier is not double-weighted.

**The defect: averaging each tier over whatever assets it happens to have makes the tier lines describe different portfolios.** `cyclone` publishes rcp26/rcp60 and **no rcp85**, so every cyclone-carrying asset is unscoreable at the high tier. The portfolio high-tier 2020s mean read **39.9 against 42.1** for low and medium — impossible on a common basis, because the shared 2020s panel is bit-identical across scenarios and the tiers *must* be equal there. The high line was 2 assets; the others were 5.

**It then recurred at the next level down.** After the portfolio fix, adding a per-location Climate Score view reproduced it exactly: Shasta read **62.3 at high vs 51.7 at low/medium** in the 2020s, because that site holds a timber asset *and* a warehouse, and the warehouse carries cyclone. Same cause, one rollup level lower, found only because the baseline-equality tell was checked again.

**Fix applied**: every chart comparing across tiers or decades restricts to a **balanced panel** — assets complete in every cell of that grid — and states what it dropped (2 of 6 assets on the example delivery). Both levels now read 39.88 / 62.3-equivalent identically across all three tiers.

**Rule created**: **the shared 2020s panel is bit-identical across scenarios, so a baseline Climate Score MUST be equal across tiers for the same hazard set.** That invariant is now asserted in `test_customer_delivery.py` at *both* the raw-row level and on the balanced panels the dashboard builds, portfolio and per-location. **Apply the balanced panel at every level of rollup, not just the one where it first bit.**

**Two further ways the score can mislead, both now guarded:**

- **Incomplete hazard coverage.** ISIMIP3b layers have no 2010s panel, so a 2010s score rests on 1 hazard where the 2020s rests on 3 — a portfolio mean of **22 → 42** that reads as risk doubling when only the hazard set changed. Charts exclude incomplete-coverage decades and say so.
- **`n_hazards` alone cannot say whether a score is complete.** A 2-hazard warehouse and a 4-hazard timber asset both read 2, one complete and one partial, and the expected count lives in `assets.csv` — so a naive consumer sorting the CSV compares unlike scores. `n_hazards_expected` was added for exactly this. The manifest's `incomplete_rows` was also comparing every row against the **global** maximum, mislabelling every warehouse row in a mixed portfolio; it is now per-asset and reads **36** (12 from the 2010s + 3 cyclone-carrying assets × 8 decades).

**Files**: `scripts/utils/delivery.py`, `scripts/generate_delivery_dashboard.py`, `scripts/test_customer_delivery.py`, `ASSET-CATALOG.md`.

---

### 2026-08-13: Delivery Dashboard — Customer Text That Could Execute, a Silent Plotly Failure, and Two Phantom Bugs From Browser Cache

**Delivered**: `scripts/generate_delivery_dashboard.py` + `scripts/utils/viz_common.py` — a single self-contained HTML page (~110 KB) over the delivery CSVs: Climate Score stat tiles, a site map with a metric toggle, bar rollups, per-location small multiples, and a filtered table. **Deliberately NOT merged with `generate_maps.py`**: that renders gridded layers at ~70,849 SVG markers per panel where marker count is the binding constraint; this renders a point portfolio that embeds as JSON. Same conventions, opposite payload profile. Shared vocabulary lives in `viz_common.py`, which also absorbed the symmetric-limit pattern that had been open-coded repeatedly across two scripts at two different percentiles (95th in `generate_maps.py`, 98th in `compare_water_index.py`).

**1. Customer-controlled text reached the page as markup — the most serious defect of the session.** `json.dumps` escapes quotes and backslashes but **not** `<`, `>` or `&`. A delivery was generated with a location literally named `Depot </script><script>alert(1)</script>`; the string appeared **verbatim inside the `const DATA` block**, terminating the script element early and executing what followed. Raw `<b>markup</b>` from another name reached the DOM through `innerHTML`. Fix: the payload is `\uXXXX`-escaped server-side (valid JSON, decodes identically) and a JS `esc()` guards all 17 places data reaches markup, including Plotly `hovertemplate` strings. Re-tested with the same input: one script element, no raw markup, payload still decodes to the exact original names.

**2. Clearing a Plotly container's `innerHTML` without `Plotly.purge()` fails silently.** It wipes the DOM but leaves Plotly's state bound to the node, so the next `Plotly.react` does nothing — **a blank panel with no console error**. Introduced while "fixing" a placeholder, and it broke the Value/Percentile views entirely. The root cause was structural: two functions drew into the same `#series` div and the mode switch between them cleared the node by hand. Fixed structurally — one render path, one lifecycle, every DOM write behind a `resetSeries()` that purges first, and the Climate Score became a *panel* in the shared grid rather than a competing mode.

**3. A browser-cached dashboard is indistinguishable from a fresh one and produced two phantom bug reports** — colours and features reported missing that were correct on disk. Fixed with a **build stamp** in the page header (8-hex hash of payload + page logic + tokens) and `--check-stamp` to read what is on disk. **When a user says a regenerated page looks wrong, check the stamp before debugging the code.**

**A second external review found four more, all confirmed and fixed:**

- **The verifier validated the rows that were present, never that the expected rows existed.** Deleting every `rcp85` row from `values.csv` passed cleanly. It now reconstructs the expected key set from assets × layers × scenarios × decades and compares against manifest counts — the same truncation now fails four ways.
- **The Climate Score check imported the writer's tier mapping** "so they cannot drift", which is exactly backwards: a wrong mapping would both produce and validate the score. It now keeps a deliberate second copy. Verified by mistyping `ssp585` as `medium`: previously silent, now **68 violations**. *A test must restate what it expects; that is what makes it a test.*
- **A partial delivery looked shippable.** If the dashboard build raised, the folder kept all its CSVs and nothing marked it. Now the run records `dashboard: failed`, writes **`DELIVERY-INCOMPLETE.md`** (mirroring the `INVALID-DO-NOT-USE.md` convention) and exits non-zero; the verifier refuses any delivery carrying that marker.
- `record_stage()` read-modify-wrote `manifest.json` non-atomically (a crash mid-write left truncated JSON that broke every later stage — now tmp + replace), and the verifier opened `glob(...)[0]` rather than the manifest's recorded source path.

**Colour decisions** (user): red = risk (`percentile`, Climate Score), blue = raw magnitude; forcing tiers blue/yellow/red for low/medium/high. The tier triple is reference-palette slots 1/4/8, **not** the validated adjacent 1/2/3, and the palette validator could not be run — **no JS runtime, no Chrome on this machine**. What pure-Python luminance math could measure was measured: the red ramp is monotonic light→dark across all 13 steps; ordinal band steps clear the 2:1 floor on their own surface (2.01 lightest on light, 2.97 darkest on dark); and **blue vs red sit 1.12:1 apart on light and 1.05:1 on dark — effectively no lightness fallback**, with yellow at 2.11:1. Marker symbol *and* line dash are therefore applied **unconditionally** as secondary encoding. Do not strip them as visual noise; they are the accessibility channel.

**Known limitation, recorded rather than glossed**: the dashboard has **never been rendered in a browser in this environment**. Static checks cover the payload, element ids, delimiter balance and script integrity; they cannot see layout, overlap or colour. Any session that generates a dashboard and does not open it must report it as *unreviewed*, exactly as the layer workflow requires for maps.

**Regeneration check after the restructure**: CSVs byte-identical to the pre-restructure snapshot, decoded dashboard payloads equal, **18/18 reference values identical to full float precision** (including some measured before the antimeridian fix), and schema changes strictly additive.

**Rule created**: `/customer-delivery` skill owns the four-stage pipeline (1 inputs → 2 CSV **and** dashboard → 3 reports *(unbuilt)* → 4 caveats *(unbuilt)*); stage status is stamped into each delivery's `manifest.json`. Stage 2 is atomic — one command produces both artifacts, and a delivery missing either fails verification.

**Files**: `scripts/generate_delivery_dashboard.py`, `scripts/utils/viz_common.py`, `scripts/generate_maps.py`, `scripts/test_customer_delivery.py`, `.claude/skills/customer-delivery/SKILL.md`, `ASSET-CATALOG.md`.

---

### 2026-08-13: Stages 3 and 4 Built — The Dangerous Failure Is a Correct Report That Reads as a Complete One

Stages 3 (reports) and 4 (caveats) were built on the IFRS S2 spine, mapped outward to CDP 3.1.1, ESRS E1-9 and California SB 261. Four findings worth keeping.

**1. Coverage had to become a declared artifact, because omission is invisible.** A report that lists the hazards it assessed and stops reads as though the rest were assessed and found immaterial — every number correct, the document still misleading, and nothing in the pipeline able to detect it. `config/hazard_taxonomy.yaml` now enumerates the **19** physical-hazard families a disclosure is expected to address and records which each registry layer evidences. **Three are covered** (tropical cyclone, wildfire, drought). Riverine flood, coastal flood and extreme heat are absent, and for many built assets one of those is the dominant risk. A "hazards not assessed" section is now mandatory in every report.

**2. `conifer-npp` is not a hazard, and its percentile hides that.** It measures the stand's own productivity, is `higher_is_better`, and was inverted at processing time — so it reads on the same 1–100 risk axis as a hazard and would silently become a fourth one in every count. `vulnerability_frame()` excludes it via the taxonomy's `non_hazard_layers`. The Climate Score still averages it in, which is now a `must_disclose` caveat rather than an unstated choice.

**3. The vulnerability threshold is the single easiest way to make a true report misleading.** The count of vulnerable assets is a monotone function of a number somebody chose. Measured on the loblolly example: at percentile ≥ 60 **all five** tracts are vulnerable; at ≥ 80, **none**. The whole portfolio sits inside one band. So the threshold is disclosed in the same sentence as any count, the neighbouring thresholds are always reported, and `source: default` raises a caveat saying the threshold is ours rather than the customer's — the same protocol as `confirmed_on` in the asset catalog.

**4. Four leaks found by looking at the rendered output, not by testing.** Each was correct code producing an inappropriate document:

- **Multi-line HTML comments would have rendered.** The narrative scaffold hands the writer their facet-profile guidance as comments, and `markdown()` skipped only lines *starting* with `<!--` — so the second line onward of every guidance block, including the profiles' own "Do not claim" lists, would have appeared in the customer's report. Comments are now stripped in `parse_narrative()` and the verifier fails any report containing `<!--`.
- **Internal vocabulary in a filing.** `UNVERIFIED`, `this pipeline`, `config/…yaml`, `viz_common.RISK_BANDS`, and an anecdote about the withdrawn sugarcane layer all reached the rendered compliance report. Fixed by splitting the taxonomy into `customer_note` (rendered) and `materiality_note`/`blocker`/`isimip_candidate` (internal), filtering caveat evidence per semicolon-separated segment against a delivered-file allowlist, and trimming `source_dataset` to its first clause. The verifier now fails on a leak list.
- **`(nan)` in an asset label.** `str(x or "")` does not catch a missing optional column: `float('nan')` is **truthy**, so it survived to `str()` and printed. 
- **"1 of 2 hazards assessed" read as a site gap.** It was `cyclone` publishing no high-forcing scenario. Both affected sections now name the layers absent at the displayed tier — without it, a coastal warehouse appears never to have been checked for cyclones when it was, at the two tiers where the layer exists.

**5. A re-extract left stale reports looking finished.** `run_delivery()` rewrites `manifest.json`, which reset every downstream stage to `not_started` while the previous run's reports sat in the folder, openable and shippable. Downstream artifacts are now marked **`stale`** with the extract timestamp, and the verifier fails on a report whose stage is anything other than `built` — caught live, then fixed and re-caught.

**Rule created**: stage order is `inputs → extract → dashboard → **caveats** → compliance_report → bespoke_report`. Caveats runs BEFORE the reports because the caveat set is an *input* to them, not a summary: each report must carry every `must_disclose` entry, enforced at build and re-checked by the verifier. Generating them last would let two documents about one delivery disagree about what is wrong with it.

**Rule created**: a report contains exactly two kinds of sentence — derived from the delivery, or researched with a citation that **resolves**. `[data:…]` must name a real CSV row, `[dossier:…]` a source in `dossier.yaml`. The guarded failure mode is not laziness but **fluency**: a confident paragraph about a customer reads identically whether researched or invented. An unresolvable citation is worse than none, because it looks like evidence.

**Rule created**: facet profiles (asset × region × persona × vertical × use case × company) **guide** the narrative and are never pasted into it. Pasting would make every report of a given asset class identical while looking bespoke — generic output nobody catches.

**Verification**: 10 injected corruptions, all caught — must-disclose caveat removed, citation stripped, citation pointing at a nonexistent decade, slot reverted to `TODO`, threshold changed without rebuilding, vulnerable count altered in the HTML, internal vocabulary injected, HTML comment left in, manifest claiming a missing report, `caveats.json` deleted. Reference-site check per §12: the coastal NC tract carries **5.4×** the cyclone exposure of the piedmont tract (0.00332 vs 0.00062). Stage 2 outputs stayed byte-identical across the schema change. **4,418 checks** on the loblolly delivery, 3,301 on the forestry one.

**Standing limitation**: no browser can be driven here, so **no report has been rendered or printed**. Layout and pagination are unreviewed and must be reported as such. PDF is Safari's Print ▸ Save as PDF; the environment has no pandoc, weasyprint, wkhtmltopdf, kaleido or headless Chrome, which is why figures are inline SVG with zero JavaScript.

**Open**: every facet profile and every asset-catalog entry is `confirmed_on: null`. The `timber land` catalog entry omits `cyclone`, which is the dominant hazard for southeastern pine and irrelevant for Pacific Northwest pine — the catalog has no region dimension, so the loblolly delivery used a row-level override rather than changing it globally.

**Files**: `config/hazard_taxonomy.yaml`, `scripts/generate_delivery_caveats.py`, `scripts/generate_compliance_report.py`, `scripts/generate_bespoke_report.py`, `scripts/utils/report_common.py`, `scripts/utils/report_figures.py`, `scripts/utils/report_profiles.py`, `scripts/utils/delivery.py`, `scripts/test_customer_delivery.py`, `docs/reporting/`, `.claude/skills/customer-delivery/SKILL.md`.

---

### 2026-08-13: First Report Review — A Regional Narrative for the Wrong Region, a Metric That Should Not Have Been Published, and Two Figures Nobody Could Read

The first generated reports were reviewed by a user. Every finding was a failure of the same kind: the pipeline produced something **plausible** where it should have produced something **true or nothing**.

**1. The bespoke report described North Carolina for a portfolio spanning California, Alabama, Virginia and Bavaria.** The facet model allowed exactly one `region` profile, and nothing checked it against the delivered coordinates — so a single region's framing was applied to a two-continent portfolio and read as confident local expertise. Fixed structurally rather than by choosing better: **every facet now accepts a list, and `region` is validated against the data.** Each region profile declares a `matches:` block of countries and states, and `assert_region_coverage()` **refuses to build** when any location is uncovered. A region profile with no `matches.countries` is a load error, because it would match everything and silently defeat the check. Two bugs surfaced while testing it — an absent CSV column arrives as the float `nan` whose `str()` is the truthy `"nan"`, so every offshore site read as uncovered everywhere; and the pre-existing `us-southeast` profile had no `matches` block and therefore matched the entire world.

**2. The vulnerable-asset count was withdrawn.** IFRS S2 29(c) asks for the amount and percentage of assets *vulnerable* to physical risk. The implementation called an asset vulnerable when its worst hazard reached the 80th percentile — but **`percentile` is a global-relative EXPOSURE rank and "vulnerable" is a claim about susceptibility to HARM**, and nothing in the pipeline connects them. The instability was visible in the output and was published anyway: on one portfolio **every asset was vulnerable at a threshold of 60 and none at 80**. The metric had a defined method, a sensitivity table and an independent verification path — everything except a reason to believe it measured vulnerability.

**Rule created — a report section may report NOTHING, and should.** Deferred decisions live in `report_common.TBD_SECTIONS` and render through `tbd_block()`, which states the requirement, why it is deferred, the decisions outstanding, and what is reported instead. Both reports render the same block so they cannot describe one gap differently. The verifier **fails any report publishing a figure whose method is still deferred**, and flips automatically to "the published counts must match an independent recomputation" when the entry is removed. The pressure runs toward filling every box a framework defines, because a complete-looking report is what a customer expects; that pressure is exactly what produces a number nobody chose, and once it is in a filing it is indistinguishable from a reasoned one. **A gap is a conversation. A wrong number is a liability with the customer's name on it.**

**3. Two figures were unreadable, and both had been "verified" by well-formedness checks.** The XML parsed, the values were correct, and nobody could read them:

- **Bar labels clipped.** Site labels run to 54 characters ("Gulf Platform Alpha — Warehouse (Offshore supply base)") in a fixed 210-unit column, truncating precisely the part that distinguishes two assets at one site. Now wrapped to two lines with ellipsis beyond. SVG cannot measure text without a layout engine and there is no browser here, so `CHAR_WIDTH_EM` is an estimate deliberately erring narrow.
- **Trend strips rendered as columns of `0.000`.** Slopes span three orders of magnitude between layers — cyclone tops out near 9e-4 per decade, wildfire reaches 3e-1 — and fixed decimal places erased whole figures. Each strip now picks its own factor so its limit lands in [1, 10) and states it on the axis (`×10⁻⁴`).

**4. The markdown converter was silently DELETING author text.** Found only because the new region narrative used bold lead-ins. A paragraph opening `**Finding.** …` starts with an asterisk, was classified as a list item, rejected as a list (no space after the marker), rejected as a paragraph, and **dropped**. Multi-line list items were separately split into an `<li>` plus a stray `<p>`. The document looked finished and was missing arguments; every remaining number was still correct, so nothing else would have noticed. Fixed both, and added `assert_narrative_rendered()` — every source line of ≥25 characters must appear in the output, compared whitespace-insensitively on a fragment BETWEEN citation markers (stripping a marker leaves a gap the rendered text fills with a number, and replacing a tag with a space inserts spaces mid-sentence; both produced false misses before the comparison was made robust).

**5. Scope prose was correct and the example was not.** `report_compliance.html` already read "6 assets at 5 locations in Germany, USA" — the country list was derived, not hardcoded. What was wrong was that a *second*, invented single-region example delivery had become the thing being reviewed. It was removed; the worked example is now the actual mixed-asset, multi-region portfolio.

**What this round cost, and what it bought**: four new region profiles (`us-west-california`, `us-gulf-coast`, `us-mid-atlantic`, `central-europe`), three new asset profiles, and a new `annual-disclosure` use case — all grounded in cited sources, because the regional claims are the ones a reader checks. The Bavarian pairing turned out to be the strongest result in the portfolio: drought exposure diverges from 33.8 at low forcing to **97.4** at high, on the one hazard-and-asset combination with a documented mortality pathway (drought → bark beetle → spruce), while the productivity layer reads **8.2** at the same site and scenario because it contains no pest module. Both numbers correct, jointly meaningless without the caveat.

**Standing limitation, unchanged**: no report has been rendered in a browser here. The clipped labels and unreadable trend axes were both caught by a human opening the file — **static checks confirmed the SVG was well-formed and told us nothing about whether it could be read.**

**Files**: `scripts/utils/report_figures.py`, `scripts/utils/report_common.py`, `scripts/utils/report_profiles.py`, `scripts/generate_compliance_report.py`, `scripts/generate_bespoke_report.py`, `scripts/generate_delivery_caveats.py`, `scripts/test_customer_delivery.py`, `docs/reporting/`.

---

### 2026-08-13: External Review of the Report Tooling — Nine Confirmed Overstatements, Including Three Regulatory Claims That Were Simply Wrong

A Codex review scoped deliberately narrowly: *is anything asserted that is not true, and does anything overpromise?* Not general code review. The output of this pipeline is a document filed with a regulator, where a plausible false statement is worse than a bug because nothing downstream catches it. Nine confirmed findings, all accepted.

**The worst one contradicted the deferral we had just built.** The `ASSET-VALUES-ABSENT` caveat told the customer "this report therefore discloses counts and percentages of assets only" — in a `must_disclose` caveat, in *both* reports, three sections after the report states it publishes no such thing. The verifier passed 3,294 checks over it, because its deferral check matched five literal phrases and that sentence was not one of them.

**Three regulatory claims were wrong, not merely loose:**

- **IFRS S2 22(b)(iii)–(iv)** was cited for the "diverse range" and Paris-alignment requirements. Those are **22(b)(i)(2)** and **22(b)(i)(4)**; 22(b)(iii) is the reporting period and **22(b)(iv) does not exist**. Found independently while checking paragraph numbers against the published standard, and confirmed by the review.
- **22(a) was marked "Supplied".** It requires an assessment of the resilience of the entity's *strategy and business model* — uncertainties, financial flexibility, ability to redeploy assets. We compare hazard exposure across scenarios, which is an *input* to that assessment and not the assessment. Now marked NOT SUPPLIED.
- **ESRS E1-9 rows were marked "Yes"** for acute/chronic disaggregation, horizon breakdown and before-adaptation. E1-9's requirement is a **monetary amount and proportion of financial-statement assets**; classifying hazards and printing horizon labels does not supply it. All now marked No or Partly.

**"19 hazard families in the standard disclosure taxonomy" was an invented provenance.** The acute/chronic split *is* standards-derived; the family list and its boundaries are **ours**. TCFD Table A1 gives broad examples, not an enumeration; ESRS E1 AR 11 publishes a substantially longer classification. So "3 of 19" was a ratio against a denominator we chose, presented as though a standard had set it. Every fraction is now removed from customer-facing surfaces — reports name the hazards covered and not covered, which is true however the list is grouped.

**"Paris-aligned" was an inference presented as a fact.** IFRS S2 deliberately prescribes no scenarios, and a scenario *code* does not establish alignment with an agreement — that depends on temperature outcomes and the entity's own commitments. `PARIS_ALIGNED` is now `LOW_FORCING_PROXY`, the column reads "Low-forcing pathway present", and the report says the judgement is the entity's.

**Two source attributions were wrong in opposite directions.** "Roughly half the trees died at named northern Bavarian sites" came from a search-result summary and could not be traced to the paper it was cited against — the full text is paywalled (HTTP 403). **Removed** from the customer narrative rather than carried on an unverified attribution, with the removal recorded in the dossier's `verified_by`. Conversely, "29 million trees lost in California in 2015" is **true** but belonged to the USDA Forest Service Aerial Detection Survey, not to either source it was cited against; the correct source was added.

**Our own guarantees overpromised.** `report_common` claimed "there is no third kind" of sentence and that `check_citations()` "enforces it mechanically"; `generate_bespoke_report` claimed it made "impossible to publish judgement that is not sourced". What the code checks is narrower: one citation **per paragraph**, a data key that **exists** (not that the row contains the number quoted), a dossier id that **exists** (the source is never fetched). A paragraph can carry one resolving citation and several unsupported assertions and pass. Both docstrings now state the gap explicitly — an overstated guarantee about a guardrail is the same category of error as the thing the guardrail guards against.

**`customer_evidence()` matched substrings**, so `internal/archive/manifest.json` passed as a delivered file. Now matches bare filename tokens.

**Citation scope: five sentences quoted a 2020s value while citing only the 2090s row.** The checker passed them because the cited row exists — it never compares the number in the sentence to the row. Fixed in the narrative by citing both endpoints. The general gap is now documented rather than papered over.

**Rule created — a measurement quoted in a filing needs a reproduction, not a memory.** "Moving a site 0.25° changed a value by 166%" and "20,000 random sites, 100% resolve to a 4-cell blend" were both genuinely measured and neither had a script, a seed, or a retained artifact. To a reviewer that is indistinguishable from invention. `scripts/measure_extraction_sensitivity.py` now reproduces both on demand: the blend result confirms exactly (20,000/20,000), and the sensitivity figure is restated from a fresh measurement on the portfolio's own sites (**44% to 569%** at ±0.25°, wildfire ssp585 2090s) rather than an orphaned 166%.

Same reason, same round: the threshold-cliff figure ("all five vulnerable at 60, none at 80") was measured on a delivery that had since been **deleted**. Restated from the shipped example, where the count falls from 4 of 5 to 1 of 5 between thresholds 60 and 90.

**Verifier**: the deferral check now tests **text nodes rather than the flattened page** — flattening replaced tags with spaces, so "Sections 2, 4, 7" running into a cell reading "Assets vulnerable" produced a false positive on a report that published nothing. Four injected corruptions caught (contradicting caveat, restored table caption, per-asset determination cell, a count in prose), clean baseline.

**What the review could not check, and said so**: per-layer global slope measurements (`sen_slope == 0` on 74.0% of active cells, etc.) are reproducible from the processed NetCDFs via `--measure-slopes` but not from a six-site delivery. Fair; the receipt exists, one level up.

**Files**: `config/hazard_taxonomy.yaml`, `scripts/generate_compliance_report.py`, `scripts/generate_delivery_caveats.py`, `scripts/generate_bespoke_report.py`, `scripts/utils/report_common.py`, `scripts/test_customer_delivery.py`, `scripts/measure_extraction_sensitivity.py`, `docs/reporting/`.

---

### 2026-08-13: Crop-Failure Ingest — A Recorded Negative That Blocked a Hazard for Weeks, a Product That Zero-Fills the Ocean, and a Contract Variable I Nearly Published as a Constant

Ingesting ISIMIP3b `cropfailure` (Zantout2025). Four findings, in descending order of how
long they would have gone unnoticed.

**1. Our own catalog understated the hazard, and the understatement was the blocker.**
`search_results.drought.exposure_lange2020.family.lec` read `impact_models: "gepic"` — one
model. `config/hazard_taxonomy.yaml` then carried that forward as the *stated reason* crop
failure was unassessable: *"`lec` is a single impact model, so its CI would carry no
structural uncertainty."* Both were wrong. Lange2020 publishes `lec` from **GEPIC and
PEPIC**, and the ISIMIP3b re-issue has **eight** crop models × 5 GCMs × 3 SSPs in a complete
120-file matrix — the **deepest ensemble of any layer in this product**. The bad figure came
from a 2026-07-24 API listing capped at 20 rows, recorded without its cap noted.

This is GUARDRAILS §11 in a new costume. The rule was written for *absence* claims; this was
an **understated positive**, which is more dangerous, because a negative at least looks like
something to re-check. `impact_models: "gepic"` reads as a measured fact. Nothing in the
catalog format distinguishes "enumerated exhaustively" from "whatever the first page of a
truncated listing showed" — and §8 already says `count=1001` means truncation, but says
nothing about a *row cap* on a listing, which is the same failure at a different limit.

**Rule reinforced**: a capped or paginated result is `UNVERIFIED` for **coverage in either
direction**, not just for absence. Record the enumeration method beside the count.

**2. The publisher zero-fills the entire globe, and `isfinite` is therefore meaningless.**
Every one of the 120 members is non-NaN over ~100% of the 259,200-cell grid. Ocean,
Antarctica, Greenland and the Sahara all read exact `0`, not NaN. Using the finite mask as
the land mask — which every other processor in this repo does, correctly, for its own
inputs — would have put **87% ocean zeros into the percentile baseline population**,
reported `n_members=40` over open ocean, and handed every ocean cell the lowest-risk
percentile. The footprint had to be derived from where the field is ever **non-zero**:
39,890 cells, 42.4% of land.

The near-miss worth recording is the *diagnosis*, not the discovery. A 100%-coverage finite
mask is also the signature of the known `floodedarea` ocean leak, and the reflex was to
call it that and mask to the generic land-sea mask. It is **not** that defect: 39,406 of
39,890 footprint cells (98.8%) fall inside the ISIMIP3b land mask, and of the 484 that do
not, 63% are directly adjacent to land and the rest are **small islands the 0.5° generic
mask does not resolve** — Lofoten (68.25 N, 13.25 E), Shetland (60.25 N, −1.75 E). They come
from `ldndc` (476) and `lpjml` (453) only. Treating them as an ocean leak would have deleted
real island cropland to satisfy a coarser product. Two failure modes, identical symptom,
opposite correct response — separated by asking *where* the off-mask cells were, not *how
many*.

**3. `n_members` and `n_models` would have shipped as constants.** Caught in the middle of
the first full build, which was killed before it wrote. Because every member is finite
everywhere, `isfinite(...).any(...)` gives **40 members and 8 models in every published
cell** — while the same file's `mask_rule` attribute tells the reader *"a 1-model cell
carries no inter-model uncertainty in its CI — filter on `n_models` if that matters."* The
file would have contradicted its own instructions, and the contract variables would have
been unusable for exactly the marginal cells they exist to flag.

Fixed by masking each member to its **own** cropland footprint before pooling, so a model
that grows nothing in a cell contributes no observation rather than a structural zero.
Measured: only **38.1%** of footprint cells have all 40 members simulating crops (1–10
members: 563 cells; 11–30: 5,471; 31–40: 33,856). Effect on the value is **1.03×** globally
(0.00535 → 0.00550) and **4.4×** on the 1–10-member cells — negligible for a headline,
decisive for a site-level query in marginal cropland, which is where such a query is
already most fragile. After the fix `n_members` spans 1–40 and `n_models` 1–8.

The general shape: **a zero-filled product silently converts "no opinion" into "zero
risk"**, and the conversion is invisible in every algebraic check because the array is
full, finite and correctly shaped. `test_shared_baseline.py` would have passed.

**4. The time axis declares no units at all.** `time` carries `long_name`, `standard_name`
and `axis` — and **no `units`** — so xarray cannot decode it and `.dt` raises
`AttributeError`. The values are bare integers `165..250` for a file named `..._2015_2100`,
i.e. `years since 1850`, undeclared. Decoded from the filename span with assertions on
length and unit stepping rather than a hardcoded epoch. Its 3b sibling Heinicke2026
(`driedarea`) has a properly declared `days since 1601-01-01`, so **time-axis handling is
per-PUBLICATION, not per-round** — the same lesson the filename grammar already taught
(Zantout2025 carries a leading `zantout2025_` token; Heinicke2026 carries none).

**Also**: Zantout2025 publishes **no `.json` sidecars**, so there is no upstream sha512.
Integrity is `Content-Length` only, and the sha512 in `download_provenance.csv` is computed
locally — a self-consistency receipt, not publisher confirmation. The provenance file now
carries a `sha512_source` column saying so, because a digest column that means two different
things in two layers is worse than no column.

**GUARDRAILS §12, both halves, passed**: all 12 named cropland reference sites are non-zero
and non-NaN in all 120 raw members *and* in the processed output, rising 2020s→2090s at 11
of 12. Contact sheets were **rendered and looked at** — 40 members and 8 decade panels — via
a new `scripts/render_contact_sheet.py`, which writes PNG directly from numpy with zlib
because this venv has neither matplotlib nor Pillow. Geography coherent, no block structure,
no seams, and the visible model ordering matches the measured 6.69× inter-model spread.

**A reading that will be challenged, and is correct**: the layer ranks **Iowa at the 99.3rd
percentile of cropland and the Sahel at 69.4**. `cropfailure` measures *unprecedentedness*
against a fixed local reference, so a reliable breadbasket with a tight historical
distribution trips its own threshold more readily than a chronically marginal region. It is
the crop analogue of the drought layer's "departure, not aridity". This layer must never be
read as a food-security or absolute-yield map, and `interpretation_caveat` says so.

**Open**: the index carries **no crop token** and its sidecar's `specifiers` block names only the variable, so
*which* crops are aggregated into it, and how they are weighted, is not readable from the
archive. That is material for site-level use — a crop-aggregated index over crops a site
does not grow is not that site's risk — and it is recorded as
`crop_composition_undocumented` rather than guessed.

**Two rules came out of the same session and are recorded with it**, because both were
"written down but not wired up":

- **A layer scoring departure from a local baseline now declares it and the machinery
  enforces it.** `relative_baseline: true` + `relative_baseline_note` in the registry →
  a `must_disclose` caveat both reports refuse to render without. `cropfailure-3b` ranks
  **Iowa at the 99.3rd percentile of world cropland and the Sahel at 69.4** — correct, and
  the first thing a reader challenges. The two drought layers were already in this class and
  their own `delivery_note` said *"Must be stated in any customer narrative"* while the
  generator filed it as `should_note`. Flagging cropfailure alone would have implied the
  drought layers were absolute, so all three were flagged and the example delivery
  regenerated (11 → **12** must-disclose caveats, 3,294 checks, text confirmed present in
  both rendered documents). Moving the substance into its own field also silently dropped it
  from the Stage 1 plan — the one place the mapping is agreed with the customer — so the plan
  now prints it as `READ AS RELATIVE:` above the ordinary `NOTE:`.
- **Per-dataset facts moved out of CLAUDE.md into [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)**
  (user instruction). CLAUDE.md had already declared *"per-hazard framing decisions live with
  their layer, not here"* and then carried a shipped-layer table, a Water Index variable
  table, and a dozen per-layer measurements. The rule now has somewhere to point. Per-layer
  slope and smoothing tables also moved out of ASSET-CATALOG.md so there is **one** copy that
  cannot drift, and `OUTPUT-SPEC.md`'s stale *"Only `let` qualifies today"* was corrected.
  Where a fact belongs: a rule that generalizes → CLAUDE.md/GUARDRAILS; a fact about one
  dataset → DATASET-ATTRIBUTES + processor docstring; a fact the processed file already
  carries → **nowhere**, it is read from the file's global attributes at delivery time.

**Files**: `scripts/download_cropfailure_isimip3b.py`, `scripts/check_cropfailure_nature.py`,
`scripts/process_cropfailure_isimip3b.py`, `scripts/render_contact_sheet.py`,
`scripts/utils/delivery.py`, `scripts/generate_delivery_caveats.py`,
`scripts/generate_customer_delivery.py`, `config/isimip_search_catalog.yaml`,
`config/hazard_taxonomy.yaml`, `config/layer_registry.yaml`, `config/asset_catalog.yaml`,
`DATASET-ATTRIBUTES.md`, `CLAUDE.md`, `ASSET-CATALOG.md`, `OUTPUT-SPEC.md`, `GUARDRAILS.md`,
`.claude/skills/isimip-search-download/SKILL.md`,
`.claude/skills/isimip-process-visualize/SKILL.md`.

---

### 2026-08-14: The CaMa-Flood Suite We Never Saw — Three Weeks of Treating a Known Path as a Search

**Outcome first**: the best river-flood data in the ISIMIP repository — hydrodynamic inundation at up to **150 arcsec (~4.6 km, 144× our grid)**, with flood **depth in metres**, protection-level variants, and **rcp85** — sat one directory level above a path this repo had already walked, for **three weeks**, while `config/hazard_taxonomy.yaml` recorded a two-item candidate list and a blocker for the family it calls *the most consequential gap in this pipeline*. It surfaced only because the user asked "is there not a higher resolution option that speaks to actual inundation risk?"

**What happened**:

1. On **2026-08-08** the Lange2020 exposure family was enumerated as `ISIMIP2b/DerivedOutputData/Lange2020/{MODEL}/{gcm}/future/` — the catalog records that exact path. We reached it from a variable-code search (`le{d,r,w,c,h,t}`) and walked *down* correctly, listing every level below.
2. **`ISIMIP2b/DerivedOutputData/` — the parent — was never listed.** It has two entries. The second is `Zimmer2023`: CaMa-Flood v3.6.2, `flddph`/`fldfrc`, 6 GHMs × 4 GCMs × rcp26/60/85 × 3 protection levels = 432 files.
3. On **2026-08-11** the *3b* `DerivedOutputData/` root **was** listed — because someone asked a newer-round question about drought. The rule in the skill read "list `DerivedOutputData/` before concluding a product has no newer-round version", so it fired for 3b and never for 2b, where we believed we already had the answer.
4. On **2026-08-13** `hazard_taxonomy.yaml` was seeded from the catalog, hardening the two-item list into a decision-grade artifact with a written blocker. Nothing in that file could express "this list has never been checked against a directory listing".
5. On **2026-08-13**, reviewing flood options, I enumerated only the two products the candidate list named. I inherited the blind spot and re-published it.
6. **The correction then repeated the mistake.** Having found `Zimmer2023`, I listed the 3a and 3b roots, read the *publication names*, and wrote into the catalog: *"No CaMa-Flood publication exists under ISIMIP3a or 3b."* Wrong twice. `Quesada-Chacon2026` (3a) **is** CaMa-Flood, and `TipESM2025` (3b) is CaMa-Flood `fldfrcmax` at 15 arcmin with **ssp126/370/585 and 32 members per scenario** — which the catalog had recorded as "water models" since 2026-08-11, from a note that identified its model directories and never opened a file.

**Impact**: no bad data shipped — flood is not a delivered layer. What was damaged is worse in kind: the *option set*. A hazard the taxonomy calls the pipeline's most consequential gap was deferred against an incomplete list of what exists, and the reason it looked complete is that every individual fact in it was true.

**Root cause — five distinct failures, all of which read green**:

- **Entering a known path is retrieval, not enumeration.** Every enumeration rule we had (§8 "list every intermediate level", §11 "a negative needs a receipt") is triggered by *doubt* — by asking whether something is absent. Nothing was triggered by *touching a directory*. So the level we started at was never listed.
- **A candidate list is an unreceipted positive.** §11 already covers absence claims and understated positives; an option list is the third form and the most decision-bearing, since work is planned off it. It had no `enumerated_on` field and no way to say "incomplete".
- **A variable-code search cannot discover a vocabulary you do not have.** `flddph`, `fldfrc`, `fldfrcmax` are unreachable from any `le*` query. Directory listing is the only discovery mechanism; code search only confirms.
- **Publications are named after people and identified by their files.** Three of the nine derived publications across all four rounds are CaMa-Flood, and **not one says so above the file level**. Two of our three mischaracterisations came from reading directory names.
- **Resolution was never a search dimension.** Every layer we ship is 0.5°, so "is there a finer product?" was not a question the workflow could ask. One hazard spans 0.5°, 15 arcmin, 150 arcsec, 300 arcsec and 900 arcsec.

**Rule created — GUARDRAILS §13, "List the Container You Entered — a Known Path Is Not a Search"**: list the parent before recording the child, regardless of the question; identify a publication by its files, never its name or model list; a candidate list needs the same receipt a negative does; resolution and filename grammar are per-publication properties; a code search is confirmation, not discovery.

**Applied, not just written**: swept all four rounds' `DerivedOutputData/` roots and opened every publication to file level (2026-08-14). The result is `repository_publication_map` in the catalog — nine publications, what each actually publishes, the non-uniform directory depths, and an explicit `known_gaps_in_this_map` list naming what the sweep did *not* cover (2a OutputData sectors, `InputData/` roots, Jaegermeyr2021's crop coverage). The skill gains the map as a table plus the refresh command; `hazard_taxonomy.yaml` gains `candidates_enumerated_on`; the false "no CaMa-Flood in 3a/3b" line is corrected in place rather than deleted.

**Also caught by the sweep**: `TipESM2025`'s directory is **five** levels deep with the scenario **above** the model (`MIT/{scenario}/{GHM}/{gcm}/`), unlike every other publication's `{MODEL}/{gcm}/{period}/`. A three-level walker returns it as empty — which under §8 is a wrong-depth signal, not a negative. And `Zimmer2023` filenames carry two extra fields (protection level, resolution), so the standard `$(NF-4)` projection reports a variable vocabulary of `150arcsec`; the fix is to read the `.json` sidecar's `specifiers` block, which names every field.

**Files**: `GUARDRAILS.md` (§13, new), `.claude/skills/isimip-search-download/SKILL.md` (publication map + refresh procedure + resolution as a search dimension + the narrow newer-round rule rewritten + the mid-filename grammar warning), `config/isimip_search_catalog.yaml` (`repository_publication_map`, `search_results.river_flood` with all three CaMa-Flood editions, corrected 3a/3b claim), `config/hazard_taxonomy.yaml` (`candidates_enumerated_on`, flood candidates), `CLAUDE.md`, `WORKFLOW-ISSUES.md`.

---

### 2026-08-14: Heatwave Ingest — A Sidecar Trap That Downgraded a Shipped Layer's Provenance, a Layer Whose Flat Trend Means the Opposite, and an Index That Changed Between Rounds

**What happened**: Three separate findings while ingesting both heatwave sources. Two are traps that will recur; one is a layer property that breaks a standing rule.

**1. `{stem}.nc.json` 404s; the sidecar is at `{stem}.json`.** The crop-failure ingest (2026-08-13, same `Zantout2025` publication) recorded *"NO `.json` SIDECARS … `{stem}.nc.json` returns 404, verified"* and fell back to `Content-Length` for integrity, stamping `sha512_source: computed-locally (no upstream sidecar)` into `download_provenance.csv`. **Zantout2025 does publish sidecars**, with `size`, sha512 and a full `netcdf_header` — the probe just used the wrong suffix. Verified 2026-08-14 on both a `heatwave` and a `cropfailure` file: `.nc.json` → 404, `.json` → 200. The directory listing shows the `.json` files plainly.

**Impact**: `cropfailure-3b`'s 120 files were verified against `Content-Length` only, which catches truncation but not silent corruption, and its provenance actively *records* that no publisher checksum existed. That is a self-inflicted downgrade on a shipped layer, and it reads to a reviewer as a property of the publisher rather than a bug in our probe.

**Root cause**: A negative was drawn from **one** URL shape without listing the directory that would have shown the answer. GUARDRAILS §11 says a negative needs its receipt; the receipt here named a path that was never going to exist.

**Correct action**: Before recording "this publication has no sidecars", list the directory and look. One request.

**Fix applied**: `download_heatwave_isimip3b.py` and `download_leh_isimip2b.py` both verify against the publisher sidecar (15/15 and 8/8 sha512-matched). The trap is recorded in `search_results.heatwave.representations.heatwave_exposure_3b.sidecars`. `download_cropfailure_isimip3b.py` was independently corrected the same day and now verifies against `{stem}.json`; the catalog records the extension-REPLACING naming under `lec_isimip2b.sidecars` and `floodedarea.sidecars`. Whether `cropfailure-3b`'s existing `download_provenance.csv` was regenerated against the publisher checksums is a separate question for that layer's owner — the rows written on 2026-08-13 still say `computed-locally`.

**2. `heatwave-3b` saturates, and a flat slope there means the opposite of what it means everywhere else.** Exposure is the annual HWMId exceeding the 97.5th percentile of *that cell's own* preindustrial control, so warming pushes cells permanently over the threshold and the binary flag pins at 1. On ssp585 the 2090s panel has **45.9% of published cells at exactly 1.0** (ssp370 32.6%, ssp126 1.2%). At 1.0 the pooled sample has no variance, so the CI collapses to zero width, **both** slopes go to ~0 and **agree** there, and the percentile ties at 100 (51.9% of cells at ≥99.5).

This breaks the standing dual-slope rule from OUTPUT-SPEC — that the estimators fail in *opposite* regimes so their disagreement flags a fragile trend. Here they fail the *same* way, and agreement near zero is not reassurance. The censoring inverts regional rankings: the Amazon's `ols_slope` **falls** +0.160 → +0.046 dec⁻¹ as it goes 0% → 100% saturated, while never-saturating Siberia **rises** +0.069 → +0.098, so the 2090s panel has Siberia out-trending the Amazon 2.1×.

**Correct action**: Do not invent a statistic to hide it. The contract fields are emitted as specified and the censoring is *declared* — `saturation_caveat` and `saturation_by_decade` in every scenario file, plus a worked example — and saturated cell-decades are identifiable from the published fields alone (`median == 1.0`, equivalently a zero-width CI at 1.0). Tie-ranking and a possible time-to-saturation companion are left as **open decisions** rather than answered unilaterally.

**3. The index changed between rounds; the layers are not versions of each other.** ISIMIP2b `HWMId-humidex` requires the relative HWMId **AND** an absolute Humidex ≥ 45, so counted events "would also adversely affect human health" (Lange et al. 2020). ISIMIP3b `HWMID-NONE` drops the absolute criterion — the "NONE" *is* the humidity term. Reading `HWMID-NONE` as a newer `HWMId-humidex` would have shipped a temperature-only relative index under copy promising workforce-safety and equipment-derating relevance.

Measured consequence, and it is large in both directions. `leh` under its absolute gate is **silent on 59.7% of land** across 2006–2099 and both RCPs — 81% of the 35–60°N temperate belt and 99% of the boreal belt never register a single exposed year in any member. Paris, Frankfurt and Yakutsk are **exactly zero in every member and every year**. That is the absolute threshold behaving correctly, not a defect, but a layer built on it reports *no heat hazard at all* for much of a Northern-Hemisphere portfolio, and "zero" there means "never crosses Humidex 45", not "low heat risk".

**4. A land-mask variable name assumed instead of resolved.** `leh` zero-fills the entire 259,200-cell grid, so `isfinite` is not a mask. The first sparsity run tried to load the ISIMIP2b land mask by the variable name `mask` — which is the **3b** name; 2b publishes `LSM` — the load threw, the check silently fell back to the finite domain, and every sparsity share came out quoted against the whole globe (89.1% "silent") instead of against land (59.7%). Same dilution family as the ocean-diluted `sen_slope == 0` misreport of 2026-08-10.

**Fix applied**: `check_leh_nature.py` resolves the mask variable name per round and prints the denominator it used.

**5. A documented status value that the loader could not load.** `config/layer_registry.yaml` has always documented `status: preferred | alternate | superseded | blocked`, but no layer had ever used `blocked`, and `LayerSpec` had no field for the reason. Registering `heatwave-3b` as blocked with a `blocked_reason` raised `TypeError: __init__() got an unexpected keyword argument` and took the whole delivery planner down until `LayerSpec` gained the field. A vocabulary documented in a header is not a vocabulary the code accepts; the first user of an unexercised branch pays for it.

**Rule created**: none new — this is GUARDRAILS §8 (never guess a specifier), §9 (measure, do not infer) and §11 (a negative needs a receipt naming what was actually enumerated) landing on three new surfaces: a sidecar URL shape, an index name, and a mask variable name. But OUTPUT-SPEC.md's dual-slope guidance IS amended: it now records that a CENSORED field breaks the "they fail in opposite regimes" premise, because at a bound both estimators fail the same way and agree at zero.

**Files**: `scripts/download_heatwave_isimip3b.py`, `scripts/check_heatwave_nature.py`, `scripts/process_heatwave_isimip3b.py`, `scripts/download_leh_isimip2b.py`, `scripts/check_leh_nature.py` (all new), `config/isimip_search_catalog.yaml` (`search_results.heatwave`), `config/hazard_taxonomy.yaml`, `DATASET-ATTRIBUTES.md`.

---

### 2026-08-14: Permafrost Layer — Four Estimator Traps, Three of Them Mine, and a Model I Nearly Dropped for the Wrong Reason

**What happened**: Building `permafrost-3b` from ISIMIP3b `thawdepth` surfaced four separate ways a correct-looking measurement produced a wrong conclusion. Each was caught, but three of them had already been written into a recommendation before the check that overturned them.

**1. A sector-only walk would have under-counted the ensemble by two thirds.** `thawdepth` is published byte-identically under `permafrost`, `biomes` and `water_global` (Content-Length **and** ETag matched across sectors). The `permafrost` sector holds two models and only **one** of them publishes the variable; JULES-ES-VN6P3 and CLASSIC publish it in `biomes` and are absent from the permafrost sector entirely. Walking the obviously-named sector answers *"1 model, 5 members"* to a question whose answer is *3 models, 12 members*.

**Root cause**: the standing "walk one sector, ingest one sector" rule (a *duplication* rule) was read as a rule about where to **look**. It is not. Recorded as `but_look_everywhere` under `repository_structure_cache.cross_sector_duplication`.

**2. An absolute tolerance on a ceiling test inverted a finding.** A permafrost-free cell is published **at the model's soil column depth** — Nairobi and Paris read exactly the column in all three models — so detecting "fully thawed" means testing equality with that ceiling. The first check used `atol=1e-4` and reported LPJmL's at-ceiling share as **0.0%** while that model's *median* sat on the ceiling: its maximum is 13.001 m and the pinned mass is at 13.000, a 1 mm gap. The same bug made its permafrost domain read as the full finite extent. Any ceiling test on a float field needs slack **proportional to the ceiling**.

**3. A split-half smoothing test that measured model composition instead of sampling noise.** The halves were taken as `members[::2]` / `members[1::2]` off a *sorted* member list, which gave one half 3 JULES + 2 LPJmL and the other 2 + 3. Different composition → different permafrost domains → Pearson **r = 0.376**, and each half read *smoother* than the full ensemble, which is the tell (fewer members should be noisier). Stratifying by model — alternating GCMs **within** each model — gave r = 0.784 with roughness stable across halves. The verdict text was also written as a conclusion before the number existed; it is now generated from the measurement with an explicit three-way branch.

**4. The one I got most wrong: recommending a model be dropped on a normalised-space artifact.** After per-model column normalisation JULES sits at 0.95 of its 3 m column while CLASSIC and LPJmL sit at 0.035–0.046, so JULES looks censored and I recommended dropping it. Checked in **raw metres** against observed active-layer thickness (~0.3–1.5 m in continuous permafrost), the ranking reverses:

| model | column | 2020s thaw p50 | p95 | Fairbanks |
|---|---|---|---|---|
| LPJmL5-7-10-fire | 13.0 m | **0.83 m** | 4.57 | 1.49 m |
| JULES-ES-VN6P3 | 3.0 m | 2.85 m | 3.00 | 3.00 m (at its column) |
| CLASSIC | 61.4 m | 2.15 m | **28.0 m** | **61.4 m — no permafrost** |

The model with the most headroom is the least physical: a 28 m p95 seasonal thaw is not an active layer, and CLASSIC reports no permafrost at Fairbanks, which sits in the discontinuous zone. JULES's baseline pattern also correlates with the other two as well as they correlate with each other (ρ 0.785 / 0.657 against 0.631). **"Has room to move" is not evidence of quality**, and a member that looks degenerate after a transform should be re-checked in the units it was published in. The user proposed keeping JULES — for a reason that was itself inverted (its column is the *shallowest*, not the deepest) — and the conclusion was right where mine was wrong.

**Impact**: no shipped artifact was affected — every one of these was caught before the layer was registered — but items 2–4 had each reached a stated recommendation before being overturned, and item 4 would have removed 5 of 12 members and the model closest to the ensemble's centre.

**5. A threshold on the central value is a property of the estimator, not of the data.** The layer reports how much of the soil column transitions from permafrost to none, so "area whose column is fully thawed" looks like a natural headline. It is not stable: under a pooled **median** the criterion is *>half the members*, under a pooled **mean** it is *effectively all of them* — **7.97 vs 0.40 M km²** for the same ensemble at ssp585 2090s. Both numbers are correct implementations of different questions. The invariant is the **member share**: 0.478 of members on average by the 2090s, 10.52 M km² where at least half agree, 0.69 M km² where all twelve do.

**6. Theil-Sen is biased LOW on a multimodal ensemble — a regime the two-slope table did not have.** `sen==0` on only 15.2% of active cells with 84.4% sign agreement reads as a healthy field, and Sen is still wrong here: the pairwise sample is dominated by *cross-cluster* member pairs whose slopes carry the level offset rather than the trend. Measured on ssp585, the 12 members' own OLS trends span +0.0104…+0.0826 dec⁻¹ with a mean of **+0.0326**; the published `ols_slope` is **+0.0326** (exact) and `sen_slope` is **+0.0069 — below every single member**. An estimator outside the range of the things it summarises is not robust. Caught only because `generate_customer_delivery.py --measure-slopes` flagged `REGISTRY DISAGREES` against a `sen_slope` entry that had been written from the *ols-bias* regime by analogy instead of measurement.

**Rules created**: OUTPUT-SPEC gains a **fourth decadal-statistic branch** (`pooled_mean_multimodel`, with the conditions distinguishing it from reaching for the mean to improve contrast) and a **fifth row in the two-slope table** (multimodal ensemble → read `ols_slope`). CLAUDE.md's "three branches" rule is updated to four and now carries the threshold-inherits-the-estimator warning.

**Files**: `scripts/download_thawdepth_isimip3b.py`, `scripts/check_thawdepth_nature.py`, `scripts/process_thawdepth_permafrost.py` (all new), `OUTPUT-SPEC.md`, `CLAUDE.md`, `DATASET-ATTRIBUTES.md`, `config/layer_registry.yaml`, `config/hazard_taxonomy.yaml`, `config/isimip_search_catalog.yaml`.

---

### 2026-08-15: `csoil` rebuilt — a lost artifact, an ensemble understated by a third, and a parallel port that OOM-killed itself twice

**What happened**: a review of soil-degradation coverage asked whether `csoil-total` had been processed to current standards. It had not, and the first finding was blunter than a standards gap: **the layer did not exist**. `data/processed/soilcarbon_csoil_annual/` was an empty directory, no raw remained, and there had never been a `download_csoil*.py` — while `config/isimip_search_catalog.yaml` still named the three output files and DATASET-ATTRIBUTES.md still said "processed but not registered". A `processed:` date is not evidence the artifact is still on disk; `data/` is local and ephemeral.

**1. The recorded ensemble was wrong by omission — 12 members where 17 exist.** The catalog listed `csoil_annual: [classic, jules-es-vn6p3, mc2-usfs]`. **LPJmL5-7-10-fire publishes `csoil-total` annually for all three SSPs** and was simply missing: +5 members and +1 structurally independent model. This is GUARDRAILS §11's understated positive, and it is the more dangerous half of that rule — a recorded negative at least looks like something to re-check, while a model list reads as a measured fact. The search skill's own notes already warned that this same list "later still missed CLM45 and VEGAS for csoil" in ISIMIP2b; the identical failure had recurred in 3b and nobody reconciled it.

**The sector trap was NOT the cause, and that is the lesson.** The obvious hypothesis after `thawdepth` was that a sibling sector held the missing model. It did not — `permafrost` holds only ELM-ECA and LPJmL, `fire` is a subset of `biomes`, `water_global` has no csoil publisher. **The omission was inside the sector that had already been walked**, which no completeness argument about sectors could ever have caught. Only projecting the variable field over a *full* listing found it. A sector-coverage claim is not a variable-coverage claim.

**2. Every framing decision from the 3-model build had to be re-taken, and one of them did not survive.** The old rationale for declining normalization was "comparable magnitudes, medians within ~2×". With LPJmL the spread is **2.92×** (5.70 / 7.55 / 10.45 / 16.82 kg m⁻²) and LPJmL carries the upper tail alone (p95 91.1 against 16.5–27.4). Normalization is still declined, but on different grounds — a gradient rather than clusters, and physically meaningful units that can be checked against soil inventories. The other decisions were re-measured rather than inherited: `pooled_median` retained because the largest adjacent-model gap is **0.59 × IQR** (branch 4 needs ~1.0), and smoothing declined on a **model-stratified split-half** — roughness 0.148, halves 0.149/0.150, **r = 0.992**.

**3. The layer's defining caveat was written where nothing would read it.** The mixed-CO₂ limitation — JULES publishes only the fixed-2015-CO₂ run, so 5 of 17 members carry no fertilisation signal — lived in a `co2_treatment` attribute. `LAYER_ATTRS_EXPORTED` in `scripts/utils/delivery.py` is a **closed allowlist** and `co2_treatment` is not on it, so the text would have been silently dropped before `layers.csv`, the caveat generator and both reports. Moved to `interpretation_caveat`, which is on the list. **Writing a caveat is not the same as delivering one**, and the allowlist is the thing to check.

**4. The mask is time-varying, by one cell.** The mask-agreement guardrail refused to write ssp585: at the 2060s panel one cell had observations early in the expanding slope window and none inside the decade, producing a finite slope against a NaN median. Fixed per OUTPUT-SPEC by masking slopes to each decade's own median mask. One cell is trivial against `npp-tempnle`'s 53→374 — and it is the same defect, on a layer nobody expected to be in that regime at all.

**5. The parallel port OOM-killed the run twice, because the memory budget was per-worker.** Single-core Theil-Sen at 17 members projected to ~7 hours, so the pool machinery from `process_thawdepth_permafrost.py` was ported. `slope_chunk_cells` sizes scratch from a fixed 400 MB constant — which **every forked worker applies independently**, so 8 workers meant 3.2 GB of slope scratch on top of a 1.24 GB cube, on a 16 GB machine with ~5 GB free. Both runs died at the identical 2030s→2040s transition. The budget is now divided by the job count so it is a total, and a resumed run frees the scenarios it is not writing: measured after the fix, processor plus 6 workers hold **0.28 GB RSS**. A `--scenarios` flag was added that narrows only the *write* loop — every scenario is still loaded and pooled into the shared baseline, because narrowing the baseline instead would silently break cross-scenario comparability.

**The first of those two kills was self-inflicted, and worth recording separately.** Waiting on the run with chained `until … done` shell loops left three stale pollers running, one for 72 minutes against a *previous, failed* log. Stopping one of them killed the processor it referenced. Background work should be waited on by ending the turn — the harness re-invokes on completion — and anything long should be detached (`nohup`/`disown`) so a watcher's death cannot reach it.

**6. A zero-carbon cell scores as maximum risk.** Measured on the processed ssp585 file: the Sahara reads 0.00 kg m⁻² at the **96th** risk percentile, the Taklamakan 0.04 at the 93rd, Greenland's ice sheet 0.00 at the **97th** — because `higher_is_better` maps a low stock to a high risk score, and a place that never held soil carbon holds the least of it. Arithmetically correct, and the conclusion inverts: a desert has nothing to lose. This is the RELATIVE-BASELINE class of trap and is recorded in the registry `delivery_note`, DATASET-ATTRIBUTES.md and the taxonomy blocker. Greenland is in the footprint at all only because two of four models write 0 over the ice sheet rather than NaN.

**7. `--measure-slopes` confirmed `sen_slope`, but the agreement figure is the story.** 70,910 active cells, sen exactly zero on only 6.8%, sign agreement **70.3%** — the lowest of any `sen_slope` layer (`conifer-npp` is 87.9%). The two estimators disagree in sign on nearly a third of active cells, and on the global mean throughout (ssp585 2090s: ols −0.0475 against sen −0.0072 kg m⁻² dec⁻¹). This is the ols-is-biased regime, not the sen-collapses regime — the only shipped layer in it.

**8. QA/QC surfaced a claim I had propagated from the old build, into a customer-facing field.** The processor docstring said JULES's fixed-CO₂ trend was **"muted"** — the deduction from the mechanism's name, carried forward unmeasured, and repeated in the `interpretation_caveat` that delivery exports. The `/isimip-process-visualize` skill already recorded the opposite for *this exact layer*. Measured on ssp585, 2090s−2020s over each member's own footprint:

| model | CO₂ | 2020s mean | change | % |
|---|---|---|---|---|
| `lpjml5-7-10-fire` | transient | 28.07 | −0.772 | −2.75% |
| `jules-es-vn6p3` | **fixed 2015** | 10.74 | −0.470 | **−4.37%** |
| `mc2-usfs-r87g5c1` | transient | 8.33 | −0.004 | −0.05% |
| `classic` | transient | 6.70 | +0.053 | **+0.79%** |

The fixed-CO₂ member is the **largest relative loser**, and one transient model *gains*. The mechanism is real — no fertilization means less litter input — but it makes the decline **stronger**, not weaker. GUARDRAILS §9 had asserted the muted version as a general rule since 2026-07-25; it is corrected there, in the catalog, the download script, the processor and `interpretation_caveat` (which required a full reprocess, since the text lives in the file). **A mechanism tells you a term exists; only a measurement tells you its sign.**

**9. The contact sheet earned its place — twice.** No `csoil_members.nc` existed, so `generate_maps.py` silently skipped the Members tab; the processor had never emitted one. Added, with a `--members-only` fast path that rebuilds it in ~2 min instead of re-running the slope stage. Rendered to PNG and **looked at**, which found two things no table could: `classic` is visibly **blocky** — confirmed at **100% adjacent-cell ties on odd column pairs against 0% on even**, a 1° field replicated 2×2 with a one-cell longitude offset, exactly as the skill's trap table warned — and `mc2-usfs` has a **truncated northern footprint** (58,919 cells against jules's 67,647; 12,417 cells ≥60°N against ~17,200), so ensemble support above 70°N falls to **12.2 of 17 members**, weighted toward the two models that read highest there. The high-latitude caveat is therefore about thin *and* biased support, not just disagreement.

**Impact**: no customer artifact was affected. The layer is registered `status: alternate`, no asset type routes to it, and the `soil-degradation` family's `covered_by` is deliberately empty, so nothing counts it as hazard coverage. Two questions are open for the user: whether this files as a hazard or as asset condition (it is `higher_is_better` with an inverted percentile, structurally identical to `conifer-npp`), and whether to annualise ELM-ECA and VISIT's monthly runs, which would take the ensemble from 17 to 27 members under GUARDRAILS §1–2.

**What it is not**: a soil-degradation layer. Soil organic carbon is one of the **ten** degradation processes in recital 4 of Directive (EU) 2025/2360, and no framework defines soil degradation as a climate hazard at all — ESRS E1 AR 11 names it and defines nothing, citing a delegated regulation whose Appendix A *is* that table. Land use is held fixed in every member, so the layer cannot see management-driven loss, which is most of observed loss.

**Rules reinforced**: a `processed:` date is not proof of an artifact — check the disk. A sector-completeness claim is not a variable-completeness claim. A caveat is delivered only if its attribute is on `LAYER_ATTRS_EXPORTED`. A per-worker memory constant is a per-worker memory *leak* once a pool is introduced. And a mechanism's name is not a measurement of its direction.

**Files**: `scripts/download_csoil_isimip3b.py`, `scripts/check_csoil_nature.py` (both new), `scripts/process_csoil_soilcarbon.py`, `scripts/utils/delivery.py` (missing `NamedTuple` import), `config/layer_registry.yaml`, `config/hazard_taxonomy.yaml`, `config/isimip_search_catalog.yaml`, `DATASET-ATTRIBUTES.md`.

---

### 2026-08-16: The threshold ladder — a censoring result that was a continent, an invented caveat name, and a `git checkout` I should not have run

**What happened**: nine threshold-exceedance layers were built from ISIMIP3b daily `tasmax`/`tasmin` (12 GCMs × 3 SSPs, ~1.34 TB streamed and deleted). Six things went wrong, all of them mine.

**1. I reported a censoring result that was Antarctica, not the product.** I measured the cold rungs as heavily censored at the calendar ceiling — `FD` 30.25% of pooled observations at ≥364 days, `ID` 28.15%, `FDm10` 21.09% — and told the user this was the `heatwave-3b` regime requiring a declared `saturation_caveat`. It was measured on the **full** ISIMIP3b land mask. On `landseamask_no-ant.nc` the same numbers are **2.01%, 1.08% and 0.00%**, and Antarctica alone accounts for 98.85% / 94.04% / 72.83%. The censored rungs are in fact the **hot** ones (`hd30` 11.8%, `TR20` 11.4% by ssp585 2090s). A mask is not a neutral denominator: 27,092 of 92,889 cells (29%) were a permanently-frozen block that both saturated the estimator and filled the top of the frost percentile with land nobody sites an asset on.

**Rule**: state the mask with every share-of-cells statistic. A percentage without its denominator named is not a measurement, and "share of land" is a different claim under two masks that both call themselves the ISIMIP3b landseamask.

**2. I read a masked array's `.data` and got the whole globe as land.** `landseamask_no-ant.nc` encodes land as `1` and everything else as `_FillValue = 1e20`, unlike the full mask which encodes `0`/`1`. `np.asarray(ds.variables["mask"][:]) > 0.5` therefore returned **259,200 "land" cells** — every cell on the grid — because the fill value is greater than 0.5. Caught only because 259,200 is obviously wrong; a subtler fill value would have passed. Fix: `.filled(np.nan)` before any comparison, never `np.asarray` around a masked array.

**3. I invented a caveat attribute name, which would have been dropped silently at delivery.** I wrote the GCM non-independence finding to `ensemble_independence_caveat`. `LAYER_ATTRS_EXPORTED` in `scripts/utils/delivery.py` is a **closed allowlist**, so it would have vanished between the file and the filing — the identical failure recorded on 2026-07-25 for `co2_treatment`. Diffing my attributes against the allowlist also found `source_dataset` expected and unwritten (no provenance at all would have reached the customer) and two applicable caveats missing (`sparsity_caveat`, `resolution_caveat`). **The incident log already contained this lesson and I reproduced it anyway** — so the check is now mechanical: assert the written attribute set against `LAYER_ATTRS_EXPORTED` before shipping, rather than trusting that a caveat written is a caveat delivered.

**4. A per-scenario decision on a per-layer field.** `recommended_slope` was computed and written as each scenario finished. `hd35` measures `sen_slope == 0` on 51–53% of active cells — straddling the 0.5 cut — so ssp126 would have recorded `ols_slope` and ssp585 `sen_slope` **for one layer**, while `layer_registry.yaml` carries exactly one. Three files would have disagreed about how to read themselves. Fixed by computing every scenario before writing any file; `emit_tasthresh_registry.py` now asserts agreement across a rung's files rather than trusting it.

**5. A caveat that said "none" would have been published as MUST-DISCLOSE.** Having added `saturation_caveat` and `sparsity_caveat`, I wrote them **unconditionally** — `"none -- the highest share at the ceiling is 0.69%"` on rungs that do not saturate. But `generate_delivery_caveats.layer_caveats()` promotes these to MUST-DISCLOSE on the test `if note and note.lower() != "nan"` — on the attribute being **non-empty**, not on it asserting anything. Seven of nine rungs would have carried a must-disclose caveat whose body reads "none", in both reports. That is precisely the defect `load_registry()` already refuses for a blank `relative_baseline_note`: a heading promising a caveat with nothing under it. The convention, visible in the shipped layers, is to **omit the attribute entirely when it does not apply** — `heatwave-3b` carries `saturation_caveat` because it saturates; `cropfailure-3b` carries none of the three. Fixed at write time, and the 27 built files patched in place (45 attributes renamed to `saturation_measured` / `sparsity_measured`, which are **not** on the allowlist, so the negative measurement stays auditable in the file without reaching a report). `resolution_caveat` remains unconditional, correctly — every rung is 0.5° support against a hazard that turns on elevation.

**Rule**: when adding a caveat attribute, read the *promotion condition*, not just the allowlist. Being on `LAYER_ATTRS_EXPORTED` decides whether a caveat can be delivered; the promotion test decides whether it *will* be — and a mechanism that fires on non-emptiness will happily publish your negative result as a warning.

**6. I ran `git checkout` on a file with uncommitted changes.** After destroying 92 comment lines by round-tripping `config/hazard_taxonomy.yaml` through `yaml.safe_load`/`yaml.dump`, I reverted with `git checkout` — a file that `git status` had shown as modified at session start. It happened to be safe (commit `07c63cc` had landed mid-session and already contained that work), but that was **luck, not diligence**: I discarded the working tree before establishing what was in it.

**Rules**: never round-trip a hand-maintained YAML through `yaml.dump` — it strips every comment; use targeted text edits. And never `git checkout` a modified file without first capturing it (`git diff > /tmp/...` or a copy), because the change you are discarding may not be the one you made.

**What went right, and is worth copying**: the **containment check**. Each higher threshold must be a strict subset of the lower one, so `hd30 ≥ hd35 ≥ hd40 ≥ hd45`, `TR20 ≥ TR25`, `FD ≥ FDm10`, and across variables `FD ≥ ID`. That is a *dimensional* check — it cannot be satisfied by a plausible-looking wrong answer — and it returned **0 violations** across all three scenarios, which no per-layer contract check could have established. It is the same class of check that caught the flood layers' "additive decomposition that was actually containment".

Also: the calendar risk the stage-1 reducer was designed around (a 360-day member gets six fewer chances a year, a ~1.5% bias that is constant per member and never averages out) **did not materialise** — all 32 members measured `proleptic_gregorian`. Measuring it was still right; the reducer's `*_days_per_year_HETEROGENEOUS` flag is an artifact of dividing each chunk's day count by its own unique-year count and must not be read as a mixed calendar.

**Files**: `scripts/download_reduce_tasthresh_isimip3b.py`, `scripts/process_tasthresh.py`, `scripts/emit_tasthresh_registry.py` (new), `config/layer_registry.yaml`, `config/hazard_taxonomy.yaml`, `config/isimip_search_catalog.yaml`, `DATASET-ATTRIBUTES.md`.

### 2026-08-18: Nature 2026 Hail Deposit -- Evidence Gate Before Any Layer Decision

**What this is**: not an incident. A dated record of the Phase 1 evidence gate run against
the deposited output of *Rising global hail damage potential in a warming world* (Nature
653, 1069-1077, 2026), before deciding whether it can support a hail layer. The gate was
designed after an adversarial review flagged four ways the obvious ingest would be wrong.
Script: `scripts/measure_hail_nature2026.py` (its docstring is the reference; this entry is
the outcome). Provenance: `data/raw/hail-nature2026/manifest.json`, 11 files, 6.81 GB,
CC-BY-4.0, sha256 per file.

**What the deposit is**: a semi-3D hail trajectory model run on 12,412 satellite-detected
severe hailstorms (2014-04-01 to 2021-03-31), historical from ERA5 soundings, future by
pseudo-global-warming on the SAME events. 20 EC-Earth3 realizations per SSP, plus
ensemble-mean runs for MPI-ESM1-2-LR and NorESM2-LM. No time axis, no frequency.

**Four traps, each of which produces a plausible wrong answer:**

1. **The row-position join is wrong and looks fine.** `hailstone_growth_radii_*` renumber
   events `0..n-1` and carry no coordinates. Joining historical row *i* to future row *i*
   gives r = -0.0002, i.e. noise. The real index is `date{type}` in the sounding profile
   files (copied into the `sample` coordinate of the growth-duration files), which is the
   row index into the catalogue's `para33_Idx{type}`. The authors join exactly this way
   (`np.isin(proff['date'], profh['date'])`). Verified four ways: per-type counts match the
   radii arrays exactly; paired Spearman 0.248 against a permutation null of 0.001; hail
   diameter correlates with the CATALOGUE's MUCAPE for the mapped event (rho 0.36 historical,
   0.45 future); and the authors' own code performs the same join.

2. **The historical and future files do not share support.** The publication seeds 330
   embryos per event (6 `w_perc` launch positions x 5 radii x 11 heights) and
   `severity_calc_n` hard-codes `ndia = 6*5*11`. The deposited FUTURE radii carry all 330;
   the deposited HISTORICAL radii carry 55 -- `w_perc = 0.5` only -- although the historical
   growth-duration files carry all six positions, so the full run exists and was not
   archived. Reducing over all six future positions against one historical position moves
   the median diameter shift from **+30% to +97%**. Every comparison must be made on the
   common `w_perc = 0.5` support, and the script refuses to run if that is violated.

3. **Three different denominators are all called "the >=30 mm change".** On ssp245 member 1,
   paired over 10,869 events: share of ALL seeded embryos landing >=30 mm **+46.0%**; share
   among embryos that produced a stone at all **+21.3%**; share of EVENTS whose largest stone
   reaches 30 mm **+19.5%**. The publication's headline is the first. A number quoted without
   its denominator is unusable here.

4. **The estimand decides the SIGN of the regional pattern.** Same events, same member
   (ssp245 m1), three statistics, 30-50N: median event diameter **-1.1%**, mean kinetic
   energy **+27.9%**, share of embryos >=30 mm **+25.1%**. US Plains: -0.2% / +32.0% /
   +26.6%. The median storm barely changes while the energy-weighted and large-stone measures
   rise ~30%, because kinetic energy goes as D^3.5 and the change is concentrated in the
   upper tail. Reporting median diameter would have contradicted the paper's own regional
   conclusion using its own data.

**The finding a customer would care about most, measured over all 60 member-runs.** The
paired median change in an event's largest stone is **+1.84 mm at ssp245, +0.52 mm at
ssp370 and -0.03 mm at ssp585**, with the share of events getting larger falling 55.0% ->
51.3% -> **49.3%** -- while over the same runs the share of embryos >=30 mm rises +51.7% ->
+54.3% -> +53.7% and mean kinetic energy +46.8% -> +50.2% -> +49.1%. The typical storm does
not change under the strongest forcing; the damaging tail grows by half. The two orderings
are not merely different in size, they run in **opposite directions across scenarios**, so
"does hail get worse here" has two defensible answers from one dataset and the answer must
name its statistic. Eligibility transitions grow with forcing too: 779 -> 1024 -> 1192
events change simulation eligibility between the historical and future climates.

**What reproduces.** Ensemble-mean runs, unpaired, common support, across 3 GCMs x 3 SSPs:
share of embryos >=30 mm **+37.5% to +52.8%** against a published **+37.9% to +51.8%** -- the
published range spans GCM x scenario, not scenario alone (MPI-ESM1-2-LR ssp245 is the floor,
EC-Earth3 ssp370 the ceiling). The deposited validation correlations reproduce **exactly**
(US r = 0.658 vs 0.66 published, China r = 0.633 vs 0.63) but only with the authors'
undocumented rule -- min-max normalise the month-then-year mean, then fill each field's NaN
with 0 where the other field is finite. Without the zero-fill: 0.573 and 0.558.

**What does not reproduce.** `<30 mm` change measures **-9.4% to -13.3%** against a published
-4.2% to -12.3%, and mean kinetic energy **+28.4% to +49.1%** against a published +36.5% to
+42.1%. Both run outside the published ranges in the direction expected if the five missing
launch positions matter. Until they are obtained, **no claim of reproducing the publication**
-- the drafted author query is `reports/maps/hail-severity/author_query.md`.

**Two facts that constrain any product built from this.** Event eligibility is
*climate-dependent*: profiles are retained only where the maximum updraft is >= 5 m/s
(`Calculation_of_Verticl_profiles_for_driving_hailstone_trajectory_model.py`, line 576),
so 11,158 events simulate historically, ~11,282 under ssp245 m1, and 10,869 pair -- events
enter and leave the population because of the treatment. And the ensemble is 20 realizations
of ONE model: across members the >=30 mm change spans +46.0 to +56.8% at ssp245, while
across the three GCMs the kinetic-energy change spans +28.4 to +49.1%. **Structural spread
is larger than internal spread**, so a member-only interval would understate uncertainty by
more than it captures.

**Status**: gate PASSED on the join, the estimand definitions, the >=30 mm headline and the
validation correlations; OPEN on the 55-vs-330 support. No layer built, no taxonomy change,
nothing written to `data/processed`. The deposit gives conditional SEVERITY only -- it cannot
express hail frequency, occurrence or return period, and at 12,412 events a regular global
grid averages 4.8 events per cell at 5 degrees.

**Files**: `scripts/measure_hail_nature2026.py`, `reports/maps/hail-severity/` (gate.csv, author_query.md),
`data/raw/hail-nature2026/manifest.json`, `isimip-pipeline/pyproject.toml` (py7zr -- the
deposit ships nested .7z and this machine has no 7z binary).

### 2026-08-18: A Spatial Average Has a WEIGHTING, and Choosing It Wrong Flipped the Sign in Six Regions

**What happened**: building the hail severity fields from the Nature 2026 deposit
(`scripts/build_hail_severity_fields.py`), each cell's mean stone size and mean kinetic
energy were computed as the mean over STORMS of each storm's own mean -- one vote per storm.
The request was for the mean size of stones at a location, which is the mean over STONES.
The first regional table went out on the wrong one.

**Impact**: the headline conclusion inverted. Under storm-weighting, mean kinetic energy per
stone read US Plains -1.0%, Argentina -12.3%, Southern Africa -1.4%, and the summary written
from it said impacts get softer almost everywhere. Pooled over stones the same events give
US Plains **+18.8%**, Argentina **+6.2%**, Southern Africa **+16.8%**. Six of eleven regions
changed sign. Both numbers are correct arithmetic on the same data.

**Root cause**: warming lets more WEAK storms produce hail at all. Those newly-producing
storms enter the per-storm average carrying small stones and drag it down, while the stones
themselves are not smaller -- pooled over landed stones the US Plains ssp585 median diameter
rises 37.3 -> 42.6 mm over the very same events. Storm-weighting answers "what does a typical
storm here drop"; stone-weighting answers "how big are the stones that fall here". The
composition shift between them is the physical signal, so the two must diverge.

**Fix applied**: both weightings are published and named in the file --
`mean_diameter`/`mean_ke`/`production`/`ake` pooled over stones and embryos (primary),
`storm_mean_diameter`/`storm_mean_ke` storm-weighted (secondary) -- with the divergence
documented in the module docstring and both rendered in the QA maps so it cannot be
overlooked again.

**A second thing the fix bought**: with pooled fields every quantity is a ratio of the same
summed numerators and denominators, so `ake == production * mean_ke` holds EXACTLY at cell
level (max residual 2.7e-16, asserted in the build). Under storm-averaging that identity
broke by 15-23% because a cell mean of a product is not the product of the cell means when
the factors covary -- a residual an earlier version of the check mislabelled as a wiring
warning. An identity that is exact by construction is worth engineering FOR, not just
asserting: it turned a diagnostic that could only ever be argued about into one that either
holds or reveals a broken denominator.

**Rule**: state the weighting with every spatial average of a per-event quantity, and where a
composition shift is part of the signal, publish both. A mean whose weighting is unstated is
not a measurement -- the same defect as a share of cells whose mask is unnamed.

---

### 2026-08-18: The precipitation layers — a mask that published as a real zero, and a slope I chose for consistency instead of measurement

**What happened**: ten layers were built from ISIMIP3b daily `pr` (14 GCMs x 3 SSPs, ~911 GB streamed and deleted). Two defects, both caught before delivery, and both instructive about *what kind of check finds what kind of error*.

**1. 599 masked cells published as `0`, not as missing — and every contract check passed on it.** The relative metrics (`R95pD`/`R99pD`) count days above each cell's own 2020s wet-day percentile, and 599 hyper-arid cells have too few wet days for a percentile to be definable. Stage 1 computes `hit = (value >= threshold) & usable`, so an unusable cell accumulates a count of **zero** — a perfectly legitimate integer. The `-1` sentinel I had built only ever marked *ocean*, where `scatter` left NaN. So those cells shipped as "0 days above the local 95th percentile": indistinguishable from a genuine zero, ranked at **percentile tier 1 — the LOWEST risk band** — in the driest interiors on Earth, when the honest answer is "no threshold could be defined here".

`test_shared_baseline.py` passes on this without complaint. `lower_ci <= median <= upper_ci` holds at 0. The percentile is well-defined. There is no ocean leak. **Nothing in the contract can distinguish a real zero from a fabricated one**, because the defect is semantic, not algebraic.

What found it was checking a *specific expected property of the one code path that was new*: "the 599 cells the mask excluded should be NaN — are they?" Answer: 0 of 599. Fixed by applying the `usable` mask on load in stage 2, plus an **assertion before writing** so a recurrence fails the build. The lesson generalises past this layer: **a sentinel that is only applied on one axis of a two-stage pipeline is not a sentinel.** Stage 1 knew the cell was unusable and encoded that as a value indistinguishable from data.

**2. I chose a slope estimator for family consistency and the measurement overruled it.** I set `recommended_slope: ols_slope` on all ten, reasoning that ols's failure mode requires uneven member coverage and coverage is 14 members everywhere — so ols cannot be biased here. That argument is correct and **establishes only that ols is SAFE, not that it is BEST**. Measured on active cells at the final decade, Theil-Sen had not collapsed at all on the three mm-valued metrics (`sen==0` on 0.0%-1.4%, sign agreement 0.92-0.95) — cleaner than either existing `sen_slope` layer in the product (`conifer-npp` 2.1%, `csoil` 6.8%) — while collapsing on 41.8%-100% of the seven day counts. `Rx1day`/`Rx5day` have a structural claim on the robust estimator besides: an annual maximum is a heavy-tailed extreme-value series, exactly where one freak year drags a least-squares fit and a median of pairwise slopes does not. Corrected to **counts -> mean + ols, mm -> median + sen**, so the slope now splits on the same boundary as the statistic and both follow measurement.

**Rule reinforced**: "which estimator can fail here" is the right question, but it has two halves. Ruling out one estimator's failure mode does not select it — check whether the other estimator actually failed before defaulting to consistency.

**What went right and is worth copying**: the containment check. `wetdays >= R10mm >= R20mm >= R50mm >= R100mm`, `R95pD >= R99pD`, and **`Rx5day >= Rx1day`** — 0 violations across all three scenarios. That last relation is the only independent test on real data of stage 1's chunk-boundary carry for the 5-day running window; a synthetic test showed that dropping the carry understates **42.5% of cells** in every year after a boundary, and no per-layer contract check would have seen it, because the values would simply be slightly too low everywhere, forever.

**Files**: `scripts/pr_baseline_percentiles.py`, `scripts/download_reduce_prthresh_isimip3b.py`, `scripts/process_prthresh.py`, `scripts/run_pr_pipeline.py`, `scripts/emit_prthresh_registry.py`, `scripts/qa_prthresh_dashboards.sh`, `scripts/check_pr_nature.py` (all new), `config/layer_registry.yaml`, `config/hazard_taxonomy.yaml`, `config/isimip_search_catalog.yaml`, `DATASET-ATTRIBUTES.md`, `scripts/README.md`.

---

### 2026-08-18: Water-stress build — five defects, three external review rounds

Layer: `waterstress-3b-*`. Record: `docs/water-stress-status-2026-08-17.md`.
Generalizable rules extracted to **GUARDRAILS §18**.

**1. Area convention settled twice, in opposite directions.** Concluded `cellarea` alone on
0.6% global agreement; reversed to `cellarea × contfrac/100` after an external review
challenged the comparison's provenance. The reference was a three-model average from a
previous simulation round; the model ran 26% low and the missing factor 27% high, and they
cancelled. Resolved by a stratified residual (structured signature) plus a global total that
only one candidate satisfies across all 15 members. → §18.1

**2. Calendar asserted from a sidecar, wrong for the model in use.** Wrote
`proleptic_gregorian`; all 45 WaterGAP files declare `365_day`. H08 genuinely declares
`proleptic_gregorian`, so it is heterogeneous per model. → §18.3

**3. `np.asarray()` dropped a mask; 1e20 fills passed `isfinite`.** Global sum returned
3.6e41. Initially misdiagnosed as a `netCDF4` defect; the library was innocent. → §18.2

**4. River mask built from the ensemble mean.** Let one member's water license another's
ratio — 686 cells passed a pooled mask while their own baseline was below it; one returned
6e9 in every year. → §18.4

**5. Three masks proposed against extreme ratios, each on a wrong hypothesis.** "Tiny
rivers" (the offending cells had a median baseline of 5.0 m³/s), then "year-specific
collapse" (excluding collapses left the max untouched). Only inspecting the single worst
cell-year found the real cause. A proposed 50 m³/s mask was rejected by the user on the
grounds that it is the Thames at London and would delete exactly the stressed dryland
basins — correct, and it had been proposed from a max value without checking how many cells
produced it (0.016%, 174 cells). → §18.6

**Process note.** Three Codex review rounds. Round 2 raised 4 blockers, round 3 confirmed the
reversal but found the precedence table still ordered the river mask after the zero rules,
and caught a mislabelled column in the evidence table (ratios of totals presented as
medians). **Every round found something the previous one had not**, and two rounds found
errors in the corrections themselves. Budget for more than one review pass on a layer whose
arithmetic is novel.

---

### 2026-08-19: Landslide — the obvious aggregation was degenerate, and the licence is the real blocker

**Context.** Asked to document the landslide options and aggregate landslide risk to 0.25° on
the standard output shape. Full options review:
[docs/landslide-data-options-2026-08-19.md](docs/landslide-data-options-2026-08-19.md).

**1. The ISIMIP negative is real but weaker than tornado's, and saying so matters.** ISIMIP
publishes no landslide output in any round *and* `InputData/geo_conditions/` carries no
slope, elevation or lithology. But landslide hazard factorises into terrain × trigger, and
the trigger IS in ISIMIP — NGI's GIRI global model is driven by ISIMIP3b precipitation
(IPSL-CM6A-LR, SSP126/585, 2061–2100). Recording this as "ISIMIP cannot express landslide",
by analogy with tornado, would have been wrong and would have foreclosed the one build that
would actually be ours. The catalog entry says so explicitly.

**2. Measuring the statistic before writing it saved a degenerate layer — again.** Each 0.25°
cell holds 900 native 1 km pixels, so within-cell quartiles look obvious. Over occupied cells
(n=89,871) they give `median == 0` in **78.48%** and `q25 == q75` in **64.73%** — three of
four published variables carrying no information over two-thirds of the map, in a file that
would still pass a structural contract check. This is the third time this failure has been
caught by measuring first (`let` 2026-08-11, tornado 2026-08-18, here). **The check is cheap
and it keeps paying.**

**3. Two further branches were tried and both failed, on numbers worth keeping.** The repo's
own `pooled_mean_zero_inflated` (mean ± SD) is non-degenerate but puts `lower_ci` below zero
in **88.55%** of occupied cells — unphysical for a frequency. Sub-block (0.05°) quartiles were
built specifically to give the areal mean a non-negative bracketing interval, and **cannot**:
the within-cell field is right-skewed, so the mean exceeds its own sub-block upper quartile in
**42.48%** of cells. That is a property of the distribution, not the block size — recorded so
nobody re-attempts it with a different block.

**4. A new departure: `percentile` ranks on a different variable than `median`.** The adopted
central value is conditional on hazard-bearing ground, which discards *extent*; the areal mean
integrates extent and intensity. Spearman between the two orderings is **0.34**, so the choice
changes essentially every customer score, and reference sites settle it (Apennines 74.2 vs
58.9; Cairo 1.4 vs 5.3 off one pixel covering 0.1% of the cell). Declared in the file, the
registry and DATASET-ATTRIBUTES, because a cell showing `median = 0` at a high percentile
looks like a bug and is not.

**5. Zero was ambiguous in the source.** The COG declares no nodata and writes exact `0.0` for
ocean *and* flat land — Pacific, Atlantic, Sahara, Amazon, Netherlands, Greenland all
`0.000000`. Masked as source extent ∩ (ISIMIP3b land ∪ any hazard-bearing cell), the union
deliberately so a coarse coastline cannot erase ground the source modelled as hazardous.

**6. The blocker is legal, not scientific, and it was nearly missed.** The one global product
with a real rate has **two publisher-side records that disagree about its licence** — the World
Bank DDH says CC BY-NC 4.0, the energydata.info mirror of the same dataset says CC-BY-4.0, and
the 113-page project report states neither. NC would forbid redistribution in a commercial
deliverable outright. The layer is `status: blocked`. **Recorded as `source_licence_status`, a
non-caveat attribute, not as a caveat** — caveat attributes are promoted verbatim into customer
reports and our licensing problem is not the customer's.

**Process note.** `rasterio` (and `pypdf`) were installed into `.venv` to do this work and are
**not declared in `isimip-pipeline/pyproject.toml`**, which had uncommitted edits from the user
at the time and was deliberately not touched. Anyone rebuilding the venv will need
`pip install rasterio` before `process_landslide_arup.py` runs.

**Follow-up, same day.** The licence was cleared by the user for our limited commercial use,
so the layer moved from `blocked:` to `layers:` with `status: preferred` and
`covered_by: [landslide-arup]` in the taxonomy. **Attribution to World Bank / GFDRR and Arup
is now a shipping requirement**, carried in the file's `attribution_required` attribute, the
registry `delivery_note` and the ingest sidecars. The publisher's own inconsistency is
retained in `source_licence` rather than deleted — it is a real defect in the source's records
and the next person to check provenance will hit it.

QA maps were built (`scripts/generate_landslide_qa.py` → `reports/maps/landslide/`), which is a
third renderer for the observational contract after `generate_maps.py` (decadal only) and
`generate_tornado_qa.py` (CONUS rungs). **`generate_layer_qa.py` and `test_shared_baseline.py`
both crash with `KeyError: 'decade'` on any observational layer** — verified identical on the
shipped tornado layer, so pre-existing and not introduced here, but it means no observational
layer has ever had a markdown QA record. Left alone; `test_shared_baseline.py` guards 23 layers
and is not mine to widen unasked.

**A wrong comment, caught by measuring the thing it asserted.** The QA page came out at 50 MB.
I attributed that to unrounded floats serialised as JSON text, added per-metric rounding, and
wrote a comment stating that decimal places were the dominant payload term. The page went to
39.8 MB — entirely from an unrelated latitude crop in the same edit. Checking the file showed
plotly 6.9 emits **base64 typed arrays** (`bdata`), so rounding changed nothing at all and the
string `"null"` appears once in the whole document. The real term is bytes-per-cell: casting
the panels to float32 took it to 21.7 MB. **The lesson is not about plotly** — it is that a
plausible explanation written into a comment as if measured is indistinguishable from a
measured one to the next reader, and this repo's own rule (a claim needs its receipt) applies
to comments about performance exactly as it does to claims about data. The corrected comment
now records what was measured, including that rounding is retained only as a guard.

**Sign-off and cleanup, 2026-08-19.** User reviewed `reports/maps/landslide/landslide-qa.html`
and signed off on the layer as an **initial high-level landslide risk mapping product**.
`qa_reviewed_on: "2026-08-19"` is set; the scope of that sign-off (screening only, still no
scenario axis) is a comment beside it. `landslide-arup` is now the first layer to go from
first ingest to signed-off delivery-ready inside one session.

Three things came out of the cleanup and all three generalise:

**1. QA/QC maps belong in `reports/maps/{hazard}/`, and a one-off renderer is exactly where
that slips.** `generate_maps.py` has always written there, but both narrow per-contract
renderers — `generate_tornado_qa.py` (2026-08-18) and `generate_landslide_qa.py` (2026-08-19)
— invented their own top-level folders, `reports/tornado-qa/` and `reports/landslide-qa/`.
Both were moved and every reference updated. The rule is now in CLAUDE.md and in the
`/isimip-process-visualize` skill under *Output Organization*, phrased as **set `OUT_DIR` when
you create the script, not after**. Note what does *not* move: `reports/qa/` (markdown
records), `reports/{var}_model_diagnostics/`, and licence/author correspondence directories are
not maps. `reports/coastal-inundation/maps/` is still non-conforming and was left alone — it
has no references anywhere in the tree, so moving it would be a guess about what it is.

**2. `config/layer_registry.yaml` has a schema, and an extra key breaks EVERY layer.**
Adding `qa_reviewed_scope:` alongside `qa_reviewed_on:` raised
`TypeError: __init__() got an unexpected keyword argument` out of `LayerSpec`, and because
`load_registry()` builds all layers in one dict comprehension, one unknown key on one layer
takes down the whole registry — not just that entry. Caught immediately because the verify step
ran, but the blast radius is the point: the registry is not a free-form notes file. **Put
anything that is not a `LayerSpec` field in a YAML comment**, which is what the landslide
sign-off scope is now.

**3. The browser-payload guidance in the skill was right for the wrong renderer.** Its
`COORD_DECIMALS` / `VALUE_SIGFIGS` knobs are real, and they work because a `Scattergeo` is
serialised as JSON text. They do nothing on a raster `go.Heatmap`, where plotly ≥ 6 emits
base64 typed arrays and the term is bytes-per-cell. Both facts now sit in the same section of
the skill so the next person does not apply the text-encoding fix to a binary-encoded page, as
I did. See the same-day entry above for how that error was made and caught.

### 2026-08-20: Tornado QA pages were never human-reviewable

**What happened**: The tornado QA renderer was orphaned 24 minutes after its last run by the `historical` → `observed` scenario-token rename — it printed "no layers found" instead of raising, so nothing flagged it. Once revived at the QA review, its uncropped global float64 Heatmap panels emitted 153.7 MB pages no browser could usefully open. The rungs shipped 2026-08-18 with evidence pages nobody could actually have read.

**Impact**: A QA sign-off would have rested on unviewable evidence; caught only because the 2026-08-20 review regenerated evidence before presenting it.

**Prevention**: A renderer that discovers zero inputs must raise, not print; regenerate evidence at review time and compare page mtimes against data mtimes; the Heatmap payload rules (crop to the data's extent, cast float32) are in the process-visualize skill. Fixed in `aef11c0` and `39db4ae`; recorded in the tornado sign-off record in DATASET-ATTRIBUTES.md.

### 2026-08-21: Date rollover split a delivery across two dated folders

**What happened**: A `--run` re-extract after midnight wrote `deliveries/storefront-test/20260821/` while the verifier and caveat commands were still pointed at the previous day's `20260820/`, which held pre-decade-policy CSVs. The verifier "failures" were real disagreements — with a stale folder.

**Impact**: ~270 phantom verification failures and a debugging detour before the path mismatch was noticed.

**Prevention**: After any `--run`, point every downstream command at the output path the run **printed**, never at a remembered path; a delivery re-extracted near midnight lands in a new dated folder by design. Delete or clearly supersede the stale folder in the same session.

### 2026-08-21: Dashboard table rebuilt itself out from under its own dropdowns

**What happened**: The Values table rebuilt its entire DOM (header, filter row, ~7,000 body rows as one innerHTML string) on every dropdown change — destroying the `<select>` the user was interacting with while its change event was still dispatching, and generating enough DOM churn to crash the tab (which auto-reloads, reading as "the filter reset"). Two fix rounds treated symptoms before the mechanism was restructured.

**Impact**: Repeated user-visible filter failures across three feedback rounds.

**Prevention**: The stable-header pattern is now the rule (ASSET-CATALOG.md, dashboard section): controls are built once and never rebuilt — only `<tbody>` re-renders; handlers are delegated to a persistent ancestor; rendered rows are capped with the truncation stated; and a global error banner surfaces any uncaught page error as text a reviewer can report verbatim, so no control can fail silently again.

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

---
