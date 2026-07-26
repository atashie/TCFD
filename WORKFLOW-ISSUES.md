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

**Processor**: `scripts/process_csoil_soilcarbon.py`. Ensemble = 3 models × CMIP6 GCMs (2+5+5) × {ssp126, ssp370, ssp585} = 12 members/scenario, raw kg C m⁻², no normalization, no spatial smoothing, `higher_is_better`, baseline-anchored trend. See memory `project_csoil_soilcarbon_layer`, `reference_isimip_enumeration_curl`, `feedback_verify_experiment_tokens`.

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
