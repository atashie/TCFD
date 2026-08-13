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

## Adding New Incidents

When documenting a new incident, include:

1. **Date and title**: Brief description of what went wrong
2. **What happened**: Detailed narrative of the error
3. **Impact**: What was the consequence (data loss, incorrect results, wasted time)
4. **Root cause**: Technical explanation of why it happened (if applicable)
5. **Correct action**: What should have been done instead
6. **Fix applied**: How was it resolved (if applicable)
7. **Rule created**: Reference to the guardrail added/updated
