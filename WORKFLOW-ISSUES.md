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

## Adding New Incidents

When documenting a new incident, include:

1. **Date and title**: Brief description of what went wrong
2. **What happened**: Detailed narrative of the error
3. **Impact**: What was the consequence (data loss, incorrect results, wasted time)
4. **Root cause**: Technical explanation of why it happened (if applicable)
5. **Correct action**: What should have been done instead
6. **Fix applied**: How was it resolved (if applicable)
7. **Rule created**: Reference to the guardrail added/updated
