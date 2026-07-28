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

---

## 4. Multi-Model Search for Vegetation Variables

**Rule**: When searching for vegetation variables (cwood, npp, gpp, cveg, cleaf, lai, burntarea), **ALWAYS enumerate all available biomes models across simulation rounds** before presenting options to the user.

**Why this matters**:
- Different models use different PFT (Plant Functional Type) classification schemes
- Generic PFT schemes (e.g., CLASSIC's `evgndltr`) combine all needleleaf into one category
- Climate-zone specific PFTs (e.g., CLM45's `needleleaf-evergreen-tree-temperate`) exist in older rounds
- Best climate-specific data may be in ISIMIP2b or ISIMIP3a, not ISIMIP3b
- Users need to make informed trade-offs between data recency and climate specificity

**Required workflow**:
1. Identify target vegetation type and climate zone (temperate, boreal, tropical, etc.)
2. Check `config/isimip_search_catalog.yaml` biomes_models section for all available models
3. List models from ALL relevant simulation rounds (3b, 3a, 2b)
4. Check each model's PFT scheme against user needs
5. Present options with explicit trade-off summary:
   - Newer scenarios (ISIMIP3b SSP) vs climate specificity (ISIMIP2b/3a)
   - Generic PFTs (larger ensemble) vs climate-zone PFTs (better proxy)
   - Global coverage vs regional detail

**Do NOT**:
- Assume the first model found is representative of all available data
- Skip older simulation rounds without explaining trade-offs to user
- Present generic PFT (e.g., `evgndltr`) when climate-specific variants exist
- Default to ISIMIP3b without checking if ISIMIP2b/3a has better PFT coverage

**Example failure**: Searching for "loblolly pine timber proxy" and returning only CLASSIC `evgndltr`, which combines temperate AND boreal conifers into one PFT, when CLM45 has `needleleaf-evergreen-tree-temperate` (temperate-specific) and MC2-USFS has `mesictemperateneedleleafforest` (biome-specific).

**Reference**: See `/isimip-search-download` skill Step 0 (Model Enumeration) and catalog `pft_equivalences` table.

---

## 5. Diverging Colorscales Must Center on Zero

**Rule**: When visualizing trend or change data, diverging colorscales (RdBu, RdYlBu, etc.) **MUST use symmetric scaling centered on zero** so that white/neutral represents no change.

**Why this matters**:
- Diverging colorscales communicate directionality: one color = positive, other color = negative
- If zero is not at the midpoint (white), the visual interpretation is misleading
- A trend of +0.1 might appear red (bad) if the scale runs from -0.05 to +0.5 with white at +0.225

**Required behavior**:
- For trend and change maps: `max_abs = np.percentile(np.abs(all_values), 98); cmin, cmax = -max_abs, max_abs`
- Do NOT use percentile-based scaling like `cmin = np.percentile(all_values, 2), cmax = np.percentile(all_values, 98)`
- **Colorscale direction is keyed to the variable's `percentile_direction` attribute** (set by the processor):
  - `higher_is_better` (or unset — vegetation/productivity): default `RdBu`, positive (increase) = blue = good, negative = red = bad
  - `higher_is_worse` (hazards — drought, mortality): `generate_maps.py` reverses `RdBu` so positive (increase) = **red = worse**, negative = blue = better
  - Zero = white = no change in both cases

**Implementation**: `scripts/generate_maps.py` uses symmetric scaling for trend and change maps, and `create_map_figure(..., reversescale=...)` flips the diverging scale when the processed file declares `percentile_direction: higher_is_worse`. Variables without that attribute keep the default `RdBu`. See WORKFLOW-ISSUES.md 2026-07-24.

---

## 6. Water Index Is a Separate Workflow From Standard TCFD

**Rule**: The water risk index (`waterIndexUnderlyingData_*.nc`, 20 value types) is a **completely independent data product** from the standard TCFD annualized pipeline (6 value classes). **NEVER apply standard TCFD pipeline concepts** (kernel smoothing, Theil-Sen trends, percentile-of-score ranking, shared 2020s baseline) to the water index workflow.

**Why this matters**:
- The standard TCFD pipeline produces 6 value classes: smoothed median, percentile score, trend, significance, lower/upper CI bounds
- The water index produces 20 value types: 12 monthly ensemble means + annual mean + 7 annual quantile breakpoints (Q05-Q95)
- These are fundamentally different statistical approaches — confusing them produces incorrect output
- The R code that generated the original water index files was NOT found in `_deprecated/` — it was a separate codebase

**Required behavior**:
- Water index scripts: `process_water_tws.py`, `process_water_variable.py`, `config_water_*.py`, `validate_water_tws.py`, `compare_water_index.py`, `download_water_*.py`, `diagnose_*_models.py`
- Standard TCFD scripts: `process_qg.py`, `process_fish_*.py`, `process_health_*.py`, `generate_maps.py`, `isimip-pipeline` CLI
- Never import functions from one workflow into the other
- Never apply `/isimip-process-visualize` skill to water index processing

---

## 7. Water Index Value Types Must Follow the 20-Type Format

**Rule**: Water index output files MUST have exactly 20 value_types with this specific structure:
- **vt 0-11**: Per-month ensemble means (Jan-Dec) within each decade
- **vt 12**: Annual mean (= mean of vt 0-11, NOT a smoothed median)
- **vt 13-19**: Annual quantile breakpoints Q05, Q15, Q25, Q50, Q75, Q85, Q95

**Why this matters**:
- The water index requires seasonal cycle information (vt 0-11) for downstream risk assessment
- Quantile breakpoints (vt 13-19) capture ensemble spread in physical units, not derived statistics
- vt 12 is simply the mean of the 12 monthly means — no smoothing, no kernel, no special processing
- Quantiles are computed by pooling raw annual values across all ensemble members and years within the decade

**Required behavior**:
- Quantiles MUST be monotonically non-decreasing: Q05 <= Q15 <= Q25 <= Q50 <= Q75 <= Q85 <= Q95
- vt 12 MUST exactly equal `np.nanmean(vt0:vt12)` — validated by `validate_water_tws.py`
- Annual values for quantiles are **always aggregated using mean** (matching vt12 units), regardless of whether the variable is a stock or flux. Summing monthly rate values (kg/m²/s) produces values ~12× the rate, breaking units consistency with vt0-12. See WORKFLOW-ISSUES.md 2026-04-07 incident.
- No kernel smoothing, no Theil-Sen trends, no percentile-of-score ranking in this workflow

---

## 8. Never Guess Specifier Codes — Enumerate; Treat count=0 / count=1001 as Triggers

**Rule**: ISIMIP specifiers (variable codes, PFTs, crops, size classes, hazard-exposure codes) are a **controlled vocabulary**. **Never conclude a variable does or doesn't exist from a guessed code.** Enumerate the actual family, and treat both an empty and a maxed-out result count as a signal to verify further — never as a final answer.

**Why this matters**:
- A plausible-looking code can be wrong: tropical-cyclone exposure is `let`, not the mnemonic `letc`. Guessing `letc` returned `count=0`, which was misread as "doesn't exist," nearly missing a ready-to-process TCFD layer. See WORKFLOW-ISSUES.md 2026-07-24.
- `count=0` for something that plausibly exists is a **red flag**, not a conclusion.
- `count=1001` is the **API maximum = truncated**; the result set is incomplete and must NOT be generalized from.

**Required behavior**:
- **Enumerate families, don't guess members.** The Lange 2020 extreme-event exposure family is exactly six: `led` (drought), `leh` (heatwave), `lew` (wildfire), `ler` (river flood), `lec` (crop failure), `let` (tropical cyclone) — all ISIMIP2b, annual (data nature varies by member: `led` binary, `let` fractional — value-check each, never assume from a sibling). Authoritative table in `config/isimip_search_catalog.yaml` → `search_results.drought.exposure_lange2020.family`. Consult it before searching TC/heatwave/flood/wildfire/crop-failure exposure.
- **A variable can also have a non-exposure representation.** "Wildfire" is both the Lange 2020 `lew` *exposure* member AND the direct `burntarea` burnt-area-fraction fire output (ISIMIP2b `biomes`, in %); enumerate all products (see `search_results.wildfire`) and present the trade-offs rather than assuming the exposure family is the only option.
- On `count=0` for a plausible variable: fall back to file-server enumeration (`https://files.isimip.org/{round}/...`) or query each candidate code before concluding non-existence.
- On `count=1001`: narrow the query (by round/product/sector/model) until under the cap before drawing any conclusion about coverage.
- **When processing one member of a known family, record the whole family in the catalog** at that time (do not leave siblings un-enumerated).
- **The catalog's per-model variable lists are opportunistic, not exhaustive.** They were compiled for prior searches, so a variable can exist in models the catalog doesn't list under it (the soil-carbon `csoil` search 2026-07-25 found **11** ISIMIP2b + **5** ISIMIP3b biomes models, vs the 5 the catalog documented). Enumerate the **file server**; do not trust the catalog's `variables:` as complete.
- **The repository API (`data.isimip.org`) is behind Anubis anti-bot and blocks WebFetch/curl** (returns an "Access Denied" challenge page). Enumerate via the **file server** `https://files.isimip.org/{round}/OutputData/{sector}/...` instead. WebFetch's page-summarizer also **silently truncates long directory listings** (it repeatedly claimed "only one file" for a many-file `csoil-total` listing) — for an exact, complete file list use `curl -s <dir-url> | grep -oE '<pattern>'`, not WebFetch.

---

## 9. Verify a Variable's Data Nature Empirically Before Choosing Statistics

**Rule**: Before selecting the decadal statistic, aggregation, or percentile method, **inspect the actual raw values**. Do not infer a variable's data nature from its name, CF long_name, or a sibling variable.

**Why this matters**:
- The Lange 2020 exposure family shares a filename pattern and a "land area ... exposed to X" long_name, but the underlying values differ by member: `led` (drought) is a **binary {0,1}** per-cell flag, while `let` (tropical cyclone) is a **continuous fraction [0,1)** (~97% exact 0, remainder smooth over (0,1), never 1). Verified 2026-07-24 after a user challenge; assuming `let` inherited `led`'s binary nature would have mis-framed the product.
- Binary vs fractional vs continuous changes the correct decadal statistic (mean vs median), the CI definition, and the percentile treatment.
- **The metadata itself can be wrong or diverge across members of the same variable.** For `burntarea` (3 fire models): `lpj-guess` mislabels its `long_name` as "Fire Return Interval" though the values are burnt-area % (and floors at 0.1%, so it never emits a true zero); `mc2-usfs` encodes its time axis as **`days since 1661`** while the others use `years since 1661`, which silently breaks a naive year parser. Trust the values, not the labels; parse time axes defensively.

**Required behavior**:
- For any new variable, print min/max, the counts of exact-0 / exact-1 / interior values, and the number of unique values before processing. Two exact values = categorical/binary; a smooth interior distribution = continuous/fractional.
- **Also verify per-member metadata**: units, `long_name`, and time-axis (`units`/`calendar`) — do NOT assume they match across members. Parse time robustly (years/days/months per calendar) and add a guard that skips/flags a member with zero in-window years rather than corrupting indices.
- **Same unit ≠ needs normalization; different unit ≠ always poolable.** If ensemble members share a unit but differ in magnitude (e.g. `burntarea`: mc2-usfs ~5–7× hotter), that spread is genuine model uncertainty — pool in raw units (no normalization) and let it widen the CI. Reserve robust z-score normalization for members on genuinely different scales (e.g. water-index TWS). Decide per variable and record it.
- **Verify the experiment tokens (soc / CO₂-sens) per model — a uniform treatment across the ensemble is NOT guaranteed.** For `csoil` (2026-07-25), `jules-es-vn6p3` publishes ONLY its fixed-2015-CO₂ run (`2015co2`) while `classic`/`mc2-usfs` publish the transient `default` run, and soc tokens differ (`2015soc-from-histsoc` vs the natural-veg `nat`). Enumerate each member's sens/soc token *before* committing to a "transient CO₂" (or any) ensemble; if you must mix treatments to preserve ensemble depth, get the user's call and record it (`co2_treatment: mixed`).
- **Do not assume which DIRECTION a fixed-CO₂ member biases the trend — measure it.** This guardrail previously asserted that a fixed-CO₂ model's stock trend is "muted (no CO₂ fertilization)". For `csoil` that is **backwards**, as measured 2026-07-27: `jules-es-vn6p3` (fixed 2015 CO₂) has the **strongest and most scenario-sensitive** trend of the five models (global-mean 2090s change +1.57% / −2.92% / −4.37% for ssp126/370/585), while the four transient-CO₂ models run flat-to-positive (+0.8% classic, ~0% mc2-usfs, +2.6% elm-eca, +7.1% visit at ssp585). Removing fertilization does not damp the signal; it removes the *offsetting gain*, so warming-driven decomposition dominates and the loss grows with forcing. Report a fixed-CO₂ member's per-model trend explicitly rather than describing it as muted.
- **Never report spatial resolution from `ds.sizes` or coordinate spacing — measure it from the values.** A model that runs natively coarser and is replicated onto the ISIMIP grid still reports the full 360×720 at 0.5°. Dimensions describe the container, not the information content. For `csoil` (2026-07-28) `classic` is natively **1°×1°**, replicated 2×2 with a **one-cell longitude offset**, while its header is indistinguishable from the four genuine 0.5° members. Test the exact-tie fraction between adjacent cells **at both block offsets** (an origin-aligned test found only 3.3% and nearly cleared it), or the inside-vs-seam gradient ratio per candidate block width. A variance-loss-under-coarsening test does **not** discriminate blockiness from ordinary geophysical smoothness. `scripts/generate_qa_report.py` now checks this automatically. Also record which cells are carried by a single model, since they inherit that model's native resolution unchanged.
- Record the verified data nature in the catalog (`data_nature`) and in the output metadata.
- Never copy a sibling variable's processor attributes verbatim — re-derive the framing from the verified values.

**Related techniques** (documented in WORKFLOW-ISSUES.md 2026-07-24; patterns, not hard rules):
- **Thin ensembles** (few members, e.g. 1 impact model × 4 GCMs): apply spatial smoothing (5×5 exponential-decay, cos(lat)-scaled, normalized over non-NaN neighbours) to per-member decadal maps to borrow strength from neighbours. Reference: `scripts/process_let_cyclone.py`.
- **Zero-inflated hazards**: two-tier percentile (zeros → 1; non-zeros ranked against the non-zero 2020s baseline → [2,100]) so the exposed gradient isn't crushed into the top quantiles.

---

## 10. Trend Must Be a Decadal Signal, Not a Within-Decade Annual Slope for Noisy Variables

**Rule**: For a variable with high interannual variability (fire, precipitation extremes, floods), do **not** report the trend as the OLS slope of the *annual* series within a single decade. Use a **baseline-anchored across-decade rate** so the trend is spatially coherent and consistent with the change map.

**Why this matters**:
- The legacy processors (`process_qg.py`, `process_led_drought.py`, `process_let_cyclone.py`) compute `trend[decade]` as the slope of the annual values *inside* that decade (10 points). For a noisy variable this slope is dominated by interannual noise → a spatially **spotty, sign-flipping** field, while the change map (a difference of two 10-year decadal means) is smooth. A user flagged exactly this for `burntarea` (2026-07-24).
- A "trend" that contradicts the change map is misleading — they should tell the same story.

**Required behavior** (for high-variance variables):
- Compute the trend from the **decadal-median series**, anchored at the baseline decade:
  `trend[decade] = (median[decade] − median[2020s]) / (elapsed decades)`, units *value* · decade⁻¹.
- This makes each decade's trend the rate **from the 2020s baseline to that decade** (2090s panel = full baseline→2090s trend), built on decadal means so it is exactly the (decade − 2020s) change map ÷ elapsed decades → coherent by construction (corr = 1.0 with change at the last decade).
- The baseline decade has no elapsed change → trend 0 (keeps the 2020s baseline bit-identical across scenarios).
- The generate_maps trend label is `[units decade⁻¹]`, so emit per-**decade** units, not per-year.
- Reference: `scripts/process_burntarea_fire.py` (`anchored_trend`). Low-variance variables may keep the within-decade slope, but confirm the trend map is coherent against the change map before finalizing.


---

## 11. Look at the Data, Not Only Its Statistics

**Rule**: Before publishing a layer, **render the field and look at it — per member, not only the ensemble.** A layer is not verified until someone has viewed an image of it.

**Why this matters**:
- Every other verification rule here is **distribution-based** (§9: min/max, exact-0/1 counts, unique values, units, `long_name`, calendar). All of those are **invariant under spatial rearrangement** — shuffle a field's cells arbitrarily and every one still passes. They are structurally incapable of detecting a spatial defect.
- That blind spot shipped a layer: `csoil` (2026-07-28) passed a full tabular value-check twice, plus 37 algebraic QA checks, and was reported as verified — while one member (`elm-eca`) was effectively **~4°×5°** on a product sold as 0.5°, rendering as large bright rectangles. A **user** caught it by looking at the maps. See WORKFLOW-ISSUES.md 2026-07-28.
- Scalar diagnostics are a weak substitute even when aimed at the right question: exact-tie fractions, origin-aligned modulo-gradient ratios and variograms *all* cleared `elm-eca`, because its blocks were smooth inside, offset from the array origin, and localized rather than a uniform global grid. One image made it obvious in seconds.
- **Per member matters.** The defect was unmistakable in a per-model panel and partly diluted in the 22-member mean. Ensemble-level review hides individual bad members.

**Required behavior**:
- Render a **per-member contact sheet** (one small global panel per member) during the §9 value-check, and view it before choosing statistics. Look for: block structure, banding/seams, unrealistic hard edges, land/ocean mask errors, hemisphere flips, and hot/cold patches unrelated to geography.
- Render the **pooled product** too (median, change, trend) and confirm the geography is plausible — mountain ranges, biome boundaries and coastlines should appear where they belong.
- **Maps are a QA gate, not a deliverable.** Do not declare a layer good, report it as verified, or run `storage.cleanup_raw` until the maps have actually been viewed. `finalize_layer` generating them is not the same as someone reviewing them.
- Say plainly whether you looked. If you have not viewed an image, report the layer as *unreviewed* rather than *verified*.
- When a visual defect is found, add an automated check for its class to `scripts/generate_qa_report.py` so the next layer fails loudly instead of relying on someone noticing.
