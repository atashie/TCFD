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

**Rule**: The water risk index (`waterIndexUnderlyingData_*.nc`, 20 value types) is a **completely independent data product** from the standard TCFD annualized pipeline. **NEVER apply standard TCFD pipeline concepts** (kernel smoothing, OLS/Theil-Sen slopes, percentile-of-score ranking, shared 2020s baseline) to the water index workflow.

**Why this matters**:
- The standard TCFD pipeline produces the variables defined in [OUTPUT-SPEC.md](OUTPUT-SPEC.md): `median`, `lower_ci`, `upper_ci`, `percentile`, `ols_slope`, `sen_slope`, `n_members`, `n_models`
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
- **Enumerate families, don't guess members.** The Lange 2020 extreme-event exposure family is **twelve, not six**: six hazards — `led` (drought), `leh` (heatwave), `lew` (wildfire), `ler` (river flood), `lec` (crop failure), `let` (tropical cyclone) — each a **land**-area-exposed code paired with a `pe{d,h,w,r,c,t}` **population**-exposed twin. All ISIMIP2b, annual (data nature varies by member: `led` binary, `let` fractional — value-check each, never assume from a sibling). Authoritative table in `config/isimip_search_catalog.yaml` → `search_results.drought.exposure_lange2020.family`. Consult it before searching TC/heatwave/flood/wildfire/crop-failure exposure.
- **A variable can also have a non-exposure representation.** "Wildfire" is both the Lange 2020 `lew` *exposure* member AND the direct `burntarea` burnt-area-fraction fire output (ISIMIP2b `biomes`, in %); enumerate all products (see `search_results.wildfire`) and present the trade-offs rather than assuming the exposure family is the only option.
- **List EVERY intermediate directory level; never path-guess past one.** An enumeration that skips a level cannot support an absence claim about anything below it. Measured 2026-08-11: the walk listed `ISIMIP3b/InputData/`, saw `climate/`, and jumped straight to `climate/atmosphere/` — skipping `climate/` itself, which is where **`tropical_cyclones/`** lives. The published inventory therefore claimed no SSP tropical-cyclone product existed, when 3b ships the newest TC hazard in the repository (MIT per-storm wind footprints, Frieler et al. 2025). The user supplied the correction.
- **Scope every absence claim to exactly the family you enumerated.** In the same review it was *correctly* established that the Lange 2020 exposure family has no ISIMIP3b re-issue — and that true finding was then allowed to stand in for "no SSP TC data exists". A family having no re-issue is not the hazard having no newer data, and newer data need not live in `DerivedOutputData/` at all.
- **Open a directory before inferring from its name.** `ISIMIP3b/DerivedOutputData/TipESM2025/MIT/` reads like the tropical-cyclone group (MIT = Emanuel's institution) but holds **water** models (CWATM, H08, JULES-W2, MIROC-INTEG-LAND, …).
- **Project the variable FIELD; never let a believed-in token narrow a harvest.** A correct, complete enumeration that is filtered by a token you have not yet observed produces a false negative *wearing the authority of the full walk* — the worst kind, because the receipt (150 directories listed, 0 empties) looks impeccable. Measured 2026-08-11: a sugarcane search grepped `sgc` — the code **our own catalog** listed under `crop_codes.isimip3b` — across every agriculture directory in all four rounds and matched zero files. Sugarcane output is `sug` (ISIMIP2a/2b); `sgc` exists only in the ISIMIP3b crop-calendar `InputData`. Only a parallel `$(NF-4)` vocabulary dump surfaced the real layer. **Rule**: reduce a harvest by projecting the filename's variable field to a distinct vocabulary, then match your target offline. A pre-filter can confirm presence; it can never establish absence, and its empty result is `UNVERIFIED` (§11), not a negative.
- **Count coverage PER GCM, never pooled — and check the structure cache before walking.** Pooling members before counting makes a variable published for one GCM of five look like a full ensemble: measured 2026-08-12, `npp-tene`/`gpp-tene`/`fapar-tene` exist in the ISIMIP3b `fire` sector for **ukesm1-0-ll only**. Report `models × GCMs × scenarios` as a matrix, not as a set union. Before enumerating, read `repository_structure_cache` in `config/isimip_search_catalog.yaml` — it holds the verified model/GCM roster and the complete per-model class vocabulary for the 3b vegetation sectors and 2b biomes, and it records that `fire` and `permafrost` are byte-identical republications of `biomes` (matching Content-Length **and** ETag), so walking all three costs hours and adds nothing.
- **Vocabulary is per (round, product) — record which one you observed.** `InputData` protocol vocabulary is a superset of what models publish: the ISIMIP3b crop calendar defines 20 crops, ISIMIP3b models publish 11. Codes also drift across rounds for the same quantity (`sug`/`sgc`, `ben`/`bea`, `whe`/`swh`+`wwh`, `ric`/`ri1`+`ri2`, `csoil`/`csoil-total`). A code inherited from another round or from `InputData` is a hypothesis, not a lookup key.
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
- **Verify the experiment tokens (soc / CO₂-sens) per model — a uniform treatment across the ensemble is NOT guaranteed.** For `csoil` (2026-07-25), `jules-es-vn6p3` publishes ONLY its fixed-2015-CO₂ run (`2015co2`) while `classic`/`mc2-usfs` publish the transient `default` run, and soc tokens differ (`2015soc-from-histsoc` vs the natural-veg `nat`). Enumerate each member's sens/soc token *before* committing to a "transient CO₂" (or any) ensemble; if you must mix treatments to preserve ensemble depth, get the user's call, record it (`co2_treatment: mixed`), and note the biased member. **Do NOT assume which DIRECTION a fixed-CO₂ member biases the trend — measure it.** This bullet asserted until 2026-08-15 that "a fixed-CO₂ model's stock trend is muted (no CO₂ fertilization)", which is the deduction from the mechanism's name and is measurably backwards on the layer it was written from. Measured on `csoil` ssp585, 2090s−2020s over each member's own footprint: `jules-es-vn6p3` (**fixed** CO₂) **−4.37%**, `lpjml5-7-10-fire` −2.75%, `mc2-usfs-r87g5c1` −0.05%, `classic` **+0.79%**. The fixed-CO₂ model is the **largest relative loser** — removing fertilization removes litter input, so the decline gets *stronger*, not weaker. A mechanism tells you a term exists; only a measurement tells you its sign and size.
- Record the verified data nature in the catalog (`data_nature`) and in the output metadata.
- Never copy a sibling variable's processor attributes verbatim — re-derive the framing from the verified values.

**Related techniques** (documented in WORKFLOW-ISSUES.md 2026-07-24 and 2026-08-11; patterns, not hard rules):
- **Thin ensembles** (few members, e.g. 1 impact model × 4 GCMs): apply spatial smoothing (5×5 exponential-decay, cos(lat)-scaled, normalized over non-NaN land neighbours) to each member's **annual** map before any pooling. Reference: `scripts/process_let_cyclone.py`.
  - **The decay length `L` is a per-layer measurement, not a constant.** `L=0.7` keeps **32.1%** of the weight on the centre cell and preserves sparse structure; `L=2.5` keeps **8.1%** and dissolves it. On `let` the sharp kernel left the raw data reading as one-cell-wide storm tracks (roughness 0.389 raw → only 0.142). Measure roughness — mean |cell − its 4-neighbour mean| over exposed land, normalized by the mean — and pick `L` from it.
  - **Prefer re-weighting inside the existing footprint over widening it.** Widening asserts that land further away is exposed: on `let`, 5×5→7×7 added 1,542 exposed cells and 9×9 added 3,940, whereas `L=0.7`→`L=2.5` at fixed 5×5 added **none** while cutting roughness to 0.044. Anchor any radius physically (2 cells = 111 km ≈ the hurricane-force wind radius) rather than by eye.
- **Zero-inflated hazards**: two-tier percentile (zeros → 1; non-zeros ranked against the non-zero 2020s baseline → [2,100]) so the exposed gradient isn't crushed into the top quantiles.
- **EXTREMELY zero-inflated hazards** (`let`: 97.84% exact-zero annually) additionally take the **third decadal-statistic branch**, `pooled_mean_zero_inflated` — mean ± 1 SD on a *continuous* field. The median branch left 2,684 exposed cells vs 15,122, erasing 93% of exposed land. This is a **declared** deviation: measure the median branch's exact-zero share, record it in `decadal_statistic_rationale`, and never take this branch merely to improve contrast (`burntarea` at 29.2% does not qualify). See OUTPUT-SPEC.md "The third branch".
- **Judge slope agreement on ACTIVE cells only** (either slope non-zero). Permanently-zero cells have a genuinely zero slope under both estimators, so counting them inflates agreement: on `let` the all-cell view reads 73% agreement / 99.2% Sen-zero while the active-cell view reads **3.0% / 97.0%**.

---

## 10. Trend Must Be a Decadal Signal, Not a Within-Decade Annual Slope for Noisy Variables

> **Superseded in part by [OUTPUT-SPEC.md](OUTPUT-SPEC.md).** The *diagnosis* below still
> holds and is why the contract exists; the *prescribed fix* has changed. Layers now emit
> **two** slopes — `ols_slope` and `sen_slope` — each fitted over an **expanding window**
> from the 2020s baseline through the target decade, via `scripts/utils/decadal_stats.py`.
> The expanding window solves the within-decade-noise problem structurally (the 2090s
> panel fits 80 years, not 10), so the baseline-anchored two-point rate is **retired**.
>
> Two consequences of the change:
> - `trend × elapsed_decades == change map` **no longer holds** — do not reinstate that
>   check. The slopes are fitted, not differenced.
> - The baseline panel is **NaN**, not 0. A finite 0 there makes the entire ocean a
>   finite zero, which QA does not catch (it only asserts that *finite* baseline trends
>   equal zero, never that the masks agree).

**Original rule (retained for context)**: For a variable with high interannual variability (fire, precipitation extremes, floods), do **not** report the trend as the OLS slope of the *annual* series within a single decade. Use an across-decade signal so the trend is spatially coherent rather than dominated by interannual noise.

**Why this matters**:
- The legacy processors (`process_qg.py`, `process_led_drought.py`, `process_let_cyclone.py`) compute `trend[decade]` as the slope of the annual values *inside* that decade (10 points). For a noisy variable this slope is dominated by interannual noise → a spatially **spotty, sign-flipping** field, while the change map (a difference of two 10-year decadal means) is smooth. A user flagged exactly this for `burntarea` (2026-07-24).
- A "trend" that contradicts the change map is misleading — they should tell the same story.

**Required behavior** (for high-variance variables):
- Compute the trend from the **decadal-median series**, anchored at the baseline decade:
  `trend[decade] = (median[decade] − median[2020s]) / (elapsed decades)`, units *value* · decade⁻¹.
- This makes each decade's trend the rate **from the 2020s baseline to that decade** (2090s panel = full baseline→2090s trend), built on decadal means so it is exactly the (decade − 2020s) change map ÷ elapsed decades → coherent by construction (corr = 1.0 with change at the last decade).
- The baseline decade has no elapsed change → trend 0 (keeps the 2020s baseline bit-identical across scenarios).
- The generate_maps trend label is `[units decade⁻¹]`, so emit per-**decade** units, not per-year.
- Reference under the CURRENT contract: `scripts/process_csoil_soilcarbon.py` calling
  `decadal_stats.expanding_slopes`. Never hand-roll a slope in a processor.
- Verify with `python scripts/test_shared_baseline.py {processed_dir}` — it asserts the
  baseline panel is NaN and that no slope is finite where `median` is NaN.


---

## 11. A Recorded Negative Is a Claim, Not Evidence — Date It or Mark It Unverified

**Rule**: An absence claim written in our own notes (`config/isimip_search_catalog.yaml`,
CLAUDE.md, a skill, a module docstring) carries **no evidential weight** on its own. It is
only as good as the enumeration behind it, and that enumeration must be **named and dated
in the same place as the claim**. A negative without a receipt is `UNVERIFIED`, and must be
re-enumerated before it is relied on or repeated to the user.

**Why this matters** — measured, 2026-08-11:

- On **2026-07-24** the drought layer shipped on rcp26/rcp60 (`edeb174`) and the catalog
  was written the same day (`8fac010`) asserting `simulation_round: "ISIMIP2b only (NO
  ISIMIP3b/SSP version of this family)"` and *"ISIMIP3b/SSP version of this family NOT
  found (0 hits)"*. **Both were false.** ISIMIP3b publishes drought exposure as
  `driedarea` under `DerivedOutputData/Heinicke2026` — a **complete** 3 models × 5 CMIP6
  GCMs × 3 SSPs = 45-file matrix with uniform `2015soc`/`default`. Nobody had listed
  `DerivedOutputData/`; the "0 hits" was a search for the *code* `led` in 3b, and the code
  is 2b-only even though the hazard is not.
- On **2026-08-07** the `/isimip-search-download` skill gained the correct fact
  (`ac524f5`), and CLAUDE.md on **2026-08-10** (`cf335fb`). From that moment the repository
  **contradicted itself** — skill says the re-issue exists, catalog says it does not — and
  the contradiction sat unresolved for **4 days**, the false negative for **18 days**.
- On **2026-08-11** the `let` rebuild (`21d8a5a`) *edited the very string* "NOT found (0
  hits)" while fixing the tropical-cyclone notes, and left the drought false negative
  standing one block away.

**Required behavior**:

- **Every negative in the catalog carries its receipt inline**: the exact path or query
  enumerated, and the date. Use `verified_absent_on: "<date> — listed <URL>"`. If you
  cannot state the path you listed, write `UNVERIFIED` instead of the negative.
- **An UNDERSTATED POSITIVE is this same failure and is harder to catch.** A capped or
  paginated result is `UNVERIFIED` for **coverage in either direction**, not only for
  absence — record the enumeration method beside the count. Measured 2026-08-13: the catalog
  read `lec: impact_models: "gepic"` — one model — taken from a 2026-07-24 API listing capped
  at **20 rows**, and `config/hazard_taxonomy.yaml` then carried that forward as the *stated
  blocker* for the entire crop-failure family: *"`lec` is a single impact model, so its CI
  would carry no structural uncertainty."* `lec` has **two** models and the ISIMIP3b re-issue
  `cropfailure` has **eight** — the deepest ensemble in the product, blocked for weeks by our
  own note. A recorded negative at least looks like something to re-check; `impact_models:
  "gepic"` reads as a measured fact and invites no scrutiny at all. §8 already treats
  `count=1001` as truncation; a **row cap on a listing is the same failure at a different
  limit**, and so is any count taken from a first page.
- **A negative about a CODE is not a negative about a HAZARD.** `led` returning 0 hits in
  ISIMIP3b is a true statement about a controlled-vocabulary token and says nothing about
  whether the hazard was re-issued under another name. Scope the claim to the token.
  (Same shape as §8's "scope every absence claim to the family you enumerated".)
- **Two of our own documents disagreeing is a STOP-AND-RESOLVE trigger, not a precedence
  question.** Precedence (skill > CLAUDE.md > catalog) tells you which to *believe right
  now*; it does not license leaving the loser in place to mislead the next session.
  Resolve by enumeration and **write the result back to both**.
- **When editing a doc region, re-read the surrounding claims about the same fact.** A
  stale negative adjacent to the line you are fixing is the one most likely to survive.
- **Never repeat a recorded negative to the user as fact without re-verifying it**, and
  when you do repeat one, say which enumeration and date it rests on.

**Related**: §8 (never guess codes; scope absence claims), and the skill's *"List
`DerivedOutputData/` before concluding a product has no newer-round version"* — that rule
already existed and was correct; what failed was that a **stale catalog negative
contradicting it was never reconciled**.

---

## 12. Verify the Layer Is Non-Trivial WHERE THE THING ACTUALLY EXISTS

**Rule**: For any layer describing something with a known real-world geography — a crop, a biome, a fishery, an industry, a hazard with a known footprint — check the field is non-trivial at a **named list of reference locations where that thing demonstrably exists**, before processing and again before shipping. Record the sites and their values in the processor docstring. §9 asks *what the values are*; this asks *where they are*, and the two failures are independent.

**Why this matters** — measured 2026-08-11, the sugarcane layers:

- ISIMIP2b LPJmL `yield-sug-*` reads exactly **0 across the entire sugarcane belt** — São Paulo, Uttar Pradesh, Guangxi, Thailand, Pakistan Punjab, Veracruz, KwaZulu-Natal, Cauca, Queensland, Florida, Louisiana — with sentinel companions (`biom-sug` pinned at exactly 0.267 t ha-1, `plantday = matyday = 1`, i.e. no season simulated) on 19.2% of land, 0% of which has a non-zero yield. The same model's ISIMIP2a run gives Florida **19.49** and São Paulo **11.69** t ha-1.
- Two layers were built from it, and **`test_shared_baseline.py` passed every check on both** — schema, shared baseline, CI ordering, percentile range and orientation, slope masks, ensemble depth. All of that was true. Contract conformance is a statement about form and says nothing about whether the input means what its name says.
- The §9 data-nature check was run properly and still missed it: the 87% zero fraction was measured, classified as *structural*, and written up as "LPJmL grows no sugarcane there" — a true sentence whose implication was backwards. Counting zeros and testing whether they move over time does not locate them.
- The defect surfaced only because a **user looked at a map** and asked why the US Midwest out-yielded Florida.

**Required behavior**:

- Keep an explicit reference-site list per layer (5–10 places, named, with lat/lon) and print the value at each during processing. For a crop, use its top producing regions; for a hazard, its known footprint.
- **A zero at a reference site is a STOP**, not a data point to rationalise. Investigate before writing any output.
- **Cross-check the same sites in a second round or second model publishing the same quantity.** A contradiction of this size (Florida 19.49 vs 0.00) costs two lookups.
- Inspect the companion variables for sentinel signatures — a constant repeated across unrelated cells (0.267 everywhere), or degenerate phenology (`plantday == matyday`), means the model did not run there, whatever the primary variable says.
- **Do not let a contract PASS stand in for a sanity check.** When a layer passes cleanly, the remaining risk is entirely in whether the input is about what you think; spend the two minutes there.

**Related**: §9 (measure data nature, never infer it) — necessary and, as shown here, not sufficient.

## 13. List the Container You Entered — a Known Path Is Not a Search

**Rule**: Before you record anything from a directory, **list its parent** and record every
sibling. Entering a path you already know the name of is retrieval, not enumeration, and it
can only ever return what you already knew. This applies the moment you touch a path,
regardless of what question you are answering — §8's "list every intermediate level" governs
walking *down*; this governs the level you *started at*.

**Why this matters** — measured 2026-08-14, the CaMa-Flood inundation suite:

- Between 2026-07-24 and 2026-08-13 the repository held that river flood was evidenced only
  by exposure flags: `ler` (2b) and `floodedarea` (3b), both 0.5°, both scoring *departure
  from a preindustrial reference*. `config/hazard_taxonomy.yaml` recorded exactly those two
  as the `isimip_candidate` list for the family it calls **the most consequential gap in
  this pipeline**, and a blocker was written explaining why neither could ship.
- **`ISIMIP2b/DerivedOutputData/` was never listed.** The catalog records
  `enumerated: "2026-08-08 (ISIMIP2b/DerivedOutputData/Lange2020/{MODEL}/{gcm}/future/)"` —
  we walked *into* the publication whose name we already had, from a variable-code search,
  and never looked at the directory containing it. That directory has **two** entries.
  The second, `Zimmer2023`, is CaMa-Flood hydrodynamic inundation at **150 arcsec (~4.6 km,
  144× our grid)** with flood **depth in metres**, three protection variants, and **rcp85** —
  a scenario the exposure family does not have at all.
- The cost is not academic: `floodedarea` reads **0.000 across the Amazon floodplain** in all
  45 members, while `Zimmer2023` puts the Amazon main stem at **6.11 m mean depth, 47.2% of
  cells flooded** — the wettest box in a six-site check. A customer shown the first would
  have rejected it on sight, correctly.
- **The same failure then repeated inside the correction.** On 2026-08-13 a listing of the
  3a and 3b `DerivedOutputData/` roots produced the claim *"no CaMa-Flood publication exists
  under ISIMIP3a or 3b"* — written from **publication names**, not contents. Both were wrong:
  `Quesada-Chacon2026` (3a) *is* CaMa-Flood, and `TipESM2025` (3b) is CaMa-Flood
  `fldfrcmax` at 15 arcmin with **ssp126/370/585** and 32 members per scenario. The catalog
  had `TipESM2025` filed as "water models" from a 2026-08-11 note that identified its models
  and never asked what it published.

**Required behavior**:

- **List the parent before recording the child.** One request. `curl -s <parent>/ | grep -oE
  'href="[^"]+/"'`. Do it even when you arrived at the child from a search hit, a citation,
  or a path in our own notes.
- **A publication is identified by its FILES, not its name or its model list.** Descend to
  the first directory containing `.nc`/`.nc4` and read two filenames. `Zimmer2023`,
  `TipESM2025` and `Quesada-Chacon2026` are all CaMa-Flood and none of them says so above
  the file level. Directory names are institutions and first authors; they tell you nothing
  about the variable.
- **A candidate/option list is a positive claim of completeness.** It needs the same receipt
  a negative does: which directory was listed, and when. Without one it is `INCOMPLETE`, not
  a list. §11 covers absence claims and understated positives; an option list is the third
  form and the most decision-bearing, because work is planned off it.
- **Resolution and filename grammar are per-PUBLICATION properties, not round constants.**
  Observed in one hazard: 0.5°, 15 arcmin, 150 arcsec, 300 arcsec and 900 arcsec, with field
  counts differing between publications so that `$(NF-4)` reads the variable in one and the
  *resolution* in another. Never carry an offset or a grid assumption across publications;
  read the `.json` sidecar's `specifiers` block, which names fields explicitly.
- **A variable-code or family search cannot discover a vocabulary you do not already have.**
  `flddph`, `fldfrc` and `fldfrcmax` were unreachable by any search keyed on `le*` codes.
  Directory listing is the only discovery mechanism; code search is confirmation.

**Related**: §8 (never guess specifier codes; list every level going down), §11 (a recorded
negative needs a receipt — as does a recorded option list).
