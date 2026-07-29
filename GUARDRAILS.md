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
| Burnt area (`burntarea`) | % of cell per interval | **sum** — accumulates over the reporting interval |
| Soil carbon (`csoil-total`) | kg C m⁻² | **mean** — a stock, not a flow |
| Species richness | count | mean or max |

**The units alone do not settle it — MEASURE the cadence semantics.** `burntarea` and
`csoil-total` are both per-cell fields on the same grid at the same monthly cadence, yet one
must be summed and the other averaged, because one is a quantity *accrued during* the
interval and the other is a *state at* the interval. Getting it backwards mis-scales the
whole layer by 12×, silently and uniformly, in a way no invariant check will catch.

**Preferred technique — cross-check against a second cadence.** When a model publishes the
same run at two cadences, the relationship between them is *observable* rather than
arguable. ISIMIP3b `classic` publishes `burntarea-total` both daily and monthly: each
published monthly value equals the **sum** of that month's daily values (agreement 1e-6,
r = 1.00000000 in all 12 months), which settles it. Look for a dual-cadence member before
reasoning from units or from a sibling variable's precedent.

**Two implementation traps when summing:**
- Sum with `skipna=False`. With `skipna=True` an all-NaN ocean cell sums to `0.0`, which is
  *finite*, so the ocean is admitted as zero-valued land and drags every spatial statistic
  down. (The same bug in analysis code produced a spurious constant ratio of
  259200/66644 = 3.889 before it was caught.)
- Drop years that lack a full 12 months rather than emitting a partial sum — a partial sum
  is a silent under-count that looks like real data.

**Do not clamp a summed quantity to its nominal per-interval ceiling.** Annual burnt area
legitimately exceeds 100% where a cell reburns within the year; clamping `upper_ci` at 100
drives it *below* the median there and breaks CI ordering.

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
- **Enumerate families, don't guess members.** The Lange 2020 extreme-event exposure family is **twelve** variables (corrected 2026-07-28), not six: `le{d,r,w,c,h,t}` land-area exposed, paired with `pe{d,r,w,c,h,t}` population exposed — `d` drought, `r` river flood, `w` wildfire, `c` crop failure, `h` heatwave, `t` tropical cyclone. All ISIMIP2b, annual (data nature varies by member: `led` binary, `let` fractional — value-check each, never assume from a sibling). Authoritative table in `config/isimip_search_catalog.yaml` → `search_results.drought.exposure_lange2020.family`. Consult it before searching TC/heatwave/flood/wildfire/crop-failure exposure.
- **Enumerate `DerivedOutputData/` publications, not just variable names across rounds.** A product can be **re-issued in a later round under a different publication directory with different variable names**, so searching the old variable name returns nothing and reads as "absent". The Lange 2020 exposure concept WAS re-issued for ISIMIP3b, split across two publications: `ISIMIP3b/DerivedOutputData/**Heinicke2026**` (`driedarea`, `floodedarea`) and `.../**Zantout2025**` (`heatwave`, `wildfire`, `cropfailure`) — hazard words, not `le*` codes. The catalog asserted "ISIMIP3b/SSP version of this family NOT found (0 hits)" for four days on that basis, and the drought layer would have been built on RCP2.6/6.0 when SSPs existed. **`curl -s https://files.isimip.org/{round}/DerivedOutputData/` is a 2–4 entry listing — always run it** before concluding a derived product has no newer-round equivalent. See WORKFLOW-ISSUES.md 2026-07-28.
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
- **Same unit ≠ needs normalization; different unit ≠ always poolable.** If ensemble members share a unit but differ in magnitude, that spread is genuine model uncertainty — pool in raw units (no normalization) and let it widen the CI. Reserve robust z-score normalization for members on genuinely different scales (e.g. water-index TWS). Decide per variable and record it.
- **The declared unit string can never gate a pooling decision — only measured magnitudes can.** For `burntarea` (2026-07-28) `clm45`, `orchidee` and `orchidee-dgvm` all declare `units='%'` while sitting **~1000× lower**, on a 0–1 fraction scale. Worse, the mis-scaling can occur **within a single model, across GCMs**: `classic`'s `2015soc-from-histsoc` burntarea files are fraction-scaled for `gfdl-esm4` (monthly mean 0.0032, max 0.32) and percent-scaled for `ukesm1-0-ll` (mean 0.31, max 35.4), with identical `units` and `long_name`. **Required**: before pooling, print each member's magnitude and compare *every GCM within each model*; a ~100× sibling gap is a unit error, not model spread. Anchor the ensemble against an external observation where one exists (GFED4 ≈ 348 Mha yr⁻¹ ≈ 3.5 Mkm² yr⁻¹ global burned area) — a member three orders of magnitude off observation is broken, not merely cold.
- **Determine a monthly variable's cadence semantics empirically: does it ACCUMULATE (annual = sum) or is it a STATE (annual = mean)?** Never carry the choice over from a sibling layer. `csoil` is a stock → annual **mean**; `burntarea` is burnt area *during* the interval → annual **SUM**, and using the mean would under-scale the layer 12×. Verified 2026-07-28 the only reliable way: `classic` publishes the same run at both **daily and monthly** cadence, and each published monthly value equals the **sum** of that month's daily values (agreement 1e-6, r = 1.00000000 in all 12 months). When a model publishes two cadences, use them as the ground truth. Guard the aggregation too: drop years lacking all 12 months rather than emitting a partial sum, and sum with `skipna=False` so an all-NaN ocean cell stays NaN instead of becoming a spurious 0.
- **Do not clamp a cumulative quantity to its nominal ceiling.** Annual burnt area legitimately exceeds 100% where a cell reburns within a year (`visit` 0.15% of annual values; `classic` 0.05%, a few to ~182%; decadal means to ~107%). Clamping `upper_ci` at 100 there drives it **below** the median and breaks CI ordering. Floor at 0, leave it unbounded above.
- **The soc / CO₂-sens variant that is sound for one variable can be broken for another in the same model.** Do not inherit a variant pinned by another layer — re-value-check it. `classic`/`2015soc-from-histsoc` is the variant the `csoil` layer uses, and it is the mixed-scale one for `burntarea`; `classic`/`2015soc` is uniformly percent across both GCMs and all three SSPs.
- **Verify the experiment tokens (soc / CO₂-sens) per model — a uniform treatment across the ensemble is NOT guaranteed.** For `csoil` (2026-07-25), `jules-es-vn6p3` publishes ONLY its fixed-2015-CO₂ run (`2015co2`) while `classic`/`mc2-usfs` publish the transient `default` run, and soc tokens differ (`2015soc-from-histsoc` vs the natural-veg `nat`). Enumerate each member's sens/soc token *before* committing to a "transient CO₂" (or any) ensemble; if you must mix treatments to preserve ensemble depth, get the user's call and record it (`co2_treatment: mixed`).
- **Do not assume which DIRECTION a fixed-CO₂ member biases the trend — measure it.** This guardrail previously asserted that a fixed-CO₂ model's stock trend is "muted (no CO₂ fertilization)". For `csoil` that is **backwards**, as measured 2026-07-27: `jules-es-vn6p3` (fixed 2015 CO₂) has the **strongest and most scenario-sensitive** trend of the five models (global-mean 2090s change +1.57% / −2.92% / −4.37% for ssp126/370/585), while the four transient-CO₂ models run flat-to-positive (+0.8% classic, ~0% mc2-usfs, +2.6% elm-eca, +7.1% visit at ssp585). Removing fertilization does not damp the signal; it removes the *offsetting gain*, so warming-driven decomposition dominates and the loss grows with forcing. Report a fixed-CO₂ member's per-model trend explicitly rather than describing it as muted.
- **Never report spatial resolution from `ds.sizes` or coordinate spacing — measure it from the values.** A model that runs natively coarser and is replicated onto the ISIMIP grid still reports the full 360×720 at 0.5°. Dimensions describe the container, not the information content. For `csoil` (2026-07-28) `classic` is natively **1°×1°**, replicated 2×2 with a **one-cell longitude offset**, while its header is indistinguishable from the four genuine 0.5° members. Test the exact-tie fraction between adjacent cells **at both block offsets** (an origin-aligned test found only 3.3% and nearly cleared it), or the inside-vs-seam gradient ratio per candidate block width. A variance-loss-under-coarsening test does **not** discriminate blockiness from ordinary geophysical smoothness. `scripts/generate_qa_report.py` now checks this automatically. Also record which cells are carried by a single model, since they inherit that model's native resolution unchanged.
- **A global or area-weighted statistic is blind to a defect confined to one latitude zone — check the zonal PROFILE, not just the total.** Polar cells carry almost no area, so a member can be badly wrong there and still produce the best-looking global number. For `burntarea` (2026-07-28) the 5 `visit` members report an **inverted** profile — 25.9 %/yr burnt above 75°N against 7.3% in the tropics, with visit-only Arctic islands near 100%/yr — while `visit`'s global area-weighted total (4.56 Mkm² yr⁻¹) is the **closest of the three models** to GFED4's ~3.5. Every aggregate check passed; only rendering the maps and then breaking the field down by latitude band exposed it. **Required**: report land-mean by latitude band for every layer (`generate_qa_report.py` now does, in JSON and HTML), and set `zonal_expectation` where a hazard has a known geography so a polar band exceeding the tropical band raises a warning.
- **A DERIVED product's metadata is no more trustworthy than a model's — often less.** The temptation is to treat `DerivedOutputData` as curated and skip the value-check. Both failure modes appeared in one publication family (2026-07-28): ISIMIP3b `Heinicke2026` `floodedarea` declares `long_name="Exposed Area Share"`, `units="1"` and `definition="Flood fraction from cama flood"` while being strictly **binary {0,1}** — an occurrence flag, so its decadal mean is flood *frequency*, not area — and its mask covers **94.7% of the globe including ocean** (mid-Pacific and Antarctica read `0.0`, while Sahara and Greenland are `NaN`), so global statistics computed on it silently include the Pacific. Meanwhile ISIMIP2b `Zimmer2023` `fldfrc` carries **no `units` and no `long_name` at all**, only `_FillValue`. Note the mask defect was **per-variable**: `driedarea`, same product and same model, is correctly land-masked. Check the mask of every variable, not of the product.
- **When a value-check is challenged, RE-VERIFY at the byte level — do not re-assert it.** A plausible counter-argument ("the sibling variable is fractional, so this one should be") deserves evidence, not repetition. The escalation that settles it: (1) confirm the download is byte-identical to the ISIMIP sidecar (`size` + `sha512`); (2) print the variable's `dtype` and check for `scale_factor` / `add_offset` packing that the reader may be applying silently; (3) re-read with `set_auto_maskandscale(False)` and print the **raw stored values**; (4) count unique values **per timestep across the whole record**, not pooled; (5) print a point time series at a known-active location. For `floodedarea` this gave raw `{0.0, 1.0, 1e20}`, `float32`, no packing, exactly 2 unique values in all 86 timesteps — binary beyond argument.
- **A per-cell FRACTION metric is diluted at coastlines, and the coastal ring is usually the thinnest-covered part of the ensemble.** Any "fraction of the grid cell affected" hazard (`burntarea`, `led`, `let`, `lew`, `fldfrc`) reports against the **whole** cell, so a cell that is 10% land can report at most ~10% even if all its land is affected — an asset on the land portion is understated. It compounds with mask heterogeneity: for `burntarea` (2026-07-28) coastal cells are 12.4% of land, their median is **18× below** the interior median (0.029% vs 0.520%), and **67.5% are backed by a single model** (vs 3.7% inland) because the widest land mask covers the ring alone. Measure the coastal-vs-interior split before shipping such a layer, and either mask the ring downstream or normalize by land fraction — do not present the raw coastal values as exposure.
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
- **The trend's finite mask MUST match the median's.** `process_burntarea_fire.py:236` and `process_csoil_soilcarbon.py:241` return a bare `np.zeros(med_stack.shape[1:])` for the baseline decade, which makes the **entire ocean a finite zero** where the median is NaN. QA does not catch it: it checks only that *finite* baseline trends equal zero, never that the two masks agree. Emit `0.0` where the baseline median is finite and `NaN` elsewhere — see `process_driedarea_drought.py` (`anchored_trend`), which takes the baseline finite-mask as an argument. **Known-unfixed in burntarea and csoil** (WORKFLOW-ISSUES.md 2026-07-28).
- **An automatic threshold rule must be checked against coverage before it is allowed to fire.** `make_pct_fn` in `process_burntarea_fire.py` auto-switches to the two-tier zero-inflated percentile above 2% exact zeros. For `driedarea` that rule would have fired on 3.59% zeros over the union land mask — but the figure is 0.18% over fully-covered cells, because the zeros live in cells that not all models cover (30 samples vs 450). The trigger would have measured **uneven model coverage, not a rare hazard**, and collapsed every never-dry cell onto percentile 1. When models disagree on the land mask, evaluate any such threshold on the fully-covered subset before trusting it, and record the measured value in the output attrs either way.
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
- Render a **per-member contact sheet** (one small global panel per member) during the §9 value-check, and view it before choosing statistics. `scripts/utils/contact_sheet.py` does this — call it from the processor while each member still exists separately, and pass the path to `finalize_layer(..., extra_maps=[...])` so it is linked from the map index and included in `maps_bundle.zip`. Panels are full 0.5°, one pixel per cell, **never downsampled** — downsampling would smooth away the block structure the sheet exists to reveal. Look for: block structure, banding/seams, unrealistic hard edges, land/ocean mask errors, hemisphere flips, and hot/cold patches unrelated to geography.
- Render the **pooled product** too (median, change, trend) and confirm the geography is plausible — mountain ranges, biome boundaries and coastlines should appear where they belong.
- **Maps are a QA gate, not a deliverable.** Do not declare a layer good, report it as verified, or run `storage.cleanup_raw` until the maps have actually been viewed. `finalize_layer` generating them is not the same as someone reviewing them.
- Say plainly whether you looked. If you have not viewed an image, report the layer as *unreviewed* rather than *verified*.
- When a visual defect is found, add an automated check for its class to `scripts/generate_qa_report.py` so the next layer fails loudly instead of relying on someone noticing.

---

## 12. Report a Coverage or Resolution Fraction Against a Denominator That Makes It Meaningful

**Rule**: When stating how much of the world a layer covers, divide by the domain the variable actually inhabits. For a **land** variable that is **land** — never the whole grid, which is ~71% ocean.

**Why this matters**:
- I described the `fldfrc` river-flood footprint as *"~76% of the grid is structurally NaN, i.e. only ~24% of land"* and escalated zero-filling to the user as a product decision (2026-07-28). Both halves were wrong: the second clause silently swapped denominators. The domain is 62,066 cells = **128.8 million km² = 95.7% of global land excluding Antarctica**, against 67,420 land cells for Lange2020 `led` (61,546 shared). 79.6% of domain cells lie ≥99% inside the model domain; the median is 100%. It is an ordinary ISIMIP land mask, and there was no coverage problem to decide about.
- The error is easy to make because 259,200 (360×720) is the number in front of you, and a plausible-sounding fraction of it invites a plausible-sounding conclusion. It cost the user a spurious decision point.
- The same trap in reverse hides a real defect: `floodedarea` looks *well* covered at 94.7% of the grid, but that number is high **because it includes ocean**. Only 13–15% of its "valid" cells ever flood. A coverage fraction against the wrong denominator can flatter a broken mask as easily as it can condemn a sound one.

**Required behavior**:
- Quote land coverage as **cells and km² against a reference land mask**, and name the reference. Compare against a known-good ISIMIP land variable (e.g. `led` at 67,420 cells) or `scripts/utils/land_mask.py`, not against `lat.size * lon.size`.
- Cross-check the implied area: global land is ~148.9 M km², ~134.6 M km² excluding Antarctica. A land variable claiming far less warrants investigation; one claiming much more is including ocean.
- Before escalating coverage to the user as a decision, verify the denominator. An unnecessary decision point costs credibility and the user's time.
- Emit a per-cell coverage diagnostic (`n_members` / `n_models`, and a domain-fraction field such as `floodplain_fraction` where a cell can be partly inside the model domain) so a partly-covered cell is auditable instead of merely looking low.

## 13. Establish a PFT Field's Value Convention Before Pooling It, and Match Baseline Composition Per Scenario

**Rule**: A PFT-resolved field carries a **per-model** reporting convention, and pooling two conventions averages two different physical quantities. Establish each member's convention from the values, never from the metadata. Separately: the shared 2020s baseline is only valid when ensemble composition is uniform across scenarios — otherwise the trend is partly manufactured.

**Why this matters** (all measured 2026-07-28 building `timber_cveg-tempnle_annual` / `timber_npp-tempnle_annual`):
- Two conventions exist. **Per-tile density** reports the value on the PFT's own tile area, so `sum_i(frac_i * value_i)` recovers the cell total. **Per-gridcell** is already cover-scaled. orchidee / orchidee-dgvm / clm45 are per-tile; lpjml / caraib are per-gridcell. Pooled raw, the apparent cross-model spread read **10.5× (cveg) and 177× (npp)**; harmonized and compared on a common mask the same models agree to **2.35× and 1.83×**. The first pair of numbers would have justified a normalization the data does not need.
- **`pft-` unit labels are wrong often enough to be useless.** CLASSIC and ORCHIDEE declare `%` but store 0–1 fractions; JULES, LPJmL and CARAIB store true percent. Decide from the values: a fraction cannot exceed 1.
- **A PFT's carbon can be undefined where the PFT is absent, and masks differ wildly** — clm45 writes ~9.7k cells (its tile extent) against orchidee's ~37k. Pooling the union mixes multi-model consensus with single-model marginal habitat sitting **12×** lower, which prints mask edges into the maps and distorts the percentile (a periphery cell ranks low for having fewer contributing models, not for low value).
- **Non-uniform composition fakes a trend.** rcp85 read **−0.72 kg m⁻²/dec** at the 2030s with a plausible "recovery" purely because orchidee-dgvm sits in the all-member baseline, is absent from rcp85, and is *higher* than orchidee on the retained cells.
- **Two sound decisions can compose into a third nobody chose.** "Require ≥2 models" plus "count orchidee + orchidee-dgvm as one family" silently became "require all models" for cveg, which has only two families — cutting coverage from 17,217 cells to 9,698 and UK cells from 275 to 134.

**Required behavior**:
- Decide the convention with the two-file test: a per-tile field **exceeds** the all-PFT total (clm45 does in 94.4% of cells, ratio p50 2.28); a per-gridcell field never does and collapses toward 0 as cover → 0. Confirm with the closure `sum(frac * value) / total == 1`.
- Harmonize onto one convention **before** measuring inter-model spread or deciding normalization, and compare models **on a common mask** — own-mask medians compare different cell populations.
- Converting per-gridcell → per-tile means dividing by cover, so **apply a cover floor** and drop floored cells from `n_members` rather than emitting an amplified value (below 1% cover the amplification exceeds 100×). Check the floor's cost in the regions the product is actually for, not just globally.
- Prefer per-tile when any member publishes **no** cover fraction (clm45), since that member can never be converted the other way.
- State the mask rule in terms of **distinct model families**, and check what the rule degenerates to at the actual family count before quoting coverage.
- **Verify baseline composition per scenario.** If any member is missing from any scenario, pool each scenario's baseline panel over that scenario's members. Accept that the baseline panel is then not bit-identical across scenarios; keep the percentile *reference* distribution global so percentiles stay comparable.
- Assert `isfinite(trend) == isfinite(median)` per decade before publishing — a bare `np.zeros()` for the baseline decade makes the whole ocean a finite zero (§10), and it propagates by copy-paste between processors.
