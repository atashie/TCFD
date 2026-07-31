# Workflow Issues Log

This document tracks workflow mistakes and their resolutions. Each incident documents what went wrong, why it matters, and how to prevent recurrence.

See [GUARDRAILS.md](GUARDRAILS.md) for the rules derived from these incidents.

---

## Incident Log

### 2026-07-30: `trend` Replaced by a Theil-Sen Slope — and the First Attempt Was Broken by Zero-Inflation — RESOLVED

**What happened**: after the p-value backfill (below), the user asked to *"serve sen slopes
on all 9 layers instead of trends"*, standard for all future layers. Two premises in the
request were factually off and were corrected before acting:

- *"instead of trends derived from linear regression"* — the trend was never a regression.
  It was the baseline-anchored **two-point rate** `(median[decade] − median[2020s]) /
  elapsed decades`.
- An earlier message assumed Sen slopes were **already** being served and offered to delete
  `theilsen()` as dead code *"if we run sen slope using a different function"*. Verified
  against S3 first: all nine layers declared the anchored rate, all six processors called
  `anchored_trend()`, and `theilsen()` was the **only** Sen implementation for these layers.
  Deleting it would have removed the sole Sen capability. **Check the premise of a cleanup
  request before executing it** — the conditional was false.

**The first implementation was wrong, and a user question caught it after three layers had
already been published.** Fitting Theil-Sen to the ensemble-mean **ANNUAL** series looked
strictly better: it uses every year and puts `trend`, `trend_tau` and `trend_pvalue` on one
series (sign agreement hit 100.000%). The user then asked whether a hazard that is
historically mostly zero could report a zero trend despite the frequency doubling, "because
most year-to-year changes are still 0 to 0". It can, and it did:

**Theil-Sen is a MEDIAN of pairwise slopes, so once more than half of all pairs are tied the
slope is EXACTLY 0** regardless of the real change. Measured share of cells with an
exactly-zero slope at the 2090s:

| layer | annual series | **decadal series** | (old anchored rate) |
|---|---|---|---|
| `driedarea` ssp126 | **91.3%** | 13.7% | 4.5% |
| `driedarea` ssp585 | 56.2% | 10.3% | 3.6% |
| `burntarea` ssp585 | 14.4% | 10.8% | 9.8% |
| `csoil` ssp585 | 3.7% | 3.7% | 0.4% |

On the annual series **25.1% of `driedarea` ssp585 cells reported `p < 0.05` beside a slope
of exactly zero**, and 81.0% of ssp126 cells reported zero while the value had materially
changed. The flood layers (~93% exact zeros in any single year) would have been worst of
all; the rollout was stopped before reaching them.

**Resolution**: `trend` is `trend_significance.theilsen_decadal()` — Theil-Sen on the
**DECADAL median** series over an expanding window from the baseline decade, in value per
decade. All nine layers republished; all six processors emit it natively. The p-value stays
on the **annual** series (rank-based, so ties cost it nothing, and n=20…80 rather than 2…8),
so `trend` and `trend_pvalue` are deliberately fitted on different series — sign agreement
85–99% on significant cells rather than 100%. **A reduction, not a cure**: 10–14% of hazard
cells still return exactly 0, so `generate_qa_report.py` warns when a zero slope coincides
with a significant p-value (0.03–1.1%, down from 25.1%).

**Verified end-to-end 2026-07-31**: `process_driedarea_drought.py` was run in full with
`publish_processed_layer`/`finalize_layer` stubbed (and `storage.push` hard-failed, so a
mismatch could not overwrite the published layer), and its staged output diffed against S3.
**Bit-identical on all 10 variables x 3 scenarios — 30/30, worst relative deviation
0.000e+00**, identical finite masks, no attr diffs. That matters because the processor and
`backfill_trend_significance.py` reach `trend`/`trend_pvalue`/`trend_tau` by different
routes — the processor accumulates the annual ensemble mean member-by-member while
streaming, the backfill rebuilds it in a separate pass. They agree to the last bit only
because both add members in **sorted filename order**; float addition is not associative, so
a different order would show up as small non-zero deviations rather than exact equality. Do
not "simplify" either path's iteration order. The other five processors share this wiring
but have not had their own end-to-end run.

**Lessons worth keeping**:

- **A statistic that is right on continuous data can be catastrophically wrong on
  zero-inflated data, and the aggregate diagnostics will not say so.** Sign agreement went
  UP to 100% under the broken version — the defect made `trend` and `tau` agree because both
  were near zero. The detector that works is "exactly-zero slope beside a significant test".
- **Never infer a method from prose that also names the alternatives.**
  `trend_definition_theilsen()` ends with *"…not the GUARDRAILS S10 baseline-anchored
  two-point rate…"*, and `generate_qa_report.py` detected the method with
  `"baseline_anchored" in trend_definition`. The disclaimer matched, and a correct
  Theil-Sen layer **failed the retired change-map identity in all three scenarios**. Detect
  from the structured `trend_method` attr, and let it win over any prose test.
- **Score an agreement metric on the SIGNIFICANT subset.** This was got wrong twice in one
  session. Over all cells, non-trending cells contribute an arbitrary sign: the check read
  74–77% at the weakest forcing and failed three *correct* layers, while their significant
  subsets read 85–99%.
- **A "replacement" needs its own gate.** `verify_append_only` demands bit-identity for
  every pre-existing variable, so `trend` had to be exempted — and a bare exemption would
  let a silent no-op pass as a successful replacement. It now *requires* `trend` to have
  changed whenever the file still declares the anchored method, while still allowing an
  idempotent re-run of an already-converted layer.
- **Units: `theilsen_decadal()` fits against the decade INDEX**, so the slope is already per
  decade. The annual variant correctly multiplies by `window_years`; doing that here would
  inflate every published trend tenfold. That is why the decadal entry point exists
  separately rather than reusing `theilsen_expanding()` with a decade axis.

---

### 2026-07-30: No Layer Carried a Trend p-value, and the First Two Methods Chosen Were Both Wrong — RESOLVED

**What happened**: Preparing the customer-delivery extractor, `scripts/utils/export_formatter.py`
was found to read a `trend_pvalue` variable (mapped to the legacy `Decadal_Trend_Significance`
column) that **no layer emitted**. It resolved to NaN, `calculate_trend_aggregated()` returned
`0.0`, and every delivered row would have read *"trend not significant"* — a silent, uniform,
wrong answer across all 9 layers. `CLAUDE.md` had advertised a "significance" value class since
the beginning; it had never existed in the new processors.

The published `trend` cannot supply one: it is a baseline-anchored rate built from two decadal
numbers (GUARDRAILS §10), so it has no residual and no degrees of freedom.

**Three method errors, each caught by measurement rather than reasoning**:

1. **Spearman was proposed, and is degenerate at these sample sizes.** A perfectly monotonic
   8-point series gives `scipy.stats.spearmanr` p = **exactly 0.0** (infinite t-statistic);
   at n=5 it gives 1.4e-24. Kendall's exact p at n=8 is 4.96e-5. Switched to Kendall.
2. **An earlier claim that MK's smallest two-sided p at n=8 was ≈0.011 was wrong by ~200×** —
   it is ≈5e-5. That error had made a pooled-decadal test look underpowered when it was not.
3. **The chosen pooling — stacking every member-year as an independent observation — destroys
   the signal on continuous layers.** It looks more powerful (n = years × members) but the
   sample is then dominated by between-model **level offsets** rather than by time. Measured on
   `csoil-total` ssp585, the between-member level SD is **68.7×** the within-member interannual
   SD. Share of land reaching p<0.05 over 2020–2039:

   | pooling | n | land p<0.05 |
   |---|---|---|
   | stack raw member values | 240 | **14.4%** |
   | de-mean each member, then stack | 240 | 78.4% |
   | ensemble-mean annual series | 20 | **83.0%** |

   On `driedarea` (binary; level ratio 0.56×) all three agree within a few points — so the
   error **hides on the layers where it is harmless and bites on the ones where it is not.**
   Inter-member rank correlation of the annual series is only +0.10, so pseudo-replication —
   the concern that had been raised — was the *lesser* problem.

**Resolution**: `scripts/utils/trend_significance.py` — two-sided tie-corrected Mann-Kendall
on the **ensemble-mean annual series**, expanding window anchored at the baseline decade
(2030s tests 2020–2039 at n=20, rising to n=80 at the 2090s). Because the window expands,
`trend_pvalue[2090]` **is** the long-term p-value, so one variable fills both legacy columns
and `export_formatter.py` needed no change. 36 unit tests check it against an independent
textbook MK implementation and tau-b against scipy. Backfilled into all 9 layers in place by
`scripts/backfill_trend_significance.py`.

**Two traps worth keeping**:

- **Exactness matters at small n and is infeasible at large n.** An exact tie-aware null is
  enumerable for n≤8 (it depends only on the tie pattern — a partition of n — so 65 tables
  cover n=2..8), and it matters: with ties `scipy.stats.kendalltau` falls back to the
  asymptotic normal and returned 0.127 where the exact value is 0.250, **2× too significant**.
  At n=20–80 enumeration is impossible and unnecessary. The method therefore changed *because*
  the sample size changed — the exactness decision was made under the old n and had to be
  revisited, not carried forward.
- **The reconstruction gate earned its place immediately.** Requiring the decade-mean of the
  rebuilt annual series to reproduce the published `median` caught a missing fraction→percent
  scale on all three `fldfrc` layers on the first run — wrong on **every cell by 100×**, and
  invisible in the p-value itself because MK is rank-based and scale-invariant. It would have
  shipped a layer whose significance was right and whose provenance claim was false. Gate on
  the *share* of cell-decades over tolerance, not the max: `csoil` (36–91 of ~495k) and
  `burntarea` (0–55 of ~472k) legitimately differ where a member's coverage varies **within**
  a decade, because mean-over-years-of-mean-over-members is not mean-over-members-of-mean-over-years.
  Gating on the max would have rejected them; gating on the median alone would have missed
  fldfrc.

**Also found**: the 2020s baseline is a **cross-scenario** construct in several processors
(`process_driedarea_drought.py:280-284` averages each member across all three scenarios), so
the published `median[2020]` does not equal any single scenario's own 2020s mean. The
significance test uses each scenario's own annual trajectory, which is the honest test of that
scenario, so the baseline decade is reported but **not** gated.

See GUARDRAILS.md §15.

---

### 2026-07-30: The `median` Variable Has Always Held a MEAN — OPEN, deferred by the user

**What happened**: While specifying the windy-days layer the user asked "I thought we were
outputting both the median and the mean?". Checking the code rather than answering from the
docs: **no processor emits a `mean` variable at all**, and every processor writes a *mean*
into the variable named `median`:

```
process_csoil_soilcarbon.py:466   median[i] = np.nanmean(layer, axis=0)
process_let_cyclone.py:319        median[i] = np.nanmean(layer, axis=0)
```

and every layer's `statistic` attribute says so plainly — `decadal_mean_soil_carbon_kg_per_m2`,
`decadal_mean_exposed_area_fraction_spatially_smoothed`, `decadal_mean_annual_burnt_area_percent`,
`decadal_mean_exposure_frequency`, `decadal_mean_annual_flooded_area_percent`.
`grep -l '"mean": ('  scripts/process_*.py` returns nothing.

**Impact**: The NetCDF variable name contradicts its content in all 9 published layers. A
consumer reading the variable name alone gets a mean while believing it has a median; a
consumer reading `statistic` is correctly informed. For `let`/`led` the mean is *deliberate*
and documented (the field is ~97% zeros, so a true median is identically 0 and useless);
`csoil`/`burntarea`/`fldfrc` simply inherited the same slot.

**Root cause**: A legacy 6-slot schema (see `generate_maps.py`: `0 median, 1 percentile,
2 trend, 3 significance, 4 lower_ci, 5 upper_ci`) fixed the *names* before the statistics were
settled per hazard. The same legacy schema documents `lower_ci`/`upper_ci` as **Q25/Q75**,
while `csoil`/`burntarea` implement them as mean ± 1 inter-member SD — a second divergence
between documented and actual semantics.

**Correct action**: Either add a real `mean` field and let `median` hold a true median (correct
but breaking for every consumer), or keep one field and fix the naming where it is documented.

**Status**: **OPEN — user decision 2026-07-30 to note it now and address it later.** Do not
"fix" it opportunistically inside an unrelated layer; it is a cross-layer schema change.
The windy-days layers ship with the existing convention (`median` holding the decadal mean,
`statistic` stating it) so they stay directly comparable with the other 9 layers.

**Second half of the same issue — `lower_ci` / `upper_ci` semantics.** The legacy 6-slot
schema in `generate_maps.py` documents them as **Q25 / Q75**; every processor implements
**ensemble mean ± 1 inter-member SD** (`csoil`'s `ci_definition` states this explicitly).
Those are different quantities: the SD form is symmetric about the mean and can be inverted
by clamping, the quantile form is not and cannot. The user initially specified Q25/Q75 for
the windy-days layers, then **decided 2026-07-30 to keep mean ± 1 inter-member SD for now,
so the five new layers stay consistent with the existing nine**, and asked that the
divergence be recorded here for fuller treatment later.

**Deferred decision, to take together**: whether `median` becomes a true median with a new
`mean` field beside it, and whether `lower_ci`/`upper_ci` become genuine quantiles. Both are
cross-layer schema changes affecting all published layers plus
`scripts/utils/export_formatter.py` and the legacy 28-column schema; neither should be done
piecemeal. Until then, **trust the `statistic` and `ci_definition` attributes, never the
variable names.**

**Rule created**: none yet — pending the schema decision.

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
- Alignment issue: **FIXED** (2026-07-28) — took the `go.Heatmap` option from the table above. `create_map_figure()` now draws one rectangle per grid cell with `zsmooth=False`, so cells are cells at every zoom instead of cell-centre point markers, and the marker/cell-size mismatch is gone. It was **not** low priority after all: markers were being drawn 2px wide at 0.83px per cell on a 300px canvas (~6× overplot), which not only looked misaligned but **actively masked the real resolution of the data** and contributed to a coarse ensemble member going unnoticed (see 2026-07-28 below). The table's stated cost — "loses interactivity" — did not materialise: hover, zoom and the colorbar all work; what was lost is the plotly `geo` basemap (coastlines/ocean fill), and the land mask outlines the continents itself.

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

### 2026-07-27: csoil Ensemble Expanded 12 → 22 Members; Two Prior Claims Refuted

**Context**: Re-reviewed the ISIMIP repository for soil-carbon datasets *without* trusting the existing documentation. Enumerated 98 directories (66,902 filenames) and read 17 file headers. The user then chose to keep the stock-only `csoil-total` layer on ISIMIP3b and expand it 12 → 22 members by annualizing the two monthly submissions.

**What the re-enumeration corrected**:

1. **A parallel harvest silently returned empty listings.** My first pass fetched 22 directories concurrently; the file server rate-limited **13 of them to zero results**, which is indistinguishable from a genuine "no data" answer. I built and reported member matrices from that truncated data before catching it. → Harvest **serially with retries and assert non-empty per directory**. An empty autoindex listing is a *failure signal*, not a finding. (Extends GUARDRAILS §8.)

2. **The variable is renamed between rounds.** ISIMIP**2b** publishes bare `csoil`; ISIMIP**3b** publishes `csoil-total`. Searching `csoil-total` in 2b returns nothing, which reads as "absent" and isn't.

3. **The 2b model list was still short by two.** The 2026-07-25 entry listed 9; the file server has **11** — it missed **CLM45** (annual) and **VEGAS** (monthly).

4. **ELM-ECA cross-publishes under two sectors, byte-identically.** `biomes/ELM-ECA/…csoil-total…` and `permafrost/ELM-ECA/…csoil-total…` share an sha512 and size (237,500,430 B). Harvesting both would have **silently double-weighted** one model. Only the biomes copy is ingested.

5. **`csoillayer-total` looks ideal for 0–30 cm SOC accounting but is unusable.** CLASSIC has 20 layers, but the depth coordinate is integer-rounded so the top ten layers all read "0 m" (no 30 cm cut is possible) and layer 0 holds ~81% of the column; JULES's copy has `depth=1` and column values ~1042 kg C m⁻², ~100× its own `csoil-total`.

6. **No `ssp245`** exists for any ISIMIP3b biomes variable — only ssp126/370/585.

**The "muted trend" claim was backwards.** Both GUARDRAILS §9 and the catalog stated that `jules-es-vn6p3`'s fixed-2015-CO₂ run makes its soil-carbon trend *muted*. Measured per-model global-mean 2090s change proves the opposite — JULES has the **strongest and most scenario-sensitive** response of the five:

| model | CO₂ | ssp126 | ssp370 | ssp585 |
|---|---|---|---|---|
| classic | transient | +0.65% | +0.77% | +0.77% |
| mc2-usfs | transient | +0.40% | −0.11% | −0.05% |
| elm-eca | transient | +2.01% | +2.43% | +2.60% |
| visit | transient | +6.70% | +7.31% | +7.14% |
| **jules-es-vn6p3** | **fixed 2015** | +1.57% | **−2.92%** | **−4.37%** |

Removing CO₂ fertilization does not damp the signal — it removes the *offsetting gain*, so warming-driven decomposition dominates and the loss scales with forcing. Corrected in GUARDRAILS §9, the catalog, CLAUDE.md and the processor metadata.

**Consequence: the headline sign flipped.** Global-mean 2090s change moved from **+1.1% / −1.4% / −2.2%** (12 members) to **+2.50% / +1.57% / +1.24%** (22 members) for ssp126/370/585. Scenario ordering stayed monotonic, but ssp370/ssp585 flipped from net loss to net gain, because JULES — the ensemble's main loss source — fell from 42% to 23% of the weight while the two added models both accumulate carbon. **Recomputing the old 3-model ensemble from the same raw files reproduced the previous published layer to within 0.05pp**, confirming this is ensemble membership and not a processing change. The sign of the global mean should be read as **contested across models** (3 of 5 gain, 1 flat, 1 loses strongly), not settled.

**Two caveats recorded rather than silently absorbed**:
- **Tail-driven pooling.** "No normalization" was justified on medians agreeing within 1.8×, but members are pooled with a **mean**, and their global means span **2.6×** (6.70 classic … 17.72 elm-eca, the latter 2.3× its own median) because elm-eca and visit carry much fatter peat tails. Those extremes are **tropical** (median latitude of the top 0.1%: 5.8 N, 31.9 N) — deep tropical peat, not boreal permafrost. The baseline global mean rose 8.78 → 10.99 kg C m⁻² purely from adding those two models. A cross-member median would be robust to this; the mean was kept for cross-layer consistency.
- **Heterogeneous land masks.** The 5 models do not share a mask (58,714–67,647 of 259,200 cells; each model's GCMs match exactly). ~81% of land carries all 22 members, 5.9% is single-model, and 2 cells are single-member. Cells are **retained, not masked** (user decision: a decade supplies 10 annual samples per member, so spread stays estimable). New `n_members` / `n_models` variables make coverage auditable. Note the CI is a spread over member *decadal means* (n = members), not over year × member samples — so on a single-model cell it is a 2-sample spread, and on the 2 single-member cells it collapses to zero width.

**Also fixed**: years are now decoded with **cftime** rather than days/365 arithmetic. The old arithmetic put December of a monthly member at ~year+0.96, where rounding pushed it into the following year — misassigning one month in twelve. Harmless for the 3 annual models, corrupting for the 2 new monthly ones.

**Verified**: 37/37 QA checks — CI ordering 0 violations, 2020s baseline bit-identical across scenarios, `trend × span == change` to 2.4e-07, percentile inverted (corr = −0.956), scenario ordering monotonic, counts ≤ 22/5 and NaN off-land. Raw: 66 files / 7.47 GB ingested to S3 with sha512 + source URLs recorded.

**Processor**: `scripts/process_csoil_soilcarbon.py`.

---

### 2026-07-28: Declared Grid ≠ Effective Grid — CLASSIC Is Natively 1°, Not 0.5°

**What happened**: Reviewing the published `soilcarbon_csoil_annual` maps, the user noticed the projections looked "blocky" — large boxes of general trend with finer detail inside them — and asked me to re-check each ingested dataset's spatial resolution. I had already "verified" the resolution during the 2026-07-27 review by reading 17 file headers, all reporting **0.5° / 360 × 720 / lat −89.75…89.75**, and reported that as uniform. That check was worthless for the question asked.

**Why it was worthless**: a model that runs natively coarser and is replicated onto the ISIMIP 0.5° grid **still reports 360 × 720**. Dimensions and coordinate spacing describe the container, not the information content. Only the values can reveal the effective resolution. This is GUARDRAILS §9 ("trust the values, not the labels") applied to geometry rather than to units — I had applied it to units and magnitudes and not thought to apply it to the grid.

**Measured result**: `classic` is **natively 1.0° × 1.0°**, replicated 2 × 2 onto 0.5° **with a one-cell longitude offset**:

- 100% of longitude pairs at offset 1 identical; 100% of latitude pairs at offset 0
- 100% of 2 × 2 blocks constant at offset (lat=0, lon=1)
- 99.35% of longitude runs are exactly 2 cells; the 0.38% longer runs are constant-value desert/ice plateaus, **not** a coarser grid — nothing coarser than 1° exists
- `elm-eca`, `jules-es-vn6p3`, `mc2-usfs`, `visit`: genuine 0.5°. No GCM contributes block structure (all 5 checked, ratios 0.94–1.08)

Note `mc2-usfs` is a *coarse biome model* with *0.5° output* and the smallest land mask — a separate property that must not be conflated with output resolution.

**Two false starts worth recording**:

1. **An origin-aligned 2 × 2 constancy test found only 3.3%** and I nearly concluded there was no blocking. CLASSIC's blocks are offset by one cell in longitude, so an aligned-only test misses them almost entirely. **Always test both offsets.** The contradiction that exposed it: 53% of adjacent cells were exact ties while "only 3%" of 2 × 2 blocks were constant — arithmetically impossible for a genuinely unblocked field, and the signal to keep digging.
2. **A variance-loss-under-coarsening test (R² of k × k block means) did not discriminate at all** — all five models produced smooth declining curves (R² ≈ 0.87–0.97 at 1°, ≈ 0.6–0.7 at 6°), because it cannot separate replication from ordinary geophysical autocorrelation. Discard this approach; use exact ties at both offsets, the inside-vs-seam gradient ratio, or an FFT of the column-gradient profile.

**Impact on the published layer**: `classic` holds 2/22 = **9.1%** of the weight where all models report, leaving a measurable global 1° signature (inside/seam longitude gradient ratio **0.838**, vs 1.0 for a native 0.5° field). More visibly, the **3,085 cells (4.36% of land) where `classic` is the only reporting model are rendered at its native 1°** — scattered across all latitude bands (30% tropics, 32% N mid-lat, 29% boreal) and low-carbon marginal land (median 2.82 vs 10.07 kg C m⁻² layer-wide). This is the blockiness the user saw. The user's estimate of ~3° was the right instinct; the measured value is 1°.

**Systematic fix**: `scripts/generate_qa_report.py` now runs an **effective-resolution check** on every layer — inside-vs-seam gradient ratio at candidate block widths 2/3/4, plus a count of cells rendered by a single model. Both are warnings, not failures, since a coarse member can be a documented trade-off for ensemble depth. On this layer it correctly reports lon = 0.838 at 1°. (The 2° warning is a harmonic of the 1° signal, not an independent finding.) Recorded in the catalog `resolution` block, the processor's `effective_resolution` attribute, and the `isimip-process-visualize` skill.

**Rule**: never report spatial resolution from `ds.sizes` or coordinate spacing. Measure it from the values, at both block offsets, per member.

**Resolution (corrected)**: my first answer to the user named `classic` (1°) as the cause. That was **wrong** — it is real but minor. The dominant artifact was **`elm-eca`, effectively ~4° lat × 5° lon**: gradient seams every 10 columns (5.0°, 62 of 75 detected) and every 8 rows (4.0°), with fine variation *inside* each block. It hid from three separate tests: declared dims are identical; an exact-tie test scores it clean (only 2.9% ties, because its blocks are smooth inside rather than constant); and my first origin-aligned modulo-gradient check stopped at 2° and at k=10 even *inverted* to 1.149, reading cleaner than clean. Compounding it, `elm-eca` had the fattest tail of the five (max ~160 vs 35–70 kg C m⁻²), so its coarse boxes were biased high and rendered as bright rectangles — it alone inflated the 2020s global mean from 9.32 to 10.99 kg C m⁻² (+18%).

What finally worked was **rendering the field and looking at it**. Global-average statistics (tie fractions, modulo gradient ratios, variograms) all missed a localized coarse patch pattern; one image made it obvious in seconds. Render before trusting a scalar diagnostic.

**Resolution taken (user decision 2026-07-28)**: `elm-eca` dropped → **17 members** (classic 2, jules 5, mc2-usfs 5, visit 5). Cost 5 land cells (0.01%); all effective-resolution checks now pass with no dominant seam spacing. Trade-off recorded: single-model cells rose 5.91% → 9.23% of land and zero-width-CI cells 0.52% → 3.6%, since fewer models cover the thin margins. `generate_qa_report.py` now tests block widths to 6° **and** measures seam spacing with no alignment assumption — the check that would have caught this unprompted. Since an automated scalar check is a backstop rather than a substitute for looking, `scripts/utils/contact_sheet.py` also renders a **per-member** contact sheet at full 0.5° (embedded PNGs, shared colour scale), linked at the top of the map index and bundled — the artifact that would have made this obvious in seconds.

**Also fixed**: `generate_maps.py` drew one 2 px marker per land cell onto a **300 px-tall** canvas — 0.83 px per 0.5° cell, markers 2.4× cell size, ~6× overplot. That exaggerated any blockiness and made real resolution unreadable. Replaced with `go.Heatmap` (one rectangle per cell, `zsmooth=False` so nothing is interpolated into fake detail) at height 460. Note the payload did **not** shrink as I expected (~57 MB vs 53.4 MB) because the z-grid carries ocean nulls the scatter omitted; the download path is the ~9.5 MB `maps_bundle.zip` regardless. Colourbar titles also moved to a `<br><sub>` line under each figure title — a long label was being laid out inside the colourbar's region, so plotly widened it and took the space from the map.

---

### 2026-07-28: Wildfire Moved to ISIMIP3b — A Mixed-Scale Member, a Wrong Annualization, and a False-FAILing QA Check

**What happened**: Re-researched wildfire from scratch (user asked for an independent review rather than trusting the catalog) and rebuilt the layer on ISIMIP3b/SSP. Four separate problems surfaced.

**1. The catalog's own heading was never honoured.** The wildfire section was titled "fire sector + biomes burntarea" but only `biomes` had ever been walked. For ISIMIP3b that was harmless — its 451 `burntarea-total` files are a **strict subset** of `biomes` (intersection 451, fire-only 0). For **ISIMIP3a it was a real gap**: the `fire` sector holds **10 models** found nowhere else (LPJ-GUESS-SPITFIRE, LPJ-GUESS-SIMFIRE-BLAZE, LPJmL5-fire, JULES-INFERNO, SSiB4-TRIFFID-Fire, ORCHIDEE-MICT, …) plus variables absent from all docs (`firesize`, `firenr`, `fireints`, `fireduration`, `fireros`). 3a is historical-only, so it cannot serve a projection layer, but it is the natural evaluation reference. **Rule**: if a catalog section names a sector, verify it was actually enumerated.

**2. A prior characterization was backwards, and it biased the shipped layer.** The catalog called `mc2-usfs` a "5–7× hotter outlier". The ratio was right; the judgement was inverted. Area-weighted against GFED4 (~3.5 Mkm² yr⁻¹): `mc2-usfs` 5.63, `visit` 4.88, `classic` 7.26 — right order; `lpj-guess` 1.09 and `lpjml` 0.94 — **~3.5× too low**. The retired 2b trio was therefore two low-biased models against one calibrated one, so its median burnt area was biased **low**. Separately, "coarse biome model" read as a grid claim and is not one: `mc2-usfs` is genuinely 0.5° (4.6% constant 2×2 blocks, clean periodicity, no seams) — the retired layer had **no** elm-eca-class spatial defect.

**3. A mixed-scale member reached a publish.** `classic`'s `2015soc-from-histsoc` burntarea files are **fraction-scaled for gfdl-esm4** (monthly mean 0.0032, max 0.32) and **percent-scaled for ukesm1-0-ll** (mean 0.31, max 35.4) — identical `units='%'` and `long_name='Burnt Area Fraction'`. That variant was chosen precisely *because* the csoil layer had pinned it for the same model, i.e. by inheritance rather than by value-check. It was caught only because the per-member load log printed one GCM at `max=1.18%` next to its sibling at `max=114%`. The build had already published `v2026-07-28_62512ee-dirty` by then; the QA gate correctly marked it **NOT verified** and the DO-NOT-USE banner fired, so nothing consumed it. **Fix**: switched to `classic`/`2015soc` (uniformly percent across both GCMs and all three SSPs, 3.58–3.74 %/yr), purged the bad version and its 6 raw objects, re-ingested with sha512 verification, re-ran. **Rule** (GUARDRAILS §9): never inherit a soc/sens variant across variables, and compare every GCM's magnitude *within* each model — a ~100× sibling gap is a unit error, not model spread.

**4. Two of my own defects, one in the layer and one in the checker.**
- **Wrong annualization would have under-scaled the layer 12×.** `visit` and `classic` publish burntarea only monthly. The csoil precedent annualizes monthly members by **mean** — correct for a *stock*, wrong here, because burnt area accumulates over its reporting interval. Settled empirically rather than by argument: `classic` publishes the same run at **daily and monthly** cadence, and each published monthly value equals the **sum** of that month's daily values (1e-6 agreement, r = 1.00000000 in all 12 months). Also surfaced that annual burnt area legitimately exceeds 100% where a cell reburns, so clamping `upper_ci` to 100 would push it **below** the median — the CI is now floored at 0 and unbounded above.
- **`generate_qa_report.py` used Pearson to test percentile direction, which false-FAILs correct layers.** `percentile` is a percentile-of-score, i.e. a deliberately **non-linear** monotone rank transform of the value. On burntarea (45% exact zeros, tail past 100%) Pearson read **+0.53** and failed all three scenarios of a correct layer. Switched to **Spearman** (rank correlation), ≈1 by construction, which still catches a genuinely inverted or scrambled percentile. This bug would equally have hit `led` and `let`. A separate FAIL — "coverage counts align with finite data" — was mine: I omitted csoil's `nmem[nmem == 0] = np.nan` masking, leaving off-land counts as a finite 0 against a NaN median.

**5. OPEN — `visit` high-latitude bias, found by rendering after every check passed.** With all four problems above fixed the layer passed QA **47/52, 0 failed**. Looking at the maps then showed bright specks on Arctic islands, which turned out to be a systematic defect: the 5 `visit` members have an **inverted zonal profile** (2090s ssp585 land-mean burnt %/yr) —

| band | visit | mc2-usfs | classic |
|---|---|---|---|
| 45–60°N | 4.76 | 0.65 | 1.15 |
| 60–70°N | 3.36 | 1.35 | 0.10 |
| 70–75°N | 5.95 | 1.45 | 0.02 |
| **>75°N** | **25.94** | 0.67 | 0.00 |

`visit` thus burns more above 75°N than in the tropics (7.3% at 23°S–0°), which is not physical; its worst cells are visit-only Arctic islands saturated near 100%/yr (81.25°N/56.75°E: 0.000% in the 2020s, 100.04% in the 2090s, single model). In the pooled layer >75°N reaches 6.5/9.0/10.0 %/yr.

**Why nothing caught it**: polar cells carry negligible area, so `visit`'s global area-weighted total (4.56 Mkm² yr⁻¹) is the **closest of the three models to GFED4's ~3.5**. Every aggregate and area-weighted statistic is blind to a defect confined to one latitude zone — the same lesson as GUARDRAILS §11, one axis over. **Fix to the process**: `generate_qa_report.py` now reports a **zonal profile** (land-mean by latitude band) for every layer, in both JSON and HTML, and warns when a polar band exceeds the tropical band for layers declaring `zonal_expectation: low_latitude_dominated`. A statistic that can see the defect now runs on every layer, not just this one.

**Status: OPEN, no decision taken.** All 12 members are retained and nothing is masked, per the user's call: *"note the potential issue with polar region anomalies, but make no decision about exclusion until the results land."* Recorded in the layer's `known_issues` attr and the catalog's `open_issues`. For reference, dropping `visit` would leave 2 models / 7 members per scenario, with `mc2-usfs` at 71% of members and 8,245 cells (12.3% of land) resting on `classic` alone at an effective 1.0° with only 2 members — which would reopen the thin-ensemble smoothing question.

**Impact**: one unverified publish (purged, never consumed); ~2.0 GB re-ingested; the QA gate and the per-member load log each caught what the other missed, and the visual review caught what both missed.

**Layer**: ISIMIP3b `biomes` `burntarea-total`, ssp126/370/585, **12 members/scenario** = `mc2-usfs` (annual, 5 GCMs, `nat/default`) + `visit` (monthly→SUM, 5 GCMs, `2015soc/default`) + `classic` (monthly→SUM, 2 GCMs, `2015soc/default`, effective **1.0°**, retained deliberately). `elm-eca` excluded (~4°×5°). Raw %, no normalization, no smoothing, shared 2020s baseline, baseline-anchored trend. Supersedes the ISIMIP2b/RCP build, whose outputs no longer exist as data anywhere.

**Rules created/updated**: GUARDRAILS §9 — declared units can't gate pooling; cross-GCM magnitude comparison within a model; cadence semantics (accumulate vs state) must be measured, ideally against a model publishing two cadences; don't clamp a cumulative quantity at its nominal ceiling; soc/sens soundness is per-variable.

---

### 2026-07-28: River Flooding — Two Products That Misreport Themselves, and a Protection Level That Decides the Answer

**What happened**: New research into severe/regional flooding, target metric "area impacted per year". Walked every round × product × sector on `files.isimip.org`, then downloaded and value-checked 10 files. Four things mattered, and three of them contradicted what the metadata or my own first reading said.

**1. There is no inundation variable in the raw water sector at all.** Confirmed across all 17 ISIMIP2b and 9 ISIMIP3b `water_global` models: no `fldfrc`, no `flddph`, no `floodedarea`. Flooded area exists **only** as `DerivedOutputData`. `maxdis` — the obvious peak-flow proxy — is in **2 of 17** 2b models and **zero** 3b models. Daily `dis` is broad but ~929 MB per 10-year chunk per member (~1 TB+ for an ensemble) and would mean re-deriving the GEV fit and hydrodynamic routing that Zimmer2023 already publishes. So the whole question reduces to which *derived* product to use — which is exactly where the metadata is least trustworthy.

**2. `floodedarea` declares an area share and is binary.** ISIMIP3b `Heinicke2026` sets `long_name="Exposed Area Share"`, `units="1"`, `definition="Flood fraction from cama flood"`. Values are strictly `{0,1}`. Verified on 6 of 15 members spanning all 3 models, 4 GCMs and all 3 SSPs. When the finding was challenged — reasonably, since `ler`, the ISIMIP2b analogue, *is* fractional — it was re-verified rather than re-asserted: byte-identical sha512 against the sidecar, `float32` with no `scale_factor`/`add_offset`, raw values read with `set_auto_maskandscale(False)` are exactly `{0.0, 1.0, 1e20}` (20,885,979 / 226,419 / 1,178,802), all 86 timesteps have exactly 2 unique values, and point series switch cleanly 0→1→0. Per cell it is an occurrence flag, so its decadal mean is flood **frequency**, not area.

**3. Binarising at 0.5° overstates flooded area 17–28×, so the binary field cannot be rescued by spatial aggregation.** Measured directly instead of argued: took Zimmer2023 `fldfrc` (true sub-grid fraction), coarsened the *same* data two ways to 0.5°, area-weighted (h08/miroc5/rcp60, 2020s decadal mean) —

| protection | true area km² yr⁻¹ | as a binary flag | inflation |
|---|---|---|---|
| none | 6,104,969 | 106,330,690 | **17.4×** |
| 100yr | 260,894 | 6,951,831 | **26.6×** |
| flopros | 743,970 | 20,693,916 | **27.8×** |

A 0.5° cell is ~2,500 km²; flagging the whole cell when any part floods counts ~71% of global land as flooded under `none` against a true ~4%. Only a sub-grid fractional product answers "how much area flooded".

**4. The protection level is the single biggest choice, and my first recommendation was wrong for the stated use case.** Same member, global flooded area, 2020s → 2090s: `none` 6.10M → 6.28M km² (**+2.9%**, saturating — it counts routine seasonal floodplain wetting that recurs every year in both decades); `100yr` 261k → 502k (**+92.4%**); `flopros` 744k → 1,130k (**+51.9%**). I initially recommended `flopros` as primary. Then the user framed the use case as *"impacts across large areas in unprotected regions"*, so I tested the tempting assumption that FLOPROS ≈ no protection where defenses are absent. **It is false**: flopros retains only 19–36% of the undefended signal even in the least-defended regions (Bangladesh 0.189, Niger/Nigeria 0.364, Mekong 0.348) and 0.8–4% in well-defended ones (Netherlands 0.008, Mississippi 0.015, Japan 0.043). FLOPROS *is* spatially real — less protective than a uniform 100yr in Niger/Bangladesh, more protective in the Mississippi/Netherlands where defenses exceed 1-in-100 — but it is not a proxy for "unprotected". Protectiveness runs `none < flopros < 100yr` globally. Recommendation changed to `none` + `100yr` as the pair, with `100yr` read as a **severity threshold** (severe events, no reliance on defense data) rather than as "protection". **All three are being built as parallel layers** at the user's request.

**5. A correction of my own framing.** I described the `fldfrc` footprint as "~76% of the grid is structurally NaN, i.e. only ~24% of land" and raised zero-filling as a decision the user had to make. That compared a **land** variable against the **whole globe including ocean**. The domain is 62,066 cells = **128.8 million km² = 95.7% of global land excluding Antarctica**, against 67,420 land cells for Lange2020 `led` (61,546 shared); 79.6% of domain cells are ≥99% inside the CaMa-Flood domain and the median is 100%. It is a normal ISIMIP land mask. There was no sparse-coverage problem and nothing needed zero-filling. **Rule**: state a coverage fraction against the denominator that makes it meaningful — for a land variable that is land, never the globe.

**Also corrected in the catalog**: Lange2020 is **twelve** variables, not six — every hazard has a `pe*` population-exposed twin (`ped/per/pew/pec/peh/pet`) that was never recorded; `lec` is published by `pepic` as well as `gepic`; `ler` is fractional (not "UNVERIFIED"); files are `.nc4` not `.nc`, which cost a false-negative listing; and the note "no ISIMIP3b/SSP version of this family exists" is true only of the *publication name* — `Heinicke2026` publishes `driedarea`, a correctly land-masked binary SSP sibling of `led`, so the drought layer **can** move to SSP. Searching a variable name across rounds misses this because the publication directory changes.

**Verified equivalence**: `ler` ≈ the `100yr` protection level of `fldfrc`, pre-aggregated to 0.5°. Matched member, 2020s decadal mean, 61,483 shared cells: mean ratio **0.99**, Spearman 0.784 — against 25.04 for `none` and 2.54 for `flopros`. So `ler` is a fixed 1-in-100 assumption with no rcp85 and no way to vary the protection.

**Layers**: `riverflood_fldfrc-{none,100yr,flopros}_annual` — ISIMIP2b `DerivedOutputData/Zimmer2023`, **rcp26/60/85**, **24 members/scenario** (6 GHMs × 4 GCMs; `clm50`/`mpi-hm`/`pcr-globwb` are historical-only). CaMa-Flood v3.6.2, GEV fit to picontrol, `discharge_threshold` 0.1 mm/d. Native **150 arcsec** coarsened area-preservingly to 0.5° at ingest (12×12 block sum, exact alignment to 1.4e-14°, area-conserving to <1e-9). Decadal **mean**, baseline-anchored trend (annual flooded area swings ~17× between adjacent decades), two-tier percentile, `higher_is_worse`, no normalization, no smoothing, CI clamped to [0,100] — safe here because flooded fraction is a bounded share, unlike burntarea's cumulative percentage.

**Deliberate deviation, recorded so it is not "corrected"**: the raw prefix holds the **0.5° coarsened** fields, not the 150 arcsec originals. The full source set is ~54 GB against 19 GB of local scratch, and raw staging is transient by contract (`STORAGE.md` — `cleanup_raw` deletes it once `source_url` + checksum are recorded). Every ingested file records the `source_url`, the **sha512 of the original**, and the exact transform, in both its own global attrs and `layer.json`. `files.isimip.org` is not behind the Anubis anti-bot that guards the `data.isimip.org` API, so re-fetching an original is routine.

**Impact**: no bad publish. One code defect caught in smoke-test (reading `t.units` after `ds.close()` → `AttributeError: NetCDF: Attribute not found`), which is why the 2-file smoke test exists before a 216-file run.

**Rules created/updated**: GUARDRAILS §9 — a challenged value-check is re-verified at the byte/dtype level (checksum, packing attrs, unscaled read), not re-asserted; §12 — report a coverage fraction against a meaningful denominator; and a new note that a *derived* product's metadata is no more trustworthy than a model's, since both `floodedarea` (binary declared as a share) and `fldfrc` (no `units`, no `long_name` at all) came from the same publication family.

---

### 2026-07-29: Flood Layers Shipped — A 60°N DEM Seam Found by Looking, and a Detector That Took Four Attempts to Calibrate

**What happened**: Built and published the three `fldfrc` protection-level layers (216 members ingested, 0 errors). QA passed and the numbers were coherent, but reviewing the contact sheet showed a faint horizontal discontinuity across northern Eurasia that no check had mentioned.

**The defect**: the zonal-mean flooded fraction **halves across a single 0.5° row at exactly 60.0°N** — 11.06% just north, 5.93% just south (1.87×), the largest single-row jump anywhere between 45° and 80°N, and identical in all six GHMs. Confirmed in the **native 150 arcsec** data, where the sharpest row step lands on exactly `60.0208 → 59.9792` in all three protection levels, so it is inherited, not introduced by the coarsening. Cause: CaMa-Flood's floodplain topography changes DEM at the **SRTM/HydroSHEDS coverage limit** (SRTM spans 60°N–56°S).

**Why nothing caught it**: the band-based zonal profile added the day before uses bands `45..60` and `60..70` — a step exactly on a band **edge** is invisible to it. Distribution statistics cannot see it either. It took a per-member image plus a row-by-row scan.

**It is a level bias, not a trend bias** — and that distinction decides whether the layer is usable. Trend and change difference the same cells against their own 2020s baseline, so the static offset largely cancels: rcp85 `none` 2020s→2090s gives 60–70°N +2.2%, 50–60°N −0.8%, 30–50°N +0.4%, tropics +8.3%, S subtropics +6.4% — coherent across the seam. Trend/change across 60°N is usable; absolute levels across it are not. Recorded in each layer's `known_issues` and `known_latitude_seams='60.0'`.

**The detector took four attempts, and each failure is worth keeping**:
1. **MAD z-score** — the *same* seam scored z = 12.2 / 11.7 / 11.7 for rcp26/60/85, straddling a threshold of 12. MAD estimated from a few hundred rows is not stable enough for a fixed cut.
2. **Ratio to the global median row-to-row change** — fired on all three *other* published layers, at −50.75°N (12–19 cells/row), 83.75°N (19–32) and 80.25°N (97–236). Fixed with a **≥150 finite cells/row** floor; the genuine seam sits on rows of 451–475.
3. **Still fired on `driedarea`** at 8–11× for a level step of only 1.25–1.33×, because that field is near-zero at high latitude so the median row variation is tiny. Fixed by additionally requiring a **≥1.5-fold level step**.
4. **Still flagged 70°N and the equator on the thresholded flood layers** as loudly as the real seam (19–49×). These are **not seams** — verified against the native data, their sharpest steps sit at 70.27° and 0.23°, i.e. the real Arctic river decline and the Congo–Amazon belt. Fixed by scoring against the **local** median (±10 rows) instead of the global one: a true seam is quiet on both sides so its ratio stays high, while a steep gradient's collapses. Final numbers: real seam **11.4–11.8×**, loudest non-seam across five other layers **6.9×**, threshold set at **9.0**.

Negative control run: the first published version, which lacks the declaration, warns on all three scenarios; the declared version passes. The check does real work rather than passing vacuously.

**A second gap I closed**: `layer.json` shipped with `inputs.files` **empty** — no `source_url`, no checksums. That is normally a nuisance; here it undercut the whole justification for coarsening at ingest, since the manifest is the *only* place the chain back to the 150 arcsec original is written down, and `cleanup_raw` refuses without it. The processor now lifts each raw file's `source_url`, original `sha512`, original byte count and transform string into the manifest (72 inputs per layer, all populated) and warns loudly if any is missing.

**Results** (24 members/scenario, global flooded area, 2020s → 2090s):

| protection | 2020s baseline | rcp26 | rcp60 | rcp85 |
|---|---|---|---|---|
| `none` | 6,157,552 km² yr⁻¹ | −0.1% | +2.2% | **+5.1%** |
| `100yr` | 332,988 km² yr⁻¹ | +12.0% | +88.9% | **+225.5%** |
| `flopros` | 872,279 km² yr⁻¹ | +1.6% | +43.7% | **+101.3%** |

Monotonic in forcing for all three, and the protection level moves the rcp85 answer from +5% to +226%. The `none` change map is **mixed** — blue (drying) across the West Siberian lowlands and parts of Europe against red elsewhere — while `100yr` is almost uniformly red. That is the saturation result made visual: routine seasonal inundation can fall while *severe* flooding rises in the same place.

**Impact**: no bad publish; three versions minted (`-dirty`, `-dirty-b`, `-dirty-c`) as the caveat and then the provenance chain were added, with the superseded ones retained as history. One processing bug caught in the 2-file smoke test before the 216-file run.

**Rules created**: GUARDRAILS **§14** — a discontinuity on a round latitude is a data-source seam; confirm against native resolution, separate seam from gradient, declare it in `known_latitude_seams`, and establish whether it biases the *level* or the *trend* before condemning the layer.

**Addendum — a proposed shortcut that was wrong, and the fourth layer it produced.** Asked whether a 1-in-100-year flood extent could be estimated as `none` + `100yr`, on the reasoning that `100yr` is the residual above the annual expectation. Tested rather than argued, and it does not hold: `100yr` is a **subset** of `none`, not its complement — `100yr > none` in **0.000%** of cells, and where `100yr > 0` it *equals* `none` in **84–88%** of them (mean ratio 0.90–0.93). The decisive check is dimensional rather than statistical: `none + 100yr` exceeds **100% of a cell** in 6,867–28,139 cells, which is impossible for a fraction of cell area. `100yr` is a *filtered copy* of `none`, kept only in years that overtopped a 1-in-100 defence; `none` is already the total undefended extent, i.e. the maximum of the three layers rather than a floor. The sum is also numerically unstable — 17% *below* the correctly measured footprint in the 2020s, drifting with scenario because the terms scale differently (+5% vs +226%). **Lesson**: when a proposed decomposition is additive, test for containment first, and prefer a *dimensional* impossibility check (a fraction exceeding 1) over a plausibility argument — the latter would have accepted it, since the sum lands within ~20% of the right answer today.

The correct construction became the fourth layer, `riverflood_fldfrc-event100yr_annual`: per member and decade, frequency = years with `100yr > 0` ÷ valid years, footprint = mean of `none` over exactly those years. **Frequency is the primary value**, because the footprint is nearly static (**+1.8%** 2020s→2090s at rcp85 — topography fixes a given flood's extent) while frequency **more than doubles (+122%)**: the preindustrial 1-in-100 flood becomes a **1-in-5.7-year** flood. Making the footprint primary would have shipped a confidently flat layer. The two change maps side by side make it unmistakable — the frequency map is deep red across every major basin, the footprint map is blank. **This retro-explains the three protection layers**: `none` saturates (+5.1%) and `100yr` explodes (+225.5%) because all four are measuring a *frequency* change; only the fourth states it directly. Caveats carried in `known_issues`: "1-in-100" is **preindustrial** (GEV fit to `picontrol`/`1860soc`, so the 2020s already sit at ~1-in-15), and the conditional mean is **selection-biased upward** (cells never exceeding in a decade drop out) → the ratio is robust, the level indicative. Published `v2026-07-29_ba7bae4-dirty`, QA PASS 52/55, 144 raw inputs from two sibling prefixes with full provenance.

**One process hazard worth not repeating**: several `until ! pgrep -f "process_fldfrc_flood"; do sleep 60; done` waiters **never terminated, because the waiter's own command line contains the pattern** — `pgrep -f` matched the waiting shell itself. Nine of them span for the rest of the session and, worse, made completed jobs look like running ones. Match on the interpreter (`pgrep -f "python.*process_x"`), use `pgrep -x`, or track the PID directly.

---

### 2026-07-28: Drought Moved to ISIMIP3b — A Product Hidden Behind a Publication Name, and a Zero-Inflation Trap That Was Really a Coverage Artefact

**What happened**: Asked to review options for tree-relevant drought, an independent re-enumeration of the ISIMIP file server found that the catalog's standing claim — *"ISIMIP3b/SSP version of this family NOT found (0 hits)"* — was wrong, and shipped a new SSP-era drought layer (`drought_driedarea_annual`) that supersedes the ISIMIP2b `led` layer.

**Issues found**:

1. **A whole product class was invisible because the search was by variable name.** The Lange 2020 exposure concept WAS re-issued for ISIMIP3b — but split across two *publication directories* with hazard-word variable names: `DerivedOutputData/Heinicke2026` (`driedarea`, `floodedarea`) and `DerivedOutputData/Zantout2025` (`heatwave`, `wildfire`, `cropfailure`). Searching `led` across rounds returns nothing, which reads as "absent". **Resolution**: enumerate `{round}/DerivedOutputData/` itself — it is a short listing (2–4 entries per round) and it is the only way to see a republished product. Added to the `isimip-search-download` skill and GUARDRAILS §8.

2. **Two other catalog assertions were also wrong.** The ISIMIP2b `forestry` sector has **10 models**, not 1 — though the conclusion survives, because every one is site-level (kroof, hyytiala, collelongo, …) and unusable globally; and the tree water-stress diagnostics (`trans`, `lai-{PFT}`, `fapar`, `pft-{tree}`, `npp-tree`) were missing from the drought variable list entirely. `trans` and `soilmoist` are **monthly-only in every model, both rounds** — there is no annual fallback.

3. **The filename field offset differs from `led` — again.** `driedarea` files are `{model}_{gcm}_w5e5_{ssp}_…`, so model/GCM are fields **[0]/[1]**; `lange2020` prefixes the publication name and puts them at [1]/[2]. This is the third time this trap has appeared (it is why `process_qg.py` could not be reused for `led`). **Resolution**: a dedicated `parse_name` plus explicit membership assertions — duplicate-member detection, expected model/GCM vocabularies, and an exact per-scenario count — because the loading pattern `dec[s][member] = maps` silently *overwrites* a duplicate rather than failing.

4. **Zero-inflation looked decisive and was measuring the wrong thing.** A single member's 2020s frequency map is **63% exact zeros**, which would trigger burntarea's automatic >2% two-tier percentile switch several times over. But on the 15-member × 3-SSP shared baseline the figure is **3.59% over the union land mask and 0.18% over fully-covered cells**. The zeros live almost entirely in the 17,431 cells (27% of the union) that not all three GHMs cover: a single-model cell rests on 30 samples against 450 for a fully covered one, so exact zeros are near-inevitable there. **The tier rule would have fired on a coverage artefact, collapsing every never-dry cell onto percentile 1 for the wrong reason.** Single-tier was applied by decision; the measured zero fraction is recorded in the output attrs. Visible consequence, left in deliberately for review: the percentile floor is 3.59 and ~3.6% of land ties there.

5. **A high cross-scenario spread did NOT mean the shared baseline was smearing signal.** A review flagged the per-member cross-SSP spread in the 2020s as 1.83× the inter-member SD, against a suggested 0.10 retain-threshold. Measuring further showed the opposite reading: the spread (mean |max−min| **0.104**) sits *below* the pure-sampling-noise floor for a 10-year mean of a Bernoulli field (`1.693·√(p(1−p)/10)` = **0.139**) and near a within-scenario 5yr-vs-5yr split (**0.079**), while the scenario global means agree to 0.005. The SSPs have barely diverged in forcing by 2020–2029; what differs is *which years happened to be dry*. Averaging across SSPs therefore suppresses weather noise rather than destroying signal. **Lesson**: for a rare-event binary field, compare a spread against its sampling-noise floor before calling it divergence.

6. **A latent mask defect found in two shipped processors.** `anchored_trend` returns `np.zeros(med_stack.shape[1:])` for the baseline decade (`process_burntarea_fire.py:236`, `process_csoil_soilcarbon.py:241`), which makes the **entire ocean a finite zero** in the baseline trend even where the median is NaN. QA does not catch it because it only checks that finite baseline trends equal zero, never that the trend and median finite-masks agree. `process_driedarea_drought.py` masks it to NaN off-land instead. **Not yet fixed in burntarea/csoil** — logged here and as a required check in GUARDRAILS §10.

7. **The Codex review harness could not start.** `codex:rescue`'s companion hard-codes a bwrap sandbox (`read-only` / `workspace-write`) and bwrap cannot create a namespace in this container, so it failed before reading a single file. **Workaround**: call the Codex CLI directly with `--sandbox danger-full-access` against a *throwaway copy* of the repo in scratchpad, so an unsandboxed agent cannot touch the live tree. The review itself was accurate — every checkable claim verified against the source, with line numbers drifting by a few.

**Impact**: none shipped wrong; all caught before or during the build. The drought layer would otherwise have been built on ISIMIP2b/RCP26+60 when an SSP product existed.

**Fix applied**: `scripts/download_driedarea_drought.py` (45 members, sha512-verified against ISIMIP sidecars) and `scripts/process_driedarea_drought.py`; catalog `drought.exposure_heinicke2026` / `exposure_zantout2025` / `forest_cover_weight` / `vegetation_impacts.tree_water_stress`; GUARDRAILS §8 and §10; the `isimip-search-download` skill.

**Rule created**: GUARDRAILS §8 (enumerate `DerivedOutputData/` publications, not just variable names across rounds) and §10 (the trend's finite mask must match the median's; an automatic tier/threshold rule must be checked against coverage before it is allowed to fire).

---

### 2026-07-28: Sitka-Spruce Timber — A Catalog Recommendation That Cannot Work, Two Value Conventions, and a Trend Manufactured by Ensemble Composition

**What happened**: Asked to find a PFT relevant to Sitka spruce speaking to general timber growth, and to check the answer against what prior research had asserted. Re-enumerated the file server from scratch. Most of the prior work held up; five recorded claims did not, and building the layers surfaced three defects that distribution statistics alone would never have shown.

**Claims refuted** (all measured, not argued):

1. **The 2026-01-21 remediation's top recommendation is not actionable.** The catalog told future readers to use MC2-USFS `mesictemperateneedleleafforest` / `maritimeevergreenneedleleafforest` for temperate-conifer productivity. MC2 publishes **no carbon or NPP resolved by biome type** — only `-total`, `-tree`, `-grass`. Its 47 biome names exist solely as `pft-{biome}` fields, and those are **one-hot presence flags valued exactly 0 or 100** despite units `%` and long_name "Plant Functional Type Grid Fraction". Worse for the maritime case: `pft-maritimeevergreenneedleleafforest` is **identically zero in every run checked** — 3b ssp370/gfdl-esm4, 3b ssp126/ukesm1-0-ll, and 3a obsclim historical (7,008,981 finite values, all 0). A remediation that fixes an incident by recording an unusable recommendation leaves the next reader worse off than no entry.
2. **ISIMIP2b forestry was undercounted tenfold.** Catalog: 1 model, 4 species, "NO oak (Quercus)". Reality: **10 models, 9 species**, including `quro` (Quercus robur) and `psme` (Douglas-fir), plus the actual forestry metrics `mai` (mean annual increment), `vol`, `harv`. ISIMIP2a adds 5 more models.
3. **CARAIB does publish conifers.** Catalog: "needleleaf PFTs defined but not in output — use for broadleaf analysis, not conifers." It publishes **seven** needleleaf classes for both `cveg` and `npp`.
4. **JULES-ES-VN6P3 was absent from the catalog entirely** — yet it is the deepest 3b source of PFT-resolved wood carbon (**5 GCMs** vs CLASSIC's 2, `cwood-ndlevg` + annual `npp-ndlevg`, and `default` transient-CO₂ runs at all three SSPs). Also, CLASSIC's entry listed 5 climate forcings; only **2** (gfdl-esm4, ukesm1-0-ll) have a `future/` directory in biomes — a real 404, not an empty listing.
5. **The "generic PFT is confounded" lesson is species-dependent, and inverts here.** It was right for loblolly (one climate zone). Sitka is hyper-oceanic and straddles zones, so the climate-zone split *fragments* its range: LPJmL temperate-NLE covers UK/Ireland (774/774 cells) but only 609/2451 in PNW-BC-SEAK; LPJ-GUESS boreal-NLE covers PNW (1866) but **6** UK cells; the generic classes cover both ends. Also worth stating precisely: a generic class is **not** an internally-weighted mixture of climate-specific conifers — `evgndltr`/`ndlevg` are each a single PFT with one global parameter set, proven by JULES splitting broadleaf evergreen into `bdlevgtemp`/`bdlevgtrop` while publishing a single `ndlevg`. The mixing therefore bites at the **percentile ranking** step (a cell ranked against a global distribution spanning boreal to subtropical), not in the raw per-cell values.

**Defects found while building** (user chose 2b for PFT specificity, accepting RCP and the loss of `cwood` — no PFT-resolved `cwood` exists in 2b at all, in `biomes` or `permafrost`):

6. **Two incompatible value conventions across models being pooled.** A PFT field is either a **per-tile density** (on the PFT's own tile area; `sum_i(frac_i * value_i)` recovers the cell total) or **per-gridcell** (already cover-scaled). Verified: orchidee / orchidee-dgvm / clm45 per-tile; lpjml / caraib per-gridcell. Pooling them raw compares different quantities — it made the apparent cross-model spread read **10.5× (cveg) and 177× (npp)**. After harmonizing and comparing on a common mask the same models agree to **2.35× and 1.83×**, at the csoil precedent, so **no normalization** was needed. The one-line test: a per-tile field *exceeds* the all-PFT total (clm45 does in 94.4% of cells); a per-gridcell field never does. Note clm45 publishes **no cover fraction at all**, so it can only ever be used per-tile.
7. **A trend manufactured entirely by ensemble composition.** The first dry run showed rcp85 losing **−0.72 kg m⁻²/dec** at the 2030s and then "recovering" — because orchidee-dgvm is in the shared cross-scenario baseline, is absent from rcp85, and sits *higher* than orchidee on the retained cells. Composition varies by scenario here (dgvm has no rcp85; clm45's GCM pair differs per RCP), so differencing an all-member baseline against per-scenario decades differences two different ensembles. Fixed by pooling each scenario's 2020s panel over **that scenario's** members, at the cost of the 2020s panel no longer being bit-identical across scenarios. **The shared-baseline convention is only safe when the ensemble is uniform across scenarios**; csoil's was, this one's is not.
8. **Level disagreement was between CELLS, not models.** The pooled mean is 6.29 vs 0.52 kg m⁻² where all vs one model reports (12×; npp 5×), because each model's periphery is marginal habitat. Pooling the union would print mask edges into the maps and distort the percentile. Also: **the ">=2 models" rule and the "treat orchidee+orchidee-dgvm as one family" rule interact** — under family counting cveg has only 2 families, so ">=2" silently becomes "all", cutting coverage from a quoted 17,217 cells to 9,698 and UK cells from 275 to 134. Two independently reasonable decisions can compose into a third one nobody chose.
9. **The `anchored_trend` ocean-zero defect (finding 6 of the 2026-07-28 drought entry) propagates by copy-paste.** Writing this layer from the csoil template reproduced it; caught here by an explicit trend-vs-median finite-mask assertion. **Still unfixed in `process_burntarea_fire.py` and `process_csoil_soilcarbon.py`.**
10. **The manifest guard fired correctly and should not have been loosened.** `utils/layer_publish.py` rejects scenario files that disagree on any non-scenario-specific attribute; `n_members` legitimately differs (8/8/6). The fix was to record a `members_by_scenario` breakdown that is byte-identical in every file — strictly more information — rather than adding `n_members` to `_SCENARIO_SPECIFIC` and weakening the check for every layer.

**Impact**: nothing shipped wrong. Findings 6–9 would each have produced a plausible-looking but wrong layer; 7 in particular would have published a fake rcp85 loss signal with a monotonic-looking recovery.

**Root cause**: the catalog is written as recommendations rather than as verified capability, so an untested suggestion (finding 1) reads exactly like a measured fact. And PFT output has more per-model degrees of freedom than any variable processed so far — convention, cover availability, mask extent, class semantics, effective grid — none of which are visible in the declared metadata.

**Fix applied**: `scripts/process_timber_tempnle.py` (both tracks, one parameterized script); catalog sections `cwood`, `pft_value_conventions`, `temperate_needleleaf_evergreen_2b`, plus corrections to `pft_codes.isimip2b_forestry`, `biomes_models` (MC2, CARAIB, CLASSIC, new JULES entry) and `pft_equivalences.temperate_needleleaf_evergreen`.

**Rule created**: GUARDRAILS §13 — establish a PFT field's value convention before pooling it; never trust a `pft-` unit label; and re-verify that the shared 2020s baseline's composition matches each scenario's before differencing.

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
