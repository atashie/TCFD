---
name: isimip-process-visualize
description: Process annualized ISIMIP NetCDF into the TCFD value-class format, then generate its QA/QC report and interactive maps. Use when processing a downloaded ISIMIP variable into a TCFD/CDP layer, reprocessing an existing layer, or regenerating QA evidence and visualizations.
---

# Process & Visualize (TCFD/CDP product)

**This skill is for the TCFD/CDP product only** — annualized decadal statistics. It is NOT
for the Water Risk Index (20 value types, monthly, standalone scripts, no `isimip-pipeline`
CLI). Never mix the two. See CLAUDE.md.

## Recommended Workflow

```bash
# 1. Search and download data
isimip-pipeline run "groundwater runoff" --name gw-runoff --keep-raw

# 2. Review QA report in browser
#    Open: ./data/processed/ outputs

# 3. Generate visualization maps
python scripts/generate_maps.py {variable} {processed_dir} {output_dir}

# 4. After verifying processed data is correct, delete raw files
isimip-pipeline cleanup ./data/raw
```

## Data Processing Parameters

- **Temporal binning**: 2010s–2090s (no data before 2010 or after 2099)
- **Adaptive windowing**: Minimum 100 data points per decade-bin
- **Percentile baseline**: Always use 2020s as reference distribution
- **Sub-annual data**: Always ask the user for the aggregation method. See GUARDRAILS.md §2.
  Do not pick silently — note the units first (a density metric wants mean, an accumulating
  one wants sum), present mean/median/sum/min/max, and use `AskUserQuestion`.

## Shared 2020s Baseline

All climate projections share **identical 2020s baseline values**, computed as the average
of ALL scenarios (projection + historical/picontrol) that have overlapping 2020s data. For
2030s–2090s, values are computed per-scenario as before. The `baseline_source` attribute in
processed NetCDF files records whether the shared baseline was used.

Implementation pattern:
1. Collect 2020s data from ALL available scenarios (projection + historical/picontrol)
2. Average across scenarios per GCM to create the shared baseline
3. Use the shared baseline for ALL scenarios' 2020s decade
4. Use scenario-specific data for 2030s–2090s
5. Set `baseline_source: "shared_across_all_scenarios"`
6. Output separate files per scenario

**Verification**: Run `python scripts/test_shared_baseline.py {output_dir}` after processing
to verify:
1. 2020s values are identical across all scenarios
2. 2030s+ values differ across scenarios (as expected)
3. Files carry `baseline_source: shared_across_all_scenarios`

> This verification step was dropped during the S3 era and no layer built since carries it.
> It is restored here as mandatory. `generate_qa_report.py` also checks baseline identity;
> run both — the script is the contract, the QA check is the safety net.

**Caveat added since (additive, keeps the rule):** the shared baseline is correct **only
when ensemble composition is uniform across scenarios.** Check first. If any member is
missing from any scenario, pool each scenario's baseline panel over *that scenario's*
members, or the change map differences two different ensembles and manufactures trend:
`timber_cveg-tempnle` read **−0.72 kg m⁻² dec⁻¹** at rcp85 2030s with a plausible later
"recovery", purely because `orchidee-dgvm` sat in the baseline and is absent from rcp85. The
2020s panel is then no longer bit-identical across scenarios — declare `members_by_scenario`
(member **identity**, not a count) so QA groups by composition instead of failing. Keep the
percentile *reference* distribution global so percentiles stay comparable.

**Reference implementations**: `scripts/process_qg.py` (groundwater runoff),
`scripts/process_timber.py` (timber/wood carbon).

## Output Organization

**Folder naming**: `{descriptive-name}_{variable}_{timestep}/`
- Example: `wildfire-burntarea_burntarea_monthly/`
- Example: `drought-severity_led_annual/`

**File naming**: `{variable}_{scenario}_processed.nc` — one file per scenario, dims
`(decade, lat, lon)`. Never a per-decade file.

> Note: the original guidance said all scenarios belong in ONE file with a `scenario`
> dimension. That was superseded by your own pre-AWS decision (`a34bafc`, and the rcp45
> incident fix which globs `{variable}_*_processed.nc`), not by the S3 work — so the
> per-scenario form stands. All scenarios still belong in ONE folder; never split folders
> per scenario (avoid `burntarea-rcp26/`, `burntarea-rcp60/`).

**CLI workflow**: Always use `isimip-pipeline run` for downloading data. It handles
multi-scenario downloads into a single output directory. Avoid manual `search` + `download`
workflows that may fragment scenarios into separate folders.

## Non-negotiables

Read [GUARDRAILS.md](../../../GUARDRAILS.md) before choosing any statistic.

- **§9 — value-check before choosing statistics.** Print min/max, median, exact-0 and
  exact-1 counts, unique-value count, units, `long_name`, time units AND calendar **per
  member**. Never infer a variable's nature from its name, its CF metadata, or a sibling
  variable. Metadata is wrong often enough that values are the only authority.
- **§9 — verify soc / CO₂ experiment tokens per model.** A uniform treatment across an
  ensemble is not guaranteed. Do not assume which *direction* a fixed-CO₂ member biases the
  trend — measure it. For `csoil`, fixed CO₂ produced the *strongest loss*, the opposite of
  the documented expectation.
- **§2 — sub-annual data: ask the user for the aggregation method.**

## Value-check every member (§9)

Print the per-member table before writing any processing logic. Record verified findings in
`config/isimip_search_catalog.yaml` under `data_nature`, and carry them into the output's
global attributes.

**Then RENDER a per-member contact sheet and look at it.** One small global panel per
member, before choosing statistics. This is not optional polish: every statistic in the
table is *invariant under spatial rearrangement*, so the table cannot see a spatial defect.
A `~4°×5°` member once passed the full table twice and 37 algebraic QA checks, and a user
caught it by looking. Look for block structure, seams/banding, hard unrealistic edges,
land-mask errors, hemisphere flips, and patches unrelated to geography. Check per member —
the pooled ensemble dilutes one bad member.

Traps that have actually occurred:

| Trap | Real example |
|---|---|
| Wrong `long_name` | `burntarea` lpj-guess labels burnt-area % as "Fire Return Interval" |
| Divergent time units | `burntarea` mc2-usfs uses `days since`, siblings use `years since` |
| Divergent calendars | `csoil` jules is `proleptic_gregorian`, the other four `365_day` |
| Missing `calendar` | two ISIMIP2b `csoil` models omit it entirely |
| Same variable, different name per round | ISIMIP2b `csoil` vs ISIMIP3b `csoil-total` |
| Byte-identical cross-sector duplicate | `elm-eca` csoil-total under both `biomes` and `permafrost` — would double-weight the model |
| Heterogeneous land masks | `csoil`: 58,714–67,647 cells across 5 models |
| **Declared grid ≠ effective grid** | `csoil` classic declares 0.5°/360×720 but is natively **1°**, replicated 2×2 with a one-cell longitude offset |
| **Declared unit ≠ actual scale** | `burntarea` clm45/orchidee declare `%` on a 0–1 **fraction** scale (~1000× low) |
| **Mis-scaling *within* one model, across GCMs** | `burntarea` classic `2015soc-from-histsoc`: gfdl-esm4 fraction-scaled, ukesm1-0-ll percent — identical `units` and `long_name`. A ~100× sibling gap is a unit error, not model spread |
| **soc/sens variant sound for one variable, broken for another** | `classic`/`2015soc-from-histsoc` is fine for `csoil`, mixed-scale for `burntarea`. Never inherit a variant |
| **Monthly cadence semantics assumed** | `burntarea` **accumulates** → annual = **SUM**; `csoil` is a stock → annual = **mean**. Copying the csoil precedent under-scales fire 12× |
| **Clamping a cumulative quantity at its nominal ceiling** | annual burnt area exceeds 100% where a cell reburns; clamping `upper_ci` to 100 drives it below the median |
| **PFT field reported per-tile vs per-gridcell** | `timber_*-tempnle`: `orchidee`/`clm45` report on the PFT's own tile area, `lpjml`/`caraib` are cover-scaled. Pooling raw made the spread read **10.5×/177×**; harmonized it is **2.35×/1.83×** |
| **`pft-` cover-fraction units lie** | `classic`/`orchidee` declare `%` but store **0–1**; `jules`/`lpjml`/`caraib` store true percent |
| **`long_name` copied from an unrelated variable** | `orchidee` `cveg-tendev` reads "crop biomass yield"; its `npp` claims "positive" yet holds negatives; `caraib` has **no `long_name` at all** |
| **PFT class name contradicts its geography** | `caraib` `ndevtecdt` has "te" in the code but peaks at 55–70°N — it is **boreal** |
| **Zero-inflation faking a coarse grid** | PFT fields are 50–72% exact zeros, so an all-zero 2×2 block counts as "constant". Run resolution tests on **strictly positive** blocks only |
| **A DERIVED product's metadata is no safer than a model's** | ISIMIP3b `Heinicke2026` `floodedarea` declares `units="1"` and is strictly **binary {0,1}** — an occurrence flag, so its decadal mean is *frequency*, not area. Do not skip §9 because a file sits under `DerivedOutputData` |
| **A mask defect that is per-VARIABLE, not per-product** | the same `floodedarea` is non-NaN over **94.7% of the globe including ocean**, while `driedarea` — same product, same model — is correctly land-masked |
| **A sharp step on a round latitude** | `fldfrc` halves across one 0.5° row at exactly **60.0°N** (1.87×) because CaMa-Flood changes DEM at the SRTM/HydroSHEDS limit. Confirm against the **native** data, then declare it in `known_latitude_seams` |
| **An additive decomposition that is actually containment** | "1-in-100 extent = `none` + `100yr`" fails: `100yr` is a **subset** of `none`. Caught by a **dimensional** check — the sum exceeds 100% of a cell — not by plausibility |

**Checking `ds.sizes` proves nothing about resolution.** Test the *values*: exact-tie
fraction between adjacent cells **at both offsets** (an aligned-only 2×2 test misses an
offset grid), or the inside-vs-seam gradient ratio per candidate block width. A
variance-loss-under-coarsening test does **not** work.

**Decode time with `cftime`**, not days-per-year arithmetic. Dividing by 365 puts December
of a monthly member at ~year+0.96, and rounding pushes it into the next year.

## Record the load-bearing choices

Each goes in the output's global attrs:

- `statistic` — median for continuous variables; **mean** for the Lange 2020 exposure family
  and for smooth stocks.
- `normalization` — `none` when members share a unit and comparable magnitudes ("model
  democracy"); robust z-score only when scales genuinely differ. **Compare the statistic you
  actually pool with**: medians can agree within 1.8× while means differ 2.6×.
- `spatial_smoothing` — 5×5 exponential-decay for thin ensembles; `none` for thick ones.
- `percentile_direction` — `higher_is_worse` (hazards) or `higher_is_better` (assets like
  stored carbon, where the risk is *loss* and the percentile is **inverted**).
- `baseline_decade`, `baseline_source` — see Shared 2020s Baseline above.
- Zero-inflated hazards → two-tier percentile (zeros → 1; non-zeros ranked against the
  non-zero baseline → [2,100]).

Emit `n_members` / `n_models` per cell whenever the ensemble's land masks differ, so the CI
is auditable. Do not silently mask thin cells — that is a product decision for the user.

**Before differencing anything, assert `isfinite(trend) == isfinite(median)` per decade.** A
bare `np.zeros()` for the baseline decade makes the entire **ocean** a finite zero, and QA
does not catch it. This defect is live in `process_burntarea_fire.py` and
`process_csoil_soilcarbon.py` and **propagates by copy-paste**.

**When members are two configurations of one model, pool by FAMILY, not by member.**
`orchidee` and `orchidee-dgvm` are the same code with/without dynamic vegetation, so a flat
member mean gave one model 6 of 8 votes. Mean within family, then across families; let
`n_models` count families. **Then re-derive any mask rule** — "≥2 models" plus family
counting silently became "≥2 *families*", and coverage fell from 17,217 cells to 9,698. Two
sound decisions can compose into a third nobody chose.

**A level step between multi-model and single-model cells is not model disagreement.** Where
model masks barely overlap, the pooled mean was 6.29 vs 0.52 kg m⁻² (**12×**) because each
model's periphery is marginal habitat. Decide a minimum-model rule with the user and state
coverage in cells for the regions the product is actually for.

## QA Report Generation

**Always use `scripts/generate_maps.py`.** Never write ad-hoc inline scripts for report
generation. The established script provides per-scenario HTML files (one per metric ×
scenario), 2020s vs 2090s comparison structure, cross-navigation, a master index with grid
layout, browser-safe file sizes (~4 MB each), anomaly detection and JSON summaries, and
auto-detection of SSP vs RCP scenarios.

```bash
python scripts/generate_maps.py {variable} {processed_dir} {output_dir}

# Example for heatwave data:
python scripts/generate_maps.py leh ./outputs/heatwave-exposure_leh-annual/processed ./reports/maps
```

### Dashboard conventions (2026-08-10 — apply to EVERY layer's dashboard)

These are settled decisions, not per-layer taste. Do not reintroduce the old shape.

**Tabs are exactly:** `Median | Percentile | Trend | Confidence | Anomaly | Members`
(`DASHBOARD_TABS`).

- **No `Lower_Ci` / `Upper_Ci` tabs** — the Confidence tab already plots both bounds side
  by side; separate tabs were pure duplication.
- **No `Ols_Slope` / `Sen_Slope` / `Change` tabs.** All three are decadal-change views and
  belong on the single **Trend** tab, so the two estimators can be read against each other
  — they fail in OPPOSITE regimes, so their *disagreement* is the diagnostic (OUTPUT-SPEC).

**Trend tab shows the 2030s and 2090s, never the 2020s** (`TREND_DECADES`). The baseline
panel's slopes are NaN by contract and its change-from-itself is identically 0, so a 2020s
panel is guaranteed blank. It also carries the change maps: 2030s − 2020s and 2090s − 2020s.

**Diverging panels are centred on zero** with equal blue/red extent: `cmin = -L`,
`cmax = +L`, where `L = SYMMETRIC_LIMIT_PERCENTILE`-th percentile of `|value|`
(**currently 95**; `None` would mean true `max|value|`). Cells beyond `L` are **clamped to
the endpoint colour by Plotly, never blanked**, and each panel's subtitle states its scale
and clamped share.

Do not go back to true max. Measured on the wildfire layer:

| | true max | 95th pct |
|---|---|---|
| `ols_slope` limit | 38.2 | **1.16** |
| cells using >10% of the colour range | 0.33% | **37.5%** |
| clamped at the ends | 0% | 5.0% |

A handful of outliers set `max|value|` (median `ols_slope` is 0.049), so the true-max scale
left ~99.7% of cells in the near-white middle.

**Slope panels do NOT share a scale.** `sen_slope`'s limit came out 14× smaller than
`ols_slope`'s (0.080 vs 1.16) because Sen collapses to zero on zero-inflated fields. Each
metric gets its own symmetric limit — a shared one would render Sen entirely white — so
equal redness across those two rows does **not** mean an equal rate. Say so in any writeup.

**No colorbar titles.** A title like "Annual burnt area (percent of grid cell burned per
year) [%]" reserves more horizontal space than the map. That text is now the figure title
directly above the map, with the panel/decade label as a smaller second line.

**Hover formats are per metric** (`HOVER_FORMATS`). Percentiles are integers on [1,100] —
`.0f`, not `.3e`; scientific notation there is noise, not precision. Counts (`n_members`,
`n_models`) likewise.

**Members tab** renders every LSM × GCM member on ONE shared colour scale, plus a table of
each member's land mean / median / P95 / max / zero% / cell count — so a member running
systematically hot, or with an unlike spatial distribution, is visible at a glance. It reads
an optional **`{variable}_members.nc`** sitting beside the processed files, dims
`(member, lat, lon)`, variable `value`, attr `member_field`. **Every processor should emit
it** (see `process_burntarea_isimip3b.py`, which also offers `--members-only` to rebuild it
from cache in seconds without redoing the slopes). It is a *diagnostic*, not part of the
OUTPUT-SPEC contract; when the file is absent the tab is skipped rather than failing, so
older layers still render.

**Reversed diverging scales.** `RdBu` runs red→blue, so a layer where increase = harm needs
`reversescale=True` to put red at the top. This now applies to `ols_slope`/`sen_slope`/
`change`, not just the legacy `trend` (previously only `trend` got it, so slope maps on
`higher_is_worse` layers rendered increases in *blue*).

### Browser payload — the binding constraint is MARKER COUNT, not file size

Every map is an SVG `Scattergeo`: **one DOM marker per land cell** (~70,849 at 0.5°). SVG
scatter degrades well before 100k nodes, so panels-per-page is what makes a dashboard slow.
The wildfire dashboard measured, before → after the 2026-08-10 reductions:

| Page | Panels | Markers | Before | After |
|---|---|---|---|---|
| Members | 22 | 1.56M → 390k | **37.6 MB** | 6.2 MB |
| Trend | 6 | 425k | 11.2 MB | 8.0 MB |
| Confidence | 4 | 283k | 7.3 MB | 5.1 MB |
| Median / Percentile | 2 | 142k | 3.8 MB | 2.6 MB |
| **Layer total** | | | **127 MB** | **69 MB** |

Three knobs, none of which change what the reader sees:

- `COORD_DECIMALS = 2` — exact for a 0.5° grid, zero loss.
- `VALUE_SIGFIGS = 4` — finer than a colour can encode.
- `MEMBERS_GRID_STRIDE = 2` — block-averages **the Members tab only** to 1°. Use
  `block_mean()`, **never** `values[::2, ::2]`: slicing samples every other cell and
  silently deletes burning cells on a sparse hazard.

**Watch the count when adding panels.** A tab with many panels (Members) is the one that
breaks first. If pages still feel slow, the real fix is switching from SVG markers to a
raster `go.Heatmap` — equirectangular is exactly linear in lon/lat so it maps 1:1 — at the
cost of the geographic coastline overlay. Not done yet; it is a larger change.

**Do NOT**:
- Write custom inline Python to generate HTML visualizations
- Create per-decade files instead of per-scenario files
- Generate monolithic single-file reports (>10 MB)
- Merge the collection into one self-contained HTML — plotly.js comes from the CDN so there
  is nothing to de-duplicate, and it would be ~57 MB with every figure instantiated at once

**Review the contact sheet FIRST** — per-member, full 0.5°. The pooled maps cannot show a
defect confined to one member.

Map values serialize at 5 significant figures, which is display-only — full precision stays
in the NetCDF. Do not "fix" this by writing full float64.

## Review the QA report before claiming success

`generate_qa_report.py` checks: value classes present and shaped; `lower_ci ≤ median ≤
upper_ci`; zero-width CIs isolated to all-zero or single-model cells; percentile in [1,100]
and oriented to match the declared direction; shared baseline bit-identical across
scenarios; coverage counts consistent; land coverage non-empty.

Read the warnings, don't just check the verdict. If a check reports itself **skipped**,
treat that as a failure to investigate — a silently skipped invariant is worse than a failed
one.

**When a check fails because it encodes an assumption your layer legitimately breaks, fix
the CHECK — never the data, and never by loosening the check for every layer.** Group on
member **identity**, not count. If a regrouping leaves no group with two scenarios, emit an
explicit `NOT TESTED` warning — a check that quietly tests nothing is worse than one that
fails.

**A green verdict is not verification.** These are algebraic self-consistency checks; a field
can satisfy all of them and still be geophysically wrong. **View the maps** and confirm the
geography is plausible — mountain ranges, biome boundaries and coastlines where they belong.
Then say plainly whether you looked: if you have not viewed an image, report the layer as
*unreviewed*, not *verified*. When you find a visual defect, add an automated check for its
class so the next layer fails loudly.

## Only then clean up raw

`isimip-pipeline cleanup ./data/raw` — never before the user has reviewed the maps, and
never before every input's `source_url` and checksum are recorded.

## The output contract

**[OUTPUT-SPEC.md](../../../OUTPUT-SPEC.md) is authoritative.** Do not restate or
re-derive it — read it, then call the shared implementation:

```python
from scripts.utils.decadal_stats import (
    is_boolean_field, pooled_decadal_stat, expanding_slopes)
```

Never hand-roll these statistics in a processor. The two divergent families that grew up
here (mean-then-median vs median-then-mean; IQR vs ±1 SD; three different trend
definitions) are exactly what the shared module exists to prevent.

`scripts/process_csoil_soilcarbon.py` is the **reference implementation** — copy its
structure. Note especially that it retains **annual** member series rather than collapsing
to decadal means on load (the expanding-window slopes need them), and repacks to land
cells to keep that affordable.

Three things it does that are easy to omit and must not be:

- **Classify the field from its VALUES** with `is_boolean_field`, and record the branch in
  `decadal_statistic` / `field_nature`. §9 — never from the name.
- **Assert the slope and median masks agree** after assembly. A bare `np.zeros()` for the
  baseline panel makes the entire ocean a finite zero, and the QA report does *not* catch
  it (it only checks that *finite* baseline trends equal zero). The reference
  implementation fails loudly instead.
- **Write compressed** (`zlib`, `complevel 4`) — ~7× smaller, nothing downstream changes.

### Verify before claiming success

```bash
python scripts/test_shared_baseline.py {processed_dir}
```

Exits non-zero on any contract violation. It checks the shared 2020s baseline is
bit-identical across scenarios *and* that later decades diverge, CI ordering, percentile
bounds and orientation, that baseline slopes are NaN rather than 0, that no slope is
finite where `median` is NaN, and ensemble-depth sanity. It prints `[SKIP] … NOT TESTED`
rather than silently passing when only one scenario exists.

It also reports `sen exactly zero on X% of cells` — when that is high the layer is
zero-inflated and `ols_slope` is the field to read.

## Tooling status

`utils/contact_sheet.py`, `utils/trend_significance.py`, `utils/finalize.py` and
`utils/layer_publish.py` exist only on `origin/main` and are S3-coupled. The *lessons*
above apply now; the *tooling* is pending port. Do not reference these modules until they
land locally. (`trend_significance.py` is superseded by the dual-slope contract.)
