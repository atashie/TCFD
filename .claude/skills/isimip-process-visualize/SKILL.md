---
name: isimip-process-visualize
description: Process annualized ISIMIP NetCDF into the TCFD 8-value-class format, then publish the layer to S3 with its QA/QC report and interactive maps. Use when processing a downloaded ISIMIP variable into a TCFD/CDP layer, reprocessing an existing layer, or regenerating QA evidence and visualizations.
---

# Process & Visualize (TCFD/CDP product)

**This skill is for the TCFD/CDP product only** — 8 value classes, annualized decadal
statistics. It is NOT for the Water Risk Index (20 value types, monthly, standalone
scripts, no `isimip-pipeline` CLI). Never mix the two. See CLAUDE.md.

## Non-negotiables

Read [GUARDRAILS.md](../../../GUARDRAILS.md) before choosing any statistic. The rules
that bite most often here:

- **§9 — value-check before choosing statistics.** Print min/max, median, exact-0 and
  exact-1 counts, unique-value count, units, `long_name`, time units AND calendar
  **per member**. Never infer a variable's nature from its name, its CF metadata, or a
  sibling variable. Metadata is wrong often enough that values are the only authority.
- **§9 — verify soc / CO₂ experiment tokens per model.** A uniform treatment across an
  ensemble is not guaranteed. And do not assume which *direction* a fixed-CO₂ member
  biases the trend — measure it. For `csoil`, fixed CO₂ turned out to produce the
  *strongest loss*, the opposite of the documented expectation.
- **§10 — `trend` is a THEIL-SEN SLOPE of the DECADAL MEDIAN series**, expanding window
  from the baseline decade, in value per decade. Use
  `trend_significance.theilsen_decadal(med_stack, DECADES, ...)` — the processor already
  holds the decadal stack, so no extra data is needed, and it needs **no rescaling**
  (fitting against the decade index already gives per-decade; a `× window_years` would
  inflate every trend 10×). **Never fit the ANNUAL series**: Theil-Sen is a *median* of
  pairwise slopes, so on a zero-inflated hazard it returns exactly 0 wherever most
  year-pairs are 0→0 — **91.3% of `driedarea` ssp126 cells**, with 25.1% of ssp585 cells
  pairing `p<0.05` with a zero slope. The decadal series cuts that to 10–14%: a
  reduction, **not a cure**, so QA warns when a zero slope meets a significant p-value.
  Never fit within a single decade, and never use the superseded two-point rate. The
  baseline panel is **NaN, not 0**, and `trend × elapsed_decades == change map` **no
  longer holds**. Declare `trend_method: theil_sen_on_decadal_median_series`.
- **§15 — trend significance is mandatory, and is computed on the ensemble MEAN annual
  series.** Emit `trend_pvalue` / `trend_tau` / `trend_n_obs` via
  `scripts/utils/trend_significance.py`. Never stack member-years as independent
  observations: between-model level offsets then dominate the sample (68.7× the
  interannual SD on `csoil-total`), collapsing detection from 83.0% of land to 14.4%.
  Never hand-roll the test. A constant series gives **p=1.0, not NaN**, and the p-value
  measures monotonicity of the ensemble mean, **not** inter-model agreement.
- **§2 — sub-annual data: ask the user for the aggregation method.** Do not pick silently.

## Workflow

### 1. Stage raw from S3

Raw lives at `TCFD/raw/isimip/{layer_id}/`. Processors call
`storage.stage_raw(LAYER_ID, RAW_PATTERN)`; if it returns nothing, the members have not
been ingested yet — see the `isimip-search-download` skill. Never hardcode an S3 key;
every key comes from `isimip_pipeline/storage.py`.

Make `RAW_PATTERN` cadence-agnostic when a variable is published at more than one
cadence (`*_{var}_global_*_{y0}_{y1}.nc` rather than `..._annual_...`), or you will
silently drop members.

### 2. Value-check every member (§9)

Print the per-member table before writing any processing logic. Record the verified
findings in `config/isimip_search_catalog.yaml` under `data_nature`, and carry them into
the output's global attributes.

**Then RENDER a per-member contact sheet and look at it (§11).** Use
`utils/contact_sheet.py` → `render_contact_sheet({member: 2D array}, ...)` from inside the
processor, where each member still exists separately, and hand the path to
`finalize_layer(..., extra_maps=[sheet])` so it is linked from the map index and
bundled. One small global panel per member, before choosing statistics. This is not optional polish: every statistic in
the table above is *invariant under spatial rearrangement*, so the table cannot see a
spatial defect. A `~4°×5°` member once passed the full table twice and 37 algebraic QA
checks, and a user caught it by looking. Look for block structure, seams/banding, hard
unrealistic edges, land-mask errors, hemisphere flips, and patches unrelated to geography.
Check per member — the pooled ensemble dilutes one bad member.

Watch for these, all of which have actually occurred:

| Trap | Real example |
|---|---|
| Wrong `long_name` | `burntarea` lpj-guess labels burnt-area % as "Fire Return Interval" |
| Divergent time units | `burntarea` mc2-usfs uses `days since`, siblings use `years since` |
| Divergent calendars | `csoil` jules is `proleptic_gregorian`, the other four `365_day` |
| Missing `calendar` | two ISIMIP2b `csoil` models omit it entirely |
| Same variable, different name per round | ISIMIP2b `csoil` vs ISIMIP3b `csoil-total` |
| Byte-identical cross-sector duplicate | `elm-eca` csoil-total under both `biomes` and `permafrost` — would double-weight the model |
| Heterogeneous land masks | `csoil`: 58,714–67,647 cells across 5 models |
| **Declared grid ≠ effective grid** | `csoil` classic declares 0.5° / 360×720 but is natively **1°**, replicated 2×2 with a one-cell longitude offset |
| **Declared unit ≠ actual scale** | `burntarea` clm45/orchidee declare `%` on a 0–1 **fraction** scale (~1000× low) |
| **Mis-scaling *within* one model, across GCMs** | `burntarea` classic `2015soc-from-histsoc`: gfdl-esm4 fraction-scaled, ukesm1-0-ll percent — identical `units` and `long_name`. Compare every GCM's magnitude inside each model; a ~100× sibling gap is a unit error, not model spread |
| **soc/sens variant sound for one variable, broken for another** | `classic`/`2015soc-from-histsoc` is fine for `csoil`, mixed-scale for `burntarea`. Never inherit a variant — re-value-check it |
| **Monthly cadence semantics assumed** | `burntarea` **accumulates** → annual = **SUM**; `csoil` is a stock → annual = **mean**. Copying the csoil precedent under-scales fire 12×. If a model publishes two cadences (classic: daily + monthly), use them as ground truth |
| **Clamping a cumulative quantity at its nominal ceiling** | annual burnt area exceeds 100% where a cell reburns; clamping `upper_ci` to 100 drives it below the median |
| **PFT field reported per-tile vs per-gridcell** | `timber_*-tempnle`: `orchidee`/`clm45` report on the PFT's own tile area, `lpjml`/`caraib` are already cover-scaled. Pooling raw made the spread read **10.5×/177×**; harmonized on a common mask it is **2.35×/1.83×** |
| **`pft-` cover-fraction units lie** | `classic` and `orchidee` declare `%` but store **0–1**; `jules`/`lpjml`/`caraib` store true percent |
| **`long_name` copied from an unrelated variable** | `orchidee` `cveg-tendev` reads "crop biomass yield" (an `ncrename` artifact in the file's own `history`); its `npp` claims "positive" yet holds negatives; `caraib` has **no `long_name` at all** |
| **PFT class name contradicts its geography** | `caraib` `ndevtecdt` has "te" in the code but peaks at 55–70°N — it is **boreal**. With no `long_name`, classes must be identified by cover-weighted latitude profile |
| **A "temperate"/"boreal" class is not a sub-mixture** | `evgndltr`/`ndlevg` are each ONE PFT with one global parameter set; only cover fraction and climate forcing vary per cell. JULES splits broadleaf evergreen into `bdlevgtemp`/`bdlevgtrop` but publishes a single `ndlevg` — proof the conifer class is not internally climate-split |
| **Zero-inflation faking a coarse grid** | these PFT fields are 50–72% exact zeros, so an all-zero 2×2 block counts as "constant": the resolution test read 49–66% for native-0.5° models. Run it on **strictly positive** blocks only |
| **A DERIVED product's metadata is no safer than a model's** | ISIMIP3b `Heinicke2026` `floodedarea` declares `long_name="Exposed Area Share"`, `units="1"`, `definition="Flood fraction from cama flood"` and is strictly **binary {0,1}** — an occurrence flag, so its decadal mean is *frequency*, not area. Its sibling `fldfrc` carries **no `units` and no `long_name` at all**. Do not skip §9 because a file sits under `DerivedOutputData` |
| **A mask defect that is per-VARIABLE, not per-product** | the same `Heinicke2026` `floodedarea` is non-NaN over **94.7% of the globe including ocean** (mid-Pacific reads `0.0`, Sahara reads `NaN`), while `driedarea` — same product, same model — is correctly land-masked. Check the mask of every variable you use |
| **A sharp step on a round latitude** | `fldfrc` halves across one 0.5° row at exactly **60.0°N** (1.87×) because CaMa-Flood changes DEM at the SRTM/HydroSHEDS limit. Confirm against the **native** data (a seam sits on the round latitude and one row; a gradient spreads over many), then declare it in `known_latitude_seams`. See GUARDRAILS §14 |
| **An additive decomposition that is actually containment** | "1-in-100 extent = `none` + `100yr`" fails: `100yr` is a **subset** of `none` (`>` in 0.000% of cells; *equal* in 84–88% where active). Caught by a **dimensional** check — the sum exceeds 100% of a cell — not by plausibility, since it lands within ~20% of the right answer. Test containment before summing |

**Checking `ds.sizes` proves nothing about resolution.** A natively coarse model
replicated onto the ISIMIP grid reports the same 360×720 as a native 0.5° model. Test the
*values*: exact-tie fraction between adjacent cells **at both offsets** (an aligned-only
2×2 test misses an offset grid — this is what let `classic` through the first review), or
the inside-vs-seam gradient ratio per candidate block width. A
variance-loss-under-coarsening test does **not** work — it cannot separate blockiness from
ordinary geophysical smoothness. `generate_qa_report.py` now runs this automatically.

**Decode time with `cftime`**, not days-per-year arithmetic. Dividing by 365 puts
December of a monthly member at ~year+0.96, and rounding pushes it into the next year —
misassigning one month in twelve.

### 3. Decide and record the load-bearing choices

Each of these goes in the output's global attrs, and `utils/layer_publish.py` lifts them
into `layer.json` so the manifest cannot drift from the data:

- `statistic` — the decadal statistic. Median for continuous variables; **mean** for the
  Lange 2020 exposure family and for smooth stocks.
- `normalization` — `none` when members share a unit and comparable magnitudes ("model
  democracy"); robust z-score only when scales genuinely differ (water-index TWS).
  **Compare the statistic you actually pool with**: medians can agree within 1.8× while
  the means differ 2.6× because of tail behaviour.
- `spatial_smoothing` — 5×5 exponential-decay for thin ensembles (e.g. 1 model × 4 GCMs);
  `none` for thick ones.
- `percentile_direction` — `higher_is_worse` (hazards) or `higher_is_better` (assets like
  stored carbon, where the risk is *loss* and the percentile is **inverted**).
- `trend_definition` + `trend_units` — see §10.
- `significance_method` + `significance_definition` + `significance_pooling` — see §15.
  `generate_qa_report.py` gates its significance checks on `significance_method` being
  present, so omitting it makes the layer WARN as "significance not carried" rather
  than fail — do not rely on that; emit the fields.
- `baseline_decade`, `baseline_source` — the shared 2020s baseline must be **identical
  across scenarios**, computed from all scenarios with overlapping 2020s data — **but only
  when ensemble composition is uniform across scenarios.** Check first. If any member is
  missing from any scenario, pool each scenario's baseline panel over **that scenario's**
  members, or the change map differences two different ensembles and manufactures trend:
  `timber_cveg-tempnle` read **−0.72 kg m⁻² dec⁻¹** at rcp85 2030s with a plausible later
  "recovery", purely because `orchidee-dgvm` sat in the baseline, is absent from rcp85, and
  is *higher* than `orchidee` on the retained cells. The 2020s panel is then no longer
  bit-identical across scenarios — declare `members_by_scenario` (member **identity**, not
  a count) so the QA check groups by composition instead of failing. Keep the percentile
  *reference* distribution global so percentiles stay comparable.
- Zero-inflated hazards → two-tier percentile (zeros → 1; non-zeros ranked against the
  non-zero baseline → [2,100]).

Emit `n_members` / `n_models` per cell whenever the ensemble's land masks differ, so the
CI is auditable. Do not silently mask thin cells — that is a product decision for the user.

**Before differencing anything, assert `isfinite(trend) == isfinite(median)` per decade.**
A bare `np.zeros()` for the baseline decade makes the entire **ocean** a finite zero, and
QA does not catch it (it only checks that finite baseline trends *equal* zero, never that
the masks agree). This defect is live in `process_burntarea_fire.py` and
`process_csoil_soilcarbon.py` and **propagates by copy-paste** to anything templated on
them.

**When members are two configurations of one model, pool by FAMILY, not by member.**
`orchidee` and `orchidee-dgvm` are the same code with/without dynamic vegetation (the
non-DGVM file's own attrs read `ORCHIDEE-MICT(w/tDGVM)` with a `nodgvm` path), so a flat
member mean gave one model 6 of 8 votes. Mean within family, then across families; let
`n_models` count families. **Then re-derive any mask rule** — "≥2 models" plus family
counting silently became "≥2 *families*", which for a 2-family track is "all models", and
coverage fell from a quoted 17,217 cells to 9,698. Two sound decisions can compose into a
third nobody chose; recompute the numbers after both are applied and re-quote them.

**A level step between multi-model and single-model cells is not model disagreement.**
Where model masks barely overlap, the pooled mean was 6.29 vs 0.52 kg m⁻² on all-model vs
one-model cells (**12×**) because each model's periphery is marginal habitat. Pooling the
union prints hard mask edges into the maps *and* distorts the percentile — a periphery cell
ranks low for having fewer contributing models, not for a low value. Decide a minimum-model
rule with the user and state coverage in cells for the regions the product is actually for.

### 4. Publish, then finalize — one call, not three

```python
version = publish_processed_layer(LAYER_ID, stage, created_by=..., notes=...)
finalize_layer(LAYER_ID, version=version)   # QA report + maps
```

`finalize_layer` (`scripts/utils/finalize.py`) is **mandatory**: every ingest-and-process
run must leave reviewable HTML behind. It never raises — the data is already published and
gated at that point, and a failed map must not look like a failed publish. It reports
failures and both artifacts are regenerable:

```bash
python scripts/generate_qa_report.py {layer_id} --version {version}
python scripts/generate_maps.py     {layer_id} --version {version}
```

A published version is **immutable**. Reprocessing mints a new
`v{YYYY-MM-DD}_{git-sha}` (`-dirty` on an uncommitted tree); publishing onto an existing
version raises unless you pass `on_exists="bump"`. Let it bump rather than overwriting —
the superseded version stays as history and `_VERSION.json` records the chain.

### 5. Hand the user something they can actually open

- **`maps/maps_bundle.zip`** (~9.5 MB) is the whole map collection in one object. The S3
  console downloads one file at a time, so never ask the user to fetch ~20 interlinked
  HTML pages individually. Unzip, open `index.html` — links are relative.
- **`maps/contact_sheet.html`** is the FIRST thing to review — per-member, full 0.5°.
  The pooled maps cannot show a defect confined to one member.
- **`{var}_trend_significant_{scenario}.html`**, not just `{var}_trend_pvalue_*`. A
  p-value field looks plausible almost regardless of what it contains — smooth, bounded,
  mostly mid-range — so reviewing only that satisfies §11 in form and not in substance.
  The masked-trend panel makes a wrong significance mask obvious, because the surviving
  pattern either follows physical geography or it does not.
- **`qa/qa_report.html`** is standalone; no bundle needed.
- Map values serialize at 5 significant figures (`_compact()` in `generate_maps.py`), which
  is display-only — full precision stays in the NetCDF. Do not "fix" this by writing full
  float64: it inflates the collection ~1.5× for precision no colour scale can show.
- Do not merge the collection into one giant self-contained HTML. plotly.js already comes
  from the CDN, so there is nothing to de-duplicate, and it would be ~57 MB with every
  figure instantiated at once.

### 6. Review the QA report before claiming success

`generate_qa_report.py` checks: value classes present and shaped; `lower_ci ≤ median ≤
upper_ci`; zero-width CIs isolated to all-zero or single-model cells; percentile in
[1,100] and oriented to match the declared direction; `trend` is **NaN** in the baseline
decade (fitted slope) or `== 0` for a legacy anchored layer, selected from the file's own
`trend_method`; `sign(trend) == sign(trend_tau)` on the **significant** subset (the
anchored `trend × elapsed_decades == change map` identity is retired); the significance
fields are present, bounded, and their exactly-zero-slope residual is reported;
shared baseline bit-identical across scenarios;
coverage counts consistent; land coverage non-empty.

Read the warnings, don't just check the verdict. And if a check reports itself
**skipped**, treat that as a failure to investigate — a silently skipped invariant is
worse than a failed one.

**When a check fails because it encodes an assumption your layer legitimately breaks, fix
the CHECK — never the data, and never by loosening the check for every layer.** The
cross-scenario baseline-identity check assumed a uniform ensemble, which held for every
layer before `timber_*-tempnle`. It now groups scenarios by the `members_by_scenario`
identity each file declares and asserts bit-identity *within* a group, so uniform layers
behave exactly as before. Two things that made this safe rather than a weakening:

- Group on member **identity**, not count. rcp26 and rcp60 both had 8 members but not the
  *same* 8 (`clm45` contributes different GCM pairs per RCP), and a count-based signature
  demanded bit-identity between panels that legitimately differ.
- If the regrouping leaves **no** group with two scenarios, the check now emits an explicit
  `NOT TESTED` warning. A check that quietly tests nothing is worse than one that fails.

Same principle at the publish gate: `utils/layer_publish.py` rejects scenario files that
disagree on any non-scenario-specific attribute, and it correctly fired on a per-scenario
`n_members`. The fix was to record a `members_by_scenario` breakdown that is byte-identical
in every file — strictly *more* information — not to add `n_members` to
`_SCENARIO_SPECIFIC` and weaken the guard for every existing layer.

**A green verdict is not verification.** These are algebraic self-consistency checks; a
field can satisfy all of them and still be geophysically wrong. **View the maps** and
confirm the geography is plausible — mountain ranges, biome boundaries and coastlines
where they belong. Then say plainly whether you looked: if you have not viewed an image,
report the layer as *unreviewed*, not *verified*. When you do find a visual defect, add an
automated check for its class to `generate_qa_report.py` so the next layer fails loudly.

### 7. Only then clean up raw

`storage.cleanup_raw` refuses unless `_COMPLETE.json` verifies and every input records a
`source_url` — the ISIMIP API is behind Anubis anti-bot, so an un-recorded source may not
be re-downloadable. Never delete raw before the user has reviewed the maps.

## Outputs

- `{variable}_{scenario}_processed.nc` — one file per scenario, dims `(decade, lat, lon)`.
  Never a per-decade file; never one monolith across scenarios.
- Published to `TCFD/tcfd/layers/{layer_id}/{version}/` with `data/` + `layer.json` under
  the `_COMPLETE.json` gate, and `qa/` + `maps/` pinned alongside but ungated.
- Consumers read `…/current/` and must call `storage.verify_complete()` first.

## Reference implementations

`process_csoil_soilcarbon.py` (mixed cadence, 17 members, `EXCLUDED_MODELS`, coverage diagnostics),
`process_burntarea_fire.py` (thick ensemble, Theil-Sen decadal trend),
`process_let_cyclone.py` (thin ensemble → spatial smoothing, zero-inflated → two-tier
percentile), `process_led_drought.py` (binary exposure flag),
`process_timber_tempnle.py` (**two tracks from one parameterized script**; mixed per-tile /
per-gridcell PFT conventions harmonized with a cover floor; model-**family** pooling;
minimum-family mask; per-scenario baseline composition).

Scripts other than these four still use pre-S3 local paths — check before copying one.
