---
name: isimip-process-visualize
description: Process annualized ISIMIP NetCDF into the TCFD 6-value-class format, then publish the layer to S3 with its QA/QC report and interactive maps. Use when processing a downloaded ISIMIP variable into a TCFD/CDP layer, reprocessing an existing layer, or regenerating QA evidence and visualizations.
---

# Process & Visualize (TCFD/CDP product)

**This skill is for the TCFD/CDP product only** — 6 value classes, annualized decadal
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
- **§10 — trend must be a decadal signal**, not a within-decade annual slope, for noisy
  variables. Prefer the baseline-anchored rate
  `(median[decade] − median[baseline]) / elapsed_decades`.
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
- `baseline_decade`, `baseline_source` — the shared 2020s baseline must be **identical
  across scenarios**, computed from all scenarios with overlapping 2020s data.
- Zero-inflated hazards → two-tier percentile (zeros → 1; non-zeros ranked against the
  non-zero baseline → [2,100]).

Emit `n_members` / `n_models` per cell whenever the ensemble's land masks differ, so the
CI is auditable. Do not silently mask thin cells — that is a product decision for the user.

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
[1,100] and oriented to match the declared direction; `trend == 0` in the baseline decade;
`trend × elapsed_decades == change map`; shared baseline bit-identical across scenarios;
coverage counts consistent; land coverage non-empty.

Read the warnings, don't just check the verdict. And if a check reports itself
**skipped**, treat that as a failure to investigate — a silently skipped invariant is
worse than a failed one.

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
`process_burntarea_fire.py` (thick ensemble, baseline-anchored trend),
`process_let_cyclone.py` (thin ensemble → spatial smoothing, zero-inflated → two-tier
percentile), `process_led_drought.py` (binary exposure flag).

Scripts other than these four still use pre-S3 local paths — check before copying one.
