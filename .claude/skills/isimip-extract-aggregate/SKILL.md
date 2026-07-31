---
name: isimip-extract-aggregate
description: Extract values from a published TCFD layer by point location or region polygon and aggregate them to CSV for customer delivery. Use when a client needs site-level or regional climate-risk values pulled out of a processed layer.
---

# Extract & Aggregate (TCFD/CDP product)

Pulls values out of an already-published layer for delivery. This does **not** process or
reprocess data — if the layer does not exist yet, use `isimip-process-visualize` first.

## Always read from a gate-verified layer

```python
version  = storage.resolve_current(layer_id)          # or pin an explicit version
vprefix  = storage.version_prefix(layer_id, version)
storage.verify_complete(vprefix, require=["layer.json"])   # REQUIRED before consuming
data_dir = storage.pull_prefix(f"{vprefix}/data")
```

A version prefix whose `_COMPLETE.json` is absent or mismatched is **in-flight or corrupt,
never data**. Skipping `verify_complete` risks delivering a half-written layer to a
customer.

**The gate proves completeness, not correctness.** It says the bytes all arrived; it says
nothing about whether the layer is geophysically sound. Before an extraction goes to a
customer, confirm the version's `qa/qa_report.html` has no unresolved warnings **and** that
its maps were actually reviewed (GUARDRAILS §11) — a layer once passed every algebraic
check while carrying a member on a ~4°×5° grid. Also read `effective_resolution` in the
layer attributes: a point value inherits the resolution of whichever member covers it, so
site-level numbers from a coarse or single-model cell should be flagged, not delivered flat.

**Pin the version in the deliverable.** Record `layer_id`, `version` and the extraction
date in the CSV (or a sidecar) so a client number is always reproducible. `current/` moves
when a layer is reprocessed; a delivered number must not.

## What the values mean — carry the semantics through

Read `layer.json` / the NetCDF global attrs and propagate them into the output. Getting
these wrong silently inverts a client's risk reading:

- **`percentile_direction`** — `higher_is_worse` (hazards) vs `higher_is_better` (assets
  like stored soil carbon). For `higher_is_better` the percentile is already **inverted**,
  so a *high* percentile means low stock / high risk. Do not re-invert it.
- **`trend_definition` / `trend_method`** — since 2026-07-30 `trend[decade]` is a
  **Theil-Sen slope of the DECADAL MEDIAN series** from the baseline decade to that decade
  (`theil_sen_on_decadal_median_series`), in value per decade. Never a within-decade slope.
  **`trend[baseline]` is NaN**, not 0 — a fitted slope has no elapsed period — so never
  report a baseline-decade trend, and do not read the blank as "no change".
  **`trend × elapsed_decades` no longer equals the change map**; if a client needs the
  change, compute `median[decade] − median[baseline]` directly.
  **An exactly-zero trend is not always "flat"**: Theil-Sen is a median of pairwise slopes,
  so on zero-inflated hazards 10–14% of cells return exactly 0 (0.03–3.7% on continuous
  layers). Where a zero slope sits beside a significant `trend_pvalue`, that is the tie
  residual, not evidence of stability — read `trend_tau` for direction.
- **`trend_pvalue` / `trend_tau` / `trend_n_obs`** — Mann-Kendall on the ensemble-mean
  **ANNUAL** series (n = 20…80 **years**, members averaged within each year), so the test
  and the slope are fitted on **different series** by design. `trend_pvalue[last decade]`
  **is** the long-term p-value because the window expands. A constant series gives
  **p=1.0, not NaN**. It measures monotonicity of the ensemble MEAN, **not** inter-model
  agreement — on a smooth stock like `csoil` it saturates near 89% at *every* scenario and
  carries no severity information. **A polygon's p-value must be RECOMPUTED** from its
  area-weighted annual series via `trend_significance.mk_pvalue()`; averaging per-cell
  p-values is meaningless.
- **`units`** and **`trend_units`** — these differ (e.g. `kg m-2` vs `kg m-2 decade-1`).
- **`baseline_decade` / `baseline_source`** — the 2020s baseline is shared across scenarios
  by design; identical 2020s values across scenarios are correct, not a bug.
- **`n_members` / `n_models`** where present — cells backed by one model have a CI that
  reflects GCM spread only. Flag or filter these for a client rather than presenting them
  at equal confidence.
- **`known_issues`** where present — an OPEN, unresolved caveat the layer ships with. Read it
  before quoting any number. `wildfire_burntarea_annual` carries one: values poleward of
  ~70°N are unreliable (one model projects up to ~100%/yr burnt on Arctic islands).

## Extraction

- **Points**: nearest-cell lookup on the 0.5° grid (all TCFD layers are 360×720,
  lat −89.75…89.75, lon −179.75…179.75). Report the actual cell centre used alongside the
  requested coordinate, and flag any point landing on a NaN (ocean / outside the model's
  land mask) rather than emitting a silent blank.
- **Coastal points need a caveat.** For any per-cell **fraction** metric (`burntarea`, `led`,
  `let`, `lew`, `fldfrc`) the value is a fraction of the *whole* cell, so a mostly-water
  coastal cell reads low no matter how exposed its land is — and the coastal ring is also the
  thinnest-covered part of the ensemble. Measured on `wildfire_burntarea_annual`: coastal
  cells are 12.4% of land, their median is **18× below** the interior median, and **67.5%**
  rest on a single model (vs 3.7% inland). Most client sites are coastal, so flag a coastal
  hit explicitly rather than reporting it as a low-risk finding.
- **Regions**: area-weight by `cos(lat)` — unweighted means over a lat/lon grid
  over-weight high latitudes. Use `scripts/extract_region_polygons.py` as the reference.
- **Aggregating across cells**: use `np.nanmean`/`np.nanmedian`, and report the count of
  contributing cells so a mostly-NaN region is visible.
- Never aggregate a **percentile** across cells by averaging it — ranks are not additive.
  Aggregate the underlying value, then rank, or report the distribution.

## Output

Exports go to `TCFD/exports/{customer}/{run_date}/*.csv` via `storage.export_prefix(...)` —
never the repo, never local disk. Include: location id, requested lat/lon, cell lat/lon,
scenario, decade, value class, value, units, `layer_id`, `version`.

Long format (one row per location × scenario × decade × value class) is safer than wide —
it survives a layer gaining a decade or a value class.

## Reference scripts

`extract_region_polygons.py`, `extract_timber_locations.py`,
`generate_extraction_report.py`. **These still use pre-S3 local paths** — they predate the
storage migration, so re-point them through `isimip_pipeline/storage.py` rather than
copying their I/O as-is.
