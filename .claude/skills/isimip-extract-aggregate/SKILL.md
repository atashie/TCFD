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
- **`trend_definition`** — for a baseline-anchored rate, `trend[decade]` is the rate **from
  the baseline decade to that decade**, not a within-decade slope. `trend[baseline] ≡ 0`,
  so never report a baseline-decade trend as a finding.
- **`units`** and **`trend_units`** — these differ (e.g. `kg m-2` vs `kg m-2 decade-1`).
- **`baseline_decade` / `baseline_source`** — the 2020s baseline is shared across scenarios
  by design; identical 2020s values across scenarios are correct, not a bug.
- **`n_members` / `n_models`** where present — cells backed by one model have a CI that
  reflects GCM spread only. Flag or filter these for a client rather than presenting them
  at equal confidence.

## Extraction

- **Points**: nearest-cell lookup on the 0.5° grid (all TCFD layers are 360×720,
  lat −89.75…89.75, lon −179.75…179.75). Report the actual cell centre used alongside the
  requested coordinate, and flag any point landing on a NaN (ocean / outside the model's
  land mask) rather than emitting a silent blank.
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
