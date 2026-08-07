---
name: isimip-extract-aggregate
description: Extract values from a processed TCFD layer by point location or region polygon and aggregate them to CSV for customer delivery. Use when a client needs site-level or regional climate-risk values pulled out of a processed layer.
---

# Extract & Aggregate (TCFD/CDP product)

Pulls values out of an already-processed layer for delivery. This does **not** process or
reprocess data — if the layer does not exist yet, use `isimip-process-visualize` first.

## Read from a verified layer

Read processed layers from `data/processed/{layer_dir}/{variable}_{scenario}_processed.nc`.

Before extracting, confirm the layer is complete and reviewed:
- every expected scenario file is present (discover them dynamically — glob
  `{variable}_*_processed.nc`, never a hardcoded scenario list; see GUARDRAILS.md §3)
- the QA report was generated and its warnings read
- the maps were actually viewed

Extracting from an unreviewed layer produces numbers that look authoritative and may not be.
If you have not confirmed the above, say so in the delivery rather than implying otherwise.

## Extraction modes

`scripts/utils/spatial_extract.py` provides:

- **`extract_by_point(ds, lat, lon)`** — Gaussian distance-weighted average over nearby cell
  centres (σ = 0.25, half a cell). Longitudes are normalized to [-180, 180]. NaN cells are
  excluded and the remaining weights re-normalized.
- **`extract_by_polygon(...)`** — area-weighted average using shapely intersection of the
  polygon with each 0.5° cell.

**Respect `percentile_direction`.** A layer declaring `higher_is_better` (stored carbon,
biomass — assets, where the risk is *loss*) carries an inverted percentile. Use
`get_percentile_direction()` / `apply_percentile_inversion()` rather than assuming higher =
worse. Getting this backwards silently reverses the risk ranking in a client deliverable.

## Known defect — polygon area weighting

`_calculate_cell_weights` computes `intersection.area / cell.area` in **planar degree
space** with no `cos(lat)` term. The per-cell coverage *fraction* is fine (the cosine
cancels within a single cell), but the **cross-cell normalization is not**: a cell at 60°N
has roughly half the true surface area of an equatorial one, yet contributes equal weight
for equal coverage.

Consequence: any polygon spanning a wide latitude band over-weights its high-latitude cells.
Small mid-latitude polygons are barely affected; large or high-latitude regions are.

Until fixed, either restrict polygon extraction to compact regions or state the caveat in
the deliverable. The fix is to multiply each weight by `cos(radians(lat))` before
normalizing.

## Output

CSV via `scripts/utils/export_formatter.py`, with a metadata sidecar recording: the layer
and its processed-file provenance, the extraction date, the mode (point vs polygon), the
`percentile_direction` applied, and the units.

**There is no significance statistic under the current contract**
([OUTPUT-SPEC.md](../../../OUTPUT-SPEC.md)). The dual-slope design replaced it: agreement
between `ols_slope` and `sen_slope` is the robustness signal, and disagreement means a
cell's trend is not robust.

The legacy 28-column schema still carries `Decadal_Trend_Significance` /
`Long_Term_Trend_Significance`. Write `export_formatter.SIGNIFICANCE_NOT_COMPUTED` into
them — **never leave them NaN.** That distinction is not cosmetic: the retired schema
declared a `significance` value class no processor ever emitted, so those columns silently
resolved to NaN and every layer ever delivered read as *"not significant"* rather than
*"not computed"*.

`Decadal_Trend_Strength` takes `sen_slope` by default. On a **zero-inflated** layer
`sen_slope` is exactly 0 almost everywhere and `ols_slope` is the correct source instead —
`test_shared_baseline.py` prints the zero fraction, so check it before exporting.

## Scenario handling

Deliver projection scenarios only — SSP `ssp126/245/370/585`, RCP `rcp26/45/60/85`.
`picontrol` and `historical` strengthen the baseline but are not client-facing.

Discover which scenarios exist from the filesystem. A hardcoded list once made 25% of
processed data invisible.
