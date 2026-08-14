# Dataset Attributes

Per-dataset facts for both products: what each shipped layer *is*, and the framing decisions
measured for it. [CLAUDE.md](CLAUDE.md) holds the rules that generalize; this file holds the
instances they were measured on.

## What is authoritative for what

This file is an **index and orientation aid**, not a source of truth. Four things outrank it,
and when one disagrees with this file, it is right and this file is stale:

| Authority | For |
|---|---|
| the processed NetCDF's **global attributes** | a layer's own framing — statistic, mask rule, smoothing, percentile direction, units, caveats. Written by the processor at build time and read at delivery time. |
| **`config/layer_registry.yaml`** | delivery routing — folder, status, which slope to read, `relative_baseline`. |
| [OUTPUT-SPEC.md](OUTPUT-SPEC.md) | the output contract every layer satisfies. **Do not restate it here.** |
| the **`/isimip-search-download` skill** + `config/isimip_search_catalog.yaml` | what exists in the ISIMIP repository. |

Deeper per-layer detail than this file carries — per-member metadata quirks, the exact
measurements behind each decision, the reference sites — lives in **the processor's module
docstring** and its dated **[WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md)** entry. Read the one
layer you are working on; see *Scope discipline* in CLAUDE.md before reading the others.

---

# Product 1 — TCFD / CDP layers

## Identity

The **`layer_id`** is the identifier the delivery workflow uses — `config/layer_registry.yaml`,
`config/asset_catalog.yaml`, `config/hazard_taxonomy.yaml` and `layers.csv` all key on it.
The **variable** is what the NetCDF calls itself. They are different vocabularies; keep both
straight.

| Hazard | `layer_id` | Variable | Round / scenarios | Processor |
|---|---|---|---|---|
| Drought (exposure), ISIMIP2b | `drought-2b` | `led` | 2b, rcp26/60 | `process_led_drought.py` (rebuilt 2026-08-11) |
| Drought (exposure), ISIMIP3b | `drought-3b` | `driedarea` | 3b `Heinicke2026`, ssp126/370/585 | `process_driedarea_isimip3b.py` (2026-08-11) |
| Tropical cyclone (exposure) | `cyclone` | `let` | 2b, rcp26/60 | `process_let_cyclone.py` (rebuilt 2026-08-11) |
| Wildfire (burnt area) | `wildfire` | `burntarea-total` | 3b `fire`+`biomes`, ssp126/370/585 | `process_burntarea_isimip3b.py` |
| Crop failure (exposure) | `cropfailure-3b` | `cropfailure` | 3b `Zantout2025`, ssp126/370/585 | `process_cropfailure_isimip3b.py` (2026-08-13) |
| Temperate conifer productivity | `conifer-npp` | `npp-tempnle` | 2b `biomes` CLM45+ORCHIDEE+LPJmL, rcp26/60/85 | `process_tempnle_npp.py` (2026-08-12) |

`conifer-npp` is **not a hazard** — it is an asset-condition layer and is excluded from every
hazard count. `config/hazard_taxonomy.yaml` records that under `non_hazard_layers`.

## Ensemble and framing, per layer

Every value below was **measured for that layer**, never inherited. Two layers agreeing on a
knob is not a reason for a third to adopt it.

| `layer_id` | Members/scenario | Decadal statistic | Mask rule | Smoothing | Percentile |
|---|---|---|---|---|---|
| `cropfailure-3b` | **40** (8 models × 5 GCMs) | `pooled_mean_zero_inflated` | full cropland footprint, no min-model cut | none | single-tier, `higher_is_worse` |
| `drought-3b` | 15 (3 × 5) | `pooled_mean_boolean` | full union, no min-model cut | none | two-tier, `higher_is_worse` |
| `drought-2b` | 31 (8 × 4) | `pooled_mean_boolean` | **≥2 impact models** | none | two-tier, `higher_is_worse` |
| `cyclone` | 4 (1 × 4) | `pooled_mean_zero_inflated` | — | **5×5, L=2.5** | two-tier, `higher_is_worse` |
| `wildfire` | 22 (5 models) | `pooled_median` | — | none | `higher_is_worse` |
| `conifer-npp` | — | `pooled_median` | 2% cover presence mask | none | **`higher_is_better`** — percentile is already inverted in the file |

### Which slope to read

Measured on active cells only (either slope non-zero) at the final decade of each layer's
first scenario. Reproduce with `python scripts/generate_customer_delivery.py --measure-slopes`.

| layer | active cells | `sen==0` | sign agreement | read |
|---|---|---|---|---|
| `conifer-npp` | 25,603 | 0.021 | 0.879 | `sen_slope` |
| `cropfailure-3b` | 39,878 | 1.000 | 0.000 | `ols_slope` |
| `cyclone` | 20,337 | 0.974 | 0.025 | `ols_slope` |
| `drought-2b` | 66,741 | 1.000 | 0.000 | `ols_slope` |
| `drought-3b` | 61,810 | 1.000 | 0.000 | `ols_slope` |
| `wildfire` | 64,039 | 0.740 | 0.226 | `ols_slope` |

Measured 2026-08-12; `cropfailure-3b` added 2026-08-13.

`wildfire` is the one worth noting: its 29.2% grid-level zero fraction suggests Sen is safe
and it is not — the collapse lives in the year-pair *differences*, not the values. That is
why this is measured per layer rather than inferred from `field_nature`.

`cropfailure-3b` reads 1.000 at ssp126 but its Sen slope does lift at higher forcing
(96.5% zero at ssp370, 94.3% at ssp585) as non-zero years become common. `ols_slope` remains
the read at every tier.

### Smoothing

Applied to every member's **annual** map before any pooling, and recorded in the file's
`spatial_smoothing` attribute. The decay length is a measurement, not a constant.

| Layer | Stage-1 smoothing |
|---|---|
| `cyclone` (`let`) | **5×5 window, w = exp(−d/2.5) in grid cells**, cos(lat)-scaled longitude distance, longitude-wrapped, normalized over non-NaN land neighbours. A declared deviation. |
| `cropfailure-3b`, `drought-3b`, `drought-2b`, `wildfire`, `conifer-npp` | **none**, each declared with a reason |

Only `cyclone` is smoothed: raw `let` is sparse storm *tracks* one cell wide and a 4-member
ensemble cannot average them into a projected impact field. `L=2.5` leaves **8.1%** of the
weight on the centre cell; the earlier `L=0.7` kernel kept **32.1%** and left the tracks
visible. The others decline it for stated reasons — `conifer-npp` because a kernel would
bleed productivity into cells where no model places a stand; `drought-2b`/`drought-3b`
because 310 and 150 Bernoulli draws per cell-decade already estimate a frequency and
smoothing would blur real aridity gradients; `wildfire` because 22 members is thick.

`cropfailure-3b` is the case where draw count alone did **not** settle it: roughness 0.347 is
close to `let`'s raw 0.389 despite 400 draws per cell-decade. A **split-half test** settled
it — two disjoint 20-member halves give roughness 0.351 and 0.359 and correlate at Pearson
0.977 / Spearman 0.991, so the roughness is real cropland structure, not sampling noise
(0.790 on footprint edges vs 0.304 in the interior). Reuse that test whenever roughness and
ensemble depth disagree.

**Read the layer's own `spatial_smoothing` attribute before interpreting its values.**

### The third decadal-statistic branch

`pooled_mean_zero_inflated` is a **declared** deviation from the contract's continuous
branch, taken only on measurement and recorded in `decadal_statistic_rationale`. Two layers
qualify today; each adoption was its own decision on its own numbers, and the count is not a
precedent.

| Layer | Annual exact-zero | What the median branch would have published |
|---|---|---|
| `cyclone` (`let`) | 97.84% | erases **93%** of exposed land — 2,684 exposed cells vs 15,122 |
| `cropfailure-3b` | 60.9% within footprint (96.14% counting zero-fill) | erases **96.6%** of exposed cropland — 1,351 exposed cells vs 39,872 |
| `wildfire` (`burntarea`) | 29.2% | **does not qualify** — took the median branch without difficulty |

### Layers that score DEPARTURE FROM A LOCAL BASELINE

Three layers define their value against a fixed historical reference, so a high score means
*"unusual for this place"* and never *"bad in absolute terms"*. This is machine-enforced —
see the rule and mechanism in CLAUDE.md.

| Layer | What a high score means | What it does **not** mean |
|---|---|---|
| `cropfailure-3b` | conditions departing furthest from what this cropland is built for | absolute agricultural risk, yield loss, or food insecurity |
| `drought-3b` | soil moisture departing furthest from a fixed preindustrial reference | that the site is dry |
| `drought-2b` | same, against preindustrial | that the site is dry |

The canonical case, and the one a reader will challenge first: **`cropfailure-3b` puts Iowa at
the 99.3rd percentile of world cropland and the Sahel at 69.4.** That is the measure working
as defined — the intensity of change relative to local norms genuinely is greater in Iowa,
and a farming system built around a stable climate has further to fall. A permanently arid
cell with a stable regime scores LOW on the drought layers for the same reason.

### Per-layer specifics worth knowing before use

- **`cropfailure-3b`** — the publisher **zero-fills the entire globe**: ocean, Antarctica,
  Greenland and the Sahara read exact `0`, not NaN, so `isfinite` carries no footprint and
  must never be used as its mask. Published on a **cropland footprint** (39,890 cells, 42.4%
  of land) derived from where the field is ever non-zero, so a non-agricultural site returns
  `OFF_LAYER_MASK` rather than zero and a percentile is a rank **among cropland**. Each
  member is additionally masked to its own footprint, because only 38.1% of cells have all 40
  members simulating crops. The index carries **no crop token** and the sidecar's
  `specifiers` block names only the variable, so which crops it aggregates is
  **undocumented** — for a specific crop, treat it as a regional signal. Time axis declares
  **no `units` at all**. Sidecars are at `{stem}.json`, **not** `{stem}.nc.json` — the
  latter 404s and reads convincingly as "this publication has no sidecars"; the ingest made
  exactly that error and ran unverified until it was caught (corrected 2026-08-14).
- **`wildfire`** — values are **not clipped at 100%**: a cell that reburns within a year can
  legitimately exceed full coverage (measured annual maxima ~575% for `elm-eca`).
- **`cyclone`** — **single impact model** (`ke-tg-meanfield`), so its CI is inter-GCM spread
  only and carries no structural model uncertainty. Only ISIMIP2b/RCP exists; there is no SSP
  re-issue of TC *exposure*.
- **`conifer-npp`** — reported **per tile** (per unit conifer-stand area), behind a 2% cover
  presence mask. Rising NPP is partly a CO₂-fertilisation response and that belongs in any
  narrative. Covers temperate needleleaf evergreen stands only.
- **`drought-3b`** — the 60–80°N signal is largely one model (`jules-w2`); thin coverage
  concentrates in the arid belt.

### Siblings, not versions

**The two drought layers are both shipped and both current.** `led` (2b, 8 models × 4 CMIP5
GCMs, rcp26/60) and `driedarea` (3b, 3 models × 5 CMIP6 GCMs, ssp126/370/585) — deeper
ensemble versus newer scenarios. Pick per delivery; neither supersedes the other. They
resolve the **minimum-model mask differently on measured evidence** (`led` ≥2, `driedarea`
full union), which is the intended behaviour.

`cropfailure-3b` has an unshipped 2b sibling in the same relationship: `lec` (Lange2020,
GEPIC + PEPIC × 4 CMIP5 GCMs, rcp26/60).

### Withdrawn, superseded and unregistered

| | Status |
|---|---|
| `yield-sug-noirr`, `yield-sug-firr` | **WITHDRAWN 2026-08-11, upstream data defect.** ISIMIP2b LPJmL does not simulate cane in the cane belt (São Paulo, UP India, Queensland, Florida all sentinel-zero). Passed the contract and meant nothing. No ISIMIP source supports a scenario-bearing sugarcane layer. Refused by the loader via `blocked:` in the registry. |
| `burntarea` (2b `biomes`, rcp26/60/85) | **SUPERSEDED 2026-08-10** by the 3b layer. |
| `csoil-total` | Processed but **not registered** for delivery. `process_csoil_soilcarbon.py` is the OUTPUT-SPEC **reference implementation**. |
| Timber, fisheries, health, … | Various `process_*.py`, not registered. |

## Source families in the ISIMIP repository

The `/isimip-search-download` skill and `config/isimip_search_catalog.yaml` are
**authoritative** here; this is orientation only.

Extreme-event exposure exists in two generations:

- **Lange 2020** (ISIMIP2b, rcp26/rcp60) — twelve members: `le{d,r,w,c,h,t}` land-area-exposed
  paired with `pe{d,r,w,c,h,t}` population-exposed twins.
- Its **ISIMIP3b/SSP re-issue** under `DerivedOutputData/`, split across `Heinicke2026`
  (`driedarea`, `floodedarea`) and `Zantout2025` (`heatwave`, `wildfire`, `cropfailure`) —
  named by **hazard word, not `le*` code**.

Never conclude a family has no newer-round version without listing `DerivedOutputData/`. And
a family having no re-issue is not the hazard having no newer data: ISIMIP3b's newest
tropical-cyclone product sits in `InputData/climate/tropical_cyclones/`, outside
`DerivedOutputData/` entirely.

---

# Product 2 — Water Risk Index

Monthly ISIMIP data for 6 water variables → per-month ensemble means plus annual quantile
breakpoints. **Standalone scripts only — not the `isimip-pipeline` CLI**, and none of the
TCFD contract applies: no trends, no percentile scoring, no kernel smoothing.

- **Output**: `C:\Cai_data\WaterIndex\waterIndexUnderlyingData_{var}_ssp.nc` — dims
  `(lat=360, lon=720, scenario=3, value_type=20, decade=9)`
- **Value types**: per-month ensemble means (vt 0–11), annual mean (vt 12), annual quantile
  breakpoints Q05–Q95 (vt 13–19). Quantile annual aggregation **always uses mean, not sum**,
  to keep vt 13–19 in the same units as vt 0–12.
- **Normalization**: robust z-score per impact model (median/IQR from a 2015–2024 reference
  period → target mean 1000, SD 200), applied before ensemble averaging **only when model
  scales diverge significantly** (e.g. TWS). Decided per variable.
- **Units**: output units match the original RCP files. Per-variable
  `unit_conversion_factor` in `config_water_variables.py`.
- **QA/QC**: `validate_water_tws.py` (quantile ordering, annual mean consistency, seasonal
  sanity, cross-scenario checks); `compare_water_index.py` (RCP vs SSP comparison with
  Theil-Sen slope maps and spatial Spearman R²).

| Variable | Aggregation | Output Units | Notes |
|---|---|---|---|
| `tws` | **mean** | kg m-2 (normalized) | Stock; 4 models, normalized to synthetic units |
| `rootmoist` | **mean** | % max capacity | Stock; WEB-DHM-SG only (÷1187.29 × 100) |
| `qr` | **sum** | kg m-2 s-1 | Flux; 4 models, raw ISIMIP units, no normalization |
| `dis` | **mean** | m3 s-1 | Stock; 5 models, no normalization, raw ISIMIP units |
| `potevap` | **sum** | kg m-2 s-1 | Flux; 4 models, h08 selectively normalized to the reference ensemble (cwatm/miroc/watergap2-2e) |
| `precip` | **sum** | TBD | TODO — climate forcing InputData, not model output |
