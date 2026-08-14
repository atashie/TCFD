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
| Heatwave (exposure), ISIMIP3b | `heatwave-3b` | `heatwave` | 3b `Zantout2025`, ssp126/370/585 | `process_heatwave_isimip3b.py` (2026-08-14) |
| Heatwave (health threshold), ISIMIP2b | `heatwave-2b` | `leh` | 2b `Lange2020`, rcp26/60 | `process_leh_isimip2b.py` (2026-08-14) |
| Permafrost thaw | `permafrost-3b` | `thawfrac` (from `thawdepth`) | 3b `permafrost`+`biomes`, ssp126/370/585 | `process_thawdepth_permafrost.py` (2026-08-14) |
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
| `heatwave-3b` | 5 (**1 index model** × 5 GCMs) | `pooled_mean_boolean` | full footprint; min-model rule **void**, not relaxed — `n_models`≡1 and all 5 GCM members share one 67,420-cell mask | none | single-tier, `higher_is_worse` |
| `heatwave-2b` | 4 (**1 index model** × 4 GCMs) | `pooled_mean_boolean` | **explicit ISIMIP2b `LSM` land mask** — the source zero-fills all 259,200 cells, so `isfinite` is not a footprint; min-model rule void | none | two-tier, `higher_is_worse` |
| `permafrost-3b` | 12 (3 models × {5,5,2} GCMs) | **`pooled_mean_multimodel`** | 2020s permafrost footprint, **≥2 of 3 models** | none | two-tier, `higher_is_worse` |
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
| `heatwave-3b` | 67,171 | 1.000 | 0.000 | `ols_slope` — **but see the saturation caveat below; on this layer a near-zero slope is ambiguous** |
| `permafrost-3b` | 14,455 | 0.152 | 0.844 | `ols_slope` — Sen has *not* collapsed here and is still wrong; see below |

Measured 2026-08-12; `cropfailure-3b` added 2026-08-13; `heatwave-3b` and `permafrost-3b`
added 2026-08-14.

**`permafrost-3b` is the counter-example to "low `sen==0` means Sen is safe".** Its Sen share
is the second-lowest of any layer (0.152, against `conifer-npp`'s 0.021) and sign agreement is
high (0.844), which reads as a well-behaved field — and Sen is still the wrong estimator,
for a reason none of the other layers exhibit. On a **multimodal** ensemble the pairwise
sample is dominated by *cross-cluster* member pairs, whose slopes carry the level offset
rather than the trend, and the median over those pairs is dragged toward the flat cluster.
Measured on ssp585 over the published footprint: the 12 members' own OLS trends span
+0.0104…+0.0826 dec⁻¹, their mean is **+0.0326**, the published `ols_slope` is **+0.0326**
(exact), and `sen_slope` is **+0.0069 — below every single member**. An estimator outside the
range of the things it summarises is not robust. Read `ols_slope`, and read it with the
Anomaly panel, because both slopes are censored where the column is fully thawed.

`wildfire` is the one worth noting: its 29.2% grid-level zero fraction suggests Sen is safe
and it is not — the collapse lives in the year-pair *differences*, not the values. That is
why this is measured per layer rather than inferred from `field_nature`.

### `heatwave-3b` saturates — the one layer where a flat slope means the opposite

Every other layer here follows the standing rule: the two slopes fail in *opposite* regimes,
so their disagreement flags a fragile trend. **`heatwave-3b` breaks that rule**, and it is the
only shipped layer that does.

Exposure is defined as the annual HWMId exceeding the 97.5th percentile of *that cell's own*
preindustrial control, so warming pushes cells permanently over the threshold and the binary
flag pins at 1. Measured on the shipped files — pooled decadal exposure frequency, and the
share of published cells sitting at exactly 1.0:

| scenario | 2020s | 2050s | 2090s | cells at 1.0, 2090s | percentile ≥ 99.5, 2090s |
|---|---|---|---|---|---|
| ssp126 | 0.314 | 0.473 | 0.462 | 1.2% | 2.1% |
| ssp370 | 0.314 | 0.607 | 0.846 | 32.6% | 39.4% |
| ssp585 | 0.314 | 0.650 | 0.903 | **45.9%** | **51.9%** |

At 1.0 the pooled sample has no variance: the CI collapses to zero width, **both** slopes go
to ~0 and **agree** there, and the percentile ties at 100. So on this layer a near-zero slope
is ambiguous between "no trend" and "pinned at the ceiling", and agreement between the
estimators is not reassurance.

The censoring **inverts trend rankings between regions**. Measured on ssp585: the Amazon
(10°S–0, 70–55°W) goes 0.601 → 1.000 and 0% → 100% saturated while its `ols_slope` *falls*
from +0.160 to +0.046 dec⁻¹; Siberia (60–70°N, 80–120°E) never saturates and its slope *rises*
from +0.069 to +0.098. On the 2090s panel Siberia out-trends the Amazon 2.1×, which reads as
"the Amazon has stabilised" and means "the Amazon is exposed every year in every member".

**Identify saturated cell-decades as `median == 1.0` (equivalently a zero-width CI at 1.0)**
and treat their slopes and percentile ranks as censored. How the tied top block should be
ranked, and whether to publish a time-to-saturation measure, are open decisions.

### `heatwave-3b` is a relative index with no humidity term

`HWMID-NONE` is HWMId **only** — the "NONE" is the humidity term. Its ISIMIP2b predecessor
`leh` (`HWMId-humidex`) additionally required Humidex ≥ 45 so that counted events "would also
adversely affect human health"; Zantout et al. 2025 drops that criterion. **This layer cannot
evidence a wet-bulb, heat-stress, workforce-safety or equipment-derating claim.**

It also carries a latitudinal artifact worth stating before anyone maps it: on the 2020s
baseline, |lat| ≤ 23.5 reads **0.561** against **0.149** for |lat| > 50 — a 3.77× ratio before
any warming has accumulated — because a relative threshold is crossed sooner where interannual
variance is low, and the tropics have the lowest temperature variance on the planet. Measured
site example: Chicago 0.688 reads *above* Delhi 0.548.

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

`permafrost-3b` reuses the test and adds a caution about **how the halves are drawn**. Split
alphabetically off a sorted member list, the two halves got 3 JULES + 2 LPJmL against 2 + 3 —
different models, therefore different permafrost domains — and the test returned Pearson
**0.376** with each half reading *smoother* than the full ensemble, which is the tell, since
fewer members should be noisier. Stratified by model (alternating GCMs **within** each model)
the same panels give roughness 0.114 with halves at 0.115 / 0.139 and **r = 0.929**. Draw the
halves so every model appears in both, or the test measures composition rather than noise.

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

### The fourth branch, and why `permafrost-3b` needed it

`pooled_mean_multimodel` exists because the median assumes the pooled sample has **one mode**.
This layer's value is thaw depth divided by *each model's own soil column*, and the three
models sit in two separated clusters on that axis in the 2020s:

| | CLASSIC (61.4 m column) | LPJmL (13.0 m) | JULES (3.0 m) |
|---|---|---|---|
| 2020s normalised thaw | 0.035 | 0.046 | **0.951** |
| members | 2 | 5 | 5 |

Seven members low, five high. Under the median branch the ssp585 spatial median went
**0.40 (2080s) → 0.93 (2090s)** — the median *crossing* between clusters as the high group
gained the majority, which reads as thaw suddenly accelerating and is an artifact of the
estimator. The mean moves smoothly and the SD carries the disagreement (2-SD width, median
**0.748** on a [0,1] field). The percentile also behaves better: `corr(median, percentile)`
rose from +0.712 to +0.973 on the same data.

**A threshold on the central value inherits that choice.** "Area whose column is fully
thawed" means *>half the members* under a median and *effectively all of them* under a mean:
**7.97 M km² vs 0.40 M km²** for the same ensemble at ssp585 2090s. The invariant quantity is
the **member share** — measured there as 0.478 of members on average, with 10.52 M km² where
at least half agree and 0.69 M km² where all twelve do. Report the agreement spread, not one
threshold.

### `permafrost-3b` — what it is, and what it is not

- **Value**: annual maximum thaw depth ÷ that model's soil column, in [0,1]. 0 = never thaws,
  1 = nothing frozen left. **The endpoints are commensurable across models; the interior is
  only ordinally so** — 0.5 of a 3 m column is not physically 0.5 of a 61.4 m column. The
  customer-facing quantity is therefore the **change against the 2020s**: the share of the
  column transitioning from permafrost to none (`transition_summary.json`, and the dashboard's
  Anomaly panel). ssp585 reaches **0.250** area-weighted by the 2090s, against 0.182 (ssp370)
  and 0.060 (ssp126).
- **A permafrost-free cell is published AT the column depth, not as NaN** — Nairobi and Paris
  read exactly the column in all three models — so `isfinite` covers 27.5% of the grid and
  carries no domain information. The footprint is the **2020s permafrost domain** where ≥2 of
  3 models find frozen ground: 15,509 cells, 18.17 M km², against Obu et al. 2019's ~14 M km²
  permafrost area and ~21 M km² permafrost region. Union of all three would be 31.07 M km²,
  above the entire observed region, because LPJmL alone claims 30.16.
- **Do not judge a member by its headroom.** JULES sits near 1 and looks censored, but in raw
  metres it is the *deepest-column* model, CLASSIC, that is unphysical: 2020s thaw p95 of
  **28.0 m**, and no permafrost at all at Fairbanks. LPJmL is closest to observed active-layer
  thickness (p50 0.83 m against an observed ~0.3–1.5 m), JULES next (p50 2.85 m, censored at
  3), CLASSIC last (p50 2.15 m but p75 4.90). JULES's baseline pattern also agrees with the
  others as well as they agree with each other (ρ 0.785 / 0.657, against 0.631 between CLASSIC
  and LPJmL). All 12 members are retained on that evidence.
- **The models agree on how much and disagree on where.** Over each model's own domain the
  2090s ssp585 loss share is comparable (CLASSIC 31.0%, JULES 58.0%, LPJmL 32.3%) — but on the
  cells CLASSIC and LPJmL *both* call permafrost, CLASSIC loses 28.6% and LPJmL 4.0%, a
  **Jaccard overlap of 2.8%**. Comparable totals, different maps. This is what the wide CI is
  reporting; a narrow CI here would mean the ensemble had been trimmed.
- **Not a ground-stability or subsidence layer.** Thaw depth carries nothing about ground ice,
  excess ice, thaw settlement or slope, and ISIMIP publishes no variable for any of them, so a
  cell losing its permafrost is not thereby a foundation-damage forecast. **Solifluction** has
  no representation anywhere in ISIMIP and remains uncovered.
- **Sector trap**: `thawdepth` is published byte-identically under `permafrost`, `biomes` and
  `water_global`, and the two models that are *not* in the permafrost sector (JULES, CLASSIC)
  are two thirds of the ensemble. Walking one sector answers "1 model, 5 members" to a
  question whose answer is 3 models and 12 members. Ingest from one sector, but look in all.

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


### The two heatwave layers are different indices, and they fail in opposite directions

`heatwave-2b` (`leh`) and `heatwave-3b` (`heatwave`) are **siblings, not versions**, in the
stronger sense than the drought pair: the drought layers measure the same thing with
different ensembles, whereas these two measure **different things**.

| | `heatwave-2b` (`leh`) | `heatwave-3b` (`heatwave`) |
|---|---|---|
| Index | HWMId **and** Humidex ≥ 45 | HWMId **only** — "NONE" is the humidity term |
| Health anchoring | yes — Humidex 45 is Environment Canada's "great discomfort; avoid exertion" band | none |
| Round / scenarios | 2b, rcp26 + rcp60 (no rcp85 exists) | 3b, ssp126/370/585 |
| Members per scenario | 4 (1 index model × 4 CMIP5 GCMs) | 5 (1 index model × 5 CMIP6 GCMs) |
| Land at exactly 0, final panel | **65.7%** (rcp60) | 0.0% |
| Cells at the ceiling (1.0), final panel | 0.0% | **45.9%** (ssp585) |
| Land active in **all** members, baseline | **8.4%** | **74.8%** |
| Active area resting on one GCM | **22%** | 0.4% |
| GCM spread (per-member land mean) | 2.16×, CV 28.5% | 1.58×, CV 17.5% |
| Failure mode | **silence** outside the humid-heat belt | **saturation** at high forcing |

**`heatwave-3b` is the selected layer** (user determination 2026-08-14, after reviewing the
dataset). It is `status: preferred` and carries `qa_reviewed_on: 2026-08-14` — the only layer
in the registry with a human QA date. `heatwave-2b` is retained, built and contract-passing,
as the only health-anchored heat source in the repository, but is not selected and stays
blocked pending its own review.

**Neither can be called "the heat layer" without the caveat that distinguishes it.** A zero
in `heatwave-2b` means "never crossed Humidex 45", not "no heat risk" — Paris, Frankfurt and
Yakutsk are exactly zero in every member and every year. A flat slope in `heatwave-3b` means
"pinned at the ceiling", not "no trend".

They also disagree about *where* heat matters, and the disagreement is structural rather
than noise: `heatwave-2b`'s footprint is a humid-heat belt (Central America, Amazon, Sahel
and West Africa, Horn of Africa, Arabian Gulf, Indus/Ganges, maritime SE Asia, northern
Australia), while `heatwave-3b` is brightest in the *same* tropics but for the opposite
reason — low interannual variance makes a relative threshold easy to cross there.


### Two must-disclose caveats now travel with these layers automatically

`LAYER_ATTRS_EXPORTED` in `scripts/utils/delivery.py` is a **closed allowlist**: a caveat
written into a processed file under any attribute name not on that list is silently dropped
before it reaches `layers.csv`, the caveat generator, or either report. Added 2026-08-14 so
the two heatwave layers' defining limitations survive the trip:

| Attribute | Caveat ID | Severity |
|---|---|---|
| `saturation_caveat` | `CENSORED-CEILING-{layer}` | **must-disclose** |
| `sparsity_caveat` | `CENSORED-FLOOR-{layer}` | **must-disclose** |

Both are promoted to must-disclose for the same reason `relative_baseline_note` was on
2026-08-13: **every number is correct and the reader's conclusion inverts.** A layer pinned
at a bound has no variance in the pooled sample, so the CI collapses to zero width, *both*
slopes go to ~0 and **agree** — the dual-slope disagreement rule gives no warning — and the
percentile ties. A flat trend then reads as "stable" and means "maximally exposed, no
headroom left"; a floor reads as "low risk" and means "never crossed the threshold this
index measures".

Any future layer with a physical or definitional bound should set one of these attributes
rather than describing the problem in prose that nothing reads.
