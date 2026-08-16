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

## Spatial support: every layer here is a REGIONAL result, not a site result

**Read this before quoting any number at a location.** The value returned for a site is the
statistic for the whole cell containing that site. It is not the value at the site, and no
amount of precision in the number changes that.

**Resolution is a per-layer property, not a product constant** — it was one until
2026-08-14, and the habit of saying "the grid is 0.5°" outlived the fact.

| Grid | Layers |
|---|---|
| **0.5°** (~55 km) — the default and preferred grid | every layer except the three below |
| **0.25°** (~28 km, 15 arcmin) | `flood-3b-flopros`, `flood-3b-40yr`, `flood-3b-none` |

Each file declares its own `spatial_resolution_degrees`, `test_shared_baseline.py` checks
that it matches the coordinates, and delivery geometry follows each layer's own cell size —
so a 0.25° layer blends a proportionally smaller neighbourhood rather than inheriting the
0.5° window. `manifest.json` records the grid and extraction geometry **per layer**, and the
delivered README states them per layer whenever a delivery spans more than one grid. Never
compare a cell count across layers on different grids without saying so: the flood layers
publish 251,890 cells where a 0.5° layer publishes ~70,000, and that is resolution, not
coverage.

For hazards that vary over comparable distances — drought, heatwave, wildfire — that is a
tolerable approximation and the cell value is a fair description of the site. For hazards
that turn on **fine terrain it is not**, and two families are in that class:

- **`sealevel-2b`** — coastal inundation depends on metres of elevation over hundreds of
  metres of ground, so a quayside and a hillside two kilometres apart sit in the same cell
  and receive the **same number**.
- **the nine threshold rungs** (`heatdays-*`, `tropicalnights-*`, `frostdays-*`, `icedays-id`,
  added 2026-08-16) — temperature falls ~6.5 °C per km of elevation, and on a day count that
  is the difference between never crossing a threshold and crossing it for months. Measured
  on each rung's own baseline, **adjacent** 0.5° cells differ by up to 59 d/yr at the 99th
  percentile on `hd35` (max 296) and 110 d/yr on `hd30`; a cell spanning a valley floor and a
  ridge varies over its own footprint by a comparable amount, and publishes the cell mean,
  which matches neither. Urban heat island is not represented either, so a city centre runs
  hotter than its cell. Visible in the `FD` contact sheet as the Andes rendering as a single
  narrow ribbon of cells across an otherwise frost-free continent.

These layers are a **first-pass screen**. They rank which regions, coastlines and sites
deserve investigation. They cannot support a site-level conclusion, a design elevation, an
asset-level financial estimate, or any statement about an individual building or berth —
those need a site-specific study at metre-scale elevation, and for coastal work that is
expected in future deliveries rather than from this product.

A layer may declare a `resolution_caveat` attribute; where present it is promoted to a
**`must_disclose`** caveat and both reports refuse to render without it (added 2026-08-14
with `sealevel-2b`; ten layers carry it as of 2026-08-16). Note what a layer's resolution
does and does **not** mean: `sealevel-2b` reads terrain at 15″ (~460 m) and uses the **full
elevation distribution inside each cell** rather than a cell mean — the coarseness is in the
*published support and the sea-level forcing*, not in the terrain that fed it. The threshold
rungs are the opposite case: their input **is** a 0.5° cell-mean temperature, so the
sub-cell variation is genuinely absent from the data rather than collapsed at publication.

**The attribute is set only where it applies, and that is enforced by convention, not by
code.** `generate_delivery_caveats.layer_caveats()` promotes `resolution_caveat`,
`saturation_caveat` and `sparsity_caveat` on the test *"is this attribute non-empty"* — not
on it asserting anything. So writing `saturation_caveat: "none — 0.69% at the ceiling"`
publishes a **must-disclose caveat whose body says "none"**, in both reports. Record a
negative measurement under a name that is **not** on `LAYER_ATTRS_EXPORTED` — the threshold
rungs use `saturation_measured` / `sparsity_measured` — so it stays auditable in the file
without reaching a filing. Seven of the nine rungs shipped with exactly that defect for a
few hours on 2026-08-16 before it was caught.

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
| Sea level rise (inundation) | `sealevel-2b` | `coastalinundation` | 2b `sealevelrise` x GEBCO_2026, rcp26/60 | `process_coastal_inundation.py` (2026-08-14) |
| Riverine flood (inundation, FLOPROS defences) | `flood-3b-flopros` | `fldfrcmax-flopros` | 3b `TipESM2025` CaMa-Flood, ssp126/370/585 | `process_fldfrcmax_isimip3b.py` (2026-08-14) |
| Riverine flood (inundation, uniform 40-yr standard) | `flood-3b-40yr` | `fldfrcmax-40yr` | 3b `TipESM2025` CaMa-Flood, ssp126/370/585 | `process_fldfrcmax_isimip3b.py` (2026-08-14) |
| Riverine flood (inundation, no defences) | `flood-3b-none` | `fldfrcmax-none` | 3b `TipESM2025` CaMa-Flood, ssp126/370/585 | `process_fldfrcmax_isimip3b.py` (2026-08-14) |
| Temperate conifer productivity | `conifer-npp` | `npp-tempnle` | 2b `biomes` CLM45+ORCHIDEE+LPJmL, rcp26/60/85 | `process_tempnle_npp.py` (2026-08-12) |
| Soil organic carbon | `csoil` | `csoil-total` | 3b `biomes`, ssp126/370/585 | `process_csoil_soilcarbon.py` (rebuilt 2026-08-15) |
| Chronic heat (day threshold) | `heatdays-hd30/hd35/hd40/hd45` | `hd30`/`hd35`/`hd40`/`hd45` | 3b bias-adjusted daily `tasmax`, ssp126/370/585 | `process_tasthresh.py` (2026-08-16) |
| Chronic heat (night threshold) | `tropicalnights-tr20/tr25` | `TR20`/`TR25` | 3b bias-adjusted daily `tasmin`, ssp126/370/585 | `process_tasthresh.py` (2026-08-16) |
| Cold / frost | `frostdays-fd`, `frostdays-fdm10` | `FD`/`FDm10` | 3b bias-adjusted daily `tasmin`, ssp126/370/585 | `process_tasthresh.py` (2026-08-16) |
| Cold / ice days | `icedays-id` | `ID` | 3b bias-adjusted daily `tasmax`, ssp126/370/585 | `process_tasthresh.py` (2026-08-16) |

`conifer-npp` is **not a hazard** — it is an asset-condition layer and is excluded from every
hazard count. `config/hazard_taxonomy.yaml` records that under `non_hazard_layers`.

`csoil` is the **same shape and its classification is an open user decision**: it is
`higher_is_better` with an inverted percentile, and declining soil carbon reads more like the
asset's condition than a hazard acting on the asset — unless the asset *is* the soil. It is
registered `status: alternate` with **no asset type routing to it** and the `soil-degradation`
family's `covered_by` deliberately left empty, so nothing counts it as a hazard yet. It is
also **not** a soil-degradation layer: soil organic carbon is one of the ten degradation
processes in recital 4 of Directive (EU) 2025/2360, and land use is held fixed in every
member, so the layer cannot see management-driven loss at all.

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
| `sealevel-2b` | 4 (**1 sea-level model** × 4 GCMs) | `pooled_median` | GEBCO `TID==0` land, ocean-connected; **per-cell consensus mask** drops a member deviating >0.5 m from the all-member median | none | **absolute calibration**, `higher_is_worse` — NOT percentile-of-score |
| `flood-3b-flopros` | 27 (6 models × {2,5,5,5,5,5} GCMs) | `pooled_mean_zero_inflated` | uniform 24.34% grid coverage in **every** member, so **no min-model rule exists to apply**; permanent water (493 always-flooded cells, mostly Caspian) excluded | none | two-tier, `higher_is_worse` — ranked against **its own** 2020s, which is also the shared reference |
| `flood-3b-40yr` | 27 (same 6 models) | `pooled_mean_zero_inflated` | as above | none | two-tier, `higher_is_worse` — ranked against the **flopros** 2020s |
| `flood-3b-none` | **32** (7 models — CWATM publishes the unprotected field only) | `pooled_mean_zero_inflated` — taken for estimator consistency, **not** because its own median failed | as above | none | two-tier, `higher_is_worse` — ranked against the **flopros** 2020s; compresses at the top |
| `conifer-npp` | — | `pooled_median` | 2% cover presence mask | none | **`higher_is_better`** — percentile is already inverted in the file |
| `csoil` | **17** (4 models × {2,5,5,5} GCMs) | `pooled_median` — branch 4 declined on measurement | union of finite cells, 71,251; `isfinite` IS a footprint here | none — split-half r=0.992 | single-tier, **`higher_is_better`** — inverted in the file |
| the nine threshold rungs (`heatdays-*`, `tropicalnights-*`, `frostdays-*`, `icedays-id`) | **12** (12 GCMs × **no impact model**) | `pooled_mean_zero_inflated` — **all nine**, on measurement | **`landseamask_no-ant.nc`, 65,797 cells — ANTARCTICA EXCLUDED**; the counts are finite over the whole globe including ocean, so `isfinite` is **not** a mask here | none — 120 draws per cell-decade on already-coherent bias-adjusted forcing | two-tier, `higher_is_worse` (frost included — user decision) |

### Which slope to read

Measured on active cells only (either slope non-zero) at the final decade of each layer's
first scenario. Reproduce with `python scripts/generate_customer_delivery.py --measure-slopes`.

| layer | active cells | `sen==0` | sign agreement | read |
|---|---|---|---|---|
| `conifer-npp` | 25,603 | 0.021 | 0.879 | `sen_slope` |
| `csoil` | 70,910 | 0.068 | **0.703** | `sen_slope` — the only layer in the **ols-is-biased** regime; see below |
| `cropfailure-3b` | 39,878 | 1.000 | 0.000 | `ols_slope` |
| `cyclone` | 20,337 | 0.974 | 0.025 | `ols_slope` |
| `drought-2b` | 66,741 | 1.000 | 0.000 | `ols_slope` |
| `drought-3b` | 61,810 | 1.000 | 0.000 | `ols_slope` |
| `wildfire` | 64,039 | 0.740 | 0.226 | `ols_slope` |
| `sealevel-2b` | 10,542 | 1.000 | 0.000 | `ols_slope` |
| `heatwave-3b` | 67,171 | 1.000 | 0.000 | `ols_slope` — **but see the saturation caveat below; on this layer a near-zero slope is ambiguous** |
| `permafrost-3b` | 14,455 | 0.152 | 0.844 | `ols_slope` — Sen has *not* collapsed here and is still wrong; see below |
| `heatdays-hd30` | 59,974 | 0.290 | 0.710 | `ols_slope` — see the ladder note below; the whole family reads ols |
| `heatdays-hd35` | 55,732 | 0.437 | 0.563 | `ols_slope` |
| `heatdays-hd40` | 48,086 | 0.699 | 0.301 | `ols_slope` |
| `heatdays-hd45` | 29,899 | 0.854 | 0.146 | `ols_slope` |
| `tropicalnights-tr20` | 56,451 | 0.439 | 0.561 | `ols_slope` |
| `tropicalnights-tr25` | 45,445 | 0.525 | 0.475 | `ols_slope` |
| `icedays-id` | 40,439 | 0.182 | 0.817 | `ols_slope` |
| `frostdays-fd` | 49,267 | 0.209 | 0.790 | `ols_slope` |
| `frostdays-fdm10` | 40,632 | 0.182 | 0.818 | `ols_slope` |

Measured 2026-08-12; `cropfailure-3b` added 2026-08-13; `heatwave-3b` and `permafrost-3b`
added 2026-08-14; `csoil` added 2026-08-15; the nine threshold rungs added 2026-08-16.

**The threshold ladder is the first family where the choice was made from the FAILURE MODES
rather than from the `sen==0` share**, and it is worth reading before applying a share
threshold to any new layer. Four of these rungs (`ID`, `FDm10`, `FD`, `hd30`) sit at
0.18–0.29 — lower than `permafrost-3b` and far below every other `ols_slope` entry — so a
mechanical "Sen has collapsed" rule would have read `sen_slope` on them. It is still wrong,
for the reason that decides it:

- `ols_slope`'s documented failure mode **requires uneven member coverage**. Measured across
  all nine rungs and all three scenarios, `n_members` is **exactly 12 in every cell**. The
  bias has no mechanism here. (Contrast `csoil`, where CLASSIC contributes 2 GCMs against 5
  and ols is genuinely biased — that is what the `sen_slope` column exists for.)
- `sen_slope`'s failure mode **is** present on every rung, because a threshold count is tied
  at zero wherever the threshold is rarely crossed.

One estimator that cannot fail here beats a mixed family where the reader must remember which
rung takes which. **Mind the denominator**: these shares are over ACTIVE cells (either slope
non-zero), the convention this table uses. Over all FINITE cells they run far higher — `ID`
reads 0.542 rather than 0.182 — because every never-crossing cell has both slopes correctly
at 0 and is not a Sen failure. Measuring on the wrong denominator initially flipped four of
these nine recommendations.

**`csoil` is the layer the other column of the two-slope table was written for.** Every other
entry above reads `ols_slope` *because Sen collapsed*; `conifer-npp` and `csoil` read
`sen_slope` because it did not. But only `csoil` sits in the regime OUTPUT-SPEC names as
ols's own failure mode — a large between-member **level offset** (~68.7× the interannual SD)
with **uneven member coverage** (CLASSIC contributes 2 GCMs against 5 for each of the other
three, and the four models' 2020s medians span 5.70–16.82 kg m⁻²). ols absorbs those offsets
as trend, and the two estimators disagree in **sign** on the global mean throughout: ssp585
reads ols −0.0119 against sen +0.0081 at the 2050s, and ols −0.0475 against sen −0.0072 by
the 2090s. Read `sen_slope`.

Its **70.3% sign agreement is the lowest of any `sen_slope` layer** (`conifer-npp` is 87.9%),
so the two estimators disagree on nearly a third of active cells — which is precisely the
standing signal that a cell's trend is not robust. Do not quote a site-level trend from this
layer without checking both.

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
| `cropfailure-3b`, `drought-3b`, `drought-2b`, `wildfire`, `conifer-npp`, `sealevel-2b` | **none**, each declared with a reason |

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

**The three flood layers are NOT in this class**, and that is the whole reason they exist.
They publish absolute inundation, so a high value means "deep or frequent flooding here",
full stop. The relative-baseline alternative for this hazard — the 0.5° `floodedarea`
exposure flag — reads **0.000 across the Amazon floodplain**, which is a correct
departure-from-preindustrial answer and an unusable flood layer. Flagging these three as
relative-baseline would be its own false claim; `relative_baseline: false` is set
deliberately on all three, not left at a default.

### Per-layer specifics worth knowing before use

- **`csoil`** — **a zero-carbon cell scores as MAXIMUM risk.** `higher_is_better` maps a low
  stock to a high risk percentile, and a place that never held soil carbon holds the least of
  it: measured on ssp585, the Sahara reads 0.00 kg m⁻² at the **96th** risk percentile, the
  Taklamakan 0.04 at the 93rd, and the **Greenland ice sheet 0.00 at the 97th**. The ranking
  is arithmetically correct and the conclusion inverts — a desert has nothing to lose. Same
  class of trap as a relative baseline. Greenland is in the footprint at all only because two
  of four models write `0` over the ice sheet rather than NaN.
  **`classic` is natively 1°**, replicated 2×2 with a one-cell longitude offset — measured
  **100%** adjacent-cell ties on odd column pairs against **0%** on even; the other three are
  genuinely 0.5°. Visible on the per-member contact sheet, invisible to every statistic in a
  value-check table. **Model footprints differ sharply at high latitude**: `mc2-usfs` covers
  58,919 cells against `jules`'s 67,647, and only 12,417 cells ≥60°N against ~17,200 — so
  above 70°N mean ensemble support falls to **12.2 of 17 members**, weighted toward the two
  models that read highest there. **Do not call the fixed-CO₂ member's trend "muted"**: it is
  the largest relative loser (−4.37%, against −2.75% `lpjml`, −0.05% `mc2-usfs`, **+0.79%**
  `classic` — which gains).
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
- **`sealevel-2b`** — **screening resolution: a regional result, not a site result** (see
  *Spatial support* above; it carries the must-disclose `resolution_caveat`). Terrain is read
  at 15″ and the within-cell elevation **distribution** is used, not a cell mean — but sea
  level is one value per 0.5° cell, connectivity is decided at ~1.85 km, and the published
  result is one value per cell. It is the only layer built from **two sources**: ISIMIP2b regional sea level
  applied to GEBCO_2026 terrain, so it is not an ISIMIP output and has no ISIMIP variable
  name. Three things about it invert a reader's expectation. (1) Its `percentile` is an
  **absolute calibration** (0 m → 100, the 10 m LECZ bound → 1, between filled by the
  empirical elevation distribution), not the contract's percentile-of-score, so it is **not
  on the same axis** as any other layer's percentile and a Climate Score averaging them must
  say so. (2) Absolute exposure is systematically **LOW**, because GEBCO is a digital
  *surface* model carrying canopy and buildings — South Florida reads `0.00000` at every
  decade despite 99.9% of its land lying below 10 m, and correcting the class of bias is
  measured to *triple* exposure estimates. A zero means "no land below projected sea level in
  this DEM", never "no coastal risk". (3) **Defences are invisible** at 460 m, so the
  Netherlands — 33% of its coastal land below projected sea level by the 2090s — reads
  exactly like an undefended coast. Highest exposure globally is **Arctic coastal tundra**,
  which is land-area exposure with almost no assets behind it. Two connectivity levels are
  used deliberately: 2 m decides what can flood, 10 m decides LECZ pool membership, because
  one level cannot do both — at 10 m the Salton Sink connects through the Colorado delta, and
  at 2 m coastal plain at 5 m would grade like a mountain.
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

**The three flood layers are variants, which is a tighter relationship than siblings.** Same
model chain, same members, same grid, same years — they differ only in the flood-protection
assumption. They are three `layer_id`s because OUTPUT-SPEC has no protection dimension, not
because they are three products.

### The flood variants: what the protection level does, and the one number it changes

`flood-3b-flopros` applies the empirical FLOPROS protection standards and is the delivery
default. `flood-3b-40yr` assumes a uniform 40-year standard **everywhere**, including places
with no defences at all — a yardstick, not a description of anywhere real. `flood-3b-none`
is the raw climate signal.

- **A zero under `flopros` means PROTECTED TO STANDARD, not NO HAZARD.** Measured: the Rhine
  (NL/DE) box falls from 0.05909 unprotected to 0.00122, the Lower Mississippi from 0.14505
  to 0.00358. Protection is held constant into the future, so the layer also embeds a
  no-further-adaptation assumption.
- **`40yr` is not an ordered midpoint between the other two.** The ranking *reverses* by
  region — Amazon 0.01077 at `40yr` against 0.01607 at `flopros`, Ganges 0.06853 against
  0.05664 — because a uniform standard exceeds the empirical one in some countries and falls
  short in others. Spatial Spearman `none` vs `flopros` is 0.367: defences **reorder** the
  map, they do not scale it.
- **The ensembles differ.** `none` carries 32 members from 7 impact models, the protected
  variants 27 from 6, because CWATM publishes only the unprotected field. A `none` − `flopros`
  difference therefore mixes a defence effect with an ensemble change unless restricted to
  the common 27.
- **All three rank against the `flopros` 2020s baseline** so a percentile means the same
  thing across them, while every raw value stays the variant's own. The cost lands on
  `none`: its median percentile is ~90 with ~11% of cells at ≥99, against 50 and 1% on its
  own baseline. Real sites show it — Manaus 99.0, Dhaka 99.9, Rotterdam 99.9, and dry Phoenix
  still 91.8. **Judge `none` on its median, not its percentile.**
- **These are ABSOLUTE inundation, so they are NOT in the departure-from-local-baseline
  class above.** That is exactly why the product was chosen: the 0.5° `floodedarea` exposure
  flag scores departure from a preindustrial reference and reads **0.000 across the Amazon
  floodplain in all 45 of its members**, where this layer reads 0.11496 unprotected.
- `sen_slope` is not collapsed on `none` (34–38% zero, ~60% sign agreement) the way it is on
  the protected variants (97–99%). `ols_slope` is still the recommended read on all three —
  on `none` that is a judgement recorded in the registry, not a reading off the numbers.

### Withdrawn, superseded and unregistered

| | Status |
|---|---|
| `yield-sug-noirr`, `yield-sug-firr` | **WITHDRAWN 2026-08-11, upstream data defect.** ISIMIP2b LPJmL does not simulate cane in the cane belt (São Paulo, UP India, Queensland, Florida all sentinel-zero). Passed the contract and meant nothing. No ISIMIP source supports a scenario-bearing sugarcane layer. Refused by the loader via `blocked:` in the registry. |
| `burntarea` (2b `biomes`, rcp26/60/85) | **SUPERSEDED 2026-08-10** by the 3b layer. |
| `csoil-total` | **Rebuilt and registered 2026-08-15** — see the `csoil` rows above. The 2026-07-25 build's output had been lost from disk and its 12-member/3-model ensemble was wrong by omission (LPJmL missing). `process_csoil_soilcarbon.py` remains the OUTPUT-SPEC **reference implementation**. |
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


### QA sign-off: the threshold ladder, 2026-08-16

**Signed off by the user on 2026-08-16**, recorded as `qa_reviewed_on: '2026-08-16'` on all
nine rungs in `config/layer_registry.yaml`. Only one other layer carries a human QA date
(`heatwave-3b`, 2026-08-14), so a null elsewhere is a real gap, not a formality.

A date is only auditable if the evidence behind it is recorded, so this is what the review had
in front of it:

- **Contract**: `test_shared_baseline.py` on all nine — 2020s panel bit-identical across
  scenarios, `percentile` within [1,100] and correctly oriented, no slope finite where the
  median is NaN, `n_members` 12 everywhere. Re-run and re-passed after the two in-place
  attribute corrections.
- **Containment across the ladder** — `hd30 ≥ hd35 ≥ hd40 ≥ hd45`, `TR20 ≥ TR25`,
  `FD ≥ FDm10`, `FD ≥ ID`: **0 violations** in all three scenarios. No per-layer check can
  establish this, and a flipped operator or swapped threshold could not survive it.
- **Reference sites, all nine rungs** (GUARDRAILS §12) — Singapore 347 `hd30` / 7 `hd35` /
  365 `TR20`; Kuwait 94 `hd45`; Yakutsk 228 `FD`, 181 `FDm10`, 178 `ID`; Nairobi 15 `hd30` and
  0 `hd35` at 1,795 m; tropics 0 on every cold rung.
- **Per-member statistics and contact sheets** for `hd35` and `FD`, viewed. 1.31× and 1.09×
  spreads with the member ordering following climate sensitivity and **inverting** between the
  hot and cold rungs; no block structure, seams, hemisphere flips or mask errors; the Andes
  and the Tibetan Plateau resolve correctly on `FD`.
- **Cross-layer**: rank correlation and top-decile agreement against `heatwave-3b`
  (+0.554 / 47.2%), establishing the two as complements rather than duplicates.

Dashboards for all nine are at `reports/maps/{rung}/`; contact sheets at
`reports/contact_sheets/`. **The review is of the DATA, not of the processor** — a future rung
built by the same script does not inherit this date.

### The threshold ladder: nine rungs, one ingest, and the absolute counterpart to `heatwave-3b`

Built 2026-08-16 from ISIMIP3b bias-adjusted daily `tasmax`/`tasmin`, 12 GCMs × 3 SSPs.
~1.34 TB was streamed, sha512-verified against publisher sidecars, reduced to annual counts
and **deleted** — provenance is `data/interim/tasthresh/download_provenance.csv`, a declared
deviation from the `data/raw/` retention convention. The whole ladder came from **one** pass:
once a day is read, testing it against nine thresholds costs nothing, so a rung we do not
ship today needs no second 1.34 TB.

**Every rung is ABSOLUTE, which is the entire point.** `heatwave-3b` scores departure from
each cell's own preindustrial distribution and ranks Chicago above Delhi; these count
crossings of a fixed physical threshold. Measured over the 65,797 shared land cells on the
2020s panel, the two are **not** redundant and **not** interchangeable:

| | value |
|---|---|
| Spearman(`hd35`, `heatwave-3b`) | **+0.554** |
| Spearman(`hd30`, `heatwave-3b`) | +0.631 |
| Spearman(`FD`, `heatwave-3b`) | −0.662 |
| **top-decile agreement, `hd35`** | **47.2%** (3,105 of 6,580 cells) |

Screening on one gives a materially different site list than screening on the other. Delhi
has 130.6 hot days a year and sits at the 53rd percentile of `heatwave-3b`; Chicago has 10.0
and sits at the 75th. **Singapore is the trap in the other direction**: 4.8 hot days in the
2020s rising to 268.8 by ssp585 2090s — a 56× rise, because its temperature distribution is
narrow and sits just below 35 °C. A threshold count reads near zero right up until the
distribution crosses it, so **a low value here is not evidence of a safe margin.**

Per-rung, from the built files (ssp585; `zero%` is the share of land cells whose 2020s value
is exactly 0, the same quantity as `percentile_zero_fraction`):

| rung | slope to read | zero% | 2020s | 2090s | saturates? | sparse? |
|---|---|---|---|---|---|---|
| `hd30` | **`sen_slope`** | 11.9% | 107.4 | 152.7 | **yes, 11.8%** | no |
| `hd35` | `ols_slope` | 19.3% | 39.9 | 102.9 | no | no |
| `hd40` | `ols_slope` | 36.3% | 10.5 | 44.0 | no | no |
| `hd45` | `ols_slope` | 70.8% | 1.0 | 13.6 | no | **yes** |
| `TR20` | `ols_slope` | 19.2% | 83.5 | 129.5 | **yes, 11.4%** | no |
| `TR25` | `ols_slope` | 38.9% | 19.7 | 81.8 | no | no |
| `ID` | `ols_slope` | 39.0% | 78.3 | 54.4 | no | no |
| `FD` | **`sen_slope`** | 25.6% | 118.8 | 85.0 | no | no |
| `FDm10` | `ols_slope` | 38.7% | 70.6 | 41.4 | no | no |

**Only `hd30` and `FD` take `sen_slope`** — the two rungs most of the land actually crosses.
Everywhere else a majority of year-pairs are 0→0 and Theil-Sen collapses (`hd45`: 94.7%
exactly zero). The split is *within one family, from one processor, on one unit*, which is
why OUTPUT-SPEC says the slope is **measured, not inferred from `field_nature`**.

**Antarctica is excluded from the mask, and that decision is load-bearing.** The full
ISIMIP3b mask has 92,889 cells, `landseamask_no-ant.nc` has 65,797, and the 27,092-cell
difference carried almost all of the cold rungs' censoring:

| | full mask | no-ant | Antarctica alone |
|---|---|---|---|
| `FD` 2020s at the 364-day ceiling | 30.25% | **2.01%** | 98.85% |
| `ID` 2020s | 28.19% | **1.08%** | 94.04% |
| `FDm10` 2020s | 21.25% | **0.00%** | 72.83% |

Keeping it would censor ~30% of the `FD` pool — the regime where both estimators go to zero
and **agree**, so the disagreement rule stops warning — and fill the top of the frost ranking
with a continent holding no assets. The rungs that genuinely saturate are the **hot** ones.

**All nine take `pooled_mean_zero_inflated`, measured against what the median would publish.**
The median branch erases cells that do cross the threshold: 29,058 cells (44.2% of land) on
`hd45` 2090s, 20,014 (30.4%) on `hd35` 2020s, each with a strictly positive expected count.
Going the other way costs 1.6% where the median works (`FD` 2020s: median 106.0, mean 104.3).
One statistic for the whole ladder, because a customer reads `hd35` beside `FD`.

**Containment was verified as a dimensional check**, not a plausibility one: `hd30 ≥ hd35 ≥
hd40 ≥ hd45`, `TR20 ≥ TR25`, `FD ≥ FDm10`, and `FD ≥ ID` (if the day's maximum is below
freezing, its minimum is too) — **0 violations** across all three scenarios. A flipped
comparison operator or a swapped threshold could not survive that.

Two further facts that travel with every file: there is **no impact model in the chain**, so
the CI carries GCM spread and interannual variability only; and the **12 GCMs are not 12
independent models** — CNRM-CM6-1/CNRM-ESM2-1 share ARPEGE/NEMO/ISBA and KACE-1-0-G runs
UKESM1-0-LL's HadGEM3 atmosphere. Family pooling was tested by correlating each member's
residual from the ensemble mean and **rejected**: on `FD` the CNRM pair ranks 1 of 66 pairs
(+0.764) but on `hd35` it ranks 14 (+0.228), behind KACE × UKESM (+0.546). The duplication is
real, rung-dependent, and does not partition, so `n_models` counts GCMs and the
non-independence is declared instead of silently resolved.

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
