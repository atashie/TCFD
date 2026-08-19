# Landslide data — what exists, at what resolution, and what each option would cost us

**Date: 2026-08-19.** Availability review plus a build record. Sections 0–9 are an
inventory: what was found, what was measured, and the decisions each option depends on.
Section 10 records what was actually built on 2026-08-19 and the measurements that forced
its design. Section 11 lists what is still open.

Read this before proposing any landslide work. It exists so the options do not have to be
re-derived, and so the ones already rejected are rejected *with their numbers*.

---

## 0. What this review established

1. **ISIMIP publishes no landslide, mass-movement, slope or geotechnical output in any
   round.** Verified by directory listing 2026-08-18; receipt and scope in §1.
2. **But the structural barrier is weaker than for tornado, and the proof is external.**
   NGI built the global GIRI landslide projection *from ISIMIP3b bias-adjusted
   precipitation*. ISIMIP lacks the impact model, not the forcing (§2).
3. **Landslide hazard factorises**, and every product that collapses the factorisation pays
   for it: `hazard(x,t) = S(x) [terrain, static, fine] × T(x,t) [rainfall, time-varying,
   coarse]` (§8).
4. **Exactly one global product publishes a RATE with physical units** — the World Bank /
   GFDRR map produced by Arup. Everything else is an ordinal class, an event point, or a
   figure in a paper (§4).
5. **Only four studies project landslides globally under a scenario, and three of them
   publish no grids at all** (§6). The one that does — GIRI — uses a single GCM.
6. **The finest published resolution is 3 arcsec (~90 m) and it buys less than it appears
   to**, because the climate factor underneath it is ~0.5–1° in every existing product (§8).
7. **The licence looked like the binding constraint and was resolved by asking.** Two
   publisher-side catalogue records for the Arup map disagree (CC BY-NC 4.0 vs CC-BY-4.0)
   and the project report states neither (§4.1). Cleared for our limited commercial use by
   user determination 2026-08-19, with **attribution required**.

---

## 1. ISIMIP — the receipt

Enumerated 2026-08-18, serially with 3 s spacing, **zero empty listings** (an empty listing
is a failure signal, not an absence — GUARDRAILS §11):

```
https://files.isimip.org/{ISIMIP2a,ISIMIP2b,ISIMIP3a,ISIMIP3b}/
https://files.isimip.org/{ISIMIP2a,ISIMIP2b,ISIMIP3a,ISIMIP3b}/OutputData/
https://files.isimip.org/{ISIMIP2a,ISIMIP2b,ISIMIP3a,ISIMIP3b}/DerivedOutputData/
https://files.isimip.org/ISIMIP{3a,3b}/InputData/  and their geo_conditions/
https://files.isimip.org/ISIMIP3b/SecondaryInputData/
```

Sector rosters found:

| Round | Sectors |
|---|---|
| 2a | agriculture biodiversity biomes forestry lakes_global marine-fishery_{global,regional} permafrost water_{global,regional} |
| 2b | as 2a **plus** health, lakes_local |
| 3a | agriculture biomes fire lakes_{global,local} peat permafrost water_global |
| 3b | agriculture biomes fire lakes_{global,local} marine-fishery_global permafrost water_global |

`DerivedOutputData/` publications: 2a `Zimmer2023`; 2b `Lange2020`, `Zimmer2023`; 3a
`Quesada-Chacon2026`, `Zantout2025`; 3b `Heinicke2026`, `Jaegermeyr2021`, `TipESM2025`,
`Zantout2025`. None is a mass-movement product.

`InputData/geo_conditions/` (3a and 3b, identical): `countrymasks fishmip_regions lakes
landseamask ocean river_routing soil`. **No slope, no elevation, no lithology.** This is the
substantive finding — the terrain factor is absent from the archive, not merely the impact
model.

**Scope of this negative.** NOT listed, so NOT covered: individual model directories under
each sector, the `InputData/climate/` subtrees, `composition_atmosphere/`, `socioeconomic/`,
and every `InputData` root in 2a/2b. No impact sector models slope stability, but that is
reasoning rather than enumeration.

---

## 2. Why ISIMIP cannot carry it — and the one thing that qualifies that

Unlike tornado, where the barrier is absolute (every tornado-environment index needs the
atmosphere in the vertical, and the bias-adjusted forcing publishes 11 surface variables and
nothing aloft), the landslide barrier is **the terrain factor, which is not climate data at
all**. Slope, lithology and soil depth do not change over a projection horizon and do not
need to come from ISIMIP.

The proof that this is a soft barrier: **the GIRI global landslide model is driven by
ISIMIP3b.** Palau et al. (2023, NGI) take 24 h rainfall intensity from W5E5 for 1979–2016
and from **IPSL-CM6A-LR under SSP126 and SSP585 for 2061–2100, via ISIMIP3b**, and combine it
with a terrain susceptibility index built from open global datasets.

So the correct statement is not "ISIMIP cannot express landslide" but: *ISIMIP supplies the
trigger and nothing else; the terrain factor must come from outside.* That is a materially
different conclusion from the tornado case and it is why §8 is worth reading.

---

## 3. Class 1 — event inventories (the observed record)

| Dataset | Coverage | Spatial | Temporal | Access |
|---|---|---|---|---|
| **NASA Global Landslide Catalog (GLC)** | global, rainfall-triggered only | points + qualitative accuracy radius (km) | **2007–2018** | open |
| **COOLR** = GLC + Landslide Reporter Catalog + contributed | global | points | 2007–ongoing | open |
| **Global Fatal Landslide Database** (Froude & Petley) | global, **non-seismic fatal only** | points | **2004–2016** (4,862 events, 55,997 deaths); online release to 2017 | open |
| **EM-DAT** (CRED) | threshold-based disasters | **country level** | 1900–present | registration |
| **DesInventar** (UNDRR) | country-by-country | **municipal level** | varies | open |
| **Tanyaş et al. / USGS EQIL database v2.0** | global, **earthquake-triggered** | 66 digital polygon inventories, 363 triggering earthquakes | 1900s–2021 | open (ScienceBase) |
| **HR-GLDD**, Sen12Landslides | 10–15 physiographic regions | **3 m** PlanetScope / 10 m Sentinel | event-based | open (Zenodo) — ML *training* sets |

**The governing property of this entire class.** Every global catalog is a *minimum* count.
NASA's own decadal analysis of the GLC finds reporting correlates positively with national
GDP and is biased toward English-language media. Absence of records is not absence of
landslides — and every global ML susceptibility model in §4 and §6 trains on exactly these
records, so the bias propagates into the gridded products.

---

## 4. Class 2 — global gridded hazard and susceptibility (present climate)

| Product | Quantity | Spatial | Temporal | Access |
|---|---|---|---|---|
| **NASA Global Landslide Susceptibility Map** (Stanley & Kirschbaum 2017; 2023 update) | 5-class susceptibility (heuristic fuzzy: slope, faults, geology, forest loss, roads) | **1 km (30″)** | static | open |
| **LHASA v1.1 nowcast** | categorical daily hazard | 1 km, 60°N–60°S | 2000-06-14 → 2020-12-31 | GES DISC |
| **LHASA v2.0** (XGBoost, probabilistic) | daily landslide probability | **1 km (0.00833°)**, 60°N–60°S | archive 2015-04-03 → 2021-02-10; operational NRT ~5 h latency | GES DISC, `10.5067/8VKQDQFFOTS3` |
| **Emberson et al. 2020** | population / road exposure | 1 km | 2001–2019 | NHESS, open |
| **Felsberg et al. 2022** | susceptibility **with explicit uncertainty** | **36 km** (satellite soil-moisture grid) | static | NHESS, open |
| **World Bank/GFDRR–Arup Global Landslide Hazard Map** | **annual frequency of significant landslides per km²** | **1 km**, land 60°S–72°N | rainfall-triggered median & mean **1980–2018**; separate earthquake layer | COG, open, licence disputed (§4.1) |
| **GIRI landslide model** (NGI/CDRI/UNDRR) | 5-class susceptibility + scenario probabilistic hazard, rainfall **and** earthquake | **3″ (~90 m)** | baseline W5E5 1979–2016 **+ SSP126/SSP585 2061–2100** | GIRI portal, registration |

### 4.1 The Arup map — the only rate, and its licence problem

This is the only global landslide product with physical units, which is the only reason it
can be aggregated to a coarser grid at all: the areal mean of a per-km² frequency is
arithmetic, whereas the mean of an ordinal susceptibility class is not a quantity.

Files (all under `https://datacatalogfiles.worldbank.org/ddh-published/0037584/`):

| File | What |
|---|---|
| `DR0045419/LS_RF_Median_1980-2018_COG.tif` | rainfall trigger, **median** 1980–2018, 124 MB — *ingested* |
| `DR0045418/LS_RF_Mean_1980-2018_COG.tif` | rainfall trigger, mean |
| `DR0045416/ls_eq_tiled.tif` | earthquake trigger |
| `DR0045417/LS_TH_COG.tif` | "TH ranks" — an ordinal |
| `DR0045411/global-landslide-hazard-map-report.pdf` | method report, 113 pp, 121 MB |

**The publisher's own licence records are inconsistent** — cleared for our use 2026-08-19,
recorded here because anyone re-deriving the provenance will hit it:

| Record | States |
|---|---|
| `datacatalog.worldbank.org` dataset 0037584 | Creative Commons Attribution-**Non Commercial** 4.0 |
| `energydata.info` CKAN mirror | **CC-BY-4.0** |
| The Arup project report itself | neither — no licence string in front or back matter (checked 2026-08-19) |

NC and BY are materially different for a commercial climate-risk product. Resolved by user
determination on 2026-08-19: cleared for our limited commercial use cases, **with attribution
to World Bank / GFDRR and Arup required wherever a value is published**. The layer built in
§10 is `status: preferred`.

---

## 5. Class 3 — national and government inventories

Highest quality, least comparable. Representative, not exhaustive:

| Source | Contents |
|---|---|
| **USGS** *Landslide Inventories across the United States* v2.0 (2022) | >570,000 points and polygons, federal + state surveys, per-feature location-confidence attribute |
| **BGS** National Landslide Database (Great Britain) | >17,000 records, continuously updated |
| **ISPRA IFFI** (Italy) | national official inventory, online since 2005 |
| **Preliminary Canadian Landslide Database** v13.0 | ~28,000 point entries (Zenodo) |

Japan (NIED), New Zealand (GNS), Australia (GA) and China hold national inventories that
this review did **not** verify — UNVERIFIED, not absent. These differ in mapping method,
epoch, minimum size and attribute schema, so they cannot be pooled into a globally
consistent frequency without the kind of explicit scoring system Tanyaş et al. built.

---

## 6. Class 4 — long-term projections (the only class that answers a scenario question)

### Global

| Study | Forcing | Scenarios | Horizons | Grid | Data published? |
|---|---|---|---|---|---|
| **GIRI / Palau et al. 2023 (NGI)** | **ISIMIP3b, IPSL-CM6A-LR — one GCM** | SSP126, SSP585 | **2061–2100** | **3″ (~90 m)** | Yes — GIRI portal (registration). *The only global projection distributed as a data product.* |
| **Duan et al. 2025, *Geoscience Frontiers*** | CMIP6 precip + static predictors, ML ensemble (AUC 0.971) | SSP126/245/370/585 | 2021–2100 vs **2001–2020** | 0.15° | **No** |
| **Duan et al. 2026, *Comms Earth & Env* 7:57** | **13 CMIP6 GCMs**, 7 ML models, predictor = max 3-day precip | 4 SSPs | vs **1981–2020**; NF 2021–2060, FF 2061–2100 | **0.15°** | **No** — statement names only ESGF + OpenStreetMap inputs |
| **Wang et al. 2023, IJDRS** | multiple GCMs | SSP2-4.5 | 2031–2060, 2066–2095 vs 1971–2000 | country-aggregated | No |

Headlines: global susceptibility **+~1%** by 2081–2100 SSP5-8.5 (Duan 2025); global road risk
**+30.6%** by 2100 SSP5-8.5, +23.9% even under SSP1-2.6 (Duan 2026); landslide frequency
**+7%/+10%** with casualty risk **+140%/+160%** (Wang 2023).

### Regional (higher fidelity, no global coverage)

- **Stanley et al. 2024, *Earth's Future*** — High Mountain Asia Landslide Hazard Indicator,
  SSP2-4.5 and SSP5-8.5 to 2100. Largest absolute increase in the Central Himalaya, largest
  relative on the Tibetan Plateau.
- **Kirschbaum et al. 2020, *GRL*** — HMA, 2061–2100 vs 1961–2000, **+30–70%** on the
  China–Nepal border.
- **Wang et al. 2025** (EGUsphere) — China, random forest, annual scale, to 2076–2100.
- **Ullah et al. 2024, *ERL*** — Pakistan Hindukush, **30 m**, 4 GCMs, SSP126/245/585,
  2040/2070/2100, *including* land-use change.
- **Roman Quintero et al. 2025** (EGUsphere) — S. Italy, RCP4.5/8.5 to 2070: **more
  landslides under a drier climate**, because rainfall re-timing keeps antecedent soil
  moisture high at the trigger moment. See §11 — this is the mechanism every global product
  is blind to.
- **Gariano & Guzzetti 2016, *Earth-Science Reviews*** — the standing review. Its finding
  still holds: the projection literature is overwhelmingly European, catchment-scale and
  physically based.

---

## 7. The Duan 0.15° products, reviewed in detail

Reviewed 2026-08-19 from the full text, Supplementary Information and the **public peer
review file**. Recorded here because "0.15° global landslide projection" sounds like exactly
what we want and is not.

- **Nothing is downloadable.** Data availability names inputs only; there is **no code
  availability statement**; no Zenodo/figshare deposit exists. Both papers are **CC
  BY-NC-ND 4.0** — the CEE licence block explicitly forbids sharing adapted material.
- **0.15° is presentation, not information.** The only time-varying predictor is ~1° CMIP6
  M3DP put through **bilinear interpolation**. No downscaling, no bias adjustment. Reviewer
  #1 said so and **recommended rejection**; the authors' response is a computational-cost
  argument. Meanwhile lithology enters at **0.5°** — coarser than the output grid.
- **AUC 0.9892 is not measuring skill at the useful task.** Negatives are drawn from
  anywhere on land ≥0.25° from a recorded event, so the classifier separates "steep, wet,
  well-reported terrain" from "random global land" — trivial with DEM, slope and lithology
  as inputs. The split is **random-stratified, not spatially blocked**, so there is direct
  autocorrelation leakage. 50/50 prevalence against a true ~1% means the output is an
  **uncalibrated relative score**, not a probability.
- **The income finding is partly a formula artifact.** V = E(1 − √C) puts GDP on both sides:
  road value grows ~linearly with wealth inside E while coping capacity enters as **√C**,
  which flattens. HIC > LIC is close to structurally guaranteed. The vulnerability term was
  **added during revision**; the submitted version had none.
- **Ensemble mean before a nonlinear model**, and no uncertainty on the output at all.
- **Two discussion claims the model cannot support**: Malta's rise to #1 is attributed to
  sea-level rise, which is **not a model input**; and Malta is ~316 km² against a ~280 km²
  equatorial cell — roughly one grid cell, 253 m high point.

**Conclusion: do not reimplement this framework.** Its defects are structural, not
incidental.

---

## 8. Resolution — what it can and cannot buy

Landslide hazard factorises into a static terrain term and a time-varying rainfall term.
Resolution behaves completely differently in each.

| | terrain factor S(x) | trigger factor T(x,t) |
|---|---|---|
| available resolution | 3″ (GIRI), 30″ (NASA, Arup) | 0.5° ISIMIP3b, ~1° raw CMIP6 |
| changes over horizon? | no | yes — this is the entire climate signal |

**Consequence: below ~0.5°, every bit of variation comes from terrain, not from climate.**
That is not automatically dishonest — it is the difference between coarse *support* and
coarse *inputs* that CLAUDE.md already draws for `sealevel-2b` (reads terrain at 15″,
publishes on a coarse grid, uses the within-cell distribution). Interpolating the *driver*
to a fine grid, as Duan does, invents information. Multiplying a coarse driver by a
genuinely fine static field does not.

**Storage, if a finer grid is ever wanted** (float32, full observational contract):

| Grid | Cells | One 2-D field | Note |
|---|---|---|---|
| 0.5° (most layers) | 259,200 | 1.0 MB | |
| **0.25° (built, §10)** | 1,036,800 | 4.1 MB | matches the tornado layer |
| 0.1° | 6.48 M | 26 MB | ~11 km; still ≫ climate resolution |
| 0.01° | 648 M | **2.6 GB** | breaks `generate_maps.py` browser limits; a full projected layer would be ~560 GB, the entire disk budget |

The recommendation on record (2026-08-18) is that resolution beyond ~0.1° is better spent on
**within-cell distribution statistics** than on more cells.

---

## 9. How each option sits against our output contract

| Option | Contract | Blocker |
|---|---|---|
| **Arup rainfall-triggered map — BUILT AND SHIPPED (§10)** | `observational-historical-v1` — no decade, no slopes, no ensemble | historical only; attribution required |
| NASA susceptibility 1 km | same, but an **ordinal** — no units, aggregation undefined | not a rate |
| GIRI 3″ + SSP126/585 | could carry a 2-period decadal axis | 1 GCM → `n_members`=1, registration, class-valued |
| Duan 0.15° | — | no data exists to ingest (§7) |
| **Build S × T ourselves** | full OUTPUT-SPEC decadal contract: 5 GCMs × 3 SSPs, real `lower_ci`/`upper_ci`/`n_members` | needs the ISIMIP3b daily `pr` harvest |

---

## 10. Build record — global historical landslide hazard at 0.25°, 2026-08-19

Built: `data/processed/landslide-arup_rf-median_hist/landslide-arup_observed_processed.nc`
Scripts: `scripts/download_landslide_arup.py`, `scripts/process_landslide_arup.py`
Verifier: `python3 scripts/test_observational_baseline.py …` → **17/17 PASS**
QA maps: `scripts/generate_landslide_qa.py` → `reports/maps/landslide/landslide-qa.html`

### 10.1 Zero is ambiguous in the source

The COG declares **no nodata value** and writes exact `0.0` for ocean *and* for flat land —
Pacific, Atlantic, Sahara, Amazon, Netherlands and Greenland all read `0.000000`. For a
landslide field, 0 on land is a legitimate result; 0 on water is not a result. Cells are
published only where inside the source extent **and** (ISIMIP3b land, upsampled 0.5→0.25°,
**OR** carrying ≥1 hazard-bearing pixel — so a coarse coastline cannot mask away ground the
source itself modelled as hazardous). Published mask: **247,271 cells**, of which **89,871
(36.35%) occupied**.

### 10.2 The statistic — four branches measured, three rejected

Over **occupied cells** (≥1 hazard-bearing pixel, the most favourable denominator, n=89,871):

| Branch | median==0 | q25==q75 | lower_ci<0 |
|---|---|---|---|
| A quantiles over all 900 native pixels | **78.48%** | **64.73%** | — |
| B areal mean ± 1 SD (`pooled_mean_zero_inflated`) | 0.00% | — | **88.55%** |
| C quartiles of 5×5 sub-block (0.05°) means | 49.38% | 33.38% | — |
| **D quantiles over hazard-bearing pixels — adopted** | **0.00%** | **5.70%** | **0.00%** |

- **A** is the documented zero-inflation failure in a new setting: three of four published
  variables identically 0 across two-thirds of the occupied map, in a file that still passes
  a structural check.
- **B** is what the shared module does for `let`, but a frequency cannot be negative and
  `mean − SD` < 0 in 88.55% of occupied cells.
- **C** was built specifically to rescue B with a non-negative bracketing interval. **It
  cannot be rescued**: the within-cell field is right-skewed, so the areal mean exceeds its
  own sub-block upper quartile in **42.48%** of cells. This is a property of the
  distribution, not of the block size — *do not re-attempt with a different block size.*
- **D** takes all three slots from **one** distribution, so the triple orders without
  clipping. (The tornado build's first attempt mixed an MLE central value with posterior
  quartiles and had to clip; that is the trap this avoids.)

### 10.3 The percentile ranks on a different variable, deliberately

`percentile` ranks on `areal_mean_rate`, not on `median`. Spearman rank correlation between
the two orderings is **0.34** over occupied cells, so the choice changes essentially every
score, and the reference sites decide it: the conditional median puts the **Apennines** at
58.9 against 74.2, and gives **Cairo** 5.3 off a single pixel covering 0.1% of the cell
against 1.4. The conditional median discards *extent*, which is first-order for exposure.

Consequence a reader must be told, and which is in the file's own `percentile` note: **a
cell can show `median` = 0 with a high percentile** — most of it is flat, the cell as a whole
still carries substantial landslide activity.

The shared verifier's two-tier check is keyed on a variable named `n_events`. This layer has
no event count, so the count is published as `n_hazard_pixels` and that check does not fire;
the equivalent assertions are made in the processor and printed (all pass).

### 10.4 The reference-site check

| Site | median | areal mean | hazard frac | **percentile** |
|---|---|---|---|---|
| Baguio PH | 0.22788 | 0.229637 | 0.999 | **99.996** |
| Wenchuan CN | 0.02203 | 0.021740 | 1.000 | 98.862 |
| Medellín CO | 0.03416 | 0.020437 | 0.726 | 98.616 |
| Shimla IN | 0.02528 | 0.019195 | 0.994 | 98.277 |
| Bergen NO | 0.00957 | 0.011918 | 0.576 | 94.325 |
| Cusco PE | 0.02000 | 0.010989 | 0.878 | 93.551 |
| Kathmandu NP | 0.02139 | 0.009102 | 0.648 | 91.550 |
| Apennines IT | 0.00426 | 0.005163 | 0.579 | 84.163 |
| Freetown SL | 0.00967 | 0.004269 | 0.270 | 81.248 |
| Rio de Janeiro BR | 0.00545 | 0.001811 | 0.158 | 67.905 |
| *Amsterdam NL* | 0.00000 | 0.000000 | 0.000 | *1.000* |
| *Des Moines IA* | 0.00000 | 0.000000 | 0.000 | *1.000* |
| *Cairo EG* | 0.00312 | 0.000003 | 0.001 | *2.437* |
| *mid-Pacific / McMurdo* | NaN | NaN | NaN | *NaN* |

Units are landslides per km² per year. Rio at 67.9 is the one to watch: a coastal cell that
is 84% water/flat, which is a resolution artifact rather than a modelling error, and exactly
what §8 predicts at 28 km.

---

## 11. Still open

1. ~~The Arup licence (§4.1).~~ **CLOSED 2026-08-19** — cleared by user determination for
   our limited commercial use. The layer is `status: preferred`. The publisher's own records
   remain inconsistent, which is recorded in the file's `source_licence`; **attribution to
   World Bank / GFDRR and Arup is required wherever a value is published**, carried in the
   file's `attribution_required` attribute and in the registry `delivery_note`.
2. **`qa_reviewed_on` is null — this is now the only thing between the layer and a delivery.**
   QA maps exist: `python3 scripts/generate_landslide_qa.py` writes
   `reports/maps/landslide/landslide-qa.html` (21.7 MB, four global panels). Read them, then set
   the date in `config/layer_registry.yaml`. The page's first review item is the deliberate
   `median`/`percentile` divergence, because that is the design decision most likely to be
   wrong.
3. **Whether to build S × T ourselves** (§9, and the 2026-08-18 recommendation). Would give
   a real scenario axis and a genuine 5-member ensemble — the one axis where we would be
   better than every published product, since GIRI uses 1 GCM and Duan discards its spread.
   Prerequisites: harvest ISIMIP3b daily `pr` (0.5°, 5 GCMs × 3 SSPs), and pick the trigger
   threshold — **compute the whole ladder in one pass** rather than pre-committing.
4. **Absolute vs locally-normalised rainfall threshold.** GIRI normalises 24 h intensity
   against the local annual-max distribution, because 100 mm/day triggers slides in a dry
   climate and nothing in a monsoon one. If we do the same the layer becomes
   `relative_baseline: true` — "unusual rainfall for this place" — and must carry the
   must-disclose caveat, with its siblings checked. An absolute mm ladder avoids that but
   under-reads arid steep terrain. Recommendation on record: compute both, decide on the
   measurements.
5. **Antecedent soil moisture — the one genuine advance available to us.** Every global
   product (Arup, GIRI, Duan, Wang) uses a precipitation index alone and is therefore blind
   to the Roman Quintero mechanism, where a *drier* climate produces *more* landslides
   through rainfall re-timing. **`soilmoist` is confirmed present in ISIMIP3b
   `water_global`** across 10 models (verified 2026-08-19). Its **cadence is not yet
   checked** — the file counts (3 per model×GCM against 30 for `qtot`/`dis`) suggest monthly,
   which would be useless for a 3-day antecedent window. Check this before promising it.
6. **Earthquake-triggered landslides are excluded** and that is a scope decision, not an
   omission: seismic triggering is stationary under every pathway. GIRI and Arup both
   publish a seismic layer, so a reader comparing totals will find ours lower.

---

## 12. Repo updates this review made

- `config/isimip_search_catalog.yaml` → `negative_results.landslide`, with the §1 receipt.
- `config/hazard_taxonomy.yaml` → the `landslide` entry, whose `blocker` said "Outside
  ISIMIP's published impact sectors. Requires terrain and geotechnical data." That is still
  true of ISIMIP and was misleading about what exists outside it.
- `config/layer_registry.yaml` → `landslide-arup`, `status: evaluation`.
- `DATASET-ATTRIBUTES.md` → index line.
- `WORKFLOW-ISSUES.md` → dated entry.

---

## Sources

ISIMIP file server (enumerated 2026-08-18) ·
[World Bank Global landslide hazard map](https://datacatalog.worldbank.org/search/dataset/0037584/global-landslide-hazard-map) ·
[energydata.info mirror](https://energydata.info/dataset/global-landslide-hazard-map) ·
[GIRI landslide model, Palau et al. 2023](https://giri.unepgrid.ch/sites/default/files/2023-06/20230615-NGI_manuscript_GIRI_landlside_hazard_model.pdf) ·
[GIRI portal](https://giri.unepgrid.ch/node/154) ·
[NASA GLC reporting-bias analysis](https://ntrs.nasa.gov/api/citations/20230001884/downloads/glc_manuscript.pdf) ·
[COOLR](https://un-spider.org/links-and-resources/data-sources/cooperative-open-online-landslide-repository-coolr-nasa) ·
[Froude & Petley 2018](https://nhess.copernicus.org/articles/18/2161/2018/) ·
[Tanyaş et al. 2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JF004236) ·
[USGS EQIL database v2.0](https://www.sciencebase.gov/catalog/item/614512b3d34e0df5fb95b5f9) ·
[Stanley & Kirschbaum 2017](https://link.springer.com/article/10.1007/s11069-017-2757-y) ·
[LHASA v2 at GES DISC](https://www.earthdata.nasa.gov/data/catalog/ges-disc-global-landslide-nowcast-2.0.0) ·
[Emberson et al. 2020](https://nhess.copernicus.org/articles/20/3413/2020/) ·
[Felsberg et al. 2022](https://nhess.copernicus.org/articles/22/3063/2022/) ·
[HR-GLDD](https://essd.copernicus.org/articles/15/3283/2023/) ·
[USGS US Landslide Inventory v2.0](https://www.usgs.gov/data/landslide-inventories-across-united-states-ver-20-june-2022) ·
[BGS National Landslide Database](https://www.bgs.ac.uk/datasets/national-landslide-database/) ·
[ISPRA IFFI](https://www.progettoiffi.isprambiente.it/?lang=en) ·
[Preliminary Canadian Landslide Database](https://zenodo.org/records/17219072) ·
[Duan et al. 2026 (CEE 7:57)](https://www.nature.com/articles/s43247-025-03073-8) ·
[its peer review file](https://static-content.springer.com/esm/art%3A10.1038%2Fs43247-025-03073-8/MediaObjects/43247_2025_3073_MOESM1_ESM.pdf) ·
[Duan et al. 2025 (Geoscience Frontiers)](https://journal.hep.com.cn/gsf/EN/10.1016/j.gsf.2025.102074) ·
[Wang et al. 2023 (IJDRS)](https://link.springer.com/article/10.1007/s13753-023-00514-w) ·
[Stanley et al. 2024 (Earth's Future)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023EF004325) ·
[Kirschbaum et al. 2020 (GRL)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019GL085347) ·
[Ullah et al. 2024 (ERL)](https://iopscience.iop.org/article/10.1088/1748-9326/ad8a72) ·
[Roman Quintero et al. 2025](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-5826/) ·
[Gariano & Guzzetti 2016](https://geomorphology.irpi.cnr.it/publications/repository/public/journals/2016/landslides-in-a-changing-climate/@@download/file/1-s2.0-S0012825216302458-main.pdf) ·
[Özturk et al. 2022 (Nature)](https://www.nature.com/articles/d41586-022-02141-9)
