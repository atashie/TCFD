# Water stress — where the work stands

**Date: 2026-08-17.** Working record covering the ISIMIP water-sector enumeration, the
metric design, the decisions taken, the external review, and the river-routing finding
that changed the plan. Companion to [water-stress-formulations.md](water-stress-formulations.md),
which is the design doc and will be rewritten to match this.


> **Revised 2026-08-17, after this record was first written.** Two scope decisions supersede
> parts of what follows; superseded passages are marked rather than deleted.
>
> - **Environmental flows are out of scope for this pass** — the requirement is set by
>   regional and legal context, not derivable from hydrology alone. Formula B is dropped,
>   taking Blockers 2 and 3 with it. Six metrics become four. This aligns the metric with
>   Aqueduct 4.0 baseline water stress and therefore ESRS E3, which adopts Aqueduct's
>   thresholds verbatim; it moves away from SDG 6.4.2, which is defined *with* the deduction.
> - **Routing is used for upstream accumulation only.** A cell's demand is the sum over every
>   cell draining through it. Basin-uniform aggregation is rejected: it would charge a
>   headwater cell with withdrawals occurring downstream of it. Downstream abstraction cannot
>   reduce upstream supply and the metric must not imply that it does. The basin identifier is
>   retained as a diagnostic grouping key only, never as the support of a published value.
> - **New risk, measurable up front:** DDM30 was built for WaterGAP. H08 and LPJmL5 route on
>   their own networks and their published discharge reflects those. Accumulate each model's
>   own `qtot` along DDM30 and compare cell-wise to its published `dis`; a model whose network
>   disagrees cannot use the routed metric.

---

## 1. Enumeration of the ISIMIP water sector

**Receipt.** 2026-08-16, serial listings with 3-retry, empty-listing-treated-as-failure.
All 109 `{MODEL}/{gcm}/future/` directories walked (47 in ISIMIP3b, 62 in ISIMIP2b),
35,977 files projected from the end of the filename. Zero EMPTY / FAILED / NO_FUTURE.
The sub-directory level was enumerated separately per model: only 2b `DBH` and 2b
`WaterGAP2` carry a `future_extended/`; DBH's holds picontrol `qtot` past 2100,
WaterGAP2's is empty.

### 1.1 Who publishes water use — ISIMIP3b

Monthly, one file per variable per scenario × soc, 2015–2100, `(time, lat, lon)` at 0.5°.

| model | GCMs | ptotww | ptotuse | atotww | atotuse | gw split |
|---|---|---|---|---|---|---|
| **H08** | 5 | yes | yes | yes | yes | — |
| **LPJmL5-7-10-fire** | 5 | yes | yes | yes | yes | — |
| **WaterGAP2-2e** | 5 | yes | yes | — | yes | yes (`…gw` for every sector) |
| **CWatM** | 5 | yes | yes | yes | — | yes (`atotusegw`, `atotwwgw`, `ptotwwgw`) |
| MIROC-INTEG-LAND | 5 | — | — | — | — | sectoral only (dom/ind/irr × a/p × ww/use) |

CLASSIC, JULES-ES-VN6P3, JULES-W2, VISIT and WEB-DHM-SG publish **no** water-use variable.

A four-model total-demand ensemble is available on `ptotww` and `ptotuse` only. `atotww`
gives three models, `atotuse` gives three — and they are **different** threes.

### 1.2 Units, read from the published `.json` sidecars

Every `.nc` has a sibling `.json` carrying the full NetCDF header. Measured 2026-08-16:

| variable | units | note |
|---|---|---|
| `ptotww`, `ptotuse`, `atotww`, `atotuse` | `kg m-2 s-1` | per-area flux, **not** a volume |
| `qtot`, `qr`, `qs`, `qsb`, `qg` | `kg m-2 s-1` | |
| `dis` | `m3 s-1` | already volumetric |
| `cellarea` | **`km2`** | published by CWatM, LPJmL5-7-10-fire, WaterGAP2-2e only |
| `contfrac` | `1` | "Continental Fraction of Grid Cell" |
| `soilmoist` | `kg m-2` | carries a `depth` dim whose length differs by model |

`qr` is *"Total groundwater recharge"*; `qg` is *"Groundwater Runoff"* — distinct
variables, both present. Note `scripts/config_water_variables.py` labels `qr` as
"Groundwater Runoff (Recharge)", which conflates the two; the file's own `long_name` is
authoritative.

`soilmoist` depth by model (gfdl-esm4): H08 1, WaterGAP2-2e 1, CWatM 3, JULES-ES-VN6P3 4,
JULES-W2 4, LPJmL5 5, MIROC-INTEG-LAND 13, CLASSIC 20; WEB-DHM-SG has no depth dim;
VISIT absent. It is not one variable across models.

### 1.3 The socioeconomic (`soc`) dimension — a catalog correction

`config/isimip_search_catalog.yaml:366` says the 3b water sector has "SIX soc variants per
scenario". That was probed on H08 and **it is H08-only**. Across all 25,375 projected 3b
rows the twelve `ssp{N}soc-{adapt,noadapt}+{image,magpie}` tokens carry 370 files each —
all H08 (5 GCM × 74 files).

| model | soc variants under SSP | demand time-varying |
|---|---|---|
| **H08** | 1850soc, 2015soc, `{ssp}soc-{adapt,noadapt}+{image,magpie}` | **yes** — 4 transient variants |
| LPJmL5-7-10-fire, WaterGAP2-2e, MIROC-INTEG-LAND | 1850soc, 2015soc, 2015soc-from-histsoc | no |
| CWatM | 2015soc-from-histsoc **only** | no |
| CLASSIC, JULES-W2 | 2015soc, 2015soc-from-histsoc | n/a |
| JULES-ES-VN6P3 | 2015soc-from-histsoc | n/a |
| VISIT, WEB-DHM-SG | 2015soc | n/a |

Consequences:

- Transient demand under SSPs exists in exactly one 3b model.
- **CWatM and H08 share no soc variant at all.** Any ensemble containing both necessarily
  mixes `2015soc` with `2015soc-from-histsoc` — the same fixed-2015 forward
  socioeconomics with a different historical spin-up. Precedent: the flood family, where
  all seven models were kept over a mixed soc (user decision 2026-08-14).
- **No `nosoc` in 3b** — verified over the complete soc token set. `1850soc` is the closest
  naturalised-demand run and exists for H08, LPJmL5, WaterGAP2-2e and MIROC — **not CWatM**.

`2015co2` sensitivity runs exist alongside `default` for CLASSIC, JULES-ES-VN6P3,
JULES-W2, LPJmL5-7-10-fire, VISIT and WEB-DHM-SG (mostly ssp585). LPJmL5 ssp585/2015soc
therefore has two sens variants of every demand file.

### 1.4 ISIMIP2b

GCMs are gfdl-esm2m, hadgem2-es, ipsl-cm5a-lr, miroc5 (4 each), except MPI-HM 3,
ORCHIDEE-DGVM 2, DBH 1. Monthly files are single-file 2006–2099.

- **Totals exist in one model only**: WaterGAP2-2c — `ptotww`, `ptotuse`, `atotuse`
  (no `atotww`), rcp26/60/85, **2005soc only**, demand fixed.
- **Transient socioeconomics** (`rcp26soc`, `rcp60soc`) exist for CWatM, H08, LPJmL and
  MATSIRO — but as **sectoral** variables only, never totals. Pairings are matched:
  rcp26↔rcp26soc, rcp60↔rcp60soc. **No `rcp85soc` anywhere** (verified 2026-08-16).
- **`nosoc`** exists for JULES-W1, ORCHIDEE, ORCHIDEE-DGVM and VIC — **none of which
  publishes any demand variable**. The naturalised run and the demand runs are disjoint
  model sets in 2b.
- **`WaterGAP2` (non-2c) is empty**: both `future/` and `future_extended/` return no files
  at all four GCMs (verified 2026-08-16 by direct listing). Superseded by WaterGAP2-2c.
  Relevant because `watergap2` appears in the 2b `led` drought ensemble list.
- 2b `dis` is monthly only for CWatM, MATSIRO and JULES-ES-55; everything else including
  WaterGAP2-2c is daily-only.

### 1.5 The demand forcing is published separately — and omits irrigation

`ISIMIP3b/InputData/socioeconomic/water_abstraction/` holds `1850soc`, `2015soc`,
`histsoc` and `ssp{126,370,585}soc-noadapt`. Annual, 0.5°, 2015–2100, already in
**m³/year**, and the SSP directories ship a `-modelstd` companion giving spread across
H08/WaterGAP2/CWatM.

But across **all six** directories the variables are `domww`, `domwc`, `indww`, `indwc`
only — **there is no irrigation field**. Irrigation is computed inside each hydrological
model from the land-use forcing, which is why `pirrww` is a model output.

This kills the tempting hybrid of pairing transient demand forcing with a deep supply
ensemble: the transient forcing covers only the two sectors that are roughly 30% of global
withdrawal and omits the dominant, most climate-sensitive one.

### 1.6 A publishing defect to remember

CWatM's 3b sidecars are misnamed. The directory holds 153 `.nc` files all carrying
`2015soc-from-histsoc`, but 123 of the 156 `.json` sidecars carry `2015soc`. A URL list
built from sidecar names 404s against the data. **Verify CWatM downloads by
Content-Length, not by fetching the sha512 at the matching path.**

### 1.7 Volumes

Measured per file: 219–266 MB for a 3b monthly 2015–2100 file (1032 × 360 × 720 float32,
deflated); 259 MB for a 2b 2006–2099 monthly file. Daily `dis` decade chunks are ~890 MB.

---

## 2. Choosing a single metric

### 2.1 The candidates

| | **A · WTA** | **B · SDG 6.4.2** | **C · Depletion** |
|---|---|---|---|
| Formula | `W / A` | `W / (A − EFR)` | `C / A` |
| Numerator | `ptotww` (or `atotww`) | same | `ptotuse` (or `atotuse`) |
| Models | 4 | 4 | 4 (`ptotuse`) / 3 (`atotuse`) |
| Members | 20 (4 × 5 GCM) | 20 | 20 |
| Bands transfer? | yes — 40/80% Aqueduct **and ESRS E3 verbatim** | yes — SDG's 25/50/75/100 | **no** — 40/80 are withdrawal-calibrated |
| Extra machinery | none | VMF from the same `dis`, no extra download | none |

**D · Standardized streamflow index (supply only).** No demand — `dis` departure from its
own baseline. Deepest ensemble available (6 models × 5 GCMs = 30 members) and it returns a
meaningful number in every cell on Earth. But it is hydrological drought, not water stress:
it sets `relative_baseline: true`, cannot carry the ESRS/CDP bands, and **cannot populate
`families.water-stress.covered_by`**.

**E · A blended 0–100 score** was rejected: the weights are unagreed and unmeasurable, it
averages an absolute field with a relative-baseline one, and it is the one construction a
reviewer cannot reproduce.

### 2.2 Why potential rather than actual withdrawal

`atot*` is what the model let the sector take *given supply*. A ratio built on it is
self-limiting by construction — it cannot exceed 1 — and would report a drying basin as
unstressed because everyone politely took less. `ptot*` is unconstrained demand. It is also
the version all four models publish.

### 2.3 Decisions taken

1. **Numerator: `ptotww`** — potential total withdrawal, all sectors.
2. **Denominator from the `1850soc` runs**, paired member-wise with `2015soc` demand.
3. ~~**The environmental-flow deduction (B) applies only to the river-flow denominator.**~~
   **Superseded** — environmental flows dropped entirely for this pass.
4. **Seasonal reduction is the highest-stress month** — the ratio computed in each of the
   12 months, maximum taken. Not the lowest-supply month.

### 2.4 Why decision 2 matters, and what it costs

River discharge under `2015soc` is flow **after** people have taken water out — reservoirs,
canals and intakes are all simulated. Dividing demand by it puts demand on top and also
subtracts it from the bottom: heavy withdrawal shrinks the denominator, which inflates the
ratio through a feedback rather than a measurement.

Taking supply from `1850soc` also made formula B *correct* rather than approximate. Pastor's
VMF method derives the environmental requirement from **naturalised** flow regimes. Computed
on `2015soc` discharge, reservoir operation would have flattened the seasonality — dams raise
dry-season flow — classifying fewer months as low-flow and reserving less water, biasing the
requirement downward exactly where it matters.

Costs, both of which belong in the file attributes:

- **CWatM drops out** (no `1850soc`), leaving three models / 15 members.
- **The denominator no longer responds to storage infrastructure.** A well-dammed basin and
  an undammed one with the same natural flow read identically. This metric cannot show the
  benefit of reservoirs.
- `1850soc` is preindustrial *conditions*, not zero-human — negligible withdrawal and no
  modern dams, but some land clearing, and preindustrial land use changes
  evapotranspiration and therefore runoff.

### 2.5 The ensemble

Verified: H08, LPJmL5-7-10-fire and WaterGAP2-2e all have `ptotww` monthly at `2015soc`,
and `qtot` + `dis` monthly at `1850soc`, for all of ssp126 / ssp370 / ssp585, five GCMs
each (gfdl-esm4, ipsl-cm6a-lr, mpi-esm1-2-hr, mri-esm2-0, ukesm1-0-ll). **15 members per
scenario.**

MIROC-INTEG-LAND can restore a fourth model (20 members) by summing
`pdomww + pindww + pirrww`; it has the `1850soc` supply and the `2015soc` sectoral demand,
5 GCMs, all three scenarios. See §4 item 10 for why the originally proposed acceptance
test for this was wrong.

---

## 3. The formulae

Cell `c`; member `m` = (impact model × GCM); calendar month `k ∈ 1…12`; year `y`.
`dt(y,k)` = actual seconds in that month — the ISIMIP time axis is `days since 1601-1-1`,
calendar `proleptic_gregorian`, so leap years are real. Baseline `B` = 2020–2029.

### 3.1 Unit conversion  *(SUPERSEDED — see §11)*

> The conversion and calendar stated in this section are **superseded by §11**: the flux is
> per continental area (`× contfrac/100`) and the WaterGAP2-2e calendar is `365_day`, not
> `proleptic_gregorian`. Retained as written because §11 documents how the first answer was
> reached and corrected. **§11 is authoritative.**


```
V [m³ s⁻¹]  =  flux [kg m⁻² s⁻¹] · Acell [km²] · 10⁶ [m² km⁻²] / 1000 [kg m⁻³]
            =  flux · Acell · 1000
```

The multiplier is **1000 · Acell(km²)**; ρ = 1000 kg m⁻³ is the only physical constant.
`dis` needs no area conversion.

`cellarea` and `contfrac` are published by CWatM, LPJmL5-7-10-fire and WaterGAP2-2e only.
H08 and MIROC publish neither, but all ten models are on the identical 0.5° grid, so one
model's grid serves the others **once the three published grids are confirmed identical**.

**Open and load-bearing:** whether the flux denominator is grid-cell area or continental
area, i.e. whether `contfrac` multiplies in. The decisive check exists — convert
WaterGAP2-2e's or CWatM's `pdomww` to m³/yr and compare cell-by-cell against
`InputData/socioeconomic/water_abstraction/2015soc/domww_2015soc_annual_2015_2100.nc`,
which is already in m³/year. Note this check is decisive for *those models and that
variable*; area normalisation must be verified per model and per variable.

### 3.2 Demand

```
W(c,y,k,m) = ptotww · 1000 · Acell(c) · dt(y,k)          [m³]
```

### 3.3 Availability

```
Q_runoff(c,y,k,m) = qtot · 1000 · Acell(c) · dt(y,k)     [m³]
Q_river(c,y,k,m)  = dis · dt(y,k)                        [m³]
```

Both from `1850soc`. **Do not add `qr` (groundwater recharge) to availability without
measuring first** — routed discharge already carries the baseflow that recharge feeds, so
`Q + recharge` double-counts. H08, LPJmL5, MIROC and WaterGAP2-2e publish `qs`, `qsb`,
`qg` and `qtot` together, so the water-balance check is directly available.

### 3.4 Environmental flow requirement — VMF  *(DROPPED — retained as a record)*

Pastor et al. 2014 (HESS 18:5041). Computed on the 2020s baseline and **frozen** — a
requirement recomputed per decade drifts with the climate and makes the index structurally
unable to detect drying.

```
MAF(c,m)   = baseline mean annual flow, duration-weighted        [m³ s⁻¹]
MMF(c,k,m) = baseline mean flow of calendar month k              [m³ s⁻¹]

           MMF < 0.40 · MAF                →  f = 0.60   (low-flow month)
f(c,k,m) = 0.40·MAF ≤ MMF ≤ 0.80·MAF       →  f = 0.45   (intermediate)
           MMF > 0.80 · MAF                →  f = 0.30   (high-flow month)

EFR(c,k,m) = f · MMF                                             [m³ s⁻¹]
```

The five constants are **0.40 and 0.80** (classification, as fractions of MAF) and
**0.60, 0.45, 0.30** (allocation, as fractions of that month's MMF). Sanity check: annual
EFR *volume* should be near **37% of annual flow volume**, with 25–46% spanning the method
family.

**EFR is per member `m`.** See §4 item 2 — the original plan pooled it and that was wrong.

### 3.5 Reference bands

Applied in the report, never baked into the values.

| scheme | thresholds |
|---|---|
| Aqueduct 4.0 / ESRS E3 verbatim | 0.40 high · 0.80 extremely high |
| SDG 6.4.2 | 0.25 · 0.50 · 0.75 · 1.00 |
| Criticality ratio (Alcamo, Vörösmarty) | 0.40 high · 0.80 extreme |

**Ratios are not capped at 1.** Capping would manufacture a censored field with the
`heatwave-3b` pathology.

---

## 4. External review (Codex, 2026-08-17)

Run against the plan plus `OUTPUT-SPEC.md`, `GUARDRAILS.md`, `decadal_stats.py` and
`layer_registry.yaml`. Full output at
`scratchpad/../tasks/bsrc19i7r.output`. The two findings that cite our own code and spec
were checked against the files and **both citations are accurate**.

### Blockers

**1 · Support mismatch — RESOLVED, see §5.** `W` was cell-local while `dis` is routed
from the whole upstream catchment. Downstream cells receive basin-scale supply against only
local demand, systematically understating stress along large rivers. Two further
consequences: the Aqueduct/ESRS/SDG bands were not merely "indicative" but unsupported;
and the claim that ratios above 1 indicate groundwater mining is **false** for a cell-local
construction — above 1 there just means the cell imports water from upstream, which is the
normal condition for most riverside cells. Our own design doc already called this "mismatched
support, not WTA in the literature's sense," and the plan proposed it anyway.

**2 · The pooled EFR is invalid.** *(MOOT — formula B dropped.)* The plan computed one `MMF`, `MAF`, class and absolute
EFR after pooling all models and GCMs, destroying member identity before a discontinuous
classification and an absolute subtraction. Worked example: two members with comparable
seasonality but baseline monthly flows of 10 and 100 m³/s give a pooled low-flow EFR of
0.6 × 55 = 33; subtracting 33 from the 10 m³/s member makes it environmentally bankrupt
when its own requirement would be 6. That is an ensemble artifact, not scarcity. EFR must
be computed **per impact-model × GCM member**, from that member's own shared,
scenario-averaged baseline climatology.

**3 · `Q − EFR` has a pole, and every treatment biases the statistics.** *(MOOT — formula B dropped. The verified code asymmetry below still applies to any infinity from a zero denominator, so item 7's precedence table is now the sole defence.)* Just above zero
the ratio becomes arbitrarily large, OLS is dominated by the near-pole year, and the
worst-month reduction almost guarantees selecting that month. Below zero, retaining the
ratio makes a more severe state look less risky. Masking removes precisely the most-stressed
observations and makes missingness state-dependent. Capping creates a censored field where
both slopes flatten together and agree.

This interacts badly with our implementation, **verified**: `decadal_stats.py:211` screens
for all-non-finite cells but the central statistic runs through `nanmedian`/`nanmean`,
which ignore NaN and **propagate ±inf**; the slope routine at line 392 filters with
`np.isfinite` and drops it. A single infinite ratio would therefore corrupt `median` and
the CI **while vanishing from both slopes**. `Q ≤ EFR` is a distinct
"environmental allocation exhausted" state and needs its own representation.

### High severity

**4 · "Demand is fixed at 2015 levels" was inferred from a filename and is probably false.**
`2015soc` fixes socioeconomic *inputs* — population, GDP, irrigated area. Irrigation
withdrawal is computed inside each model and responds to precipitation, evapotranspiration
and growing-season conditions, and irrigation is the majority of withdrawal. The assertion
that every trend is supply-driven does not follow from the soc token. This is a GUARDRAILS
§9 violation. Measure per member and sector: temporal variance and slope of `ptotww`,
whether non-irrigation sectors are fixed, and how much variability comes from irrigation.
Until then describe *socioeconomic conditions* as fixed, not demand.

**5 · The max-over-months reduction is not a stable seasonal metric.** After pooling it is
the median of 150 member-year maxima per decade — not the maximum month of the decadal
climatology, not one consistent calendar month, not a trend in a particular season. The
argmax month can switch between years and members; changes in monthly variance or one
unusually small denominator can create a trend even if typical monthly stress is unchanged.
It also biases low if any month is missing unless complete years are required. Measure the
argmax-month distribution, month-switch rate, the denominator at the maximum, and
sensitivity to requiring complete years.

**6 · The shared-baseline ordering was underspecified.** `OUTPUT-SPEC.md:226` requires each
member's 2020s window to be averaged across scenarios *first*, then pooled. Required order:
(i) construct the shared scenario-averaged baseline per member; (ii) derive that member's
MAF, MMF, class and EFR; (iii) compute each scenario's annual ratios using that member's
EFR; (iv) build the shared 2020s panel with the contract's cross-scenario averaging.

**7 · Zero denominators are unhandled.** The plan covered only the river mask and
`Q_river − EFR ≤ 0`. It omitted `qtot = 0`; `dis = 0`; `W = 0, Q = 0` which gives `0/0 = NaN`
rather than the zero the Sahara reference-site check demands; `W > 0, Q = 0` which gives
infinity; future monthly zeros inside a cell that passed a *baseline* river mask; and
negative or fill-value leakage. A precedence table is needed, with `W == 0` mapping
explicitly to zero **before** denominator handling.

### Medium severity

**8 · Arithmetic tightening.** For monthly A1 and A2 `dt` cancels exactly; for A1 area
cancels too, so `contfrac` cannot affect A1 provided `ptotww` and `qtot` share an area
basis. Using `dt_bar(k)` for EFR while using actual `dt(y,k)` for discharge introduces a
needless leap-year inconsistency — store EFR as m³/s and multiply both terms by the actual
month duration. `MAF` must be duration-weighted; a plain mean of 12 monthly rates is not an
annual mean. The "37% of mean annual flow" check must compare annual *volume* to annual
*volume*. And unit attributes alone do not establish that monthly values are monthly means —
verify `cell_methods` and time bounds.

**9 · The zero-inflated mean branch is dangerous for a ratio.** That branch was justified
for bounded expected-exposure fractions. A water-stress ratio has a heavy right tail and a
potential pole; mean ± SD can make the central estimate depend on a few near-zero
denominators and can produce a negative lower CI. Zero share alone is not sufficient
evidence — measure the finite positive tail, per-model modes, mean stability under
denominator perturbation, and the counterfactual median output.

**10 · The MIROC acceptance test was wrong.** Showing that `dom + ind + irr ≈ ptotww` in
CWatM establishes CWatM's accounting identity and nothing about MIROC's sector coverage,
definitions, source-water treatment or omissions. Inclusion must be justified from MIROC's
own metadata and values; the cross-model comparison is a plausibility check, not the test.

### Redundancy

Monthly maximum is mechanically ≥ the annual ratio whenever denominators are positive, so
those two are not independent information until the gap is measured. Frozen-EFR B2 is an
amplified transform of A2 and the annual pair may be nearly rank-identical. Before retaining
six layers, measure spatial Spearman correlation, threshold reclassification, trend-sign
agreement, reference-site ranking, and the fraction of maxima caused by near-zero denominators.

---

## 5. River routing — ISIMIP publishes it, no HydroBASINS needed

`ISIMIP3b/InputData/geo_conditions/river_routing/` holds three files, ~1.7 MB total:

| file | contents |
|---|---|
| `ddm30_basins_cru_neva.nc` | basin ID per cell — **11,558 basins** |
| `ddm30_flowdir_cru_neva.nc` | D8 flow direction, codes 1–8 plus 0 and −1 |
| `ddm30_slopes_cru_neva.nc` | channel slope |

DDM30 is the 30-arcminute drainage direction map, **natively on our 0.5° grid** — no
regridding. Longitude is the full 720 cells (−179.75 … 179.75); latitude ascends
south-to-north across 280 rows (−55.75 … 83.75), trimming only Antarctica and the polar
ocean. Land mask **67,424 cells**, against the 67,420 that H08, LPJmL and WaterGAP2 use in
the drought layers.

### 5.1 Direction convention — determined, not assumed

A first guess (1=N, clockwise) left 39,773 of 67,424 cells trapped in cycles. A systematic
search over all 16 candidate mappings (8 rotations × 2 handednesses × 2 row orders) found
**exactly one** with zero cycles; the runner-up leaves 2,893.

```
1 = E    2 = SE    3 = S    4 = SW    5 = W    6 = NW    7 = N    8 = NE
```

That is the ESRI ordering renumbered 1–8 instead of powers of two (ESRI: 1, 2, 4, 8, 16,
32, 64, 128 → log₂ + 1). Codes 0 and −1 are sinks and outlets.

Code counts across the 67,424 land cells: `-1`: 782, `0`: 10,776, `1`: 13,920, `2`: 3,466,
`3`: 8,140, `4`: 2,811, `5`: 12,505, `6`: 3,018, `7`: 8,747, `8`: 3,259.

### 5.2 Validation

- **Acyclic and complete** — a topological pass reached all 67,424 land cells.
- **Largest accumulation at 0.25°N, 50.25°W** — the mouth of the Amazon, 1,946 upstream
  cells. At ~3,090 km² per cell at the equator that is **6.0 million km²** against the
  Amazon's actual ~6.1 M km².
- **11,474 of 11,558 basins (99.3%)** have every cell draining to a single terminus. There
  are 14,095 terminal cells for 11,558 basins, so some basins reach the sea through more
  than one mouth — the likely explanation for the 84 exceptions, to be characterised before
  use rather than assumed.

### 5.3 What it changes

Blocker 1 becomes a bug fix rather than a choice. The two denominators survive; the
**numerator** on the river-flow versions becomes upstream-accumulated.

| | before | after |
|---|---|---|
| **L · local** | local demand ÷ local runoff | unchanged — genuinely local self-sufficiency |
| **R · routed** | local demand ÷ routed flow ✗ | **upstream-accumulated demand** ÷ routed flow ✓ |
| ~~B2~~ | — | dropped — environmental flows out of scope |

**The settled metric set is four:** L and R, each as an annual value and as the
highest-stress month.

```
L(c,y,k,m) = W(c) / Q_runoff(c)                      local self-sufficiency
R(c,y,k,m) = W_hat(c) / Q_river(c)                   water stress
  where W_hat(c) = sum of W over every cell draining through c, c included
```

The bands become defensible rather than indicative: the 0.4/0.8 criticality thresholds
originate with Alcamo and Vörösmarty, whose work is grid-scale with upstream-accumulated
demand, so the routed construction is closer to their origin than basin aggregation is.

### 5.4 Upstream accumulation, not basin aggregation

Demand at a cell is summed over every cell draining through it. A headwater sees only its
own small demand against its own small flow; a river mouth sees the whole catchment against
the whole catchment. Numerator and denominator cover the same ground at every point.

The basin-uniform alternative is **rejected**. It assigns every cell in a basin the basin's
total demand, which would charge a Rocky Mountain headwater with Las Vegas's withdrawals —
a claim on water that has already passed the headwater by. Downstream abstraction cannot
reduce upstream supply and the metric must not imply that it does.

The constraint is doubly satisfied: upstream accumulation is one-directional by
construction, and the `1850soc` denominator carries no withdrawal effect anywhere.
`basinnumber` is retained as a diagnostic grouping key only, never as the support of a
published value.

### 5.5 A risk this creates: each model has its own routing

DDM30 was developed for WaterGAP. H08 and LPJmL5 route on their own networks and their
published `dis` reflects those. Accumulating demand on DDM30 while taking discharge from a
model that routes differently would put numerator and denominator back on different
catchments — the very defect the routing was brought in to fix.

Measurable before anything is built: accumulate each model's own `qtot` along DDM30 and
compare cell-wise to that model's published `dis`. Close agreement means the networks concur
and the pairing is sound; divergence localises exactly where that model's drainage differs,
and a model that disagrees badly cannot use the routed metric. This alone justifies
downloading `qtot` for every model, independent of the local metric.

Computationally cheap: the topological order is computed once, then all 86 years × 12
months accumulate in a single pass with the time series as a vector payload.

### 5.4 Two properties to keep visible

- **The network is static.** One flow-direction grid, no time or scenario dimension. No
  channel migration, no new reservoir rerouting flow. Standard practice, and an assumption
  that belongs in the file attributes.
- **DDM30 spans 280 of our 360 latitude rows.** Its land mask is 67,424 against the water
  models' 67,420 — but "within four" is not "identical", and any water-model land cell with
  no routing entry would silently receive no upstream demand. The intersection must be
  checked explicitly during ingest.

---

## 6. Inter-model comparison, and the decision to narrow to one model

**Measured 2026-08-17** on `gfdl-esm4`, `ssp370`, `1850soc`, 2015–2100 time means, all
three models that publish both an all-sector total withdrawal and a naturalised supply
(H08, LPJmL5-7-10-fire, WaterGAP2-2e). 1.43 GB, six files.

### 6.1 Routing agreement with DDM30

Accumulated `qtot` along DDM30 versus each model's own published `dis`. Because
evaporation can only ever make accumulated runoff *exceed* discharge, the fraction of
cells where accumulated is **less** than published is a one-sided network-mismatch
signal, uncontaminated by losses.

| model | acc/dis < 0.9 | median acc/dis by catchment size (1–2 cells → >1000) |
|---|---|---|
| H08 | 0.13% | 1.003 · 1.002 · 1.002 · 1.002 · 1.002 |
| LPJmL5-7-10-fire | 0.10% | 1.005 · 1.005 · 1.007 · 1.014 · 1.023 |
| WaterGAP2-2e | 3.45% | 1.004 · 1.003 · 1.009 · 1.055 · 0.988 |

All three networks are compatible with DDM30. LPJmL5's apparent disagreement is **not** a
network problem — its deficit tail is the smallest of the three and its ratio climbs
monotonically with catchment size, the signature of water lost in transit. WaterGAP2-2e's
deficit tail is 26× H08's and roughly flat across scales, which is a genuine partial
topology difference and is recorded as a known limitation.

**This does not require accumulated runoff in the metric.** The denominator is the model's
own published `dis`, carrying its own channel physics; DDM30 is used only to accumulate
demand, where no such physics applies.

### 6.2 The apparent ranking inverted on inspection

WaterGAP2-2e scored worst on naive agreement **because it is the only one of the three that
simulates lake, wetland and floodplain evaporation** — and it is correspondingly the only
one whose loss-dominated rivers are near reality:

| river | WaterGAP2-2e | H08 | LPJmL5 | observed natural |
|---|---|---|---|---|
| Nile | **3,053** | 8,215 | 16,303 | ~2,700–2,800 m³/s |
| Murray | **232** | 795 | 2,136 | ~380 m³/s |

H08 and LPJmL5 essentially route runoff to the ocean, which is why they agree beautifully
with a naive accumulation and overshoot the Nile by 3× and 6×. **The first agreement score
was measuring model physics, not model plumbing**, and using it to exclude a model would
have removed the most realistic one.

### 6.3 The spread is not a scale offset — the ranking reverses

| ordering (driest → wettest), published `dis` | % cells | % weighted by flow |
|---|---|---|
| H08 < WGAP < LPJmL | 31.6 | 30.3 |
| H08 < LPJmL < WGAP | 23.7 | 28.1 |
| WGAP < H08 < LPJmL | 17.5 | 22.8 |
| WGAP < LPJmL < H08 | 13.4 | 5.2 |
| LPJmL < WGAP < H08 | 7.5 | 6.9 |
| LPJmL < H08 < WGAP | 6.3 | 6.6 |

**No ordering holds in even a third of cells; the ranking reverses in 68.4%.** Pairwise sign
flips: H08/LPJmL 27.2%, H08/WGAP 38.3%, LPJmL/WGAP 37.5%. Pairwise ratios span 0.11–2.22
between the 10th and 90th percentiles — a twentyfold range across space.

There is a systematic tendency (H08 driest in 55% of cells, LPJmL wettest in 49%) but it is
a majority habit, not a rule.

**The disagreement is in the land surface, not the routing.** Reversal rate on *accumulated
runoff*, before any channel physics: **62.9%**. On discharge, after routing: 68.4%. Routing
adds 5.5 points; the other 63 were already present in how much water each model generates
and where.

**And it is regime-structured, not random:**

| basin type | driest model | max/min spread |
|---|---|---|
| arid / loss-dominated | WaterGAP2-2e, by a wide margin | Colorado **29.7×**, Orange 14.9×, Murray 9.2×, Nile 5.3× |
| Arctic | LPJmL5 | Lena 1.3×, Mackenzie 1.9×, Yenisei 1.2× |
| wet tropical / temperate | H08 | Amazon 1.1×, Ganges 1.1×, Congo 1.3× |

**The spread is concentrated precisely in the basins that will score as stressed.** Wet
basins agree within 10–30%; arid basins disagree by 5–30×. Since the 40% and 80% bands are
*absolute* thresholds, a denominator uncertain by 30× makes the band assignment in arid
basins closer to a model choice than a measurement. The relative *ranking* of basins
survives this spread; the absolute class assignment does not.

### 6.4 DECISION — WaterGAP2-2e only, for this pass. TEMPORARY.

**User decision, 2026-08-17.** The initial layer uses **WaterGAP2-2e alone**: 5 GCMs ×
3 scenarios = **5 members per scenario, `n_models` = 1**.

Rationale: of the three candidates it is the only one that simulates open-water evaporation,
which is the dominant term in exactly the arid basins this layer exists to identify, and it
is the only one whose naturalised Nile and Murray are near observed. It also publishes its
own `cellarea`, so the ~0.9% grid-borrowing bias (§3.1) does not arise.

**This is explicitly provisional and is to be revisited.** What it costs:

- **The confidence interval carries no structural uncertainty whatsoever** — only GCM
  spread. `generate_delivery_caveats.py` already emits `SINGLE-MODEL-{layer}` automatically
  when `impact_models` holds one entry, so the caveat fires without new machinery. But its
  standing text says the interval is narrower than the truth "by an unknown amount — there
  is no second model to measure the gap against", and **that clause is now false for this
  layer**: §6.3 measured the gap. The caveat text must carry the measured numbers, and
  whether it is promoted from `should` to `must` is a separate decision, because changing
  the severity globally would affect every other single-model layer.
- Choosing the model that best reproduces observed loss-dominated flow is a defensible
  selection, but it is a **selection on the outcome**, and the layer should say so rather
  than present WaterGAP2-2e as an arbitrary pick.
- The 3.45% deficit-tail cells from §6.1 are now the only routing check we have, with no
  second model to cross-read against.

**Revisit triggers:** a second model gaining an `1850soc` naturalised run with total
withdrawal; a decision to accept the mixed-soc denominator and recover the 3–4 model
ensemble; or a customer requirement for structural uncertainty in the interval.

---

## 7. Downloads

Three-model set, 135 files, ~31 GB:

| variable | soc | files |
|---|---|---|
| `ptotww` | 2015soc | 45 (3 models × 5 GCM × 3 SSP) |
| `qtot` | 1850soc | 45 |
| `dis` | 1850soc | 45 |

Plus ~17 GB if MIROC is included, plus `cellarea`/`contfrac` (three 32 KB files), plus the
1.7 MB of DDM30 routing. The 2020s baseline sits inside the 2015–2100 files, so no
historical download.

---

## 8. Open items before any layer is written

1. `contfrac` in or out — settled by the `pdomww` vs InputData `domww` [m³/yr] comparison,
   verified per model and per variable rather than generalised from one.
2. Confirm the three published `cellarea` grids are identical so one can serve H08 and MIROC.
3. `qtot ≈ qs + qsb`, and where `qg`/`qr` sit — decides whether recharge may enter
   availability at all.
4. Whether `ptotww` actually is flat in time (item 4 above) — measured, not inferred.
5. `ptotww` vs Σ sectors per model, and which sectors each model's "tot" spans; MIROC
   justified from its own metadata.
6. River mask threshold for the routed metrics, proposed from the measured baseline flow
   distribution.
7. Precedence table for every zero, `0/0` and `x/0` case — now the **sole** defence against
   an infinity reaching the pooled statistics.
8. Exact-zero share of `W`, plus the tail diagnostics in item 9 above, before choosing a
   decadal branch.
9. Argmax-month stability diagnostics for the worst-month reduction.
10. Redundancy measurements across the six metrics.
11. Which slope to recommend per metric — measured after the run, never inferred.
12. Reference sites, before and after processing: Indus, Nile, Colorado, Murray–Darling,
    North China Plain and Central Valley should read high; Amazon, Congo and boreal Siberia
    low; a Sahara cell with no demand must read **zero**, not high.
13. Characterise the 84 basins with multiple termini.
15. **Routing agreement per model** — `qtot` accumulated along DDM30 versus that model's own
    published `dis`. A model whose network disagrees cannot use the routed metric.
14. Verify the DDM30 land mask against each water model's land mask.

---

## 9. Taxonomy consequences

- Only the demand-based routed metric may populate `families.water-stress.covered_by`, and
  only **in the same change that routes it in `config/asset_catalog.yaml`**, per that
  file's own rule.
- A supply-only standardized index would not qualify; filing it under `water-stress` would
  assert that the demand-versus-supply question had been assessed when it had not.
- The existing `water-stress` blocker text — "SEPARATE PRODUCT, different output contract …
  do not silently substitute" — describes the Water Risk Index. These layers are the
  alternative route it anticipated, on the TCFD contract, and the blocker needs rewriting
  when one ships.

---

## 10. Next step

Rewrite [water-stress-formulations.md](water-stress-formulations.md) as the settled
specification with the routed construction and all of the review fixes folded in, then send
that revised version back through Codex before any download starts.

---

## 11. Pre-flight measurement 1 — the area convention

> **CORRECTED 2026-08-18, same day.** The first pass concluded `cellarea` alone and was
> **wrong**. Two independent tests reverse it. The error is retained here because the way
> it happened — two mistakes cancelling — is the reusable lesson.

### 11.1 The convention (corrected)

```
V [m3 s-1] = flux [kg m-2 s-1] * cellarea [km2] * (contfrac / 100) * 1000
```

The flux is per **continental** area. `contfrac` **does** multiply in, after division by
100. This is what the file's own `unit_conversion_info` attribute says, and the first pass
dismissed it on bad evidence.

### 11.2 The two tests that settle it

> **Corrected 2026-08-18 after review round 3.** The first write-up of this section
> labelled a column "median ratio" while printing the *ratio of stratum totals*, and its
> global totals were not calendar-weighted. Both are fixed below. The conclusion is
> unchanged; the reasoning that reached it was partly misstated.

**Test (a) — global total withdrawal**, calendar-weighted with `365_day` month lengths,
across **all 15 members** (5 GCMs × 3 scenarios), calendar year 2015:

| candidate | range across members | mean | vs FAO AQUASTAT ~3,900–4,000 km³/yr |
|---|---|---|---|
| A `cellarea` only | 4,995.0 – 5,237.6 | 5,102.1 | **~27% high, every member** |
| B `cellarea × contfrac/100` | **3,944.6 – 4,163.7** | **4,045.7** | **brackets the reference** |

Every one of the 15 members lands inside or beside the reference range under B, and none
does under A. This is **strong corroboration, not precision validation** — AQUASTAT
aggregates country reporting across differing periods and is not a single observed 2015
total.

**Test (b) — residual stratified by `contfrac`.** If `contfrac` were the missing factor,
the residual must track 1/contfrac and applying B must flatten it. Both hold:

| contfrac band | n cells | A: median cell ratio | A: ratio of totals | **B: ratio of totals** |
|---|---|---|---|---|
| 0.000–0.500 | 3,047 | 3.03 | 4.08 | **0.84** |
| 0.500–0.800 | 1,744 | 1.05 | 1.36 | **0.85** |
| 0.800–0.950 | 1,213 | 0.71 | 0.93 | **0.82** |
| 0.950–0.999 | 337 | 0.63 | 0.81 | **0.78** |
| 0.999–1.010 | 45,415 | 0.48 | 0.71 | **0.71** |

Under A the residual spans 4.08 → 0.71, a factor of 5.7 tracking land fraction. Under B it
spans 0.84 → 0.71, a factor of 1.2. **The flattening is the evidence**, not the raw
magnitude in any one band.

**What test (b) does not establish on its own.** The reference product's spatial allocation
could itself correlate with coastal settlement and therefore with `contfrac`, so (b) is
persuasive rather than independently decisive. It is (a) and (b) together — an absolute
total that only B satisfies, and a spatial residual that only B flattens — that settle it.

### 11.3 How the first pass got it wrong — two errors cancelling

The first test compared WaterGAP2-2e's `pdomww` against
`InputData/.../water_abstraction/2015soc/domww` and found candidate A agreeing to 0.6%
globally. Both halves of that were flawed:

1. **The reference is not the same quantity.** That input file's own `comment` attribute
   reads *"prepared from ISIMIP2a data (H08, PCR-GLOBWB, WaterGAP) for ISIMIP2b"* — it is a
   **three-model average from a previous simulation round**, not the ISIMIP3b WaterGAP2-2e
   forcing. Its `history` even records the `cdo setmisstoc,0` that zero-filled it.
2. **WaterGAP's domestic withdrawal is genuinely ~26% below that 2b-era average**, and
   candidate A overestimates by ~27% for want of `contfrac`. **The two cancelled**, producing
   a spurious 0.6% agreement.

The claim that "a units error would appear as a constant multiplicative offset" was also
**false**: `contfrac` varies spatially and correlates with coastal settlement, so a missing
`contfrac` is a *spatially structured* error, which is precisely what §11.2 detects.

**Lesson, and it generalises:** a global-total agreement between two products of different
provenance is weak evidence. It cannot distinguish "correct" from "two errors cancelling".
Prefer a test whose signature is *structured* — here, the residual tracking 1/contfrac —
over one that yields a single scalar.

### 11.4 Findings that survive the correction

1. **`contfrac` is in PERCENT (0–100) while its `units` attribute says `1`.** Measured max
   exactly 100.00000, median 100.00000, min 1.8e-5. Now load-bearing: the conversion divides
   by 100, and taking the attribute at face value is a **100× error**.
2. **`np.asarray()` silently discards a masked array's mask.** *(Diagnosis corrected in
   round 3 — the original text blamed `netCDF4` and that was wrong.)* `set_auto_mask(True)`
   **does** mask correctly: it returns a MaskedArray with 191,781 of 259,200 cells masked.
   The fills reappear the moment that array goes through `np.asarray()`, which drops the
   mask and restores the raw 1e20 — which is *finite* and passes `isfinite`. An unguarded
   sum returned 3.6e41. The defence stands (**mask `|x| >= 1e19` explicitly**) but the
   reason is the conversion, not the library.
3. **The InputData `domww` zero-fills the globe** — 138,264 of 201,600 cells are exact zeros,
   all finite. Same trap as `cropfailure-3b`.
4. **`contfrac < 0.99` on 10,556 of 67,419 land cells — 15.7%.** *(Corrected in round 3;
   the earlier 4.1% wrongly used all 259,200 grid cells, ocean included, as the
   denominator.)* The correction therefore touches **one land cell in six**, not one in
   twenty-five — and it concentrates on coastlines, where population and demand concentrate.

### 11.5 CALENDAR — corrected

The plan asserted `proleptic_gregorian` with real leap years. **All 45 downloaded
WaterGAP2-2e monthly files declare `calendar = 365_day`.** Annual weighting must use no-leap
month lengths and a 365.0-day year. Monthly *ratios* are unaffected because the duration
cancels; annual sums are not.

**This is heterogeneous across models** — H08's sidecar declares `proleptic_gregorian` for
the same variable and period. Read the calendar per file; never assume it.

The files carry **no `time` bounds and no `cell_methods`**, so "monthly values are monthly
means" is an assumption, not a declaration, and the processor must record it as such.

## 12. Plan for completing the work

Data is in hand: 47 files, 11 GB, `data/raw/water_stress/`, WaterGAP2-2e, 5 GCMs ×
3 scenarios × {`ptotww`@2015soc, `qtot`@1850soc, `dis`@1850soc} + `cellarea` + `contfrac`,
plus DDM30 routing (1.7 MB) and the two files used for the area test.

### 12.1 Stage A — routed demand

Build the DDM30 topological order once (verified acyclic, 67,424 cells, mapping
`1=E 2=SE 3=S 4=SW 5=W 6=NW 7=N 8=NE` for ascending latitude). Note DDM30 latitude
**ascends** (−55.75 → 83.75) while the model output **descends** (89.75 → −89.75) and
covers 360 rows against DDM30's 280 — align by exact coordinate match, never by
assuming orientation.

Per member (15) and month (1032):

```
W(c,t)     = ptotww(c,t) * cellarea(c) * (contfrac(c)/100) * 1000 * dt(t)   [m3]
W_hat(c,t) = sum of W over every cell draining through c, c included
```
`dt` uses **365_day** month lengths (§11.5). All inputs fill-masked at `|x| >= 1e19`
BEFORE any arithmetic (§11.4).

Accumulation is one topological pass with the 1032-month series as a vector payload.

### 12.2 Stage B — the four metrics

```
Q_run(c,t) = qtot(c,t) * cellarea(c) * (contfrac(c)/100) * 1000 * dt(t)   [m3]
Q_riv(c,t) = dis(c,t)  * dt(t)                             [m3]

L_ann(c,y) = SUM_k W(c,y,k)     / SUM_k Q_run(c,y,k)
L_max(c,y) = MAX_k [ W(c,y,k)     / Q_run(c,y,k) ]
R_ann(c,y) = SUM_k W_hat(c,y,k) / SUM_k Q_riv(c,y,k)
R_max(c,y) = MAX_k [ W_hat(c,y,k) / Q_riv(c,y,k) ]
```

No cap at 1. Values above 1 are a real state and must survive.

### 12.3 Precedence for degenerate cells — REBUILT twice (rounds 2 and 3)

Two phases. Conflating them was a defect: a monthly rule applied to a year, or an annual
rule applied to a month, changes which cells survive.

**Phase 1 — per month, per metric.**

| # | condition | result | note |
|---|---|---|---|
| 1 | **domain**: routed metrics require a DDM30 entry; local metrics do not | NaN | applying this to local metrics discards valid model land |
| 2 | any input `abs >= 1e19`, or NaN | NaN | fill normalisation, before any arithmetic |
| 3 | **any input negative** | NaN, **counted** | physically impossible; must surface, and must be tested before the zero and mask rules or it hides inside them |
| 4 | **routed metrics: cell fails the FIXED river mask** | NaN | **moved ahead of the zero rules in round 3.** Previously an off-river cell with zero demand returned `0` instead of being masked, and an off-river cell with demand and no flow could trigger the rule-6 hard stop before ever being excluded. The mask is derived once from the shared 2020s baseline climatology and held — thresholding each year's own denominator would create state-dependent missingness exactly when rivers dry |
| 5 | numerator `== 0` — `W` for local, **`W_hat`** for routed | **0** | no demand is not stress. Routed metrics test the *accumulated* numerator |
| 6 | numerator `> 0` and denominator `== 0` | NaN, **counted and mapped** | the most-stressed state, deleted by necessity. If material anywhere on the eligible domain, **stop and design an explicit representation** — a caveat is not sufficient |
| 7 | otherwise | monthly ratio, float64 | |

**Phase 2 — per year.**

| # | condition | result |
|---|---|---|
| 8 | fewer than 12 valid months in the year | NaN — annual sums and maxima must never be built from partial years; `nansum`/`nanmax` would publish them silently |
| 9 | otherwise | `SUM_k` for the annual metrics, `MAX_k` for the worst-month metrics |
| 10 | **postcondition** | assert no `inf` anywhere in either phase; equality-with-zero is not an overflow defence against a near-zero denominator |

**Accumulation needs its own missing-input policy.** Ignoring a missing upstream land cell
understates every downstream numerator; propagating NaN without first fixing the valid land
mask blanks an entire river. Establish the land mask first, then treat within-mask gaps as
a declared, counted decision.

**On what NaN actually does:** it does not numerically poison `nanmedian`/`nanmean` — it
silently changes the sample and the per-cell member coverage, which is why rules 3 and 6
must be counted and mapped rather than merely masked. An *infinity* does break the mean
branch; whether it moves a median depends on how many there are.

### 12.4 Stage C — the output contract

Per `OUTPUT-SPEC.md`, with the baseline ordering from review finding 6:

1. Shared 2020s baseline: average **each member** across scenarios first, then pool.
2. Decadal statistic: pool (year × member) within the decade. **Branch to be measured,
   not assumed** — see 12.5. `n_models` = **1**. `n_members` must be emitted **per cell and per decade**, not as a
   constant — rules 5, 6 and 7 reduce coverage unevenly. There are five member *identities*
   repeated across three scenarios, not fifteen members.
3. Percentile against the shared 2020s baseline; two-tier if the baseline exceeds 2%
   exact zeros, which it almost certainly will.
4. `ols_slope` and `sen_slope` on the expanding window; recommended slope measured after
   the run.
5. `percentile_direction: higher_is_worse`. `relative_baseline: false` — this is an
   absolute ratio.

### 12.5 Measurements that gate the run

| # | measurement | decides |
|---|---|---|
| 1 | exact-zero share of `W`, and of each metric | decadal branch; two-tier percentile |
| 2 | shape of the pooled per-cell distribution (tail, modes) | whether `pooled_mean_zero_inflated` is safe on an unbounded ratio, or the median branch holds |
| 3 | baseline `dis` distribution | the `river_mask` threshold, chosen from the histogram and declared |
| 4 | count of rule-4 and rule-6 cells | whether either is material enough to caveat |
| 5 | argmax-month distribution and switch rate for `L_max`/`R_max` | whether the worst-month metrics are stable enough to publish |
| 6 | `ols` vs `sen` agreement on active cells | recommended slope per metric |
| 7 | pairwise rank correlation across the four metrics | whether four layers carry four pieces of information |
| 8 | temporal variance and slope of `ptotww` per member | **whether demand is actually flat** — review finding 4, still open |
| 9 | reference sites | Indus, Nile, Colorado, Murray–Darling, North China Plain, Central Valley high; Amazon, Congo, boreal Siberia low; a Sahara cell with no demand **exactly 0** |

### 12.6 Stage D — registry and delivery

- `config/layer_registry.yaml`: four entries. **Each metric gets its OWN 2020s percentile
  reference**, shared only across that metric's three scenarios. The flood precedent does
  NOT apply — those were aligned protection variants of one quantity; annual, worst-month,
  local and routed are four different estimands with different distributions and different
  masks, and ranking maxima against an annual reference would mechanically inflate them. `impact_models: watergap2-2e` — one entry, which is what makes
  `generate_delivery_caveats.py` emit `SINGLE-MODEL-{layer}` automatically.
- The `SINGLE-MODEL` caveat text says the interval is narrower "by an unknown amount —
  there is no second model to measure the gap against". **That is now false**; §6.3
  measured it. The text must carry the numbers. Promotion from `should` to `must` is a
  separate decision because the severity is global.
- `config/hazard_taxonomy.yaml`: only the routed metric may populate
  `families.water-stress.covered_by`, and only in the same change that routes it in
  `config/asset_catalog.yaml`.
- Rewrite the `water-stress` family blocker, which currently describes the Water Risk
  Index and warns against silent substitution.
- `python scripts/test_shared_baseline.py {processed_dir}` must pass.

### 12.7 Known limitations to carry into the file attributes

1. One impact model; the interval carries **no structural uncertainty**, and §6.3 measured
   the missing component as large (68.4% ranking reversal; 29.7× on the Colorado).
2. The model was **selected on the outcome** — best reproduction of observed loss-dominated
   flow. Defensible, not neutral, and the layer should say so.
3. Demand is `2015soc`; socioeconomic *conditions* are fixed. Whether *demand* is fixed is
   measurement 8 and is **not yet established**.
4. Supply is `1850soc`: preindustrial conditions, so the metric **cannot show the benefit
   of reservoir storage**, and preindustrial land use also changes runoff.
5. 3.45% of cells show a DDM30/WaterGAP routing deficit, now with no second model to
   cross-read against.
6. The 40%/80% bands are absolute thresholds; §6.3 showed arid-basin denominators varying
   5–30× across models, so **band assignment there is fragile** even though basin ranking
   is robust.
7. Environmental flows are out of scope; this is Aqueduct-style baseline water stress, not
   SDG 6.4.2.
