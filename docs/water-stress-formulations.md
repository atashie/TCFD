# Water stress — candidate formulations

**Status: design. Nothing is built, no data has been enumerated, and no formulation is
chosen.** Written 2026-08-15 from the literature review of the same date. Every threshold
and rule below carries its source; every property of *our* data marked UNVERIFIED has not
been measured yet and must be before any of this is coded (GUARDRAILS §9).

This document exists because three separate layers are on the table and only one of them
may be called "water stress" in a filing. Keeping them in one file is what stops the other
two from drifting into that name.

---

## 0. What the data constraint removed, and what it left

Declared available: **total water demand** (withdrawal), **total consumptive use**,
**natural groundwater recharge**, **soil moisture**, **natural volumetric streamflow**.
Declared *not* available: total soil matrix storage potential — so no field capacity, no
wilting point, no plant-available water.

That single absence closes an entire family. SWDI, relative extractable water and the
FAO-56 stress coefficient all divide by `θ_FC − θ_WP`:

| Removed | Formula | Why it is out |
|---|---|---|
| SWDI (Martínez-Fernández 2015) | `10·(θ − θ_FC)/(θ_FC − θ_WP)` | both constants unavailable |
| REW (Granier 1999) | `(θ − θ_WP)/(θ_FC − θ_WP)` | same denominator |
| FAO-56 `Ks` (Allen 1998) | `TAW = 1000(θ_FC − θ_WP)Z_r`; `Ks = (TAW − D_r)/((1−p)TAW)` | TAW is that denominator times rooting depth |

There is no substitute that preserves their meaning. An *empirical* range normalisation
(§1.3) is available and is a different quantity — the sample extremes of a simulation are
not soil hydraulic constants, and treating them as such is the kind of inference GUARDRAILS
§9 forbids.

**The consequence that shapes everything below:** soil moisture can now only be expressed
**relative to its own history**, while the blue-water fields (demand, consumption, flow)
carry their own units and can be expressed **absolutely**. So W1 is necessarily a
departure-from-baseline layer and W3 is necessarily not — and those two things read in
opposite directions, which is the distinction `relative_baseline` exists to enforce.

---

## 1. Layer W1 — soil moisture deficit (green water)

**This layer may not be named "water stress".** ESRS E3 defines that term as a
withdrawal-to-availability ratio (§3). W1 answers "is the soil drier than this place is
used to", which is agricultural/ecosystem water stress. Candidate names: *soil moisture
deficit*, *agricultural water stress*. It is a sibling of `drought-3b`, not of W3.

### 1.0 Source field — an open question before anything else

Two candidates and they are not interchangeable: `soilmoist` (total column) and
`rootmoist` (root zone). The Water Risk Index took `rootmoist` and found **one usable
model** — WEB-DHM-SG, after excluding miroc-integ-land for glacier mass in the field
(`scripts/config_water_variables.py`). A one-model ensemble carries no structural
uncertainty and would make the CI meaningless. `soilmoist`'s ensemble depth is UNVERIFIED.

Whichever is used, **column depth differs by model**, so raw `kg m⁻²` is not poolable
across an ensemble — the structural problem `permafrost-3b` solved by normalising each
model against its own column. All three formulations below normalise per model by
construction, which is the reason to prefer them over any raw-value treatment.

### 1.1 Standardized Soil Moisture Index, nonparametric (SSI)

Per member *m*, cell *c*, year *y*, reduce the monthly series to one annual value
`A(m,c,y)` (choice in §1.4), then rank it against the cell's own reference sample pooled
over members and years:

```
R(c)   = { A(m,c,y) : y ∈ reference period, m ∈ members }      N = |R(c)|
p      = (i − 0.44) / (N + 0.12)        Gringorten plotting position, i = rank of A in R(c)
SSI    = Φ⁻¹(p)
D      = −SSI                            published value; positive = drier than normal
```

Publishing the **deficit** `D` rather than `SSI` keeps the raw value and the contract
percentile pointing the same way (higher = worse) and avoids `percentile_direction:
higher_is_better`, which every downstream reader has to remember to un-invert.

**The reference sample size is a hard ceiling on the index, and it censors.** Gringorten's
maximum plotting position is `(N−0.44)/(N+0.12)`, so:

| Reference period | N (10 members) | max \|SSI\| |
|---|---|---|
| 2020s decade only | 100 | **2.54** |
| 2020s decade, 5 members | 50 | **2.28** |
| historical 1971–2014 | 440 | **3.02** |

By the 2090s a drying cell parks at that bound, the CI collapses, and both slopes go to
zero and *agree* — the exact `heatwave-3b` pathology, where agreement is ambiguous between
"no trend" and "no headroom left" (OUTPUT-SPEC, censored fields). Three mitigations, all of
which must be decided rather than defaulted:

1. **Use the historical runs as the reference period**, not the 2020s decade. Raises the
   bound to ~3.0σ and costs a download of `historical` members.
2. **Fit a parametric distribution** instead (gamma/log-logistic per cell), which
   extrapolates past the sample. Buys range, spends it on distribution-choice
   artifacts — Tijdeman et al. (2020) measured those at large sample scale.
3. **Publish the unstandardized anomaly** — `(A − mean_baseline)/SD_baseline` per model,
   then pool. No ceiling at all, no plotting position, still per-model normalised. Loses
   the probability interpretation the standardized indices are chosen for.

Whichever is taken, the share of cells at the bound must be measured per panel and, if
material, declared in a `saturation_caveat`.

### 1.2 Soil moisture percentile (SMP)

`SMP = 100 · p` from the same empirical CDF — the US Drought Monitor / Sheffield & Wood
construction, no distribution assumption at all.

**Do not adopt this without resolving the name collision.** The output contract already
publishes a variable called `percentile`, which is a rank against the *global 2020s
baseline distribution*. Publishing a per-cell historical percentile as the `median` value
would produce a file whose `percentile` variable is the global rank of a local percentile.
Not wrong, but two different percentiles in one file is a mis-reading waiting to happen,
and it would reach the customer CSV as two columns with the same word in them.

### 1.3 Empirical saturation index — the only quasi-absolute route left

```
S = (θ − θ_min) / (θ_max − θ_min)        per model, per cell, over the reference period
```

This is the closest available analogue to plant-available water and it is **a proxy, not a
substitute**: `θ_min`/`θ_max` are simulation sample extremes, so they move with record
length, member count and model, and they are not soil hydraulic properties. The Water Risk
Index already does a version of this — `rootmoist` normalised by a single model's observed
maximum of 1187.29 kg m⁻² — and that constant is model-specific, which is precisely the
weakness to declare.

If taken, it needs: per-model `θ_min`/`θ_max` (never pooled), the record length stated in
the file's attributes, and an `interpretation_caveat` saying it is not field capacity.

### 1.4 The annual reduction (all three formulations need one)

| Reduction | Reads as | Property |
|---|---|---|
| annual mean | average wetness | smooth, insensitive to timing, weakest signal |
| annual **minimum** monthly | the driest month | matches the seasonal framing the scarcity literature insists on (Liu 2017; Brauman 2016) |
| months below the cell's baseline 20th percentile | duration of deficit | bounded 0–12, **censors** at 12, zero-inflated at 0 |

The count is attractive and carries two contract hazards at once (a ceiling and
zero-inflation). Measure all three on one member before choosing.

---

## 2. Layer W2 — streamflow stress (blue water supply)

Input: natural volumetric streamflow, `m³ s⁻¹`, on the 0.5° grid.

**A cell's discharge is the routed accumulation of everything upstream.** Read as "the
river at this location" it is fine; read as "this cell's water" it is wrong. Cells with no
significant river carry a tiny discharge dominated by local runoff, so a published layer
should either carry a `river_mask` derived from baseline mean discharge (threshold
declared) or state that low-flow cells are meaningless. This is a masking decision, not a
caveat to write later.

### 2.1 Standardized Streamflow Index (SSI-dis)

The §1.1 construction applied to annual or accumulated-monthly discharge — Vicente-Serrano
et al. (2012), with the nonparametric form as in Tijdeman et al. (2020). Same censoring
arithmetic, same reference-period decision.

Two additions specific to flow:

- **Timescale must be declared and justified.** The HESS 2025 technical note shows the
  SSI's meaning changes with accumulation window and flow regime — snowmelt basins track
  multi-month storage, rainfall basins track recent anomalies — and that the apparent
  SPI→SSMI→SSI propagation order can reverse on timescale choice alone.
- **Intermittent rivers break the fit.** Any cell with exact zeros in the reference sample
  needs an explicit rule (mask, or a two-tier treatment like the zero-inflated percentile).
  Measure the zero share first.

### 2.2 Low-flow exceedance (fixed threshold from the baseline)

```
T(c)          = Q90 of the baseline monthly discharge sample, pooled over members    [m³ s⁻¹]
months_below  = Σ_months 1[ Q(c,month) < T(c) ]                                      [0–12]
deficit_vol   = Σ_months max(0, T(c) − Q(c,month)) · Δt_month                        [m³]
```

Anchored to the regulatory family — Q90/Q95 from the flow-duration curve, MAM7 (UK
abstraction licensing), 7Q10 (US discharge permits) — and absolute once `T` is fixed.

**Prefer the deficit volume over the month count.** The count has a ceiling at 12 and will
saturate in drying basins, reproducing the censoring problem; the volume has none. If the
count is published anyway, treat it as a censored field from the start.

### 2.3 Environmental-flow shortfall (VMF)

Computable from natural flow alone, and it produces the EFR term W3 needs, so building it
once serves both layers.

```
MAF(c)              = baseline mean annual flow
MMF(c,month)        = baseline mean flow of that calendar month
class(c,month)      = low          if MMF < 0.40·MAF        →  f = 0.60
                      intermediate if 0.40 ≤ MMF ≤ 0.80·MAF →  f = 0.45
                      high         if MMF > 0.80·MAF        →  f = 0.30
EFR(c,month)        = f · MMF(c,month)
shortfall_vol       = Σ_months max(0, EFR(c,month) − Q(c,month)) · Δt_month          [m³]
```

Pastor et al. (2014); VMF and Tessmann correlated best with locally derived requirements
(R² = 0.91). Global mean EFR ≈ 37% of mean annual flow, range 25–46% across methods — that
spread is the size of the choice.

**Fix `MMF`/`MAF` on the baseline period and hold them.** If the requirement is recomputed
per decade it drifts with the climate and the index becomes structurally unable to detect
drying — the moving-baseline trap.

---

## 3. Layer W3 — water stress (the framework definition)

This is the only one of the three that answers what ESRS E3 and CDP mean by water stress,
and the only one whose `covered_by` may populate the `water-stress` family in
`config/hazard_taxonomy.yaml`.

### 3.1 Unit reconciliation, first

Demand and recharge arrive as fluxes (`kg m⁻² s⁻¹`), flow as a volume rate (`m³ s⁻¹`):

```
V [m³ s⁻¹] = flux [kg m⁻² s⁻¹] · A_cell [m²] / 1000
```

`A_cell` varies with latitude and must come from the same grid definition the layer
publishes on. UNVERIFIED: the actual units on the demand and consumption files — read them
from the values and the attributes, not from this document.

### 3.2 The denominator, where the double-counting trap is

| Construction | Availability | Property |
|---|---|---|
| **D1 river** | `Q(c)` routed discharge | denominator is upstream-accumulated while the numerator is local — mismatched support; the ratio is not WTA in the literature's sense |
| **D2 local** | locally generated runoff | supports match; a city on a big river with no local runoff reads infinite stress — the classic failure of grid-cell WTA |
| **D3 basin** | `Q` at basin outlet, demand summed over the basin | the Aqueduct / SDG 6.4.2 construction, and the only one the thresholds were calibrated on. Needs a basin or flow-direction asset (HydroBASINS, DDM30) we do not have |

D3 makes every cell in a basin share one number. That is a **support** change, not an
accuracy loss, and it is exactly the situation `resolution_caveat` was built for on
`sealevel-2b`. It is also what Aqueduct publishes, so a customer comparing our number to
Aqueduct's is comparing like with like only under D3.

**Groundwater recharge is not a free addition.** In most global hydrological models the
routed discharge already contains the baseflow that recharge feeds, so `Q + recharge`
double counts. Aqueduct includes rechargeable groundwater because its supply is assembled
from runoff *components*, not from routed discharge. Before recharge enters any
denominator, run the water-balance check per model — does total runoff already include the
groundwater term? — and record the answer. Model conventions differ; this is a measurement,
not a convention to assume (GUARDRAILS §9).

### 3.3 The indices

```
WTA        = W / A                          withdrawal-to-availability  (Alcamo, Vörösmarty)
Depletion  = C / A                          consumption-to-availability (Brauman 2016)
WS_SDG     = W / (A − EFR)                  SDG 6.4.2 / Smakhtin 2004, EFR from §2.3
```

Reference bands, for the report to apply — not to bake into the values:

| Scheme | Bands |
|---|---|
| Aqueduct 4.0 / **ESRS E3** | high **40–80%**, extremely high **>80%** |
| SDG 6.4.2 | <25 none · 25–50 low · 50–75 medium · 75–100 high · **>100 critical** |
| Criticality ratio | >0.4 high, >0.8 extreme |

**Withdrawal and consumption are two variants of one layer, not two layers with one
number.** The flood family is the precedent — three protection variants sharing one
percentile reference (`config/layer_registry.yaml`). Same structure fits here:
`waterstress-wta` and `waterstress-depletion`, both ranked against the same 2020s baseline
so a percentile means the same thing across them.

### 3.4 Rules that must be decided, not defaulted

- **Do not clip the ratio at 1.** Withdrawal exceeding renewable supply is a real state
  (groundwater mining) and is SDG's ">100% critical" class. Clipping manufactures a
  censored field with the `heatwave-3b` pathology.
- **`A − EFR ≤ 0`** happens in arid basins. Mask, or publish a declared cap — either way
  the rule goes in the file's attributes.
- **`W = 0` is not stress.** A desert with no demand scores zero and that is the indicator
  working as defined. Expect a strongly zero-inflated field, expect the third decadal
  branch to be a candidate, and measure the exact-zero share before choosing (OUTPUT-SPEC
  §"third branch").
- **Percentile vs band.** The contract's `percentile` is a global relative rank; the ESRS
  test is an absolute band. Both can travel, but the file must say which is which, or a
  reader will run the ESRS test on the percentile. `sealevel-2b` has the same duality and
  declares it.

### 3.5 Demand is socioeconomic, and this decides what the layer means

If the available runs hold the socioeconomic forcing fixed, `W` and `C` are constant in
time and **every trend in the ratio is supply-side** — the slopes measure discharge, not
stress evolution. That is a defensible layer, and a completely different one from what
"water stress projections" implies to a reader. If demand *does* vary, the layer combines a
socioeconomic pathway with a climate pathway and the report must name both.

UNVERIFIED and pivotal: which `soc` forcing the demand files carry, and whether demand is
time-varying. Nothing else in this section can be finalised before that is read off the
files.

---

## 4. How each plugs into the output contract

| | W1 soil moisture | W2 streamflow | W3 water stress |
|---|---|---|---|
| Annual value | §1.4 reduction | shortfall / deficit volume, or SSI | ratio (dimensionless) |
| Field nature | continuous, ~symmetric | zero-inflated (volumes) | strongly zero-inflated + heavy right tail |
| Expected branch | `pooled_median` | **measure** — `pooled_mean_zero_inflated` likely | **measure** — same |
| `percentile_direction` | higher_is_worse (publishing the deficit) | higher_is_worse | higher_is_worse |
| `relative_baseline` | **true** — must-disclose caveat | true for SSI, false for 2.2/2.3 | **false** |
| Censoring risk | high (§1.1) | high if a month-count is published | only if the ratio is clipped |
| Mask | model soil column, arid low-variance cells | `river_mask` | `A − EFR ≤ 0`, no-demand cells |
| Support | cell | cell = river reach | **basin** under D3 |

Every "expected branch" above is a hypothesis. Both slopes get measured with
`generate_customer_delivery.py --measure-slopes` on active cells after the run, and the
recommended slope follows the measurement, never the table.

---

## 5. Pre-flight measurements — none of this is codeable before them

1. **Units and nature from the values** for all five fields, per model (GUARDRAILS §9).
2. **Soil moisture ensemble depth**: how many models publish `soilmoist` annually, and what
   column depth each uses. `rootmoist` is one usable model in the Water Index today.
3. **Does total runoff already include groundwater recharge?** Per-model water balance.
4. **`soc` forcing on the demand and consumption files; is demand time-varying?**
5. **Ensemble intersection.** The models publishing demand are probably a subset of those
   publishing discharge, and `n_members` is the intersection, not the maximum of the three.
6. **Zero shares**: exact-zero discharge cells, zero-demand cells, and the exact-zero share
   of each candidate annual index.
7. **Reference sites (GUARDRAILS §12)**, before and after processing. High stress must
   appear in the Indus, Nile, Colorado, Murray-Darling, North China Plain, Central Valley;
   low in the Amazon, Congo, boreal Siberia. A Sahara cell with no demand must read zero,
   not high — if it reads high, the denominator is broken.

---

## 6. Taxonomy consequences

- W3 populates `families.water-stress.covered_by` — **in the same change that routes it in
  `config/asset_catalog.yaml`**, per that file's own rule.
- W1 and W2 do **not**. Filing either under `water-stress` would state that the
  demand-versus-supply question had been assessed when it had not — the failure mode
  `hazard_taxonomy.yaml` exists to prevent. W1 belongs beside drought (agricultural
  moisture deficit); W2 is hydrological drought / low flow and has no family yet.
- The existing `water-stress` blocker — "SEPARATE PRODUCT, different output contract … do
  not silently substitute" — describes the Water Risk Index. These three layers are the
  alternative route it anticipated, on the TCFD contract, and the blocker text needs
  rewriting when one ships.

---

## Sources

Martínez-Fernández et al. 2015 (SWDI) · Granier et al. 1999 (REW) · Allen et al. 1998
FAO-56 ch.8 (TAW/RAW/Ks) · Hao & AghaKouchak 2013, Farahmand & AghaKouchak 2015
(nonparametric standardized indices) · Sheffield & Wood 2004 (soil moisture percentile) ·
Vicente-Serrano et al. 2012 (SSI) · Tijdeman et al. 2020 WRR (parametric vs nonparametric) ·
HESS 29:1981 2025 (what SSI reflects) · Pastor et al. 2014 HESS (VMF, EFR comparison) ·
Smakhtin et al. 2004 (WSI with environmental water requirements) · Richter et al. 2012
(presumptive standard) · Alcamo et al. 1997/2003, Vörösmarty et al. 2000 (criticality
ratio) · Brauman et al. 2016 (depletion, seasonal and dry-year) · Hoekstra & Mekonnen 2012,
Mekonnen & Hoekstra 2016 (monthly blue water scarcity) · Liu et al. 2017 Earth's Future
(review) · Kuzma et al. 2023 (Aqueduct 4.0) · FAO/UN-Water SDG 6.4.2 metadata · EFRAG ESRS
E3 · Boulay et al. 2018 (AWARE) · Schewe et al. 2014 PNAS, Gosling & Arnell 2016 (ISIMIP
precedent).
