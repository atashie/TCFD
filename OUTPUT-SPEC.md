# TCFD/CDP Output Specification

The contract every hazard layer must satisfy, so layers are mutually comparable.
Implemented by `scripts/utils/decadal_stats.py`; tested in
`isimip-pipeline/tests/test_decadal_stats.py`.

Supersedes the two divergent families that grew up in this codebase (see
[Appendix: what changed](#appendix-what-changed)).

## File

`{variable}_{scenario}_processed.nc` — **one file per scenario**, dims
`(decade, lat, lon)` = `(9, 360, 720)` at 0.5°. Decades `2010…2090` (layers whose source
starts in 2015 begin at the `2020` baseline instead). All scenarios of a layer live in
**one folder**; never split folders per scenario.

## The sample

Two aggregation samples are used, and the difference is load-bearing:

| Name | Definition | Used by |
|---|---|---|
| **decade pool** | every `(year, member)` observation inside the 10-year window | `median`, `lower_ci`, `upper_ci` |
| **expanding stack** | every `(year, member)` observation from the START of the baseline decade through the END of the target decade | `ols_slope`, `sen_slope` |

A *member* is one impact-model × GCM combination. Members are pooled as independent
observations — see [Known consequences](#known-consequences).

## Variables

| Variable | Continuous field | Boolean {0,1} field | Units | Baseline panel |
|---|---|---|---|---|
| `median` | **median** of the decade pool | **mean** of the decade pool | native | shared across scenarios |
| `lower_ci` | **25th percentile** of the decade pool | **mean − 1 SD** of the decade pool | native | shared |
| `upper_ci` | **75th percentile** of the decade pool | **mean + 1 SD** of the decade pool | native | shared |
| `percentile` | percentile-of-score vs the shared 2020s baseline distribution | same | 1–100 | shared |
| `ols_slope` | least-squares slope over the expanding stack | same | native **per year** | **NaN** |
| `sen_slope` | Theil-Sen slope over the expanding stack | same | native **per year** | **NaN** |
| `n_members` | count of contributing members per cell | same | count | per scenario |
| `n_models` | count of contributing model *families* per cell | same | count | per scenario |

`median` keeps its name for backward compatibility even though it holds a **mean** on
boolean layers. The layer must declare which branch was taken via the
`decadal_statistic` global attribute (`pooled_median`, `pooled_mean_boolean`,
`pooled_mean_zero_inflated`, or `pooled_mean_multimodel`).

### Boolean detection

Taken from the **values**, never the name (GUARDRAILS §9). A field is boolean only if
every finite value is exactly 0 or 1. `led` is binary; `let` is a continuous fraction in
[0,1) despite the sibling naming. A *zero-inflated* continuous field is **not** boolean
and takes the median/IQR branch — unless it is zero-inflated to the degenerate degree
described next. The scan is exhaustive — a sampled scan can step over a rare continuous
tail and silently switch the published statistic.

### The third branch: extreme zero-inflation

The boolean/continuous split is a *proxy* for the real question: **is the decade pool's
distribution degenerate at zero?** The proxy holds for almost every field. It fails when
a continuous field is so zero-inflated that its median is 0 nearly everywhere, and there
the median branch does not lose precision — it erases the layer.

Measured on `let` (tropical-cyclone exposure), 2020s panel after smoothing:

| branch | exposed cells | exact-zero land |
|---|---|---|
| pooled median + IQR | 2,684 | 96.07% |
| pooled mean ± 1 SD | 15,122 | 77.84% |

93% of exposed land reads exactly 0 under the median, and the two-tier percentile then
assigns tier 1 to 96% of land. `let` is 97.84% exact-zero at annual resolution; for
contrast `burntarea` is 29.2% and takes the median branch without difficulty.

A layer in this regime may take **mean ± 1 SD** — the same treatment the spec already
grants boolean fields, for the same reason: `let` and `led` both estimate an *expected
annual exposed fraction*, and a median of near-Bernoulli draws is meaningless. Such a
layer passes `central="mean"` to `pooled_decadal_stat` and **must** declare
`decadal_statistic: pooled_mean_zero_inflated` plus a
`decadal_statistic_rationale` attribute carrying the measured numbers.

This is a **declared** deviation, never a silent one. Do not reach for it to improve
contrast on an ordinary field: measure the median branch's exact-zero share first, and
record it. Each adoption is its own decision on its own measurement, and the count is not a
precedent — **which layers qualify today, and on what numbers, is in
[DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)** (two, as of 2026-08-13).

### The fourth branch: a multimodal ensemble

The median assumes the pooled sample has **one** mode. When the members separate into
clusters, the median does not summarise them — it *selects* the cluster holding more
members, and it moves discontinuously when the balance tips. A layer in that regime may
take **mean ± 1 SD**, declaring `decadal_statistic: pooled_mean_multimodel` plus a
`decadal_statistic_rationale` carrying the measured numbers.

Measured on `permafrost-3b` (added 2026-08-14, the only layer in this regime today). Its
value is thaw depth normalised by each model's own soil column, and the three models sit at
different fractions of their columns in the 2020s:

| | CLASSIC (61.4 m column) | LPJmL (13.0 m) | JULES (3.0 m) |
|---|---|---|---|
| 2020s normalised thaw | 0.035 | 0.046 | **0.951** |

Seven members in the low cluster, five in the high one. Under the median branch the ssp585
spatial median went **0.40 (2080s) → 0.93 (2090s)** — not thaw accelerating, but the median
crossing between clusters as the high cluster gained the majority. The mean moves smoothly
through the same panels and the SD carries the disagreement, which on that layer is the
honest content: the 2-SD width has a median of 0.748 on a [0,1] field.

Two conditions distinguish this from reaching for the mean to improve contrast:

- the separation is **between models**, not within one, so it is structural uncertainty
  rather than a heavy tail — verify by measuring the per-model distributions separately;
- the median is **discontinuous in the decade series**, which is visible as a jump no
  physical process explains. Record the counterfactual at run time.

**A threshold applied to the central value inherits this choice, so do not report one as if
it were a property of the data.** On `permafrost-3b`, "area whose column is fully thawed"
means *>half the members* under the median and *effectively all of them* under the mean —
7.97 vs 0.40 M km² for the same ensemble and scenario. Where a layer needs such a count,
compute it from the **member share** (the fraction of members over the threshold), which is
invariant to the central statistic, and publish the agreement spread rather than one number.

### `percentile`

Each cell's decadal `median` is ranked against the shared 2020s baseline land
distribution — one global reference, so values are comparable across decades and
scenarios. Clipped to [1, 100].

- **Inverted** (`101 − raw`) when the layer declares `percentile_direction:
  higher_is_better` — for assets like stored carbon the risk is *loss*, so a low stock
  earns a high risk percentile.
- **Two-tier** when the baseline is >2% exact zeros: zeros → 1, non-zeros ranked against
  the non-zero baseline → [2, 100].

### The two slopes

Both are emitted because **they fail in opposite regimes**, and neither alone is safe:

| Regime | `ols_slope` | `sen_slope` | Read |
|---|---|---|---|
| Well-behaved continuous field | correct | correct | either |
| **Zero-inflated hazard** (`led`, `driedarea`, `let`) | correct | **collapses to exactly 0** — most year-pairs are 0→0; measured at 91.3% of `driedarea` ssp126 cells and **96.9% of `let`'s exposed land** | `ols_slope` |
| **Unbalanced members with level offsets** | **biased** — measured +40% (0.70 vs true 0.50) when masking is uneven | correct | `sen_slope` |
| **Multimodal ensemble** (`permafrost-3b`) | correct — reproduces the ensemble-mean member trend exactly (+0.0326 dec⁻¹ against a member mean of +0.0326) | **biased LOW** — +0.0069 dec⁻¹, *below every one of the 12 members* (range +0.0104…+0.0826), because the pairwise sample is dominated by cross-cluster pairs carrying the level offset rather than the trend | `ols_slope` |
| Outliers / a single wild year | pulled off | robust | `sen_slope` |

Units are **per year**. Multiply by 10 to declare per-decade, and record which in
`slope_units` — do not leave it implicit. (Fitting against a decade index and *also*
multiplying by the window length inflates every trend 10×; that mistake has been made
here before.)

**The baseline panel is NaN, not 0.** Writing 0 makes the entire ocean a finite zero,
and the QA report will not catch it because it only asserts that *finite* baseline
trends equal zero.

A slope is fitted from each cell's own finite observations — a cell need not be present
in every member. Cells with fewer than 3 finite observations, or with all observations
in a single year, get NaN rather than a spurious trend.

**A layer whose mask is time-varying must mask its slopes to each decade's median mask.**
The slope window *expands* from the baseline while the median window is the decade alone,
so a cell with observations early in the window but none inside the decade gets a finite
slope against a NaN median — which the mask-agreement check rejects, correctly: a trend
over a decade the subject was absent from is not a trend. Measured on `npp-tempnle`
(a 2%-cover conifer presence mask): 53 such cells at the 2030s rising to 374 by the 2090s,
and they were the artifacts — dropping 53 of 25,821 moved the mean slope from **−1.89 to
+0.64**, because a stand vanishing mid-window produces a wild trend.

**The "opposite regimes" premise fails on a CENSORED field, and then agreement means
nothing.** The table above is the reason both slopes are emitted: each is wrong where the
other is right, so disagreement flags a fragile cell. A field pinned at a **bound** breaks
that. Where every observation in the window is the maximum, the OLS slope and the Theil-Sen
slope are both ~0 — not because the trend is zero but because there is no headroom left to
measure — and the two **agree**. Agreement near zero is therefore ambiguous on such a layer
between "no trend" and "maximally exposed, permanently".

Measured on `heatwave-3b` (2026-08-14), a binary exposure flag defined against each cell's
own preindustrial distribution: 45.9% of ssp585 2090s cells sit at exactly 1.0, and there
the CI collapses to zero width and the percentile ties at 100 (51.9% of cells at ≥99.5).
The censoring **inverts trend rankings between regions**: the Amazon's `ols_slope` *falls*
+0.160 → +0.046 dec⁻¹ as it saturates 0% → 100%, while never-saturating Siberia *rises*
+0.069 → +0.098, so the final panel has Siberia out-trending the Amazon 2.1×.

A layer in this regime must declare it. Emit the contract fields unchanged — do **not**
invent a statistic to hide the ceiling — and add a `saturation_caveat` naming how a
saturated cell is identified from the published fields alone (`median` at the bound with a
zero-width CI), plus per-panel shares. The bound need not be 1: any layer with a physical or
definitional ceiling is a candidate, so check the share at the bound before reading a slope.

**Judge slope agreement on ACTIVE cells only.** A cell that is permanently 0 — never
burns, never sees a cyclone — has a genuinely zero slope under *both* estimators, so
including it inflates apparent agreement and dilutes the Sen zero-fraction. On `let` the
all-cell view reads 73% sign agreement / 99.2% Sen-zero while the active-cell view
(either slope non-zero) reads **3.0% / 97.0%** — opposite conclusions from the same
array. `test_shared_baseline.py` reports the active-cell figure and the all-cell figure
in parentheses. This is the same dilution that once turned a `sen_slope == 0` share of
66–76% into a reported 20.7% by counting ocean.

## Shared 2020s baseline

The 2020s panel of `median` / `lower_ci` / `upper_ci` / `percentile` is **bit-identical
across scenarios**: each member's 2020s window is averaged across scenarios first, then
pooled. Valid only when ensemble composition is uniform across scenarios — if a member
is missing from a scenario, pool that scenario's baseline over *its own* members and
declare `members_by_scenario` (member **identity**, not a count).

## Known consequences

Stacking `(year, member)` as independent observations is the agreed estimand. What it
implies, so nobody rediscovers it as a bug:

- Members with more valid observations carry more weight.
- Interannual variability and inter-model spread are conflated in the CI.
- Cross-member pairs mix model level offsets with temporal change. A symmetric, complete
  ensemble largely cancels this in the Sen median; **uneven missingness does not**, and
  that is what biases `ols_slope`.
- On `csoil-total` the between-member level offset measured **68.7×** the interannual SD,
  so this is a live regime, not a hypothetical.

## This contract is about FORM, not meaning

`test_shared_baseline.py` passing means the file is *shaped* right — schema, shared
baseline, CI ordering, percentile range and orientation, slope masks, ensemble depth. It
cannot tell you the input is about what its name says.

Measured 2026-08-11: both sugarcane layers passed every check and were meaningless —
ISIMIP2b LPJmL simulates no sugarcane in the sugarcane belt, so the layer was zero across
São Paulo, Uttar Pradesh, Queensland and Florida while marginal temperate cells read as the
maximum. Once a layer passes this contract, the entire remaining risk sits in the input.
Check it per **GUARDRAILS §12** — verify the field is non-trivial at named reference
locations where the subject demonstrably exists, and record those sites and values in the
processor docstring.

## Performance

Theil-Sen is the only expensive term (quadratic in stacked points); everything else is
linear. A 12-member × 80-year expanding window is 455,040 valid pairs per cell.

| | single-core | 10 cores |
|---|---|---|
| one layer (~67k land cells × 24 panels) | 219 min | **21.9 min** |
| 21-layer backfill | — | 7.7 h |

Pair subsampling is available (`max_pairs=`) but **off by default**: at csoil's 68×
offset ratio a 100k cap costs ~15% of the slope, and the only tolerable cap buys just
1.14×. Use it for iteration, not production.

## Appendix: what changed

| | family A (`process_qg.py`) | family B (2026 layers) | this spec |
|---|---|---|---|
| `median` | mean over models, median over years | mean over years, mean over members | median over the pooled sample |
| CI | IQR pooled over member × year | mean ± 1 SD across members | IQR (mean ± 1 SD if boolean) |
| trend | OLS within each decade | baseline-anchored two-point rate | `ols_slope` + `sen_slope`, expanding window |
| significance | declared in `CLAUDE.md`, never emitted | never emitted | superseded — read the two slopes |

`CLAUDE.md` previously declared a 6th value class, `significance`, that no processor
ever wrote; `export_formatter.py` resolved `Decadal_Trend_Significance` to NaN, which
read as "not significant" rather than "not computed." Carrying both slopes replaces it:
disagreement between `ols_slope` and `sen_slope` is the signal that a cell's trend is
not robust.
