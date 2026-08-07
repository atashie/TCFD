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
`decadal_statistic` global attribute (`pooled_median` or `pooled_mean_boolean`).

### Boolean detection

Taken from the **values**, never the name (GUARDRAILS §9). A field is boolean only if
every finite value is exactly 0 or 1. `led` is binary; `let` is a continuous fraction in
[0,1) despite the sibling naming. A *zero-inflated* continuous field is **not** boolean
and takes the median/IQR branch. The scan is exhaustive — a sampled scan can step over a
rare continuous tail and silently switch the published statistic.

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
| **Zero-inflated hazard** (`led`, `driedarea`) | correct | **collapses to exactly 0** — most year-pairs are 0→0; measured at 91.3% of `driedarea` ssp126 cells | `ols_slope` |
| **Unbalanced members with level offsets** | **biased** — measured +40% (0.70 vs true 0.50) when masking is uneven | correct | `sen_slope` |
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
