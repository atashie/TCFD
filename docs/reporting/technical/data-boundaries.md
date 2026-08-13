# Data boundaries — what a report may and may not say

Read this before writing any sentence about what a number means. Everything here is a
property of the data, not a house style, and most of it has already been got wrong once.

## The single most misread quantity: `percentile`

**`percentile` ranks a site against the global 2020s land distribution.** It is a LEVEL, not
a change. Percentile 80 means "more exposed than 80% of global land was in the 2020s". It
does not mean the site has worsened by 80%, and it does not mean the site has worsened at
all.

| Say | Never say |
|---|---|
| "more exposed than 80% of global land in the 2020s" | "an 80% risk" |
| "at the 80th percentile, and rising" | "risk has increased to 80" |
| "ranks 3rd of 5 sites on wildfire exposure" | "3rd most likely to burn" |

A site can sit at the 95th percentile and be almost flat; another can climb steeply and
remain less exposed in absolute terms. **Level and trend are separate axes and must be
reported separately.** The compliance report has separate sections for them for this reason.

Percentile is also **already oriented for risk** on every layer: 100 is worst everywhere,
because `higher_is_better` layers were inverted at processing time. Never invert again.

## There is no p-value, and disagreement is the signal

Two slope estimators are fitted because they fail in opposite regimes. `sen_slope` collapses
to exactly zero on zero-inflated hazards; `ols_slope` absorbs member-level offsets as trend
when coverage is uneven. Which one to read is **measured per layer** and recorded in
`recommended_slope`.

Where the two disagree, the trend is not robust. Say that. Do not substitute a confidence
statement the data does not support — there is no significance test in this contract, and the
retired schema's significance columns silently resolved to "not significant" for every layer
ever delivered, which is exactly the failure this design removes.

## Resolution: what "at this site" actually means

- The model grid is **0.5°** (~55 km at the equator).
- Extraction blends the **four surrounding cell centres** with a Gaussian weighting —
  20,000 random sites, 100% resolve to a 4-cell blend. Footprint ≈ **1° × 1°, 111 km
  north–south**.
- `cyclone` is *additionally* smoothed at processing time, so its two stages compound. Its
  values are the least site-specific numbers in any delivery.
- Moving a site **0.25°** changed 2090s burnt-area values at the example portfolio's own
  sites by **44% to 569%**.

**Both figures are reproducible**, and that is deliberate — an unreproducible measurement in
a filing is indistinguishable from an invented one:

```bash
python scripts/measure_extraction_sensitivity.py
```

Therefore: this is a screening instrument for comparing sites, horizons and scenarios. It is
not a site-level engineering assessment, it cannot resolve topography, drainage, defensible
space or flood defences, and a coordinate is an input worth confirming before anything is
decided on it.

## The claim list

### Not currently reported at all

- **A count or percentage of assets vulnerable to physical risk.** The method is deferred —
  see [../compliance/vulnerability-definition.md](../compliance/vulnerability-definition.md).
  `percentile` is a global-relative exposure rank and "vulnerable" is a statement about
  susceptibility to harm; nothing here connects them. The verifier fails any report that
  publishes the figure while the entry sits in `TBD_SECTIONS`.

### Never claim

- An expected financial loss, a loss probability, a return impact, or a valuation effect.
  We have hazard exposure. Converting it to money needs vulnerability functions and asset
  values we do not have.
- That a hazard is immaterial because it does not appear in the report. Sixteen of nineteen
  hazard families are not assessed at all — see `config/hazard_taxonomy.yaml`.
- That a site is low-risk because a value is missing. `NOT_ASSESSED` is not
  `NOT_VULNERABLE`, and the two are separate statuses in `vulnerability_frame()` precisely so
  they cannot be conflated.
- That a wildfire number reflects this site's fire risk. The layer contains no suppression
  capacity, no firebreaks and no fuel management.
- That a drought number reflects aridity. It measures departure from a fixed historical
  reference, so a permanently dry site can read low.
- That rising productivity is a benefit, without the CO₂-fertilisation caveat.
- That the portfolio is diversified against climate risk. Site-level results treat each
  asset independently and cannot see correlated regional loss.

### Claim only with the qualification attached

| Claim | Required qualification |
|---|---|
| any portfolio percentage | that it is count-weighted, unless asset values were supplied |
| a cross-tier comparison | the balanced panel it was computed over |
| a comparison of two sites | that both were assessed on the same hazard set |
| a cyclone result | single impact model, wind only, no high-forcing scenario |
| a trend | which estimator, and whether the two agree |

## Composition traps, both measured

**Decade.** ISIMIP3b layers publish no 2010s panel, so a 2010s Climate Score can rest on one
hazard where the 2020s rests on three — measured 30.4 → 39.9 on the example portfolio, a 31%
jump entirely manufactured by hazards arriving. `portfolio_score_series()` drops incomplete
decades for this reason.

**Tier.** `cyclone` publishes no `rcp85`, so cyclone-carrying assets are unscoreable at high
forcing. Averaging each tier over whatever assets it has made the high-tier baseline read
39.9 against 42.1 — impossible on a common basis.

The invariant that catches both: **the 2020s panel is bit-identical across scenarios, so a
baseline Climate Score must be equal across tiers on a balanced panel.**
`assert_baseline_tier_equality()` enforces it and both report builders call it.

Apply the balanced panel at **every** level of rollup. Fixing it at the portfolio level did
not fix it at the location level — the same artifact reappeared one level down.

## Where the boundaries are enforced rather than remembered

- `vulnerability_frame()` — hazard layers only, `NOT_ASSESSED` distinct, max not mean
- `assert_baseline_tier_equality()` — refuses to publish a composition artifact
- `assert_caveats_present()` — refuses a report missing a required disclosure
- `check_narrative()` — refuses an uncited claim or an unresolvable citation
- `customer_evidence()` — drops evidence pointing at files the customer does not hold
- `test_customer_delivery.py` — recomputes the headline metric and compares it to the
  printed table
