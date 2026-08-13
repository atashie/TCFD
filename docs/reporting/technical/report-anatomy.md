# Report anatomy

## Compliance report — `generate_compliance_report.py`

IFRS S2 spine, mapped outward. **Fully deterministic**: every sentence is either fixed text
about method or a number computed from the delivery. No narrative slots, no researched
claims, so an auditor can rebuild it and get the same document.

| § | Section | Content |
|---|---|---|
| 1 | Scope and basis | portfolio size, layers, and an explicit is/is-not callout |
| 2 | Risks identified | hazards assessed, acute/chronic class, ensemble, scenarios |
| 3 | Time horizons | the three configured horizons; a callout if they are ours |
| 4 | Where risk concentrates | ranked bar + band stack + per-asset table |
| 5 | Resilience and scenario analysis | tiers, Paris-aligned member, provenance, Climate Score series |
| 6 | Assets vulnerable | **DEFERRED** — states the requirement, why, and what is outstanding |
| 7 | Direction of change | one trend strip per hazard — units differ, so never one axis |
| 8 | Sections the entity must complete | governance, risk management, financial effects, targets |
| 9 | Limitations | every `must_disclose` caveat, then the `should_note` table |
| 10 | Hazards not assessed | the 16 uncovered families with customer-facing notes |
| A | Framework mapping | IFRS S2, CDP 3.1.1, ESRS E1-9, SB 261 |
| B | Per-asset results | the full table at three horizons × three tiers |

**Section 8 is not padding.** IFRS S2 has four pillars and this report can populate parts of
two. Publishing sections 1–7 without naming the gap would present hazard screening as an S2
filing. A missing pillar is invisible unless it is named.

**Section 6 currently reports nothing, on purpose.** See
[../compliance/vulnerability-definition.md](../compliance/vulnerability-definition.md). It is
an instance of the general rule: a section whose method has not been agreed states that
instead of publishing a provisional figure.

**Section 10 is not padding either.** A report that lists what it assessed and stops reads as
though the rest was assessed and found immaterial — every number correct, the document still
misleading. This is the most dangerous thing this pipeline could do to a customer.

### Adding a section

Add to `SECTIONS` and write `sec_<anchor>()`. It must take the delivery and return HTML, and
it must not compute anything a `report_common` function could compute — the bespoke report
will need the same number.

## Bespoke report — `generate_bespoke_report.py`

Same shell, different spine: the reader's decision rather than the standard's structure.

| Slot | Heading | Cited | Guided by profile section |
|---|---|---|---|
| `executive_summary` | What we found | yes | What this reader needs decided |
| `portfolio_context` | Your portfolio in context | yes | Framing and vocabulary |
| `hazard_reading` | What these hazards mean for this asset class | yes | Transmission channels |
| `site_findings` | Site-by-site | yes | Metrics that lead |
| `decision` | What this changes for the decision | yes | What this reader needs decided |
| `confounders` | What would change this conclusion | yes | Known confounders |
| `next_steps` | What we recommend doing next | no | Questions for the customer |

Then three generated sections: **The numbers behind this**, **What we could not assess**,
**Sources and provenance** (numbered references first, then the source and climate-data
tables).

`next_steps` is the one uncited slot, because it is a request rather than a claim.

### Citations

```
[data:values.csv#AST-001/wildfire/ssp370/2050]
[data:climate_score.csv#AST-001/high/2090]
[data:caveats.json#SPATIAL-FOOTPRINT]
[data:layers.csv#cyclone]   [data:assets.csv#AST-003]   [data:locations.csv#LOC-005]
[data:manifest.json#counts.value_rows]
[dossier:usda-srs-gtr075-pine-silviculture]
```

Each must resolve. They are numbered in document order across all slots, rendered as
superscript links, and listed in a reference table — hover titles are not citations in a
document that will be printed and filed.

### Writing the narrative

1. `--scaffold` writes `narrative.md` and `dossier.yaml`. It surfaces every selected
   profile's guidance for each slot as HTML comments, plus the combined "Do not claim" list.
2. Write from the guidance. **Do not paste it** — see the composition rule in the
   [README](../README.md).
3. Rebuild. The build refuses on an unfilled slot, a `TODO`, an uncited paragraph in a cited
   slot, or an unresolvable citation.

Both files are never overwritten by `--scaffold`, so re-running it is safe.

## What both reports share

Loaded once in `report_common.py`: the delivery, the report config, the caveat set, the
coverage summary, the vulnerability frame, the balanced-panel Climate Score series, the
document shell, the print CSS and the build stamp. Both call
`assert_caveats_present()` before returning, and `assert_baseline_tier_equality()` before
drawing any cross-tier series.

## `report_config.yaml`

Written with defaults by the extract stage and **never overwritten afterwards**, so edits
survive regeneration.

```yaml
horizons:      {short/medium/long → decade, label}  + source: default|customer
vulnerability: {metric, threshold, sensitivity, basis} + source: default|customer
asset_values:  {supplied, n_assets_valued, currency, basis}
facets:        {asset, region, persona, vertical, use_case, company}
frameworks:    {spine, mapped}
```

`source: default` on the `horizons` block raises a `must_disclose` caveat. That is the point
of the field: IFRS S2 10(d) requires the *entity* to define its horizons. Change the value and
the source together — filling in `customer` without a customer actually saying so is the one
failure the field exists to prevent, exactly as with `confirmed_on` in the asset catalog.

**The `vulnerability` block is currently INERT.** The metric is deferred
([../compliance/vulnerability-definition.md](../compliance/vulnerability-definition.md)), so a
threshold sitting in that file is not an endorsement of one and nothing reads it for
publication. It stays because the machinery on both sides of the decision is already written:
remove the `vulnerability_metric` entry from `TBD_SECTIONS` and the block goes live, with the
verifier switching from "must not publish" to "the published counts must match an independent
recomputation".
