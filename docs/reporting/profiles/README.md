# Facet profiles

One short profile per **value** of each facet. A bespoke report selects one from each and
composes them; six facets with ten values each covers a million combinations from sixty
files. See the composition model in [../README.md](../README.md).

```
profiles/
  asset/       what the thing is            data-center, timber-land-conifer, timber-land-loblolly, warehouse-logistics
  region/      where it is                  central-europe, us-gulf-coast, us-mid-atlantic, us-southeast, us-west-california
  persona/     who reads the report         sustainability-team
  vertical/    what business they are in    investment-manager
  use-case/    what decision they face      acquisition-screening, annual-disclosure
  company/     the specific customer        (none yet — seeded from a delivery dossier)
```

Loaded by `scripts/utils/report_profiles.py`; selected in a delivery's `report_config.yaml`
under `facets:`.

## A profile is written for the writer, never for the customer

**Nothing in a profile is rendered into a report.** The scaffold surfaces it as HTML comments
and the writer works from it; `parse_narrative()` strips comments and the verifier fails any
report containing `<!--`.

This is the rule the whole library depends on. If profile prose were pasted into reports,
every loblolly report would carry the same paragraphs — generic output that *looks* tailored,
which is worse than obviously generic output, because nobody catches it.

So do not write customer-facing sentences here. Write what someone needs in order to think:
which hazards actually reach this asset and how, what decision this reader is making, what
vocabulary lands and what misfires, what will make them distrust the numbers, what to ask
them, and what must never be claimed.

## Every facet accepts a list

`report_config.yaml` may name one profile or several per facet, and most real portfolios need
several. The worked example carries three asset profiles and four region profiles.

## Region is determined by the data, and is validated

A region profile MUST declare a `matches:` block. `assert_region_coverage()` checks every
delivered location against the selected profiles and **refuses to build** if any is
uncovered — because a report whose regional context does not match its sites is confidently
wrong about the first thing a reader checks. A profile without `matches.countries` is a load
error: it would match everything and silently defeat the check.

```yaml
matches:
  countries: [USA, United States, US]
  states: [California, Oregon, Washington]   # empty list = the whole country
```

Matching is: country must match; if states are listed the location's state must be among
them; a location with no state recorded matches on country alone, so an offshore or
country-level site can still be covered. Absent CSV columns arrive as the float `nan`, whose
`str()` is the truthy string `"nan"` — `_blank_safe()` handles that, and without it every
offshore site reads as uncovered everywhere.

## Schema

```markdown
---
facet: asset            # asset | region | persona | vertical | use_case | company
id: timber-land-conifer    # must equal the filename stem
name: Temperate conifer timberland
aliases: [...]          # optional
matches: {...}          # REQUIRED for region profiles, ignored elsewhere
confirmed_on: null      # null until a human signs it off
sources: [...]          # required for any factual claim; see ../research/method.md
---

## Transmission channels
## What this reader needs decided
## Framing and vocabulary
## Metrics that lead
## Known confounders
## Questions for the customer
## Do not claim
```

Use the subset that fits — an asset profile has transmission channels and no reader
decisions; a persona profile is the reverse. **Headings outside this list are a load error**,
not a warning: a typo'd heading is guidance nobody will ever see.

| Section | What belongs in it |
|---|---|
| Transmission channels | how each hazard actually reaches value, ranked by materiality, not by how alarming it sounds |
| What this reader needs decided | the question they are actually answering; write for their second reading, under pressure |
| Framing and vocabulary | their words; terms that misfire; what to define in place |
| Metrics that lead | which of our columns goes first, and which to keep out of the headline |
| Known confounders | what will make this reader distrust the report, and what genuinely limits it |
| Questions for the customer | the inputs that would most improve the next version |
| Do not claim | facet-specific overclaim traps. The combined list from every selected profile goes into the scaffold. |

## `confirmed_on`

`null` means **nobody has reviewed this profile**. Everything seeded 2026-08-13 is still
null, and both the scaffold and the build print which. Filling the date in to quiet the
warning, without someone actually approving it, is the one failure this field exists to
prevent — the same protocol as `confirmed_on` in `config/asset_catalog.yaml`.

## Writing a new one

1. Copy the closest existing profile; keep the section order.
2. Cut every sentence that could appear in any profile of that facet. Generic guidance
   produces generic reports.
3. Attach a source to every factual claim. Structural guidance ("this reader reports upward")
   needs none; "hurricane wind is the dominant acute forest risk in the southeast" does.
4. Write **Do not claim** last and make it specific. It is the section that does the most
   work, because it names the mistakes a fluent writer makes without noticing.
5. Leave `confirmed_on: null` and say so when you hand it over.

## General and specific profiles

Both are legitimate. `timber-land-conifer` covers managed temperate conifer stands generally;
`timber-land-loblolly` covers the specific species-and-region system with its own documented
damage pathway. Select the general one, or both, or the specific one alone.

Where a hazard ranking depends on region — wind dominates southeastern pine and is irrelevant
in the Pacific Northwest — the **region profile owns the ranking** and the asset profile says
so rather than assuming. This is why the two are separate facets.

## Open catalog question

The `timber land` entry in `config/asset_catalog.yaml` omits `cyclone`. That is right for the
Pacific Northwest and wrong for the southeastern US, and the catalog has no region dimension
to express the difference. It has not been changed globally, because doing so would attach a
hazard with no transmission channel to inland holdings. **Raise it before the next timber
delivery.**
