# Research method

Everything in a bespoke report that is not computed from the delivery is researched, and
every researched claim carries a citation that resolves to `dossier.yaml`. This file is how
the dossier gets filled.

## The rule

**Never assert anything you have not read.** Not "it is well known that", not "typically",
not a plausible number recalled rather than looked up.

The failure mode is not laziness, it is **fluency**. A confident paragraph about a customer's
operations reads identically whether it was researched or invented. By the time it reaches a
disclosure nobody can tell which, and the customer is the one who signed it. That is why the
citation check is machinery rather than a habit: `check_narrative()` refuses to render an
uncited paragraph in a cited slot, and refuses a citation that does not resolve.

An unresolvable citation is **worse** than no citation, because it looks like evidence.

## Source precedence

Prefer the highest tier that answers the question, and say in the narrative when you are
relying on a lower one.

| Tier | `source_type` | Notes |
|---|---|---|
| 1 | `customer-direct` | Something the customer told us. Record who, when, and in what forum. |
| 2 | `customer-published` | Their annual report, website, regulatory filing. |
| 3 | `regulator` | A filing, register, or official dataset (USDA, EPA, a national statistics office). |
| 4 | `peer-reviewed` | Published research. |
| 5 | `trade-press` | Industry reporting. Useful for orientation, weak for facts. Never the only support for a load-bearing claim. |

**Pull from the customer wherever possible.** Their own account of their operations beats any
inference, and asking is cheap. Where a conversation is the source, write down what was said,
by whom and when — an unrecorded conversation is not a source.

## Temporal lag is a field, not a feeling

Two dates, always:

- `as_of` — when the fact was true
- `retrieved_on` — when we read it

A 2019 figure retrieved today is a 2019 figure. If the lag matters to the claim, say so in
the narrative. The example dossier records a 2004 USDA publication and notes explicitly that
it describes practice as of the early 2000s and does not reflect later change — that
sentence is the difference between a citation and an alibi.

## Double-check before you write

`verified_by` is a required habit, not an optional field. Record **how** the fact was
checked and against what. Where a fact could not be independently verified, say that in the
field rather than leaving it blank — the example dossier does exactly this for a relative-
mortality claim, and consequently the narrative states the *ranking* (weather exceeds fire,
insects and disease) without restating proportions that rest on a single source.

Assess provenance actively:

- Does the source actually say what you are citing it for, or does it cite someone else?
- Is it about the same geography, species, asset class and period as your claim?
- Is there a more recent version?
- Would the customer's own team recognise it as authoritative?

## The dossier

`deliveries/<customer>/<date>/dossier.yaml`, scaffolded by
`generate_bespoke_report.py --scaffold`.

```yaml
sources:
  - id: usda-srs-gtr075-pine-silviculture   # cited as [dossier:usda-srs-gtr075-pine-silviculture]
    title: "..."
    publisher: "..."
    url: "https://..."
    source_type: regulator
    retrieved_on: "2026-08-13"
    as_of: "2004"
    what_it_supports: >-
      The specific claim this source backs.
    verified_by: >-
      How it was double-checked, and against what.
```

`id` is stable and semantic — it appears in `narrative.md` and in the rendered reference
list, so renaming one breaks citations.

**One source, one claim.** A source listed without `what_it_supports` cannot be checked by
anyone else, which defeats the purpose of listing it.

A delivery with no external sources is legitimate — it means every statement derives from the
delivered climate data, and the report says so.

## Facet profiles

Profiles carry their own `sources:` front matter under the same rules. A profile making a
factual claim about an asset class or region needs the source attached; a profile giving
structural guidance ("this reader is deciding X") does not.

To cite a profile's source in a narrative, copy it into that delivery's `dossier.yaml` —
citations resolve against the dossier only, deliberately, so the report's reference list is a
single list a reader can check rather than a union across files they cannot see.

See [../profiles/README.md](../profiles/README.md) for the profile schema.

## Researching a new customer

1. **Ask them first.** What do they own, where, on what horizon, for what decision, and what
   do they already file?
2. **Read what they publish.** Annual report, sustainability report, website, filings. Record
   what you use.
3. **Fill the asset and region context** from tier 3–4 sources, not from memory.
4. **Write down what you could not establish.** A gap you name is a question for the
   customer; a gap you paper over is a defect in a document with their name on it.
5. **Seed a company profile** under `docs/reporting/profiles/company/` if this is likely to
   be a repeat engagement.
