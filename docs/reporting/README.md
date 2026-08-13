# Reporting — the map

Stages 3 and 4 of the customer delivery. **The `/customer-delivery` skill is the entry
point**; this tree is what it points at. Read the skill first, then whichever file below
covers what you are about to do.

| Read | Before you |
|---|---|
| [technical/report-anatomy.md](technical/report-anatomy.md) | change what a report contains or add a section |
| [technical/figures-and-output.md](technical/figures-and-output.md) | add a figure, touch the print CSS, or ask about PDF |
| [technical/data-boundaries.md](technical/data-boundaries.md) | write any sentence about what a number means |
| [compliance/framework-map.md](compliance/framework-map.md) | claim the report satisfies a disclosure requirement |
| [compliance/vulnerability-definition.md](compliance/vulnerability-definition.md) | quote or change the vulnerable-asset count |
| [research/method.md](research/method.md) | look anything up about a customer, a region or an asset class |
| [profiles/README.md](profiles/README.md) | write or edit a facet profile |
| `config/hazard_taxonomy.yaml` | say what this assessment does or does not cover |

## Two documents, one set of numbers

| | Compliance report | Bespoke report |
|---|---|---|
| Spine | IFRS S2 paragraph order | the reader's decision |
| Author | fully generated | generated data + a written narrative |
| Reproducible byte-for-byte | yes | no, and it should not be |
| Command | `generate_compliance_report.py` | `generate_bespoke_report.py` |

Both read the same delivery through `scripts/utils/report_common.py`, so **they cannot
disagree about a number**. Every count, figure and caveat comes from one implementation. If
you find yourself computing something in a report script, it belongs in `report_common.py`
instead — that is not a style preference, it is the only thing keeping two documents about
one delivery consistent.

## The stage order, and why caveats comes first

```
inputs → extract → dashboard → caveats → compliance_report → bespoke_report
```

**Stage 4 runs before Stage 3.** The caveat set is a mechanical derivation from the manifest,
the CSVs and the hazard taxonomy, and it is an *input* to both reports: each is required to
carry every `must_disclose` entry, `report_common.assert_caveats_present()` refuses to render
without them, and `test_customer_delivery.py` re-checks afterwards. Generating caveats last
would mean each report derived its own list, and the one thing a caveat list must not do is
differ between two documents describing one delivery.

## Say "not yet decided" rather than fill the box

**A report section is allowed to report nothing, and should, wherever the method for
translating our data into a framework's requirement has not been thought through and agreed
with the user.**

The pressure runs the other way. A complete-looking report is what a customer expects to
receive, and that pressure is exactly what produces a number nobody chose deliberately —
which, once it is in a filing, is indistinguishable from a reasoned one.

Deferred decisions live in `report_common.TBD_SECTIONS` and render through `tbd_block()`,
which states the requirement, why it is deferred, the decisions outstanding, and what is
reported instead. Both reports render the same block, so they cannot describe the same gap
differently. The verifier enforces the deferral: while an entry is in `TBD_SECTIONS`, a
report that publishes the corresponding figure **fails**.

Currently deferred: **the vulnerable-asset count** (IFRS S2 29(c) / ESRS E1-9). See
[compliance/vulnerability-definition.md](compliance/vulnerability-definition.md) for what was
tried, why it was wrong, and what has to be settled.

## The composition model

A report should be specific to one **asset × region × persona × vertical × company × use
case**. That is combinatorial: five asset types, ten regions, four personas, five verticals
and four use cases is four thousand documents.

So the library holds **one short profile per value of each facet** — `timber-land-loblolly`,
`us-southeast`, `sustainability-team` — and a report selects one from each. Sixty files
cover a million combinations. Facet selections live in the delivery's `report_config.yaml`.

**Every facet accepts a list**, and most portfolios need one. The worked example carries
three asset profiles and four region profiles, because it holds timber, warehouses and a data
centre across California, the Gulf coast, the mid-Atlantic and Bavaria.

**Region is the one facet the DATA determines, and it is validated.** Each region profile
declares a `matches:` block of countries and states; `assert_region_coverage()` checks every
delivered location against the selected profiles and **refuses to build** if any is
uncovered. A report whose regional context does not match the sites it describes is
confidently wrong about the thing a reader checks first — and describing a portfolio spanning
two continents with one region's framing is the specific failure this prevents. A region
profile without a `matches:` block is a load error, because it would match everything and
silently defeat the check.

`company` is different in kind. It is researched per engagement and lives in the delivery's
`dossier.yaml`; a repeat customer earns a company profile seeded from that dossier.

### Profiles guide the narrative. They are never pasted into it.

This is the load-bearing rule. If profile prose were rendered into the report, every loblolly
report would contain the same paragraphs — generic output that *looks* tailored, which is
worse than obviously generic output. So a profile contains no sentences meant for a customer.
It contains what the writer needs in order to write, and the scaffold surfaces it as HTML
comments that are stripped before rendering.

`report_common.parse_narrative()` strips comments from every slot body, and the verifier
fails any report containing `<!--`. That check exists because a line-by-line comment filter
catches only the opening `<!--` and renders the rest — which would put the profiles' own
"Do not claim" lists into the customer's document.

## The two-author boundary, enforced in code

A report may contain exactly two kinds of sentence:

1. **Derived from the delivery** — computed in `report_common.py`, never retyped.
2. **Researched** — carrying a citation that *resolves*.

There is no third kind. `[data:values.csv#AST-001/wildfire/ssp370/2050]` must name a row that
exists; `[dossier:some-source]` must name a source in `dossier.yaml`. An unresolvable
citation fails the build.

The failure mode this guards against is not laziness, it is **fluency**. A confident
paragraph about a customer's operations reads identically whether it was researched or
invented, and by the time it reaches a disclosure nobody can tell which. An unresolvable
citation is worse than no citation, because it looks like evidence.

## Running the stages

```bash
# Stage 2 + 4 + 3a in one command. Caveats runs before the compliance report.
python scripts/generate_customer_delivery.py --customer "<name>" --input <sites.csv> --run --reports

# or separately
python scripts/generate_delivery_caveats.py      deliveries/<c>/<d>
python scripts/generate_compliance_report.py     deliveries/<c>/<d>

# Stage 3b: pick facets in report_config.yaml first, then
python scripts/generate_bespoke_report.py deliveries/<c>/<d> --scaffold   # writes narrative.md + dossier.yaml
python scripts/generate_bespoke_report.py deliveries/<c>/<d>              # renders, once they are filled

# always
python scripts/test_customer_delivery.py deliveries/<c>/<d>
```

The bespoke report is deliberately **not** chained onto `--run`: it needs facets chosen and a
narrative written, neither of which a batch run can supply.

Re-running the extract marks any existing caveats or reports **stale** in the manifest, and
the verifier refuses a stale artifact — a report built against an earlier extract no longer
describes the data sitting next to it.

Two measurements quoted in reports are reproducible on demand, and should be re-run after any
change to the extraction parameters:

```bash
python scripts/measure_extraction_sensitivity.py
```

## What the verifier checks about reports

Every `must_disclose` caveat appears in every built report · no unfilled narrative slot · no
uncited paragraph in a cited slot · every citation resolves · no HTML comment survived · no
internal vocabulary leaked · the build stamp is present · a stage claiming to be built has
its file · **and the vulnerable-asset count printed in the compliance report is recomputed
independently and compared**. That last one is the end-to-end proof for the headline metric;
everything else about a report can be right while the table says something else.

All of it has been shown to fail on injected corruption. A check nobody has watched fail is a
check nobody should trust.
