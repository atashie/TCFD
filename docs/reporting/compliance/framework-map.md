# Framework map

What the assessment can evidence, per framework, and what it cannot. Appendix A of every
compliance report renders a customer-facing version of this; the extra material here is the
reasoning and the dates.

## Why IFRS S2 is the spine

TCFD was **disbanded in 2023** and the IFRS Foundation took over monitoring. IFRS S2 carries
its four pillars forward and is what the successor regimes point at:

- **California SB 261** (Climate-Related Financial Risk Act) names the IFRS Sustainability
  Disclosure Standards as an accepted framework. Companies over USD 500m revenue doing
  business in California file a climate-related financial risk report biennially. *Status as
  of 2026-08-13: the original 1 January 2026 deadline has been subject to litigation and
  CARB has indicated an alternate date would follow — **confirm the current deadline before
  relying on it**.*
- **CDP** has aligned its corporate questionnaire with IFRS S2.
- **ESRS E1** asks many of the same questions in its own vocabulary.
- Mandatory IFRS S2 adoption began in several jurisdictions from 1 January 2026 (Brazil for
  publicly accountable entities, Chile, Mexico, Qatar among them).

One document on the S2 spine with a mapping appendix therefore serves several filings. Four
documents on four spines would drift apart at the first data refresh.

## IFRS S2 — what we supply

| Paragraph | Requirement | Status |
|---|---|---|
| 10, 10(a) | Risks that could affect prospects; physical or transition | **Supplied** — physical only |
| 10(d) | Definitions of short, medium, long term | **Supplied**, but entity-owned; ours are defaults and disclosed as such |
| 13(a)–(b) | Effects on business model; where risks concentrate | **Partly** — concentration by site; effects are entity-owned |
| 14 | Response and adaptation | **Not supplied** — entity-owned |
| 15–21 | Effects on financial position, performance, cash flows | **Not supplied** — needs vulnerability functions and asset values |
| 22(a) | Assessment of the resilience of the entity's **strategy and business model** — implications, significant uncertainties, financial flexibility, ability to redeploy assets | **NOT SUPPLIED.** We compare hazard exposure across scenarios; that is an input to a resilience assessment, not the assessment |
| 22(b) | Which scenarios and sources, over which horizons, on which assumptions | **Supplied** |
| 22(b)(i) | Inputs, incl. whether a diverse range was covered and whether one scenario was aligned with the latest international agreement | **Partly** — three forcing tiers covers the diverse range. Alignment is offered as a **proxy** only: a scenario code does not establish alignment, and the judgement is the entity's |
| 5–7 | Governance | **Not supplied** — entity-owned |
| 24–26 | Risk management process | **Not supplied** — entity-owned |
| 28–29 | Metrics; basis and limitations | **Supplied** |
| 29(c) | Amount **and** percentage of assets vulnerable to physical risk | **NOT REPORTED** — method deferred, see [vulnerability-definition.md](vulnerability-definition.md) |
| 33–37 | Targets | **Not supplied** — entity-owned |

29(c) is the paragraph that gets scrutinised, and it is the one we currently do not answer.
Two separate obstacles: no climate model contains a carrying amount, and — the larger one —
converting exposure into vulnerability is a judgement we have not made. The report states
that position rather than publishing a provisional count. Section 8 likewise prints the whole
"not supplied" column rather than omitting it.

## CDP question 3.1.1

CDP's per-risk table is understood to ask around ten fields per physical risk, with
management credit contingent on time horizon, likelihood and magnitude being populated, and
leadership credit on scenario analysis under at least two pathways. **All of that is
unverified against the live form** — see the warning below. On the scenario count this
assessment would comfortably satisfy a two-pathway requirement, since it uses three.

| Field | Supplied | Note |
|---|---|---|
| Primary climate risk driver | yes | from the hazard list |
| Value-chain location | partly | **direct operations only**; upstream and downstream are not assessed |
| Country / area | yes | from `locations.csv` |
| Time horizon | yes | the configured horizons |
| Likelihood | **proxy only** | percentile is an exposure level, not an event probability |
| Magnitude of impact | **proxy only** | a percentile band, not a financial magnitude |
| Financial impact figure | no | entity-owned |
| Primary financial effect | no | entity-owned |
| Explanation | partly | method here; business consequence entity-owned |
| Cost of response / management method | no | entity-owned |

**The likelihood and magnitude rows are the honest weak point.** CDP wants a probability and a
financial magnitude; we have a spatial exposure ranking. Mapping one onto the other is a
judgement the customer must make and own, and the report should say so rather than
manufacturing a correspondence.

**THE WHOLE CDP MAPPING IS UNVERIFIED, not just the hazard labels.** As of 2026-08-13 the
live 2026 questionnaire had not been opened: neither the hazard dropdown NOR the per-risk
field list above NOR the scoring rules ("ten fields", "management credit", "leadership tier")
were transcribed from it. They are inferences from published guidance and secondary summaries.
An earlier version of this file marked only the labels as unverified, which was too narrow and
implied the field list had been checked. Open the questionnaire and reconcile the whole
appendix before filing. GUARDRAILS §11 — a recorded claim needs its receipt or the word
UNVERIFIED.

## ESRS E1-9

| Datapoint | Supported |
|---|---|
| Monetary amount **and proportion** of assets at material physical risk (para 66(a)) | **No.** This is the datapoint, and it is financial. We supply neither half: the definition of "at material risk" is not agreed, and no asset values were supplied |
| Disaggregation of those monetary amounts into acute and chronic | **No.** We classify the *hazards* acute/chronic, which is not the same thing as the disaggregated monetary figure |
| Breakdown over short / medium / long term | **Partly** — exposure is reported at three horizons; the monetary datapoint those horizons would qualify is absent |
| Before adaptation actions | n/a until the above exists. No adaptation is modelled anywhere, so any future figure would be a pre-adaptation one |
| Location of significant assets *at material physical risk*, by NUTS 3 | **No** — applies to assets already determined to be at material risk. Coordinates are delivered; no NUTS lookup exists here |
| Anticipated financial effects | **No** — entity-owned |

**E1-9 is a financial disclosure, and this assessment is not one.** Its core requirement is a
monetary amount and proportion of financial-statement assets, disaggregated acute/chronic and
broken down by horizon. Classifying hazards and printing horizon labels does not supply it,
and an earlier version of this file marked three of those rows "Yes" on that basis — which
overstated coverage of the one datapoint E1-9 actually asks for.

The NUTS 3 gap is worth flagging early to any EU customer: it is a lookup we do not have, and
inventing a mapping from decimal degrees would be a fabricated regulatory identifier. Note it
applies to assets already determined to be at material risk, so it sits downstream of the
definition we have not yet agreed.

## Hazard coverage

**Not restated here.** `config/hazard_taxonomy.yaml` is authoritative: it enumerates the
hazard families a physical-risk disclosure is expected to address, which of them each
registry layer evidences, and for each gap what would close it. Coverage is read from the
file — any count written here would go stale the moment a layer ships.

That file also carries the customer-facing / internal split. `customer_note` is rendered into
reports; `materiality_note`, `blocker` and `isimip_candidate` are ours and must never be —
they carry repository paths, dataset defects and the word UNVERIFIED, which mean something
precise to us and something alarming in a filing.

## Keeping this current

Disclosure requirements move faster than the code. Before any delivery is filed:

1. Re-check the framework version and deadline for the customer's jurisdiction. The SB 261
   date in particular has already moved.
2. Enumerate the CDP dropdown if filing against CDP, and update the taxonomy labels.
3. Re-read 29(c) and E1-9 for changes to the vulnerability metric — that definition is where
   an amendment would hit hardest.

Date anything you verify, here and in the taxonomy. An undated framework claim is the same
kind of liability as an undated absence claim.
