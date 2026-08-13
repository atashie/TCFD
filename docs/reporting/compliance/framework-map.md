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
| 22(a) | Climate resilience assessed via scenario analysis | **Supplied** |
| 22(b) | Which scenarios, sources, horizons, assumptions | **Supplied** |
| 22(b)(iii)–(iv) | Diverse range; a Paris-aligned scenario among them | **Supplied** — low/medium/high tiers, `rcp26`/`ssp126` as the Paris-aligned member |
| 5–7 | Governance | **Not supplied** — entity-owned |
| 24–26 | Risk management process | **Not supplied** — entity-owned |
| 28–29 | Metrics; basis and limitations | **Supplied** |
| 29(c) | Amount **and** percentage of assets vulnerable to physical risk | **Percentage always; amount only if asset values are supplied** |
| 33–37 | Targets | **Not supplied** — entity-owned |

29(c) is the paragraph that gets scrutinised, and it is the one our data cannot fully meet on
its own: no climate model contains a carrying amount. Section 8 of the report prints the
whole "not supplied" column rather than omitting it.

## CDP question 3.1.1

CDP's per-risk table asks ten fields per physical risk and awards management credit only when
time horizon, likelihood and magnitude are all populated with something other than "Unknown".
Leadership-tier credit requires scenario analysis under at least two pathways, which this
assessment satisfies (it uses three).

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

**The hazard labels are UNVERIFIED.** `config/hazard_taxonomy.yaml` maps our families onto
CDP's dropdown by name; as of 2026-08-13 the live 2026 questionnaire dropdown had not been
enumerated option by option. Open the questionnaire and reconcile before filing. This is a
GUARDRAILS §11 situation — a recorded claim needs its receipt or the word UNVERIFIED, and it
carries the word.

## ESRS E1-9

| Datapoint | Supported |
|---|---|
| Assets at material physical risk — count and share | yes |
| Assets at material physical risk — monetary amount | only with supplied asset values |
| Disaggregation by acute and chronic | yes — `esrs_class` in the hazard taxonomy |
| Location of significant assets by NUTS 3 | **no** — coordinates are delivered; no NUTS lookup exists here |
| Short / medium / long term breakdown | yes |
| Before adaptation actions | yes — no adaptation is modelled anywhere |
| Anticipated financial effects | no — entity-owned |

The NUTS 3 gap is the one worth flagging early to any EU customer: it is a lookup we do not
have, and inventing a mapping from decimal degrees would be a fabricated regulatory
identifier.

## Hazard coverage

**Not restated here.** `config/hazard_taxonomy.yaml` is authoritative: it enumerates the
nineteen families a physical-risk disclosure is expected to address, which of them each
registry layer evidences, and for each gap what would close it. Three families are covered
today.

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
