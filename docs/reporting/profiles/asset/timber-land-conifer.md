---
facet: asset
id: timber-land-conifer
name: Temperate conifer timberland
aliases: [timberland, timber land, conifer plantation, softwood, Douglas-fir, Norway spruce, managed forest]
catalog_entry: timber land
confirmed_on: null
sources:
  - id: usda-srs-gtr075-pine-silviculture
    title: "The Evolution of Pine Plantation Silviculture in the Southern United States (GTR-SRS-075, ch. 8)"
    publisher: USDA Forest Service Southern Research Station
    url: "https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs075/gtr_srs075-fox002.pdf"
    source_type: regulator
    retrieved_on: "2026-08-13"
    as_of: "2004"
    what_it_supports: >-
      Managed conifer rotation lengths are measured in decades, which sets the horizon
      against which any projection must be read. Figures are for southern pine; other
      species and regions differ and should be taken from the owner.
---

The general profile for managed temperate conifer stands. Where the holding is a specific
species-and-region system with its own documented damage pathway — loblolly in the US
Southeast, Norway spruce in Central Europe — use that profile as well, or instead.

## Transmission channels

Rank by how loss actually reaches the balance sheet, not by how alarming a hazard sounds.
The ranking is **region-dependent and species-dependent**, and getting it from the region
profile rather than assuming it is the point of keeping these separate:

- **Wind** — stem breakage and windthrow. Dominant where storms reach the stand; converts
  standing sawtimber into salvage in a single event, and salvage depresses the local price at
  the moment the owner is forced to sell. Our coverage is tropical cyclone wind only.
- **Fire** — total loss of a stand where it occurs, but mediated by suppression, access and
  fuel management, none of which our layer models.
- **Drought** — two channels on different timescales: growth loss within a rotation, and
  predisposition to pest and fire damage. The second is usually the larger and is the one our
  drought metric speaks to most directly.
- **Productivity** (`conifer-npp`) — not a hazard. The yield side, and the only layer that
  addresses return rather than risk of loss.

**Rotation length is the clock.** Managed conifer rotations run in decades
[usda-srs-gtr075-pine-silviculture], so a hazard that bites in the 2090s but not the 2050s is
a hazard for a later owner. Establish the rotation and the hold period before deciding which
horizon leads.

## Metrics that lead

- Level (`percentile`) and direction (the recommended slope) reported separately, never
  merged.
- The hazard the region profile identifies as dominant goes first, even if another ranks
  higher on the global percentile scale.
- The Climate Score is a screening aid for ranking many holdings. For a single stand it
  hides more than it shows, because it averages productivity in with the hazards.

## Known confounders

- **CO₂ fertilisation sits inside the productivity number.** Rising modelled NPP is partly a
  response to elevated CO₂, its real-world magnitude is contested, and the models treat
  nutrient and water limitation differently. Never present rising productivity as a plain
  benefit.
- **No pest or pathogen module exists anywhere in this assessment.** For several conifer
  systems the dominant mortality pathway is drought-mediated insect attack, and modelled
  productivity can rise while stands are dying.
- **`conifer-npp` is reported per unit stand area behind a 2% cover mask.** A site with no
  modelled conifer stand returns "off layer mask", which is a statement about the model's
  vegetation map, not about whether trees grow there.
- **Management dominates at stand scale.** Species and provenance choice, stocking, thinning
  and rotation move outcomes more than between-scenario spread over a single rotation. A
  report implying climate is the controlling variable will be dismissed by anyone who manages
  stands.

## Framing and vocabulary

Rotation, stand, stocking, basal area, thinning, salvage. Use the owner's volume units —
cubic metres in Europe, board feet or tons in North America. Say "growth" rather than "net
primary productivity" once past the technical section.

Avoid "risk score" for the percentile; it reads as an actuarial probability of loss and is
not one.

## Questions for the customer

- Species, provenance, rotation length and current age-class distribution, per stand.
- Product mix, and any carbon, conservation or recreation income.
- Insurance and salvage arrangements, and mill capacity within haul distance.
- Recent loss events and how they were handled.
- Whether the holding runs to rotation or to a fund term, and when that term ends.

## Do not claim

- Do not state or imply an expected financial loss, loss probability or return impact.
- Do not describe wildfire exposure as this stand's fire risk.
- Do not present rising productivity as a net benefit without the CO₂ and pest caveats.
- Do not carry one region's hazard ranking into another. Read the region profile.
