---
facet: asset
id: timber-land-loblolly
name: Loblolly pine timberland
aliases: [loblolly, loblolly pine, southern yellow pine, southern pine plantation, pine plantation]
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
      Scale and intensification of southern pine plantations, and the halving of rotation
      lengths as mean annual increment roughly doubled.
  - id: jof-hurricane-wind-risk-forests
    title: "Potential Hurricane Wind Risk to US Rural and Urban Forests"
    publisher: "Journal of Forestry 119(4)"
    url: "https://academic.oup.com/jof/article/119/4/393/6260883"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2021"
    what_it_supports: >-
      Hurricane wind is the dominant acute forest risk in the southeastern United States.
  - id: frontiers-loblolly-stem-breakage-mortality
    title: "Long-term tree mortality prediction following stem breakage in loblolly pine plantations"
    publisher: "Frontiers in Forests and Global Change"
    url: "https://www.frontiersin.org/journals/forests-and-global-change/articles/10.3389/ffgc.2025.1708548/full"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2025"
    what_it_supports: >-
      Weather-event mortality in southeastern US forests exceeds mortality from fire,
      insects and disease.
---

## Transmission channels

Rank these by how the loss actually reaches the balance sheet, not by how alarming the
hazard sounds.

- **Wind (tropical cyclone).** Stem breakage and windthrow. This is the dominant acute loss
  mechanism for southeastern pine and the evidence is unambiguous: weather-event mortality
  exceeds fire, insect and disease mortality combined in this region
  [frontiers-loblolly-stem-breakage-mortality], and hurricane wind is the largest forest
  risk in the southeast [jof-hurricane-wind-risk-forests]. Loss is not gradual — a single
  storm converts standing sawtimber into salvage, and salvage floods the local market and
  depresses the price at the moment the owner is forced to sell.
- **Wildfire.** Real, but for a *managed* plantation it is mediated by suppression, firebreak
  and prescribed-burn practice, none of which our layer models. Treat our wildfire number as
  a landscape hazard signal that a stand sits in a burning region, never as this stand's
  probability of burning.
- **Drought.** Two distinct channels and they operate on different timescales: growth loss
  within a rotation, and fire preconditioning. Our drought layer measures departure from a
  historical reference, so it speaks to the second more directly than the first.
- **Productivity (`conifer-npp`).** Not a hazard. It is the yield side, and it is the only
  layer that speaks to the asset's return rather than its risk of loss.

Rotation length is the clock everything is measured against: roughly 15–20 years to
pulpwood, 25–35 years to sawtimber [usda-srs-gtr075-pine-silviculture]. A hazard that
matters in the 2090s and not in the 2040s is a hazard for the *next* owner. Ask which it is.

## Metrics that lead

- Lead with **cyclone exposure and its trend** where the site is within reach of the coast.
  Everything else is secondary for this asset class.
- Lead with **productivity level and trend** where the question is yield rather than loss.
- Report **percentile for the level** and the recommended **slope for the direction**, and
  never merge them: a stand at the 90th percentile that is flat is a known operating
  condition, and one climbing steeply from the 40th is a changing one.
- The Climate Score is a screening aid for ranking many sites. For a single stand it hides
  more than it shows, because it averages the productivity axis into the hazard axes.

## Known confounders

- **CO₂ fertilisation is inside the productivity number.** Rising NPP in these projections is
  partly a modelled response to elevated CO₂, not purely a climate response. Its real-world
  magnitude is contested and depends on nutrient and water limitation, which these models
  treat differently from one another. Never present rising NPP as a straightforward benefit.
- **The productivity layer is per unit stand area, behind a 2% cover mask.** A site with no
  modelled conifer stand returns "off layer mask", not zero. That is a statement about the
  model's vegetation map, not about whether trees grow there.
- **Wind is not fully covered.** The cyclone layer measures tropical cyclone wind exposure
  and stops there: no storm surge, no rainfall flooding, no mid-latitude windstorm, and no
  high-forcing scenario. For the hazard that matters most to this asset class, we have the
  thinnest scenario coverage in the assessment. Say so plainly.
- **Management dominates at the site scale.** Species and provenance choice, stocking,
  thinning regime and rotation length move outcomes more than the between-scenario spread
  over a single rotation. A report that implies climate is the controlling variable will be
  dismissed by anyone who manages stands, and rightly.

## Framing and vocabulary

Use the owner's units: rotation, stocking, stand, basal area, salvage, sawtimber and
pulpwood. Say "growth" rather than "net primary productivity" once the technical section is
past.

Avoid "risk score" as a noun for the percentile — foresters will read it as an actuarial
probability of loss and it is not one. "More exposed than 80% of global land" is longer and
survives challenge.

## Questions for the customer

- Rotation length and current stand age, per tract. This decides which horizon matters.
- Product mix — pulpwood, sawtimber, both, and any carbon or conservation-easement income.
- Insurance and salvage arrangements, and whether a windthrow event is insured or absorbed.
- Species and provenance actually planted, and whether any deployment is already
  climate-informed.
- Whether the tract is held to rotation or to a fund term, and when that term ends.

## Do not claim

- Do not state or imply an expected financial loss, an annual loss probability, or a return
  impact. We have hazard exposure only.
- Do not describe wildfire exposure as this stand's fire risk — suppression and fuel
  management are absent from the model.
- Do not present rising productivity as a net benefit without the CO₂ caveat.
- Do not imply the wind assessment is complete. It covers tropical cyclone wind only, and
  only under low and medium forcing.
