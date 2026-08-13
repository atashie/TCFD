---
facet: region
id: central-europe
name: Central Europe
matches:
  countries: [Germany, Austria, Switzerland, Czechia, Czech Republic, Poland, Slovakia]
  states: []
confirmed_on: null
sources:
  - id: sciencedirect-drought-bark-beetle-central-europe
    title: "Drought initialised bark beetle outbreak in Central Europe: meteorological factors and infestation dynamics"
    publisher: "Forest Ecology and Management"
    url: "https://www.sciencedirect.com/science/article/abs/pii/S0378112723009003"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2023"
    what_it_supports: >-
      The 2018 drought initiated the most severe spruce bark beetle outbreak in Central
      Europe, and drought is the inciting factor rather than a coincident one.
  - id: pmc-drought-bark-beetle-transformation
    title: "The increasing role of drought as an inciting factor of bark beetle outbreaks can cause large-scale transformation of Central European forests"
    publisher: "Global Change Biology (PMC)"
    url: "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12098194/"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2025"
    what_it_supports: >-
      Drought-driven bark beetle outbreaks are capable of transforming Central European
      forest composition at landscape scale, not only reducing a single rotation's yield.
  - id: sciencedirect-beech-spruce-mortality-2018-19
    title: "Tree mortality of European beech and Norway spruce induced by 2018–2019 hot droughts in central Germany"
    publisher: "Agricultural and Forest Meteorology"
    url: "https://www.sciencedirect.com/science/article/abs/pii/S0168192321001659"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2021"
    what_it_supports: >-
      That the 2018-2019 hot droughts induced tree mortality of Norway spruce in central
      Germany. Title and subject only -- the full text is paywalled, so no site-level
      mortality percentage is taken from it.
---

## Transmission channels

**Drought is the hazard that matters here, and it does its damage indirectly.** The dominant
mechanism for Norway spruce is drought → physiological stress → loss of resin defence →
*Ips typographus* bark beetle outbreak → stand mortality
[sciencedirect-drought-bark-beetle-central-europe]. The 2018 drought initiated the most
severe outbreak on record; the 2018–2019 hot droughts induced severe Norway spruce mortality
in central Germany [sciencedirect-beech-spruce-mortality-2018-19], and the effect is capable
of transforming
forest composition at landscape scale rather than merely reducing one rotation's yield
[pmc-drought-bark-beetle-transformation].

This matters for how our drought layer is read. It measures departure from a fixed
historical reference — which is *exactly* the quantity the beetle mechanism responds to,
because spruce defence fails relative to the conditions the stand established under, not
relative to a global aridity scale. **Of all our layers and all our regions, this is the
pairing where the hazard metric most directly matches the documented damage pathway.**

**Wildfire is a weak signal here** compared with North America. Central European fire regimes
are largely ignition-limited rather than fuel- or climate-limited, and our layer models
neither ignition nor suppression. Do not lead with it.

**Tropical cyclone does not apply.** The layer's values are effectively nil at these
latitudes, and reporting it invites the reader to conclude the assessment is generic.

## Known confounders

- **The beetle is not in the model.** Our drought layer captures the *inciting* condition,
  not the outbreak. Actual mortality depends on stand age, density, thinning history,
  sanitation felling, and the beetle population already present. A moderate drought signal
  in a dense mature spruce monoculture is a far worse prospect than a severe one in a mixed
  stand.
- **Species conversion is already under way.** Much of the German spruce estate is being
  converted to mixed and broadleaf stands after 2018–2020. A projection for *spruce*
  productivity may be describing a species the owner is actively replacing — ask before
  reporting a productivity trend as the asset's outlook.
- **`conifer-npp` covers temperate needleleaf evergreen stands**, which includes Norway
  spruce, but it is a productivity model with no pest module at all. Rising modelled
  productivity is fully compatible with catastrophic beetle mortality.
- **Timber markets here are salvage-sensitive.** Large mortality events flood the regional
  market with damaged wood and depress prices, so loss reaches the owner through price as
  well as volume.

## Framing and vocabulary

Hectares, not acres. Cubic metres, not board feet. "Calamity wood" / *Kalamitätsholz* is the
term for storm- and beetle-damaged salvage volume and it will be understood immediately.
Reference the 2018–2020 period directly — every forest owner in this region dates their
thinking from it, and a report that does not mention it reads as though it were written
about somewhere else.

## Questions for the customer

- Species composition and age-class distribution per stand, and any conversion programme
  already committed.
- Losses during 2018–2020, and how salvage was handled.
- Current sanitation-felling and monitoring practice.
- Whether the holding is certified (FSC/PEFC) in a way that constrains response options.

## Do not claim

- Do not present the drought result as a forecast of beetle mortality. It is the inciting
  condition; the outbreak depends on stand and management factors we do not model.
- Do not report tropical cyclone exposure for these sites as a finding.
- Do not present rising modelled spruce productivity as a favourable outlook without the
  pest caveat — the two are entirely compatible.
- Do not import North American fire framing into this region.
