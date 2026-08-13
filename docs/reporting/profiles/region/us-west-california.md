---
facet: region
id: us-west-california
name: California and the US West
matches:
  countries: [USA, United States, US]
  states: [California, Oregon, Washington, Nevada, Idaho, Montana, Arizona, Utah, Colorado]
confirmed_on: null
sources:
  - id: usgs-sierra-catastrophic-wildfire
    title: "The Risks and Consequences of Catastrophic Wildfires on the Conifer Forests of the Sierra Nevada"
    publisher: "U.S. Geological Survey, Climate Adaptation Science Centers"
    url: "https://www.usgs.gov/programs/climate-adaptation-science-centers/science/risks-and-consequences-catastrophic-wildfires"
    source_type: regulator
    retrieved_on: "2026-08-13"
    as_of: "2024"
    what_it_supports: >-
      Sierra Nevada conifer forests have become increasingly vulnerable to wildfire through
      the combination of fire suppression, elevated tree mortality and warming.
  - id: stephens-drought-mortality-wildfire
    title: "Drought, Tree Mortality, and Wildfire in Forests Adapted to Frequent Fire"
    publisher: "BioScience / USDA Forest Service"
    url: "https://research.fs.usda.gov/treesearch/55621"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2018"
    what_it_supports: >-
      Acute drought compounded by a century of fire exclusion drove mass conifer mortality in
      the Sierra Nevada, and the resulting dead fuel raises the potential for very large,
      severe fires.
---

## Transmission channels

**Drought and wildfire are one coupled mechanism here, not two hazards.** Acute drought
compounded by a century of fire exclusion produced mass conifer mortality — the Forest
Service estimated 29 million trees lost in California in 2015 alone — and the standing dead
fuel that follows raises the potential for very large, severe fire
[stephens-drought-mortality-wildfire] [usgs-sierra-catastrophic-wildfire]. A drought result
and a wildfire result at the same western site are therefore **not independent**, and adding
or averaging them double-counts a single causal chain.

This is the region where our two layers are most likely to be misread as corroborating each
other. They are closer to being two measurements of the same process.

**Fire suppression history is the dominant variable and is absent from the model.** Our
wildfire layer is a climate-and-vegetation burnt-area model. It does not know that a stand
has missed a century of low-intensity fire, which is the single strongest predictor of
catastrophic loss here. Treat the number as a landscape signal and say so.

**Tropical cyclone does not apply.** Values are effectively nil; reporting them as a finding
invites the reader to distrust the rest.

## Known confounders

- **Fuel load and treatment history dominate.** Two adjacent stands with identical climate
  exposure can differ by an order of magnitude in outcome depending on thinning and
  prescribed-burn history.
- **Elevation and aspect matter at a scale we cannot resolve.** A ~1° footprint spans
  several thousand metres of relief in the Sierra Nevada and the northern Coast Ranges. Site
  values here are more of an approximation than in flat terrain.
- **Mortality events are episodic, not gradual.** A decadal mean smooths over the two- or
  three-year die-offs that actually move a holding's value.
- **Regulatory and insurance conditions move faster than the climate signal.** Availability
  and cost of fire insurance, and harvest-permit conditions, may bind before any projected
  change does.

## Framing and vocabulary

Acres and board feet. Name the region precisely — Cascades, Klamath, northern Sierra, Coast
Range — rather than "California", which spans very different fire and drought regimes.
"Fuel treatment", "defensible space", "salvage" and "mortality event" are the operating
vocabulary.

## Questions for the customer

- Fire history and fuel-treatment record for each holding.
- Insurance status, and whether cover has been repriced or withdrawn.
- Elevation band and aspect, which matter more here than the coordinate alone conveys.
- Access and egress — a single-road holding has a different loss profile.
- Any losses in the 2012–2016 drought or subsequent fire seasons.

## Do not claim

- Do not treat the drought and wildfire results as independent corroboration. They are
  substantially the same mechanism.
- Do not present the wildfire percentile as this stand's probability of burning; suppression
  and fuel management are not in the model.
- Do not report tropical cyclone exposure for these sites as a finding.
- Do not describe a value as site-specific in mountainous terrain.
