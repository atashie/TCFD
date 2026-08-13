---
facet: region
id: us-southeast
name: US Southeast
covers: [North Carolina, South Carolina, Georgia, Florida, Alabama, Mississippi, Louisiana, Arkansas, Tennessee, Virginia, East Texas]
confirmed_on: null
sources:
  - id: jof-hurricane-wind-risk-forests
    title: "Potential Hurricane Wind Risk to US Rural and Urban Forests"
    publisher: "Journal of Forestry 119(4)"
    url: "https://academic.oup.com/jof/article/119/4/393/6260883"
    source_type: peer-reviewed
    retrieved_on: "2026-08-13"
    as_of: "2021"
    what_it_supports: >-
      Hurricane wind risk to forests is concentrated in the southeastern United States.
  - id: srs-hurricane-timber-markets-2025
    title: "Evaluating hurricane impacts on timber markets in the Southeastern United States"
    publisher: USDA Forest Service Southern Research Station
    url: "https://www.srs.fs.usda.gov/pubs/ja/2025/ja_2025_brandeis_002.pdf"
    source_type: regulator
    retrieved_on: "2026-08-13"
    as_of: "2025"
    what_it_supports: >-
      Hurricane damage propagates into regional timber markets, not only into the damaged
      stand — the salvage and price channel.
---

## Transmission channels

The regional signature is **compound coastal events**: wind, rainfall flooding and storm
surge arriving together. Our assessment resolves the first of those three and neither of the
others, which is the single most important thing to say honestly about any coastal site here
[jof-hurricane-wind-risk-forests].

The second regional signature is **market propagation**. A storm large enough to damage one
holding damages the region's holdings, and the salvage volume that follows moves the local
price [srs-hurricane-timber-markets-2025]. A site-by-site exposure report understates
correlated regional loss by construction, because it treats each site independently.

Inland sites in this region have a genuinely different profile from coastal ones: much
lower cyclone exposure, and drought and fire doing proportionately more of the work.

## Known confounders

- **Our grid is coarse relative to this coastline.** A ~1° footprint spans the distance over
  which cyclone wind exposure decays inland. A site tens of kilometres from the coast is
  blended with cells that are not, in both directions.
- **The `cyclone` layer is additionally smoothed** at processing time because raw storm
  tracks are one cell wide. Values in this region are the most spatially smoothed numbers in
  the whole assessment, and should never be described as site-specific.
- **No high-forcing cyclone scenario exists.** In the region where tropical cyclone is the
  leading hazard, the high tier is the one tier that cannot include it. Any high-tier
  statement about a southeastern portfolio is a statement about fire and drought only.
- **Sea level rise and coastal flooding are absent entirely.** For low-lying coastal tracts
  this is the dominant omission, and it is invisible in the results.

## Framing and vocabulary

"Hurricane" in conversation, "tropical cyclone" where it must match the data label. Name the
states; "the Southeast" is too coarse for an audience that operates there. Distinguish
coastal plain from piedmont — readers here do, and a report that does not sounds remote.

## Questions for the customer

- Distance from the coast, per tract, and elevation for anything low-lying.
- Recent storm history on the holding and how salvage was handled.
- Whether nearby mill capacity constrains salvage — this determines whether a wind event is
  a loss of volume or a loss of value.
- Whether any tract sits in a floodplain or is served by a single access road.

## Do not claim

- Do not describe a coastal site as fully assessed for storm risk. Wind only.
- Do not compare a coastal and an inland site on the Climate Score alone: the coastal site's
  dominant hazard is thinly covered and the inland one's is not, so the comparison is
  between differently complete assessments.
- Do not attribute a regional price effect to a site-level exposure result.
