---
facet: region
id: us-mid-atlantic
name: US Mid-Atlantic
matches:
  countries: [USA, United States, US]
  states: [Virginia, Maryland, Delaware, Pennsylvania, New Jersey, New York, West Virginia, District of Columbia]
confirmed_on: null
sources: []
---

## Transmission channels

**This is the region where our hazard set fits the assets worst, and saying so is the most
useful thing a report can do here.**

The mid-Atlantic carries a heavy concentration of data centres, logistics and commercial
property. The hazards that actually reach those assets are:

- **extreme heat** — cooling load, equipment derating, workforce limits. **Not covered.**
- **water availability** — evaporative cooling and the local grid. Our drought layer is the
  nearest thing we have, but it measures departure from a historical reference rather than
  supply against demand, which is the operative question. **Partially covered, and not by
  the right metric.**
- **grid reliability** under peak load. **Not covered, and not a climate layer at all.**
- **riverine and pluvial flooding.** **Not covered.**
- **remnant tropical cyclone** — systems reaching this far north deliver rainfall flooding
  more often than damaging wind, and our cyclone layer measures wind. **Poorly matched.**

Wildfire and drought will produce numbers for a mid-Atlantic site. Those numbers are real
and correctly computed, and they are **not the hazards that threaten the asset**. A report
that leads with them because they are what we have is the exact failure the coverage section
exists to prevent.

## Known confounders

- **Water stress is not drought.** Where a facility's water question is supply against
  competing demand, the drought layer does not answer it. The Water Risk Index is a separate
  product with a different output contract and cannot simply be substituted.
- **Grid and utility dependence is the real transmission channel** for a data centre, and it
  is an infrastructure question rather than a climate-layer question.
- **Urban and peri-urban siting** means the surrounding land cover in the model bears little
  relation to the site itself.

## Framing and vocabulary

Facility, campus, capacity, uptime, PUE, cooling load. For data centres, "availability" and
"interruption" are the operative risk words, not "damage" — the asset rarely fails
structurally, it fails to operate.

## Questions for the customer

- Cooling technology and water dependence. This decides whether water is a process input or
  a background condition.
- Grid supply arrangements, on-site generation, and contracted uptime commitments.
- Flood zone status and stormwater arrangements.
- Design temperature assumptions and what happens when they are exceeded.
- Whether an existing heat or water study exists that we should reconcile against rather
  than duplicate.

## Do not claim

- Do not lead with wildfire or drought for a mid-Atlantic built asset simply because they
  are the layers available.
- Do not present drought exposure as a water-availability assessment.
- Do not imply the assessment addresses heat, cooling load or grid reliability.
- Do not describe a mid-Atlantic facility as assessed for its material hazards. On current
  coverage it is not, and the report should say which are missing.
