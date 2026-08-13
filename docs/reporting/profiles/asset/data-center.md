---
facet: asset
id: data-center
name: Data centre
aliases: [datacenter, data centre, colocation facility, server farm, hyperscale]
catalog_entry: data center
confirmed_on: null
sources: []
---

## Transmission channels

**The honest headline for this asset class is that our coverage is poorly matched to it**,
and a report that does not open with that is misleading by omission.

The hazards that actually reach a data centre are, in rough order:

- **Extreme heat** — cooling load, equipment derating, and the design temperature the
  facility was built to. **Not covered by this assessment.**
- **Water availability** — for evaporatively-cooled facilities, and for the thermal plant on
  the local grid. The operative question is supply against competing demand. **Our drought
  layer is not that metric**; it measures departure from a fixed historical reference.
- **Grid reliability** under peak load. Not a climate layer, and not covered.
- **Flooding** — the classic total-loss mechanism for a facility with plant at or below
  grade. **Not covered.**
- **Wind** — real for coastal siting, but rarely the binding constraint for a hardened
  facility.

The asset almost never fails structurally. It fails to **operate**, and the loss is
contractual: uptime commitments, service credits, and customers with alternatives.

## Metrics that lead

- Lead with what is **not** assessed. For this asset class that is the finding.
- Where water is a genuine process input, report the drought result with an explicit
  statement of what it does and does not measure.
- Wind and fire exposure where the region profile makes them material — and say plainly when
  it does not.

## Known confounders

- **Drought is not water stress.** Reporting one as the other is the single most likely
  error for this asset class, because the words are close and the metrics are not.
- **Design temperature is the threshold that matters**, and it is a facility property we do
  not hold. A modest shift in extreme heat can cross it while a large shift in mean
  conditions does not.
- **Urban and peri-urban siting** means the land cover in the model bears little relation to
  the site.
- **Redundancy and tier rating** dominate outcomes and are invisible here.

## Framing and vocabulary

Facility, campus, capacity, load, uptime, availability, interruption, PUE, WUE. "Availability"
and "interruption", not "damage".

## Questions for the customer

- Cooling technology and water dependence — this decides whether water is a process input.
- Design temperature assumptions and the behaviour when exceeded.
- Grid arrangements, on-site generation, and contracted uptime commitments.
- Flood zone status and the elevation of critical plant.
- Any existing heat or water study we should reconcile against rather than duplicate.

## Do not claim

- Do not present the drought result as a water-availability or water-stress assessment.
- Do not imply the assessment addresses heat, cooling load, grid reliability or flooding.
- Do not describe a data centre as assessed for its material hazards on current coverage.
- Do not report a favourable result for this asset class without stating that the hazards
  most likely to affect it were not examined.
