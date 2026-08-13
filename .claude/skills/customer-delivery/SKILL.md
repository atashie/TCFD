---
name: customer-delivery
description: End-to-end customer climate-risk delivery — assemble the location/asset inputs, generate the deterministic CSV extract and its QA dashboard, then the reports and caveat documentation. Use whenever a client needs their sites assessed, or when the asset→layer catalog needs updating.
---

# Customer Delivery

**This skill is the entry point for the whole delivery.** Everything a delivery needs is
either here or pointed at from here. Supplementary detail lives in the files below; do not
restate them, read them.

| Read | For |
|---|---|
| [ASSET-CATALOG.md](../../../ASSET-CATALOG.md) | The reference doc — schemas, the Climate Score, the dashboard, spatial averaging, registry maintenance |
| [OUTPUT-SPEC.md](../../../OUTPUT-SPEC.md) | What the underlying layers guarantee |
| [GUARDRAILS.md](../../../GUARDRAILS.md) | §3 dynamic scenario discovery, §12 reference sites |
| `config/asset_catalog.yaml` | asset type → hazard layers |
| `config/layer_registry.yaml` | layer → disk location, status, which slope to read |

## The four stages

A delivery is **not** a single command. It is four stages, and each one's status is stamped
into `manifest.json` so a half-finished delivery is visible in its own folder rather than in
someone's memory of what they ran.

| Stage | Produces | Built? |
|---|---|---|
| **1. Inputs** | a confirmed location/asset list | yes — manual, guided below |
| **2. Extract** | the CSV star schema **and** `dashboard.html` | yes |
| **3. Reports** | compliance report + bespoke customer report | **not built** |
| **4. Caveats** | per-delivery anomaly and caveat documentation | **not built** |

Check where a delivery stands:

```bash
python -c "import json;print(json.load(open('deliveries/<c>/<d>/manifest.json'))['stages'])"
```

---

## Stage 1 — Assemble and confirm the inputs

The customer provides sites; you turn them into a valid input file. **One row per
location-asset combination** — a site with three assets is three rows.

| Column | Required | Notes |
|---|---|---|
| `Location` | yes | Site name. Same name + coordinates ⇒ one `location_id`. |
| `Lat`, `Lon` | yes | Decimal degrees. |
| `Asset_Type` | yes | Matched case-insensitively against catalog names and aliases. |
| `Sub_Asset_Unit` | no | Species, facility class, whatever the customer calls it. |
| `Country`, `State`, `City`, `Region`, `Subregion` | no | Carried through to `locations.csv`. |
| `Layers` | no | `;`-separated layer_ids overriding the catalog for that row. |
| `Coord_Source` | no | `supplied` (default) or how you derived it. See below. |

### When information is missing, produce it — then prove it

You are expected to fill gaps rather than bounce the request back. What is **not** allowed
is filling them silently.

**Missing coordinates.** Derive them from the place name, then do all three of:

1. Set `Coord_Source` to how you got it (`derived: city centroid`, `derived: customer
   address`, …). It lands in `locations.csv`, so a derived coordinate is never
   indistinguishable from a surveyed one in the delivered file.
2. Cross-check against any `Country`/`State`/`City` the customer gave. A coordinate that
   disagrees with the stated country is wrong, not approximately right.
3. After the run, check `data_status`. `OUTSIDE_DOMAIN` on a site the customer believes is
   on land means the coordinate is wrong or the pair is reversed.

Remember the footprint: extraction is a ~1° × 1° blend, and moving a site 0.25° changed one
measured burnt-area value by 166%. A city-centroid coordinate for a large estate is a real
approximation — say so.

**Missing or vague asset type.** Ask. Do not guess. "Facility" is not an asset type; what
matters is which hazards have a transmission channel to it, which is a question about what
the site *does*. An asset type absent from the catalog is a hard error by design — work it
out with the user and add it (Stage 1 also owns catalog maintenance, below).

**Sites that may be offshore.** Customer locations are expected to be masked to land
upstream. Extraction will still return a value for a coastal site by drawing on its land
neighbours, and will not warn you — that is accepted behaviour, not a bug. If upstream
masking is uncertain, check `data_status` and the coastal caveat in ASSET-CATALOG.md.

### Catalog maintenance is part of this stage

- An unknown asset type → work out which hazards genuinely reach it, add the entry to
  `config/asset_catalog.yaml` with a `rationale`, set `confirmed_on` to today.
- A bundle that proved wrong in use → fix `layers:` and date it.
- A user approving an existing entry → set its `confirmed_on`.

`confirmed_on: null` means **nobody has ever reviewed that entry**. Everything seeded
2026-08-12 is still null. Filling the date in to quiet the warning, without the user
actually approving it, is the one failure this field exists to prevent.

### Then show the plan and get agreement

```bash
python scripts/generate_customer_delivery.py --customer "<name>" --input <sites.csv>
```

Planning is the default and touches no data. **Show the resolved asset → layer mapping to
the user and get agreement before adding `--run`.** The mapping encodes a claim about which
hazards reach which asset, and the customer is the one who knows their business.

Surface the per-layer `NOTE:` lines from the plan in that conversation — they carry the
things that change how a number should be read (`drought-3b` measures departure from a fixed
reference rather than aridity; `cyclone` has a single impact model and so no structural
uncertainty in its CI; `conifer-npp` covers temperate needleleaf evergreen stands only).

---

## Stage 2 — Extract and build the dashboard

```bash
python scripts/generate_customer_delivery.py --customer "<name>" --input <sites.csv> --run
```

**One command produces both the CSVs and the dashboard.** A delivery shipped as bare CSVs
gives nobody a way to look at what was extracted, and the reference-site check below then
depends on someone remembering a second command.

If the dashboard is skipped (`--no-dashboard`, for debugging) or its build raises, the run
**exits non-zero** and writes `DELIVERY-INCOMPLETE.md` into the folder. The verifier refuses
any delivery carrying that marker. The CSVs may well be correct; the delivery is not
finished, and the folder says so rather than looking shippable.

Output in `deliveries/{customer-slug}/{YYYYMMDD}/` — schemas are in ASSET-CATALOG.md.

### Verify

```bash
python scripts/test_customer_delivery.py deliveries/{customer}/{date}
```

Recomputes every metric from the source NetCDF with an implementation independent of the one
that wrote it, then checks referential integrity, source hashes, percentile orientation, the
`slopes_agree` rule, NaN baseline slopes, `data_status`, the Climate Score (including the
balanced-panel baseline invariant), the dashboard payload, and the stage record. Exits
non-zero on any violation. A missing dashboard is a failure, not a note.

### Then LOOK at it

```bash
open deliveries/{customer}/{date}/dashboard.html
```

The verifier checks payloads, not layout. **Open the page.** If you have not, report the
dashboard as *unreviewed*, exactly as the layer workflow requires for maps.

If it looks stale, check the build stamp before debugging anything:

```bash
python scripts/generate_delivery_dashboard.py deliveries/{customer}/{date} --check-stamp
```

A cached page is indistinguishable from a fresh one by eye and has produced two phantom bug
reports. The header shows the stamp; a mismatch means hard-reload.

### Reference-site sanity, per GUARDRAILS §12

A schema-valid, passthrough-verified delivery can still be meaningless — the verifier proves
the numbers were carried faithfully, not that the layer means what its name says. Check two
or three sites where you know the answer:

- a coastal Gulf site should carry non-zero cyclone exposure; an inland western site should
  read exactly 0 (percentile 1)
- a northern-California site should outrank a central-European one on wildfire
- on `conifer-npp`, high productivity should read as *low* risk

That last one is an orientation check, not proof of passthrough — value and percentile are
independently spatially averaged, so a wrong transformation can preserve the ordering. The
verifier is what proves it.

Empty dashboard panels annotate themselves with the reason (`OUTSIDE_DOMAIN`, etc.), so an
empty panel is information, not a bug.

---

## Stage 3 — Reports *(not built)*

Two distinct outputs, not one:

- **Physical hazard compliance report** — TCFD/CDP disclosure framing, standardised
  structure, auditable against the CSV.
- **Bespoke customer report** — narrative for this customer's portfolio and decisions.

Both are downstream of Stage 2 and must read `values.csv` / `climate_score.csv` rather than
recomputing anything. When built, they record themselves via
`record_stage(out_dir, "compliance_report" | "bespoke_report", "built")`.

Decisions already made that constrain them: RCP↔SSP harmonization belongs here, not in the
CSV; the Climate Score is the portfolio-level headline; percentile is the only cross-hazard
comparable axis.

---

## Stage 4 — Caveats and anomalies *(not built)*

Per-delivery documentation of what is odd about *this* customer's data. Most of the raw
material already exists and is currently scattered — the intent is to collect it:

- layers with `qa_reviewed_on: null` (all of them today)
- `OFF_LAYER_MASK` / `OUTSIDE_DOMAIN` sites, and derived coordinates
- assets dropped from a balanced panel because a layer lacks a forcing tier
- per-layer `delivery_note` and `interpretation_caveat`
- incomplete `n_hazards` on any Climate Score row
- the coastal-renormalization and ~1° footprint caveats

Records itself via `record_stage(out_dir, "caveats", "built")`.

---

## After reprocessing any layer

```bash
python scripts/generate_customer_delivery.py --measure-slopes
```

Flags any layer whose measured sen-zero share disagrees with the registry's
`recommended_slope`. The correct slope is a property of the data, not the hazard —
`wildfire` proved a layer's zero fraction does not predict it.

## Scope

Point extraction only; polygon/region carry a documented `cos(lat)` defect and are not wired
in. Stage 2 produces numbers and provenance only — scoring, ranking and narrative belong to
Stage 3.
