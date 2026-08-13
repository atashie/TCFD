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
| [ASSET-CATALOG.md](../../../ASSET-CATALOG.md) | Stages 1–2 — schemas, the Climate Score, the dashboard, spatial averaging, registry maintenance |
| [docs/reporting/README.md](../../../docs/reporting/README.md) | **Stages 3–4** — the map to everything below |
| [docs/reporting/technical/data-boundaries.md](../../../docs/reporting/technical/data-boundaries.md) | before writing any sentence about what a number means |
| [docs/reporting/compliance/framework-map.md](../../../docs/reporting/compliance/framework-map.md) | before claiming a disclosure requirement is met |
| [docs/reporting/research/method.md](../../../docs/reporting/research/method.md) | before looking anything up about a customer |
| [docs/reporting/profiles/README.md](../../../docs/reporting/profiles/README.md) | before writing or editing a facet profile |
| [OUTPUT-SPEC.md](../../../OUTPUT-SPEC.md) | What the underlying layers guarantee |
| [GUARDRAILS.md](../../../GUARDRAILS.md) | §3 dynamic scenario discovery, §12 reference sites |
| `config/asset_catalog.yaml` | asset type → hazard layers |
| `config/layer_registry.yaml` | layer → disk location, status, which slope to read |
| `config/hazard_taxonomy.yaml` | which hazard families exist, which we cover, and what each gap would take |

## The four stages

A delivery is **not** a single command. It is four stages, and each one's status is stamped
into `manifest.json` so a half-finished delivery is visible in its own folder rather than in
someone's memory of what they ran.

| Stage | Produces | Built? |
|---|---|---|
| **1. Inputs** | a confirmed location/asset list | yes — manual, guided below |
| **2. Extract** | the CSV star schema **and** `dashboard.html` | yes |
| **4. Caveats** | `caveats.json` + `caveats.md` | yes |
| **3a. Compliance report** | `report_compliance.html` — IFRS S2 spine | yes |
| **3b. Bespoke report** | `report_bespoke.html` — composed for this reader | yes |

**Stage 4 runs before Stage 3, and the row order above is not a typo.** The caveat set is a
mechanical derivation from the manifest, the CSVs and `config/hazard_taxonomy.yaml`, and it
is an *input* to both reports rather than a summary of them: each report is required to carry
every `must_disclose` entry, the builder refuses without them, and the verifier re-checks
afterwards. Generating them last would let two documents about one delivery disagree about
what is wrong with it.

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

Remember the footprint: extraction is a ~1° × 1° blend, and shifting sites by 0.25° changed
measured burnt-area values by 44–569% on one example portfolio
(`python scripts/measure_extraction_sensitivity.py` reproduces it). A city-centroid
coordinate for a large estate is a real approximation — say so.

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

## Stage 4 — Caveats

```bash
python scripts/generate_delivery_caveats.py deliveries/{customer}/{date}
```

Derives this delivery's caveat set mechanically and writes `caveats.json` + `caveats.md`.
Nothing here is authored: it comes from the manifest, `layers.csv`, `values.csv`,
`climate_score.csv`, `report_config.yaml` and `config/hazard_taxonomy.yaml`.

Ids are **semantic and stable** (`HAZARD-COVERAGE`, not `CAV-001`) because narratives cite
them, and severity decides enforcement:

| Severity | Meaning |
|---|---|
| `must_disclose` | must appear in **every** report built from this delivery; the builder refuses and the verifier re-checks |
| `should_note` | changes how one figure reads, rather than what the delivery can claim |
| `informational` | context for whoever maintains the delivery |

It does **not** judge materiality. Severity says how certainly something must be disclosed,
never how much it matters to this customer — that needs business context and belongs in the
bespoke narrative, with a citation.

---

## Stage 3a — Compliance report

```bash
python scripts/generate_compliance_report.py deliveries/{customer}/{date}
```

IFRS S2 spine with a CDP / ESRS E1-9 / California SB 261 mapping appendix.
**Fully deterministic** — no narrative, no researched claims, so an auditor can rebuild it
byte-for-byte. Structure and rationale: `docs/reporting/technical/report-anatomy.md`.

**Every regulatory claim needs its paragraph checked, and every measurement needs a
reproduction.** An external review found three wrong IFRS S2/ESRS citations, an invented
"standard taxonomy" provenance, and two quoted measurements with no retained receipt. Before
asserting that a section satisfies a requirement, read the requirement; before quoting a
number as measured, make sure something in the repo re-measures it.

**A section may report nothing.** Where the method for turning our data into a framework's
requirement has not been agreed with the user, the section states that and publishes no
figure. Entries live in `report_common.TBD_SECTIONS`; both reports render the same block, and
the verifier **fails** any report that publishes a figure whose method is still deferred.
Currently deferred: the IFRS S2 29(c) vulnerable-asset count. Err toward an explicit gap over
a fast answer — a gap is a conversation, a wrong number is a liability with the customer's
name on it.

Two things to check before it goes anywhere:

- **`report_config.yaml`.** While `horizons.source` or `vulnerability.source` reads
  `default`, the report discloses that we chose them. IFRS S2 10(d) wants the *entity's*
  horizons and the vulnerability threshold is a risk-appetite decision — get both from the
  customer and set the value and the source together.
- **Asset values.** Without `Asset_Value` in the input, 29(c)'s monetary half cannot be
  reported and a `must_disclose` caveat says so.

---

## Stage 3b — Bespoke report

**Facets accept lists, and `region` is determined by the data.** Each region profile declares
a `matches:` block; `assert_region_coverage()` refuses to build if any delivered location is
uncovered. A portfolio spanning California, the Gulf coast, the mid-Atlantic and Bavaria needs
four region profiles, not one — describing it with a single region's framing is confidently
wrong about the first thing a reader checks.

```bash
# 1. choose facets in report_config.yaml, then
python scripts/generate_bespoke_report.py deliveries/{customer}/{date} --scaffold
# 2. write narrative.md and dossier.yaml, then
python scripts/generate_bespoke_report.py deliveries/{customer}/{date}
```

Not chained onto `--run`: it needs facet profiles chosen and a narrative written.

**The generator owns the numbers; a person owns the narrative.** The boundary is enforced —
an unfilled slot, a `TODO`, an uncited paragraph in a cited slot, or a citation that does not
resolve all fail the build. Read `docs/reporting/research/method.md` before researching
anything and `docs/reporting/profiles/README.md` before touching a profile.

**Profiles guide the narrative and are never pasted into it.** The scaffold surfaces them as
comments; comments are stripped before rendering and the verifier fails a report containing
any. Pasting them would make every report of a given asset class identical while looking
bespoke.

---

## Verify, then look

```bash
python scripts/test_customer_delivery.py deliveries/{customer}/{date}
open deliveries/{customer}/{date}/report_compliance.html
```

The verifier recomputes the vulnerable-asset count independently and compares it to the
number the report printed, checks every `must_disclose` caveat appears, re-validates the
narrative and its citations, and refuses a report leaking internal vocabulary or an HTML
comment.

It does not check layout. **No browser can be driven on this machine**, so a report that has
not been opened is *unreviewed* and must be reported as such — the same rule the layer
workflow applies to maps. PDF is Safari's File ▸ Print ▸ Save as PDF; pagination has never
been verified here.

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
