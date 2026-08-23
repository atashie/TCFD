# Ingest licensing audit — commercial use and the always-serve-free list

**Date:** 2026-08-22. **Scope:** every data source ingested for the TCFD/CDP product —
all 43 registry layers plus the evaluated-but-not-shipped ingests. Water Risk Index is out
of scope (out of storefront v1; its sources are ISIMIP 2b/3b only, so the round rules
below cover it anyway). This is build item 1 of the storefront plan
([docs/storefront/DECISIONS.md](storefront/DECISIONS.md) §3: licensing audit precedes
anything publicly served).

**Framing (user decision 2026-08-22):** a non-commercial restriction is not a blocker —
an NC dataset can be served in the free lane alongside the for-sale data. The audit's job
is therefore to sort every source into *commercially clean* vs *must always serve free*,
and to collect the attribution obligations that ride with serving at all.

Every licence below was re-verified against the live publisher record on 2026-08-22
unless noted. Receipts inline.

## The answer

**No ingested source forbids us from serving it. One layer has a genuinely disputed
licence and is the only candidate for the always-serve-free list: `landslide-arup`.**
Everything else — every ISIMIP layer, the tornado ladder, `hail-vlh`, GEBCO terrain,
Natural Earth boundaries — is commercially clean under CC0, CC BY 4.0, or US-government
public domain.

| Source | Feeds | Licence (verified 2026-08-22) | Commercial use | Free-lane flag |
|---|---|---|---|---|
| ISIMIP3b (all: OutputData, DerivedOutputData, bias-adjusted forcing) | 34 of 39 hazard layers | CC0 1.0, per-dataset | Yes, unrestricted | — |
| ISIMIP2b DerivedOutputData Lange2020 | `drought-2b`, `cyclone`, `heatwave-2b` | CC0 1.0, per-dataset | Yes, unrestricted | — |
| ISIMIP2b InputData sealevelrise | `sealevel-2b` | CC BY 4.0 | Yes, with attribution | — |
| ISIMIP2b OutputData biomes | `conifer-npp` | CC BY 4.0 | Yes, with attribution | — |
| NOAA SPC severe weather database | `tornado-*` (4 rungs) | US Gov public domain | Yes | — |
| Natural Earth admin-0 (CONUS mask) | `tornado-*` | Public domain | Yes | — |
| GEBCO_2026 grid (SRTM15+ land base) | `sealevel-2b` | Public domain, citation requested | Yes | — |
| Battaglioli 2026 article Source Data | `hail-vlh` | CC BY 4.0 (article route — see guard) | Yes, with attribution | — |
| **World Bank/GFDRR–Arup landslide map** | **`landslide-arup`** | **Disputed at publisher: CC BY-NC 4.0 vs CC BY 4.0** | **Publisher records disagree** | **Recommended — open user call** |
| Hail severity, Nature 653 (figshare) | evaluation only, nothing served | CC BY 4.0 (data and code) | Yes | n/a |

## The one open call: `landslide-arup`

The publisher's own records still disagree, re-verified today:

- `datacatalog.worldbank.org/search/dataset/0037584` — **"Creative Commons
  Attribution-Non Commercial 4.0"**, classification "Public".
- `energydata.info/dataset/global-landslide-hazard-map` (World Bank Group's own CKAN
  mirror) — **"Creative Commons Attribution 4.0"**.
- The 113-page Arup project report — no licence statement at all (checked 2026-08-19).

The standing determination (2026-08-19) cleared the layer for our limited commercial use
with mandatory attribution — made when the alternatives were ship-or-drop. Today's
framing adds a third option that did not exist then: **serve it, but only in the free
lane.** That posture costs nothing (the free lane exists regardless), removes the only
licence ambiguity in the paid roster, and keeps the layer in every delivery. The
alternatives remain: (a) keep the 2026-08-19 determination, leaning on the publisher's
own CC BY mirror record and "Public" classification; (b) ask the World Bank data help
desk to state which record is authoritative — the 2026-08-19 review resolved a licence
question by asking once already, and either answer here is actionable.

Recommendation: free-lane flag until (b) resolves it. Not applied anywhere yet — the
registry has no licence/free-lane field to carry it (see *Wiring*, below).

Whatever the call, **attribution is required under both candidate licences** (both are
BY): "World Bank / GFDRR Global landslide hazard map, produced by Arup (2021)" — already
carried in the file's `attribution_required` attribute and the registry `delivery_note`.

## Non-ISIMIP sources, in detail

### NOAA SPC tornado database → `tornado-f2plus/-all/-f1plus/-f3plus`

US Government work, public domain. Re-verified via the NWS disclaimer
(`weather.gov/disclaimer`): *"The information on National Weather Service (NWS) Web pages
are in the public domain, unless specifically noted otherwise, and may be used without
charge for any lawful purpose."* Two real conditions ride with it: do not imply NOAA/NWS
endorsement, and do not present modified data as official government material — met by
ordinary source attribution ("derived from NOAA SPC tornado reports, 1950–2025"), which
is good practice though not legally required. The CONUS mask input, Natural Earth
admin-0 boundaries, is likewise public domain (*"You may use the maps in any manner…
commercial purposes… no permission is needed"* — naturalearthdata.com terms, verified
today).

### Battaglioli et al. 2026 Source Data → `hail-vlh`

The dual-licence split documented in DATASET-ATTRIBUTES.md and
`reports/maps/hail-vlh/essl_licence_query.md` re-verified today and holds:

- The *Nature Geoscience* article (doi 10.1038/s41561-025-01868-0) is **CC BY 4.0** —
  *"permits use, sharing, adaptation, distribution and reproduction in any medium or
  format"* — and its Source Data files carry the full 0.25° grids. This is the route the
  layer is built on: `process_hail_vlh_battaglioli.py` reads `MOESM2.csv`/`MOESM7.csv`
  (article supplementary files), and the processed file's `source_licence` attribute
  states the Zenodo deposit is not used.
- The Zenodo deposit (10.5281/zenodo.17064885) is **CC BY-NC-ND 4.0**: NC blocks
  commercial use, ND separately blocks the regridding/aggregation any ingest performs.

**Standing guard:** the article's own data-availability statement points readers to
Zenodo — the trap is that "the same data" sits there under the incompatible licence.
Never ingest the Zenodo annual fields, seasonal series, or trend fields into this layer;
that forecloses a 2014–2023 window and a self-computed trend (both documented in the
registry entry). If a recent window, a ≥2 cm threshold, or scenarios ever matter, the
route is the parked ESSL licence query, not the deposit.

Attribution is legally required (CC BY). **Gap:** unlike landslide, the hail registry
`delivery_note` says nothing about attribution and the file has no `attribution_required`
attribute — `source_dataset` (which is exported to layers.csv) carries the citation, but
nothing forces it into a rendered report. See *Wiring*.

Upstream chain note: AR-CHaMo is run on ERA5, and the lightning/loss inputs are
proprietary to Earth Networks, the Met Office, Vaisala, BoM and Munich Re — those
obligations bind the authors' publication, not our reuse of their CC BY output. Our
obligation is the article credit alone.

### GEBCO_2026 terrain → `sealevel-2b`

*"The GEBCO Grid is placed in the public domain and may be used free of charge"* —
commercial use explicitly permitted (gebco.net, verified today). Citation is requested,
not required: "GEBCO Bathymetric Compilation Group 2026 (2026). The GEBCO_2026 Grid.
doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa". Cheap to honor; add it wherever the
coastal layer is served.

### Hail severity dataset (Nature 653, 2026) — evaluation only

Figshare deposit 10.6084/m9.figshare.30103471.v3, **"CC-BY-4.0 (data and code)"** per the
retrieved manifest (`data/raw/hail-nature2026/manifest.json`, 2026-08-18). Not a layer,
nothing served; no obligation until something is. If its findings are ever quoted in a
report, cite the paper (doi 10.1038/s41586-026-10543-2).

## ISIMIP — the baseline under everything else

The repository licences per dataset, and the split is by simulation round
(isimip.org licences page, verified today): **ISIMIP3 is CC0 1.0** (public-domain
dedication — *"anyone [may] use the data for derivative work in any way and for any
purpose, including commercial purposes"*); **ISIMIP2 is CC BY 4.0 preferred, with CC
BY-NC 4.0 and CC BY-SA 4.0 available to modelling groups on request.** Per-dataset
`rights` fields checked today via the repository API for every product we ship:

| Product | Datasets checked | Rights |
|---|---|---|
| 3b OutputData `fire` (wildfire) | listing | CC0 1.0 |
| 3b OutputData `biomes` incl. LPJmL5-7-10-fire thawdepth (permafrost, csoil) | listing | CC0 1.0 |
| 3b InputData bias-adjusted daily forcing (heat/cold/pluvial ladders) | listing | CC0 1.0 |
| 3b DerivedOutputData TipESM2025 CaMa-Flood (`flood-3b-*`) | listing | CC0 1.0 |
| 3b DerivedOutputData Heinicke2026 (`drought-3b`) | listing | CC0 1.0 |
| Zantout2025 (`cropfailure-3b`, `heatwave-3b`) | listing (query surfaced the 3a twins; round rule and every 3b sibling check say CC0) | CC0 1.0 |
| 2b DerivedOutputData Lange2020, KE-TG-meanfield and exposure family (`cyclone`, `drought-2b`, `heatwave-2b`) | listing, n=28 + n=219 | **CC0 1.0** |
| 2b InputData sealevelrise (`sealevel-2b` deltas) | listing, n=72 | **CC BY 4.0** |
| 2b OutputData biomes (`conifer-npp`: CLM45/ORCHIDEE/LPJmL) | listing + licences-page table | **CC BY 4.0** |

So nothing we ship from ISIMIP is NC or SA. Two obligations and one trap remain:

1. **The terms-of-use page carries a blanket "sale of the data is strictly forbidden"
   sentence** alongside the per-dataset licences. Read against CC0/CC BY (which expressly
   permit selling copies and derivatives), the coherent reading is: do not sell the
   repository's files as such. Our paid lane sells derived decadal statistics and
   site-level analysis; our bulk lane redistributes *processed derivatives* free. Both
   are consistent with even the strict reading. **Standing rule: never put raw ISIMIP
   file mirrors behind the paywall** — we have no plan to, and now it's written down.
2. **Attribution.** CC BY datasets require it; for CC0 ISIMIP still *requests* citation
   of dataset DOIs and credit to the modelling groups, plus notification of publications.
   A credits page on the storefront listing ISIMIP, the impact models and GCMs per layer
   satisfies all of it at once, and the bulk-release format already reserves an
   attribution file (DECISIONS.md decision 4).
3. **Trap for future ingests: an ISIMIP2 dataset can individually be CC BY-NC or CC
   BY-SA** ("on request" per model). None of ours is, but any *new* 2b ingest must have
   its `rights` field checked at download time — the field is in the repository API
   response and in the file sidecars the downloader saves. Cheap to check, expensive to
   discover after processing.

## Attribution obligations that ride with serving (free or paid)

| Layer(s) | Credit | Legally required? | Carried today |
|---|---|---|---|
| `landslide-arup` | World Bank / GFDRR Global landslide hazard map, produced by Arup (2021) | Yes (BY under either record) | `attribution_required` attr + registry `delivery_note` |
| `hail-vlh` | Battaglioli et al., *Nature Geoscience* 19, 52–58 (2026), doi 10.1038/s41561-025-01868-0 | Yes (CC BY) | `source_dataset` attr only — **not forced into reports** |
| `sealevel-2b` | ISIMIP2b sea-level input (Nauels et al. 2016; Bamber & Riva 2010); GEBCO_2026 citation (requested) | Yes (CC BY) / requested | `source_dataset` attr only |
| `conifer-npp` | ISIMIP2b biomes: CLM4.5, ORCHIDEE, LPJmL + GCMs | Yes (CC BY) | `source_dataset` attr only |
| All ISIMIP3b layers | ISIMIP + modelling groups + dataset DOIs | Requested (CC0) | `source_dataset` attr |
| `tornado-*` | NOAA SPC tornado reports 1950–2025 | No (PD; no-endorsement rule) | `source_dataset` attr |

## Wiring — what this audit found missing

1. **The registry has no `license` / `attribution` / free-lane fields** (`LayerSpec`
   rejects unknown keys), so the free-lane flag and the credits above live only in prose.
   The storefront build plan already calls for these fields; this audit supplies their
   contents. Adding them touches `LayerSpec`/`load_registry()` and its tests.
2. **`attribution_required` is not on `LAYER_ATTRS_EXPORTED`** (`scripts/utils/delivery.py`)
   — the closed allowlist silently drops it at delivery, so even landslide's credit
   travels as `delivery_note` prose rather than as a structured field a report builder
   can enforce. When the registry fields land, the promotion should work like
   `relative_baseline`: non-empty attribution → rendered credit, verifier fails without it.
3. **`hail-vlh` has no attribution wording anywhere a customer sees** despite CC BY. Until
   the structured wiring lands, the one-line fix is a sentence in its registry
   `delivery_note`, same as landslide's.

None of these blocks serving today's deliveries (reports carry `source_dataset` for every
layer); they block *machine-enforced* compliance, which is the standard the rest of the
caveat system holds.
