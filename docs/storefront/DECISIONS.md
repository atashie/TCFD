# Storefront — decision record

*Created 2026-08-21 to give the 2026-08-20 working-session decisions a durable home in
version control (until now they existed only in session context, and the page below only as
a gitignored file under `reports/`). Reconstructed after the fact — correct or extend as
needed.*

## Settled 2026-08-20

1. **Serving architecture: Cloudflare R2 + Cloudflare Pages + Google Cloud Run.** R2 holds
   the data artifacts, Pages serves the static front end, Cloud Run carries anything
   dynamic.
2. **The QA backlog gates launch.** Nothing ships to the storefront ahead of the layer QA
   sign-offs (`qa_reviewed_on` in `config/layer_registry.yaml`).
3. **License audit first.** Layer licensing and attribution obligations are reviewed before
   anything is publicly served — e.g. `landslide-arup` carries `attribution_required`
   (World Bank / GFDRR / Arup) wherever a value from it is published. **Audit completed
   2026-08-22 — see the section below.**

## Settled 2026-08-22 — licensing

The audit: [docs/licensing-audit-2026-08-22.md](../licensing-audit-2026-08-22.md)
(published artifact:
https://claude.ai/code/artifact/25b52d8c-228f-4ad6-bc15-947d8a7b44ad). Every licence was
re-verified against the live publisher record on 2026-08-22; per-dataset for ISIMIP via
the repository API `rights` field.

4. **A non-commercial licence is not a blocker — it routes the layer to the free lane.**
   NC-restricted data is served free alongside the for-sale data, never in a paid lane.
5. **Audit outcome: nothing ingested forbids serving, and exactly one licence is
   disputed.** ISIMIP3 products are CC0; 2b Lange2020 is CC0; 2b sealevelrise and 2b
   biomes are CC BY 4.0; NOAA SPC, GEBCO_2026 and Natural Earth are public domain;
   `hail-vlh` is CC BY 4.0 strictly via the article Source-Data route (its Zenodo twin is
   CC BY-NC-ND — never ingest it). `landslide-arup` remains publisher-disputed (WB DDH:
   CC BY-NC 4.0 vs energydata.info mirror: CC BY 4.0) — **recommended free-lane-only until
   the World Bank help desk resolves the record; not yet ruled.**
6. **Raw ISIMIP file mirrors are never paywalled.** ISIMIP's terms carry a blanket "sale
   of the data is strictly forbidden" sentence beside per-dataset CC0/CC BY licences;
   selling derived statistics and giving processed layers away free are both consistent
   with even the strict reading — reselling the repository's files as such is not.

## Artifacts

- [`climate-hazard-storefront.html`](climate-hazard-storefront.html) — the hand-built page
  prototype (first authored 2026-08-21 under gitignored `reports/`; tracked here so it can
  no longer be lost). It is hand-edited: **there is no generator script**.
- `location-analyses/storefront-test-sites.csv` — three synthetic warehouse test sites
  (Rotterdam, Memphis, Singapore) used for the storefront test deliveries under
  `deliveries/storefront-test/` (gitignored).

## Not yet decided / not yet built

No deploy configuration, R2 bucket, Pages project, or Cloud Run service exists in this
repo yet. UI conventions carry over from the delivery dashboard preferences (plain English
plus FAQ, opt-in heavy sections, cascading filters, stable controls, the established color
system, absent-vs-unobserved rendering).

Licensing wiring from the 2026-08-22 audit, not yet built: registry
`license` / `attribution` / free-lane fields (`LayerSpec` rejects unknown keys today);
promotion of `attribution_required` into `LAYER_ATTRS_EXPORTED` and the report builders on
the `relative_baseline` pattern; a customer-visible credit for `hail-vlh` (CC BY —
currently carried nowhere a customer sees). The `landslide-arup` free-lane call is open
(decision 5 above).
