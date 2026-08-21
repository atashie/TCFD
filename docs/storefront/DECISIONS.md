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
   (World Bank / GFDRR / Arup) wherever a value from it is published.

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
