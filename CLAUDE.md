# ISIMIP Climate Data Pipeline - Development Guide

This codebase pulls, processes, and validates ISIMIP climate model data for two distinct downstream products. Each product has its own processing logic, output format, and tooling. Workflows are consistent within each product class but differ across them. **Never confuse them.**

## Two Data Products

### 1. TCFD / CDP Reporting (annualized decadal statistics)

Processes ISIMIP data into annualized decadal statistics for physical climate risk assessment.

- **Purpose**: TCFD and CDP climate risk disclosures (timber, fisheries, health, agriculture, drought, tropical cyclones, wildfire, soil/subsurface carbon, etc.)
- **Output contract**: [OUTPUT-SPEC.md](OUTPUT-SPEC.md) is authoritative — do not restate it elsewhere. `{variable}_{scenario}_processed.nc` on `(decade, lat, lon)`, carrying `median`, `lower_ci`, `upper_ci`, `percentile`, `ols_slope`, `sen_slope`, `n_members`, `n_models`.
- **Shared statistics**: `scripts/utils/decadal_stats.py` implements the contract; every processor calls it rather than rolling its own. `scripts/process_csoil_soilcarbon.py` is the **reference implementation**.
- **Verification**: `python scripts/test_shared_baseline.py {processed_dir}` after every run — exits non-zero on any contract violation.
- **Customer delivery**: the **`/customer-delivery` skill is the entry point** — invoke it before improvising a client deliverable. It owns the pipeline (1 inputs → 2 CSV + dashboard → **4 caveats → 3a compliance report → 3b bespoke report**), and stage status is stamped into each delivery's `manifest.json`. **Stage 4 runs before Stage 3 deliberately**: the caveat set is an *input* to both reports (each must carry every `must_disclose` entry, enforced in the builder and re-checked by the verifier), not a summary of them. **Stage 2 always produces both the CSVs and `dashboard.html`** from one command; a delivery missing either fails `scripts/test_customer_delivery.py`. [ASSET-CATALOG.md](ASSET-CATALOG.md) covers stages 1–2 and **[docs/reporting/](docs/reporting/README.md) covers stages 3–4** — both are reference docs the skill points at; do not restate them here.
- **Reports are two documents from one set of numbers.** `generate_compliance_report.py` is an IFRS S2-spined, **fully deterministic** artifact (mapped outward to CDP 3.1.1, ESRS E1-9, CA SB 261); `generate_bespoke_report.py` composes **facet profiles** (asset × region × persona × vertical × use case × company) with a written narrative. Both read the delivery through `scripts/utils/report_common.py` so they cannot disagree about a number, and both refuse to render without every `must_disclose` caveat. The narrative's claims are enforced, not trusted: an unfilled slot, an uncited paragraph, or a citation that does not resolve to a real CSV row or a `dossier.yaml` source fails the build. **Profiles guide the narrative and are never pasted into it** — comments are stripped before rendering and a report containing `<!--` fails the verifier. Figures are inline SVG with zero JavaScript (`scripts/utils/report_figures.py`), because this machine has no PDF toolchain at all.
- **A report section may report NOTHING, and should, where the method has not been agreed with the user.** Deferred decisions live in `report_common.TBD_SECTIONS`, render identically into both reports via `tbd_block()`, and the verifier **fails** any report that publishes a figure whose method is still deferred. **Currently deferred: the IFRS S2 29(c) vulnerable-asset count** — `percentile` is a global-relative *exposure* rank and "vulnerable" is a claim about susceptibility to *harm*; nothing connects them, and on the worked example the first attempt's count fell from 4 of 5 assets to 1 of 5 as the threshold moved from 60 to 90. Err toward an explicit gap over a fast answer. **A measurement quoted in a filing needs a reproduction, not a memory** — `scripts/measure_extraction_sensitivity.py` exists because two footprint figures had no retained receipt, which to a reviewer is indistinguishable from invention. **Facets accept lists and `region` is validated against the delivered locations** — `assert_region_coverage()` refuses to build if a site is uncovered, because a report framed on the wrong region is wrong about the first thing a reader checks.
- **A layer scoring DEPARTURE FROM A LOCAL BASELINE must say so in every report, and the registry enforces it.** Such a layer's high score means *"unusual for this place"*, never *"bad in absolute terms"* — and the difference **reverses the ranking a reader expects**. Because a reader can be wrong about that while every number in front of them is correct, it is machine-enforced rather than left to prose: a layer sets `relative_baseline: true` + `relative_baseline_note` in `config/layer_registry.yaml`, `generate_delivery_caveats.py` promotes it to a **`must_disclose`** caveat, both reports refuse to render without it, and `load_registry()` fails on a blank note. It is printed as `READ AS RELATIVE:` **above** the ordinary `NOTE:` in the Stage 1 plan, because that is where the mapping is agreed with the customer. When flagging a new one, **check whether its siblings are flagged too** — flagging one alone implies the others are absolute. Which layers, and the worked example, are in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md). (Added 2026-08-13; before it, two layers' own notes said "Must be stated in any customer narrative" while the machinery filed it as optional — the rule was written down and not wired up.)
- **Every layer is a REGIONAL result at 0.5° (~55 km), never a site result — and where that gap matters it is machine-enforced.** A site's value is its cell's value. That is fair for hazards varying over comparable distances and **wrong for hazards that turn on fine terrain**: on `sealevel-2b` a quayside and a hillside 2 km apart share a number, because coastal inundation depends on metres of elevation over hundreds of metres of ground. Such a layer sets `resolution_caveat` in its processed file, `generate_delivery_caveats.py` promotes it to a **`must_disclose`** caveat, and both reports refuse to render without it — the same mechanism as `relative_baseline`, for the same reason: the number is correct and the reader's conclusion does not follow from it. These layers **screen**; they rank which sites deserve investigation and cannot support a design elevation or an asset-level estimate. Do not confuse coarse *support* with coarse *inputs* — `sealevel-2b` reads terrain at 15″ and uses the full within-cell elevation distribution. (Added 2026-08-14.)
- **Hazard coverage is declared, positively and negatively.** `config/hazard_taxonomy.yaml` names the physical hazards a disclosure is expected to address and records which are covered — read it rather than restating a count here, which goes stale the moment a layer ships. A report that lists what it assessed and stops reads as though the rest was found immaterial, so a "hazards not assessed" section is mandatory in every report. **The family list is OURS, not a standard's** — only the acute/chronic split is standards-derived — so reports name hazards covered and not covered and **never quote a fraction**, which would be a ratio against a denominator we chose. The file splits `customer_note` (rendered) from `materiality_note` / `blocker` / `isimip_candidate` (internal — they carry repo paths, dataset defects and the word UNVERIFIED). Do not swap them. Shared chart vocabulary (validated palette, scenario→forcing-tier colour mapping, symmetric-limit helper) lives in `scripts/utils/viz_common.py`; `generate_maps.py` stays the tool for **gridded** layers and the two are deliberately not merged — opposite payload profiles. `scripts/generate_customer_delivery.py` + `config/asset_catalog.yaml` (asset type → layers) + `config/layer_registry.yaml` (layer → disk location and which slope to read). Output is a normalized star schema in gitignored `deliveries/{customer}/{date}/`; the Looker 28-column `Export-Key.csv` contract is **retired** (user decision 2026-08-12). Planning is the default and extraction needs `--run`, because the resolved asset→layer mapping must be shown to the user before any run. An asset type absent from the catalog is an **error, never a default**.
- **Tooling**: `isimip-pipeline` CLI + `scripts/process_*.py`
- **Skills**: `/isimip-search-download`, `/isimip-process-visualize`, `/isimip-extract-aggregate`. **Invoke the skill before improvising.** Any question of the form "what data exists for X / what could we process / is X available" is `/isimip-search-download` — load it *first*. It carries the enumeration mechanics (serial harvest, `.nc4?` matching, `DerivedOutputData/`) and coverage facts that this file does not, and it is **authoritative over this file** where they disagree. Hand-rolled `curl` sweeps that skip it have re-derived documented findings and drawn absence conclusions from rate-limited empty listings. **Start from the skill's publication map, not from a variable name** — derived publications are named after first authors and identified only by their files, and entering a path we already know the name of is retrieval, not a search. That distinction hid the CaMa-Flood inundation suite for three weeks behind a path we had already walked (GUARDRAILS §13).
- **Visualization**: `scripts/generate_maps.py` — six tabs (`Median | Percentile | Trend | Confidence | Anomaly | Members`); diverging slope/change panels are zero-centred on the 95th percentile of |value| and auto-reverse to red=worse when a file sets `percentile_direction: higher_is_worse`. Conventions and the browser-payload limits live in the `/isimip-process-visualize` skill — read it before changing a dashboard.
- **Key concepts**: Shared 2020s baseline, adaptive windowing, pooled (year × member) decadal statistics, dual OLS + Theil-Sen slopes on an expanding window, percentile-of-score ranking
- **Two slopes, not one** — they fail in **opposite** regimes, so neither is safe alone. `sen_slope` collapses to exactly 0 on zero-inflated hazards; `ols_slope` absorbs between-member level offsets as trend when member coverage is uneven (measured +40% bias). Read `ols_slope` on zero-inflated fields, `sen_slope` otherwise, and treat **disagreement between them as the signal that a cell's trend is not robust**. Which slope each shipped layer takes is **measured, not inferred from `field_nature`** — the table and the re-measurement command are in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md). This supersedes the retired single `trend` variable and the never-emitted `significance` class. **The "opposite regimes" premise fails on a CENSORED field** — where a layer is pinned at a bound both estimators go to ~0 and *agree*, so agreement there is ambiguous between "no trend" and "maximally exposed, permanently", and the disagreement rule gives no warning. Check the share of cells sitting at the bound before reading any slope; `heatwave-3b` is the shipped instance (45.9% at 1.0 by ssp585 2090s) and the treatment — declare, never re-estimate — is in [OUTPUT-SPEC.md](OUTPUT-SPEC.md).
- **Per-layer facts are NOT in this file.** [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md) is the index of what each shipped layer is and how it was framed; the full detail — per-member metadata quirks, the measurements behind each decision, reference sites — lives in **the processor's module docstring** and its dated **WORKFLOW-ISSUES.md** entry, and the file's own **global attributes** outrank all of them. Read the one layer you are working on. Do **not** read the others for comparison unless the task genuinely depends on it — see *Scope discipline* below.
- **A recorded negative is a claim, not evidence — and so is an understated positive (GUARDRAILS §11).** The catalog asserted "no ISIMIP3b/SSP drought" for 18 days and it was false. It separately recorded a one-model ensemble that had eight, taken from a row-capped listing, and that number became the stated reason a whole hazard was unassessable. Every negative needs `verified_absent_on: "<date> — listed <URL>"` or the word `UNVERIFIED`; every **capped or paginated** result is `UNVERIFIED` for coverage *in either direction*. A negative about a *code* is not a negative about a *hazard* — re-issues get renamed.
- **Availability questions stop at the inventory.** Report the model × GCM × scenario matrix with cadence, volume and the open decisions, then **end the turn**. No recommendation, no ranking, no `AskUserQuestion` to force a pick — dataset choice is a discussion the user opens, and we decide together from the table. Give your read on the trade-offs when asked. (Rule added 2026-08-08; see the skill's *Answering "what data exists for {hazard}?"*.)
- **Scope discipline.** When the user asks about **one** hazard — what data exists, what to process, how to frame it — answer within that hazard. Enumerate its options from the repository and judge each on **its own** properties: cadence, ensemble depth, units, data nature, scenario coverage, volume. Do **not** justify or rank a dataset by reference to an unrelated shipped layer ("aligns with `csoil`", "same shape as `let`"). Cross-layer comparison is warranted only when the user asks for it, or when it is a real constraint (a shared processor, the OUTPUT-SPEC contract, a guardrail). Precedent from another hazard is a *hypothesis to re-verify*, never a reason to prefer a dataset.
- **Cross-cutting framing rules** (they generalize; the per-layer instances do not — those are in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)): data nature is **measured, never inferred** — from the variable name, its `long_name`, a sibling member, or another layer (GUARDRAILS §9). The measured nature then determines the statistic, the CI, and the percentile treatment per [OUTPUT-SPEC.md](OUTPUT-SPEC.md). Thin ensembles get 5×5 exponential-decay spatial smoothing — but **the decay length is a per-layer measurement, not a constant**, and so is whether to smooth at all: when roughness and ensemble depth disagree, a **split-half test** settles it (halve the ensemble; if roughness barely moves and the halves correlate, the roughness is real structure, not sampling noise). Zero-inflated hazards get a two-tier percentile (zeros→1, non-zeros ranked vs the non-zero 2020s baseline). Variables where a decline is the risk set `higher_is_better` and invert the percentile.
- **Four decadal-statistic branches, not two.** Continuous → pooled median/IQR; boolean {0,1} → pooled mean±SD; **extreme zero-inflation** → `pooled_mean_zero_inflated`, mean±SD on a *continuous* field; and a **multimodal ensemble** → `pooled_mean_multimodel`. The third and fourth exist because the median assumes something the pool does not always satisfy: the third, that the pool is not degenerate at zero (it failed hard — the median branch erased 93% of exposed land on one layer); the fourth, that the pool has *one mode* — when members separate into clusters the median selects the larger cluster and jumps when the balance tips (measured: a spatial median moving 0.40 → 0.93 in one decade with no physical process behind it). Both are **declared** deviations: measure what the median branch *would* have published, record it in `decadal_statistic_rationale`, and get the user's decision. Never take either to improve contrast. **A threshold on the central value inherits the branch** — the same ensemble reported 7.97 vs 0.40 M km² of "fully thawed" area under median vs mean — so derive such counts from the **member share**, which is invariant. Each adoption is its own decision on its own measurement and **the count is not a precedent**; which layers qualify, and on what numbers, is in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md).
- **Extreme-event exposure families** — `config/isimip_search_catalog.yaml` and the `/isimip-search-download` skill are **authoritative over this file** for repository coverage; [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md) orients. Two generations exist, and the re-issue is named by **hazard word, not `le*` code**. Never conclude a family has no newer-round version without listing `DerivedOutputData/`. **And a family having no re-issue is not the hazard having no newer data** — ISIMIP3b's newest tropical-cyclone product sits outside `DerivedOutputData/` entirely. List every directory level; scope each absence claim to exactly what was enumerated (2026-08-11).

### 2. Water Risk Index (20 value types, monthly)

Processes monthly ISIMIP data for 6 water variables into per-month ensemble means plus annual quantile breakpoints, feeding a dedicated Water Risk Index product.

- **Purpose**: Water risk scoring (total water storage, discharge, runoff, evapotranspiration, soil moisture, precipitation)
- **Tooling**: Standalone scripts only (NOT the `isimip-pipeline` CLI)
- **None of the TCFD contract applies here** — no trends, no percentile scoring, no kernel smoothing, no shared baseline. Never carry a decision across the two products.
- **Per-variable attributes** — output path and dims, the 20 value types, the normalization rule, aggregation and units per variable — are in **[DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)**, with `config_water_variables.py` authoritative for `unit_conversion_factor`.
- **QA/QC**: `validate_water_tws.py` (quantile ordering, annual mean consistency, seasonal sanity, cross-scenario checks); `compare_water_index.py` (trend-focused RCP vs SSP HTML comparison with Theil-Sen slope maps and spatial Spearman R²)

## Project Structure

```
TCFD/
├── isimip-pipeline/          # Python package (CLI) — TCFD/CDP workflow
│   ├── src/isimip_pipeline/  # cli, search, download, processing, visualization
│   └── tests/                # 299 tests
│
├── scripts/                  # Standalone scripts — both workflows
│   ├── utils/                # Shared utilities (land_mask, water_index_compare)
│   ├── process_qg.py         # TCFD: example annualized processor
│   ├── generate_maps.py      # TCFD: interactive Plotly maps
│   ├── config_water_*.py     # Water: variable configuration
│   ├── process_water_*.py    # Water: processing scripts
│   ├── validate_water_tws.py # Water: QA/QC validation
│   ├── compare_water_index.py # Water: RCP vs SSP comparison report
│   ├── diagnose_tws_models.py # Water: model distribution diagnostics
│   └── download_water_*.py    # Water: per-variable ISIMIP download scripts
│
├── config/                   # ISIMIP search catalog cache
├── data/                     # Raw + processed data (gitignored)
├── reports/                  # Generated reports (gitignored)
└── _deprecated/              # Archived legacy R code
```

## CLI Quick Reference (TCFD/CDP only)

```bash
isimip-pipeline search         # Search ISIMIP repository
isimip-pipeline download       # Batch download
isimip-pipeline process        # Process raw data
isimip-pipeline report         # Generate QA report
isimip-pipeline run            # Complete pipeline
isimip-pipeline find           # Search local processed datasets
isimip-pipeline catalog        # Manage ISIMIP catalog
isimip-pipeline cleanup        # Delete raw data after verification
```

## Skills

| Skill | When to Use | Product |
|-------|-------------|---------|
| `/isimip-search-download` | Searching ISIMIP repository, downloading data | TCFD/CDP |
| `/isimip-process-visualize` | Processing annualized NetCDF, generating QA reports | TCFD/CDP |
| `/isimip-extract-aggregate` | Extracting data by location/region for CSV export | TCFD/CDP |

## Key Documentation

| Document | Purpose |
|----------|---------|
| [GUARDRAILS.md](GUARDRAILS.md) | Critical rules that must never be violated |
| [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md) | **What each dataset is** — shipped layers, their measured framing, the Water Index variables. Per-dataset facts belong here, not in this file. |
| [OUTPUT-SPEC.md](OUTPUT-SPEC.md) | The TCFD output contract (authoritative; do not restate) |
| [ASSET-CATALOG.md](ASSET-CATALOG.md) | Customer delivery, stages 1–2 |
| [docs/reporting/](docs/reporting/README.md) | Customer delivery, stages 3–4 |
| [WORKFLOW-ISSUES.md](WORKFLOW-ISSUES.md) | Incident log and resolutions |

**Where a fact belongs.** A *rule that generalizes* goes in this file or GUARDRAILS.md. A
*fact about one dataset* goes in DATASET-ATTRIBUTES.md, and its full detail in the
processor docstring plus a WORKFLOW-ISSUES entry. A *fact a processed file already carries*
goes nowhere — it is read from the file's global attributes at delivery time, so restating
it creates a second source of truth that drifts the moment a layer is reprocessed.
