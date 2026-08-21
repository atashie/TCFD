# TCFD Utility Scripts

Standalone Python scripts for data processing, visualization and customer delivery.

Two groups, and they belong to different stages of the same product:

- **Layer processing** — `process_*.py`, `download_*.py`, `generate_maps.py`,
  `test_shared_baseline.py`. One processor per shipped layer; each writes the
  [OUTPUT-SPEC.md](../OUTPUT-SPEC.md) contract.
- **Customer delivery** — everything under *Customer delivery* below. Reads processed layers,
  never reprocesses them.

`utils/` holds the shared modules: `delivery.py` (extraction and star schema),
`report_common.py` (report loading, citations, exposure, document shell),
`report_figures.py` (inline-SVG figures), `report_profiles.py` (facet library),
`viz_common.py` (validated palette, scenario→tier colour), `spatial_extract.py`.

---

## Customer delivery

The **`/customer-delivery` skill is the entry point.** Reference documentation:
[ASSET-CATALOG.md](../ASSET-CATALOG.md) for stages 1–2, [docs/reporting/](../docs/reporting/README.md)
for stages 3–4.

Stage order is `inputs → extract → dashboard → caveats → compliance report → bespoke report`.
**Caveats runs before the reports**: the caveat set is an input to them, each report must
carry every must-disclose entry, and the verifier re-checks it.

### generate_customer_delivery.py

The driver. **Planning is the default and extraction requires `--run`**, because the resolved
asset → layer mapping has to be shown to the user before any data is touched.

```bash
python scripts/generate_customer_delivery.py --customer "Acme" --input sites.csv            # plan
python scripts/generate_customer_delivery.py --customer "Acme" --input sites.csv --run --reports
python scripts/generate_customer_delivery.py --measure-slopes    # after reprocessing a layer
```

`--run` always builds the dashboard too; skipping it marks the delivery incomplete on disk
and exits non-zero. `--reports` chains stage 4 and stage 3a.

**Output**: `deliveries/{customer}/{date}/` — a normalized star schema, `manifest.json`,
`report_config.yaml`, `dashboard.html`.

### generate_delivery_dashboard.py

Per-delivery QA dashboard (Plotly, interactive). Also importable as `build_dashboard()`.
`--check-stamp` reads the build stamp on disk — **check it before debugging a page that looks
stale**, since a cached page is indistinguishable from a fresh one.

### generate_delivery_caveats.py — stage 4

Derives the delivery's caveat set mechanically into `caveats.json` + `caveats.md`. Nothing
here is authored; it comes from the manifest, the CSVs, `report_config.yaml` and
`config/hazard_taxonomy.yaml`. Ids are semantic and stable because narratives cite them.

### generate_compliance_report.py — stage 3a

`report_compliance.html` on the IFRS S2 spine, mapped outward to CDP, ESRS E1-9 and CA SB 261.
**Fully deterministic** — no narrative, so an auditor can rebuild it byte-for-byte.

### generate_bespoke_report.py — stage 3b

`report_bespoke.html`, composed from facet profiles plus a written narrative. `--scaffold`
writes `narrative.md` and `dossier.yaml` first. The build refuses on an unfilled slot, an
uncited paragraph, or a citation that does not resolve.

### test_customer_delivery.py

**The verifier. Run it on every delivery.** Recomputes every metric from the source NetCDF
with a Gaussian weighting written independently of `spatial_extract`, then checks referential
integrity, source hashes, percentile orientation, the Climate Score, the dashboard payload,
the caveat coverage of both reports, and narrative citations. Exits non-zero on any violation.

It also checks one thing **upstream** of the delivery: that every `covered_by` entry in
`config/hazard_taxonomy.yaml` names a layer `config/layer_registry.yaml` actually has. A
dangling entry can never match the layers a delivery carries, so the family reports as NOT
ASSESSED in every report while the taxonomy claims it is covered — and both files read
correctly on their own. `permafrost-3b` sat in exactly that state for two days.

```bash
python scripts/test_customer_delivery.py deliveries/acme/20260813
```

### measure_extraction_sensitivity.py

Reproduces the two extraction-footprint measurements that reports quote as fact — the 4-cell
blend and the coordinate sensitivity. Exists because both were once quoted without a retained
receipt, which to a reviewer is indistinguishable from invention. Re-run after any change to
the extraction parameters.

---

## Layer processing

### Every processor, and where its facts live

**This index is the complete list; the long-form sections below it are not.** They cover 8
of the 27 `process_*.py` scripts, for historical reasons rather than by design, and a reader
who treats them as the inventory will conclude that layers shipped since 2026-08-12 do not
exist. **Nothing per-layer is restated here on purpose** — per
[CLAUDE.md](../CLAUDE.md) *Where a fact belongs*, a fact about one dataset lives in
[DATASET-ATTRIBUTES.md](../DATASET-ATTRIBUTES.md), its full detail in the processor's own
module docstring plus a dated [WORKFLOW-ISSUES.md](../WORKFLOW-ISSUES.md) entry, and the
processed file's global attributes outrank both. Restating any of it in this file would
create a source of truth that drifts the moment a layer is reprocessed.

Shipped layers — every one is registered in `config/layer_registry.yaml`:

| Processor | Registry layer(s) |
|---|---|
| `process_led_drought.py` | `drought-2b` |
| `process_driedarea_isimip3b.py` | `drought-3b` |
| `process_let_cyclone.py` | `cyclone` |
| `process_burntarea_isimip3b.py` | `wildfire` |
| `process_cropfailure_isimip3b.py` | `cropfailure-3b` |
| `process_heatwave_isimip3b.py` | `heatwave-3b` |
| `process_leh_isimip2b.py` | `heatwave-2b` (`status: blocked`) |
| `process_thawdepth_permafrost.py` | `permafrost-3b` |
| `process_coastal_inundation.py` | `sealevel-2b` |
| `process_fldfrcmax_isimip3b.py` | `flood-3b-flopros`, `flood-3b-40yr`, `flood-3b-none` |
| `process_tempnle_npp.py` | `conifer-npp` (asset condition, **not** a hazard) |
| `process_csoil_soilcarbon.py` | `csoil` |
| `process_prthresh.py` | `pluvial-r10mm/r20mm/r50mm/r100mm`, `pluvial-r95pd/r99pd`, `pluvial-rx1day/rx5day`, `pluvial-wetdays`, `pluvial-prcptot` — ten metrics from one ingest |
| `process_tasthresh.py` | `heatdays-hd30/hd35/hd40/hd45`, `tropicalnights-tr20/tr25`, `icedays-id`, `frostdays-fd/fdm10` — nine rungs from one ladder |

Not registered — earlier or withdrawn work. **These are not shipped layers**, and passing the
output contract never made one of them meaningful:

| Processor | State |
|---|---|
| `process_burntarea_fire.py` | **superseded 2026-08-10** by `process_burntarea_isimip3b.py` |
| `process_sug_sugarcane.py` | **withdrawn** — LPJmL does not simulate cane in the cane belt; both layers passed every check and were meaningless |
| 10 exploratory processors (qg, timber/loblolly, evgndltr, tebrsu, fish, health) | **removed 2026-08-21** — predated the current contract, never registered, so none could reach a delivery; recoverable via git history (see the dated WORKFLOW-ISSUES.md entry) |

Water Risk Index — **a different product**; none of the TCFD contract applies (see
[CLAUDE.md](../CLAUDE.md)):

| Processor | Scope |
|---|---|
| `process_water_tws.py` | total water storage |
| `process_water_variable.py` | the generic monthly water variable processor |

### process_tornado_spc.py — **the shipped CONUS tornado ladder** (NOT ISIMIP)

Builds four rungs (`all`, `f1plus`, `f2plus`, `f3plus`) from the NOAA SPC tornado database.
The only non-ISIMIP layer in the product, because ISIMIP forcing publishes 11 SURFACE
variables and nothing aloft, so the shear term every tornado index needs cannot be formed.

Answers `observational-historical-v1`, not the OUTPUT-SPEC decadal contract: no decade axis,
no slopes, scenario `observed`. **Global 0.25° grid, CONUS data, everything else NaN** — a
regional grid is not deliverable (GUARDRAILS §17). Value is a Gamma(k+½, T) posterior median
of the track-crossing rate; the pooled-median branch was measured first and is degenerate
here. `--estimator mle` reproduces the superseded hybrid; `--geometry touchdown` and
`--start-year` are the sensitivity switches.

Ingest is `download_tornado_spc.py`. Verify with `test_observational_baseline.py`, review with
`generate_tornado_qa.py`, and read
[docs/tornado-data-options-2026-08-18.md](../docs/tornado-data-options-2026-08-18.md) for the
options review, the adversarial review and every correction it forced.

### download_tornado_spc.py

Fetches the SPC database and the Natural Earth boundary used for the CONUS mask, each with a
`.json` sidecar carrying `source_url`, `sha256`, size and retrieval date — `data/` is
gitignored and ephemeral, so the sidecar is the only record of what was ingested. Resumable.
(`scripts/utils/natural_earth.py` is unusable here: it imports geopandas at module level and
geopandas is not installed in this venv.)

### test_observational_baseline.py

Contract check for `observational-historical-v1`. Separate from `test_shared_baseline.py` on
purpose — that one guards 23 projected layers and its strictness is the product, so an
observational layer failing it is correct behaviour, not something to relax.

Beyond shape, it enforces three rules this repo learned the hard way: decadal-contract
variables must be **absent, not faked to 1**; two-tier consistency is checked in **both**
directions and keyed on the observed count (keying on `median <= 0` would pass vacuously
under the posterior estimator); and a caveat attribute may not open with a negation, because
a non-empty caveat is published verbatim and `resolution_caveat: "none"` would print "none"
into a filing.

### measure_tornado_smoothing.py

Whether to smooth, and at what decay length — measured, never assumed. Splits the record
**odd vs even years** (early-vs-late would measure the reporting trend instead of sampling
noise). Primary criterion is a floor-free negative-binomial predictive; the plug-in Poisson
column is kept for comparison and systematically over-smooths. The naive "smooth both halves
and correlate" criterion is documented in the docstring as **degenerate** — it rises
monotonically with σ and always returns the largest value offered. Writes nothing.

### generate_tornado_qa.py

QA maps for the observational contract, so `qa_reviewed_on` can stop being null. Separate
from `generate_maps.py`, which selects on `ds.decade` and builds Trend and Members tabs from
variables this contract does not have. Colourbars are asserted clear of the map panels — the
first version drew them on top of the first row, which is invisible in code and obvious only
once a human opens the file.

### generate_maps.py

Generates interactive Plotly maps from processed climate data.

```bash
python scripts/generate_maps.py
```

**Input**: `data/processed/*.nc` (processed NetCDF files)
**Output**: `reports/maps/{variable}/` (HTML map files)

Features:
- Per-scenario map files (SSP126, SSP370, SSP585)
- Multiple metrics: median, percentile, trend, confidence, change, anomaly
- Scientific notation on colorbars
- Index page with navigation grid

### process_led_drought.py

Processes `led` drought exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. `led` is a **binary** per-cell annual flag, so the decadal statistic is the **mean** (drought frequency), not the median. Handles the `years since 1661` / `360_day` calendar, `(model, gcm)` ensemble members, shared 2020s baseline, and zero-inflation-aware percentile/CI. See WORKFLOW-ISSUES.md 2026-07-24.

```bash
python scripts/process_led_drought.py
```

**Input**: `data/raw/drought_led_annual/*_2006_2099.nc4`
**Output**: `data/processed/drought_led_annual/led_{rcp26,rcp60}_processed.nc`

### process_let_cyclone.py — **the shipped tropical-cyclone layer**

Processes `let` tropical-cyclone exposure (Lange 2020, ISIMIP2b `DerivedOutputData/KE-TG-meanfield`) into the TCFD output contract. Ensemble = **4 members/scenario**: 1 impact model (`ke-tg-meanfield`, an Emanuel kinetic-energy wind model) × 4 CMIP5 GCMs × {rcp26, rcp60}. **Rebuilt 2026-08-11** onto `utils/decadal_stats.py`; the 2026-07-24 version was a family-B processor (retired single `trend`, no `ols_slope`/`sen_slope`/`n_members`/`n_models`, no `members.nc`) and is fully superseded. Value-checked (see WORKFLOW-ISSUES.md 2026-08-11):

- `let` is a **continuous fraction** of each cell's land area exposed to TC winds — measured [0, 0.952], never exactly 1, `is_boolean_field` → False. Land mask identical across all 8 files (68,249 cells) and time-invariant.
- **The decadal statistic is the pooled MEAN, a declared deviation from OUTPUT-SPEC's median default.** At **97.84% annual exact-zero**, the median branch leaves 2,684 exposed cells vs 15,122 under the mean — 93% of exposed land erased — and the two-tier percentile then puts 96% of land in tier 1. `let` is the continuous analogue of `led`, which the spec already grants the mean. Declared as `decadal_statistic: pooled_mean_zero_inflated`; CI = mean ± 1 SD, clamped to [0,1]. Contrast `burntarea` at 29.2% zeros, which takes the median branch fine.
- **Aggressive 5×5 smoothing, L=2.5**, applied to each member's **annual** map before any pooling. Raw `let` is sparse one-cell-wide storm tracks and 4 members cannot average them out. The previous L=0.7 kernel kept 32.1% of weight on the centre and left tracks visible; L=2.5 puts 8.1% there, cutting roughness 0.389 raw / 0.142 previous → **0.044**, while adding **no** newly exposed cell. Radius 2 = 111 km ≈ the hurricane-force wind radius; do not widen past 3 cells. Mass conserved to +3.1%.
- **Read `ols_slope`, not `sen_slope`** — on baseline-exposed land Sen is exactly 0 on 96.9% of cells vs OLS's 1.8%, and they agree in sign only 4.9% of the time.
- Single impact model → `n_models` is 1 everywhere and the CI is **inter-GCM spread only**; there is no structural impact-model uncertainty in this layer. soc/sens uniformly `nosoc`/`co2`.
- Layer starts at the **2010s** (source spans 2006–2099), so the baseline sits at decade index **1**, not 0.

```bash
python scripts/download_let_cyclone.py    # 8 files, 58.9 MiB, resumable + sha512-verified
python scripts/process_let_cyclone.py     # ~3 min including exact Theil-Sen
python scripts/test_shared_baseline.py data/processed/cyclone_let_annual
python scripts/generate_maps.py let data/processed/cyclone_let_annual reports/maps
```

**Input**: `data/raw/cyclone_let_annual/*_let_global_annual_2006_2099.nc4`
**Output**: `data/processed/cyclone_let_annual/let_{rcp26,rcp60}_processed.nc` (+ `let_members.nc`)

### download_let_cyclone.py

Fetches the 8 `let` files (4 GCMs × rcp26/rcp60, 58.9 MiB). Resumable via Range, and **verified by sha512** against each file's `.json` sidecar — ISIMIP2b publishes `size` *and* `checksum` per file, unlike ISIMIP3b `biomes`, so this checks more than byte count. Writes `download_provenance.csv`. Serial by design: parallel requests to `files.isimip.org` get rate-limited into empty responses.

### process_burntarea_isimip3b.py — **the shipped wildfire layer**

Processes `burntarea-total` (ISIMIP3b `fire` + `biomes`) into the TCFD output contract. Ensemble = **22 members/scenario**: 5 impact models × their CMIP6 GCMs — `classic` (2), `elm-eca` (5), `lpjml5-7-10-fire` (5), `visit` (5), `mc2-usfs-r87g5c1` (5, the only **annual** burntarea publisher in the whole round) × {ssp126, ssp370, ssp585}. Value-checked 2026-08-08 (see WORKFLOW-ISSUES.md 2026-08-10):
- **Monthly → annual = SUM.** Burnt area accumulates; a mean would under-scale fire 12×. NaN-preserving — `np.nansum` maps all-NaN ocean to 0.0, which reads as "land that never burns" and tripled the apparent land mask.
- All 5 report **%** with comparable land means (1.89–3.72) → **no normalization**; no spatial smoothing (22 members). **`higher_is_worse`** → percentile not inverted. Two-tier percentile (29.2% baseline zeros): zeros → 1, non-zeros ranked against the **non-zero** baseline → [2,100].
- Values are **not clipped at 100%** — a cell that reburns within a year legitimately accumulates more (`elm-eca` reaches ~575%).
- **soc is mixed by necessity** (no token spans the ensemble); `classic` pinned to `2015soc` to avoid its cross-GCM mis-scaled `2015soc-from-histsoc` run. CO₂ treatment is uniform transient.
- **Read `ols_slope`, not `sen_slope`** — Sen is exactly 0 on 66–76% of finite cells (zero-inflation).
- `--members-only` rebuilds the `{variable}_members.nc` dashboard diagnostic from cache in seconds. Exact Theil-Sen at 22 members costs ~4.5 h for all three scenarios; **do not run the scenarios as parallel processes** (see the 2026-08-10 incident).

```bash
python scripts/download_burntarea_isimip3b.py     # 66 files, 6.03 GiB, resumable + size-verified
python scripts/process_burntarea_isimip3b.py
python scripts/test_shared_baseline.py data/processed/wildfire-isimip3b_burntarea-total_annual
```

**Input**: `data/raw/wildfire-isimip3b_burntarea-total_annual/*_2015_2100.nc`
**Output**: `data/processed/wildfire-isimip3b_burntarea-total_annual/burntarea_{ssp126,ssp370,ssp585}_processed.nc` (+ `burntarea_members.nc`)

### process_burntarea_fire.py — **SUPERSEDED 2026-08-10**

> Superseded by `process_burntarea_isimip3b.py` above: newer experiment generation (SSP over
> RCP), deeper ensemble in both dimensions (22 vs 12 members, 5 vs 3 impact models), and a
> high-forcing scenario. Kept for provenance — it documents the ISIMIP2b generation's
> per-member metadata defects, which are real findings about *that* data. Do not extend it,
> and do not treat its framing decisions as precedent (GUARDRAILS §9: measure, never inherit).

Processes `burntarea` (wildfire burnt-area fraction, ISIMIP2b `biomes`) into the TCFD 6-value-class format. This is the **direct biophysical fire signal** — the fraction of each cell burned per year, in **%** — and is distinct from the Lange 2020 `lew` *exposure* family. Ensemble = 3 annual fire models (`lpj-guess`, `lpjml`, `mc2-usfs`) × 4 GCMs × {rcp26, rcp60, rcp85} = 12 members/scenario. Key points, all value-checked (see WORKFLOW-ISSUES.md 2026-07-24):
- All 3 models report **%** [0,100] → **no normalization** (equal-weight "model democracy" in raw %); inter-member spread = CI. (Contrast the water-index TWS, which needed normalization.)
- Per-member metadata **diverges**: `mc2-usfs` uses `days since 1661` (others `years since`); `lpj-guess` mislabels its `long_name` as "Fire Return Interval" (data are burnt %) and floors at 0.1%; `mc2-usfs` is a zero-inflated, ~5–7× hotter outlier. The time parser handles days/months/years per calendar.
- Slopes follow [OUTPUT-SPEC.md](../OUTPUT-SPEC.md): `ols_slope` + `sen_slope`, fitted over an expanding window from the 2020s baseline (baseline panel = NaN). The baseline-anchored two-point rate is retired, and `trend x elapsed_decades == change map` no longer holds. No spatial smoothing (12-member ensemble is thick).

```bash
python scripts/process_burntarea_fire.py
```

**Input**: `data/raw/wildfire_burntarea_annual/*_2006_2099.nc4`
**Output**: `data/processed/wildfire_burntarea_annual/burntarea_{rcp26,rcp60,rcp85}_processed.nc`

### process_driedarea_isimip3b.py — **the shipped ISIMIP3b drought layer**

`driedarea` from ISIMIP3b `DerivedOutputData/Heinicke2026`, ssp126/370/585. A **sibling** of
the ISIMIP2b `led` layer, not a replacement — deeper ensemble versus newer scenarios, and
both are current. Resolves the minimum-model mask differently from `led` (full union rather
than ≥2 models) on measured evidence; re-measure that rule per layer, never inherit it.

**Output**: `data/processed/drought-isimip3b_driedarea_annual/driedarea_{ssp126,ssp370,ssp585}_processed.nc`

### process_tempnle_npp.py — **the shipped conifer productivity layer**

Temperate needleleaf-evergreen stand NPP from ISIMIP2b `biomes` (CLM45 + ORCHIDEE + LPJmL),
rcp26/60/85. Reported on a **measured per-tile denominator** behind a 2% cover presence mask.
`higher_is_better`, so its percentile is inverted at processing time and rising productivity
reads as falling risk — state the CO₂-fertilisation caveat wherever it is reported.

**Output**: `data/processed/conifer-temperate_npp-tempnle_annual/npp-tempnle_{rcp26,rcp60,rcp85}_processed.nc`

### process_csoil_soilcarbon.py

Processes `csoil-total` (soil organic carbon stock, ISIMIP3b `biomes`) into the TCFD output contract — the direct **subsurface carbon-storage** signal (distinct from the vegetation pools `cveg`/`croot`/`cvegbg` and the net-sink flux `nbp`). **Rebuilt 2026-08-15**; the 2026-07-25 build is superseded, its output had been lost from disk, and its ensemble was understated. Ensemble = **4** models (`classic`, `jules-es-vn6p3`, `lpjml5-7-10-fire`, `mc2-usfs-r87g5c1`) × their CMIP6 GCMs (2 + 5 + 5 + 5) × {ssp126, ssp370, ssp585} = **17 members/scenario**. Value-checked 2026-08-15 by `check_csoil_nature.py` (see WORKFLOW-ISSUES.md):
- All 4 report **kg C m⁻²**, but magnitudes are **not** as comparable as the 3-model build recorded: 2020s medians 5.70 / 7.55 / 10.45 / 16.82, a **2.92×** spread with `lpjml` carrying the upper tail alone (p95 91.1 vs 16.5–27.4). **No normalization** anyway — the models form a gradient rather than clusters, the disagreement about stock size *is* the structural uncertainty the CI carries, and the units are checkable against soil inventories. Layer starts at the 2020s baseline (ISIMIP3b csoil begins 2015 — no full 2010s).
- **Direction is `higher_is_better`** (stored carbon is an asset; the risk is **loss**) → percentile **inverted** (low stock / decline → high risk; red = carbon loss in maps). Note the consequence: a **zero-carbon cell scores as maximum risk** — the Sahara sits at the 96th risk percentile and Greenland's ice sheet at the 97th. Correct arithmetic, inverted conclusion; see the registry `delivery_note`.
- **Mixed CO₂**: `jules-es-vn6p3` publishes only its fixed-2015-CO₂ run. **Do not call that trend muted** — measured, those members are the *largest relative losers* (−4.37% vs −2.75% `lpjml`, −0.05% `mc2-usfs`, +0.79% `classic`), because removing fertilization removes litter input. All members hold land use fixed, so the layer cannot see management-driven loss. Retained as 17 members per user decision.
- **`pooled_median`** retained on measurement (largest adjacent-model gap 0.59 × IQR, below the ~1.0 that triggers OUTPUT-SPEC's fourth branch); **no spatial smoothing** on a model-stratified split-half (roughness 0.148, halves 0.149/0.150, r = 0.992); **read `sen_slope`** (measured: 6.8% Sen-zero, but only 70.3% sign agreement — the lowest of any Sen layer).
- **`classic` is natively 1°**, replicated 2×2 with a one-cell longitude offset — measured 100% adjacent-cell ties on odd column pairs, 0% on even. The other three are genuinely 0.5°. Visible only on the per-member contact sheet.

```bash
python scripts/download_csoil_isimip3b.py          # 51 files, publisher-sha512 verified
python scripts/check_csoil_nature.py               # GUARDRAILS §9 + §12, before processing
python scripts/process_csoil_soilcarbon.py --jobs 6
python scripts/process_csoil_soilcarbon.py --members-only   # rebuild the Members diagnostic only
```

**Input**: `data/raw/soilcarbon_csoil_annual/*_csoil-total_global_annual_2015_2100.nc`
**Output**: `data/processed/soilcarbon_csoil_annual/csoil_{ssp126,ssp370,ssp585}_processed.nc`

### download_reduce_tasthresh_isimip3b.py + process_tasthresh.py — **the threshold ladder (9 layers)**

Two stages, built 2026-08-16. Stage 1 streams ISIMIP3b bias-adjusted daily `tasmax`/`tasmin` and reduces it to annual threshold-exceedance counts; stage 2 turns those counts into nine OUTPUT-SPEC layers. **12 GCMs × ssp126/370/585 = 36 members, uniform across scenarios.**

**THE RAW DATA IS NEVER RETAINED, and that is not an optimisation.** The ingest moves **~1.34 TB** against ~600 GB of usable disk, so each chunk is downloaded, verified against the publisher sidecar's sha512, reduced, and **deleted** — peak disk ~10 GB. Provenance (url, size, sha512, UTC) is written to `data/interim/tasthresh/download_provenance.csv` **before** deletion, so the ingest stays auditable against files we no longer hold. A declared deviation from the `data/raw/{layer_id}/` retention convention, recorded in the interim file's attributes. Because storage never became binding, the GCM pool did **not** have to be cut.

**The whole ladder comes from ONE pass.** Once a day is read, testing it against nine thresholds costs nothing, so a rung we do not ship today needs no second 1.34 TB. `hd30/hd35/hd40/hd45` (tasmax >30/35/40/45 °C), `TR20/TR25` (tasmin >20/25 °C, ETCCDI tropical nights), `ID` (tasmax <0 °C, ETCCDI ice days), `FD` (tasmin <0 °C, ETCCDI frost days), `FDm10` (tasmin <−10 °C, **not** a standard index).

- **HDF5 IS NOT THREAD-SAFE.** Four workers in a `ThreadPoolExecutor` crashed the interpreter twice with SIGBUS inside `H5S_select_iter_init`; the single-worker smoke test passed, which is why it survived to the full run. Parallelism comes from `ProcessPoolExecutor`. Do not "simplify" it back to threads.
- **Chunk-level `.npz` checkpoints**, not member-level: a member dying on its ninth chunk would otherwise re-download ~30 GB. Three members timed out on the first pass and resumed for the cost of the ~38 chunks they were missing.
- **`landseamask_no-ant.nc` (65,797 cells), Antarctica EXCLUDED** — measured, not conventional. See GUARDRAILS §14: Antarctica is 29% of the full mask and carries 98.85% of `FD`'s ceiling censoring.
- **The counts are finite over the whole globe including ocean**, so `isfinite` is **not** a mask here — the landseamask does all the work.
- **All nine take `pooled_mean_zero_inflated`** (the median branch would publish exactly 0 for 44% of land on `hd45`) and **all nine read `ols_slope`** — chosen from the failure modes, since `n_members` is exactly 12 in every cell so ols's bias has no mechanism. See DATASET-ATTRIBUTES.md.
- **The calendar risk did not materialise**: all 32 members measured `proleptic_gregorian`, so no days-per-year correction is applied. Stage 1's `*_days_per_year_HETEROGENEOUS` attribute is an artifact of per-chunk arithmetic (a 6-year chunk gives 365.33, a 10-year one 365.30) and is **not** evidence of a mixed calendar — read `*_calendar`.
- **Containment is the check that matters** and no per-layer verifier can do it: `hd30 ≥ hd35 ≥ hd40 ≥ hd45`, `TR20 ≥ TR25`, `FD ≥ FDm10`, `FD ≥ ID`. **0 violations** across all three scenarios.

```bash
.venv/bin/python3 scripts/check_tas_thresholds_nature.py            # GUARDRAILS §9 + §12, before anything
.venv/bin/python3 scripts/download_reduce_tasthresh_isimip3b.py --plan
.venv/bin/python3 scripts/download_reduce_tasthresh_isimip3b.py --run --workers 4
.venv/bin/python3 scripts/process_tasthresh.py --plan --rung hd35   # reference sites, no slope stage
.venv/bin/python3 scripts/process_tasthresh.py --run --jobs 8       # all nine rungs, ~1.2 h
.venv/bin/python3 scripts/emit_tasthresh_registry.py                # PRINTS the registry block
```

**Input**: `data/interim/tasthresh/{GCM}_{scenario}_counts.nc` (36 files, ~1.6 GB total)
**Output**: `data/processed/{heatdays|tropicalnights|icedays|frostdays}-isimip3b_{rung}_annual/{rung}_{scenario}_processed.nc`

`emit_tasthresh_registry.py` **prints** the `layer_registry.yaml` block read from the built files and does **not** write the registry — that file carries human decisions (`status`, `delivery_note`) alongside measured ones, and a script that rewrote it in place would be one bad run from clobbering the other fourteen layers. It also asserts that a rung's three scenario files agree on `recommended_slope`, which is the defect that forced a rebuild.

### The precipitation pipeline — **the ten pluvial layers**

Three stages plus a driver, built 2026-08-18. **14 GCMs x ssp126/370/585**, ~911 GB moved and
none retained (peak disk ~10 GB).

```bash
.venv/bin/python3 scripts/check_pr_nature.py                 # GUARDRAILS 9 + 12, first
.venv/bin/python3 scripts/pr_baseline_percentiles.py --run    # stage 0: per-cell wet-day p95/p99
.venv/bin/python3 scripts/download_reduce_prthresh_isimip3b.py --run   # stage 1: annual statistics
.venv/bin/python3 scripts/process_prthresh.py --run --jobs 8  # stage 2: OUTPUT-SPEC layers
.venv/bin/python3 scripts/run_pr_pipeline.py --attach         # drives 0 -> 1 with bounded retries
bash scripts/qa_prthresh_dashboards.sh                        # contract + dashboards + containment
.venv/bin/python3 scripts/emit_prthresh_registry.py           # PRINTS the registry block
```

**STAGE 0 EXISTS BECAUSE THE RELATIVE RUNGS CANNOT BE COUNTED IN ONE PASS.** `R95pD`/`R99pD`
count days above each cell's OWN 2020s wet-day percentile, and that threshold is not knowable
from the chunk in front of you. Stage 0 reads only the two baseline chunks per member (~147 GB),
accumulates a per-cell log-spaced histogram (256 bins over 1-1000 mm/day, ~67 MB), and writes
`baseline_percentiles.nc`. Stage 1 refuses to start without it — running anyway would silently
emit empty relative rungs. Percentile reconstruction from the histogram was validated against
exact percentiles: worst error 2.63% against a 2.7% bin width, i.e. at the resolution floor.

- **The percentile is POOLED across members, not per member.** A per-member threshold would pin
  every member at exactly 5% baseline exceedance and collapse the inter-member spread the CI
  carries. It is over **WET days (>= 1 mm)**, not calendar days — on an all-days basis the
  metric is exactly 0 mm wherever wet days fall below 5% of the year.
- **Rx5day windows cross chunk boundaries.** Each chunk carries its last 4 days forward and the
  carry is stored INSIDE the checkpoint, because a resume that lost it would reintroduce the bug
  invisibly. Without the carry a synthetic test understated **42.5% of cells** in every year
  following a boundary; the real data is verified by `Rx5day >= Rx1day`, 0 violations.
- **Two statistic branches, split by unit family**: seven day counts take
  `pooled_mean_zero_inflated` (the median erases 53.5% of land on R50mm); three mm metrics keep
  `pooled_median` (it erases nothing). The slope follows the same boundary: counts read
  `ols_slope` (sen collapses on 41.8-100% of active cells), mm metrics read `sen_slope`
  (0.0-1.5%, and an annual maximum is heavy-tailed where robustness matters most).
- **599 cells are masked for the relative rungs only** and must be NaN, never 0 — stage 1 counts
  `hit & usable`, so an unusable cell accumulates a legitimate-looking count of ZERO. Stage 2
  applies the mask on load and **asserts** it before writing.

**Input**: `data/interim/prthresh/{GCM}_{scenario}_pr.nc` (42 files, 7.1 GB)
**Output**: `data/processed/pluvial-isimip3b_{metric}_annual/{metric}_{scenario}_processed.nc`

### test_shared_baseline.py

The layer contract verifier — the counterpart to `test_customer_delivery.py` one level up.
Run after every processing run; exits non-zero on any OUTPUT-SPEC violation.

```bash
python scripts/test_shared_baseline.py data/processed/{layer_dir}
```

## Dependencies

Same as the `isimip-pipeline` package: xarray, plotly, numpy, pandas, pyyaml.
Install with `pip install -e isimip-pipeline/`.

**Not present on this machine, and the report tooling is built around their absence**: no
pandoc, weasyprint, wkhtmltopdf, kaleido, node or headless Chrome. Report figures are
therefore inline SVG with zero JavaScript, and PDF is produced by printing from a browser.
