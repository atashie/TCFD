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

### process_qg.py

Processes raw groundwater runoff (qg) data from ISIMIP.

```bash
python scripts/process_qg.py
```

**Input**: `data/raw/*.nc` (raw NetCDF files from ISIMIP)
**Output**: `data/processed/qg_*.nc` (processed files by scenario)

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

Processes `csoil-total` (soil organic carbon stock, ISIMIP3b `biomes`) into the TCFD 6-value-class format — the direct **subsurface carbon-storage** signal (distinct from the vegetation pools `cveg`/`croot`/`cvegbg` and the net-sink flux `nbp`). Ensemble = 3 models (`classic`, `jules-es-vn6p3`, `mc2-usfs`) × their CMIP6 GCMs (2 + 5 + 5) × {ssp126, ssp370, ssp585} = 12 members/scenario. Value-checked 2026-07-25 (see WORKFLOW-ISSUES.md):
- All 3 report **kg C m⁻²** with **comparable magnitudes** (2020s medians ~5.8/7.7/10.3) → **no normalization** (equal-weight "model democracy"); inter-member spread = CI. Layer starts at the 2020s baseline (ISIMIP3b csoil begins 2015 — no full 2010s).
- **Direction is `higher_is_better`** (stored carbon is an asset; the risk is **loss**) → percentile **inverted** (low stock / decline → high risk; red = carbon loss in maps).
- **Mixed CO₂**: `jules-es-vn6p3` publishes only its fixed-2015-CO₂ run for csoil, so its trend is muted; `classic`/`mc2-usfs` are transient. All hold land use fixed. Retained as 12 members per user decision. Baseline-anchored trend (kg C m⁻² decade⁻¹); no spatial smoothing (thick ensemble).

```bash
python scripts/process_csoil_soilcarbon.py
```

**Input**: `data/raw/soilcarbon_csoil_annual/*_csoil-total_global_annual_2015_2100.nc`
**Output**: `data/processed/soilcarbon_csoil_annual/csoil_{ssp126,ssp370,ssp585}_processed.nc`

### generate_qa_report.py

Generates QA reports for processed data.

```bash
python scripts/generate_qa_report.py
```

**Output**: `reports/qg_qa_report.json`

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
