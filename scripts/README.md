# TCFD Utility Scripts

Standalone Python scripts for data processing and visualization.

## Scripts

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

### process_let_cyclone.py

Processes `let` tropical-cyclone exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. Unlike `led`, `let` is a **continuous fraction** [0,1) per cell (area of the cell exposed to TCs), so the decadal **mean** is the expected annual exposed-area fraction. Because the ensemble is thin (1 impact model `ke-tg-meanfield` × 4 GCMs = 4 members), it applies **5×5 exponential-decay spatial smoothing** (L=0.7, cos(lat)-scaled, non-NaN-normalized) to each per-member decadal map, and a **two-tier percentile** for the zero-inflated field (zeros→1; non-zeros ranked vs the non-zero 2020s baseline → [2,100]). See WORKFLOW-ISSUES.md 2026-07-24.

```bash
python scripts/process_let_cyclone.py
```

**Input**: `data/raw/drought_let_annual/*_2006_2099.nc4`
**Output**: `data/processed/cyclone_let_annual/let_{rcp26,rcp60}_processed.nc`

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

## Dependencies

These scripts require the same dependencies as the `isimip-pipeline` package:
- xarray
- plotly
- numpy
- pandas

Install with: `pip install -e isimip-pipeline/`
