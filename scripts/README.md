# TCFD Utility Scripts

Outputs go to S3, not local disk — see [STORAGE.md](../STORAGE.md). Scripts still referencing `data/` are pending the Python rebuild.

Standalone Python scripts for data processing and visualization.

## Scripts

### generate_maps.py

Generates interactive Plotly maps from processed climate data.

```bash
python scripts/generate_maps.py
```

**Input**: a published layer in S3 (`…/layers/{layer_id}/current/`, gate-verified)
**Output**: `…/layers/{layer_id}/{version}/maps/` (HTML map files, ungated/regenerable)

```bash
python scripts/generate_maps.py wildfire_burntarea_annual
python scripts/generate_maps.py drought_led_annual --version v2026-07-27_3412446
python scripts/generate_maps.py cyclone_let_annual --local-only
```

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

**Input**: `s3://…/TCFD/raw/isimip/drought_led_annual/*_2006_2099.nc4`
**Output**: published layer `drought_led_annual` — `led_{rcp26,rcp60}_processed.nc`

### process_let_cyclone.py

Processes `let` tropical-cyclone exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. Unlike `led`, `let` is a **continuous fraction** [0,1) per cell (area of the cell exposed to TCs), so the decadal **mean** is the expected annual exposed-area fraction. Because the ensemble is thin (1 impact model `ke-tg-meanfield` × 4 GCMs = 4 members), it applies **5×5 exponential-decay spatial smoothing** (L=0.7, cos(lat)-scaled, non-NaN-normalized) to each per-member decadal map, and a **two-tier percentile** for the zero-inflated field (zeros→1; non-zeros ranked vs the non-zero 2020s baseline → [2,100]). See WORKFLOW-ISSUES.md 2026-07-24.

```bash
python scripts/process_let_cyclone.py
```

**Input**: `s3://…/TCFD/raw/isimip/cyclone_let_annual/*_2006_2099.nc4`
  (raw prefix renamed from `drought_let_annual`, which held cyclone data)
**Output**: published layer `cyclone_let_annual` — `let_{rcp26,rcp60}_processed.nc`

### process_burntarea_fire.py

Processes `burntarea` (wildfire burnt-area fraction, ISIMIP2b `biomes`) into the TCFD 6-value-class format. This is the **direct biophysical fire signal** — the fraction of each cell burned per year, in **%** — and is distinct from the Lange 2020 `lew` *exposure* family. Ensemble = 3 annual fire models (`lpj-guess`, `lpjml`, `mc2-usfs`) × 4 GCMs × {rcp26, rcp60, rcp85} = 12 members/scenario. Key points, all value-checked (see WORKFLOW-ISSUES.md 2026-07-24):
- All 3 models report **%** [0,100] → **no normalization** (equal-weight "model democracy" in raw %); inter-member spread = CI. (Contrast the water-index TWS, which needed normalization.)
- Per-member metadata **diverges**: `mc2-usfs` uses `days since 1661` (others `years since`); `lpj-guess` mislabels its `long_name` as "Fire Return Interval" (data are burnt %) and floors at 0.1%; `mc2-usfs` is a zero-inflated, ~5–7× hotter outlier. The time parser handles days/months/years per calendar.
- `trend` is a **baseline-anchored rate**: `(median[decade] − median[2020s]) / elapsed decades` (% decade⁻¹), so each panel is the trend *from the 2020s baseline to that decade* — spatially coherent (∝ the change map, corr = 1.0 at 2090s), unlike a within-decade annual slope (fire is too noisy year-to-year). No spatial smoothing (12-member ensemble is thick).

```bash
python scripts/process_burntarea_fire.py
```

**Input**: `s3://…/TCFD/raw/isimip/wildfire_burntarea_annual/*_2006_2099.nc4`
**Output**: published layer `wildfire_burntarea_annual` — `burntarea_{rcp26,rcp60,rcp85}_processed.nc`

### process_csoil_soilcarbon.py

Processes `csoil-total` (soil organic carbon stock, ISIMIP3b `biomes`) into the TCFD 6-value-class format — the direct **subsurface carbon-storage** signal (distinct from the vegetation pools `cveg`/`croot`/`cvegbg` and the net-sink flux `nbp`). Ensemble = 3 models (`classic`, `jules-es-vn6p3`, `mc2-usfs`) × their CMIP6 GCMs (2 + 5 + 5) × {ssp126, ssp370, ssp585} = 12 members/scenario. Value-checked 2026-07-25 (see WORKFLOW-ISSUES.md):
- All 3 report **kg C m⁻²** with **comparable magnitudes** (2020s medians ~5.8/7.7/10.3) → **no normalization** (equal-weight "model democracy"); inter-member spread = CI. Layer starts at the 2020s baseline (ISIMIP3b csoil begins 2015 — no full 2010s).
- **Direction is `higher_is_better`** (stored carbon is an asset; the risk is **loss**) → percentile **inverted** (low stock / decline → high risk; red = carbon loss in maps).
- **Mixed CO₂**: `jules-es-vn6p3` publishes only its fixed-2015-CO₂ run for csoil, so its trend is muted; `classic`/`mc2-usfs` are transient. All hold land use fixed. Retained as 12 members per user decision. Baseline-anchored trend (kg C m⁻² decade⁻¹); no spatial smoothing (thick ensemble).

```bash
python scripts/process_csoil_soilcarbon.py
```

**Input**: `s3://…/TCFD/raw/isimip/soilcarbon_csoil_annual/*_csoil-total_global_annual_2015_2100.nc`
**Output**: published layer `soilcarbon_csoil_annual` — `csoil_{ssp126,ssp370,ssp585}_processed.nc`

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
