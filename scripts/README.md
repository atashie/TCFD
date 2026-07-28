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
- Trend panels compare the **2030s vs 2090s** (the 2020s trend is identically 0, since trends are anchored to the 2020s baseline)
- **`contact_sheet.html`** — per-member panels at full 0.5° (built by the processor via `utils/contact_sheet.py`, adopted through `finalize_layer(extra_maps=...)`, linked at the top of `index.html`). Embedded base64 PNGs, one pixel per grid cell, so 17 members cost ~1.8 MB instead of ~4.4M JSON numbers.
- Colourbars carry **no title** — the quantity/units label is a `<sub>` line under each figure title, so the colourbar stays as narrow as its tick labels instead of stealing width from the map.
- **`maps_bundle.zip`** — the whole collection in one object (~9.5 MB, 21 files), because the S3 console downloads one file at a time and the ~20 HTML pages are interlinked. Unzip and open `index.html`; all links are relative so navigation still works.
- Values are serialized at **5 significant figures** (`_compact()`) — a ~35% payload cut at a max relative error of 1.3e-5, far below what a colour scale can render. Display only; the NetCDF in `data/` keeps full precision. The collection is ~57 MB uncompressed / ~9.5 MB zipped (the `go.Heatmap` z-grid carries ocean nulls the old scatter omitted, so the raster is slightly larger but pixel-exact).

### generate_qa_report.py

Runs the invariant checks the 6-value-class contract depends on, plus per-scenario statistics, and publishes `qa_report.json` + `qa_report.html`. **Layer-generic**: checks are driven by each file's own declared attributes (`trend_definition`, `percentile_direction`, `baseline_source`), so the same tool serves every annualized layer.

```bash
python scripts/generate_qa_report.py soilcarbon_csoil_annual
python scripts/generate_qa_report.py drought_led_annual --version v2026-07-24_abc1234
python scripts/generate_qa_report.py cyclone_let_annual --local-only
```

**Input**: a published layer in S3 (gate-verified before any check runs)
**Output**: `…/layers/{layer_id}/{version}/qa/` (ungated/regenerable). Exit code is non-zero if any check fails.

Checks: all value classes present and correctly shaped; `lower_ci ≤ median ≤ upper_ci`; zero-width CIs isolated to all-zero or single-model cells; percentile within [1,100]; percentile orientation matches the declared `percentile_direction`; `trend == 0` in the baseline decade; `trend × elapsed_decades == change map` (GUARDRAILS §10); shared baseline bit-identical across scenarios; `n_members`/`n_models` consistent with the finite data; land coverage non-empty. An unrecognized `trend_definition` is reported as a **skipped** check rather than silently dropped.

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

Processes `burntarea-total` (wildfire burnt-area fraction, **ISIMIP3b** `biomes`) into the TCFD 6-value-class format. This is the **direct biophysical fire signal** — the fraction of each cell burned per year, in **%** — and is distinct from the Lange 2020 `lew` *exposure* family and the `ffire` emissions flux. **Moved from ISIMIP2b/RCP to ISIMIP3b/SSP on 2026-07-28** (newer round and scenario family win where the newer data is viable); the 2b build is retired. Ensemble = 12 members/scenario:

| model | cadence | GCMs | soc/sens | effective grid |
|---|---|---|---|---|
| `mc2-usfs` | annual | 5 | `nat/default` | 0.5° clean |
| `visit` | monthly | 5 | `2015soc/default` | 0.5° clean |
| `classic` | monthly | 2 | `2015soc/default` | **1.0°** (kept deliberately) |

Key points, all value-checked (see WORKFLOW-ISSUES.md 2026-07-28):
- **Monthly members annualize by SUM, not mean** — burnt area *accumulates* over its reporting interval. Verified against `classic`, which publishes the same run daily *and* monthly: each monthly value equals the sum of that month's daily values (1e-6, r = 1.0). A mean would under-scale the layer **12×**. Years lacking all 12 months are dropped, not partially summed.
- **Annual values legitimately exceed 100%** where a cell reburns (decadal means to ~107%), so the CI is floored at 0 and **unbounded above** — clamping to 100 would push `upper_ci` below the median.
- All 3 models report **%** and agree within **1.6×** → **no normalization** (equal-weight "model democracy" in raw %); inter-member spread = CI. But note the declared unit is *not* sufficient grounds: 2b `clm45`/`orchidee` declare `%` on a 0–1 fraction scale, and `classic`'s `2015soc-from-histsoc` variant is fraction-scaled for one GCM and percent for the other — hence the `2015soc` pin here, **not** the `2015soc-from-histsoc` that `csoil` uses for the same model.
- `elm-eca` is **excluded** (~4°×5°); `classic` at an effective 1.0° is retained. Emits `n_members`/`n_models` per cell — the 3 models do not share a land mask.
- `trend` is a **baseline-anchored rate**: `(median[decade] − median[2020s]) / elapsed decades` (% decade⁻¹), so each panel is the trend *from the 2020s baseline to that decade* — spatially coherent (∝ the change map), unlike a within-decade annual slope (fire is too noisy year-to-year). No spatial smoothing (12-member ensemble is thick). ISIMIP3b starts 2015 → the layer begins at the 2020s baseline.

```bash
python scripts/process_burntarea_fire.py
```

**Input**: `s3://…/TCFD/raw/isimip/wildfire_burntarea_annual/*_burntarea-total_global_*_2015_2100.nc`
**Output**: published layer `wildfire_burntarea_annual` — `burntarea_{ssp126,ssp370,ssp585}_processed.nc`

### process_csoil_soilcarbon.py

Processes `csoil-total` (soil organic carbon stock, ISIMIP3b `biomes`) into the TCFD 6-value-class format — the direct **subsurface carbon-storage** signal (distinct from the vegetation pools `cveg`/`croot`/`cvegbg` and the net-sink flux `nbp`). Ensemble = 4 of the 5 models that publish it (`classic`, `jules-es-vn6p3`, `mc2-usfs`, `visit`) × their CMIP6 GCMs (2 + 5 + 5 + 5) × {ssp126, ssp370, ssp585} = **17 members/scenario**. **`elm-eca` is excluded** — declares 0.5° but is effectively ~4°×5° (see `EXCLUDED_MODELS`). Note ISIMIP**2b** names this variable bare `csoil`; only **3b** uses `csoil-total`. Fully value-checked 2026-07-27 (see WORKFLOW-ISSUES.md):
- **Mixed cadence**: `visit` publishes csoil-total **only monthly** and is annualized by the **mean of each year's 12 months** — verified immaterial (within-year CV 0.108%; |Dec − annual mean|/mean 0.103%).
- All 4 report **kg C m⁻²** with medians within **1.8×** → **no normalization** (equal-weight "model democracy"); inter-member spread = CI. Caveat: members are pooled with a *mean* and `visit` carries a much fatter (tropical peat) tail, so pooling is more tail-sensitive than the median agreement suggests. `classic` is natively **1°**, 2/17 = 12% of the weight. Layer starts at the 2020s baseline (ISIMIP3b csoil begins 2015 — no full 2010s).
- **Direction is `higher_is_better`** (stored carbon is an asset; the risk is **loss**) → percentile **inverted** (low stock / decline → high risk; red = carbon loss in maps).
- **Mixed CO₂**: `jules-es-vn6p3` publishes only its fixed-2015-CO₂ run. Contrary to this file's earlier note, that does **not** mute its trend — it gives `jules` the **strongest loss** of the five (−4.4% at ssp585) while the four transient models run flat-to-positive. Baseline-anchored trend (kg C m⁻² decade⁻¹); no spatial smoothing (thick ensemble).
- **Also emits `n_members` / `n_models`** per cell: the 4 models do **not** share a land mask, so ~81% of land carries all 17 members and 9.2% is single-model. Those cells are **kept, not masked**; the counts let consumers filter.
- Years are decoded with **cftime**, not days/365 arithmetic (which misassigns December of a monthly member).

```bash
python scripts/process_csoil_soilcarbon.py
```

**Input**: `s3://…/TCFD/raw/isimip/soilcarbon_csoil_annual/*_csoil-total_global_*_2015_2100.nc` (66 files staged, 7.47 GB, both `annual` and `monthly`; 51 used after excluding `elm-eca`)
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
