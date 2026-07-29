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

Checks: all value classes present and correctly shaped; `lower_ci ≤ median ≤ upper_ci`; zero-width CIs isolated to all-zero or single-model cells; percentile within [1,100]; percentile orientation matches the declared `percentile_direction`; `trend == 0` in the baseline decade; `trend × elapsed_decades == change map` (GUARDRAILS §10); shared baseline bit-identical across scenarios; `n_members`/`n_models` consistent with the finite data; land coverage non-empty; effective spatial resolution (seam spacing + inside-vs-seam gradient ratio, GUARDRAILS §11). An unrecognized `trend_definition` is reported as a **skipped** check rather than silently dropped.

Two details worth knowing:

- **The percentile-orientation check uses Spearman, not Pearson.** `percentile` is a percentile-of-score — a monotone but deliberately **non-linear** rank transform of the value — so on a heavy-tailed, zero-inflated variable Pearson reads far below 1 even when the ordering is perfect. It scored `burntarea` at +0.53 and failed a correct layer; Spearman gives +1.0000 and still catches a genuinely inverted or scrambled percentile.
- **An undeclared sharp latitude seam raises a warning.** A band profile cannot see a step that falls on a band *edge* — the 60.0°N discontinuity in `fldfrc` sits exactly between the `45..60` and `60..70` bands and was invisible to it. The check scans **row-to-row** jumps in the zonal mean and flags the largest when it is ≥9× the **local** median change (±10 rows) *and* a ≥1.5-fold level step, ignoring rows with <150 finite cells. Declare a confirmed seam in `known_latitude_seams` and it passes with the evidence still printed. Calibrated against five other layers: real seam 11.4–11.8×, loudest non-seam 6.9×. Three earlier formulations were rejected — a MAD z-score straddled a fixed cut between scenarios (12.2/11.7/11.7 for the *same* seam), a global-median denominator flagged steep gradients as loudly as seams, and without the row-count floor it fired on 12–32-cell polar rows of three existing layers. See GUARDRAILS §14.
- **A zonal profile (land-mean by latitude band) is reported for every layer**, in both JSON and HTML. An aggregate or area-weighted statistic is structurally blind to a defect confined to one latitude zone, because polar cells carry almost no area. A layer that declares `zonal_expectation: low_latitude_dominated` additionally **warns** when a polar band exceeds its tropical band — which is how `visit`'s inverted fire profile became visible (GUARDRAILS §9).

### process_qg.py

Processes raw groundwater runoff (qg) data from ISIMIP.

```bash
python scripts/process_qg.py
```

**Input**: `data/raw/*.nc` (raw NetCDF files from ISIMIP)
**Output**: `data/processed/qg_*.nc` (processed files by scenario)

### download_driedarea_drought.py

Ingests ISIMIP3b `driedarea` (Heinicke2026 drought exposure) into S3 raw. 45 files = 3 GHMs (`h08`, `jules-w2`, `watergap2-2e`) × 5 CMIP6 GCMs × ssp126/370/585, ~170 MB total. Each file is sha512-verified against its ISIMIP `.json` sidecar before upload, and the run is resumable — a cached copy matching the sidecar is reused, and members already in S3 are skipped. `historical` is deliberately not ingested (the 2020s baseline lives inside the future files).

```bash
python scripts/download_driedarea_drought.py [--dry-run] [--force]
```

**Output**: `s3://…/TCFD/raw/isimip/drought_driedarea_annual/`

### process_driedarea_drought.py

Processes `driedarea` into the TCFD 6-value-class format — **the current drought layer**, superseding `led` (newer round, SSP scenario family). Like `led` it is a **binary** per-cell annual flag, so the decadal statistic is the **mean** (drought frequency), never the median — and it is an *occurrence frequency*, not a within-cell area share.

Three things differ from `process_led_drought.py` and must not be copied across:
- **Filename fields `[0]`/`[1]`** hold model/GCM, not `[1]`/`[2]` — the `led` parser mis-keys every member. Membership is asserted (duplicates, vocabularies, exact per-scenario count) because the load pattern silently overwrites a duplicate.
- **Union land mask, nothing filtered.** The three GHMs disagree (union 63,455 / intersection 46,024 cells); partially covered cells are kept, with per-cell `n_members`/`n_models` emitted.
- **Single-tier percentile by decision, not auto-selected.** The baseline's exact-zero mass is 3.59% over the union but 0.18% over fully-covered cells, so burntarea's >2% two-tier trigger would fire on a coverage artefact. Baseline-anchored trend (masked to NaN off-land), no smoothing.

```bash
python scripts/process_driedarea_drought.py
```

**Input**: `s3://…/TCFD/raw/isimip/drought_driedarea_annual/*_driedarea_global_annual_2015_2100.nc`
**Output**: published layer `drought_driedarea_annual` — `driedarea_{ssp126,ssp370,ssp585}_processed.nc`

### process_led_drought.py

Superseded by `process_driedarea_drought.py` (ISIMIP3b/SSP); retained for the ISIMIP2b/RCP build. Processes `led` drought exposure (Lange 2020, ISIMIP2b) into the TCFD 6-value-class format. `led` is a **binary** per-cell annual flag, so the decadal statistic is the **mean** (drought frequency), not the median. Handles the `years since 1661` / `360_day` calendar, `(model, gcm)` ensemble members, shared 2020s baseline, and zero-inflation-aware percentile/CI. See WORKFLOW-ISSUES.md 2026-07-24.

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

### process_timber_tempnle.py

Processes the ISIMIP2b **temperate needleleaf evergreen** conifer PFT into two TCFD layers from **one parameterized script** — `timber_cveg-tempnle_annual` (**kg m⁻²**) and `timber_npp-tempnle_annual` (**kg m⁻² yr⁻¹**, converted from `kg m-2 s-1` on the members' own 365-day calendar). The timber-growth-productivity signal for maritime/temperate conifers (the Sitka-spruce use case). Fully value-checked 2026-07-28 (see WORKFLOW-ISSUES.md and GUARDRAILS §13):
- **Why 2b, not 3b**: there is **no PFT-resolved `cwood` anywhere in ISIMIP2b** (all 11 `biomes` models *and* `permafrost` checked), and **no climate-zone-resolved conifer PFT carries `cwood` in any round**. 3b has `cwood` only for the all-climate classes (`classic` `evgndltr` 2 GCMs, `jules-es-vn6p3` `ndlevg` 5), which pool boreal with temperate conifers — **34.6% / 29.3%** of their global PFT wood sits at 55–70°N. Climate specificity was chosen over the wood-only pool; wood is p50 **77–90%** of conifer `cveg`, so `cveg` is a damped wood surrogate, not a different signal.
- **Ensemble** = `orchidee` `tendev` (4 GCMs) + `orchidee-dgvm` `tendev` (2, no rcp85) + `clm45` `needleleaf-evergreen-tree-temperate` (ragged: 2/1+1 fixed-CO₂/2) + `lpjml` `temperate-needleleaved-evergreen-tree` (4, **npp track only** — no per-PFT `cveg`) → **cveg 8/8/6, npp 12/12/10** at rcp26/60/85. rcp85 is thinner by construction.
- **`caraib` `ndevteclt` is excluded** (~1–2° effective grid, per-gridcell convention, no rcp85) but stays ingested as a sensitivity panel. Its codes carry **no `long_name`** and were identified by cover-weighted biogeography — `ndevteclt` is temperate; **`ndevtecdt` is boreal despite "te"** in the code.
- **Two value conventions harmonized to per-tile density**: `orchidee`/`orchidee-dgvm`/`clm45` report on the PFT's own tile area, `lpjml`/`caraib` are already cover-scaled and are divided by cover with a **1% floor** (floored cells drop out of `n_members`). Per-tile was required because `clm45` publishes **no cover fraction at all**. Harmonizing collapsed the apparent spread from **10.5×/177×** to **2.35×/1.83×** on a common mask → **no normalization**.
- **Pooled by model FAMILY, not member** (`orchidee` + `orchidee-dgvm` are one model in two configurations), with the CI as ±1 SD **across family means** — a deliberate departure from the flat inter-member SD used by led/let/burntarea/csoil, since `orchidee` supplies 6 of 8 cveg members.
- **Mask rule: ≥2 model families.** The pooled level differs **12×** between multi-model and single-model cells (marginal habitat), so the union would print hard mask edges and distort the percentile. Note this degenerates to "all models" for cveg, which has only 2 families.
- **Baseline composition is per scenario.** Composition varies by scenario, so an all-member shared baseline manufactured a **−0.72 kg m⁻² dec⁻¹** rcp85 loss with a fake recovery. The 2020s panel is therefore *not* bit-identical across scenarios; `members_by_scenario` (member **identity**) is declared so the QA check groups by composition.
- **Direction is `higher_is_better`** (risk = declining growth potential); baseline-anchored trend; npp CI **not floored** (npp is legitimately negative — net carbon loss). Also emits `pft_cover` as a context field (incomplete: `clm45` publishes none).

```bash
python scripts/process_timber_tempnle.py cveg
python scripts/process_timber_tempnle.py npp
```

**Input**: `s3://…/TCFD/raw/isimip/timber_{cveg,npp}-tempnle_annual/` (54 + 78 objects, 1.22 GB incl. `pft-` cover fractions)
**Output**: published layers `timber_cveg-tempnle_annual` / `timber_npp-tempnle_annual` — `{cveg,npp}-tempnle_{rcp26,rcp60,rcp85}_processed.nc`

### download_fldfrc_flood.py

Ingests ISIMIP2b `Zimmer2023` `fldfrc` (CaMa-Flood annual flooded-area fraction) — **216 files** = 6 GHMs × 4 GCMs × {rcp26, rcp60, rcp85} × {none, 100yr, flopros} — into **three** raw prefixes, one per protection level.

**Coarsens at ingest, deliberately.** The source is 150 arcsec (4320 × 8640) and the full set is ~54 GB against 19 GB of local scratch, so each member is streamed: downloaded, sha512-verified against its ISIMIP sidecar, coarsened **12×12 to 0.5°**, then the original deleted. What lands in raw is the 0.5° field, ~36× smaller. That is a departure from "raw is byte-for-byte what ISIMIP served", so it is made auditable rather than silent — every file records the `source_url`, the **sha512 of the original**, and the exact transform, in both its own global attrs and `layer.json`. Raw staging is transient by contract anyway (STORAGE.md — `cleanup_raw` deletes it once `source_url` + checksum are recorded), and `files.isimip.org` is not behind the Anubis anti-bot that guards the API, so re-fetching an original is routine.

- The aggregation is a **block sum, not interpolation**: 4320/12 = 360 and 8640/12 = 720 exactly, and each block centre coincides with the ISIMIP 0.5° centre to 1.4e-14°. Sub-cells are weighted by true spherical area (they differ across a block's 12 latitude rows), and the run **asserts area conservation** per file (observed rel. err ≤ 1e-15).
- The denominator is the **full 0.5° cell area**, not the valid sub-cells. So `sum(value × cell_area)` over any region is exactly the flooded km² — the field is directly aggregatable. A valid-only denominator would report "fraction of the modelled floodplain that flooded" and inflate partly-covered cells.
- Emits `floodplain_fraction` (area share of each cell inside the CaMa-Flood domain) so a partly-covered cell is auditable rather than merely looking low.
- Resumable (skips members already in S3 raw) and parallel; `--protection`, `--workers`, `--limit`, `--dry-run`, `--force`.

```bash
python scripts/download_fldfrc_flood.py                      # all 216
python scripts/download_fldfrc_flood.py --protection flopros --workers 6
python scripts/download_fldfrc_flood.py --protection 100yr --limit 2   # smoke test first
```

**Input**: `https://files.isimip.org/ISIMIP2b/DerivedOutputData/Zimmer2023/{MODEL}/{gcm}/future/`
**Output**: `s3://…/TCFD/raw/isimip/riverflood_fldfrc-{none,100yr,flopros}_annual/` (~3.2 MB per member; ~700 MB total)

### process_fldfrc_flood.py

Processes `fldfrc` into the TCFD 6-value-class format as **three parallel layers, one per flood-protection level** — 24 members/scenario (6 GHMs × 4 GCMs), rcp26/60/85. This is the only ISIMIP flood product that is a genuine **area** share; see `config/isimip_search_catalog.yaml` → `search_results.flooding` for the options review that rejected the alternatives.

**The protection level is the layer's biggest choice, not a parameter.** Same member, global flooded area, 2020s → 2090s: `none` 6.10M → 6.28M km² (**+2.9%** — it counts routine seasonal floodplain wetting that recurs every year, so the metric *saturates*), `100yr` 261k → 502k (**+92.4%**), `flopros` 744k → 1,130k (**+51.9%**). Protectiveness runs `none < flopros < 100yr`. **`flopros` is not a proxy for "unprotected"** — measured, not assumed: it retains only 19–36% of the undefended signal even in the least-defended regions (Bangladesh 0.189, Niger 0.364) and 0.8–4% in well-defended ones (Netherlands 0.008, Mississippi 0.015). It *is* spatially real, being less protective than a uniform 100yr in Niger/Bangladesh and more so in the Mississippi/Netherlands. Read `100yr` as a **severity threshold** that needs no defense database.

- Decadal **mean**, never median — an exposed-area frequency measure (the Lange 2020 rule), and at 0.5° a single year is ~93% exact zeros, so a decadal median would be 0 nearly everywhere.
- `trend` is **baseline-anchored** (% of cell area decade⁻¹): annual flooded area swings ~**17×** between adjacent decades (h08/miroc5/flopros: 2080 = 1,207,332 km² vs 2100 = 69,751 km²), so a within-decade annual slope is noise (GUARDRAILS §10).
- Two-tier percentile (zero-inflated); `higher_is_worse`; **no normalization** (one hydrodynamic model, one unit — spread is genuine GHM+GCM uncertainty); **no spatial smoothing** (24 members is thick).
- CI clamped to **[0, 100] %**. The upper clamp is safe *here* and is not the `burntarea` mistake: flooded fraction is a **bounded** share of cell area, so `median ≤ 100` and `min(median+sd, 100) ≥ median`. Burnt area is cumulative and legitimately exceeds 100%, which is why that layer leaves `upper_ci` unbounded.
- Coverage **62,066 cells = 128.8 M km² = 95.7% of land excluding Antarctica** (Lange2020 `led`: 67,420; 61,546 shared). A normal ISIMIP land mask, not a sparse floodplain subset — 79.6% of domain cells are ≥99% inside the domain, median 100%. Cells outside stay NaN; nothing is zero-filled. Logs whether member domains differ and takes the union.

```bash
python scripts/process_fldfrc_flood.py                        # all three layers
python scripts/process_fldfrc_flood.py --protection flopros
python scripts/process_fldfrc_flood.py --protection none --no-publish
```

**Input**: `s3://…/TCFD/raw/isimip/riverflood_fldfrc-{protection}_annual/*_fldfrc_*_halfdeg_global_annual_*.nc`
**Output**: published layers `riverflood_fldfrc-{none,100yr,flopros}_annual` — `fldfrc_{rcp26,rcp60,rcp85}_processed.nc`

### process_fldfrc_event100yr.py

The **fourth** flood layer, and the only one that is not an expected-annual-area. The three protection layers answer "how much area floods in an average year, given protection P"; this one answers the two questions a 1-in-100-year flood product needs: **how often** is the (preindustrial) 1-in-100 flood exceeded, and **how much area** does it cover when it is.

**It is not `none` + `100yr`.** That sum double-counts: `100yr` is a **subset** of `none`, not its complement. Measured per cell per year at native 150 arcsec — `100yr > none` in **0.000%** of cells; where `100yr > 0` it **equals** `none` in **84–88%** (mean ratio 0.90–0.93); and `none + 100yr` would exceed **100% of a cell** in 6,867–28,139 cells, impossible for an area fraction. `100yr` is a *filtered copy* of `none`, kept only in years that overtopped a 1-in-100 defence. The sum is numerically unstable too — 17% *below* the measured footprint in the 2020s, drifting with scenario because the terms scale differently (`none` +5%, `100yr` +226%).

Correct construction, per member and decade, from the paired annual series:
`exceedance_frequency = (years with 100yr > 0) / valid years` and `event_footprint = mean of none over exactly those years`. A year counts as valid only where both fields are present, so numerator and denominator share one year set.

- **`median` carries the FREQUENCY (% of years), not the footprint.** The footprint grows only **+1.8%** from the 2020s to the 2090s at rcp85 — topography fixes the extent of a given-magnitude flood — while the frequency **more than doubles (+122%)**: the preindustrial 1-in-100 flood becomes a **1-in-5.7-year** flood. Making the footprint primary would ship a confidently flat layer, the trap the `none` level falls into. The footprint rides along as `event_footprint` with its own CI.
- This also explains the three siblings: `none` saturates (+5.1%) and `100yr` explodes (+225.5%) because **all four are measuring a frequency change**; only this layer states it directly.
- **Two caveats in `known_issues`**: "1-in-100" is **preindustrial** (GEV fit to `picontrol`, `fit_soc=1860soc`), so the 2020s already sit at ~1-in-15 years; and the conditional mean is **selection-biased upward** because cells that never exceed within a decade drop out — read the 2020s→2090s *ratio* as robust and the level as indicative.
- `event_footprint` is undefined where a member never exceeds in a decade, so its ensemble depth is thinner and more variable than the frequency's (median 13 of 24 members). Qualify it with **`n_members_footprint`**, not `n_members`.
- Inputs come from **two sibling raw prefixes**; this layer's own raw prefix stays empty by design, and `raw_entries` records all 144 inputs so `layer.json` still carries the chain back to the 150 arcsec originals.

```bash
python scripts/process_fldfrc_event100yr.py [--no-publish]
```

**Input**: `s3://…/TCFD/raw/isimip/riverflood_fldfrc-{none,100yr}_annual/*_halfdeg_global_annual_*.nc` (matched on model × GCM × scenario)
**Output**: published layer `riverflood_fldfrc-event100yr_annual` — `fldfrc_{rcp26,rcp60,rcp85}_processed.nc`

## Dependencies

These scripts require the same dependencies as the `isimip-pipeline` package:
- xarray
- plotly
- numpy
- pandas

Install with: `pip install -e isimip-pipeline/`
