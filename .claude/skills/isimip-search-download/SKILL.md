---
name: isimip-search-download
description: Search the ISIMIP repository for climate/impact variables and download members into S3 raw staging. Use when looking for a new variable, enumerating which models/GCMs/scenarios exist, checking data availability before processing, or ingesting raw ISIMIP files.
---

# Search & Download (ISIMIP)

## Read the catalog first

Before any search, consult `config/isimip_search_catalog.yaml` — it records prior
empirical findings (variable definitions, units, PFT codes, which models and forcings
have data, coverage, scenario availability). It avoids redundant queries and contextualizes
new results.

**But the catalog is opportunistic, not authoritative for coverage.** Its per-model
`variables:` lists were compiled for earlier searches. It has been wrong more than once:
it under-listed ISIMIP2b biomes models (5 documented vs 11 real), and later still missed
`CLM45` and `VEGAS` for `csoil`. **Enumerate the file server; treat the catalog as a
starting hypothesis.** When you find something new, write it back with the search date.

## Enumerating — the mechanics that actually work

**The API at `data.isimip.org` is behind Anubis anti-bot** and returns an "Access Denied"
challenge to WebFetch/curl. Use the file server instead:

```
https://files.isimip.org/{round}/OutputData/{sector}/{MODEL}/{gcm}/{future|historical}/
```

- **Use `curl … | grep -oE`, not WebFetch.** WebFetch's summarizer silently truncates long
  autoindex listings — it repeatedly claimed "only one file" for a many-file `csoil-total`
  directory. Single directories here can hold 3,000+ files.
- **Harvest SERIALLY with retries, and assert each listing is non-empty.** Parallel curls
  get rate-limited and return **empty** listings that are indistinguishable from a genuine
  "no data" answer. This has produced wrong availability matrices: 13 of 22 directories
  came back empty on one parallel pass. An empty listing is a *failure signal*.
- Confirm a genuinely absent path with its HTTP status (a real 404) rather than an empty body.
- Each file has a `.json` sidecar with `size`, `sha512` and parsed `specifiers` — use it to
  verify downloads and to read the soc/sens tokens without downloading the NetCDF.
- To inspect a header without downloading a large file, open it lazily over HTTP
  (`fsspec` + `xarray`/`h5netcdf`) instead of pulling hundreds of MB.

Filename grammar (parse from the END; model/GCM names contain hyphens, never underscores):

```
{model}_{gcm}_{forcing}_{scenario}_{soc}_{sens}_{variable}_{region}_{timestep}_{y0}_{y1}.nc
                                                     [-5]     [-4]     [-3]   [-2] [-1]
```

**Parsing from the END is not optional for PFT-qualified variables.** The variable field
itself contains hyphens (`cveg-needleleaf-evergreen-tree-temperate`), so any fixed forward
index breaks on exactly the datasets that need it most.

**Anchor the extension in harvest regexes: `\.nc4?$`.** An unanchored `[\w.-]+\.nc` matches
the `.nc` *inside* `.nc4` and truncates every ISIMIP2b filename by one character, so every
URL built from the harvest 404s — which reads as "no data" rather than as a parsing bug.
(See also the round/extension note under §8.) 2b files do carry `.json` sidecars with
`size` + `sha512`; 3b `biomes` largely does not, so verify 3b by HTTP status instead.

## GUARDRAILS §8 — never guess specifier codes

Specifiers are a controlled vocabulary. **`count=0` for something that plausibly exists is
a red flag, not a conclusion** (tropical-cyclone exposure is `let`, not the mnemonic
`letc`). `count=1001` is the API maximum — the result set is truncated and must not be
generalized from.

- **Enumerate families, don't guess members.** The Lange 2020 exposure family is **twelve,
  not six** (corrected 2026-07-28): `le{d,r,w,c,h,t}` = **land-area** fraction exposed,
  each paired with a `pe{d,r,w,c,h,t}` = **population** fraction exposed twin. Hazards are
  `d` drought, `r` river flood, `w` wildfire, `c` crop failure, `h` heatwave, `t` tropical
  cyclone. Documenting only the `le*` half and calling the family "six" stood uncorrected
  for days. Also: `lew` has **no rcp85** (rcp26/rcp60 only).
- **List `DerivedOutputData/` before concluding a product has no newer-round version.**
  A product can be re-issued in a later round under a **different publication directory
  with different variable names**, so searching the old variable name returns nothing and
  reads as "absent". Lange 2020's exposure concept *was* re-issued for ISIMIP3b, split
  across `Heinicke2026` (`driedarea`, `floodedarea`) and `Zantout2025` (`heatwave`,
  `wildfire`, `cropfailure`) — hazard words, not `le*` codes. The catalog asserted
  "ISIMIP3b/SSP version of this family NOT found (0 hits)" on that basis, and the drought
  layer nearly shipped on rcp26/rcp60 when ssp126/370/585 existed. The listing is 2–4
  entries and costs one call:

  ```bash
  curl -s https://files.isimip.org/ISIMIP3b/DerivedOutputData/ | grep -oE 'href="[^"]+/"'
  ```
- **A variable can have more than one representation.** "Wildfire" is the `lew` exposure
  member, the `ffire` emissions flux, the ISIMIP3a-only `fire`-sector diagnostics
  (`firesize`/`firenr`/`fireints`/…), *and* the direct `burntarea` burnt-area fraction.
  Enumerate all products and present the trade-offs.
- **If a catalog entry names a sector, verify that sector was actually walked.** The
  wildfire section was titled "fire sector + biomes burntarea" but only `biomes` had ever
  been enumerated. In ISIMIP3b that was harmless (the `fire` sector's burntarea files are a
  strict subset of `biomes`), but ISIMIP3a's `fire` sector holds **10 models found nowhere
  else**. A section heading is not evidence of coverage.
- **File extensions differ by round: ISIMIP2b publishes `.nc4`, ISIMIP3a/3b publish `.nc`.**
  A `*.nc` filter silently drops the **entire 2b round** and reads as "no data". This
  produced a false negative twice in one session — once on a Lange2020 listing, once on a
  `burntarea` inventory. Match both extensions.
- **The same quantity can be named differently per round** — ISIMIP2b `csoil` vs ISIMIP3b
  `csoil-total`. Searching the 3b name in 2b returns nothing, which reads as "absent".
- **Watch for cross-sector duplicates.** `elm-eca` publishes csoil-total under both
  `biomes` and `permafrost`, byte-identical (same sha512). Ingesting both double-weights
  the model. Compare checksums.
- When processing one member of a family, record the **whole** family in the catalog then.

## Vegetation / PFT searches — treat variable × PFT × round as one lookup

A PFT-resolved variable exists per **(round, model, variable, PFT)**, and the combination
that a product needs frequently does not exist at all. Establish the intersection before
promising anything (all verified 2026-07-28):

- **PFT-resolved `cwood` exists ONLY in ISIMIP3b, from two models** — `classic` `evgndltr`
  (2 GCMs) and `jules-es-vn6p3` `ndlevg` (5 GCMs). There is **no PFT-resolved `cwood`
  anywhere in ISIMIP2b** — not in `biomes` (11 models) and not in `permafrost`. 2b PFT
  output is `cveg` / `npp` / `gpp` / `pft` only.
- **No climate-zone-resolved conifer PFT carries `cwood` in ANY round.** The climate-zone
  classes (`clm45` `needleleaf-evergreen-tree-temperate`, `lpjml`
  `temperate-needleleaved-evergreen-tree`, `orchidee` `tendev`, `caraib` `ndev*`) are 2b, so
  wanting *both* climate specificity and the wood-only pool is unsatisfiable. Surface that
  trade-off explicitly rather than picking one silently. Substituting `cveg` costs the
  root+leaf fraction: wood is p50 **77–90%** of conifer `cveg`.
- **`cveg` is not universal either** — `lpjml` publishes PFT-resolved `npp`/`gpp`/`pft` but
  **no** per-PFT `cveg`, so a cveg track and an npp track have different member counts.
- **Check the model's own PFT scheme before mapping a species onto it.** `jules` splits
  broadleaf evergreen into `bdlevgtemp`/`bdlevgtrop` but publishes a single `ndlevg`; the
  only 3b `cwood` token containing "temp" is `cwood-bdlevgtemp`, which is a temperate
  **hardwood**, not a conifer. Reading the code as "temperate evergreen ⇒ what I want" is
  how a broadleaf class gets proposed for a spruce question.
- **A generic class is one PFT with one global parameter set**, not an internally weighted
  mixture of climate-specific sub-types. Its cost is a shared parameterisation, and the
  climate-zone mixing bites at the **percentile ranking** step (cells ranked against a
  global distribution spanning boreal to subtropical), not in the raw per-cell values.
- **Climate-zone specificity is not automatically better — it depends on the species.** It
  helped loblolly (one zone). For a hyper-oceanic species spanning zones (Sitka spruce) it
  *fragments* the range: `lpjml` temperate-NLE covers UK/Ireland (774/774 cells) but only
  609/2451 in PNW-BC-SEAK, while `lpj-guess` boreal-NLE covers PNW (1866) and just **6** UK
  cells. Check the candidate PFT's coverage **in the regions the product is for**.
- **Species-level ISIMIP is 8–9 European stand sites, and there is no Sitka spruce.** The
  forestry sector (2b: 10 models; 2a: 15) is the only place individual species and the real
  forestry metrics (`mai` mean annual increment, `vol`, `harv`) appear, but it is not
  griddable. Species vocabulary is exactly
  `{acpl, bepe, fasy, lade, piab, pipi, pisy, psme, quro}`.
- **`pft-{code}` cover fractions are published alongside** by most models and are how you
  answer "does this model actually place this plant here". `clm45` publishes none.

**The catalog's PFT guidance is written as recommendations, and an untested suggestion reads
exactly like a measured fact.** Its top recommendation for temperate conifers pointed at
MC2-USFS biome classes for years; MC2 publishes **no carbon or NPP resolved by biome type**
(only `-total`/`-tree`/`-grass`), its 47 biome names are **one-hot 0/100 presence flags**,
and `pft-maritimeevergreenneedleleafforest` is **identically zero in every run** (3b
ssp370/gfdl, 3b ssp126/ukesm, 3a obsclim — 7,008,981 finite values, all 0). Before relaying
a catalog recommendation, confirm the variable it names actually exists at the resolution it
claims — and when a recommendation turns out to be unusable, **withdraw it in the catalog**
rather than leaving it to be found again.

## Choosing what to download

- **Ask the user for temporal resolution** when both exist (§1) — monthly files are often
  ~12× larger. Do not download monthly by default. Note that monthly data for a smooth
  *stock* can be annualized cheaply, which may be the way to deepen a thin ensemble.
- **Maximize ensemble depth.** Prefer multi-model ensembles; never settle for a single
  model when more exist. More members improve the CI.
- **Prefer ISIMIP3b (CMIP6/SSP) over 2b (CMIP5/RCP)** for new layers, but check depth:
  2b is sometimes far deeper (2b `csoil` has 9 annual models vs 3b's 3), which can be worth
  a cross-round sensitivity check even when 3b is the primary.
- **Verify soc/sens tokens per model before committing** to a uniform ensemble (§9) — they
  frequently differ, and a strict filter can silently collapse the ensemble.
- **Prefer global over regional** coverage for methodological consistency.
- **ISIMIP first.** Push back on non-ISIMIP sources unless ISIMIP genuinely lacks the data
  after a thorough search.
- No direct match? List proxy variables and relevant PFTs, suggest sectors/models that may
  approximate it, and only then discuss external data with the user.

Report availability as a matrix — model × GCM × scenario, with cadence, year span, soc/sens
tokens and dataset counts — plus the temporal record and the spatial resolution.

**Report resolution as DECLARED, not verified, until you have measured it from the values.**
Every ISIMIP file reports the round's nominal grid (0.5° / 360×720), including members that
ran natively coarser and were replicated onto it — `csoil`'s `elm-eca` is effectively
~4°×5° with an identical-looking header. Reading `ds.sizes` describes the container, not
the information content. Effective resolution is measured at processing time (§9, §11 and
the `isimip-process-visualize` skill); don't promise it here.

## Ingesting into S3

Raw goes to `TCFD/raw/isimip/{layer_id}/`, where `layer_id` is the canonical
`{hazard}_{variable}_{cadence}` and **must match the processed layer's id**.

```python
entries = storage.ingest_raw(files, layer_id, source_urls={name: url, ...})
```

- **Always pass `source_urls`.** It is what later makes `storage.cleanup_raw` safe; without
  it, raw cannot be deleted, because the API being behind anti-bot means it may not be
  re-downloadable.
- Verify each download against its sidecar `size` and `sha512` before ingesting, and make
  the downloader resumable — skip a file already present and matching.
- Download into `storage.cache_path(raw_key)` so a later `storage.stage_raw()` finds the
  size-matching cached copy and skips re-downloading.
- Never hardcode an S3 key — everything comes from `isimip_pipeline/storage.py`. Never
  write data to the repo or EFS; S3 is canonical and `/tmp/tcfd-cache` is scratch.

Then hand off to the `isimip-process-visualize` skill.

## Scenarios

- **Projections** (appear in reports): SSP `ssp126/245/370/585`; RCP `rcp26/45/60/85`.
  Availability is per-variable — ISIMIP3b `biomes` has **no ssp245** at all.
- **Non-projections**: `picontrol`, `historical`. Worth downloading to strengthen the 2020s
  baseline, but they are excluded from map generation automatically.
- **Discover scenarios dynamically from the filesystem** — never hardcode the list (§3).
