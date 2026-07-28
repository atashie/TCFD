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

## GUARDRAILS §8 — never guess specifier codes

Specifiers are a controlled vocabulary. **`count=0` for something that plausibly exists is
a red flag, not a conclusion** (tropical-cyclone exposure is `let`, not the mnemonic
`letc`). `count=1001` is the API maximum — the result set is truncated and must not be
generalized from.

- **Enumerate families, don't guess members.** The Lange 2020 exposure family is exactly
  six: `led` drought, `leh` heatwave, `lew` wildfire, `ler` river flood, `lec` crop
  failure, `let` tropical cyclone.
- **A variable can have more than one representation.** "Wildfire" is both the `lew`
  exposure member and the direct `burntarea` burnt-area fraction. Enumerate all products
  and present the trade-offs.
- **The same quantity can be named differently per round** — ISIMIP2b `csoil` vs ISIMIP3b
  `csoil-total`. Searching the 3b name in 2b returns nothing, which reads as "absent".
- **Watch for cross-sector duplicates.** `elm-eca` publishes csoil-total under both
  `biomes` and `permafrost`, byte-identical (same sha512). Ingesting both double-weights
  the model. Compare checksums.
- When processing one member of a family, record the **whole** family in the catalog then.

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
