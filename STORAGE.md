# Storage Layout — S3 as Canonical Store

All pipeline data lives in S3. Local disk is an **ephemeral cache**, never a source of truth.

This follows the org-wide house rule (`wildfireModeling/CLAUDE.md:30-32`):

> **Never save data locally. Always write to S3.** We work in AWS SageMaker Studio; EFS/local filesystem is for code only. Data — raw, intermediate, outputs — lives in S3.

- **Bucket**: `climate-ai-data-science-shiny-app-data` (region `us-east-2`)
- **Root prefix**: `TCFD/`
- **Encryption**: bucket default SSE-S3 (AES256) with bucket keys. SSE-C is *blocked* — never pass SSE-C parameters on a write.
- **Bucket versioning is OFF.** An overwrite is permanent loss. This is why versions live in the key (see [Versions](#versions)).
- **Bucket lifecycle configuration is off-limits.** A bucket-wide `IntelligentTiering` rule (empty filter, Days 0) already applies to every object here, including ours. We never call `PutBucketLifecycleConfiguration` — it would replace the whole config and affect other teams' prefixes. Raw-data cost control is handled in code instead (see [Raw staging](#raw-staging)).

## Principles

1. **A layer version is immutable.** Once published, its objects are never rewritten. Reprocessing creates a new version.
2. **Adding a layer is additive.** A new hazard layer creates one new prefix and appends one registry entry. It never edits an existing layer's objects.
3. **Evidence travels with the data.** A version's QA/QC reports, maps, and provenance manifest sit *inside* that version prefix, so a report can never be read against the wrong data.
4. **Products never mix.** TCFD/CDP (8 value classes, annualized) and Water Risk Index (20 value types, monthly) are separated at the top level, per `CLAUDE.md`.
5. **One place builds keys.** All key construction goes through `isimip_pipeline/storage.py`. No script hardcodes an S3 key.
6. **Nothing publishes half-written.** A consumer reads a version only if its `_COMPLETE.json` verifies (see [Publish protocol](#publish-protocol)).

## Prefix map

```
TCFD/
├── tcfd/                                    # Product 1: annualized, 8 value classes
│   └── layers/
│       └── {layer_id}/
│           ├── v2026-07-27_3412446/         # immutable version
│           │   ├── data/{var}_{scenario}_processed.nc
│           │   ├── qa/qa_report.html, qa_report.json
│           │   ├── maps/index.html, {var}_{metric}_{scenario}.html,
│           │   │        contact_sheet.html, maps_bundle.zip
│           │   ├── layer.json               # science manifest: provenance + decisions
│           │   └── _COMPLETE.json           # sha256 gate — written LAST
│           ├── current/{var}_{scenario}_processed.nc     # copy of active version's data/
│           └── _VERSION.json
│
├── water-index/                             # Product 2: monthly, 20 value types
│   └── variables/{tws,dis,qr,potevap,rootmoist,precip}/
│       ├── v{date}_{sha}/  data/ qa/ variable.json _COMPLETE.json
│       ├── current/
│       └── _VERSION.json
│
├── exports/{customer}/{run_date}/*.csv      # Export-Key deliverables + *.meta.json sidecars
├── reference/                               # shared inputs: land mask, region polygons, hazard tables
├── raw/isimip/{layer_id}/                   # staging, deleted by `cleanup` (see below)
└── _registry/
    ├── layers.json                          # index of every layer × version
    └── processing_log.yaml                  # migrated from ./outputs/
```

Underscore-prefixed names are control/metadata objects, per house convention (`_manifest.json`, `_COMPLETE.json`, `_dev/` elsewhere in the org).

## Deliberate deviations from house convention

Both are intentional decisions, recorded so they aren't "corrected" later:

1. **Raw ISIMIP members are staged in this bucket, not the datasets bucket.** House convention (`operational/parcels/CLAUDE.md:55-57`) reserves `climate-ai-data-science-shiny-app-data` for light dashboard-facing products and sends heavy data to `climate-ai-data-science-datasets`. Keeping everything under one `TCFD/` prefix was chosen deliberately (2026-07-27) for a single mental model and one place to look.
   - Related: `s3://climate-ai-data-science-datasets/arrakis-data/ISIMIP_data/` already holds ISIMIP source data (13,168 objects, 2.2 GiB, June 2024, zarr, laid out `{variable}/{Historical|Future}/{Biomes}/{annual|monthly}/`) — currently `csoil/` and `npp/`. **`csoil` overlaps `soilcarbon_csoil_annual`; check it before re-downloading.**
2. **Gridded outputs stay NetCDF, not zarr.** House convention is zarr for gridded data. `{variable}_{scenario}_processed.nc` is the documented TCFD product contract in `CLAUDE.md`, and every processor, `generate_maps.py`, and the extraction utilities assume it. Changing format is a separate project from the storage migration.
3. **The `riverflood_fldfrc-*` raw prefixes hold a COARSENED field, not what ISIMIP served** (2026-07-29). Their source is 150 arcsec (4320 × 8640) and the full 216-file set is ~54 GB against 19 GB of local scratch, so `download_fldfrc_flood.py` streams each member: download → sha512-verify against the ISIMIP sidecar → area-preserving 12×12 aggregation to 0.5° → **delete the original**. What lands in raw is the 0.5° annual field, ~36× smaller (~700 MB total instead of 54 GB).
   - This is the one place in the bucket where "raw" is not byte-for-byte upstream, so it is made **auditable rather than silent**: every ingested file records `source_url`, the **sha512 of the 150 arcsec original**, its byte count and the exact transform, in both its own NetCDF global attrs *and* `layer.json`'s `inputs.files[]`. Verified populated: 72/72 per protection layer, 144/144 for the derived event layer.
   - Justified by the fact that **raw staging is transient by contract anyway** (see *Raw staging* below — `cleanup_raw` deletes it once `source_url` + checksum are recorded), so uploading 54 GB only to delete it later is pure churn. And `files.isimip.org` is **not** behind the Anubis anti-bot that guards the `data.isimip.org` API, so re-fetching an original is routine rather than precarious.
   - **Do not "fix" this by re-ingesting the 150 arcsec files.** If a native-resolution product is ever wanted, that is a different layer with a different id, and it is a rewrite rather than a re-ingest: the current processor's member × decade stack would need **86 GB of RAM** at native resolution, and `generate_maps.py` cannot serialize 8.1 M cells per panel. See `CLAUDE.md` and the cost table in the session notes.

### Layers created after the migration

The table above is a **migration** record (pre-S3 layout → canonical id); layers created since were born canonical, so they are listed here instead rather than given a meaningless "migration fix" column:

| `layer_id` | notes |
|---|---|
| `drought_driedarea_annual` | ISIMIP3b/SSP; supersedes `drought_led_annual` |
| `timber_cveg-tempnle_annual`, `timber_npp-tempnle_annual` | ISIMIP2b PFT-resolved conifer pair |
| `riverflood_fldfrc-none_annual` | flood, no protection assumed |
| `riverflood_fldfrc-100yr_annual` | flood, uniform 1-in-100 severity threshold |
| `riverflood_fldfrc-flopros_annual` | flood, actual FLOPROS defence standards |
| `riverflood_fldfrc-event100yr_annual` | derived from the two above it; **its own raw prefix stays empty by design** — inputs come from the `-none` and `-100yr` prefixes, and `raw_entries` records all 144 |

Note the last row is the first layer whose raw prefix is legitimately empty. `stage_raw` on it returns nothing, which is correct and not a fault; `layer.json` still carries the full provenance chain.

## Layer IDs

Grammar: `{hazard}_{variable}_{cadence}` — hazard/asset name, ISIMIP variable code, `annual`|`monthly`.

**The raw prefix uses the same `layer_id` as the processed prefix.** Today they diverge in several places; migration normalizes them.

| Canonical `layer_id` | current processed dir | current raw dir | migration fix |
|---|---|---|---|
| `drought_led_annual` | same | same | — |
| `cyclone_let_annual` | same | `drought_let_annual` | **raw renamed** (holds cyclone data) |
| `wildfire_burntarea_annual` | same | same | — |
| `soilcarbon_csoil_annual` | same | same | — |
| `timber-cwood_cwood_annual` | same | `timber-cwood-future` | raw renamed |
| `oak-timber_cwood_annual` | same | — | — |
| `fish-tcb_tcb_monthly` | same | `fish-tcb` | raw renamed |
| `fish-b30cm_b30cm_monthly` | same | `fish-b30cm` | raw renamed |
| `health-mortality_an-tot-heat_annual` | same | same | — |
| `evgndltr_npp-cveg_annual` | `evgndltr_npp-cveg` | `evgndltr_npp-cveg` | cadence appended |
| `loblolly-pine-proxy_cwood-evgndltr_annual` | same | same | — |
| `loblolly-temperate_cveg-needleleaf-evergreen-tree-temperate_annual` | same | same | — |
| `loblolly-temperate_npp-needleleaf-evergreen-tree-temperate_annual` | same | same | — |
| `groundwater_qg_annual` | flat `data/processed/qg_*.nc` | flat `data/raw/` | given a layer dir |
| _tebrsu_ | CLI arg (`{var}-tebrsu_{scen}_processed.nc`) | CLI arg | **ID to assign at migration** |

## Versions

`v{YYYY-MM-DD}_{short-git-sha}` — e.g. `v2026-07-27_3412446`. UTC publish date plus the 7-char commit that produced it, matching the house run-id convention (`experiments/CLAUDE.md:202`). Chronologically sortable and commit-traceable.

Two refinements keep immutability real rather than aspirational:

- **A dirty working tree appends `-dirty`** (`v2026-07-27_3412446-dirty`). The SHA alone does not identify the code that ran if there are uncommitted changes, and pretending otherwise makes a version unreproducible.
- **A same-day, same-commit republish collides by construction** — the common case when iterating on a processor without committing. `publish_layer_version` therefore **refuses by default** (`on_exists="error"`) rather than overwriting. Pass `on_exists="bump"` to publish alongside as `-b`, `-c`, … , or `on_exists="overwrite"` only to replace a known-bad publish.

Without that guard, a rerun silently overwrites a published version — including "repairing" a corrupted one, which would hide the very problem `_COMPLETE.json` exists to catch.

`current/` is a **server-side copy** of the active version's `data/` objects (`CopyObject`, no re-upload). It is the stable read path for downstream consumers; ~10s of MB per layer, negligible against the reproducibility benefit.

## `layer.json` — science manifest

Records provenance, ensemble composition, and the **processing decisions that `GUARDRAILS.md` makes load-bearing** — decadal statistic, normalization, smoothing, percentile mode and direction, trend definition.

Most fields are already emitted by the processors as NetCDF global attributes (`percentile_direction`, `baseline_decade`, `baseline_source`, `window_years`, `n_members`, `impact_models`, `gcms`, `trend_units`, `ci_definition`, `source_dataset`, …), so the manifest is **derived from the processed file plus S3/git context** rather than hand-maintained.

```json
{
  "schema_version": 1,
  "layer_id": "wildfire_burntarea_annual",
  "product": "tcfd",
  "version": "v2026-07-27_3412446",
  "created_utc": "2026-07-27T18:22:05Z",
  "created_by": "scripts/process_burntarea_fire.py",
  "git_commit": "3412446",
  "supersedes": "v2026-07-20_cdc892b",

  "variable": "burntarea",
  "units": "%",
  "cadence": "annual",
  "simulation_round": "ISIMIP2b",
  "sector": "biomes",
  "scenarios": ["rcp26", "rcp60", "rcp85"],
  "decades": [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090],

  "decisions": {
    "decadal_statistic": "median",
    "normalization": "none",
    "spatial_smoothing": "none",
    "percentile_mode": "single_tier",
    "percentile_direction": "higher_is_worse",
    "trend_definition": "Theil-Sen slope of the DECADAL MEDIAN series over an expanding window anchored at 2020s ...",
    "trend_method": "theil_sen_on_decadal_median_series",
    "trend_units": "% decade-1",
    "significance_method": "mann_kendall_tie_corrected_asymptotic_two_sided_on_ensemble_mean_annual_series",
    "baseline_decade": 2020,
    "baseline_source": "shared_across_all_scenarios"
  },

  "ensemble": {
    "n_members_per_scenario": 12,
    "impact_models": ["lpj-guess", "lpjml", "mc2-usfs"],
    "gcms": ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
  },

  "inputs": {
    "raw_prefix": "TCFD/raw/isimip/wildfire_burntarea_annual/",
    "files": [
      {"name": "…nc4", "source_url": "https://files.isimip.org/ISIMIP2b/OutputData/biomes/…",
       "bytes": 123456789, "sha256": "…"}
    ]
  },

  "notes": "See WORKFLOW-ISSUES.md 2026-07-24; GUARDRAILS §9–§10."
}
```

Recording `inputs.files[].source_url` + `sha256` is what makes raw deletion safe: a deleted raw prefix can be re-fetched exactly. This matters because the ISIMIP repository API is behind Anubis anti-bot (`GUARDRAILS.md` §8) — re-downloading is possible but not casual.

## `_COMPLETE.json` — publication gate

Adopted from the house implementation (`operational/parcels/parcel_lib.py:666-706`, `write_completion_marker` / `verify_completion`). Holds the sha256 and byte count of **every** artifact in the version prefix, including `layer.json`, and is written **last**.

```json
{"completed_utc": "2026-07-27T18:22:11Z",
 "artifacts": {"data/burntarea_rcp26_processed.nc": {"sha256": "…", "bytes": 0},
               "layer.json": {"sha256": "…", "bytes": 0}}}
```

**Consumers must verify it and refuse to proceed on absence or mismatch.** A version prefix without a verifying `_COMPLETE.json` is in-flight or corrupt, never data.

## `_VERSION.json` and the registry

```json
{"layer_id": "wildfire_burntarea_annual", "current": "v2026-07-27_3412446",
 "updated_utc": "…", "history": ["v2026-07-20_cdc892b", "v2026-07-27_3412446"]}
```

`_registry/layers.json` indexes every layer so a consumer can discover what exists without listing S3:

```json
{"schema_version": 1, "updated_utc": "…",
 "layers": {
   "wildfire_burntarea_annual": {
     "product": "tcfd", "current": "v2026-07-27_3412446",
     "versions": ["v2026-07-20_cdc892b", "v2026-07-27_3412446"],
     "variable": "burntarea", "units": "%",
     "scenarios": ["rcp26", "rcp60", "rcp85"],
     "percentile_direction": "higher_is_worse"}}}
```

## Publish protocol

Write order matters — a consumer must never see a half-published version:

1. Upload `v{ver}/data/`, `v{ver}/qa/`, `v{ver}/maps/` — each object written to a `…​.writing-{uuid}` temp key then `fs.mv`'d into place (house atomic-publish idiom).
2. Upload `v{ver}/layer.json`.
3. Upload `v{ver}/_COMPLETE.json` — sha256 of every artifact above. **This marks the version consumable.**
4. Server-side copy `v{ver}/data/*` → `current/`.
5. Write `_VERSION.json`.
6. Update `_registry/layers.json`.

A failure before step 3 leaves an orphaned incomplete version that consumers ignore and `cleanup` can delete.

## Raw staging

`TCFD/raw/isimip/{layer_id}/` holds downloaded ISIMIP members. Since we cannot add a lifecycle expiry rule, cost is controlled in code: a `cleanup` step deletes a layer's raw prefix **only after** the published version's `_COMPLETE.json` verifies and `layer.json` records every input's `source_url` + `sha256`. Deletion is opt-in and dry-run by default.

## Credentials — read this before writing S3 code

SageMaker task-role tokens are temporary (~1 h). The org's costliest documented incident came from pinning a **static** copy of them into `AWS_*` env vars: the copy cannot refresh, so long jobs died ~90 minutes in with `OSError: Bad Request` (`operational/scaling-issues-and-solutions.md` §2).

**Do this** — drop static creds so botocore uses the auto-refreshing container-credential provider (`operational/test/finetune_togo.py:56-69`):

```python
for k in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(k, None)          # keep AWS_CONTAINER_CREDENTIALS_RELATIVE_URI
s3fs.S3FileSystem.clear_instance_cache()
```

Re-assert it before each fresh round of S3 I/O in a long job. If the provider is absent this fails at the first read (fast), not 90 minutes in.

**Do not** hand-refresh creds from the container metadata endpoint into env vars. **Do not** use named profiles or static keys — there is no `~/.aws/config` or `credentials` on these instances; everything is the SageMaker IAM role.

Fallback for stubborn `ExpiredToken` cases: shell out to `aws s3 cp` with an `AWS_*`-stripped env, which re-resolves IAM per invocation. House rule is that *every* such subprocess — read or write — must use the stripped-env helper; violating it once produced a silent empty-result corruption (`experiments/landfire-lessons-learned.md:143`).

## Access idioms

`s3fs` is the primary interface, matching the org's `s3_filesystem` / `s3_mapper` / `s3_object_exists` trio (`data-collection/v1-derived-features/lib.py:132-148`). Note the **scheme-less path convention** — s3fs takes `bucket/key`, so URIs are passed through `.replace("s3://", "")`.

## Local cache

Per the house rule, the local filesystem holds **no persistent data**. The cache is ephemeral scratch under `/tmp`, mirroring the S3 key path so the mapping is trivial to reason about:

```
/tmp/tcfd-cache/tcfd/layers/wildfire_burntarea_annual/v2026-07-27_3412446/data/burntarea_rcp26_processed.nc
     ↕
TCFD/tcfd/layers/wildfire_burntarea_annual/v2026-07-27_3412446/data/burntarea_rcp26_processed.nc
```

Default root is `/tmp/tcfd-cache`; override with `TCFD_CACHE_ROOT`. Deleting the cache must never lose data. NetCDF and GeoTIFF writes stage here before upload, since those libraries cannot write to S3 directly.
