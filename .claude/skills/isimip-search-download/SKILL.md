---
name: isimip-search-download
description: Search the ISIMIP repository for climate/impact variables and download members into local raw storage. Use when looking for a new variable, enumerating which models/GCMs/scenarios exist, checking data availability before processing, or ingesting raw ISIMIP files.
---

# Search & Download (ISIMIP)

## Search Principles

These are the standing rules. They predate every mechanical lesson below and outrank them.

**ISIMIP First.** Always search the ISIMIP repository before considering external data
sources. Push back if the user requests non-ISIMIP data unless there's a compelling
reason (e.g. ISIMIP has no relevant data after a thorough search).

**Global Over Regional.** Prefer global-coverage datasets over regional/local data when
available — global biome models (LPJmL, ORCHIDEE, VISIT, CLM) over site-specific models.
This ensures consistent methodology across geographies.

**No Direct Match Strategy.** When exact data isn't available in ISIMIP:
1. List close-enough proxy variables (e.g. NPP/GPP as proxy for growth, vegetation carbon
   for biomass)
2. Identify relevant plant functional types (e.g. "temperate broadleaf deciduous" for oak)
3. Suggest which ISIMIP sectors/models might contain usable approximations
4. Only after exhausting ISIMIP options, discuss external alternatives with user approval

**Favor data abundance over restrictiveness — when in doubt, include more datasets.**

**Maximize dataset count.** Always prioritize including the maximum number of
datasets/models/GCMs available. More data points improve statistical robustness and
uncertainty quantification. Never settle for single-model analysis when multi-model
ensembles are available.

## Answering "what data exists for {hazard}?"

This is the most common request, and the easiest one to answer beside the point. The user
is asking **what the repository holds for that hazard**. They are not asking how it compares
to whatever shipped last.

**Stay inside the hazard.** Enumerate that hazard's representations and judge each on **its
own** properties — cadence, ensemble depth (models × GCMs), units, measured data nature,
scenario and year coverage, download volume, and what the variable physically means for an
asset. Those properties decide it.

**Do not rank datasets by resemblance to an already-processed layer.** "Aligns with the
`csoil` round", "same shape as the `let` processor", "matches the units of the existing 2b
layer" — none of these are properties of the candidate data, and none of them make it the
right answer for *this* hazard. A previously shipped layer is precedent to **re-verify**,
never a reason to prefer. Mention another layer only when the user asks, or when it is a
genuine constraint: an existing processor that can be reused, the OUTPUT-SPEC contract, or a
guardrail. Where the current product already ships this same hazard, say so in one line —
round, scenarios, processor — so the user knows what is being extended or replaced, and stop
there.

**STOP AT THE INVENTORY. This is a hard stop, not a section break.** Report *everything
found* as the matrix described under "Reporting availability", state the open decisions it
depends on (sub-annual aggregation per GUARDRAILS §1–§2, soc/sens harmonization, the
data-nature value check per §9), and **end the turn there**. Choosing the dataset is a
discussion the user opens, not a question you close.

In that turn, do **not**:
- attach a recommendation, a "primary/runner-up" ranking, or a preferred option;
- call `AskUserQuestion` to make them pick;
- pre-answer the open decisions or start downloading.

The user's product judgement is the input the choice waits on, and an inventory delivered
with a recommendation stapled to it forecloses the discussion it was supposed to open. Give
your reading of the trade-offs **when asked** — then the choice is made together. Do not
silently drop an option because you would not pick it.

**Framing choices are not transferable.** How another hazard resolved normalization,
smoothing, percentile tiering, or `higher_is_better` tells you nothing about this one. Those
follow from *this* variable's measured values, and the measurement has not happened yet at
search time. Flag them as pending, do not pre-answer them from precedent.

**If you name an alternative dataset, ENUMERATE it before offering the choice.** Naming one
and quantifying the other is not a fair fork — it is a recommendation wearing a question
mark, because the option with numbers is the one that sounds real. Measured 2026-08-11: the
`led` rebuild was offered against "an SSP `driedarea` layer… a different product with a
different ensemble", with no matrix attached, and the enumeration that came *afterwards*
showed a complete 3 models × 5 CMIP6 GCMs × 3 SSPs = 45-file set, uniform `2015soc`, 168 MB
— cheap to enumerate (15 listings, ~2 min) and materially different from the prose. Offering
to enumerate "if you'd rather" pushes the cost of an informed choice onto the user. Do the
listings first, put both matrices on the table, then stop.

## Search Catalog Reference

**Before starting any ISIMIP search**, consult `config/isimip_search_catalog.yaml` to check
if the variable, model, or scenario has already been investigated. It records variable
definitions and units, available PFT codes, which impact models and climate forcings have
data, temporal and spatial coverage, and scenario availability. Using it avoids redundant
API queries and contextualizes new results. When a new search yields useful information not
already there, update the file with the search date, dataset counts, and any notes on data
quality or limitations.

**But the catalog is opportunistic, not authoritative for coverage.** Its per-model
`variables:` lists were compiled for earlier searches. It has been wrong more than once: it
under-listed ISIMIP2b biomes models (5 documented vs 11 real), and later still missed
`CLM45` and `VEGAS` for `csoil`. **Enumerate the file server; treat the catalog as a
starting hypothesis.**

**Its NEGATIVES are the most dangerous entries in it** (GUARDRAILS §11). A positive claim
in the catalog is self-correcting — you go to the path and either the files are there or
they are not. A negative claim ends the search before it starts, so it is never tested. The
drought section asserted *"ISIMIP3b/SSP version of this family NOT found (0 hits)"* from
**2026-07-24 to 2026-08-11**, and it was false: 3b publishes `driedarea` under
`DerivedOutputData/Heinicke2026`, a complete 3 × 5 × 3 = 45-file SSP matrix. The `led`
drought layer shipped on rcp26/rcp60 while ssp126/370/585 sat unlisted. Worse, from
2026-08-07 **this skill already said the re-issue existed** — the repository contradicted
itself for four days and nobody reconciled it.

So, before you rely on any negative in the catalog:

- **Demand its receipt.** A negative is trustworthy only with `verified_absent_on: "<date>
  — listed <URL>"` naming the directory actually enumerated. No receipt → treat as
  `UNVERIFIED` and enumerate it yourself.
- **Check whether the negative is about a CODE or about the HAZARD.** "`led` returns 0 hits
  in 3b" is true and nearly useless: the token is 2b-only, the hazard is not. Re-issues get
  **new names** (`driedarea`, `floodedarea`, `heatwave`, `wildfire`, `cropfailure`).
- **If the catalog and this skill disagree, stop and enumerate**, then write the answer back
  to *both*. Believing the skill and leaving the catalog stale is how this persisted.
- **Never hand the user a recorded negative as fact** without saying which enumeration and
  date it rests on.

**The catalog's PFT guidance is written as recommendations, and an untested suggestion
reads exactly like a measured fact.** Its top recommendation for temperate conifers pointed
at MC2-USFS biome classes for years; MC2 publishes **no carbon or NPP resolved by biome
type** (only `-total`/`-tree`/`-grass`), its 47 biome names are **one-hot 0/100 presence
flags**, and `pft-maritimeevergreenneedleleafforest` is **identically zero in every run**.
Before relaying a catalog recommendation, confirm the variable it names actually exists at
the resolution it claims — and when one turns out unusable, **withdraw it in the catalog**
rather than leaving it to be found again.

Summary tables should include: variable code and description, simulation round
(ISIMIP2a/2b/3a/3b), time step (annual/monthly/daily), land surface models available,
climate scenarios, and dataset counts.

Reporting principles:
1. Show all data sources (include all simulation rounds)
2. Prefer newer data (ISIMIP3b over ISIMIP2a) — **but check depth**: 2b is sometimes far
   deeper (2b `csoil` has 9 annual models vs 3b's 3), which can be worth a cross-round
   sensitivity check even when 3b is primary
3. Emphasize multi-model ensemble availability
4. Maximize dataset count (see Search Principles)

## ISIMIP Data Context

**Simulation Rounds:**
- ISIMIP3a/3b: SSP scenarios (ssp126, ssp370, ssp585)
- ISIMIP2a/2b: RCP scenarios (rcp26, rcp60, rcp85)

**Key Variables:** `led` (drought), `leh` (heatwave), `lew` (wildfire), `qg` (groundwater
runoff), `burntarea`, `potevap`, `npp`, `gpp`, `cveg`, `b30cm` (large fish biomass)

**Documentation Updates:** When new simulation rounds, variables, or data sources are
discovered in ISIMIP, work with the user to update this documentation. The ISIMIP
repository evolves; keep this guide current.

## Scenario Handling

**Projection scenarios** (included in reports):
- SSP: ssp126, ssp245, ssp370, ssp585
- RCP: rcp26, rcp45, rcp60, rcp85

Availability is per-variable — ISIMIP3b `biomes` has **no ssp245** at all.

**Scenario discovery**: Scripts must discover scenarios dynamically from the filesystem,
not hardcoded lists. See GUARDRAILS.md §3. (A hardcoded `RCP_SCENARIO_LABELS` dict once
made 25% of processed data invisible — rcp45 was silently dropped from every map.)

**Non-projection scenarios** (excluded from reports): `picontrol`, `historical`. These may
be downloaded and processed alongside projection data to strengthen the 2020s baseline, but
are **automatically excluded** from `generate_maps.py` report generation.

**New Scenarios:** If you encounter scenarios not listed above, discuss with the user
whether to add them to this documentation and update the processing/reporting code
accordingly.

## Experiment & SOC Preferences

**Sensitivity experiments:** Prefer `default` (transient CO2). Ask the user if `2015co2` or
others should be included.

**SOC scenarios:** Prefer `2015soc` (isolates climate signal). Ask the user if `histsoc` or
`2015soc-from-histsoc` should be included.

**Verify soc/sens tokens per model before committing** to a uniform ensemble (GUARDRAILS
§9) — they frequently differ, and a strict filter can silently collapse the ensemble. A
variant sound for one variable can be broken for another: `classic`/`2015soc-from-histsoc`
is fine for `csoil` but mixed-scale for `burntarea`. Never inherit a variant — re-check it.

## Choosing what to download

- **Ask the user for temporal resolution** when both exist (GUARDRAILS §1) — monthly files
  are often ~12× larger. Do not download monthly by default. Use `AskUserQuestion` and
  present the options with file counts and sizes. (28 monthly fish files / 3.62 GB were
  once downloaded when 4 annual files / ~40 MB would have done.) Note that monthly data for
  a smooth *stock* can be annualized cheaply, which may be the way to deepen a thin ensemble.
- **Maximize ensemble depth** — see Search Principles.
- **Prefer global over regional** coverage for methodological consistency.

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
  get rate-limited and return **empty** listings indistinguishable from a genuine "no data".
  13 of 22 directories came back empty on one parallel pass. An empty listing is a *failure
  signal*.
- Confirm a genuinely absent path with its HTTP status (a real 404), not an empty body.
- Each file has a `.json` sidecar with `size`, `sha512` and parsed `specifiers` — use it to
  verify downloads and read soc/sens tokens without downloading the NetCDF. 2b files carry
  these; 3b `biomes` largely does not, so verify 3b by HTTP status instead.
- To inspect a header without downloading a large file, open it lazily over HTTP
  (`fsspec` + `xarray`/`h5netcdf`). **Check the interpreter first — in this repo neither is
  installed.** `python` is not on PATH at all (only `python3`), and `.venv` has `xarray` but
  **no `fsspec`, no `h5netcdf`, no `yaml`**. So there is no lazy remote open and no YAML
  validator available; budget for that instead of discovering it mid-search.
- **Do not chase CF attributes with HTTP range reads.** NetCDF4 is HDF5 and places attribute
  object headers at unpredictable offsets: a 1 MB *and* a 4 MB prefix read of a
  `Zantout2025` file surfaced only the dims (`lat`/`lon`/`time`/`bnds`), the variable name
  and `_FillValue` — **no `units`, no `long_name`**. Units are a download-time measurement
  (GUARDRAILS §9); record them as unverified and move on.

### Harvest only what the question needs (2026-08-08: this cost ~75 min of wall clock)

Enumeration is ~1 HTTP request per directory plus the serial sleeps, so every directory you
list and every file you store is wall time. One wildfire inventory spent 75 minutes; the
measured waste:

- **Reduce during the harvest — but project the variable FIELD, never grep a token you
  believe in.** Storing all 86,065 filenames across three sectors left **87–93% irrelevant**
  (`qtot`, `snd`, `lai`, `tsl`, `trans-*`…) for a question about `burntarea`. The fix is to
  emit `$(NF-4)` (plus scenario `$(NF-7)` / sens `$(NF-5)` when you need them) and reduce to
  a distinct **vocabulary**, then match your target against that vocabulary offline:

  ```bash
  # RIGHT — one pass, tiny output, and the vocabulary is auditable
  curl -s "$URL" | grep -oE 'href="[^"]+\.nc4?"' | sed -E 's/href="([^"]+)"/\1/' \
    | awk -F'_' '{print $(NF-7), $(NF-5), $(NF-4)}' | sort -u

  # WRONG — a filter keyed to a token you have not yet observed
  curl -s "$URL" | grep -E '_yield-sgc-'      # returns nothing; means NOTHING about sugarcane
  ```

  **A pre-filter can only confirm presence, never absence.** Measured 2026-08-11: a
  sugarcane search grepped `sgc` — the code our own catalog listed for ISIMIP3b — across
  150 directories and matched **zero** files. Published sugarcane output is `sug`, and `sgc`
  is real only in the crop-calendar `InputData` vocabulary. The vocabulary projection running
  alongside is the only reason the layer was found; the targeted grep alone would have
  reported "sugarcane does not exist in ISIMIP" off a correct, complete enumeration. If you
  ever do pre-filter, an empty result is **`UNVERIFIED`**, not a negative (GUARDRAILS §8, §11).
- **List `future/` first; touch `historical/`/`pre-industrial/` only when the baseline or a
  specific check needs them.** 65% of 142 directory listings (93 of them) went to
  `historical/` + `pre-industrial/` and were never used by the inventory.
- **Before harvesting a SECOND sector, test whether it is a duplicate of the first — 2
  requests, not 28 minutes.** ISIMIP3b `fire` and `biomes` publish the *same* `burntarea`
  files (2001/2031 basenames shared; a spot-checked pair matched on `Content-Length` **and**
  `ETag`). Compare basenames against the sector you already hold, then `curl -sI` one shared
  pair. A full second-sector walk cost 28 min and added nothing. (Cross-sector duplicates are
  a known ISIMIP pattern — `elm-eca` csoil-total in `biomes`+`permafrost`.)
- **Probe the directory depth with ONE listing before writing a walker.** Depth is not
  uniform: `OutputData` is `{MODEL}/{gcm}/{future|historical}/`, but publications under
  `DerivedOutputData/` vary. Guessing Lange2020's depth twice — a 3-level walk that returned
  empty, then a 2-level walk that found 0 files because the level really is
  `{MODEL}/{gcm}/{future,historical,pre-industrial}/` — burned **~35 min** before one correct
  14-min run. An `http 200` on a directory with 0 matched files means **wrong depth**, not
  "no data".
- **Size by sampling, not exhaustively.** File size is near-uniform within a
  (model, variable, cadence, span) group — all 36 `Zantout2025` wildfire files were *exactly*
  425 MB. HEAD 2–3 per group and multiply; 36 sequential HEADs told us nothing extra.

Filename grammar (parse from the END; model/GCM names contain hyphens, never underscores):

```
{model}_{gcm}_{forcing}_{scenario}_{soc}_{sens}_{variable}_{region}_{timestep}_{y0}_{y1}.nc
                                                     [-5]     [-4]     [-3]   [-2] [-1]
```

**Parsing from the END is not optional for PFT-qualified variables.** The variable field
itself contains hyphens (`cveg-needleleaf-evergreen-tree-temperate`), so any fixed forward
index breaks on exactly the datasets that need it most.

**`DerivedOutputData/` filenames carry a LEADING publication token**, which shifts every
from-the-start field by one and makes a forward index read the wrong column *silently* —
`$4` lands on the forcing (`w5e5`/`ewembi`) and gets reported as the scenario:

```
zantout2025_classic_gfdl-esm4_w5e5_ssp126_2015soc_default_wildfire_global_annual_2015_2100.nc
lange2020_caraib_gfdl-esm2m_ewembi_rcp26_2005soc_co2_lew_global_annual_2006_2099.nc4
^^^^^^^^^^ publication    model=$2  gcm=$3  forcing=$4  scenario=$5  soc=$6  sens=$7
```

This happened twice in one search and both times produced a plausible-looking matrix that
had to be re-derived. From the end it is invariant: `variable=$(NF-4)`, `timestep=$(NF-2)`.

**Mind token ORDER when grepping a filename.** The scenario precedes the variable, so
`grep '_lew_.*rcp26'` matches **nothing** — it must be `'rcp26.*_lew_'` (or two `grep`s).
A zero-match grep here looks exactly like absent data.

**Anchor the extension in harvest regexes: `\.nc4?$`.** An unanchored `[\w.-]+\.nc` matches
the `.nc` *inside* `.nc4` and truncates every ISIMIP2b filename by one character, so every
URL built from the harvest 404s — which reads as "no data" rather than a parsing bug.

## GUARDRAILS §8 — never guess specifier codes

Specifiers are a controlled vocabulary. **`count=0` for something that plausibly exists is a
red flag, not a conclusion** (tropical-cyclone exposure is `let`, not the mnemonic `letc`).
`count=1001` is the API maximum — the result set is truncated and must not be generalized.

- **Enumerate families, don't guess members.** The Lange 2020 exposure family is **twelve,
  not six**: `le{d,r,w,c,h,t}` = **land-area** fraction exposed, each paired with a
  `pe{d,r,w,c,h,t}` = **population** fraction exposed twin. Hazards are `d` drought, `r`
  river flood, `w` wildfire, `c` crop failure, `h` heatwave, `t` tropical cyclone. Also:
  `lew` has **no rcp85** (rcp26/rcp60 only).
- **List EVERY intermediate directory level. Never path-guess past one.** An enumeration
  that skips a level cannot support an absence claim about anything below it. Measured
  cost (2026-08-11): listing `ISIMIP3b/InputData/` returned `climate/`, then the walk
  jumped straight to `climate/atmosphere/` because that was the obvious next hop —
  skipping `climate/` itself, which is where **`tropical_cyclones/`** lives. That directory
  holds the newest TC hazard in the repository (MIT per-storm wind footprints, Frieler et
  al. 2025, ssp126/370/585), and the inventory went out saying no SSP tropical-cyclone
  product existed. The user had to supply the correction.
- **A verified negative about one product family is NOT a negative about the hazard.** In
  the same review it was correctly established that the Lange 2020 exposure family has no
  ISIMIP3b re-issue — and that true finding was then allowed to stand in for "no SSP TC
  data exists". Scope every absence claim to exactly the family you enumerated.
- **Beware publication-directory name collisions.** `ISIMIP3b/DerivedOutputData/TipESM2025/MIT/`
  looks like the tropical-cyclone group (MIT = Emanuel's institution) and is **water models**
  (CWATM, H08, JULES-W2, MIROC-INTEG-LAND, …). The real TC data is under
  `InputData/climate/tropical_cyclones/MIT/`. Open the directory; do not infer from the name.
- **List `DerivedOutputData/` before concluding a product has no newer-round version.** A
  product can be re-issued in a later round under a different publication directory with
  different variable names. Lange 2020's exposure concept *was* re-issued for ISIMIP3b,
  split across `Heinicke2026` (`driedarea`, `floodedarea`) and `Zantout2025` (`heatwave`,
  `wildfire`, `cropfailure`) — hazard words, not `le*` codes. The drought layer nearly
  shipped on rcp26/rcp60 when ssp126/370/585 existed. One call:

  ```bash
  curl -s https://files.isimip.org/ISIMIP3b/DerivedOutputData/ | grep -oE 'href="[^"]+/"'
  ```
- **A variable can have more than one representation.** "Wildfire" is the `lew` exposure
  member, the `ffire` emissions flux, the ISIMIP3a-only `fire`-sector diagnostics, *and* the
  direct `burntarea` burnt-area fraction. Enumerate all products and present the trade-offs.
- **If a catalog entry names a sector, verify that sector was actually walked.** The wildfire
  section was titled "fire sector + biomes burntarea" but only `biomes` had been enumerated.
  ISIMIP3a's `fire` sector holds **10 models found nowhere else**. A section heading is not
  evidence of coverage.
- **File extensions differ by round: ISIMIP2b publishes `.nc4`, ISIMIP3a/3b publish `.nc`.**
  A `*.nc` filter silently drops the **entire 2b round**. Match both.
- **The same quantity can be named differently per round** — ISIMIP2b `csoil` vs ISIMIP3b
  `csoil-total`. Crop codes drift the same way and more often: sugarcane is `sug` (2a/2b
  output) but `sgc` in the 3b crop calendar; beans are `ben` (2a) vs `bea` (3); wheat is one
  code `whe` (2a/2b) but splits into `swh`/`wwh` (3a/3b); rice is `ric` (2a/2b) vs
  `ri1`/`ri2` (3a/3b). **Never carry a crop code across rounds** — project the vocabulary
  for the round you are in.
- **`InputData` vocabulary ≠ `OutputData` vocabulary.** The ISIMIP3b crop calendar defines
  20 crops (`bar bea cas cot mai mil nut pea pot rap ri1 ri2 rye sgb sgc sor soy sun swh
  wwh`); ISIMIP3b models actually publish **11** (`bea cas mai mil pot ri1 ri2 sor soy swh
  wwh`). A code being in the protocol says nothing about a model having run it. When
  recording any vocabulary, record **which product you observed it in**.
- **Watch for cross-sector duplicates.** `elm-eca` publishes csoil-total under both `biomes`
  and `permafrost`, byte-identical (same sha512). Ingesting both double-weights the model.
- When processing one member of a family, record the **whole** family in the catalog then.

## Vegetation / PFT searches — treat variable × PFT × round as one lookup

A PFT-resolved variable exists per **(round, model, variable, PFT)**, and the combination a
product needs frequently does not exist at all. Establish the intersection before promising
anything:

- **PFT-resolved `cwood` exists ONLY in ISIMIP3b, from two models** — `classic` `evgndltr`
  (2 GCMs) and `jules-es-vn6p3` `ndlevg` (5 GCMs). There is **no PFT-resolved `cwood`
  anywhere in ISIMIP2b**. 2b PFT output is `cveg` / `npp` / `gpp` / `pft` only.
- **No climate-zone-resolved conifer PFT carries `cwood` in ANY round.** Wanting both
  climate specificity and the wood-only pool is unsatisfiable — surface that trade-off
  explicitly rather than picking one silently. Substituting `cveg` costs the root+leaf
  fraction: wood is p50 **77–90%** of conifer `cveg`.
- **`cveg` is not universal either** — `lpjml` publishes PFT-resolved `npp`/`gpp`/`pft` but
  **no** per-PFT `cveg`, so a cveg track and an npp track have different member counts.
- **Check the model's own PFT scheme before mapping a species onto it.** `jules` splits
  broadleaf evergreen into `bdlevgtemp`/`bdlevgtrop` but publishes a single `ndlevg`; the
  only 3b `cwood` token containing "temp" is `cwood-bdlevgtemp`, a temperate **hardwood**,
  not a conifer.
- **A generic class is one PFT with one global parameter set**, not an internally weighted
  mixture of climate-specific sub-types. Its cost bites at the **percentile ranking** step,
  not in the raw per-cell values.
- **Climate-zone specificity is not automatically better — it depends on the species.** It
  helped loblolly (one zone). For a hyper-oceanic species spanning zones (Sitka spruce) it
  *fragments* the range: `lpjml` temperate-NLE covers UK/Ireland (774/774 cells) but only
  609/2451 in PNW-BC-SEAK, while `lpj-guess` boreal-NLE covers PNW (1866) and just **6** UK
  cells. Check coverage **in the regions the product is for**.
- **Species-level ISIMIP is 8–9 European stand sites, and there is no Sitka spruce.** The
  forestry sector (2b: 10 models; 2a: 15) is the only place individual species and real
  forestry metrics (`mai`, `vol`, `harv`) appear, but it is not griddable. Species
  vocabulary is exactly `{acpl, bepe, fasy, lade, piab, pipi, pisy, psme, quro}`.
- **`pft-{code}` cover fractions are published alongside** by most models and answer "does
  this model actually place this plant here". `clm45` publishes none.

## Downloading

Raw goes to `data/raw/{layer_id}/`, where `layer_id` is `{hazard}_{variable}_{cadence}` and
**must match the processed layer's name**.

**Prefer the CLI**: `isimip-pipeline run` handles multi-scenario downloads into a single
output directory. Avoid manual `search` + `download` workflows that fragment scenarios into
separate folders.

- Verify each download against its sidecar `size` and `sha512` before use.
- Make the downloader resumable — skip a file already present and matching.
- Record each file's `source_url` and checksum. The API being behind anti-bot means an
  unrecorded source may not be casually re-downloadable.
- Use ISIMIP's server-side masking (`mask_landonly`, `select_bbox`, `mask_country`) to cut
  download volume where the product only needs land or a region.

## Reporting availability

Report as a matrix — model × GCM × scenario, with cadence, year span, soc/sens tokens and
dataset counts — plus the temporal record and the spatial resolution.

**Report resolution as DECLARED, not verified, until you have measured it from the values.**
Every ISIMIP file reports the round's nominal grid (0.5° / 360×720), including members that
ran natively coarser and were replicated onto it — `csoil`'s `elm-eca` is effectively
~4°×5° with an identical-looking header. Reading `ds.sizes` describes the container, not the
information content. Effective resolution is measured at processing time — don't promise it
here.

Then hand off to the `isimip-process-visualize` skill. The output contract it must satisfy
is [OUTPUT-SPEC.md](../../../OUTPUT-SPEC.md).
