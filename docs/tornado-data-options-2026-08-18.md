# Tornado occurrence data — what exists, and what each option would cost us

**Date: 2026-08-18.** Availability review, opened by the question "is there tornado data in
ISIMIP, and if not what are our options elsewhere". This is an inventory, not a
recommendation: it reports what was found, what was measured, and the decisions the choice
depends on. No option is ranked.

---

## 0. What the enumeration established

1. **ISIMIP publishes no tornado data in any round.** Verified by directory listing today,
   receipt in §1. Scoped precisely: what was listed, and what was not.
2. **A tornado proxy cannot be constructed from ISIMIP forcing either.** The bias-adjusted
   atmosphere carries **11 surface variables and nothing aloft**. Every tornado-environment
   index needs vertical wind shear, which needs winds at height. Those are not published.
   This is a hard stop, not a gap that effort closes (§2).
3. **A convective product we had never enumerated does exist** — `InputData/climate/lightning/`,
   flagged as *unswept* in our own catalog since 2026-08-14. It is a lightning flash-rate
   field from **one GCM**, and it is a thunderstorm-activity proxy, not a tornado proxy (§1.3).
4. **Outside ISIMIP there are three classes of option, and they answer different questions**:
   occurrence records (§3), environment proxies (§4), loss records (§5). Only the environment
   proxies are projectable under a scenario.
5. **The obstacle shared by every occurrence record is measured, not theoretical.** In the US
   database — the best one in the world — reports of the weakest tornadoes rose **7.7×** since
   the 1950s while reports of strong tornadoes **fell**. The raw count series is a record of
   detection capability, not of climate (§3.2). That confound points at population density,
   which is the exposure we are supposed to score independently.
6. **Tornado has no hazard family in `config/hazard_taxonomy.yaml`** — 20 families, none of
   them this one — so it is currently invisible in the mandatory "hazards not assessed"
   section of every report (§7).

---

## 1. ISIMIP — the receipt

Serial listings, 3 s spacing, every level enumerated rather than path-guessed
(GUARDRAILS §8, §13). No empty listings.

### 1.1 `DerivedOutputData/` refreshed, all four rounds

| Round | Publications found |
|---|---|
| ISIMIP2a | `Zimmer2023` |
| ISIMIP2b | `Lange2020`, `Zimmer2023` |
| ISIMIP3a | `Quesada-Chacon2026`, `Zantout2025` |
| ISIMIP3b | `Heinicke2026`, `Jaegermeyr2021`, `TipESM2025`, `Zantout2025` |

Identical to the standing map in the `/isimip-search-download` skill. No new publication has
appeared, and none of the nine is convective. The Lange 2020 exposure family is six hazards
(`d` drought, `r` river flood, `w` wildfire, `c` crop failure, `h` heatwave, `t` tropical
cyclone) — no wind hazard other than TC, and no convective member.

### 1.2 `InputData/climate/` walked level by level

```
ISIMIP3b/InputData/           climate/  composition_atmosphere/  geo_conditions/  socioeconomic/
ISIMIP3b/InputData/climate/   atmosphere/  atmosphere_composition/  lightning/  ocean/  tropical_cyclones/
ISIMIP3a/InputData/climate/   atmosphere/  atmosphere_composition/  lightning/  ocean/  sealevel/
                              tropical_cyclones/  tropical_cyclones_flooding/
ISIMIP3b/SecondaryInputData/  climate/  geo_conditions/  socioeconomic/
```

`tropical_cyclones/` is the storm hazard ISIMIP does carry, and it is a different hazard —
TC and tornado share the word "wind" and nothing else physically.

### 1.3 `lightning/` — enumerated for the first time

Our catalog (line 3787) lists `3b InputData/climate/lightning` among directories **"NOT
swept"**. It is swept now:

| Round | Files | GCM | Scenarios | Cadence | Span | Size/file |
|---|---|---|---|---|---|---|
| ISIMIP3b | 4 + `lightning_fixed.nc` | `ukesm1-0-ll` **only** | ssp126, ssp245, ssp370, ssp585 | monthly | 2015–2100 | ~505 MB |
| ISIMIP3a | `lightning_fixed.nc` only | — | — | fixed | — | 6.6 MB |

Read from the `.json` sidecar, not inferred:

- variable `lightning`, **"Monthly Flash Rate"**, units **`km-2 d-1`**, 720×360 (0.5°), 1032 months
- `actual_range` 0 – 6.29
- title: *"estimated lightning stroke density using LPJ_cape and WGLC log(CAPE):log(WGLC)
  linear model"*
- source: Kaplan, Koch & Lau (2023), Zenodo `10.5281/zenodo.7511843`, prepared for ISIMIP3b
  by the ISIMIP data team

Two things follow. First, **one GCM is not an ensemble** — per the skill's per-GCM counting
rule this is a 1-of-5 fragment, and OUTPUT-SPEC's `n_models`/`n_members` and the CI would be
degenerate. Second, and more useful: this file is **derived from CAPE**, which ISIMIP itself
had to import from outside the forcing archive. That is the tell for §2.

Note also `ssp245` exists here — it does not exist in the bias-adjusted atmosphere.

### 1.4 The negative, scoped

> **No tornado, hail, thunderstorm, derecho or severe-convective variable is published in
> ISIMIP2a, 2b, 3a or 3b.** Verified 2026-08-18 by listing: `DerivedOutputData/` root in all
> four rounds; `ISIMIP3b/` round root; `ISIMIP3a|3b/InputData/` and their `climate/`
> subtrees; `ISIMIP3b/SecondaryInputData/`; `climate/lightning/` in 3a and 3b; and the full
> `atmosphere/bias-adjusted/global/daily/` depth in 3b.
>
> **Not listed, so not covered by this claim:** every `OutputData/` impact sector (no impact
> sector models tornadoes, but that is reasoning, not enumeration); `DerivedInputData/`;
> `atmosphere_composition/`, `ocean/`, `geo_conditions/`, `socioeconomic/`; and every
> `InputData` root in 2a/2b.

The unlisted directories are where a surprise would have to live, and on the evidence of the
TC incident — where skipping exactly one level hid the newest storm product in the
repository — "implausible" is not the same as "checked". The claim above is worth what its
listing is worth and no more.

---

## 2. Why we cannot build the index ourselves from ISIMIP

Every tornado-environment index in the literature — STP (significant tornado parameter),
Tippett's TEI, the supercell composite — is built from the same ingredients:

| Ingredient | What it needs | In ISIMIP forcing? |
|---|---|---|
| CAPE | temperature + humidity **profile** through the depth of the troposphere | **no** |
| 0–6 km bulk shear | wind vector at ~500 hPa **and** at the surface | **no** (surface only) |
| Storm-relative helicity (0–1 km, 0–3 km) | wind **profile** in the lowest kilometres | **no** |
| LCL height | surface T + dewpoint | derivable |
| CIN | temperature profile | **no** |

Measured today, ISIMIP3b `atmosphere/bias-adjusted/global/daily/ssp370/GFDL-ESM4/` — 99
files, projecting the variable field gives exactly **11 variables, 9 files each**:

```
hurs  huss  pr  prsn  ps  rlds  rsds  sfcwind  tas  tasmax  tasmin
```

All surface. `sfcwind` is a **daily mean** 10 m wind — not a gust, and not a wind aloft.
There is no pressure-level data, no geopotential, no upper-air temperature or humidity, in
the bias-adjusted product at all.

So the shear term — which is *the* discriminator between a thunderstorm and a tornadic
thunderstorm, and the reason CAPE alone is not a tornado predictor — is simply absent. One
could compute a CAPE-like instability proxy from surface temperature and humidity, and it
would be a **thunderstorm-favourability** index that says nothing about rotation. Publishing
that as a tornado layer would be the kind of error that is invisible because every number in
it is correct.

**This is the finding that shapes the rest of the review**: any tornado-relevant projection
has to come from a source that publishes the atmosphere in the vertical, which means raw
CMIP6, a reanalysis, or a convection-permitting model — never ISIMIP.

---

## 3. Class 1 — occurrence records (observed tornadoes)

These answer *"where have tornadoes been reported, and how hard did they hit"*. All are
historical; none carries a scenario.

| Source | Region | Period | Volume | Access | Licence |
|---|---|---|---|---|---|
| **NOAA SPC severe weather database** | US (+ AK, HI, PR) | 1950–2025 | **73,458 tracks** (measured) | direct CSV, no auth | US Gov public domain |
| **NOAA NCEI Storm Events** | US | 1950 – Feb 2026 | all event types, incl. hail/wind | bulk CSV over HTTP/FTP | US Gov public domain |
| **Northern Tornadoes Project** (Western Univ.) | Canada | 1980–present (1980–2009 and 1991–2020 revised sets) | several hundred previously unreported events added | ArcGIS open-data portal, KML/API | **non-commercial free; commercial requires contacting `ntp@uwo.ca`** |
| **ESWD** (ESSL) | Europe + Mediterranean | multi-decadal | "many tens of thousands" of reports, all severe types | web UI; bulk by application | **signed User's Agreement; commercial use not allowed on public data; fee for funded/organisational use** |
| **TORRO** | UK | ~1000 yr historical, modern era usable | >4,000 events (UNVERIFIED) | by request | research org, terms unstated |
| **BoM Severe Storms Archive** | Australia | 18th C – present | — | web archive | Australian Gov |
| **JMA tornado database** | Japan | 1961–2016, ~980 events (UNVERIFIED) | — | web | public domain |
| **Tornado Archive** | **global** | varies by country | **>100,000 tornadoes** | web Data Explorer; bulk terms unclear | **no commercial use without written permission** |

Two access notes. The **Tornado Archive** (Boyd et al., *BAMS* 105(7), 2024) is the only
genuinely global compilation and would otherwise be the obvious starting point; its sources
page and the BAMS article both returned **HTTP 403** to automated fetching today, so its
per-country coverage table and its bulk-download route are **UNVERIFIED** — that needs a
human visit or an email before anyone plans work on it. It also inherits ESWD's restriction,
so it is not licence-clean for a commercial product by default. The **ESWD** is the
authoritative European source and is explicitly **not free for commercial use** — for our
product that is a procurement conversation with ESSL, not a download.

### 3.1 What the SPC file actually contains

Measured today: `1950-2025_actual_tornadoes.csv`, 9.0 MB, 73,458 records, 29 fields
(`om,yr,mo,dy,date,time,tz,st,stf,stn,mag,inj,fat,loss,closs,slat,slon,elat,elon,len,wid,ns,sn,sg,f1..f4,fc,edat,etime`).
Start/end lat-lon, path length and width, magnitude, casualties, loss. 123 records have zero
track length; 1,863 carry `fc=1`, meaning the F-scale was **estimated after the fact** by
NCEI rather than assigned at survey.

SPC's own caveat, quoted from the WCM page:

> "these data are used by the NWS for verification purposes and may not accurately reflect
> all storm events. Monetary loss information is highly suspect and should be used with
> caution, if at all."

### 3.2 The measurement that governs this entire class

Per-year rates by decade, computed from the file today (2020s = 6 years, normalised):

| period | yrs | all/yr | **F0/yr** | F1/yr | **F2+/yr** | unrated/yr |
|---|---|---|---|---|---|---|
| 1950–1959 | 10 | 479.3 | 105.3 | 185.2 | **188.8** | 0.0 |
| 1960–1969 | 10 | 681.1 | 182.4 | 255.0 | **243.7** | 0.0 |
| 1970–1979 | 10 | 857.9 | 244.8 | 367.0 | **246.1** | 0.0 |
| 1980–1989 | 10 | 819.5 | 330.1 | 330.9 | **158.5** | 0.0 |
| 1990–1999 | 10 | 1213.7 | 737.1 | 327.2 | **149.4** | 0.0 |
| 2000–2009 | 10 | 1277.9 | 813.5 | 337.7 | **124.6** | 2.1 |
| 2010–2019 | 10 | 1207.0 | 631.1 | 408.8 | **137.2** | 29.9 |
| 2020–2025 | 6 | 1349.0 | 481.0 | 510.2 | **153.5** | 204.3 |

Reported tornadoes rose **2.8×** (479 → 1349 per year). Weak (F0) reports rose **7.7×** at
peak. Strong (F2+) reports **fell 19%**, from 189/yr to 154/yr, and the F2+ series has no
trend in the modern era at all — annual values 2006–2025 run 85 to 283 with no direction.

Strong tornadoes destroy buildings and have always been recorded. Weak tornadoes are found
only if somebody is there to see one, and since 1950 the US has added the WSR-88D radar
network, organised spotters, a much larger population and a camera in every pocket. The
growth is in the category that detection explains and absent from the category it doesn't.

**A third break, and a recent one**: unrated reports (`mag = -9`) were unheard of before 2007
and now run **204–268 per year**, roughly 15% of all reports (2019: 185, 2023: 268, 2025: 234).
So the naive fix — "just use F2+, it's homogeneous" — is itself degrading, because a growing
slice of the record now has no rating to filter on. And the F/EF scale changed in
February 2007, which is a rating-system discontinuity in the middle of the series.

For our product the consequence is specific and worse than a normal data-quality caveat.
A hazard layer built on report density would score risk highest where reports are densest,
and report density tracks population. That is not a caveat one discloses one's way out of —
it is a **confound with the exposure term**, and it would systematically under-score risk at
exactly the remote, low-population industrial sites (pipelines, solar farms, substations,
rail) that a physical-risk product is most often asked about.

---

## 4. Class 2 — convective environment proxies (the only projectable class)

These answer *"how often does the atmosphere here support a tornado"*. They are what the
literature actually uses for climate projection, precisely because of §3.2.

### 4.1 ERA5 (reanalysis, observed period)

CAPE and CIN are **native ERA5 single-level hourly fields**; 0–6 km shear and storm-relative
helicity are derived from pressure-level `u`/`v`. This is the substrate for Taszarek, Allen,
Marchio & Brooks, *Global climatology and trends in convective environments from ERA5 and
rawinsonde data* (npj Clim. Atmos. Sci., 2021), and for Taszarek et al. (J. Climate 33(23),
2020) which is the reference mapping ERA5 environments onto observed lightning, hail, wind
and **tornadoes** across both Europe and the US.

**Volume — computed from grid dimensions, not measured.** ERA5 is 0.25° (1440 × 721 =
1,038,240 cells). Hourly for 1940–2025 is ~753,900 steps; one 4-byte field is therefore
**~3.1 TB uncompressed**, before any of the four-to-five fields an index needs. Against the
~600 GB budget that is a non-starter at native cadence. Reductions available:

- CDS **`derived-era5-single-levels-daily-statistics`** aggregates to daily max/mean/min at
  retrieval time (it is computed on request, not a stored archive) — ÷24
- restrict to a 30-year window (e.g. 1991–2020) — ÷2.9
- land-only / bbox masking — ÷~3.4 globally, far more for a regional product

A daily-max, 30-year, land-only pull of ~5 fields lands in the tens of GB. That is a
tractable ingest, and it is stream-and-reduce shaped. But it is one realization with no
scenario — see §6.

### 4.2 CMIP6 directly (projections)

Lepore, Abernathey, Henderson, Allen & Tippett, *Future Global Convective Environments in
CMIP6 Models* (Earth's Future, 2021) is the direct precedent: **7 CMIP6 models**, computing
CAPE, CIN, S06, SRH and combined severe-weather proxies. Their model count was set by data
availability — the vertical fields needed live in the **`6hrLev`** table (6-hourly, native
model levels), which few models publish and which is large. Headline result: severe-weather
proxy frequency rises **5–20% per °C** of global warming, driven principally by CAPE.

This is the only class that can produce a scenario-resolved, multi-model layer.

### 4.3 Calibrated hazard models (proxy → actual tornado probability)

The gap between "the environment supports tornadoes" and "tornadoes occur" is closed by
statistical models fitted to occurrence reports:

- **AR-CHaMo** (Additive Regressive Convective Hazard Model) — ESSL; Rädler, Groenemeijer,
  Púčik, Battaglioli, with Taszarek. Predicts hazard probability as *P(storm) × P(hazard |
  storm)*, fitted on ERA5 + ESWD. Developed under STEPCLIM, ARCS (**co-funded by Munich Re**)
  and CHECC/ClimXtreme. Already applied to a 14-member EURO-CORDEX ensemble (Rädler et al.,
  npj Clim. Atmos. Sci., 2019); **CHECC II plans to extend it to CMIP6 for large hail and
  tornadoes in Central Europe.** This is the closest thing that exists to a purpose-built
  tornado hazard layer under climate change — Europe-focused, and not a dataset we can simply
  download.
- **TEI** (tornado environment index, Tippett) and **HEI** (hail, Allen) — Poisson
  regressions calibrated to **CONUS** monthly climatological frequency. Used inside the
  Lepore CMIP6 work. Calibration is regional; applying them outside CONUS is an unresolved
  question, not a configuration flag.

Note what both imply: the calibration is fitted against occurrence reports, so §3.2's
reporting inhomogeneity enters the proxy through the back door. It is handled (fitting on
monthly climatology and on stable subsets rather than raw counts) but it does not disappear.

### 4.4 Convection-permitting simulation

**CONUS404** (NCAR + USGS) — WRF at **4 km**, hourly, 1979/1980–2024, ~200 2-D variables plus
3-D fields, over CONUS. Its companion **CONUS404-PGW** re-runs the same 45 years perturbed by
CESM2-LENS2 warming signals (pseudo-global-warming), giving a matched present/future pair.
Available via USGS ScienceBase and the NCAR RDA/GDEX (`d559000`).

At 4 km this resolves the storms themselves rather than only their environment — this is the
approach behind Woods et al. (GRL, 2023) on future tornado intensity. Two costs: it is
**CONUS-only**, and PGW is a single warming perturbation, not a scenario ensemble — so no
`n_models` in the OUTPUT-SPEC sense. Volume is large (multi-hundred TB for the full archive;
a variable subset is the only sane ingest).

---

## 5. Class 3 — loss records

These answer *"what has it cost"*, and are exposure×vulnerability records, not hazard.

- **NOAA Billion-Dollar Weather and Climate Disasters** — **retired by NOAA in May 2025**,
  frozen at December 2024. **Relaunched outside government by Climate Central as of
  2025-07-28**, which now maintains it. Severe convective storms are ~203 of the events since
  1980, about half of all billion-dollar disasters — this is the single most-cited US SCS loss
  series and its custody changed hands last year; any pipeline pointed at the old NCEI
  endpoint is pointed at a dead dataset.
- **EM-DAT** (CRED, UCLouvain) — global, free for non-commercial with registration.
- **Commercial cat models** (Verisk, Moody's RMS, Karen Clark) — the industry standard for
  SCS, licensed, not ingestible here.
- SPC's own `loss`/`closs` fields — explicitly disclaimed by SPC as "highly suspect".

---

## 6. How each class sits against our output contract

OUTPUT-SPEC requires `(decade, lat, lon)` carrying `median`, `lower_ci`, `upper_ci`,
`percentile`, `ols_slope`, `sen_slope`, `n_members`, `n_models`, per scenario.

| Class | Scenario dim | `n_models` | Fits the contract? |
|---|---|---|---|
| Occurrence records (§3) | none | n/a | **No.** Historical only; no ensemble, so CI and both slopes are undefined as specified. Could only be a baseline/validation layer, never a projection layer. |
| ERA5 environments (§4.1) | none | 1 | **No.** Single realization, observed period. Same limitation. |
| **CMIP6 environments (§4.2)** | ssp126/370/585 | 7 (data-limited) | **Yes** — this is the only candidate that can populate the contract as written. |
| AR-CHaMo (§4.3) | inherits its driver | inherits | Conditionally — it is a method applied to a driver, not a dataset to ingest. Europe-calibrated. |
| CONUS404-PGW (§4.4) | 1 perturbation | 1 | **No** as an ensemble; yes as a high-resolution sensitivity check. CONUS only. |
| Loss records (§5) | none | n/a | No. Different quantity entirely. |

**A resolution caveat harder than any we currently carry.** Our layers are 0.5° (~55 km) and
`sealevel-2b` already carries a `resolution_caveat` because coastal inundation turns on metres
of elevation. A tornado damage path is **~100 m wide and a few km long** — six orders of
magnitude below a 0.5° cell in area. Any tornado layer we ship would be a statement about the
*environment over a 55 km cell*, and the gap between that and "will this building be hit" is
not a caveat, it is the whole quantity. Under the rule in CLAUDE.md this would set
`resolution_caveat` and promote to `must_disclose` in both reports — and honestly it is worth
asking whether a screening layer at that gap is useful to a customer at all, which is a
product question, not a data question.

---

## 7. A gap in the taxonomy this review exposed

`config/hazard_taxonomy.yaml` has **20 families**. None of them is tornado or severe
convective storm. The two nearest:

- `storm-wind-extratropical` — its own note says *"The `cyclone` layer is TROPICAL cyclone
  only. Mid-latitude windstorm is a different hazard"*. Extratropical cyclone, not convective.
- `heavy-precipitation` — "Heavy precipitation, hail, snow and ice". Hail is convective, but
  the family is precipitation-typed and its blocker reads "sub-daily hazard".

So a tornado sits between two families and inside neither. Since the caveat generator emits
the "hazards not assessed" section **from this file**, a hazard with no family is not
disclosed as unassessed — it is simply absent. That is the exact failure the file's own header
warns about: *"a report that lists the hazards we assessed, and stops there, reads as though
the hazards we did NOT assess were assessed and found immaterial."*

Worth noting the family list is **ours**, so this is a boundary choice, not a compliance
defect — but CDP's storm dropdown and the EU taxonomy's wind-hazard column both distinguish
convective storms, and our `cdp_label` values are marked UNVERIFIED against the live 2026
questionnaire in any case.

---

## 8. Open decisions

Listed, not answered.

1. **What question is the layer for?** Occurrence density, environmental favourability, or
   expected loss. §3, §4 and §5 answer different ones and are not substitutes.
2. **Global or CONUS?** The only good occurrence record, the only calibrated US index (TEI)
   and the only convection-permitting archive are all US. Europe is reachable but through a
   licensed database (ESWD) and a model we would have to reimplement (AR-CHaMo). Everywhere
   else is thin.
3. **Whether a 0.5° tornado layer is a product at all**, given §6. This one probably wants
   answering before any of the others.
4. **Whether the reporting-bias confound (§3.2) is disclosable or disqualifying** for any
   occurrence-based layer, given it correlates with exposure.
5. **Licence posture** — ESWD and Tornado Archive both restrict commercial use; NTP requires
   contact for commercial. Whether this product counts as commercial use is a question for
   you, and it gates Europe and the global compilation.
6. **Whether to record the lightning product** (§1.3) as a hazard candidate in its own right —
   it is single-GCM and a thunderstorm not tornado proxy, but it is real, cheap, and currently
   unrecorded.

## 9. Repo updates this review implies, not yet made

- `config/isimip_search_catalog.yaml` — record today's listings as a receipt: the `lightning`
  sweep (closing the "NOT swept" note at line 3787) and the scoped tornado negative from §1.4
  with `verified_absent_on`. Per GUARDRAILS §11 an unrecorded negative gets re-derived; per
  the same rule this one must carry exactly what was listed.
- `config/hazard_taxonomy.yaml` — a decision on §7. Adding a family renders a new paragraph
  into every customer report, so it is a product decision, not a cleanup.

## Sources

ISIMIP file server listings, 2026-08-18 (see §1). NOAA SPC
[WCM severe weather database](https://www.spc.noaa.gov/wcm/) and
[SVRGIS](https://www.spc.noaa.gov/gis/svrgis/); NOAA NCEI
[Storm Events Database](https://www.ncei.noaa.gov/stormevents/details.jsp?type=collection).
Western University [Northern Tornadoes Project](https://www.uwo.ca/ntp/index.html) and
[NTP Open Data](https://ntpopendata-westernu.opendata.arcgis.com/). ESSL
[European Severe Weather Database](https://www.essl.org/cms/european-severe-weather-database/)
and [CHECC](https://www.essl.org/cms/checc/).
[Tornado Archive](https://tornadoarchive.com/home/) and Boyd et al.,
[BAMS 105(7), 2024](https://journals.ametsoc.org/view/journals/bams/105/7/BAMS-D-23-0123.1.xml).
[BoM Severe Storms Archive](https://www.bom.gov.au/australia/stormarchive/). Taszarek et al.,
[npj Clim. Atmos. Sci. 2021](https://www.nature.com/articles/s41612-021-00190-x) and
[J. Climate 33(23), 2020](https://journals.ametsoc.org/view/journals/clim/33/23/jcliD200346.xml).
Lepore et al., [Earth's Future 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021EF002277).
Rädler et al., [npj Clim. Atmos. Sci. 2019](https://www.nature.com/articles/s41612-019-0083-7);
[AR-CHaMo (JAMC 57(3), 2018)](https://journals.ametsoc.org/view/journals/apme/57/3/jamc-d-17-0132.1.xml).
Woods et al., [GRL 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL104796).
[CONUS404 (BAMS 104(8), 2023)](https://journals.ametsoc.org/view/journals/bams/104/8/BAMS-D-21-0326.1.xml),
[CONUS404-PGW](https://www.sciencebase.gov/catalog/item/65ff28d3d34e64ff1548df1b),
[NCAR GDEX d559000](https://gdex.ucar.edu/datasets/d559000/). ECMWF
[ERA5 daily statistics](https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=overview)
and [ERA5 documentation](https://confluence.ecmwf.int/spaces/CKB/pages/76414402/ERA5+data+documentation).
Climate Central [Billion-Dollar Disasters](https://www.climatecentral.org/climate-services/billion-dollar-disasters);
NOAA NESDIS [notice of retirement](https://www.nesdis.noaa.gov/about/documents-reports/notice-of-changes/2025-notice-of-changes/billion-dollar-weather-and-climate-disasters).

---

# 10. Review of the proposed historical-only design

**Added 2026-08-18**, second pass. Proposal under review: historical occurrence data only;
spatial aggregation to 0.1–0.25°; outputs `median`, `q25`/`q75`, and a two-tier `percentile`
(zeros → 1, non-zeros ranked 2–100); no decadal values, no trends.

All numbers below are measured from `1950-2025_actual_tornadoes.csv` (73,458 records; 73,378
inside a 20–50 °N, 60–130 °W analysis box; 80 dropped as out-of-box or null).

## 10.1 Historical-only — agreed, but the stated reason is not the binding one

Correct call, and it removes the §6 blocker: no scenario means no `n_models`, and there was
no honest way to publish the contract's ensemble fields from a single observational record.

One correction for the record, because it will be re-litigated otherwise. *"We can't calculate
decadal values or trends using only historical data"* is not quite right — 76 years is ample
for a trend, and decadal counts are literally computable (§10.3 uses them). The binding reason
is §3.2: **the trend would be a reporting-capability trend.** F0 reports up 7.7×, F2+ down 19%,
unrated from 0 to ~204/yr. Any slope fitted to this series measures the observing system.

That belongs in the file's global attributes as the reason no slope is emitted. "Not
computable" invites someone to notice it *is* computable and add it.

## 10.2 Resolution — 0.25° survives, 0.1° does not

Three independent measurements, any one of which is disqualifying for 0.1°.

**(a) Positional accuracy of the source.** Latitude quantization by era:

| era | n | slat on an exact 0.1° multiple | distinct slat values |
|---|---|---|---|
| 1950–1975 | 16,830 | **25.6%** | 1,373 |
| 1976–1999 | 23,645 | 16.8% | 1,424 |
| 2000–2025 | 32,903 | 6.0% | 20,139 |

Pre-1976 coordinates were recorded to roughly 0.1° — 1,373 distinct latitudes for 16,830
events, i.e. a county-centroid or nearest-town geocode, not a survey position. A 0.1° grid
**resolves that quantization into the map**: artificial hot spots and stripes on round
coordinates, strongest in exactly the early decades. At 0.25° the quantization is absorbed
inside the cell. This argument is independent of anything to do with tornado climatology.

**(b) At 0.1° the map shows the record, not the hazard.** 46.1% of occupied cells contain
**exactly one** tornado in 76 years (79.0% if restricted to F2+). A cell whose entire evidence
base is one event in 1963 is not a hazard estimate.

**(c) Track geometry.** Median path length is 1.0 mi, so most tornadoes are genuinely sub-cell
and point-assignment is defensible at 0.25°. But 14.3% of tracks exceed a 0.1° cell versus
4.2% at 0.25° (max 234.7 mi), and long-track tornadoes are disproportionately the violent
ones. Assigning them wholly to the touchdown cell displaces risk systematically upstream. At
0.1° that needs track rasterization — and **only 64.5% of records carry a usable end point**,
so a third of the data has no geometry to rasterize.

**0.25° also costs nothing against our existing grid**: it is a clean 2× refinement of the
0.5° ISIMIP grid every other layer uses, so it downsamples exactly for cross-layer work.

## 10.3 `median` + `q25`/`q75` — degenerate as specified. This is the finding.

Within-cell quantiles of the annual count, computed over occupied cells only (cells that have
**ever** recorded a tornado — the most favourable possible denominator):

| resolution | pooling block | q25 | **median** | q75 | q90 |
|---|---|---|---|---|---|
| 0.25° | annual | 0 | **0** | 0 | 0 |
| 0.25° | 5-year | 0 | **0** | 1 | 2 |
| 0.25° | decadal | 0 | **0** | 1 | 3 |
| 0.5° | annual | 0 | **0** | 0 | 1 |
| 0.5° | 5-year | 0 | **1** | 2 | 5 |
| 0.5° | decadal | 0 | **2** | 5 | 8 |

Zero fraction of (cell, year) pairs among occupied cells: **97.3%** at 0.1°, **91.6%** at 0.25°,
**79.9%** at 0.5°.

So with historical-only and no decadal blocking, the pool is the annual count series, and
`median`, `lower_ci` and `upper_ci` would be published as **identically (0, 0, 0) over the
entire map**. Not sparse — uniformly empty. Three of the four output variables would carry no
information, and the file would pass a structural contract test while saying nothing.

This is the documented `pooled_mean_zero_inflated` failure (CLAUDE.md: the median branch
"erased 93% of exposed land on one layer"), in its most extreme form yet measured here.

**The two proposals are in direct tension.** The first combination that yields a non-degenerate
within-cell median is **0.5° with decadal blocks** — coarser than proposed, using the blocking
that was excluded. Note the "full record" row is deliberately absent from the table above:
pooling the whole record leaves *one value per cell*, so a within-cell IQR does not exist at
all. (Quantiles across cells of the per-cell total are a different quantity — a spatial
distribution, not an uncertainty — and must not be written into `lower_ci`/`upper_ci`.)

**What the quantity probably wants to be instead.** The per-cell estimand for an event process
is a **rate**, not a count: events per 10⁴ km² per year over the full record, which is
continuous, non-degenerate, and area-normalised so it survives the latitude convergence of a
regular lat-lon grid. Then:

- `median` carries the rate (the `pooled_mean_zero_inflated` branch, declared as a deviation);
- `lower_ci`/`upper_ci` carry a **Poisson or bootstrap interval on that rate**, which is the
  quantity a customer actually needs — whether 3 events in 76 years is distinguishable from 6.
  Empirical q25/q75 cannot express this: they are 0 regardless of whether the cell saw 1 event
  or 40. Worth being explicit that **q25/q75 across years is variability, not uncertainty**,
  and on this field it is not even that.

Per CLAUDE.md this is a declared deviation: measure what the median branch would have
published (done — identically zero, table above), record it in `decadal_statistic_rationale`,
and get sign-off. The receipt here is unusually clean.

**Standard tornado climatology does something further** — Gaussian kernel smoothing of
touchdowns or tracks to a practically-significant radius (the Brooks-style approach), which
converts a point process into a continuous intensity field and is what makes published US
tornado risk maps look like fields rather than confetti. That is a smoothing-length decision
of exactly the kind CLAUDE.md says is a per-layer measurement, not a constant.

## 10.4 The two-tier percentile — right scheme, but the DOMAIN is unstated and decides everything

Zeros → 1, non-zeros ranked 2–100 is the repo's existing zero-inflated treatment and is the
correct shape. Two problems, both about the population being ranked.

**(a) The zero bucket is the majority of the map, so the percentile is nearly binary.**

| scoring domain | cells | non-zero | **zero bucket** |
|---|---|---|---|
| CONUS box @ 0.25° | 33,600 | 9,259 | **72.4%** |
| global land @ 0.25° | ~371,600 | 9,259 | **97.5%** |

On the CONUS domain, percentiles 2–100 are compressed into 27.6% of cells, so "50th
percentile" means *median among cells that have ever had a report* — not median risk. A reader
will read it the other way. This is disclosable, but it must actually be disclosed.

**(b) On a global grid the percentile ranks reporting systems, not hazard.** SPC is US-only.
Every non-US land cell would score 1 — not because tornadoes are absent but because SPC does
not observe there. Bangladesh, which has the deadliest tornadoes on record anywhere, would
render at the 1st percentile beside the Sahara. That is not a caveat one discloses; it is a
wrong answer at the first thing a reader checks, and it is the same class of error as the
`relative_baseline` reversal but without the defence that the number is correct.

So the domain has to be pinned before anything else: either

- the layer is **explicitly CONUS-only**, `assert_region_coverage()` refuses non-US sites, and
  the asset catalog never maps a non-US asset to it; or
- it is a **multi-source merge** (Tornado Archive / ESWD / NTP), which reopens the licensing
  problems in §3 and adds cross-country rating heterogeneity (F vs EF vs IF vs T scales) that
  the merge would have to reconcile.

There is no third option in which a global grid is populated from SPC.

## 10.5 The reporting-bias fix collides with the sparsity problem

The standard mitigation for §3.2 is to restrict to F2+, the reporting-stable subset. Measured
cost of doing so:

| subset | n | res | occupied cells | once-only cells | full-record count/cell (median) |
|---|---|---|---|---|---|
| all reports | 73,378 | 0.25° | 9,259 (27.6%) | 16.0% | 6 |
| **F2+ only** | **13,400** | 0.25° | 5,430 (16.2%) | **38.5%** | 2 |
| all reports | 73,378 | 0.5° | 2,911 (34.7%) | 10.5% | 19 |
| **F2+ only** | **13,400** | 0.5° | 2,038 (24.3%) | 17.1% | 5 |

F2+ discards **82% of the records**. At 0.25° it leaves 38.5% of occupied cells resting on a
single event. So the fix for the confound that makes the layer misleading is the same knob
that makes the layer statistically thin — they pull opposite ways, and the resolution choice
sits between them. This is the real design trade, and it is worth deciding deliberately rather
than discovering after processing.

Note also that F2+ is itself degrading as a filter: unrated (`mag = -9`) reports now run
~204/yr, ~15% of the modern record, and cannot be filtered on magnitude at all (§3.2).

## 10.6 The output contract needs an explicit decision, not a discovery at test time

`scripts/test_shared_baseline.py` line 24 requires
`["median", "lower_ci", "upper_ci", "percentile", ...]` on a `decade` dimension, and CLAUDE.md
adds `ols_slope`, `sen_slope`, `n_members`, `n_models`. The proposal drops the decade
dimension, both slopes and both ensemble counts — so this is a **new output class**, and the
verifier will reject it.

Two ways, both needing your call: a second contract for observational layers (with its own
verifier path), or the existing shape with a single degenerate `decade` coordinate
(e.g. `1950-2025`) and slopes emitted as NaN. The second keeps one code path but publishes
variables that exist only to be empty, which is the pattern CLAUDE.md warns against elsewhere.

Whichever it is, `n_members`/`n_models` have no meaning here and should not be faked to 1.

## 10.7 Summary of the review

| Proposal | Verdict |
|---|---|
| Historical only | **Sound.** Fix the stated reason (reporting-capability trend, not "not computable"). |
| 0.25° | **Sound**, and a clean 2× refinement of the existing 0.5° grid. |
| 0.1° | **Not supported** — below the source's own positional accuracy (25.6% of pre-1976 records on exact 0.1° multiples), 46% of cells single-event, 14.3% of tracks longer than a cell. |
| `median` + `q25`/`q75` | **Degenerate as specified — identically (0,0,0) everywhere.** Needs a rate-based quantity with a Poisson/bootstrap interval, declared as a `pooled_mean_zero_inflated` deviation. |
| Two-tier percentile | **Right scheme, undecided domain.** CONUS-only or a licensed multi-source merge; a global SPC-only grid ranks reporting systems. |
| No trends | **Sound**, for the reason in 10.1. |

Open, in the order they gate each other: (1) scoring domain — CONUS-only or merge; (2) all
reports vs F2+; (3) rate vs count and the uncertainty estimator; (4) smoothing length, if any;
(5) which output contract this layer answers to.

---

# 11. Build record — CONUS historical tornado hazard, 2026-08-18

Decision taken: **CONUS-only, historical/current-only**, after establishing that no
licence-clean global occurrence alternative exists (§11.1). Built in this session.

## 11.1 The global question, closed

| Candidate | Coverage | Record depth | Blocker |
|---|---|---|---|
| **ESWD** (ESSL) | Europe + Mediterranean | **9,563 tornadoes, 1800–2014**; ~242/yr recently | Commercial use not allowed; fee for organisations. Severe underreporting in the Mediterranean and eastern Europe; reports rose sharply 1995–2006 from collection effort, not weather |
| **Tornado Archive** | global, >100,000 | varies by country | "No commercial use without written permission"; inherits ESWD's terms. Site returned **HTTP 403** to every automated access attempt on 2026-08-18 (4 URLs) — bulk route remains **UNVERIFIED** |
| NTP (Canada), TORRO (UK), JMA, BoM | single-country | thin | Separate ingests, incompatible rating scales (F / EF / IF / T); NTP needs contact for commercial use |

Europe's entire 214-year pan-continental record is smaller than **seven years** of US
reporting (73,458 records since 1950). There is no global occurrence dataset with detection
comparable to SPC's, and the two that approach global coverage both restrict commercial use.

The only genuinely global, licence-clean, historical option is a **reanalysis environment
climatology** (ERA5 CAPE × shear × SRH, §4.1) — a different quantity, not occurrence, and
recorded here so the option is not lost.

## 11.2 What was built

`scripts/download_tornado_spc.py` → `scripts/process_tornado_spc.py`.
Raw 12 MB, processed **1.4 MB across 8 layers**. Every raw file carries a `.json` sidecar
with `source_url`, `sha256`, byte size and retrieval date, because `data/` is ephemeral.

Grid **0.25°** on exact quarter-degree edges (nests 4-to-1 inside the 0.5° ISIMIP grid),
CONUS box 24–50 °N / 66–125 °W, **13,377 in-mask cells** of 24,544 (54.5%). Mask is
Natural Earth 1:50m USA by cell-centre containment — `scripts/utils/natural_earth.py` could
not be used because it imports **geopandas at module level and geopandas is not installed**
in this venv (shapely is). Outside the mask every field is **NaN, never 0**.

Counting is by **track**, not touchdown: each tornado increments every cell its damage path
crosses, so `n_events` sums to more than the tornado count by design. 35.5% of records carry
no end point and fall back to touchdown; `--geometry touchdown` forces that everywhere.

Published per cell: `median` (crossings per 10⁴ km² per year), `lower_ci` / `upper_ci`
(quartiles of the Gamma(k+½, T) rate posterior), `percentile` (two-tier), plus `n_events`
and `conus_mask` as diagnostics. Cell area varies **42%** across the box, so the rate is
area-normalised — without it Texas would outrank Iowa for being further south.

**Contract**: `observational-historical-v1`, declared in the file's own attributes as **not**
the OUTPUT-SPEC decadal contract. No `decade`, no slopes, and `n_members`/`n_models` are
**absent rather than set to 1** — a 1 would read as a thin ensemble instead of a different
kind of product. `scripts/test_shared_baseline.py` will reject these files; that is correct
and needs a verifier decision, not a workaround.

**One defect found and fixed during the build.** For a zero-count cell the Jeffreys posterior
puts mass above zero, so `lower_ci` came out at 0.010 against a `median` of exactly 0.000 —
violating the `lower_ci <= median <= upper_ci` invariant that `test_shared_baseline.py:172`
enforces. Observing nothing bounds a rate from above, not below, so the lower bound is now
clipped to the point estimate, giving `[0, upper]` on unoccupied cells with `upper_ci`
carrying the information. Only k=0 is affected; for k≥1 the posterior quartile already sits
below k/T. An assertion now guards it.

## 11.3 The ladder — the threshold trade, measured across 8 layers

Cells occupied of 13,377 in-mask, and the resulting zero-tier share of the percentile:

| rung | records (full) | occupied (full) | zero tier | records (1996–) | occupied (1996–) | zero tier |
|---|---|---|---|---|---|---|
| `all` | 73,378 | 9,235 (69.0%) | 31.0% | 37,983 | 8,135 (60.8%) | 39.2% |
| `f1plus` | 38,565 | 7,974 (59.6%) | 40.4% | 15,948 | 6,224 (46.5%) | 53.5% |
| `f2plus` | 13,400 | 6,390 (47.8%) | 52.2% | 4,125 | 3,866 (28.9%) | 71.1% |
| `f3plus` | 3,325 | 4,061 (30.4%) | 69.6% | 1,002 | 1,826 (13.7%) | 86.4% |

Every step that buys reporting homogeneity costs occupancy. `f3plus` on the modern window
leaves 86% of CONUS in the zero tier; `all` on the full record fills 69% of it with a series
that is partly an artefact of observation. The choice lives on this table.

## 11.4 The reference-site check did its job

Checked at seven sites before accepting the output. Phoenix, AZ read the **64.5th percentile**
on `all` / full record — implausible for the Sonoran desert. Inspection of the raw rows shows
the 12 records are **real**: all (E)F0–F1, path lengths 0.2–1.2 mi, and Arizona statewide is
196 F0 / 76 F1 / 16 F2 / 4 F3. So it is not a bug — it is the reporting-density confound
made visible, a large metro accumulating weak reports in a low-tornado state.

Phoenix across the ladder: **64.5th → 45.1st → 1st → 1st** percentile
(`all` → `f1plus` → `f2plus` → `f3plus`), and 4.4th on `all` restricted to 1996–2025.
Moore, OK holds 99.9–100th on every rung. That is the clearest single demonstration that the
`all` rung is contaminated and the F-rungs are not.

Second finding from the same check: Caribou, ME reads the **51.0th percentile** on `all` —
arithmetically correct (7 crossings, against many occupied cells holding 1–2) and a perfect
illustration of why *percentile 50 is not median risk*. It is the median among cells that
have ever recorded a tornado. That inversion is why the percentile note is a must-disclose
caveat rather than a footnote.

## 11.5 Still open

1. **Which rung and which window ship** — §11.3 is the decision table; no default was chosen.
2. **The verifier**: a second contract for observational layers, or an exemption path.
3. **`config/layer_registry.yaml` and `config/hazard_taxonomy.yaml` are untouched.** The
   taxonomy gap from §7 is now concrete — a shipped tornado layer with no hazard family.
4. **`qa_reviewed_on` is null in all 8 files.** Nobody has viewed the maps; the reference-site
   check is not that review.
5. **Smoothing not applied.** Kernel smoothing is standard in tornado climatology and the
   decay length is a per-layer measurement, not a constant — unmeasured here, so left off.

---

# 12. Follow-through, 2026-08-18

## 12.1 Enumeration receipts recorded

`config/isimip_search_catalog.yaml`:

- **`negative_results.tornado`** — the scoped absence with `verified_absent_on` naming every
  URL listed, an explicit `scope_of_this_negative` block for what was *not* listed, and a
  `structural_reason` recording that this negative will not be overturned by a new
  publication: the shear term needs the atmosphere in the vertical, and ISIMIP forcing has
  11 surface variables and nothing aloft.
- **`search_results.lightning`** — the first-ever enumeration of `InputData/climate/lightning/`
  in 3a and 3b, with the sidecar-read metadata, the single-GCM blocker, and an explicit
  `what_it_is_and_is_not` stating it is a convective-**activity** proxy and must not stand in
  for the severe-convective family.
- The `repository_structure_cache` "NOT swept" line was corrected — `lightning` moved from
  unswept to swept with its date.

## 12.2 The taxonomy gap closed

`config/hazard_taxonomy.yaml` gains **`storm-convective`** ("Tornado and severe convective
storm"), acute — the 21st family. Before this, a tornado was the one hazard that could not
even be *disclosed* as unassessed: it fell between `storm-wind-extratropical` (explicitly
extratropical only) and `heavy-precipitation` (precipitation-typed), and a family absent from
that file is absent from the mandatory "hazards not assessed" section.

`covered_by: []` **on purpose**, following the `cold-frost` precedent — no rung is chosen,
`qa_reviewed_on` is null everywhere, and the files answer a different contract. The
`materiality_note` says so, and says what must be true before coverage is claimed.
Verified: `customer_note` carries no repo paths, no `UNVERIFIED`, no internal vocabulary;
taxonomy→registry linkage has zero violations.

## 12.3 A second contract, with its own verifier

`scripts/test_observational_baseline.py` — **152 checks across the 8 layers, all passing.**
The shared verifier was left untouched: it guards 23 shipped layers and its strictness is the
product. Beyond the ordinary shape checks it encodes three rules this repository learned the
hard way:

- **Decadal-contract variables must be ABSENT, not faked.** `n_members: 1` reads as a thin
  ensemble; absence reads as a different kind of product.
- **Two-tier consistency both ways** — zero cells must be at percentile 1, *and* percentile-1
  cells must be zero. Checking one direction only would miss a mis-tiered field.
- **No negative caveats.** A caveat attribute is promoted on being non-empty, so
  `resolution_caveat: "none"` publishes a must-disclose caveat whose body reads "none". The
  verifier fails any caveat opening with a negation, catching that at build time rather than
  in a filing.

## 12.4 QA maps exist, so the review can happen

`scripts/generate_tornado_qa.py` → `reports/tornado-qa/tornado-qa-{full,from1996}.html`,
3.6 MB each, 12 panels per page (4 rungs × rate / percentile / crossing count), plus the rung
summary and the caveats rendered from each file's own attributes.

Deliberately **not** built into `generate_maps.py`: that tool selects on `ds.decade` and
builds Trend and Members tabs from variables this contract does not have, and teaching it a
second contract would put a branch through the renderer for 23 shipped layers to serve 8
without a decade axis.

The page opens with what to look for — whether the maximum *moves* between rungs, round-
coordinate stripes, whether the percentile panel reads as a hazard field or a coastline, and
state-border discontinuities (an edge on a state line is an NWS-office artefact, never
weather). It supports the review; it is not the review, and `qa_reviewed_on` stays null.

## 12.5 Smoothing — measured, and a first criterion discarded

`scripts/measure_tornado_smoothing.py`. No ensemble exists to halve, so the *record* is
halved: **odd years vs even years**, deliberately not early-vs-late, which would confound the
test with the reporting trend and measure the observing system instead of sampling noise.

**The obvious criterion is degenerate and was discarded after measuring it.** "Smooth both
halves, correlate, take the maximum" rose *monotonically* with σ on all four rungs out to the
largest length tested (ρ 0.81 → 0.99 on `all`; 0.38 → 0.98 on `f3plus`). Blurring any two
fields toward a constant makes them agree, so that rule always returns the largest σ offered
and selects nothing. Replaced with **held-out Poisson predictive likelihood** — smooth one
half, score how well it predicts the *other half's actual counts* — which is cross-validation
and has a genuine interior optimum.

| rung | ρ between halves, **unsmoothed** | σ maximising held-out likelihood |
|---|---|---|
| `all` | 0.812 | **25 km** |
| `f1plus` | 0.765 | **25 km** |
| `f2plus` | 0.612 | **60 km** |
| `f3plus` | **0.378** | **90 km** |

Two findings. **Smoothing is warranted on every rung** — σ=0 scores worse than every smoothed
option on all four, and on the three F-rungs its held-out likelihood is outright negative.
And **the sparser the rung, the longer the decay length the data supports**, which is what
should happen: fewer events need more spatial pooling before the field is estimable. The
`f3plus` unsmoothed halves correlate at 0.378 — that field barely reproduces itself at 0.25°.

**Nothing was smoothed.** The decay length is rung-dependent, which is itself a decision (a
single σ across rungs would be wrong for at least two of them), and adopting it changes every
published number.

## 12.6 Deliberately not done

- **`config/layer_registry.yaml` untouched.** `LayerSpec(**spec)` puts an entry straight into
  the delivery path, and `generate_customer_delivery.py` reads decades and slopes these files
  do not have. Registering them would break delivery for everything, not just tornado.
- **`qa_reviewed_on` still null in all 8 files.** It is a human's date. The maps now exist;
  the reference-site check is not that review.
- **No rung and no window chosen** — §11.3 is the table.
- **No smoothing applied** — §12.5 is the measurement.

---

# 13. External review and the corrections it forced

An adversarial review was run against §10–§12 and all five scripts. It found real defects.
Below: what was accepted and fixed, what was tested and did **not** hold up, and what is
correct in principle but currently inert. Everything here was re-measured, not conceded.

## 13.1 The central criticism, accepted

**"The product calls rasterised line intersections an area-exposed Poisson point process."**

Half right, and the half that is right matters. Two separate claims were tangled together:

- **The Poisson *form* survives track counting.** Thinning a Poisson process of tornado
  objects down to those intersecting one cell is still Poisson, so a per-cell marginal
  interval is not invalidated. Measured independently: the variance-to-mean ratio of annual
  per-cell counts is **0.99 at the median cell** (p90 1.70; 16.4% of cells above 1.5 for all
  reports, 6.2% for F2+). The marginal Poisson assumption holds well at the typical cell and
  understates the interval in the clustered tail. Now recorded as `overdispersion_measured`.
- **The area *normalisation* does not survive.** Expected crossings scale with roughly
  (cell size + path length) × cell size, not with area, so "per 10⁴ km²" is not a
  resolution-invariant intensity. **Measured: aggregating the four 0.25° children onto their
  0.5° parent reads +7.9% at the median cell (p90 +23%), and total crossings inflate 1.123×
  from 0.5° to 0.25°.**

That measurement also **falsifies a claim in §10.2 of this document**: that 0.25° "downsamples
exactly" onto the 0.5° ISIMIP grid. The *grid* nests exactly; the *values* do not aggregate.
Corrected in the text and in a new `rate_is_resolution_dependent` attribute that says: never
resample this layer, recompute it.

## 13.2 The estimator was incoherent — fixed

`median` was the MLE k/T while `lower_ci`/`upper_ci` were Gamma(k+½, T) posterior quartiles,
with the lower bound **clipped** to the MLE so the triple would order. The review named that
correctly as a fudge over a real incoherence: an MLE on the boundary at k=0 cannot be
bracketed by quantiles of a continuous posterior, and the published interval was therefore
not the Gamma quartile interval it was documented as.

Now all three are percentiles of **one** posterior — no clipping, ordering holds by
construction, and `--estimator mle` reproduces the old behaviour for comparison. A cell with
no report no longer claims a rate of exactly zero; it gets the small positive rate that
seeing nothing in 76 years implies, which is the more honest statement.

**This changed the zero tier's key.** With a positive rate everywhere, a `rate == 0` test
selects nothing — so the tier is now keyed on the **observed count**. The verifier check for
it would have passed *vacuously* under the new estimator, which is worse than no check; it
now keys on `n_events` too, and gained a monotonicity check and a `[2, 100]` range check.

## 13.3 Tested and did NOT hold up

**"The smoothing selection is substantially controlled by an arbitrary floor."** It is not.
Sweeping the floor across four orders of magnitude (mean/10 → mean/100,000) moves the selected
σ **not at all** — 25 km, 60 km, 90 km on `all`, `f2plus`, `f3plus` at every floor. The
concern was reasonable and the measurement refutes it.

But testing it surfaced something the review did not: the *criterion* matters even though the
floor does not. A **floor-free negative-binomial predictive** — smooth the training counts and
exposure, integrate the Gamma posterior out — is the principled version, and it selects
**shorter** lengths, because the plug-in ignores estimation uncertainty and compensates with
extra smoothing. It is now the primary criterion in the script.

| rung | σ, plug-in Poisson | σ, **NB predictive (primary)** |
|---|---|---|
| `all` | 25 km | **15 km** |
| `f2plus` | 60 km | **40 km** |
| `f3plus` | 90 km | **40 km** |

**Two corrections to §12.5 follow.** The σ values there were the plug-in ones and are revised
above. And the claim that *"the sparser the rung, the longer the decay length the data
supports"* **does not survive** — under the floor-free criterion `f2plus` and `f3plus` both
land at 40 km. That apparent pattern was an artefact of the plug-in criterion.

What does survive: **σ = 0 loses on every rung under both criteria**, so smoothing is
warranted. That conclusion was challenged as floor-dependent and is now established without a
floor at all.

## 13.4 Correct in principle, currently inert — fixed anyway

Three findings were real code defects with no effect on any published number. Verified before
fixing rather than assumed:

- **Exposure years from `yr.nunique()`** counted only years containing a qualifying report.
  Measured: all four rungs have 38/38 in both halves, so it was inert — but it would inflate
  the rate for any thinner subset. Now calendar-derived.
- **NaN endpoints would crash** `int(np.ceil(span/…))`. Measured: zero non-finite coordinates
  in the current file. Now guarded.
- **Percentile tie endpoints.** Measured: minimum occupied score 2.005 on `all`, 2.000 on the
  F-rungs; maximum 100.000 everywhere. Negligible, now documented rather than claimed exact.

Also fixed: `all(genexp)` in the downloader short-circuited on first failure, contradicting
its own "report and continue" comment; `nan_aware_gaussian` was mask-aware but not NaN-aware,
so its name was false; the verifier's out-of-mask check was labelled "every field" but tested
only `median`; and the verifier could excuse a NaN `lower_ci` under a finite `median` by
intersecting the masks before checking ordering.

One suspicion of mine was checked and **dismissed**: the `roughness()` mask/diff alignment is
correct as written.

## 13.5 Overclaims withdrawn

- *"No crossed cell is skipped"* at res/8 sampling — **false**. Dense point sampling is not
  raster traversal; a line can clip a cell corner between samples. Sampling tightened to
  res/16 and the claim replaced with an accurate one, naming supercover traversal as the
  exact fix.
- *"The probability a building is struck ≈ footprint / cell area"* — **wrong**, and it was in
  a `resolution_caveat` bound for customer reports. Strike probability depends on swath width,
  path length and orientation, none of which this layer models. Replaced with an explicit
  instruction not to convert a cell value into a strike probability.
- *"Odd/even is immune to the reporting trend"* — too strong, and it contradicted the same
  file's own note that shared bias is invisible to the test. It controls the secular trend;
  both halves still carry identical spatial bias.
- *"The F-rungs are not contaminated"* — one Phoenix/Moore comparison cannot establish that
  CONUS-wide. Narrowed to what was measured: at that site, on that ladder.
- *"Reporting-stable subset"* for F2+ — the literature's term, but a falling aggregate count
  does not establish spatially stable detection. Flagged as inherited, not demonstrated.

## 13.6 Open, and now better specified

The estimand itself is the real outstanding question, and no amount of interval arithmetic
settles it. The honest options: publish **touchdown counts**, for which area × years is a
sound exposure; publish **track intersections per cell per year without area normalisation**
and let it be explicitly resolution-specific; or model a **line/swath process** with path
length, width and orientation, which is the only route to a defensible strike probability.
The current layer is the second dressed as the first, and `rate_is_resolution_dependent` now
says so in the file.

---

# 14. Routing into customer deliveries

**Added 2026-08-18.** The layers now resolve for US sites through the normal delivery path,
and are masked — not scored low — everywhere else.

## 14.1 Proven end to end, not asserted

A delivery was run over seven warehouse sites spanning the tornado maximum, two quiet US
regions, two non-CONUS US states and one European site:

| site | value (per 10⁴ km²/yr) | percentile | `data_status` |
|---|---|---|---|
| Moore, OK | 4.89 | **99.88** | OK |
| Tuscaloosa, AL | 1.99 | 94.08 | OK |
| Seattle, WA | 0.24 | 18.99 | OK |
| Phoenix, AZ | 0.14 | 8.13 | OK |
| Honolulu, HI | — | — | **OFF_LAYER_MASK** |
| Anchorage, AK | — | — | **OFF_LAYER_MASK** |
| Frankfurt, DE | — | — | **OFF_LAYER_MASK** |

`OFF_LAYER_MASK` is the delivery's existing phrase for *"your site is on modelled land and
this layer does not model it"* — exactly right here, and structurally different from a low
score. `scripts/test_customer_delivery.py` passes: **2,948 checks, 2,380 metric values
independently recomputed.**

## 14.2 The grid had to become global — a regional grid is not deliverable

The first build published a regional 24–50 °N / 66–125 °W grid. It crashed the delivery:
`spatial_extract.extract_by_point` **raises** when a point falls outside the grid extent, so
the Honolulu site aborted the whole run rather than being reported as unmodelled.

That guard is correct and was left strict — on a global layer an off-grid point means bad
coordinates, which should fail loudly. The fix is that the layer is now published on the
**global 0.25° grid**, with 13,377 CONUS cells carrying data and everything else NaN. Cost
on disk is nil (**192 KB** per layer — an almost-all-NaN array compresses away), and the
mask now carries the meaning instead of the grid extent. Values are unchanged: Moore reads
13.608 before and after.

## 14.3 Registry and catalog

Four rungs registered. **`tornado-f2plus` is `preferred`**; `all`, `f1plus`, `f3plus` are
`alternate`, reversible by moving one `status` field. The reason is §11.4: report density
tracks population, roads and radar — the same thing the exposure term measures — and Phoenix
runs 64.4th → 45.0 → 1st percentile across the ladder while Moore holds 99.9–100 on every
rung. The weak-tornado record flatters a dense metro and under-scores a remote industrial
site, which is the wrong error for this product.

Added to `warehouse`, `office`, `data center` and `agricultural land`. **Not** added to
`timber land`, and that omission is recorded in the catalog as a decision rather than left
to look like an oversight: a ~100 m path against a holding measured in thousands of hectares
is immaterial beside wildfire and drought.

## 14.4 Five shared-code changes, each with a reason

- **`spatial_extract.as_period_dataset()`** — lifts a `(lat, lon)` observational layer onto a
  single-period `decade` axis **at read time**. The file stays 2-D on disk: writing a fake
  time axis onto data with no time information is the drift this repo forbids elsewhere.
  Called from three sites (extraction, domain mask, delivery plan). Projected layers are
  returned untouched.
- **The Climate Score now excludes non-forcing scenarios.** It is a cross-*forcing*
  comparison; `observed` has no pathway, and letting it default into "medium" raised that one
  tier and would have read as a forcing effect. `viz_common.is_forcing_scenario()` declares
  it; the delivery and its independent verifier apply it separately.
- **The dashboard says "observational", not "missing scenario".** It previously warned that
  the layer "has no high/low-tier scenario", which sends someone looking for files that do
  not exist.
- **`--measure-slopes` skips slope-less layers** instead of raising `KeyError`.
- **The verifier's independent recompute now derives geometry per layer.** It hardcoded
  `SIGMA = 0.25` / `SEARCH_RADIUS = 0.5` — the values the rule gives on a 0.5° grid, correct
  while every layer was 0.5°. Against the 0.25° tornado layer it reported a passthrough
  "mismatch" that was really the *test* using a window twice the layer's own. It now applies
  the rule (radius = 1 cell, σ = ½ cell) reimplemented, not imported, so it stays independent.
  **This was a latent bug that would have mis-verified any future fine-resolution layer.**

Two checks were made stricter rather than merely tolerant of the new shape: a metric whose
source variable is absent must be delivered as **NaN** (a number there would mean the
pipeline invented it), and a layer with no `baseline_decade` must carry **no finite slopes**
(so a projected layer cannot pass by omitting its baseline).

## 14.5 Unrelated finding

`deliveries/example-forestry-co/20260813` now fails its verifier with **11 × "source changed
since delivery"** — its conifer-npp, cyclone, drought and wildfire sources were reprocessed
after it was built. All 2,736 of its passthrough values and 3,304 checks still pass, so this
is checksum drift, not a defect; the delivery is simply stale against its own manifest. Not
caused by this work and not fixed here.
