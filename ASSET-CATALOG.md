# Asset Catalog & Customer Delivery — reference

> **The `/customer-delivery` skill is the entry point, not this file.** It owns the
> four-stage pipeline and the conversational protocol; this document is the reference it
> points at. Read the skill first, come here for detail.

How a customer's list of sites becomes a deterministic CSV extract and its QA dashboard:

| File | Holds |
|---|---|
| [`config/asset_catalog.yaml`](config/asset_catalog.yaml) | asset type → which hazard layers it is exposed to |
| [`config/layer_registry.yaml`](config/layer_registry.yaml) | layer id → where it lives on disk, and which slope to read |
| `scripts/generate_customer_delivery.py` | the driver — extraction **and** dashboard |
| `scripts/utils/delivery.py` | extraction, star-schema writing, stage recording |
| `scripts/generate_delivery_dashboard.py` | the dashboard, also importable as `build_dashboard()` |
| `scripts/test_customer_delivery.py` | the verifier |

This is **stage 2**. It produces numbers, provenance and a QA dashboard — never scores,
rankings or narrative. Stage 1 (inputs) is owned by the skill; **stages 4 (caveats) and 3
(compliance and bespoke reports) are built**, and their reference documentation is
[docs/reporting/](docs/reporting/README.md). Note the ordering: caveats runs BEFORE the
reports, because the caveat set is an input to them.

---

## The three rules

### 1. A run never extracts silently

Planning is the default. `generate_customer_delivery.py` resolves every asset to its layers,
prints the mapping with each layer's units, scenarios, row count and caveats, and **stops**.
Extraction requires `--run`.

This is enforced in code rather than documented as a habit, because the requirement is that
the specific variables resolved for an asset are shown to the user *before* the run — and a
flag is the only way to make that a property of the tool instead of a property of whoever is
driving it.

```bash
# resolve and show — touches no data
python scripts/generate_customer_delivery.py \
    --customer "Acme Timber" --input location-analyses/acme-sites.csv

# after the mapping is confirmed
python scripts/generate_customer_delivery.py \
    --customer "Acme Timber" --input location-analyses/acme-sites.csv --run
```

### 2. An unknown asset type is an error, not a default

If `Asset_Type` is not in the catalog (by name or alias), the run fails and names the known
types. It does not guess, and it does not fall back to a generic hazard bundle. Work the
right layers out with the user, then add the entry.

Getting this wrong is worse than failing: a silently defaulted bundle delivers a hazard with
no transmission channel to the asset, and the customer has no way to know it was a guess.

### 3. What a run teaches goes back into the catalog

A layer that should not have been in a bundle, one that was missing, a customer's stated
preference, a new asset class — all of it is written back to `config/asset_catalog.yaml` in
the same session it was learned, with `confirmed_on` set to the date the user approved it.

An entry with `confirmed_on: null` has **not** been reviewed by anyone. `--list-layers`
prints that status explicitly. The catalog seeded on 2026-08-12 is entirely unconfirmed.

---

## Input format

One row per **location-asset combination**. A location with three assets is three rows;
they collapse to one `location_id` and get three `asset_id`s.

| Column | Required | Notes |
|---|---|---|
| `Location` | yes | Site name. Same name + coordinates ⇒ same `location_id`. |
| `Lat`, `Lon` | yes | Decimal degrees. Longitude is normalized to [-180, 180]. |
| `Asset_Type` | yes | Matched case-insensitively against catalog names and aliases. |
| `Sub_Asset_Unit` | no | Free text, e.g. species or facility class. |
| `Country`, `State`, `City`, `Region`, `Subregion` | no | Carried into `locations.csv`. |
| `Layers` | no | Semicolon-separated layer ids overriding the catalog for that row. |
| `Coord_Source` | no | `supplied` (default), or how a missing coordinate was derived. |
| `Asset_Value` | no | Carrying amount, a bare number. Unlocks the monetary half of the disclosure metrics. |
| `Currency` | no | Required whenever `Asset_Value` is given. Mixed currencies are refused — an FX rate is a financial assumption with its own date and source. |
| `Valuation_Date` | no | When the value was struck. |
| `Value_Basis` | no | book / insured / market / replacement. **Matters as much as the figure** — the four give four different answers. |

`Coord_Source` lands in `locations.csv`. Stage 1 of the workflow permits deriving a missing
lat/lon, but a derived coordinate and a surveyed one must never be indistinguishable in the
delivered file — and the ~1° extraction footprint makes a city centroid for a large estate a
real approximation, not a rounding.

`--layers "a;b"` overrides every row. Overrides are recorded per asset in
`assets.csv:layer_source` (`catalog`, `row-override`, `cli-override`) so a delivery always
says where its layer selection came from.

## Output — a normalized star schema

`deliveries/{customer-slug}/{YYYYMMDD}/`

| File | Grain |
|---|---|
| `locations.csv` | one row per distinct site |
| `assets.csv` | one row per location-asset combination |
| `layers.csv` | one row per hazard layer |
| `values.csv` | `asset_id` × `layer_id` × `scenario` × `decade` |
| `climate_score.csv` | `asset_id` × `scenario_tier` × `decade` |
| `manifest.json` | source files with mtimes and sizes, extraction parameters, row counts |
| `README.md` | how to read the delivery, generated per delivery |
| `dashboard.html` | the QA dashboard — built by the same command; skipping it marks the delivery incomplete |
| `input_locations.csv` | verbatim copy of the submitted input |

`manifest.json` also carries a `stages` block recording which of the four delivery stages
have run, so a half-finished delivery is visible in its own folder.

Location metadata is written once, layer metadata once; `values.csv` is keys plus
measurements. Nothing derivable from another column is stored: no hazard bands (compute
from `percentile`), no long-term-trend column (it is the final decade's slope), no
significance columns (no such statistic exists under this contract).

`values.csv` columns: `value`, `lower_ci`, `upper_ci`, `percentile`, `ols_slope`,
`sen_slope`, `slopes_agree`, `n_members`, `n_models`, `data_status`.

`data_status` is `OK`, `OFF_LAYER_MASK` (on modelled land, but this layer does not cover the
site — e.g. no conifer stand present), or `OUTSIDE_DOMAIN` (offshore or off-grid). Rows are
never dropped; a site with no data says so.

The domain that separates those two is the union of finite cells across **every registry
layer**, not across the delivery's layers. Scoping it to the delivery would make a status
depend on what else the customer ordered: a conifer-only delivery for an Amazon site would
report `OUTSIDE_DOMAIN` — "your site is offshore" — when the truth is that no conifer stand
is modelled on perfectly good land.

## The Climate Score

`climate_score.csv` — the unweighted mean of `percentile` across an asset's hazards, for one
forcing tier and one decade. Higher = higher aggregate physical climate risk.

Percentile is the **only** cross-hazard comparable axis: `value` is in native units that
differ per layer (g C m⁻² yr⁻¹, a dimensionless fraction, a percent), so averaging values
across hazards is meaningless arithmetic. Percentiles are also already oriented so 100 is
worst on *every* layer, `higher_is_better` ones included — which is what makes the mean
legitimate.

### It is keyed on a forcing tier, not a scenario code

This is forced by the data, not chosen. No native ISIMIP code spans both rounds. Measured on
a timber asset carrying three hazards:

| key | hazards it can see |
|---|---|
| `rcp26` | 1 of 3 (conifer-npp) |
| `ssp126` | 2 of 3 (drought-3b, wildfire) |
| **tier `low`** | **3 of 3** |

A score keyed on `ssp126` would average two hazards and be labelled "across all hazards."
The `scenarios` column records exactly which native codes contributed, so the harmonization
is visible rather than implied. **RCP and SSP tiers are only approximately comparable** —
different scenario families from different CMIP generations — and any narrative must say so.

This is the one place harmonization enters the CSV; `values.csv` still carries native codes
only.

### Two ways it can lie, and the guards

**Incomplete hazard coverage.** A hazard that does not cover a site is excluded, not counted
as zero, and ISIMIP3b layers have no 2010s panel. So the 2010s score rests on 1 hazard where
the 2020s rests on 3 — a 22 → 42 jump that reads as risk doubling when only the hazard set
changed. `n_hazards` exists for this; **two scores with different `n_hazards` are not like
for like.**

**Uneven tier coverage across layers.** `cyclone` publishes rcp26/rcp60 and *no rcp85*, so
every cyclone-carrying asset is unscoreable at the high tier. Averaging each tier over
whatever assets it happens to have makes the tier lines describe different portfolios.
Measured before the fix: the high-tier 2020s portfolio mean read **39.9 against 42.1** for
low and medium, which is impossible on a common basis — the shared 2020s panel is
bit-identical across scenarios, so all three tiers *must* be equal there.

Any chart comparing across tiers or decades therefore restricts to a **balanced panel** —
assets complete in every cell of that grid — and states what it dropped. On the example
delivery that is 2 of 6 assets, and the baseline then reads 39.88 in all three tiers.

This applies at **every** level of rollup, not just the portfolio. It bit a second time at
location level: Shasta holds a timber asset and a warehouse, and the warehouse carries
cyclone, so at the high tier only the timber asset survived — the 2020s read 62.3 against
51.7 for low and medium. Same fix, same invariant.

`test_customer_delivery.py` enforces both: it recomputes every score from `values.csv` with
an independent two-stage mean, and asserts the baseline-decade score is identical across
tiers for any asset with the same hazard set.

### Hazards are weighted equally

There is no materiality weighting. What keeps an irrelevant hazard out of an asset's average
is the **asset catalog**, not the arithmetic — which is another reason entry 2 of the three
rules (an unknown asset is an error, never a default) matters.

## The QA/QC dashboard

```bash
python scripts/generate_delivery_dashboard.py deliveries/{customer}/{date}
```

Writes `dashboard.html` into the delivery folder — one self-contained page (~110 KB for a
6-asset delivery) with five views over the same filtered slice:

The page header carries a **build stamp** (an 8-hex hash of the payload plus the page
logic). Confirm what is on disk with:

```bash
python scripts/generate_delivery_dashboard.py deliveries/{customer}/{date} --check-stamp
```

If the stamp in the browser differs, the page is cached — hard-reload. A cached dashboard is
indistinguishable from a correctly regenerated one by eye, and that ambiguity sent us
chasing two phantom bugs before the stamp existed.

| View | Shows |
|---|---|
| **Overview** | Climate Score stat tiles — portfolio score, change vs baseline, highest-risk location, coverage |
| **Map** | one marker per site; metric toggle, **Climate Score selected by default** |
| **Summaries** | mean risk percentile by hazard and by asset class, a portfolio band histogram, and Climate Score over time beside it — toggle between the Score and any single hazard's percentile, on the same 1–100 axis |
| **Time series** | one location at a time; the Climate Score panel first, then one panel per hazard, in a single grid. Toggle switches the hazard panels between Percentile and Value |
| **Table** | every filtered row, sortable, with a text search |

Climate Score leads because it is the portfolio-level answer and every other metric is a
component of it. The overview is stat tiles rather than a chart — the story is one number.

**One filter row, pinned to the top** so it stays reachable while reading the charts below.
It scopes everything — forcing tier, hazard, asset type, decade — with three stated
exceptions, each because varying that dimension *is* the view's purpose:

| View | Ignores | Why |
|---|---|---|
| **Every time series** | tier **and** decade | seeing the full scenario spread across time is the point of the chart; filtering to one tier collapses it to a single line |
| Summaries | hazard | "impacts by hazard" is the chart |

**No time series is ever filtered by the forcing-tier selector** (user decision 2026-08-12).
That applies to all four: the per-location hazard panels, the per-location Climate Score,
the portfolio Climate Score, and the per-hazard portfolio percentile. Only the asset-type
filter reaches them.

Climate Score views also ignore the hazard filter — the metric spans hazards by definition.

### Why it is a separate script from `generate_maps.py`

`generate_maps.py` renders **gridded** layers: ~70,849 SVG markers per panel, where marker
count is the binding constraint and the output is a multi-file collection. This renders a
**point** portfolio: hundreds of sites, where the whole delivery embeds as JSON and
interactivity is the point. Same conventions, opposite payload profile. They share vocabulary
through `scripts/utils/viz_common.py` and should not be merged.

### Decisions embedded in it

**Scenario colour follows the forcing tier, not the code.** `rcp26` and `ssp126` are the same
blue, so a mixed-round portfolio (cyclone is RCP-only, wildfire is SSP-only) shares one legend
and one filter. Native codes are what get *displayed*; the tier only drives colour and
filtering. RCP and SSP tiers are approximately comparable and any narrative must say so.
`check_tier_collisions()` refuses to let two scenarios in one layer silently share a hue.

**Values are deduped to (location, layer, scenario, decade).** Extraction is purely
coordinate-driven, so two assets at one site read identically for a given layer — verified,
not assumed: `test_customer_delivery.py` asserts no such group ever disagrees before trusting
the dedup. On the example delivery that is 24 shared groups, 0 disagreements, 342 rows → 318.

**Red always means worse.** On a `higher_is_better` layer (`conifer-npp`) an increase is an
improvement, so the diverging slope scale reverses. The map subtitle says when it has.

**Diverging limits use max |value| below 40 sites, the 95th percentile above.** The
95th-percentile rule exists because a 0.5° global grid has a heavy tail that leaves 99.7% of
cells near-white. A site portfolio has no such tail, so a percentile limit would clamp real
sites for no readability gain. Whichever rule applied, and the clamped count, is printed in
the map subtitle.

**Risk bands are display-only.** The five bands (Very Low … Very High) are derived from
`percentile` and therefore deliberately absent from `values.csv` under the no-derived-columns
rule. Thresholds carry over from the retired `export_formatter.RELATIVE_HAZARD_THRESHOLDS` so
a band means what it did in prior deliveries. Being *ordered* categories, they get the ordinal
ramp; the nominal bar charts get a single colour, never a value-ramp.

### Colour

**Red means risk; blue means magnitude.** `percentile` and Climate Score use a single-hue
**red** ramp light→dark, because 100 is worst by construction and red carries that without a
legend lookup (user decision 2026-08-12). Raw `value` keeps the blue ramp — large is not
inherently bad — which also separates "how much" from "how bad" at a glance. Signed slopes
keep the blue↔red diverging scale with a neutral gray midpoint, reversed on
`higher_is_better` layers so red is always worse.

**Forcing tiers are blue / yellow / red** for low / medium / high — the conventional
climate-scenario reading (user decision 2026-08-12). These are reference-palette slots 1, 4
and 8, **not** the validated adjacent triple (1, 2, 3), so the documented all-pairs CVD
result does not transfer.

What could be measured without the JS validator (unavailable on this machine) was measured:

| Check | Result |
|---|---|
| Red ramp luminance monotonic light→dark | yes, all 13 steps |
| Ordinal band steps vs their surface | 2.01:1 lightest on light, 2.97:1 darkest on dark — both clear the 2:1 floor |
| Blue vs red luminance separation | **1.12:1 light, 1.05:1 dark — effectively none** |
| Yellow on the light surface | 2.11:1, below the 3:1 bar (the reference palette already flags this hue) |

Blue and red therefore have no lightness fallback if hue fails, and yellow↔red converge
under deuteranopia. **Secondary encoding is applied unconditionally rather than treated as
optional**: every tier carries a distinct marker symbol (circle / square / diamond) *and*
line dash (solid / dot / dash) alongside its hue, the legend is always present, and the table
view carries every value. Identity is never colour-alone.

Light and dark modes use different steps for the red ramps — the near-black end vanishes on
`#1a1a19`, so the dark variant is truncated at a step measured at 2.42:1.

**If a JS runtime ever becomes available, re-run `scripts/validate_palette.js` on the tier
triple and the red ramps.** Everything above is the subset of the checks that pure-Python
luminance math can reach; it does not cover CVD ΔE.

## Verifying a delivery

```bash
python scripts/test_customer_delivery.py deliveries/{customer}/{date}
```

`test_shared_baseline.py` asserts a processed layer is shaped right; this asserts a delivery
faithfully carries what the layer said. The central check recomputes **every** metric from
the source NetCDF using a Gaussian weighting implemented independently of
`spatial_extract.extract_by_point`, and requires a bit-for-bit match. Calling the same
function the delivery called would prove only that it is deterministic; an independent
reimplementation is what catches a percentile inverted twice, a slope scaled by 10, a wrong
cell, or a broken longitude wrap. It also checks referential integrity, source SHA-256s,
percentile orientation and range, the `slopes_agree` rule, NaN baseline slopes, and
`data_status` consistency.

Confirmed 2026-08-12 to fail on injected corruption: a double-inverted `conifer-npp`
percentile, a ×10 `wildfire` slope, and a single wrong value were all caught (161
violations across 3,106 checks).

---

## Five things that are easy to get wrong

**`value` is not always a median.** Three of five shipped layers are boolean or extreme
zero-inflated and carry a pooled *mean*. The column is named `value` rather than `median` for
exactly this reason; `layers.csv:decadal_statistic` says which branch produced it.

**Do not re-invert the percentile.** Layers declaring `higher_is_better` applied the
inversion at processing time. `spatial_extract.apply_percentile_inversion()` exists for
pre-contract layers and would double-invert a current one — silently reversing the risk
ranking of every conifer-NPP row. The guarantee is the passthrough check in
`test_customer_delivery.py`, not an orientation spot-check: value and percentile are
independently spatially averaged, so "the highest-NPP site has a low percentile" can hold
under transformations that are still wrong.

**Slopes are per decade on every shipped layer, and that is read from the file.**
`OUTPUT-SPEC.md` fits per year and requires the layer to declare what it stored in
`slope_units`. All five shipped layers declare `decade-1`, so no conversion is applied here.
Multiplying by 10 would inflate every trend tenfold — a mistake this codebase has made.

**`slopes_agree` has three states, not two.** True when both slopes are non-zero and share a
sign; false when they disagree or exactly one has collapsed to zero; **blank when both are
zero**. That last case is an inactive site — never burns, never sees a cyclone — and
OUTPUT-SPEC requires agreement to be judged on active cells only. Writing `false` there would
make a downstream "unreliable trend" filter flag every quiet site.

**Scenarios are globbed, never listed** (GUARDRAILS §3). This is not theoretical: on
2026-08-12 the conifer layer's `rcp85` file was written mid-session and the next planning run
picked it up with no code change. `picontrol` and `historical` are excluded as not
client-facing.

## Ensemble composition can vary by scenario

`layers.csv` deliberately carries **no flat `n_members` column**. It is a per-scenario global
attribute, and `conifer-npp` reads 10/9/10 across rcp26/rcp60/rcp85 — publishing one
scenario's count as the layer's would be false for rcp60. `n_members_by_scenario` and
`members_by_scenario` are emitted instead, and `values.csv` already carries the per-cell
per-scenario count.

Attributes that *must* be identical across a layer's scenarios (`units`, `slope_units`,
`decadal_statistic`, `percentile_direction`, `variable`) are checked, and a layer that
disagrees with itself stops the delivery rather than having one scenario speak for the rest.

## QA review status travels with the delivery

Each registry layer carries `qa_reviewed_on` — the date a human read the QA report warnings
and viewed the maps. It is `null` on every layer today, which surfaces as `NOT CONFIRMED` in
`layers.csv` and as an explicit warning block in the delivered README. A layer nobody has
looked at and a reviewed one must not produce indistinguishable deliverables.

---

## Registered but unrouted — layers no asset type reaches yet

A layer in `config/layer_registry.yaml` is **not** thereby in any delivery. It reaches a
customer only when an asset type in `config/asset_catalog.yaml` names it, and that mapping is
a user decision. Layers can sit registered, contract-passing and QA-signed-off while
delivering nothing — deliberately, because declaring hazard coverage before the mapping exists
produces a report that lists a hazard as assessed and carries no data for it.

As of 2026-08-16 that applies to **eleven** layers: the nine threshold rungs
(`heatdays-hd30/hd35/hd40/hd45`, `tropicalnights-tr20/tr25`, `icedays-id`,
`frostdays-fd/fdm10` — chronic heat and cold/frost), plus `heatwave-3b` and `permafrost-3b`.
`hazard_taxonomy.yaml` keeps `covered_by: []` for all three families to match.

**When routing any of them, three things need deciding with the user, not inferred:**

1. **Which rung.** The ladder ships `hd35` and `FD` as `status: preferred` and the other seven
   as `alternate`. A data centre may want `TR20` (overnight heat drives cooling load and is
   the better mortality predictor) rather than the daytime headline; a desert asset may want
   `hd45`, which is 70.8% zeros globally and *should* be, and reads all-zero for most
   portfolios by design.
2. **Whether the absolute or the relative heat layer, or both.** `heatdays-*` and `heatwave-3b`
   answer different questions and agree on only **47.2%** of their worst-decile cells. Routing
   one is a choice about the question, not a choice of source.
3. **The direction of the cold rungs.** They are scored `higher_is_worse` on the frost COUNT,
   so they report risk *falling*. For assets harmed by losing frost — vernalisation chill
   hours, ice roads, pest overwintering kill — the risk runs the other way and these layers do
   not carry it. See `hazard_taxonomy.yaml` `families.cold-frost.materiality_note`.

## Maintaining the layer registry

The registry deliberately holds **only what a processed NetCDF cannot say about itself**:
folder, file prefix, customer-facing hazard label, status, and which slope to read.
Units, statistics, ensemble composition, percentile direction and caveats are read out of the
file at delivery time. Duplicating them in YAML would create a second source of truth that
drifts the moment a layer is reprocessed.

### `recommended_slope` is measured, not assumed

`values.csv` carries both slopes, so nothing is hidden — but a reader needs to know which one
carries signal on a given layer. Re-measure after any reprocessing:

```bash
python scripts/generate_customer_delivery.py --measure-slopes
```

It prints the sen-zero share and sign agreement on **active cells only** (either slope
non-zero) at the final decade, and flags any layer where the data disagrees with the
registry. Judging on all cells dilutes both figures with ocean and permanently-zero land and
produces the opposite conclusion.

**The measured table lives in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)**, not here — one
copy, so a reprocessed layer cannot leave two tables disagreeing. Re-run the command above
after reprocessing anything and update it there.

The point worth carrying in your head: `wildfire`'s 29.2% grid-level zero fraction suggests
Sen is safe, and it is not — the collapse lives in the year-pair *differences*, not the
values. That is why this is measured per layer rather than inferred from `field_nature`.

### Preferring the more recent or more robust option

Where a hazard has more than one shipped layer, the standing preference is the more recent or
more robust source, expressed as `status: preferred` vs `alternate`. Drought has two
siblings and `drought-3b` (ISIMIP3b/SSP, CMIP6) is preferred; `drought-2b` (`led`, 31 members
across 8 impact models, RCP) stays selectable and is the right pick when RCP scenarios or
ensemble depth is the actual requirement. `alternate` means *not chosen by default*, never
*deprecated*.

### Blocked layers

A layer folder containing `INVALID-DO-NOT-USE.md` is refused by the loader regardless of what
the registry says, and the marker's first line is quoted in the error. The `blocked:` section
is a second, independent guard. Both sugarcane layers are blocked (withdrawn 2026-08-11,
upstream data defect: ISIMIP2b LPJmL simulates no sugarcane in the sugarcane belt).

---

## Spatial averaging — the complete picture

A delivered number has been spatially averaged **twice** on some layers, and the two stages
are set by different people for different reasons. Everything measured below is dated
2026-08-12 and reproducible from the processed files.

### Stage 1 — processing-time smoothing (per layer, in the processor)

Applied to every member's annual map *before* any pooling. It is a per-layer decision
recorded in the file's `spatial_smoothing` attribute, and the decay length is a **measurement,
not a constant**.

**The per-layer table lives in [DATASET-ATTRIBUTES.md](DATASET-ATTRIBUTES.md)**, kept in one
place so it cannot drift from what the files declare. Only `cyclone` is smoothed today, and
each other layer declines it for a stated, measured reason.

**Read the layer's own `spatial_smoothing` attribute before interpreting its values.** Never
assume a decay length carries across layers — it is a measurement, not a constant, and so is
whether to smooth at all.

### Stage 2 — extraction-time point averaging (here, uniform across layers)

`spatial_extract.extract_by_point`: take every grid-cell **centre** within 0.5° of the site
(a square window, applied per axis), weight by `exp(−½(d/σ)²)` with **σ = 0.25°**, normalize,
sum. NaN cells are dropped and the surviving weights renormalized.

**In practice this is always a 4-cell blend.** The 3×3 stencil requires the site to sit
*exactly* on a cell centre; drift 0.0005° and it collapses to 2×2. Across 20,000 random
sites, 100% got 2×2. So a delivered value is a weighted blend of the four cell centres
bracketing the site, over a **1° × 1° footprint ≈ 111 km north–south**.

How lopsided that blend is depends on sub-cell position:

| Site position | Weights | Behaviour |
|---|---|---|
| On a grid vertex | `0.25 / 0.25 / 0.25 / 0.25` | flat 4-cell mean |
| Median random site | dominant ≈ `0.50` | — |
| Near a cell centre | `0.78 / 0.10 / 0.10 / 0.01` | effectively nearest-neighbour |

Two sites 20 km apart can therefore receive structurally different kernels. That is inherent
to a fixed-σ kernel on a fixed grid, not a defect, but it is why "the value at your site" is
a loose phrase.

**Four properties of stage 2, all deliberate:**

1. **Truncated, not tapered.** The 0.5° radius is exactly 2σ, where the Gaussian is still at
   13.5% of peak — the kernel is chopped off. A complete Gaussian would need radius ≥ 3σ =
   0.75°.
2. **Degree space, not physical distance.** On a regular 0.5° grid that is grid-cell space,
   which is defensible, but the ground footprint stretches with latitude: 0.5° of longitude
   is 55.7 km at the equator and 21.7 km at 67° N, so the kernel reaches **2.6× further
   east–west than north–south at 67° N**. Note stage 1's `let` kernel *does* cos(lat)-scale
   longitude; stage 2 does not. The two stages are inconsistent on this point, deliberately
   left so rather than silently changing every extraction the pipeline has ever produced.
3. **Longitude wraps the antimeridian.** It did not until 2026-08-12: the window was
   one-sided at the seam, so 180° E and 180° W returned different numbers — 0.775 against
   0.962 on `burntarea` ssp585 at 17° S, and 62× apart at 67° N. This was a defect in the
   shared `extract_by_point`, so `extract_timber_locations.py` carried it too.
4. **NaN cells are dropped and weights renormalized**, so a site whose own cell is masked
   still returns a value built from its surviving neighbours. Marching west from the Oregon
   coast at 44° N the land weight goes 100% → 100% → 69% → 0%, and `data_status` reads `OK`
   throughout the first three. **This is accepted, not a gap**: customer locations are masked
   to land cells upstream before they reach this workflow, so a genuinely offshore site
   should never arrive here (user decision 2026-08-12). If that upstream masking is ever
   removed, revisit it — the mechanism will silently source a coastal value from a minority
   of land cells.

### The two stages compound

On `cyclone` a delivered value has passed through the L=2.5 5×5 kernel **and** the extraction
blend, so its effective footprint is materially wider than the 1° × 1° box above. On the
other four layers stage 2 is the only averaging applied. Do not describe `cyclone` values as
site-specific in a customer narrative.

### Coordinate precision matters

Because the footprint is ~1°, a wrong coordinate is not a rounding error. Measured
2026-08-12, moving a site by ±0.25° — one grid cell — changed 2090s burnt area at Shasta from
1.248 to 3.979, a swing of **166% of the centre value**; Mobile swung 115%.

**Reproduce it rather than quoting it from here**, and re-measure after any reprocessing:

```bash
python scripts/measure_extraction_sensitivity.py
```

That script also reproduces the 4-cell-blend result below. Both figures used to live only in
a session transcript, which for a number that reaches a customer document is indistinguishable
from having been invented. Validate customer coordinates before a run.

### Not wired in

Polygon and region extraction exist in `scripts/utils/spatial_extract.py` but are **not**
used by this workflow. `_calculate_cell_weights` normalizes across cells in planar degree
space with no `cos(lat)` term, so any polygon spanning a wide latitude band over-weights its
high-latitude cells. Fix that before offering polygon delivery.

Ensemble counts (`n_members`, `n_models`) are weighted like every other field, so a site
straddling a mask edge yields a fractional depth; they are rounded to integers for delivery.
