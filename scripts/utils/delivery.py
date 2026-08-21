"""Customer-delivery extraction: processed TCFD layers -> a normalized CSV star schema.

This is the DETERMINISTIC half of the customer workflow. It takes a list of
location-asset rows, resolves each asset to its hazard layers via the asset catalog, pulls
every contract metric out of each layer at each point, and writes a normalized set of CSVs.
No scoring, no ranking, no narrative -- those belong to the report workflow downstream.

OUTPUT (star schema, chosen by the user 2026-08-12 -- "no redundant columns")

    locations.csv   location_id, name, lat, lon, country, state, city, region, subregion
    assets.csv      asset_id, location_id, asset_type, sub_asset_unit
    layers.csv      layer_id + everything read back OUT of the processed NetCDF
    values.csv      asset_id, layer_id, scenario, decade + the ten metrics
    manifest.json   provenance: source files, mtimes, sizes, extraction parameters
    README.md       how to read the delivery

Location metadata is written once in locations.csv, layer metadata once in layers.csv, and
values.csv is pure keys plus measurements.

THREE THINGS THAT ARE EASY TO GET WRONG HERE
--------------------------------------------

1. DO NOT RE-INVERT THE PERCENTILE. Every layer following OUTPUT-SPEC.md has already
   applied the `higher_is_better` inversion when it wrote the file -- its
   `percentile_baseline` attribute says so explicitly. Calling
   `spatial_extract.apply_percentile_inversion()` on a current-contract layer would invert a
   second time and silently reverse the risk ranking of every conifer-NPP row. That function
   exists for pre-contract layers only. `percentile_direction` travels into layers.csv as
   documentation of what was ALREADY done, never as an instruction to do it here.

2. THE SLOPES ARE ALREADY PER DECADE ON EVERY SHIPPED LAYER. OUTPUT-SPEC.md fits per YEAR
   and requires the layer to declare which it stored in `slope_units`; all five shipped
   layers declare `decade-1`. So `slope_units` is READ from the file and carried into
   layers.csv rather than assumed. Multiplying by 10 here would inflate every trend tenfold,
   which is a mistake this codebase has made before.

3. SCENARIOS ARE GLOBBED, NEVER LISTED (GUARDRAILS.md §3). A hardcoded scenario list once
   made 25% of processed data invisible. `picontrol` and `historical` are dropped -- they
   strengthen a baseline but are not client-facing.

The baseline decade legitimately carries NaN slopes (the expanding window has no span yet),
so a NaN slope in the 2020s row is the contract working, not missing data.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from .spatial_extract import (
    as_period_dataset,
    extract_by_point,
    grid_cell_size,
    normalize_longitude,
    search_radius_for,
    sigma_for,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
LAYER_REGISTRY_PATH = PROJECT_ROOT / "config" / "layer_registry.yaml"
ASSET_CATALOG_PATH = PROJECT_ROOT / "config" / "asset_catalog.yaml"
DELIVERIES_ROOT = PROJECT_ROOT / "deliveries"

#: A folder carrying this file is refused by the loader regardless of the registry.
INVALID_MARKER = "INVALID-DO-NOT-USE.md"

#: Not client-facing (see the isimip-extract-aggregate skill).
NON_DELIVERABLE_SCENARIOS = frozenset({"picontrol", "historical", "obsclim", "counterclim"})

#: Metrics pulled from every layer, in values.csv column order. `median` is renamed to
#: `value` on the way out: the contract itself notes that `median` holds a MEAN on boolean
#: and extreme-zero-inflated layers, and three of five shipped layers are in that regime.
#: The layer's `decadal_statistic` in layers.csv says which statistic produced it.
CONTRACT_METRICS = (
    "median",
    "lower_ci",
    "upper_ci",
    "percentile",
    "ols_slope",
    "sen_slope",
    "n_members",
    "n_models",
)

METRIC_OUTPUT_NAMES = {"median": "value"}

#: Counts are Gaussian-weighted like everything else, so a point straddling a mask edge
#: yields a fractional depth. Rounded to an integer for delivery.
COUNT_METRICS = frozenset({"n_members", "n_models"})

VALUES_COLUMNS = (
    "asset_id",
    "layer_id",
    "scenario",
    "decade",
    "value",
    "lower_ci",
    "upper_ci",
    "percentile",
    "ols_slope",
    "sen_slope",
    "slopes_agree",
    "n_members",
    "n_models",
    "data_status",
)

#: Global attributes lifted from the processed file into layers.csv. Everything here is a
#: fact the NetCDF asserts about itself; the registry deliberately does not restate any of
#: it (see config/layer_registry.yaml).
#:
#: `n_members` is DELIBERATELY ABSENT. It is a per-scenario global attribute, and a layer
#: whose ensemble composition varies by scenario would have one scenario's count published
#: as if it were the layer's -- `conifer-npp` reads 10/9/10 across rcp26/60/85, so a flat
#: "10" is false for rcp60. `n_members_by_scenario` below is emitted instead, and
#: values.csv already carries the per-cell per-scenario count.
LAYER_ATTRS_EXPORTED = (
    "variable",
    "long_name",
    "units",
    "slope_units",
    "decadal_statistic",
    "field_nature",
    "percentile_direction",
    "percentile_zero_fraction",
    "baseline_decade",
    "window_years",
    "ensemble_uniform_across_scenarios",
    "members_by_scenario",
    "impact_models",
    "gcms",
    "source_dataset",
    "spatial_smoothing",
    "interpretation_caveat",
    #: A layer whose value is pinned at a BOUND -- the ceiling or the floor -- for a large
    #: share of cells. There the pooled sample has no variance, so BOTH slopes go to ~0 and
    #: agree, the CI collapses, and the percentile ties. Every number is correct and the
    #: conclusion inverts, which is the same class of hazard as `relative_baseline` and is
    #: promoted to MUST-DISCLOSE in `generate_delivery_caveats.layer_caveats()`.
    #: Added 2026-08-14: this allowlist is closed, so a caveat written into a processed file
    #: under any other attribute name is silently dropped at delivery. `heatwave-3b` carries
    #: `saturation_caveat` (45.9% of ssp585 2090s cells at 1.0) and `heatwave-2b` carries
    #: `sparsity_caveat` (65.7% of land at 0); before this both were invisible downstream.
    "saturation_caveat",
    "sparsity_caveat",
    #: A layer whose SPATIAL SUPPORT is coarser than the decision it will be used for. Every
    #: layer here is 0.5 deg (~55 km), so a site value is the statistic for the cell
    #: containing the site, never for the site. Tolerable where the hazard varies over
    #: comparable distances; NOT tolerable where it turns on fine terrain. Promoted to
    #: MUST-DISCLOSE in `generate_delivery_caveats.layer_caveats()` -- added 2026-08-14 with
    #: `sealevel-2b`, where metres of elevation over hundreds of metres of ground decide the
    #: answer and two assets in one cell can differ completely while sharing a number.
    "resolution_caveat",
)

#: Attributes that MUST be identical across a layer's scenarios. If they are not, the layer
#: cannot be described by a single layers.csv row and the delivery stops rather than
#: publishing one scenario's value as the layer's.
LAYER_ATTRS_MUST_MATCH = (
    "units",
    "slope_units",
    "decadal_statistic",
    "percentile_direction",
    "variable",
)

# data_status values -------------------------------------------------------------------
STATUS_OK = "OK"
#: NaN in this layer, but the site has data in at least one other registry layer -- so it is
#: on land and this layer simply does not model the site (e.g. no conifer stand present).
STATUS_OFF_LAYER_MASK = "OFF_LAYER_MASK"
#: NaN in every registry layer -- offshore, or outside the modelled domain entirely.
STATUS_OUTSIDE_DOMAIN = "OUTSIDE_DOMAIN"

# Point-extraction parameters, recorded in the manifest so a delivery is reproducible.
#
# ASSET-CATALOG.md "Spatial averaging -- the complete picture" is authoritative for what
# these mean and what they were measured to do; do not restate it here. The three facts a
# reader of THIS file needs:
#
#   - sigma 0.25 with radius 0.5 is a 2-sigma TRUNCATED Gaussian, and on this grid it
#     resolves to a 4-cell blend for every site that is not exactly on a cell centre
#     (measured: 100% of 20,000 random sites) -- a 1 deg x 1 deg footprint, ~111 km N-S.
#   - a delivered `cyclone` value has ALSO been through that layer's processing-time L=2.5
#     5x5 kernel. The two stages compound; cyclone values are not site-specific.
#   - NaN neighbours are dropped and weights renormalized, so a masked cell still returns a
#     value. Accepted because customer locations are masked to land upstream.
EXTRACT_SIGMA = 0.25
EXTRACT_SEARCH_RADIUS = 0.5


class DeliveryError(RuntimeError):
    """A delivery cannot proceed. Always carries what the operator must fix."""


# ---------------------------------------------------------------------------------------
# Registry and catalog
# ---------------------------------------------------------------------------------------

@dataclass(frozen=True)
class LayerSpec:
    """One shipped layer, as the registry knows it.

    Deliberately thin: this holds only what a processed NetCDF cannot say about itself.
    Units, statistics, ensemble composition and caveats are read from the file at delivery
    time so the registry can never drift from the data.
    """

    layer_id: str
    folder: str
    file_prefix: str
    hazard: str
    hazard_measure: str
    status: str
    recommended_slope: str
    recommended_slope_rationale: str = ""
    delivery_note: str = ""
    #: Why a layer sits at ``status: blocked`` -- built and contract-passing, but not
    #: cleared for delivery. The registry header has always documented `blocked` as a legal
    #: status while the loader had no field to carry the reason, so the first layer to use
    #: it (``heatwave-3b``, 2026-08-14) raised a TypeError on load. A blocked layer is left
    #: in the registry rather than deleted so it is discoverable instead of invisible; it
    #: reaches no delivery because no asset type references it.
    blocked_reason: str = ""
    #: True when the layer's value is defined against a FIXED HISTORICAL REFERENCE, so a
    #: high score means "unusual for this place" and NOT "bad in absolute terms". Three
    #: shipped layers are in this class and it is the single most misreadable property any
    #: of them has: `cropfailure-3b` ranks Iowa at the 99.3rd percentile of cropland and the
    #: Sahel at 69.4, and `drought-3b`/`drought-2b` score a permanently arid cell LOW.
    #: Declared here rather than left to prose because `generate_delivery_caveats.py`
    #: promotes it to a MUST-DISCLOSE caveat, which both reports are then required to carry.
    #: Before this existed, the drought layers' own delivery_note said "Must be stated in
    #: any customer narrative" while the machinery filed it as optional.
    relative_baseline: bool = False
    #: Customer-facing wording for that caveat. Required when relative_baseline is true.
    relative_baseline_note: str = ""
    #: Date a human read this layer's QA report warnings and viewed its maps. Null until
    #: that actually happened. A contract PASS means the file is SHAPED right, not that the
    #: input is about what its name says -- both sugarcane layers passed every check and
    #: were meaningless. Absence is surfaced as "NOT CONFIRMED" in the delivered
    #: layers.csv and README rather than silently omitted.
    qa_reviewed_on: Optional[str] = None


@dataclass
class Registry:
    layers: Dict[str, LayerSpec]
    blocked: Dict[str, str]
    processed_root: Path
    measured_on: str

    def get(self, layer_id: str) -> LayerSpec:
        if layer_id not in self.layers:
            known = ", ".join(sorted(self.layers))
            raise DeliveryError(
                f"Unknown layer_id {layer_id!r}. Registered layers: {known}"
            )
        return self.layers[layer_id]


def load_registry(path: Path = LAYER_REGISTRY_PATH) -> Registry:
    raw = yaml.safe_load(path.read_text())
    layers = {
        lid: LayerSpec(layer_id=lid, **spec) for lid, spec in raw.get("layers", {}).items()
    }
    # A relative-baseline layer with no wording would emit an empty MUST-DISCLOSE caveat,
    # which is worse than none: the report would carry a heading promising a caveat and say
    # nothing under it. Fail at load rather than at render.
    missing = sorted(
        lid for lid, s in layers.items()
        if s.relative_baseline and not s.relative_baseline_note.strip()
    )
    if missing:
        raise DeliveryError(
            f"{path}: layer(s) {missing} declare relative_baseline: true but carry no "
            "relative_baseline_note. That note is rendered into every report as a "
            "must-disclose caveat and cannot be blank."
        )
    processed_root = Path(raw.get("processed_root", "data/processed"))
    if not processed_root.is_absolute():
        processed_root = PROJECT_ROOT / processed_root
    return Registry(
        layers=layers,
        blocked=raw.get("blocked") or {},
        processed_root=processed_root,
        measured_on=str(raw.get("measured_on", "")),
    )


@dataclass
class AssetCatalog:
    """Asset types, hazard families and Climate Score weights. Absence is an error.

    catalog_version 2 (model change 2026-08-20): entries carry `family_weights`
    (family -> 0/1) instead of a `layers` extraction filter; `families` holds the
    score-family definitions (family -> {"layers": [...], "observed": bool}); and
    `standard_excluded` maps layer_ids excluded from the standard set to their reasons.
    """

    entries: Dict[str, dict]
    families: Dict[str, dict] = field(default_factory=dict)
    standard_excluded: Dict[str, str] = field(default_factory=dict)
    _index: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name, entry in self.entries.items():
            for key in [name] + list(entry.get("aliases") or []):
                k = str(key).strip().lower()
                if k in self._index and self._index[k] != name:
                    raise DeliveryError(
                        f"Asset catalog alias {key!r} resolves to both "
                        f"{self._index[k]!r} and {name!r} -- resolution would be "
                        f"order-dependent. Remove one of them."
                    )
                self._index[k] = name

    def weights_for(self, canonical: str) -> Dict[str, int]:
        """The 0/1 family weights for a catalog entry; error if never walked through."""
        weights = self.entries[canonical].get("family_weights")
        if not isinstance(weights, dict):
            raise DeliveryError(
                f"Catalog entry {canonical!r} has no family_weights. The delivery model "
                f"changed 2026-08-20 (deliver-all + 0/1 family weights in the Climate "
                f"Score) and this entry has not been through its weights walk-through "
                f"with the user. Set family_weights -- every family in score_families, "
                f"each 0 or 1 -- and confirmed_on, then re-run."
            )
        return {str(k): int(v) for k, v in weights.items()}

    def resolve(self, asset_type: str) -> Tuple[str, dict]:
        key = str(asset_type).strip().lower()
        if key not in self._index:
            known = ", ".join(sorted(self.entries))
            raise DeliveryError(
                f"Asset type {asset_type!r} is not in the catalog.\n"
                f"  Known asset types: {known}\n"
                f"  This is deliberate -- an unknown asset must be worked out with the "
                f"user and added to config/asset_catalog.yaml, never defaulted."
            )
        canonical = self._index[key]
        return canonical, self.entries[canonical]


def load_asset_catalog(path: Path = ASSET_CATALOG_PATH) -> AssetCatalog:
    raw = yaml.safe_load(path.read_text())

    families: Dict[str, dict] = {}
    for name, spec in (raw.get("score_families") or {}).items():
        layers = list((spec or {}).get("layers") or [])
        if not layers:
            raise DeliveryError(f"score_families entry {name!r} lists no layers")
        families[name] = {"layers": layers, "observed": bool((spec or {}).get("observed"))}
    seen: Dict[str, str] = {}
    for fam, spec in families.items():
        for lid in spec["layers"]:
            if lid in seen:
                raise DeliveryError(
                    f"layer {lid!r} appears in two score families ({seen[lid]!r} and "
                    f"{fam!r}) -- it would be double-counted in the Climate Score")
            seen[lid] = fam

    excluded = {str(k): str(v) for k, v in
                ((raw.get("standard_set") or {}).get("excluded") or {}).items()}

    catalog = AssetCatalog(entries=raw.get("assets") or {}, families=families,
                           standard_excluded=excluded)

    # An entry WITH weights must name every family explicitly, each 0 or 1 -- a newly
    # added family must never default into or out of anyone's score. `null` weights mark
    # a legacy entry pending its walk-through; using one fails in weights_for().
    for name, entry in catalog.entries.items():
        w = entry.get("family_weights")
        if w is None:
            continue
        if not isinstance(w, dict):
            raise DeliveryError(f"asset {name!r}: family_weights must be a mapping or null")
        missing = set(families) - {str(k) for k in w}
        extra = {str(k) for k in w} - set(families)
        bad = {str(k): v for k, v in w.items() if v not in (0, 1)}
        if missing or extra or bad:
            raise DeliveryError(
                f"asset {name!r} family_weights invalid -- missing {sorted(missing)}, "
                f"unknown {sorted(extra)}, non-0/1 {bad}")
    return catalog


# ---------------------------------------------------------------------------------------
# Layer files on disk
# ---------------------------------------------------------------------------------------

def layer_dir(registry: Registry, spec: LayerSpec) -> Path:
    path = registry.processed_root / spec.folder
    if not path.is_dir():
        raise DeliveryError(
            f"Layer {spec.layer_id!r} is registered but not on disk: {path}\n"
            f"  data/ is local and ephemeral -- the layer may need reprocessing."
        )
    if (path / INVALID_MARKER).exists():
        marker = (path / INVALID_MARKER).read_text().strip().splitlines()
        reason = marker[0] if marker else "no reason recorded"
        raise DeliveryError(
            f"Layer {spec.layer_id!r} is marked invalid and cannot be delivered.\n"
            f"  {path / INVALID_MARKER}: {reason}"
        )
    if spec.folder in registry.blocked:
        raise DeliveryError(
            f"Layer {spec.layer_id!r} is blocked by the registry.\n"
            f"  {registry.blocked[spec.folder]}"
        )
    return path


def discover_scenarios(registry: Registry, spec: LayerSpec) -> List[str]:
    """Glob the scenarios that actually exist. GUARDRAILS §3 -- never a hardcoded list."""
    path = layer_dir(registry, spec)
    pattern = re.compile(rf"^{re.escape(spec.file_prefix)}_(.+)_processed\.nc$")
    scenarios = []
    for f in sorted(path.glob(f"{spec.file_prefix}_*_processed.nc")):
        m = pattern.match(f.name)
        if m and m.group(1) not in NON_DELIVERABLE_SCENARIOS:
            scenarios.append(m.group(1))
    if not scenarios:
        raise DeliveryError(
            f"Layer {spec.layer_id!r} has no deliverable scenario files in {path}"
        )
    return scenarios


def scenario_path(registry: Registry, spec: LayerSpec, scenario: str) -> Path:
    return layer_dir(registry, spec) / f"{spec.file_prefix}_{scenario}_processed.nc"


# ---------------------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------------------

class Domain(NamedTuple):
    """The modelled domain, held as ONE MASK PER LAYER rather than one unioned grid.

    WHY NOT A SINGLE UNIONED ARRAY -- this is not a style choice, it is a correctness fix.
    The previous implementation built `union = finite_a | finite_b | ...` across every
    registry layer. That is an xarray binary op, so it ALIGNS ON COORDINATES, and layers on
    different grids share no coordinate values at all: 0.5 deg centres sit at +/-(k+0.25),
    0.25 deg centres at +/-(k+0.125). Measured 2026-08-14:

        half | quarter  ->  shape (0, 0), 0 True cells

    So registering a single 0.25 deg layer would have collapsed the domain to EMPTY and
    returned OUTSIDE_DOMAIN for every customer site in every delivery -- silently, with no
    exception raised. Keeping the masks separate and OR-ing the per-layer ANSWERS makes the
    domain grid-agnostic.

    On an all-0.5-deg registry this is arithmetically identical to the union: every mask
    shares one grid, so "any finite cell within the window of any layer" is the same set
    either way.
    """

    masks: List[Tuple[str, xr.DataArray]]
    consulted: List[str]


def _domain_mask(registry: Registry) -> Domain:
    """Per-layer finite-cell masks across EVERY available registry layer.

    Used only to tell OFF_LAYER_MASK ("on land, this layer does not model your site") apart
    from OUTSIDE_DOMAIN ("offshore or off-grid"). Built from the layers rather than a
    downloaded land mask, so a delivery needs no network and stays reproducible.

    THE DOMAIN IS THE WHOLE REGISTRY, NOT THE DELIVERY'S LAYERS. Scoping it to the delivery
    makes the meaning of a status depend on what else the customer happened to order: a
    conifer-only delivery for an Amazon site would report OUTSIDE_DOMAIN -- "your site is
    offshore" -- when the truth is that no conifer stand is modelled on perfectly good land.
    Layers registered but absent from disk are skipped, so a delivery still runs on a
    partial data directory; the layers actually consulted are recorded in the manifest.
    """
    masks: List[Tuple[str, xr.DataArray]] = []
    consulted: List[str] = []
    for layer_id, spec in sorted(registry.layers.items()):
        try:
            scenario = discover_scenarios(registry, spec)[0]
        except DeliveryError:
            continue  # not on disk / withdrawn -- it just does not widen the domain
        with xr.open_dataset(scenario_path(registry, spec, scenario)) as raw:
            # Observational layers are (lat, lon); lift them onto the single-period axis
            # so `.any(dim="decade")` is valid for every layer without a special case.
            ds = as_period_dataset(raw)
            finite = np.isfinite(ds["median"]).any(dim="decade").compute()
        masks.append((layer_id, finite))
        consulted.append(layer_id)
    if not masks:
        raise DeliveryError(
            "No registry layer is available on disk, so no modelled domain can be built."
        )
    return Domain(masks=masks, consulted=consulted)


def _extraction_geometry_sentence(manifest: dict) -> str:
    """Customer-facing description of the extraction geometry.

    A single "0.5 deg grid" sentence was true of every delivery until 2026-08-14 and is
    false the moment a delivery carries layers on different grids. When the grids agree it
    reads exactly as before; when they do not it names each layer rather than picking one.
    """
    ex = manifest["extraction"]
    if ex.get("sigma_degrees") is not None:
        return (f"sigma = {ex['sigma_degrees']} deg, search radius "
                f"{ex['search_radius_degrees']} deg on a {ex['grid_degrees']} deg grid")
    parts = [
        f"{layer_id}: sigma {g['sigma_degrees']} deg, radius {g['search_radius_degrees']} "
        f"deg on a {g['cell_size_degrees']} deg grid"
        for layer_id, g in sorted(ex.get("per_layer_geometry", {}).items())
    ]
    return ("geometry follows each layer's own grid -- " + "; ".join(parts))


def _layer_geometry(domain: Domain) -> Dict[str, Dict[str, float]]:
    """Cell size and derived extraction geometry per registry layer, read from its grid."""
    out: Dict[str, Dict[str, float]] = {}
    for layer_id, mask in domain.masks:
        cell = grid_cell_size(mask)
        out[layer_id] = {
            "cell_size_degrees": round(cell, 6),
            "search_radius_degrees": round(search_radius_for(mask), 6),
            "sigma_degrees": round(sigma_for(mask), 6),
        }
    return out


def _uniform_or_none(geometry: Dict[str, Dict[str, float]], key: str) -> Optional[float]:
    """The shared value of `key` across layers, or None when they disagree.

    None is the honest answer for a mixed-resolution delivery: a single "search radius =
    0.5 deg" line in a manifest or a customer README would be false for half the layers.
    """
    values = {round(g[key], 6) for g in geometry.values()}
    return values.pop() if len(values) == 1 else None


def _point_in_domain(domain: Domain, lat: float, lon: float) -> bool:
    """True if ANY registry layer models a cell within its own search window of the site.

    The window is a multiple of each layer's OWN cell size, so a finer layer searches a
    proportionally smaller degree window and does not silently widen the domain. At 0.5 deg
    that multiple reproduces the historical constant exactly (EXTRACT_SEARCH_RADIUS = 0.5 =
    one 0.5 deg cell), so the answer for every currently shipped layer is unchanged.
    """
    lon = normalize_longitude(lon)
    for _layer_id, mask in domain.masks:
        radius = search_radius_for(mask)
        sub = mask.sel(
            lat=slice(lat - radius, lat + radius),
            lon=slice(lon - radius, lon + radius),
        )
        if sub.size == 0:
            # lat is stored descending on some layers; fall back to boolean selection.
            sub = mask.where(
                (np.abs(mask.lat - lat) <= radius) & (np.abs(mask.lon - lon) <= radius),
                drop=True,
            )
        if bool(np.asarray(sub).any()):
            return True
    return False


def _slopes_agree(ols: float, sen: float) -> Optional[bool]:
    """The dual-slope robustness signal that replaced the retired p-value.

    OUTPUT-SPEC.md requires agreement to be judged on ACTIVE cells only -- either slope
    non-zero -- so this has three outcomes, not two:

        both slopes 0      -> None. INACTIVE, not applicable. A site that never burns and
                              never sees a cyclone has a genuinely zero trend under both
                              estimators; calling that "disagreement" would make a
                              downstream "unreliable trend" filter flag every quiet site.
        exactly one is 0   -> False. This is the zero-collapse regime and the disagreement
                              is real and informative.
        both non-zero      -> do the signs match?

    None also covers a NaN slope, which is what the baseline decade legitimately carries.
    """
    if not (np.isfinite(ols) and np.isfinite(sen)):
        return None
    if ols == 0 and sen == 0:
        return None
    if ols == 0 or sen == 0:
        return False
    return (ols > 0) == (sen > 0)


def extract_layer_for_points(
    registry: Registry,
    layer_id: str,
    points: Sequence[Tuple[str, float, float]],
    domain: xr.DataArray,
) -> Tuple[List[dict], dict]:
    """Extract every scenario x decade of one layer at each (asset_id, lat, lon).

    Returns (value rows, layer metadata read back out of the NetCDF).
    """
    spec = registry.get(layer_id)
    scenarios = discover_scenarios(registry, spec)

    rows: List[dict] = []
    layer_meta: dict = {}
    per_scenario_attrs: Dict[str, dict] = {}

    for scenario in scenarios:
        path = scenario_path(registry, spec, scenario)
        with xr.open_dataset(path) as raw_ds:
            # An observational layer has no decade axis; lift it onto a single-period one
            # so the (decade x scenario) extraction below is unchanged for every layer.
            ds = as_period_dataset(raw_ds)
            per_scenario_attrs[scenario] = dict(ds.attrs)
            if not layer_meta:
                layer_meta = _layer_metadata(spec, ds, scenarios)

            for asset_id, lat, lon in points:
                data = extract_by_point(
                    ds,
                    lat=lat,
                    lon=normalize_longitude(lon),
                    variables=list(CONTRACT_METRICS),
                    # Per-LAYER geometry, from that layer's own grid. Resolves to the
                    # historical 0.5 / 0.25 on every 0.5 deg layer; a finer layer gets a
                    # proportionally smaller window instead of silently blending cells
                    # ~55 km away while claiming ~28 km resolution.
                    search_radius=search_radius_for(ds),
                    sigma=sigma_for(ds),
                )
                for decade in [int(d) for d in ds.decade.values]:
                    row = {
                        "asset_id": asset_id,
                        "layer_id": layer_id,
                        "scenario": scenario,
                        "decade": decade,
                    }
                    for metric in CONTRACT_METRICS:
                        val = data.get(metric, {}).get(decade, np.nan)
                        if metric in COUNT_METRICS and np.isfinite(val):
                            val = int(round(val))
                        row[METRIC_OUTPUT_NAMES.get(metric, metric)] = val

                    row["slopes_agree"] = _slopes_agree(
                        row["ols_slope"], row["sen_slope"]
                    )
                    if np.isfinite(row["value"]):
                        row["data_status"] = STATUS_OK
                    elif _point_in_domain(domain, lat, lon):
                        row["data_status"] = STATUS_OFF_LAYER_MASK
                    else:
                        row["data_status"] = STATUS_OUTSIDE_DOMAIN
                    rows.append(row)

    _assert_scenarios_describable(spec, per_scenario_attrs)
    layer_meta["n_members_by_scenario"] = ";".join(
        f"{s}:{per_scenario_attrs[s].get('n_members', '?')}" for s in scenarios
    )
    return rows, layer_meta


def _assert_scenarios_describable(spec: LayerSpec, attrs_by_scenario: Dict[str, dict]) -> None:
    """Refuse to publish one scenario's metadata as the whole layer's.

    layers.csv carries one row per layer, which is only honest while the load-bearing
    attributes agree across scenarios. Units or a decadal statistic that differ between
    scenarios mean the files are not the same measurement and cannot share a row.
    """
    for attr in LAYER_ATTRS_MUST_MATCH:
        seen = {s: a.get(attr, "") for s, a in attrs_by_scenario.items()}
        if len(set(seen.values())) > 1:
            detail = "; ".join(f"{s}={v!r}" for s, v in sorted(seen.items()))
            raise DeliveryError(
                f"Layer {spec.layer_id!r} disagrees with itself on {attr!r} across "
                f"scenarios and cannot be described by one layers.csv row: {detail}\n"
                f"  The scenario files are not the same measurement. Reprocess the layer."
            )


def _layer_metadata(spec: LayerSpec, ds: xr.Dataset, scenarios: Sequence[str]) -> dict:
    """Build the layers.csv row -- registry labels plus facts read from the NetCDF.

    Attributes come from the layer's first scenario. That is safe only for attributes
    verified identical across scenarios by `_assert_scenarios_describable`, plus attributes
    that are scenario-resolved by construction (`members_by_scenario`). Anything genuinely
    per-scenario is emitted as its own `*_by_scenario` column instead -- see the note on
    `n_members` at LAYER_ATTRS_EXPORTED.
    """
    meta = {
        "layer_id": spec.layer_id,
        "hazard": spec.hazard,
        "hazard_measure": spec.hazard_measure,
        "registry_status": spec.status,
        "recommended_slope": spec.recommended_slope,
        "recommended_slope_rationale": " ".join(spec.recommended_slope_rationale.split()),
        "delivery_note": " ".join(spec.delivery_note.split()),
        "relative_baseline": "yes" if spec.relative_baseline else "no",
        "relative_baseline_note": " ".join(spec.relative_baseline_note.split()),
        "qa_reviewed_on": spec.qa_reviewed_on or "NOT CONFIRMED",
        "scenarios": ";".join(scenarios),
        "decades": ";".join(str(int(d)) for d in ds.decade.values),
        "source_folder": spec.folder,
    }
    for attr in LAYER_ATTRS_EXPORTED:
        value = ds.attrs.get(attr, "")
        meta[attr] = " ".join(str(value).split()) if value != "" else ""
    return meta


# ---------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------

REQUIRED_INPUT_COLUMNS = ("Location", "Lat", "Lon", "Asset_Type")

#: Optional per-asset carrying amount. IFRS S2 paragraph 29(c) and ESRS E1-9 both require
#: the *monetary amount* -- not just the count -- of assets vulnerable to physical risk, and
#: nothing in a NetCDF can supply it. So it is a customer input, and its ABSENCE is a
#: disclosed limitation rather than a silent one: with values the compliance report meets
#: 29(c) directly, without them it reports counts and percentages and states plainly that
#: the monetary figure is customer-owned.
#:
#: `Value_Basis` matters as much as the number. A book value, an insured value and a market
#: valuation give three different answers to "amount of assets vulnerable", and a report
#: that does not say which was used is not auditable.
ASSET_VALUE_COLUMNS = (
    "Asset_Value",
    "Currency",
    "Valuation_Date",
    "Value_Basis",
)

OPTIONAL_INPUT_COLUMNS = (
    "Sub_Asset_Unit",
    "Country",
    "State",
    "City",
    "Region",
    "Subregion",
    "Layers",
    "Coord_Source",
) + ASSET_VALUE_COLUMNS


def load_input(path: Path) -> pd.DataFrame:
    """Read the customer's location-asset list.

    One row per location-asset combination. `Layers` (semicolon-separated layer_ids) is an
    optional per-row override of the asset catalog.
    """
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.strip() for c in df.columns]

    missing = [c for c in REQUIRED_INPUT_COLUMNS if c not in df.columns]
    if missing:
        raise DeliveryError(
            f"Input {path} is missing required column(s): {', '.join(missing)}\n"
            f"  Required: {', '.join(REQUIRED_INPUT_COLUMNS)}\n"
            f"  Optional: {', '.join(OPTIONAL_INPUT_COLUMNS)}"
        )
    for col in OPTIONAL_INPUT_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan

    # Coerce BEFORE any arithmetic. A customer CSV is untrusted text: a stray "N/A", a
    # degrees-minutes string or an Excel artifact in Lat makes the column dtype object, and
    # `.abs()` on it raises a bare TypeError that escapes the CLI's DeliveryError handler
    # and shows the customer a traceback.
    for col in ("Lat", "Lon"):
        coerced = pd.to_numeric(df[col], errors="coerce")
        unparseable = df[coerced.isna() & df[col].notna()]
        if len(unparseable):
            examples = ", ".join(
                f"{r['Location']!s}={r[col]!r}" for _, r in unparseable.head(5).iterrows()
            )
            raise DeliveryError(
                f"{len(unparseable)} row(s) have a non-numeric {col}: {examples}\n"
                f"  Coordinates must be decimal degrees."
            )
        df[col] = coerced

    bad = df[df["Lat"].isna() | df["Lon"].isna()]
    if len(bad):
        names = ", ".join(str(n) for n in bad["Location"].tolist()[:5])
        raise DeliveryError(
            f"{len(bad)} row(s) have no coordinates and cannot be extracted: {names}"
        )
    out_of_range = df[(df["Lat"].abs() > 90) | (df["Lon"].abs() > 360)]
    if len(out_of_range):
        names = ", ".join(str(n) for n in out_of_range["Location"].tolist()[:5])
        raise DeliveryError(f"Row(s) with out-of-range coordinates: {names}")

    # Same untrusted-text treatment as the coordinates. "$4.2M", "4,200,000" and "TBD" are
    # all things customers put in a value column, and a monetary figure that silently
    # becomes NaN would drop an asset out of the 29(c) denominator without anyone noticing.
    value_raw = df["Asset_Value"]
    stripped = value_raw.astype(str).str.replace(r"[,\s]", "", regex=True)
    coerced = pd.to_numeric(stripped.where(value_raw.notna()), errors="coerce")
    unparseable = df[coerced.isna() & value_raw.notna()]
    if len(unparseable):
        examples = ", ".join(
            f"{r['Location']!s}={r['Asset_Value']!r}" for _, r in unparseable.head(5).iterrows()
        )
        raise DeliveryError(
            f"{len(unparseable)} row(s) have a non-numeric Asset_Value: {examples}\n"
            f"  Asset_Value must be a bare number -- put the unit in Currency and the\n"
            f"  measurement basis (book / insured / market / replacement) in Value_Basis.\n"
            f"  Leave the cell EMPTY if the value is unknown; an empty cell is reported as\n"
            f"  a disclosed limitation, a wrong cell is reported as a fact."
        )
    df["Asset_Value"] = coerced

    negative = df[df["Asset_Value"] < 0]
    if len(negative):
        names = ", ".join(str(n) for n in negative["Location"].tolist()[:5])
        raise DeliveryError(f"Row(s) with a negative Asset_Value: {names}")

    # A number with no currency cannot be aggregated or disclosed. Catch it here rather
    # than letting the report print a bare total.
    valued = df["Asset_Value"].notna()
    no_currency = df[valued & df["Currency"].isna()]
    if len(no_currency):
        names = ", ".join(str(n) for n in no_currency["Location"].tolist()[:5])
        raise DeliveryError(
            f"{len(no_currency)} row(s) have an Asset_Value but no Currency: {names}"
        )
    currencies = sorted(df.loc[valued, "Currency"].astype(str).str.strip().str.upper().unique())
    if len(currencies) > 1:
        raise DeliveryError(
            f"Asset values are given in {len(currencies)} currencies: {', '.join(currencies)}\n"
            f"  This pipeline does not convert currency -- an FX rate is a financial\n"
            f"  assumption with its own date and source, and inventing one would put an\n"
            f"  unsourced number into a disclosure. Convert upstream to a single currency\n"
            f"  and record the rate and date in Value_Basis."
        )
    return df


#: Statuses a layer may carry and still be delivered. `blocked`, `superseded`,
#: `development` and `diagnostic` are all outside the standard set by construction.
DELIVERABLE_STATUSES = ("preferred", "alternate")


def standard_layer_ids(registry: Registry, catalog: AssetCatalog) -> List[str]:
    """The standard delivered set (delivery model of 2026-08-20).

    Every layer carrying a human QA date with a deliverable status, minus the catalog's
    explicit exclusions. COMPUTED, never enumerated, so a newly signed layer joins the
    standard set without a catalog edit -- and an unsigned or blocked layer can never
    slip into a delivery through this path.
    """
    out = sorted(
        lid for lid, spec in registry.layers.items()
        if spec.qa_reviewed_on
        and spec.status in DELIVERABLE_STATUSES
        and lid not in catalog.standard_excluded
    )
    if not out:
        raise DeliveryError("standard set is empty -- no QA-signed deliverable layers")
    return out


def build_plan(
    df: pd.DataFrame,
    catalog: AssetCatalog,
    registry: Registry,
    layer_override: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[dict]]:
    """Resolve locations, assets and each asset's layers WITHOUT touching the data.

    This is what `--plan` prints for user confirmation before any extraction runs. IDs are
    assigned in input order and are stable for a given input file.
    """
    locations: List[dict] = []
    loc_ids: Dict[Tuple[str, float, float], str] = {}
    assets: List[dict] = []

    def _clean(v, default=""):
        return default if pd.isna(v) else str(v).strip()

    for _, row in df.iterrows():
        name = str(row["Location"]).strip()
        lat, lon = float(row["Lat"]), normalize_longitude(float(row["Lon"]))
        key = (name, round(lat, 6), round(lon, 6))
        if key not in loc_ids:
            loc_ids[key] = f"LOC-{len(loc_ids) + 1:03d}"
            locations.append(
                {
                    "location_id": loc_ids[key],
                    "name": name,
                    "lat": lat,
                    "lon": lon,
                    # Provenance for a coordinate we supplied rather than the customer.
                    # Stage 1 of the workflow permits deriving a missing lat/lon, but a
                    # derived coordinate and a surveyed one must never be indistinguishable
                    # in the delivered file.
                    "coord_source": _clean(row.get("Coord_Source"), "supplied"),
                    "country": _clean(row.get("Country")),
                    "state": _clean(row.get("State")),
                    "city": _clean(row.get("City")),
                    "region": _clean(row.get("Region")),
                    "subregion": _clean(row.get("Subregion")),
                }
            )

        asset_type_raw = _clean(row["Asset_Type"])
        if layer_override:
            canonical, layer_ids = asset_type_raw, list(layer_override)
            source = "cli-override"
        elif _clean(row.get("Layers")):
            canonical = asset_type_raw
            layer_ids = [s.strip() for s in _clean(row["Layers"]).split(";") if s.strip()]
            source = "row-override"
        else:
            canonical, entry = catalog.resolve(asset_type_raw)
            # v2 model (2026-08-20): the catalog no longer filters extraction -- every
            # asset receives the standard set. Resolving still validates the asset type
            # and that its family weights exist (weights_for raises on a legacy entry
            # that has not been through its walk-through).
            catalog.weights_for(canonical)
            layer_ids = standard_layer_ids(registry, catalog)
            source = "standard-set"

        for layer_id in layer_ids:
            registry.get(layer_id)  # fail fast on a typo, before any I/O

        asset_value = row.get("Asset_Value")
        assets.append(
            {
                "asset_id": f"AST-{len(assets) + 1:03d}",
                "location_id": loc_ids[key],
                "asset_type": asset_type_raw,
                "catalog_entry": canonical,
                "sub_asset_unit": _clean(row.get("Sub_Asset_Unit")),
                "layer_ids": layer_ids,
                "layer_source": source,
                # Optional and usually absent. Emitted as columns either way so the schema
                # is fixed: a consumer should not have to discover whether this delivery
                # happens to carry values.
                "asset_value": None if pd.isna(asset_value) else float(asset_value),
                "currency": _clean(row.get("Currency")).upper(),
                "valuation_date": _clean(row.get("Valuation_Date")),
                "value_basis": _clean(row.get("Value_Basis")),
            }
        )

    assets_df = pd.DataFrame(assets)
    locations_df = pd.DataFrame(locations)

    # One work item per layer, carrying the assets that need it.
    by_layer: Dict[str, List[str]] = {}
    for a in assets:
        for layer_id in a["layer_ids"]:
            by_layer.setdefault(layer_id, []).append(a["asset_id"])
    work = [
        {"layer_id": lid, "asset_ids": ids, "n_assets": len(ids)}
        for lid, ids in sorted(by_layer.items())
    ]

    # Fail fast on a typo in the score-family definitions, before any I/O.
    for fam, spec in catalog.families.items():
        for lid in spec["layers"]:
            registry.get(lid)

    return locations_df, assets_df, work


# ---------------------------------------------------------------------------------------
# Delivery
# ---------------------------------------------------------------------------------------

#: Columns of climate_score.csv. Grain is (asset_id, scenario_tier, decade) -- a DIFFERENT
#: grain from values.csv, which is why it is its own table rather than extra columns there.
#:
#: v2 (2026-08-20): the score aggregates hazard FAMILIES, not layers. `hazards` names the
#: families that contributed, `layers` the member layers behind them, and the two n_hazards
#: columns count families (present vs weighted). The v1 column NAMES survive so downstream
#: readers keep working; ASSET-CATALOG.md documents the semantics.
CLIMATE_SCORE_COLUMNS = (
    "asset_id",
    "scenario_tier",
    "decade",
    "climate_score",
    "n_hazards",
    "n_hazards_expected",
    "hazards",
    "layers",
    "scenarios",
)


def compute_climate_score(
    values_df: pd.DataFrame, assets_df: pd.DataFrame, catalog: AssetCatalog
) -> pd.DataFrame:
    """Family-weighted mean risk percentile, per forcing tier per decade.

    THE MODEL (user determinations, 2026-08-20)
    -------------------------------------------
    A THREE-stage mean. Native scenario codes are averaged WITHIN a layer for a tier;
    layers are averaged WITHIN their hazard family; the score is the mean over the asset
    type's weight-1 families PRESENT at the site. Weights are 0/1 per family per asset
    type from config/asset_catalog.yaml -- an irrelevant hazard is still delivered, it
    just carries no weight. Families absent at a site (off-footprint, or no member
    publishing the tier) renormalize out: `n_hazards` / `n_hazards_expected` count
    families present vs weighted, and two scores with different n_hazards are not like
    for like, exactly as under v1.

    OBSERVED FAMILIES ENTER AS CONSTANTS -- superseding the 2026-08-18 exclusion (user
    decision 2026-08-20). A family flagged `observed` has no forcing pathway; its value,
    from the single observed period, is folded unchanged into EVERY tier x decade cell of
    the asset's score grid. A constant term damps the score's scenario and time contrast
    in proportion to its weight share, and the caveat machinery discloses that wherever
    the score travels. Off-domain observed members (tornado outside CONUS) drop out of
    the family mean like any absent member -- unobserved is never zero.

    WHY PERCENTILE, AND WHY A FORCING TIER, ARE UNCHANGED FROM v1
    -------------------------------------------------------------
    `percentile` is the only cross-hazard comparable axis and is already oriented so 100
    is worst on every layer. The tier keying is forced by the data: no native code spans
    both ISIMIP rounds, so a score keyed on a native code would average a subset while
    being labelled "across all hazards". RCP and SSP tiers are only APPROXIMATELY
    comparable; any narrative must say so. The native codes that contributed are recorded
    per row in `scenarios` (including `observed` where an observed family contributed).

    EDGES
    -----
    A layer in no score family (the flood counterfactuals 40yr/none) is delivered
    unscored. An asset whose catalog entry has no family_weights (a legacy entry, or a
    row override on an unknown type) produces NO score rows rather than a wrongly
    weighted one. An asset with no forcing family present anywhere has no (tier, decade)
    grid and likewise produces no rows -- data_status in values.csv tells that story.
    """
    from .viz_common import TIER_ORDER, is_forcing_scenario, tier_of

    layer_family = {lid: fam for fam, spec in catalog.families.items()
                    for lid in spec["layers"]}

    weights: Dict[str, Optional[Dict[str, int]]] = {}
    for _, a in assets_df.iterrows():
        entry = catalog.entries.get(str(a.get("catalog_entry")))
        w = entry.get("family_weights") if entry else None
        weights[a["asset_id"]] = (
            {str(k): int(v) for k, v in w.items()} if isinstance(w, dict) else None
        )

    df = values_df[["asset_id", "layer_id", "scenario", "decade", "percentile"]].copy()
    df = df[df["percentile"].notna()]
    df["family"] = df["layer_id"].map(layer_family)
    df = df[df["family"].notna()]
    keep = [
        (weights.get(aid) or {}).get(fam, 0) == 1
        for aid, fam in zip(df["asset_id"], df["family"])
    ]
    df = df[keep]
    if df.empty:
        return pd.DataFrame(columns=list(CLIMATE_SCORE_COLUMNS))

    forcing = df["scenario"].map(is_forcing_scenario)
    fc = df[forcing].copy()
    ob = df[~forcing].copy()

    # Observed constants: codes -> layer -> family, one value per (asset, family).
    obs_by_asset: Dict[str, List[Tuple[str, float]]] = {}
    obs_layers: Dict[Tuple[str, str], List[str]] = {}
    obs_codes: Dict[str, set] = {}
    if not ob.empty:
        o1 = (ob.groupby(["asset_id", "family", "layer_id"])
                .agg(pct=("percentile", "mean"),
                     codes=("scenario", lambda s: set(s))).reset_index())
        for (aid, fam), grp in o1.groupby(["asset_id", "family"]):
            obs_by_asset.setdefault(aid, []).append((fam, float(grp["pct"].mean())))
            obs_layers[(aid, fam)] = sorted(grp["layer_id"])
            obs_codes.setdefault(aid, set()).update(*grp["codes"])

    rows: List[dict] = []
    if not fc.empty:
        fc["scenario_tier"] = fc["scenario"].map(tier_of)
        # Stage 1: native codes within a layer.
        l1 = (fc.groupby(["asset_id", "family", "layer_id", "scenario_tier", "decade"])
                .agg(pct=("percentile", "mean"),
                     codes=("scenario", lambda s: set(s))).reset_index())
        # Stage 2: layers within a family.
        f1 = (l1.groupby(["asset_id", "family", "scenario_tier", "decade"])
                .agg(fam_pct=("pct", "mean"),
                     layers=("layer_id", lambda s: sorted(set(s))),
                     codes=("codes", lambda s: set().union(*s))).reset_index())
        # Stage 3: families within the asset; observed constants join every cell.
        for (aid, tier, dec), grp in f1.groupby(["asset_id", "scenario_tier", "decade"]):
            fams = list(zip(grp["family"], grp["fam_pct"]))
            layer_ids = sorted({x for lst in grp["layers"] for x in lst})
            codes = set().union(*grp["codes"])
            for fam, pct in obs_by_asset.get(aid, []):
                fams.append((fam, float(pct)))
                layer_ids = sorted(set(layer_ids) | set(obs_layers.get((aid, fam), [])))
                codes |= obs_codes.get(aid, set())
            w = weights.get(aid) or {}
            rows.append({
                "asset_id": aid,
                "scenario_tier": tier,
                "decade": int(dec),
                "climate_score": round(sum(p for _, p in fams) / len(fams), 2),
                "n_hazards": int(len(fams)),
                # WITHOUT THIS COLUMN THE CSV IS UNSAFE TO SORT: n_hazards alone cannot
                # say whether a score is complete. complete == (n_hazards == expected).
                "n_hazards_expected": int(sum(1 for v in w.values() if v == 1)),
                "hazards": ";".join(sorted(f for f, _ in fams)),
                "layers": ";".join(layer_ids),
                "scenarios": ";".join(sorted(codes)),
            })

    out = pd.DataFrame(rows, columns=list(CLIMATE_SCORE_COLUMNS))
    if out.empty:
        return out
    out["_t"] = out["scenario_tier"].map(lambda t: TIER_ORDER.index(t)
                                         if t in TIER_ORDER else 99)
    out = out.sort_values(["asset_id", "_t", "decade"]).drop(columns="_t")
    return out


#: The customer-delivery pipeline. Stages 3 and 4 are not built yet; they are declared here
#: so a half-finished delivery is visible in its own manifest rather than looking complete.
#: See .claude/skills/customer-delivery/SKILL.md, which is the entry point for all of them.
#:
#: ORDER MATTERS AND `caveats` COMES BEFORE THE REPORTS. The caveat set is a mechanical
#: derivation from this manifest, the CSVs and config/hazard_taxonomy.yaml, and it is an
#: INPUT to both reports rather than a summary of them: each report is required to carry
#: every `must_disclose` caveat, and the verifier checks that it does. Building it after the
#: reports would mean each report derived its own list and the two drifted apart -- and the
#: one thing a caveat list must not do is differ between two documents describing one
#: delivery.
DELIVERY_STAGES = (
    ("inputs", "Location/asset list assembled and confirmed with the user"),
    ("extract", "values.csv + climate_score.csv + the star schema"),
    ("dashboard", "dashboard.html"),
    ("caveats", "caveats.json + caveats.md -- the disclosure input for both reports"),
    ("compliance_report", "report_compliance.html -- IFRS S2 spine, mapped outward"),
    ("bespoke_report", "report_bespoke.html -- composed for this reader"),
)


#: (stage, the file it produces). Used to detect artifacts left over from a previous extract.
DOWNSTREAM_ARTIFACTS = (
    ("caveats", "caveats.json"),
    ("compliance_report", "report_compliance.html"),
    ("bespoke_report", "report_bespoke.html"),
)


def record_stage(out_dir: Path, stage: str, status: str, detail: str = "") -> None:
    """Stamp a stage's status into the delivery manifest.

    Deliveries are built in stages that may run in separate sessions, so "what has actually
    been produced for this customer?" has to be answerable from the folder itself rather
    than from anyone's memory of what they ran.
    """
    known = {s for s, _ in DELIVERY_STAGES}
    if stage not in known:
        raise DeliveryError(f"Unknown delivery stage {stage!r}; known: {sorted(known)}")
    path = out_dir / "manifest.json"
    manifest = json.loads(path.read_text()) if path.exists() else {}
    stages = manifest.setdefault("stages", {})
    stages[stage] = {
        "status": status,
        "at": datetime.now(timezone.utc).isoformat(),
        **({"detail": detail} if detail else {}),
    }
    for name, desc in DELIVERY_STAGES:
        stages.setdefault(name, {"status": "not_started", "description": desc})
    # Atomic replace: a crash midway through a plain write_text leaves truncated JSON and
    # every later stage then fails to read the manifest at all.
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    tmp.replace(path)


# ---------------------------------------------------------------------------------------
# Report configuration
# ---------------------------------------------------------------------------------------

REPORT_CONFIG_FILENAME = "report_config.yaml"

#: IFRS S2 paragraph 10(d) requires the ENTITY to define short/medium/long term and to
#: explain how those definitions link to its own planning horizons. We cannot know a
#: customer's investment horizon, so these are placeholders that map our decades onto the
#: three buckets, and `source: default` on them raises a `must_disclose` caveat. Filling the
#: customer's real horizons in is a Stage 1 conversation, not a Stage 3 guess.
DEFAULT_HORIZONS = {
    "short": {"decade": 2030, "label": "Short term (2030s)"},
    "medium": {"decade": 2050, "label": "Medium term (2050s)"},
    "long": {"decade": 2090, "label": "Long term (2090s)"},
}

#: The percentile at or above which an asset-hazard pair is called "vulnerable" for
#: IFRS S2 29(c). 80 is the lower bound of the "Very High" band already used across this
#: codebase (viz_common.RISK_BANDS), so a band label and a vulnerability call cannot
#: disagree.
#:
#: There is no natural threshold here and pretending otherwise is how this report would
#: become dishonest: the count of vulnerable assets is a monotone function of a number
#: somebody chose. So the threshold is ALWAYS disclosed, and always disclosed alongside the
#: counts it would produce at the two neighbouring bands.
DEFAULT_VULNERABILITY_THRESHOLD = 80
VULNERABILITY_SENSITIVITY = (60, 90)

#: The facets a bespoke report composes from. `company` is per-engagement and lives in the
#: delivery's dossier; the rest are library profiles under docs/reporting/profiles/.
REPORT_FACETS = ("asset", "region", "persona", "vertical", "use_case", "company")


def default_report_config(customer: str, assets_df: pd.DataFrame) -> dict:
    """The starting report configuration for a delivery. Every default is marked as one."""
    if "asset_value" in assets_df.columns:
        values = pd.to_numeric(assets_df["asset_value"], errors="coerce")
        supplied = bool(values.notna().any())
        currencies = sorted(
            {c for c in assets_df.get("currency", pd.Series(dtype=str)).fillna("") if c}
        )
        bases = sorted(
            {b for b in assets_df.get("value_basis", pd.Series(dtype=str)).fillna("") if b}
        )
        n_valued = int(values.notna().sum())
    else:  # pragma: no cover - a pre-schema delivery
        supplied, currencies, bases, n_valued = False, [], [], 0

    return {
        "customer": customer,
        "config_version": 1,
        "_readme": (
            "Written with DEFAULTS by the extract stage and never overwritten afterwards, "
            "so edits survive regeneration. Every block carries `source:`; while it reads "
            "`default` the reports disclose that we chose it rather than the customer. "
            "Change the value AND the source together."
        ),
        "horizons": {
            **{k: dict(v) for k, v in DEFAULT_HORIZONS.items()},
            "source": "default",
            "_note": (
                "IFRS S2 10(d): the entity defines these and explains the link to its own "
                "planning horizons. Replace with the customer's and set source: customer."
            ),
        },
        "vulnerability": {
            "metric": "percentile",
            "threshold": DEFAULT_VULNERABILITY_THRESHOLD,
            "sensitivity": list(VULNERABILITY_SENSITIVITY),
            # Customer-facing: this string is printed verbatim in the compliance report's
            # threshold disclosure, so it names the concept rather than the module.
            "basis": "the lower bound of the 'Very High' risk band used throughout the assessment",
            "source": "default",
            "_note": (
                "percentile is ranked against the shared 2020s GLOBAL land distribution, so "
                "this threshold means 'worse than N% of global land in the 2020s' -- a "
                "level, not a change. Trend is reported separately and the two are never "
                "merged."
            ),
        },
        "asset_values": {
            "supplied": supplied,
            "n_assets_valued": n_valued,
            "n_assets_total": int(len(assets_df)),
            "currency": currencies[0] if len(currencies) == 1 else None,
            "basis": bases[0] if len(bases) == 1 else (bases or None),
        },
        "facets": {f: None for f in REPORT_FACETS},
        "frameworks": {
            "spine": "IFRS S2",
            "mapped": ["CDP 3.1.1", "ESRS E1-9", "California SB 261"],
        },
    }


def write_report_config(out_dir: Path, customer: str, assets_df: pd.DataFrame) -> Path:
    """Seed `report_config.yaml` if it is absent. NEVER overwrites an existing one.

    Re-running the extract must not silently revert a customer's stated time horizons or an
    agreed vulnerability threshold back to our placeholders -- that would change what the
    report claims without changing anything a reviewer would look at.
    """
    path = out_dir / REPORT_CONFIG_FILENAME
    if path.exists():
        return path
    path.write_text(
        yaml.safe_dump(
            default_report_config(customer, assets_df), sort_keys=False, allow_unicode=True
        )
    )
    return path


def load_report_config(out_dir: Path) -> dict:
    path = out_dir / REPORT_CONFIG_FILENAME
    if not path.exists():
        raise DeliveryError(
            f"{path} not found. It is written by the extract stage; re-run it, or copy the "
            f"defaults from delivery.default_report_config()."
        )
    cfg = yaml.safe_load(path.read_text()) or {}
    for block in ("horizons", "vulnerability"):
        if block not in cfg:
            raise DeliveryError(f"{path} is missing the `{block}` block.")
        if "source" not in cfg[block]:
            raise DeliveryError(
                f"{path}: `{block}` has no `source:`. It must say whether the value came "
                f"from the customer or from our defaults -- the reports disclose which."
            )
    return cfg


def horizon_decades(cfg: dict) -> Dict[str, int]:
    """{'short': 2030, 'medium': 2050, 'long': 2090} from a report config."""
    return {
        k: int(cfg["horizons"][k]["decade"])
        for k in ("short", "medium", "long")
        if k in cfg["horizons"]
    }


def _sha256(path: Path) -> str:
    """Content hash of a source file.

    Path, size and mtime do not identify content -- a layer reprocessed in place can keep
    all three and change every number. The hash is what lets a delivery be re-derived and
    audited later, and the same reasoning applies to the registry and catalog YAML, whose
    hand-maintained `*_version` integers nobody reliably increments.
    """
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(text: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", str(text).strip().lower()).strip("-")
    if not slug:
        raise DeliveryError(f"Customer name {text!r} does not produce a usable folder name")
    return slug


def delivery_dir(customer: str, run_date: Optional[str] = None) -> Path:
    run_date = run_date or datetime.now().strftime("%Y%m%d")
    return DELIVERIES_ROOT / slugify(customer) / run_date


def run_delivery(
    customer: str,
    input_path: Path,
    locations_df: pd.DataFrame,
    assets_df: pd.DataFrame,
    work: List[dict],
    registry: Registry,
    catalog: AssetCatalog,
    out_dir: Path,
) -> dict:
    """Extract every work item and write the star schema. Returns the manifest."""
    domain = _domain_mask(registry)

    # Index once, not once per asset -- the naive form re-indexes inside the comprehension
    # and is quadratic in site count.
    loc_by_id = locations_df.set_index("location_id")
    coords = {
        a["asset_id"]: (
            float(loc_by_id.at[a["location_id"], "lat"]),
            float(loc_by_id.at[a["location_id"], "lon"]),
        )
        for _, a in assets_df.iterrows()
    }

    all_rows: List[dict] = []
    layer_rows: List[dict] = []
    sources: List[dict] = []

    for item in work:
        layer_id = item["layer_id"]
        spec = registry.get(layer_id)
        points = [(aid, *coords[aid]) for aid in item["asset_ids"]]
        print(f"  {layer_id}: {len(points)} asset(s) x "
              f"{len(discover_scenarios(registry, spec))} scenario(s)")
        rows, meta = extract_layer_for_points(registry, layer_id, points, domain)
        all_rows.extend(rows)
        layer_rows.append(meta)
        for scenario in discover_scenarios(registry, spec):
            p = scenario_path(registry, spec, scenario)
            stat = p.stat()
            sources.append(
                {
                    "layer_id": layer_id,
                    "scenario": scenario,
                    "path": str(p.relative_to(PROJECT_ROOT)),
                    "bytes": stat.st_size,
                    "modified_utc": datetime.fromtimestamp(
                        stat.st_mtime, tz=timezone.utc
                    ).isoformat(),
                    "sha256": _sha256(p),
                }
            )

    out_dir.mkdir(parents=True, exist_ok=True)

    values_df = pd.DataFrame(all_rows)
    if not values_df.empty:
        values_df = values_df[list(VALUES_COLUMNS)].sort_values(
            ["asset_id", "layer_id", "scenario", "decade"]
        )
        # Nullable ints: an off-mask row has no ensemble depth, and a count must not be
        # delivered as "10.0" just because some other row is NaN.
        for col in COUNT_METRICS:
            values_df[col] = values_df[col].astype("Int64")

    locations_df.to_csv(out_dir / "locations.csv", index=False)
    assets_out = assets_df.copy()
    assets_out["layer_ids"] = assets_out["layer_ids"].apply(";".join)
    assets_out.to_csv(out_dir / "assets.csv", index=False)
    pd.DataFrame(layer_rows).to_csv(out_dir / "layers.csv", index=False)
    values_df.to_csv(out_dir / "values.csv", index=False)

    scores_df = compute_climate_score(values_df, assets_df, catalog)
    scores_df.to_csv(out_dir / "climate_score.csv", index=False)

    # The score's full configuration, self-contained in the delivery: the verifier and the
    # caveat generator read THIS, not the live catalog, so a delivery stays checkable after
    # the catalog moves on. `asset_weights: null` marks an asset scored on no weights (a
    # legacy entry or an override on an unknown type) -- it has no climate_score rows.
    score_config = {
        "model": "family-weighted mean, v2 (user determinations 2026-08-20)",
        "families": catalog.families,
        "standard_excluded": catalog.standard_excluded,
        "standard_layer_ids": standard_layer_ids(registry, catalog),
        "asset_weights": {
            a["asset_id"]: (
                (catalog.entries.get(str(a.get("catalog_entry"))) or {}).get("family_weights")
            )
            for _, a in assets_df.iterrows()
        },
        "observed_families_note": (
            "Families flagged observed enter the score as constants across every decade "
            "and tier (user decision 2026-08-20, superseding the 2026-08-18 exclusion). "
            "A constant term damps the score's scenario and time contrast in proportion "
            "to its weight share."
        ),
    }
    (out_dir / "score_config.json").write_text(json.dumps(score_config, indent=2) + "\n")

    status_counts = (
        values_df["data_status"].value_counts().to_dict() if not values_df.empty else {}
    )

    manifest = {
        "customer": customer,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path),
        "registry_version": yaml.safe_load(LAYER_REGISTRY_PATH.read_text()).get(
            "registry_version"
        ),
        "registry_measured_on": registry.measured_on,
        "catalog_version": yaml.safe_load(ASSET_CATALOG_PATH.read_text()).get(
            "catalog_version"
        ),
        "config_sha256": {
            "layer_registry.yaml": _sha256(LAYER_REGISTRY_PATH),
            "asset_catalog.yaml": _sha256(ASSET_CATALOG_PATH),
            "delivery.py": _sha256(Path(__file__)),
            "spatial_extract.py": _sha256(Path(__file__).with_name("spatial_extract.py")),
        },
        "qa_review": {
            layer["layer_id"]: layer["qa_reviewed_on"] for layer in layer_rows
        },
        "extraction": {
            "mode": "point",
            "weighting": "gaussian distance-weighted over cell centres",
            # Geometry is PER LAYER because layers may sit on different grids. The scalars
            # below stay populated whenever every delivered layer shares one cell size --
            # which is every delivery of only 0.5 deg layers -- and are null otherwise, so
            # nothing can quote a single radius for a mixed-resolution delivery.
            "sigma_degrees": _uniform_or_none(_layer_geometry(domain), "sigma_degrees"),
            "search_radius_degrees": _uniform_or_none(
                _layer_geometry(domain), "search_radius_degrees"),
            "grid_degrees": _uniform_or_none(_layer_geometry(domain), "cell_size_degrees"),
            "per_layer_geometry": _layer_geometry(domain),
            "nan_handling": "NaN cells excluded, remaining weights renormalized",
            "longitude_wrapping": "search window wraps the antimeridian",
            "domain_mask_layers": ";".join(domain.consulted),
            "percentile_inversion_applied_here": False,
            "percentile_note": (
                "Percentiles are delivered exactly as stored. Layers declaring "
                "higher_is_better already applied the inversion at processing time; "
                "inverting again here would reverse the risk ranking."
            ),
            "slope_units_note": (
                "slope_units is read per layer from the processed file and reported in "
                "layers.csv. No unit conversion is applied."
            ),
        },
        "counts": {
            "locations": int(len(locations_df)),
            "assets": int(len(assets_df)),
            "layers": int(len(layer_rows)),
            "value_rows": int(len(values_df)),
            "climate_score_rows": int(len(scores_df)),
            "data_status": {k: int(v) for k, v in status_counts.items()},
        },
        "climate_score": {
            "definition": (
                "Family-weighted mean of `percentile`, per forcing tier per decade "
                "(model v2, user determinations 2026-08-20): codes average within a "
                "layer, layers within their hazard family, and the score is the mean "
                "over the asset type's weight-1 families present at the site. Observed "
                "families enter as constants across every tier and decade. Percentile "
                "is the only cross-hazard comparable axis and is already oriented for "
                "risk (100 = worst on every layer). Full configuration in "
                "score_config.json."
            ),
            "keyed_on": "forcing tier, NOT a native scenario code",
            "why": (
                "No native code spans both ISIMIP rounds -- an rcp code sees only the 2b "
                "layers and an ssp code only the 3b ones -- so a score keyed on a native "
                "code would average a subset and be labelled 'across all hazards'. "
                "RCP and SSP tiers are only approximately comparable; say so in narrative."
            ),
            "hazard_weighting": (
                "0/1 per hazard family per asset type, set with the user in the "
                "2026-08-20 walk-through; families weight equally among themselves and "
                "absent families renormalize out (n_hazards discloses the count)"
            ),
            # Compared per row against ITS OWN asset's expected count. Comparing against
            # the global maximum called every warehouse row incomplete in a portfolio that
            # also held timber.
            "incomplete_rows": int(
                (scores_df["n_hazards"] < scores_df["n_hazards_expected"]).sum()
            ) if len(scores_df) else 0,
        },
        "asset_values": {
            "supplied": bool(pd.to_numeric(assets_df.get("asset_value"), errors="coerce").notna().any())
            if "asset_value" in assets_df.columns else False,
            "n_assets_valued": int(
                pd.to_numeric(assets_df.get("asset_value"), errors="coerce").notna().sum()
            ) if "asset_value" in assets_df.columns else 0,
            "note": (
                "IFRS S2 29(c) and ESRS E1-9 require the monetary AMOUNT of assets "
                "vulnerable to physical risk. No processed layer can supply it. Where "
                "values are absent the reports disclose counts and percentages only and "
                "name the monetary figure as customer-owned."
            ),
        },
        "source_files": sources,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_report_config(out_dir, customer, assets_out)

    # A re-extract rewrites the manifest, which resets every downstream stage to
    # not_started -- while the previous run's caveats and reports are still sitting in the
    # folder, still openable, still shippable. Mark them STALE rather than letting a
    # not_started manifest sit beside a finished-looking report; the verifier refuses a
    # stale artifact, so the only way forward is to rebuild them against the new extract.
    for stage, filename in DOWNSTREAM_ARTIFACTS:
        if (out_dir / filename).exists():
            record_stage(
                out_dir, stage, "stale",
                f"{filename} predates this extract ({manifest['generated_utc']}); rebuild it",
            )
    record_stage(out_dir, "inputs", "confirmed",
                 f"{len(locations_df)} locations, {len(assets_df)} assets from "
                 f"{Path(input_path).name}")
    record_stage(out_dir, "extract", "built",
                 f"{len(values_df)} value rows, {len(scores_df)} climate-score rows")
    _write_readme(out_dir, manifest, layer_rows)

    # Echo the exact input alongside the outputs so the delivery is self-contained.
    (out_dir / "input_locations.csv").write_text(Path(input_path).read_text())

    return manifest


def _write_readme(out_dir: Path, manifest: dict, layer_rows: List[dict]) -> None:
    layers_md = "\n".join(
        f"| `{r['layer_id']}` | {r['hazard']} | {r['hazard_measure']} | "
        f"{r.get('units','')} | {r.get('slope_units','')} | `{r['recommended_slope']}` | "
        f"{r['qa_reviewed_on']} |"
        for r in layer_rows
    )
    unreviewed = [r["layer_id"] for r in layer_rows if r["qa_reviewed_on"] == "NOT CONFIRMED"]
    review_md = (
        "> **QA review not confirmed for: "
        + ", ".join(f"`{lid}`" for lid in unreviewed)
        + ".** For these layers this delivery cannot state that the QA report warnings were\n"
        "> read and the maps were viewed. Passing the output contract means a file is\n"
        "> *shaped* right, not that its input is about what its name says -- two layers have\n"
        "> passed every check in this pipeline and been meaningless. Treat the numbers as\n"
        "> provisional until the layers are reviewed.\n"
        if unreviewed
        else "All layers in this delivery have a recorded QA review date.\n"
    )
    text = f"""# Climate hazard extract -- {manifest['customer']}

Generated {manifest['generated_utc']} from ISIMIP processed layers.
{manifest['counts']['locations']} location(s), {manifest['counts']['assets']} asset(s),
{manifest['counts']['value_rows']} value rows.

## Files

| File | Contents |
|---|---|
| `locations.csv` | One row per distinct site. Join key `location_id`. |
| `assets.csv` | One row per location-asset combination. Join key `asset_id`. |
| `layers.csv` | One row per hazard layer: units, ensemble, statistic, caveats. |
| `values.csv` | The measurements. `asset_id` x `layer_id` x `scenario` x `decade`. |
| `climate_score.csv` | Cross-hazard Climate Score. `asset_id` x `scenario_tier` x `decade`. |
| `manifest.json` | Provenance: source files, mtimes, extraction parameters. |
| `dashboard.html` | Interactive QA dashboard — open this first. |
| `input_locations.csv` | Verbatim copy of the submitted input. |

## Layers in this delivery

| layer_id | Hazard | Measure | Units | Slope units | Read this slope | QA reviewed |
|---|---|---|---|---|---|---|
{layers_md}

{review_md}
Ensemble composition can vary by scenario, so `layers.csv` reports
`n_members_by_scenario` and `members_by_scenario` rather than a single count.

## Reading `values.csv`

- `value` is the decadal central statistic. **Which** statistic differs by layer -- see
  `decadal_statistic` in `layers.csv` (`pooled_median`, `pooled_mean_boolean`, or
  `pooled_mean_zero_inflated`). It is a mean, not a median, on the boolean and
  zero-inflated layers.
- `lower_ci` / `upper_ci` bound the same decade pool. On median layers they are the 25th
  and 75th percentiles; on mean layers they are mean -/+ 1 SD. `ci_definition` in the
  source file records which.
- `percentile` is 1-100 against the shared 2020s baseline and is **already oriented for
  risk**: on a `higher_is_better` layer it has been inverted at processing time, so 100
  always means highest risk. Do not invert it again.
- `ols_slope` and `sen_slope` are both reported because they fail in opposite regimes.
  Read the one named in `recommended_slope`; that choice is measured per layer and the
  measurement is in `recommended_slope_rationale`.
- `slopes_agree` is the robustness signal, judged on ACTIVE cells only. True when both
  slopes are non-zero and share a sign; false when they disagree or when one has collapsed
  to zero; **blank when both are zero**, which means the site is inactive for this hazard
  (never burns, never sees a cyclone) rather than that the trends disagree. There is no
  p-value under this contract -- disagreement is what tells you a trend is not robust.
- Slopes are **NaN in the baseline decade** by design; the expanding window has no span
  there. That is the contract working, not missing data.
- `data_status`: `OK`; `OFF_LAYER_MASK` = the site is on modelled land but this layer does
  not cover it (e.g. no conifer stand present); `OUTSIDE_DOMAIN` = offshore or off-grid.

## Reading `climate_score.csv`

The **Climate Score** is a family-weighted mean of `percentile` (model v2, 2026-08-20):
native codes average within a layer, layers average within their hazard family, and the
score is the mean over the asset type's weight-1 families present at the site (weights are
0/1 per family per asset type; full configuration in `score_config.json`). Families flagged
*observed* (landslide, storm-convective) enter as constants across every decade and tier,
which damps the score's scenario and time contrast in proportion to their weight share.
Percentile is the only cross-hazard comparable axis, and it is already oriented so that 100
is worst on every layer — which is what makes the average meaningful. Higher score = higher
aggregate physical climate risk. `n_hazards` / `n_hazards_expected` count families present
vs weighted — two scores with different `n_hazards` are not like for like.

- **Keyed on a forcing tier (`low` / `medium` / `high`), not a scenario code.** No native
  ISIMIP code spans both rounds: an `rcp*` code sees only the ISIMIP2b layers and an `ssp*`
  code only the ISIMIP3b ones. A score keyed on a native code would therefore average a
  *subset* of an asset's hazards while claiming to cover all of them. The `scenarios` column
  lists exactly which native codes contributed to each row.
- **RCP and SSP tiers are only approximately comparable.** They are different scenario
  families from different CMIP generations. State that in any narrative built on this score.
- **Check `n_hazards` before comparing two scores.** A hazard that does not cover a site is
  excluded rather than counted as zero, and ISIMIP3b layers have no 2010s panel — so an
  early decade or an off-mask site legitimately scores on fewer hazards. Two scores with
  different `n_hazards` are not like for like.
- **Hazards are weighted equally.** There is no materiality weighting. What keeps an
  irrelevant hazard out of an asset's average is the asset catalog, not the arithmetic.

## Scenarios

Native ISIMIP codes are delivered as-is in `values.csv`. ISIMIP2b layers carry
`rcp26`/`rcp60`/`rcp85` and ISIMIP3b layers carry `ssp126`/`ssp370`/`ssp585`; the two
families are **not** the same scenarios and are not harmonized there. `climate_score.csv` is
the one exception, by necessity — see above.

## Extraction

Point extraction, Gaussian distance weighting over grid-cell centres
({_extraction_geometry_sentence(manifest)}). The search window
wraps the antimeridian. NaN cells are excluded and the remaining weights renormalized, so a
coastal site can draw on its land neighbours.

`manifest.json` records a SHA-256 for every source file and for the code and configuration
that produced this delivery, so it can be re-derived and audited.
"""
    (out_dir / "README.md").write_text(text)
