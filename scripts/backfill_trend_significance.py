#!/usr/bin/env python
"""Append trend significance (`trend_pvalue`, `trend_tau`, `trend_n_obs`) to a
published TCFD layer, deriving the test from the ingested RAW annual data.

Why this exists
---------------
The published ``trend`` is a baseline-anchored rate built from two decadal
numbers (GUARDRAILS S10), so it carries no p-value. Both the legacy customer
schema (``Decadal_Trend_Significance``, ``Long_Term_Trend_Significance``) and
``scripts/utils/export_formatter.py`` already expect one; until now it resolved
to NaN and was silently reported as "not significant".

What it does
------------
For each scenario it rebuilds the **ensemble-mean annual series** from the raw
member files — reusing each processor's own ``parse_name`` / ``load_member`` so
the annualization rule, unit conversion and model exclusions cannot drift — then
runs an expanding-window Mann-Kendall (see ``utils.trend_significance``). The
three new variables are APPENDED to the existing processed files; every
pre-existing variable is copied through and checked bit-identical before publish.

Members are pooled as a FLAT MEAN per year. For most layers that is exactly how
the published ``median`` is built, which gives a strong reconstruction gate: the
decade-mean of the annual series must reproduce ``median``. The two timber layers
publish a family-mean-of-family-means instead, so there the gate is reported
rather than enforced, and ``significance_pooling`` records the difference.

Usage::

    python scripts/backfill_trend_significance.py <layer_id> [--check-only]
    python scripts/backfill_trend_significance.py --all [--check-only]

``--check-only`` computes everything and prints the reconstruction diagnostics
without staging or publishing, which is the right first move on a new layer.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import sys
import warnings
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "isimip-pipeline" / "src"))

from isimip_pipeline import storage                      # noqa: E402
from utils import trend_significance as tsig             # noqa: E402
from utils.layer_publish import manifest_from_processed  # noqa: E402

#: Year grid the significance test runs on: the baseline decade's first year
#: through the last decade's last year. Every layer's decade axis is
#: 2020..2090 with window_years=10 (verified from the published attrs).
YEAR0, YEAR1 = 2020, 2099

#: Reconstruction tolerance, relative to each layer's own p99 magnitude so one
#: number fits kg m-2, % and a bare frequency. The decade-mean of the annual
#: ensemble-mean series is algebraically IDENTICAL to the published `median`
#: wherever every member covers every year of the decade, so a residual above
#: float noise means the two estimators saw different member sets.
RECON_RTOL = 1e-4

#: Share of cell-decades allowed to exceed RECON_RTOL. Not zero, because the two
#: estimators legitimately differ where a member's coverage varies WITHIN a
#: decade: mean-over-years-of-mean-over-members is only the same as
#: mean-over-members-of-mean-over-years when the member set is constant across
#: the decade's years. csoil and burntarea each have a handful of such cells
#: (a monthly member dropping an incomplete year, or a model's mask moving).
#: Gating on the max instead would reject those layers, while gating on the
#: median alone would have MISSED the fldfrc fraction-vs-percent bug, which was
#: wrong on every cell by a constant factor. This bounds the exceptions instead.
RECON_MAX_EXCEED_FRAC = 1e-3


def log(msg: str = "") -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Per-layer adapters
# ---------------------------------------------------------------------------
class Adapter:
    """Binds a layer to the processor that produced it.

    The point is to REUSE the processor's member construction rather than
    restate it: annualization (burntarea sums months, csoil means them), unit
    conversion, calendar handling and model exclusions all live in the processor
    and are imported here, so the significance pass cannot silently diverge from
    the layer it describes.
    """

    def __init__(self, layer_id: str, module: str, pattern: str,
                 to_annual: Optional[Callable] = None,
                 file_filter: Optional[Callable] = None,
                 strict_recon: bool = True,
                 value_scale: float = 1.0,
                 pooling_note: str = "flat mean across members within each year",
                 setup: Optional[Callable] = None):
        self.layer_id = layer_id
        self.module_name = module
        self.pattern = pattern
        self._to_annual = to_annual
        self._file_filter = file_filter
        self.strict_recon = strict_recon
        #: Multiplier the PROCESSOR applies after loading but that its
        #: `load_member` does not, e.g. fldfrc's fraction -> % of cell area. MK is
        #: rank-based so this cannot change p or tau, but the reconstruction gate
        #: compares against the published median and would otherwise reject a
        #: correct layer for a units mismatch (it did, which is how this was found).
        self.value_scale = value_scale
        self.pooling_note = pooling_note
        self._setup = setup
        self.mod = importlib.import_module(module)
        self.ctx: Dict = {}

    # -- discovery ---------------------------------------------------------
    def files(self) -> List[str]:
        paths = [str(p) for p in storage.stage_raw(self.layer_id, self.pattern)]
        if not paths:
            raise FileNotFoundError(
                f"{self.layer_id}: no raw files matching {self.pattern!r} under "
                f"s3://{storage.BUCKET}/{storage.raw_prefix(self.layer_id)}")
        excluded = getattr(self.mod, "EXCLUDED_MODELS", {})
        kept = []
        dropped: Dict[str, int] = {}
        for f in paths:
            model = self.mod.parse_name(f)["model"]
            if model in excluded:
                dropped[model] = dropped.get(model, 0) + 1
                continue
            if self._file_filter and not self._file_filter(self, f):
                continue
            kept.append(f)
        for model, n in sorted(dropped.items()):
            log(f"  EXCLUDING {n} {model} files -- {excluded[model][:70]}")
        if self._setup:
            self._setup(self, paths)
        return sorted(kept)

    def scenario_of(self, fpath: str) -> str:
        return self.mod.parse_name(fpath)["scenario"]

    def member_of(self, fpath: str) -> str:
        return self.mod.parse_name(fpath)["member"]

    # -- member values -----------------------------------------------------
    def annual(self, fpath: str) -> xr.DataArray:
        """One member as an annual DataArray with a `year` coord, in layer units."""
        da = self._to_annual(self, fpath) if self._to_annual \
            else self.mod.load_member(fpath)
        return da if self.value_scale == 1.0 else da * self.value_scale


def _annual_first_of_tuple(ad: Adapter, fpath: str) -> xr.DataArray:
    """For loaders that return (da, *extras) — drought, fldfrc."""
    out = ad.mod.load_member(fpath)
    return out[0] if isinstance(out, tuple) else out


def _fldfrc_filter(ad: Adapter, fpath: str) -> bool:
    """Keep only this layer's protection level.

    Each protection level has its own raw prefix, but the glob is shared, so a
    mis-ingested file would otherwise be pooled into the wrong layer.
    """
    want = ad.layer_id.split("fldfrc-")[1].split("_")[0]
    return ad.mod.parse_name(fpath)["protection"] == want


def _timber_setup(ad: Adapter, paths: List[str]) -> None:
    """Index the `pft-*` cover files that per-gridcell members are divided by."""
    cover = {}
    for f in paths:
        info = ad.mod.parse_name(f)
        if info["var"].startswith("pft-"):
            cover[(info["model"], info["gcm"], info["scenario"])] = f
    ad.ctx["cover"] = cover


def _timber_value_filter(ad: Adapter, fpath: str) -> bool:
    """Value files only — the `pft-*` cover files are inputs, not members."""
    return not ad.mod.parse_name(fpath)["var"].startswith("pft-")


def _timber_annual(ad: Adapter, fpath: str) -> xr.DataArray:
    """Harmonize one timber member to per-tile density in the layer's units.

    Mirrors process_timber_tempnle.main()'s pass 1, but per YEAR instead of per
    decade: per-gridcell models are divided by their own annual cover with the
    1% floor (GUARDRAILS S13 — the two conventions are not comparable), then the
    unit scale is applied. Ranks are invariant to the scale factor, but the cover
    division is not: cover moves year to year, so skipping it would test a
    different signal.
    """
    info = ad.mod.parse_name(fpath)
    da = ad.mod.load_raw(fpath)
    scale = ad.ctx["unit_scale"]
    if info["model"] in ad.mod.PER_GRIDCELL:
        key = (info["model"], info["gcm"], info["scenario"])
        cov_path = ad.ctx["cover"].get(key)
        if cov_path is None:
            raise FileNotFoundError(
                f"{os.path.basename(fpath)}: {info['model']} is per-gridcell but no "
                f"pft- cover file was staged for {key}; refusing to pool it uncorrected")
        cov = ad.mod.cover_fraction(ad.mod.load_raw(cov_path))
        cov = cov.reindex(year=da.year.values)
        floor = ad.mod.COVER_FLOOR
        ok = np.isfinite(cov.values) & (cov.values >= floor)
        vals = np.where(ok, da.values / np.maximum(cov.values, floor), np.nan)
        da = da.copy(data=vals)
    return da * scale


def build_adapters() -> Dict[str, Adapter]:
    """One adapter per raw-backed layer."""
    ads = {}
    for lid, module, pattern, to_annual in [
        ("drought_driedarea_annual", "process_driedarea_drought",
         "*_driedarea_global_annual_2015_2100.nc", _annual_first_of_tuple),
        ("soilcarbon_csoil_annual", "process_csoil_soilcarbon",
         "*_csoil-total_global_*_2015_2100.nc", None),
        ("wildfire_burntarea_annual", "process_burntarea_fire",
         "*_burntarea-total_global_*_2015_2100.nc", None),
    ]:
        ads[lid] = Adapter(lid, module, pattern, to_annual=to_annual)

    for protection in ("none", "100yr", "flopros"):
        lid = f"riverflood_fldfrc-{protection}_annual"
        ads[lid] = Adapter(lid, "process_fldfrc_flood",
                           "*_fldfrc_*_halfdeg_global_annual_*.nc",
                           to_annual=_annual_first_of_tuple,
                           file_filter=_fldfrc_filter,
                           # load_member returns the raw FRACTION; the processor
                           # multiplies by 100 to reach % of cell area.
                           value_scale=100.0)

    for track in ("cveg", "npp"):
        lid = f"timber_{track}-tempnle_annual"
        ad = Adapter(
            lid, "process_timber_tempnle",
            f"*_global_annual_2006_2099.nc4",
            to_annual=_timber_annual,
            file_filter=_timber_value_filter,
            setup=_timber_setup,
            strict_recon=False,
            pooling_note=(
                "flat mean across members within each year. NOTE: the published "
                "median pools family-mean-of-family-means with a >=2-family mask, "
                "so the significance series is not the identical estimator -- it "
                "is the same members, equally weighted. Coverage is still taken "
                "from median, so no cell gains a p-value that median masks out."),
        )
        ad.ctx["unit_scale"] = importlib.import_module(
            "process_timber_tempnle").TRACKS[track]["unit_scale"]
        ad.ctx["track"] = track
        ads[lid] = ad
    return ads


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------
def annual_ensemble_mean(ad: Adapter, files: List[str], lats: np.ndarray,
                         lons: np.ndarray) -> tuple:
    """Flat member mean of the annual values, on the YEAR0..YEAR1 grid.

    Accumulated member by member so peak memory stays at two grids rather than
    the whole (member, year, lat, lon) stack, which would be 1.4 GB for csoil.

    Returns:
        ``(years, mean, n_members)`` where ``mean`` is ``(n_year, lat, lon)``
        float32 and ``n_members`` counts contributing members per year and cell.
    """
    years = np.arange(YEAR0, YEAR1 + 1)
    shape = (len(years), len(lats), len(lons))
    total = np.zeros(shape, np.float64)
    count = np.zeros(shape, np.float32)

    for f in files:
        da = ad.annual(f)
        if da.year.size == 0:
            log(f"    WARNING: {os.path.basename(f)} has no usable years; skipped")
            continue
        if not (np.array_equal(da.lat.values, lats)
                and np.array_equal(da.lon.values, lons)):
            raise ValueError(
                f"{os.path.basename(f)}: grid does not match the published layer "
                f"({da.lat.size}x{da.lon.size} vs {len(lats)}x{len(lons)})")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vals = da.reindex(year=years).values.astype(np.float32)
        good = np.isfinite(vals)
        total[good] += vals[good]
        count += good
        del da, vals, good

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / np.maximum(count, 1), np.nan)
    return years, mean.astype(np.float32), count


def reconstruction_report(mean_annual: np.ndarray, years: np.ndarray,
                          published_median: np.ndarray, decades: List[int],
                          window: int, n_mem: Optional[np.ndarray] = None) -> Dict:
    """Compare the decade-mean of the annual series with the published median.

    Wherever every member covers every year of a decade this is an algebraic
    identity, so a material residual means the reconstruction is NOT the same
    ensemble the layer was built from — which would make the p-value describe a
    different quantity than the trend beside it. This is the gate that catches
    that, and it is why the member construction is imported rather than rewritten.

    Reported PER DECADE, because the baseline decade is expected to differ: most
    processors publish a SHARED 2020s baseline averaged across scenarios so that
    change and trend maps are comparable, whereas the significance test uses each
    scenario's own annual trajectory (the honest test of that scenario). Only the
    post-baseline decades are therefore gated.
    """
    scale = max(float(np.nanpercentile(np.abs(published_median), 99)), 1e-12)
    tol = RECON_RTOL * scale
    per_decade = {}
    for i, d in enumerate(decades):
        sel = (years >= d) & (years <= d + window - 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mine = np.nanmean(mean_annual[sel], axis=0)
        both = np.isfinite(mine) & np.isfinite(published_median[i])
        if not both.any():
            continue
        ad = np.abs(mine[both] - published_median[i][both])
        exceed = ad > tol
        entry = {
            "n": int(ad.size),
            "median_abs": float(np.median(ad)),
            "max_abs": float(ad.max()),
            "max_rel": float(ad.max() / scale),
            "n_exceed": int(exceed.sum()),
            "frac_exceed": float(exceed.mean()),
        }
        if n_mem is not None and exceed.any():
            # Direct test of the ragged-coverage explanation rather than an
            # assertion of it: does the member count actually move across the
            # decade's years on the cells that disagree?
            cnt = n_mem[sel][:, both]
            varying = cnt.min(axis=0) != cnt.max(axis=0)
            entry["exceed_ragged_share"] = float(np.mean(varying[exceed]))
        per_decade[d] = entry

    gated = [v for d, v in per_decade.items() if d != decades[0]]
    n_gated = sum(v["n"] for v in gated)
    return {
        "per_decade": per_decade,
        "scale_p99": scale,
        "n": sum(v["n"] for v in per_decade.values()),
        "baseline": per_decade.get(decades[0]),
        "max_rel": max((v["max_rel"] for v in gated), default=0.0),
        "max_abs": max((v["max_abs"] for v in gated), default=0.0),
        "median_abs": float(np.median([v["median_abs"] for v in gated])) if gated else 0.0,
        "n_exceed": sum(v["n_exceed"] for v in gated),
        "frac_exceed": (sum(v["n_exceed"] for v in gated) / n_gated) if n_gated else 0.0,
        "ragged_share": float(np.mean([v["exceed_ragged_share"] for v in gated
                                       if "exceed_ragged_share" in v]))
        if any("exceed_ragged_share" in v for v in gated) else float("nan"),
    }


def significance_for_scenario(ad: Adapter, files: List[str], ds: xr.Dataset,
                              decades: List[int], window: int) -> Dict:
    """Build the three new fields for one scenario, masked to `median`'s coverage."""
    lats, lons = ds.lat.values, ds.lon.values
    years, mean_annual, n_mem = annual_ensemble_mean(ad, files, lats, lons)

    covered = int(np.sum(np.isfinite(mean_annual).any(axis=0)))
    log(f"    annual ensemble mean: {len(years)} years, {covered:,} cells with data, "
        f"max members/cell={int(n_mem.max())}")

    published = ds["median"].values
    recon = reconstruction_report(mean_annual, years, published, decades, window,
                                  n_mem=n_mem)
    if recon.get("n"):
        log("    reconstruction vs published median, per decade "
            f"(scale={recon['scale_p99']:.4g}, tol={RECON_RTOL:.0e}):")
        for d, v in sorted(recon["per_decade"].items()):
            tag = "  <- baseline, shared across scenarios; not gated" \
                if d == decades[0] else ""
            rag = ""
            if "exceed_ragged_share" in v:
                rag = f", {100*v['exceed_ragged_share']:.0f}% of them ragged-coverage"
            log(f"       {d}s median|d|={v['median_abs']:.3e}  max={v['max_abs']:.3e}"
                f"  over-tol {v['n_exceed']:,}/{v['n']:,} "
                f"({100*v['frac_exceed']:.4f}%{rag}){tag}")
        log(f"       gated decades: {recon['n_exceed']:,} of "
            f"{recon['n'] - (recon['baseline'] or {}).get('n', 0):,} cell-decades "
            f"over tolerance ({100*recon['frac_exceed']:.4f}%, "
            f"limit {100*RECON_MAX_EXCEED_FRAC:.4f}%)")

    p, tau, nobs = tsig.mk_expanding(
        years, mean_annual, decades, window_years=window,
        baseline_decade=decades[0])
    # Sen slope on the DECADAL median series, NOT the annual one. On a zero-inflated
    # hazard the annual version returns exactly 0 wherever most year-pairs are
    # 0-to-0 -- 91% of driedarea ssp126 cells, 25% of ssp585 cells pairing a
    # significant p-value with a zero slope. See trend_significance's rationale block.
    sen = tsig.theilsen_decadal(
        published, decades, window_years=window, baseline_decade=decades[0])

    # Coverage is the layer's, not the reconstruction's: never emit a p-value on
    # a cell whose median is masked, and never leave one blank where it is not.
    off = ~np.isfinite(published)
    p[off] = np.nan
    tau[off] = np.nan
    sen[off] = np.nan
    nobs[off] = 0
    return {"p": p, "tau": tau, "sen": sen, "n_obs": nobs, "recon": recon,
            "n_years": len(years)}


def annotate(ds: xr.Dataset, out: Dict, ad: Adapter, decades: List[int],
             window: int) -> xr.Dataset:
    """Append the three variables and their provenance attrs to one file."""
    ds = ds.copy()
    ds = add_significance_vars(ds, out["p"], out["tau"], out["n_obs"],
                               of="the layer's median value")
    ds = replace_trend_with_sen(ds, out["sen"], decades, window)
    ds.attrs["significance_method"] = tsig.METHOD
    ds.attrs["significance_definition"] = tsig.significance_definition(
        decades, window_years=window, baseline_decade=decades[0])
    ds.attrs["significance_pooling"] = ad.pooling_note
    ds.attrs["significance_source"] = (
        f"recomputed from the ingested raw members under "
        f"s3://{storage.BUCKET}/{storage.raw_prefix(ad.layer_id)} by "
        f"scripts/backfill_trend_significance.py")
    if out["recon"].get("n"):
        r = out["recon"]
        rag = "" if np.isnan(r["ragged_share"]) else \
            f"; {100*r['ragged_share']:.0f}% of the exceptions are cells whose " \
            f"member count varies within the decade, which makes the two " \
            f"estimators legitimately differ"
        ds.attrs["significance_reconstruction_check"] = (
            f"decade-mean of the annual ensemble-mean series vs published median "
            f"(post-baseline decades): median|diff|={r['median_abs']:.3e}, "
            f"max|diff|={r['max_abs']:.3e}, {r['n_exceed']:,} of {r['n']:,} "
            f"cell-decades over {RECON_RTOL:.0e} of the p99 scale "
            f"({100*r['frac_exceed']:.4f}%){rag}")
    return ds


NEW_VARS = ("trend_pvalue", "trend_tau", "trend_n_obs")

#: Variables this pass deliberately REPLACES rather than appends. `trend` moves
#: from the GUARDRAILS S10 baseline-anchored two-point rate to a Theil-Sen slope on
#: the ensemble-mean annual series (user decision 2026-07-30), so it is exempted
#: from the bit-identity gate — but only this one, and the exemption is explicit so
#: no other published value class can drift through it unnoticed.
REPLACED_VARS = ("trend",)


def replace_trend_with_sen(ds: xr.Dataset, sen: np.ndarray,
                           decades: List[int], window: int) -> xr.Dataset:
    """Overwrite `trend` with the Theil-Sen slope and restate its provenance.

    Two behaviour changes a consumer must be told about, both recorded in the attrs:

    * The baseline decade is now **NaN**, not exactly 0. A fitted slope needs an
      elapsed period; the old two-point rate was identically zero there by
      construction. This matches `trend_pvalue`/`trend_tau`.
    * ``trend x elapsed_decades == median[decade] - median[baseline]`` no longer
      holds. That identity was specific to the anchored rate. The replacement
      invariant is that ``sign(trend)`` agrees with ``sign(trend_tau)``, which is
      strictly stronger: both are now computed from the same series and window, so
      they can only disagree where the trajectory is genuinely non-monotonic.
    """
    ds["trend"] = (("decade", "lat", "lon"), sen.astype(np.float32))
    ds["trend"].attrs = {
        "long_name": "Theil-Sen slope of the decadal median series",
        "units": str(ds.attrs.get("trend_units", "")),
        "note": ("FITTED SLOPE over an expanding window anchored at the baseline "
                 "decade, in value per decade. The baseline panel is NaN (no "
                 "elapsed period), and trend x elapsed_decades no longer equals "
                 "the change map -- that identity belonged to the superseded "
                 "baseline-anchored rate. Fitted on the DECADAL series while "
                 "trend_pvalue is tested on the ANNUAL series, so the two describe "
                 "different series; see trend_definition."),
    }
    ds.attrs["trend_definition"] = tsig.trend_definition_decadal(
        decades, window_years=window, baseline_decade=decades[0])
    ds.attrs["trend_method"] = tsig.TREND_METHOD
    ds.attrs["trend_supersedes"] = (
        "baseline-anchored two-point rate (median[decade] - median[baseline]) / "
        "elapsed decades, used by this layer before 2026-07-30")
    return ds


def add_significance_vars(ds: xr.Dataset, p: np.ndarray, tau: np.ndarray,
                          nobs: np.ndarray, prefix: str = "",
                          of: str = "the layer's primary value") -> xr.Dataset:
    """Attach one p/tau/n_obs triple, optionally under a name prefix.

    `event100yr` publishes two values (event frequency and event footprint) and
    both need their own test, so the variable names are parameterized rather than
    fixed.
    """
    dims = ("decade", "lat", "lon")
    ds[f"{prefix}trend_pvalue"] = (dims, p.astype(np.float32))
    ds[f"{prefix}trend_tau"] = (dims, tau.astype(np.float32))
    ds[f"{prefix}trend_n_obs"] = (dims, nobs.astype(np.float32))
    ds[f"{prefix}trend_pvalue"].attrs = {
        "long_name": f"Two-sided Mann-Kendall p-value of the trend in {of}",
        "units": "1",
        "note": ("Expanding window anchored at the baseline decade. Tests the "
                 "ENSEMBLE-MEAN annual trajectory, NOT inter-model agreement -- "
                 "read with the CI and n_models. A constant series gives exactly 1.0."),
    }
    ds[f"{prefix}trend_tau"].attrs = {
        "long_name": f"Kendall tau-b of the trend in {of}",
        "units": "1",
        "note": ("Direction and rank-strength of the VALUE, not of the risk. For "
                 "risk direction combine sign(tau) with percentile_direction."),
    }
    ds[f"{prefix}trend_n_obs"].attrs = {
        "long_name": f"Number of years contributing to {prefix}trend_pvalue",
        "units": "1",
        "note": ("Sample size of the Mann-Kendall test. The baseline decade is 0 "
                 "(no elapsed period); a full record gives 20 at the first "
                 "post-baseline decade rising to 80 at the last."),
    }
    return ds


def publish_backfilled(layer_id: str, stage_root: Path, version: str,
                       prefix: str, on_exists: str) -> str:
    """Rebuild the manifest from the staged files and republish the version."""
    manifest = manifest_from_processed(
        stage_root, layer_id,
        created_by="scripts/backfill_trend_significance.py",
        notes=("Trend significance appended in place; see WORKFLOW-ISSUES.md "
               "2026-07-30 and GUARDRAILS.md S15."))
    try:
        old = storage.read_json(f"{prefix}layer.json")
        for carry in ("inputs", "supersedes"):
            if carry in old and carry not in manifest:
                manifest[carry] = old[carry]
    except FileNotFoundError:
        log("  WARNING: no existing layer.json to carry `inputs` forward from")
    return storage.publish_layer_version(
        layer_id, stage_root, manifest, version=version, on_exists=on_exists)


def carry_encoding(before: xr.Dataset, new_vars: Sequence[str]) -> dict:
    """Compression settings for the written file.

    Deliberately uniform rather than copying each source variable's own encoding.
    An earlier version preserved per-variable encoding to protect
    `scale_factor`/`add_offset` packing, but no layer uses packing — all nine store
    plain float32 with no scale/offset and a NaN `_FillValue` (verified 2026-07-30)
    — so it guarded a case that does not exist while making a re-run non-idempotent.
    Revisit only if a layer is ever published with packed variables.
    """
    return {v: {"zlib": True, "complevel": 4}
            for v in list(before.data_vars) + list(new_vars)}


def verify_append_only(before: xr.Dataset, after_path: Path,
                       new_vars: Sequence[str] = NEW_VARS,
                       replaced_vars: Sequence[str] = REPLACED_VARS) -> None:
    """Every pre-existing variable must survive bit-identical, bar named exceptions.

    The bit-identity gate is the contract with already-delivered numbers, so it is
    asserted rather than assumed — a silent dtype change or a reordered coordinate
    would rewrite published values while looking like a no-op.

    ``replaced_vars`` is the deliberate escape hatch (currently ``trend`` only).
    Those are still checked for shape and dtype, and required to have actually
    CHANGED — a "replacement" that silently no-ops would otherwise pass unnoticed
    and leave a layer claiming a Theil-Sen trend while holding the old rate.

    VALUES, shapes, dtypes and coordinates only. Attrs are deliberately not
    compared: on a re-run the significance variables are themselves present in
    ``before``, so comparing their attrs makes the pass non-idempotent the moment
    any label is reworded — which is exactly what happened on 2026-07-30.
    """
    with xr.open_dataset(after_path) as after:
        missing = [v for v in new_vars if v not in after.data_vars]
        if missing:
            raise ValueError(f"{after_path.name}: new variables absent: {missing}")
        for v in before.data_vars:
            a, b = before[v].values, after[v].values
            if a.shape != b.shape or a.dtype != b.dtype:
                raise ValueError(
                    f"{after_path.name}: {v} changed shape/dtype "
                    f"{a.shape}/{a.dtype} -> {b.shape}/{b.dtype}")
            if v in replaced_vars:
                # Idempotent re-runs legitimately produce an identical `trend`, so
                # only require a change when the layer still declares the old
                # method. Otherwise a no-op rerun would fail for being correct.
                was_anchored = "anchored" in str(
                    before.attrs.get("trend_definition", "")).lower()
                if was_anchored and np.array_equal(a, b, equal_nan=True):
                    raise ValueError(
                        f"{after_path.name}: {v} was meant to be REPLACED but is "
                        f"unchanged while the layer still declares the anchored "
                        f"rate -- the new trend was not written")
                continue
            if not np.array_equal(a, b, equal_nan=True):
                n = int(np.sum(a != b))
                raise ValueError(
                    f"{after_path.name}: {v} is NOT bit-identical "
                    f"({n:,} cells changed) -- refusing to publish")
        for c in before.coords:
            if not np.array_equal(before[c].values, after[c].values):
                raise ValueError(f"{after_path.name}: coordinate {c} changed")


def backfill(layer_id: str, ad: Adapter, check_only: bool = False,
             on_exists: str = "overwrite") -> int:
    log("=" * 78)
    log(f"{layer_id}")
    log("=" * 78)

    version = storage.resolve_current(layer_id)
    prefix = storage.version_prefix(layer_id, version).rstrip("/") + "/"
    marker = storage.verify_complete(prefix)
    n_art = len(marker.get("artifacts", {}))
    log(f"  version {version} (gate verified, {n_art} artifacts re-hashed)")

    local = storage.pull_prefix(f"{prefix}data/")
    published = sorted(Path(local).glob("*_processed.nc"))
    if not published:
        raise FileNotFoundError(f"{prefix}data/: no *_processed.nc")

    files = ad.files()
    by_scen: Dict[str, List[str]] = {}
    for f in files:
        by_scen.setdefault(ad.scenario_of(f), []).append(f)
    log(f"  raw members: {len(files)} across scenarios "
        f"{ {s: len(v) for s, v in sorted(by_scen.items())} }")

    stage = storage.staging_dir(layer_id) / "data" if not check_only else None
    if stage:
        stage.mkdir(parents=True, exist_ok=True)

    summary = []
    for path in published:
        scenario = path.name.split("_")[-2]
        with xr.open_dataset(path) as ds:
            ds.load()
        decades = [int(d) for d in ds.decade.values]
        window = int(ds.attrs.get("window_years", 10))
        scen_files = by_scen.get(scenario, [])
        if not scen_files:
            raise FileNotFoundError(
                f"{layer_id}: published scenario {scenario} has no raw members "
                f"(have {sorted(by_scen)}) -- cannot compute significance")
        log(f"\n  -- {scenario}: {len(scen_files)} members, decades "
            f"{decades[0]}s-{decades[-1]}s, window {window}y")

        out = significance_for_scenario(ad, scen_files, ds, decades, window)

        fin = np.isfinite(out["p"])
        old_trend = ds["trend"].values
        for i, d in enumerate(decades[1:], start=1):
            m = np.isfinite(out["p"][i])
            if not m.any():
                continue
            sig = float(np.mean(out["p"][i][m] < 0.05) * 100)
            up = float(np.mean(out["tau"][i][m] > 0) * 100)
            # Report how far the fitted slope moves the published rate, and whether
            # it agrees in sign — the check that replaces the retired
            # trend x elapsed == change identity.
            g = m & np.isfinite(out["sen"][i]) & np.isfinite(old_trend[i]) \
                & (out["sen"][i] != 0) & (old_trend[i] != 0)
            agree = float(np.mean(np.sign(out["sen"][i][g])
                                  == np.sign(old_trend[i][g])) * 100) if g.any() else float("nan")
            log(f"     {d}s: n={int(np.nanmax(out['n_obs'][i])):2d}  "
                f"p<0.05 {sig:5.1f}%  tau>0 {up:5.1f}%  cells {int(m.sum()):,}  "
                f"| sen vs old-anchored sign agree {agree:5.1f}%")
        summary.append((scenario, out["recon"],
                        float(np.mean(out["p"][fin] < 0.05) * 100) if fin.any() else 0.0))

        if ad.strict_recon and out["recon"].get("n"):
            r = out["recon"]
            if r["frac_exceed"] > RECON_MAX_EXCEED_FRAC:
                raise ValueError(
                    f"{layer_id}/{scenario}: the annual ensemble-mean series does "
                    f"not reproduce the published median -- {r['n_exceed']:,} "
                    f"cell-decades ({100*r['frac_exceed']:.3f}%) exceed "
                    f"{RECON_RTOL:.0e} of the p99 scale, limit "
                    f"{100*RECON_MAX_EXCEED_FRAC:.3f}%. Worst |diff| = "
                    f"{r['max_abs']:.3e}. The significance test would describe a "
                    f"different ensemble than the trend beside it -- check the "
                    f"adapter's value_scale and annualization before publishing.")

        if check_only:
            continue

        new = annotate(ds, out, ad, decades, window)
        target = stage / path.name
        new.to_netcdf(target, encoding=carry_encoding(ds, NEW_VARS))
        verify_append_only(ds, target)
        log(f"     staged {target.name} (+3 vars, trend REPLACED with Sen slope, "
            f"{len(ds.data_vars) - len(REPLACED_VARS)} carried through bit-identical)")

    if check_only:
        log("\n  CHECK-ONLY: nothing staged, nothing published.")
        return 0

    stage_root = stage.parent
    published_version = publish_backfilled(layer_id, stage_root, version, prefix,
                                           on_exists)
    log(f"\n  published {published_version} (overwrote in place)")
    log(f"  s3://{storage.BUCKET}/{storage.version_prefix(layer_id, published_version)}")

    for scenario, recon, sig in summary:
        log(f"    {scenario}: {sig:.1f}% of finite cell-decades significant at p<0.05")
    log("\n  NEXT: regenerate qa/ and maps/ and LOOK AT THEM (GUARDRAILS S11):")
    log(f"    python scripts/generate_qa_report.py {layer_id}")
    log(f"    python scripts/generate_maps.py {layer_id}")
    return 0


EVENT_LAYER = "riverflood_fldfrc-event100yr_annual"


def backfill_event100yr(check_only: bool = False,
                        on_exists: str = "overwrite") -> int:
    """Backfill the 1-in-100-year EVENT layer, which publishes TWO values.

    This layer has no raw of its own — it is derived from the paired `none` and
    `100yr` annual series — and its two values need separate tests because they
    move differently: the frequency more than doubles by the 2090s while the
    footprint is nearly static (+1.8%), which is the layer's headline claim. A
    single p-value would leave that claim unevidenced, and attaching the
    frequency's p-value to the footprint would assert the opposite of the truth.

    The footprint's series is ragged by construction: a year contributes only
    where the 1-in-100 threshold was exceeded, so `event_footprint_trend_n_obs`
    is much thinner than `trend_n_obs` and must be read alongside it.
    """
    mod = importlib.import_module("process_fldfrc_event100yr")
    log("=" * 78)
    log(f"{EVENT_LAYER}")
    log("=" * 78)

    version = storage.resolve_current(EVENT_LAYER)
    prefix = storage.version_prefix(EVENT_LAYER, version).rstrip("/") + "/"
    marker = storage.verify_complete(prefix)
    log(f"  version {version} (gate verified, "
        f"{len(marker.get('artifacts', {}))} artifacts re-hashed)")
    local = storage.pull_prefix(f"{prefix}data/")
    published = sorted(Path(local).glob("*_processed.nc"))

    files_n = {mod.member_key(str(p)): str(p)
               for p in storage.stage_raw(mod.SRC_UNPROTECTED, mod.RAW_PATTERN)}
    files_h = {mod.member_key(str(p)): str(p)
               for p in storage.stage_raw(mod.SRC_THRESHOLD, mod.RAW_PATTERN)}
    keys = sorted(set(files_n) & set(files_h))
    if not keys:
        raise FileNotFoundError(
            f"{EVENT_LAYER}: no matched none/100yr raw pairs "
            f"({len(files_n)} / {len(files_h)} found)")
    by_scen: Dict[str, List] = {}
    for k in keys:
        by_scen.setdefault(k[2], []).append(k)
    log(f"  matched pairs: {len(keys)} across "
        f"{ {s: len(v) for s, v in sorted(by_scen.items())} }")

    stage = (storage.staging_dir(EVENT_LAYER) / "data") if not check_only else None
    if stage:
        stage.mkdir(parents=True, exist_ok=True)
    years = np.arange(YEAR0, YEAR1 + 1)
    new_vars = list(NEW_VARS) + [f"event_footprint_{v}" for v in NEW_VARS]

    for path in published:
        scenario = path.name.split("_")[-2]
        with xr.open_dataset(path) as ds:
            ds.load()
        decades = [int(d) for d in ds.decade.values]
        window = int(ds.attrs.get("window_years", 10))
        scen_keys = by_scen.get(scenario)
        if not scen_keys:
            raise FileNotFoundError(
                f"{EVENT_LAYER}: published scenario {scenario} has no matched pairs")
        log(f"\n  -- {scenario}: {len(scen_keys)} pairs")

        shape = (len(years), ds.lat.size, ds.lon.size)
        f_sum = np.zeros(shape, np.float64)     # frequency indicator, %
        f_cnt = np.zeros(shape, np.float32)
        p_sum = np.zeros(shape, np.float64)     # footprint over exceeding years, %
        p_cnt = np.zeros(shape, np.float32)

        for k in scen_keys:
            na, exceeded, valid, yrs, _ = mod.load_pair(files_n[k], files_h[k])
            idx = {y: i for i, y in enumerate(yrs)}
            rows = [(t, idx[y]) for t, y in enumerate(years) if y in idx]
            for t, i in rows:
                v = valid[i]
                f_sum[t][v] += 100.0 * exceeded[i][v]
                f_cnt[t] += v
                e = exceeded[i]
                p_sum[t][e] += 100.0 * na[i][e]
                p_cnt[t] += e
            del na, exceeded, valid

        with np.errstate(invalid="ignore", divide="ignore"):
            freq = np.where(f_cnt > 0, f_sum / np.maximum(f_cnt, 1), np.nan)
            foot = np.where(p_cnt > 0, p_sum / np.maximum(p_cnt, 1), np.nan)
        freq = freq.astype(np.float32)
        foot = foot.astype(np.float32)

        recon = reconstruction_report(freq, years, ds["median"].values, decades,
                                     window, n_mem=f_cnt)
        log(f"    frequency reconstruction (scale={recon['scale_p99']:.4g}): "
            f"median|d|={recon['median_abs']:.3e} max={recon['max_abs']:.3e} "
            f"over-tol {recon['n_exceed']:,} ({100*recon['frac_exceed']:.4f}%, "
            f"limit {100*RECON_MAX_EXCEED_FRAC:.4f}%)")
        if recon["frac_exceed"] > RECON_MAX_EXCEED_FRAC:
            raise ValueError(
                f"{EVENT_LAYER}/{scenario}: the annual ensemble-mean FREQUENCY "
                f"series does not reproduce the published median "
                f"({recon['n_exceed']:,} cell-decades, "
                f"{100*recon['frac_exceed']:.3f}% over tolerance). Refusing to publish.")

        pv, tv, nv = tsig.mk_expanding(years, freq, decades, window_years=window,
                                       baseline_decade=decades[0])
        pf, tf, nf = tsig.mk_expanding(years, foot, decades, window_years=window,
                                       baseline_decade=decades[0])
        # `trend` on this layer describes the PRIMARY value (event frequency), and is
        # fitted on the published DECADAL medians of that frequency -- not the annual
        # series, which is the most zero-inflated of any layer (~93% exact zeros per
        # year) and would return an exactly-zero slope almost everywhere.
        sen = tsig.theilsen_decadal(ds["median"].values, decades,
                                    window_years=window,
                                    baseline_decade=decades[0])
        off = ~np.isfinite(ds["median"].values)
        for a in (pv, tv, sen):
            a[off] = np.nan
        nv[off] = 0
        off_f = ~np.isfinite(ds["event_footprint"].values)
        for a in (pf, tf):
            a[off_f] = np.nan
        nf[off_f] = 0

        for i, d in enumerate(decades[1:], start=1):
            mk, mf = np.isfinite(pv[i]), np.isfinite(pf[i])
            if not mk.any():
                continue
            log(f"     {d}s freq: n={int(np.nanmax(nv[i])):2d} "
                f"p<0.05 {100*np.mean(pv[i][mk] < 0.05):5.1f}% "
                f"tau>0 {100*np.mean(tv[i][mk] > 0):5.1f}%   |  footprint: "
                f"n(med)={int(np.nanmedian(nf[i][mf])) if mf.any() else 0:2d} "
                f"p<0.05 {100*np.mean(pf[i][mf] < 0.05) if mf.any() else 0:5.1f}% "
                f"tau>0 {100*np.mean(tf[i][mf] > 0) if mf.any() else 0:5.1f}%")

        if check_only:
            continue

        out = ds.copy()
        out = add_significance_vars(out, pv, tv, nv,
                                    of="event frequency (the primary value)")
        out = add_significance_vars(out, pf, tf, nf, prefix="event_footprint_",
                                    of="the conditional event footprint")
        out = replace_trend_with_sen(out, sen, decades, window)
        out["event_footprint_trend_n_obs"].attrs["note"] += (
            " The footprint series is RAGGED by construction -- a year contributes "
            "only where the 1-in-100 threshold was exceeded -- so this is much "
            "thinner than trend_n_obs and must be read alongside it.")
        out.attrs["significance_method"] = tsig.METHOD
        out.attrs["significance_definition"] = tsig.significance_definition(
            decades, window_years=window, baseline_decade=decades[0])
        out.attrs["significance_pooling"] = (
            "flat mean across members within each year, computed separately for the "
            "two published values. FREQUENCY: the per-year exceedance indicator, "
            "which reproduces the published median. FOOTPRINT: the undefended extent "
            "averaged over exceeding years only -- a CONDITIONAL mean whose "
            "contributing year set differs per member and per year, so it is NOT "
            "gated against the published event_footprint (the two estimators "
            "legitimately differ) and its n_obs is much thinner.")
        out.attrs["significance_source"] = (
            f"recomputed from the paired {mod.SRC_UNPROTECTED} / "
            f"{mod.SRC_THRESHOLD} raw members by "
            f"scripts/backfill_trend_significance.py")
        out.attrs["significance_reconstruction_check"] = (
            f"frequency: median|diff|={recon['median_abs']:.3e}, "
            f"max|diff|={recon['max_abs']:.3e}, {recon['n_exceed']:,} of "
            f"{recon['n']:,} cell-decades over {RECON_RTOL:.0e} of the p99 scale")

        target = stage / path.name
        out.to_netcdf(target, encoding=carry_encoding(ds, new_vars))
        verify_append_only(ds, target, new_vars=new_vars)
        log(f"     staged {target.name} (+6 vars, {len(ds.data_vars)} carried "
            f"through bit-identical)")

    if check_only:
        log("\n  CHECK-ONLY: nothing staged, nothing published.")
        return 0
    v = publish_backfilled(EVENT_LAYER, stage.parent, version, prefix, on_exists)
    log(f"\n  published {v} (overwrote in place)")
    return 0


def main() -> int:
    ads = build_adapters()
    choices = sorted(list(ads) + [EVENT_LAYER])
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("layer_id", nargs="?", choices=choices,
                    help="layer to backfill")
    ap.add_argument("--all", action="store_true", help="backfill every layer")
    ap.add_argument("--check-only", action="store_true",
                    help="compute and report, stage nothing")
    ap.add_argument("--on-exists", default="overwrite",
                    choices=("overwrite", "bump", "error"),
                    help="what to do about the existing version (default overwrite)")
    args = ap.parse_args()

    if not args.all and not args.layer_id:
        ap.error("give a layer_id or --all")
    targets = choices if args.all else [args.layer_id]

    rc, failed = 0, []
    for lid in targets:
        try:
            if lid == EVENT_LAYER:
                rc |= backfill_event100yr(check_only=args.check_only,
                                          on_exists=args.on_exists)
            else:
                rc |= backfill(lid, ads[lid], check_only=args.check_only,
                               on_exists=args.on_exists)
        except Exception as exc:                      # keep going; report at the end
            log(f"\n  FAILED {lid}: {type(exc).__name__}: {exc}")
            failed.append(lid)
            rc = 1
        log("")
    if failed:
        # Print the roll-up last: with 9 layers the individual failures scroll off,
        # and a run that half-succeeded must not look like a success.
        log(f"FAILED {len(failed)} of {len(targets)}: {', '.join(failed)}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
