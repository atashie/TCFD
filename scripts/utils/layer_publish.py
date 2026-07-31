"""Derive a layer manifest from processed NetCDF files and publish to S3.

The processors already record every processing decision as NetCDF global
attributes (`statistic`, `percentile_direction`, `normalization`,
`spatial_smoothing`, `trend_definition`, `baseline_*`, `n_members`, …). This
module reads those back rather than asking each processor to restate them, so
`layer.json` cannot drift out of sync with the data it describes.

All global attributes are carried through verbatim under `source_attrs`; the
subset that GUARDRAILS.md treats as load-bearing is additionally projected into
a structured `decisions` block. Carrying everything means a processor can add a
new attribute without silently losing provenance here.

Usage from a processor::

    from utils.layer_publish import publish_processed_layer

    version = publish_processed_layer(
        layer_id="wildfire_burntarea_annual",
        stage_dir=stage,
        created_by="scripts/process_burntarea_fire.py",
        raw_entries=raw_entries,
    )
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402

#: Global attrs projected into the structured `decisions` block.
#: Keys are manifest field names; values are the NetCDF attribute they come from.
DECISION_ATTRS = {
    "decadal_statistic": "statistic",
    "normalization": "normalization",
    "spatial_smoothing": "spatial_smoothing",
    "percentile_direction": "percentile_direction",
    "trend_definition": "trend_definition",
    "trend_units": "trend_units",
    "significance_method": "significance_method",
    "significance_definition": "significance_definition",
    "significance_pooling": "significance_pooling",
    "baseline_decade": "baseline_decade",
    "baseline_source": "baseline_source",
    "window_years": "window_years",
}

#: Attributes that legitimately differ between scenario files.
_SCENARIO_SPECIFIC = {"scenario", "history", "description", "percentile_zero_fraction",
                      # Each scenario's reconstruction residual is its own diagnostic;
                      # requiring agreement would reject a correctly built layer.
                      "significance_reconstruction_check"}

_PROCESSED_RE = re.compile(r"^(?P<var>.+)_(?P<scenario>[a-z0-9]+)_processed\.nc$")


def _percentile_mode(attrs: Dict[str, Any]) -> Optional[str]:
    """Read the two-tier / single-tier percentile mode out of `percentile_baseline`."""
    baseline = str(attrs.get("percentile_baseline", ""))
    for mode in ("two_tier", "single_tier"):
        if mode in baseline:
            return mode
    return None


def _split_list(value: Any) -> Any:
    """Turn a comma-joined attribute back into a list, leaving other types alone."""
    if isinstance(value, str) and "," in value:
        return [part.strip() for part in value.split(",") if part.strip()]
    return value


def manifest_from_processed(
    stage_dir: Path,
    layer_id: str,
    created_by: str,
    raw_entries: Optional[Sequence[Dict[str, Any]]] = None,
    notes: Optional[str] = None,
    supersedes: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a `layer.json` manifest from the processed files in `stage_dir/data`.

    Args:
        stage_dir: Staging directory from :func:`storage.staging_dir`.
        layer_id: Canonical layer id.
        created_by: Script path that produced the outputs.
        raw_entries: Entries from :func:`storage.ingest_raw`, recording each
            input's name, size, checksum and origin URL.
        notes: Free-text pointer to the relevant WORKFLOW-ISSUES / GUARDRAILS entries.
        supersedes: Version id this one replaces, if any.

    Returns:
        A manifest dict ready for :func:`storage.publish_layer_version`.

    Raises:
        FileNotFoundError: If no processed files are staged.
        ValueError: If scenario files disagree on a non-scenario-specific attribute,
            which would mean the ensemble was not processed uniformly.
    """
    data_dir = Path(stage_dir) / "data"
    files = sorted(data_dir.glob("*_processed.nc"))
    if not files:
        raise FileNotFoundError(f"{data_dir}: no *_processed.nc files staged")

    scenarios: List[str] = []
    variable: Optional[str] = None
    merged: Dict[str, Any] = {}
    conflicts: Dict[str, set] = {}

    for path in files:
        match = _PROCESSED_RE.match(path.name)
        if match:
            scenarios.append(match.group("scenario"))
            variable = variable or match.group("var")

        with xr.open_dataset(path) as ds:
            attrs = dict(ds.attrs)

        for key, value in attrs.items():
            if key in _SCENARIO_SPECIFIC:
                continue
            if key in merged and merged[key] != value:
                conflicts.setdefault(key, set()).add(str(merged[key]))
                conflicts[key].add(str(value))
            merged.setdefault(key, value)

    if conflicts:
        detail = "; ".join(f"{k}: {sorted(v)}" for k, v in conflicts.items())
        raise ValueError(
            f"{layer_id}: scenario files disagree on {detail} — the ensemble was "
            f"not processed uniformly, so one manifest cannot describe it"
        )

    decisions = {
        field: merged[attr] for field, attr in DECISION_ATTRS.items() if attr in merged
    }

    # GUARDRAILS S10 (from 2026-07-30): `trend` must be a Theil-Sen slope of the
    # ensemble-mean annual series. Warn rather than raise, so a layer already in
    # flight under the old definition can still be published — but make the gap
    # loud and record it in the manifest, because the alternative is a layer
    # shipping a trend that no longer matches what the docs and the p-value claim.
    tm = str(merged.get("trend_method", ""))
    if "theil_sen" not in tm.lower():
        decisions["trend_method"] = tm or "MISSING"
        print(f"  WARNING: {layer_id} declares trend_method={tm or 'MISSING'!r}, not "
              f"theil_sen_on_ensemble_mean_annual_series. GUARDRAILS S10 requires the "
              f"Theil-Sen slope of the ensemble-mean annual series; use "
              f"utils.trend_significance.theilsen_expanding(). Publishing anyway.",
              flush=True)
    mode = _percentile_mode(merged)
    if mode:
        decisions["percentile_mode"] = mode

    manifest: Dict[str, Any] = {
        "layer_id": layer_id,
        "created_by": created_by,
        "variable": merged.get("variable", variable),
        "units": merged.get("units"),
        "cadence": "monthly" if layer_id.endswith("_monthly") else "annual",
        "scenarios": sorted(set(scenarios)),
        "source_dataset": merged.get("source_dataset"),
        "decisions": decisions,
        "ensemble": {
            "n_members_per_scenario": merged.get("n_members"),
            "impact_models": _split_list(merged.get("impact_models")),
            "gcms": _split_list(merged.get("gcms")),
        },
        "source_attrs": merged,
    }

    if raw_entries is not None:
        manifest["inputs"] = {
            "raw_prefix": storage.raw_prefix(layer_id),
            "files": list(raw_entries),
        }
    if notes:
        manifest["notes"] = notes
    if supersedes:
        manifest["supersedes"] = supersedes

    return manifest


def publish_processed_layer(
    layer_id: str,
    stage_dir: Path,
    created_by: str,
    raw_entries: Optional[Sequence[Dict[str, Any]]] = None,
    notes: Optional[str] = None,
    on_exists: str = "error",
) -> str:
    """Derive the manifest and publish the staged layer version to S3.

    Returns:
        The published version id.
    """
    previous: Optional[str] = None
    try:
        previous = storage.resolve_current(layer_id)
    except FileNotFoundError:
        pass

    manifest = manifest_from_processed(
        stage_dir,
        layer_id,
        created_by=created_by,
        raw_entries=raw_entries,
        notes=notes,
        supersedes=previous,
    )
    version = storage.publish_layer_version(
        layer_id, stage_dir, manifest, on_exists=on_exists
    )
    print(f"  published s3://{storage.BUCKET}/"
          f"{storage.version_prefix(layer_id, version)}", flush=True)
    print(f"  current -> s3://{storage.BUCKET}/{storage.current_prefix(layer_id)}",
          flush=True)
    return version
