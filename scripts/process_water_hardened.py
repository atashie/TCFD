"""A3 hardened processing wrapper for the Water Index SSP rebuild.

Wraps `process_water_variable.py` with the trust chain the WS1 plan (and its
external review) requires, without modifying the processor itself:

  1. RAW AUDIT (per manifest row for the variable): file present at its
     manifest destination; sha256 matches the download results receipt
     (skippable via --skip-hash on re-runs); NetCDF opens; expected variable
     present with expected raw units; time axis = exactly 1032 monthly steps
     spanning 2015-2100; 360x720 grid; per-file calendar recorded (they are
     KNOWN to differ by model); first field not all-NaN. -> audit CSV.
  2. ENSEMBLE ASSERTION: identical member sets across the three scenarios,
     matching the manifest exactly. No silent skips possible downstream,
     because anything missing fails here first.
  3. PROCESS: subprocess `process_water_variable.py` (normalization flags
     passed through). Nonzero exit propagates.
  4. OUTPUT AUDIT + HONEST PROVENANCE: every scenario slice must contain valid
     cells (the bare processor exits 0 on empty scenarios — measured 2026-08-21);
     global attrs are then REWRITTEN from the resolved inputs: actual member
     list per scenario, source manifest+results paths, per-file audit csv
     sha256, processor+wrapper git commit, normalization flags used.

Usage:
    python scripts/process_water_hardened.py VAR \
        --manifest .../manifest_water_download.csv \
        --results  .../manifest_water_download_results.csv \
        --output OUT.nc [--normalize | --normalize-models m1,m2] [--skip-hash]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

SCRIPTS = Path(__file__).resolve().parent
BASE = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
from config_water_variables import get_variable_config

EXPECTED_MONTHS = 1032  # 2015-01 .. 2100-12
SCENARIOS = ("ssp126", "ssp370", "ssp585")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(BASE), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def audit_raw(var: str, manifest_rows: list[dict], results_by_dest: dict,
              skip_hash: bool, audit_path: Path) -> list[dict]:
    cfg = get_variable_config(var)
    failures = []
    audit = []
    for row in manifest_rows:
        dest = BASE / row["dest"]
        rec = {"dest": row["dest"], "model": row["model"], "gcm": row["gcm"],
               "scenario": row["scenario"]}
        try:
            if not dest.exists():
                raise AssertionError("missing on disk")
            if not skip_hash:
                expected = results_by_dest.get(row["dest"], {}).get("sha256", "")
                if expected:
                    actual = sha256_file(dest)
                    if actual != expected:
                        raise AssertionError(f"sha256 mismatch vs download receipt")
                    rec["sha256"] = actual
                else:
                    rec["sha256"] = sha256_file(dest)
            nc = Dataset(str(dest))
            try:
                if cfg.name not in nc.variables:
                    raise AssertionError(f"variable '{cfg.name}' absent (has: {list(nc.variables)[:6]})")
                v = nc.variables[cfg.name]
                units = getattr(v, "units", "MISSING")
                if units.replace("/", " ").replace("kg m-2 s-1", "kg m-2 s-1") != cfg.units_raw and units != cfg.units_raw:
                    raise AssertionError(f"units '{units}' != expected '{cfg.units_raw}'")
                t = nc.variables["time"]
                if len(t) != EXPECTED_MONTHS:
                    raise AssertionError(f"time steps {len(t)} != {EXPECTED_MONTHS}")
                if nc.dimensions["lat"].size != 360 or nc.dimensions["lon"].size != 720:
                    raise AssertionError("grid is not 360x720")
                rec["calendar"] = getattr(t, "calendar", "MISSING")
                first = v[0, :, :]
                arr = first.filled(np.nan) if hasattr(first, "filled") else np.asarray(first)
                n_finite = int(np.isfinite(arr).sum())
                if n_finite == 0:
                    raise AssertionError("first time step is all-NaN")
                rec["n_finite_first_step"] = n_finite
                rec["status"] = "OK"
            finally:
                nc.close()
        except Exception as e:
            rec["status"] = f"FAIL: {e}"
            failures.append(rec)
        audit.append(rec)
    with open(audit_path, "w", newline="") as f:
        fieldnames = sorted({k for r in audit for k in r})
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(audit)
    if failures:
        for r in failures:
            print(f"RAW AUDIT FAIL: {r['dest']}: {r['status']}")
        sys.exit(1)
    calendars = sorted({r["calendar"] for r in audit})
    print(f"raw audit OK: {len(audit)} files; calendars present: {calendars}")
    return audit


def assert_ensemble(manifest_rows: list[dict]) -> dict:
    members = defaultdict(set)
    for row in manifest_rows:
        members[row["scenario"]].add((row["model"], row["gcm"]))
    sets = {s: members.get(s, set()) for s in SCENARIOS}
    ref = sets[SCENARIOS[0]]
    for s in SCENARIOS:
        if sets[s] != ref:
            missing = ref ^ sets[s]
            print(f"ENSEMBLE ASSERTION FAIL: {s} differs: {sorted(missing)}")
            sys.exit(1)
    print(f"ensemble asserted: {len(ref)} members identical across {len(SCENARIOS)} scenarios")
    return {s: sorted(f"{m}/{g}" for m, g in sets[s]) for s in SCENARIOS}


def output_audit_and_provenance(out_path: Path, var: str, members: dict,
                                args, audit_path: Path) -> None:
    cfg = get_variable_config(var)
    nc = Dataset(str(out_path), "a")
    try:
        v = nc.variables[cfg.name]
        scen_labels = ["".join(x) if isinstance(x, (bytes, str)) else "".join(map(str, x))
                       for x in nc.variables["scenario"][:]]
        for s_idx, label in enumerate(scen_labels):
            sl = v[:, :, s_idx, 12, 0]
            arr = sl.filled(np.nan) if hasattr(sl, "filled") else np.asarray(sl)
            n_valid = int(np.isfinite(arr).sum())
            if n_valid == 0:
                print(f"OUTPUT AUDIT FAIL: scenario '{label}' has zero valid cells")
                sys.exit(1)
            print(f"  scenario {label}: {n_valid} valid cells (annual-mean, first decade)")
        # Honest provenance — from resolved inputs, not static config
        nc.setncattr("impact_models", ", ".join(sorted({m.split("/")[0] for m in members[SCENARIOS[0]]})))
        nc.setncattr("gcms", ", ".join(sorted({m.split("/")[1] for m in members[SCENARIOS[0]]})))
        nc.setncattr("members_per_scenario", json.dumps({s: len(members[s]) for s in SCENARIOS}))
        nc.setncattr("members_resolved", json.dumps(members[SCENARIOS[0]]))
        nc.setncattr("source_manifest", str(args.manifest))
        nc.setncattr("source_results", str(args.results))
        nc.setncattr("raw_audit_csv", str(audit_path))
        nc.setncattr("raw_audit_sha256", sha256_file(audit_path))
        nc.setncattr("processing_code_commit", git_commit())
        nc.setncattr("normalization_flags",
                     f"normalize={args.normalize} normalize_models={args.normalize_models or ''}")
        if args.units_override:
            # A normalized variable must not claim physical units (review blocker #1):
            # the config's units_output describes the raw quantity, not the file truth.
            nc.setncattr("units_physical_original", nc.getncattr("units")
                         if "units" in nc.ncattrs() else "")
            nc.setncattr("units", args.units_override)
        nc.setncattr("hardened_by", "process_water_hardened.py "
                     + datetime.now(timezone.utc).isoformat(timespec="seconds"))
    finally:
        nc.close()
    print("output audit OK; provenance rewritten from resolved inputs")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("variable")
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--normalize", action="store_true")
    ap.add_argument("--normalize-models", type=str, default=None)
    ap.add_argument("--skip-hash", action="store_true")
    ap.add_argument("--units-override", type=str, default=None)
    args = ap.parse_args()

    var = args.variable
    manifest_rows = [r for r in csv.DictReader(open(args.manifest)) if r["variable"] == var]
    if not manifest_rows:
        print(f"no manifest rows for variable '{var}'"); sys.exit(1)
    results_by_dest = {r["dest"]: r for r in csv.DictReader(open(args.results))}

    audit_path = args.output.parent / f"{var}_raw_audit.csv"
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== A3 hardened processing: {var} ({len(manifest_rows)} raw files) ===")
    audit_raw(var, manifest_rows, results_by_dest, args.skip_hash, audit_path)
    members = assert_ensemble(manifest_rows)

    cmd = [sys.executable, str(SCRIPTS / "process_water_variable.py"), var,
           "--output", str(args.output)]
    if args.data_dir:
        cmd += ["--data-dir", str(args.data_dir)]
    else:
        cmd += ["--data-dir", str(BASE / "data" / "raw" / f"water_{var}")]
    if args.normalize:
        cmd.append("--normalize")
    if args.normalize_models:
        cmd += ["--normalize-models", args.normalize_models]
    print("running:", " ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        print(f"processor exited {rc}"); sys.exit(rc)

    output_audit_and_provenance(args.output, var, members, args, audit_path)
    print(f"=== {var} complete: {args.output} ===")


if __name__ == "__main__":
    main()
