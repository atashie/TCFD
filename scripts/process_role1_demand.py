"""Role-1 demand-stack processor (WRI WS2 Role 1; Gate G10-E mechanics).

Turns the twelve raw ISIMIP potential withdrawal/consumption variables
(manifest: waterRiskIndex_beta/v2/enumeration/ROLE1_20260831/) into the
per-variable demand product the WRI Engine-A stack build consumes:

    data/processed/ws2_role1/demand_{var}_2015soc.nc
    dims (member, scenario, decade, value_type, lat, lon)

Value types (13, per Gate G10-E — decided 2026-08-28: monthly climatology +
annual mean per decade, NO demand quantiles):
    vt 0-11  climatological monthly demand VOLUME, km3/month (mean over the
             decade's years of that calendar month's volume)
    vt 12    annual demand volume, km3/yr (mean over years of the annual sum;
             equals the sum of vt 0-11 by linearity — asserted at validation)

Volumetric conversion follows the pelec pattern (build_pelec_baseline.py):
flux [kg m-2 s-1] x seconds-in-month x cell_area [m2] x 1e-12 -> km3/month,
with seconds-in-month taken from EACH FILE'S OWN CALENDAR (GUARDRAILS
§18.3: calendars are heterogeneous across models; these files carry no time
bounds and no cell_methods, so "monthly values are monthly means" is a
recorded assumption, not a read fact).

This is NOT a water-index family file: GUARDRAILS §6/§7 (20-value-type
contract) govern `waterIndexUnderlyingData_*` only. This product follows the
ws2_role3 pelec precedent — a WRI Product-2 dataset, never wired into the
TCFD delivery layer registry (two-product separation, GUARDRAILS §6).

Trust chain per variable (process_water_hardened.py discipline, integrated):
  1. RAW AUDIT: file at manifest dest; sha256 vs download receipt
     (--skip-hash to skip); NetCDF opens; expected variable present AND the
     sole data variable; units exactly 'kg m-2 s-1'; 1032 monthly steps;
     360x720 grid; per-file calendar recorded; first step not all-NaN;
     time axis asserted to be exactly 2015-01..2100-12 in order.
  2. ENSEMBLE ASSERTION: identical member sets across the three scenarios,
     matching the manifest exactly.
  3. PROCESS: per member x scenario, decade climatologies (2010s = 2015-19
     stub, then 2020s..2090s; 2100 dropped — family convention).
  4. VALIDATE + HONEST PROVENANCE: per panel — finite cells > 0, all values
     >= 0, vt12 == sum(vt0-11) to float tolerance; per-member valid-cell
     census and global totals printed; provenance rewritten from resolved
     inputs (member roster, calendars, soc tokens, manifest+results+audit
     sha256, code commit).

Usage:
    .venv/bin/python scripts/process_role1_demand.py \
        --manifest .../manifest_role1_demand_download.csv \
        --results  .../manifest_role1_demand_download_results.csv \
        [--variables pdomww,pdomuse] [--skip-hash]
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
from netCDF4 import Dataset, num2date

SCRIPTS = Path(__file__).resolve().parent
BASE = SCRIPTS.parent

SCENARIOS = ("ssp126", "ssp370", "ssp585")
DECADES = (2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090)
MIN_YEAR, MAX_YEAR = 2015, 2099
EXPECTED_MONTHS = 1032  # 2015-01 .. 2100-12
EXPECTED_UNITS = "kg m-2 s-1"
N_VT = 13
VALUE_TYPE_NAMES = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
    "Nov", "Dec", "Annual_Total",
]
DAYS_365 = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def decade_years(decade: int) -> tuple[int, int]:
    return (2015, 2019) if decade == 2010 else (decade, decade + 9)


def days_in_month(year: int, month: int, calendar: str) -> float:
    cal = calendar.lower()
    if cal == "360_day":
        return 30.0
    if cal in ("365_day", "noleap"):
        return float(DAYS_365[month - 1])
    if cal in ("standard", "gregorian", "proleptic_gregorian"):
        if month == 2:
            leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
            return 29.0 if leap else 28.0
        return float(DAYS_365[month - 1])
    # 'julian' deliberately unhandled (review 2026-09-02 L11): its century
    # leap rule differs from Gregorian; refuse rather than silently misdate.
    raise AssertionError(f"unhandled calendar '{calendar}'")


def cell_area_m2(lat: np.ndarray, n_lon: int, dlat: float = 0.5,
                 dlon: float = 0.5) -> np.ndarray:
    """Spherical-band cell area, (lat, lon) in m^2 (pelec pattern)."""
    r = 6_371_000.0
    lat_rad = np.deg2rad(lat)
    band = (np.sin(lat_rad + np.deg2rad(dlat / 2))
            - np.sin(lat_rad - np.deg2rad(dlat / 2)))
    area_band = 2 * np.pi * r * r * band / (360.0 / dlon)
    return np.repeat(area_band[:, None], n_lon, axis=1)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(BASE), "rev-parse", "--short", "HEAD"],
            text=True).strip()
    except Exception:
        return "unknown"


def audit_raw(var: str, rows: list[dict], results_by_dest: dict,
              skip_hash: bool, audit_path: Path) -> dict:
    """Stage 1. Returns {dest: calendar}; exits nonzero on any failure."""
    failures, audit, calendars = [], [], {}
    for row in rows:
        dest = BASE / row["dest"]
        rec = {"dest": row["dest"], "model": row["model"], "gcm": row["gcm"],
               "scenario": row["scenario"]}
        try:
            if not dest.exists():
                raise AssertionError("missing on disk")
            if not skip_hash:
                expected = results_by_dest.get(row["dest"], {}).get("sha256", "")
                actual = sha256_file(dest)
                if expected and actual != expected:
                    raise AssertionError("sha256 mismatch vs download receipt")
                rec["sha256"] = actual
                # Review 2026-09-02 H2: an empty receipt digest is NOT a
                # verification — record it loudly so the audit CSV and log
                # show exactly which files are self-hashed only.
                rec["receipt_digest"] = "present" if expected else "ABSENT"
            nc = Dataset(str(dest))
            try:
                if var not in nc.variables:
                    raise AssertionError(
                        f"variable '{var}' absent (has: {list(nc.variables)[:6]})")
                dims = set(nc.dimensions)
                data_vars = [n for n in nc.variables if n not in dims
                             and n not in ("lat_bnds", "lon_bnds", "time_bnds")]
                if set(data_vars) != {var}:
                    raise AssertionError(
                        f"expected sole data variable '{var}', file has {data_vars}")
                v = nc.variables[var]
                units = getattr(v, "units", "MISSING").strip()
                if units != EXPECTED_UNITS:
                    raise AssertionError(
                        f"units '{units}' != expected '{EXPECTED_UNITS}'")
                t = nc.variables["time"]
                if len(t) != EXPECTED_MONTHS:
                    raise AssertionError(f"time steps {len(t)} != {EXPECTED_MONTHS}")
                if nc.dimensions["lat"].size != 360 or nc.dimensions["lon"].size != 720:
                    raise AssertionError("grid is not 360x720")
                lat0 = float(nc.variables["lat"][0])
                if abs(lat0 - 89.75) > 1e-6:
                    raise AssertionError(f"lat[0]={lat0} != 89.75 (orientation)")
                cal = getattr(t, "calendar", "standard")
                # Time-axis identity: exactly 2015-01..2100-12 in order.
                dates = num2date(t[:], t.units, calendar=cal)
                ym = [(d.year, d.month) for d in dates]
                expected_ym = [(y, m) for y in range(2015, 2101)
                               for m in range(1, 13)]
                if ym != expected_ym:
                    raise AssertionError(
                        f"time axis is not 2015-01..2100-12 monthly "
                        f"(first={ym[0]}, last={ym[-1]})")
                rec["calendar"] = cal
                calendars[row["dest"]] = cal
                first = v[0, :, :]
                arr = first.filled(np.nan) if hasattr(first, "filled") \
                    else np.asarray(first)
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
        w.writeheader()
        w.writerows(audit)
    if failures:
        for r in failures:
            print(f"RAW AUDIT FAIL: {r['dest']}: {r['status']}")
        sys.exit(1)
    n_absent = sum(1 for r in audit if r.get("receipt_digest") == "ABSENT")
    if n_absent:
        print(f"  WARNING: {n_absent}/{len(audit)} files have NO download-"
              f"receipt sha256 (verified-existing path) — self-hashed only; "
              f"see the audit CSV's receipt_digest column")
    print(f"  raw audit OK: {len(audit)} files; calendars: "
          f"{sorted(set(calendars.values()))}")
    return calendars


def assert_ensemble(rows: list[dict]) -> list[tuple[str, str]]:
    """Stage 2. Returns the sorted member list; exits nonzero on mismatch."""
    members = defaultdict(set)
    for row in rows:
        members[row["scenario"]].add((row["model"], row["gcm"]))
    ref = members[SCENARIOS[0]]
    for s in SCENARIOS:
        if members[s] != ref:
            print(f"ENSEMBLE ASSERTION FAIL: {s} differs: "
                  f"{sorted(ref ^ members[s])}")
            sys.exit(1)
    print(f"  ensemble asserted: {len(ref)} members identical across "
          f"{len(SCENARIOS)} scenarios")
    return sorted(ref)


def process_member_file(path: Path, var: str, area_m2: np.ndarray,
                        calendar: str) -> np.ndarray:
    """One raw file -> (9 decades, 13 vt, 360, 720) float32 volumes."""
    out = np.full((len(DECADES), N_VT, 360, 720), np.nan, dtype=np.float32)
    nc = Dataset(str(path))
    try:
        v = nc.variables[var]
        for d_idx, dec in enumerate(DECADES):
            y0, y1 = decade_years(dec)
            i0, i1 = (y0 - 2015) * 12, (y1 - 2015) * 12 + 12
            block = v[i0:i1, :, :]
            block = block.filled(np.nan) if hasattr(block, "filled") \
                else np.asarray(block, dtype=np.float64)
            n_years = y1 - y0 + 1
            block = block.reshape(n_years, 12, 360, 720)
            sec = np.array([[days_in_month(y, m, calendar) * 86400.0
                             for m in range(1, 13)]
                            for y in range(y0, y1 + 1)])
            vol = block * sec[:, :, None, None] * area_m2[None, None, :, :] * 1e-12
            out[d_idx, 0:12] = vol.mean(axis=0)          # km3/month climatology
            out[d_idx, 12] = vol.sum(axis=1).mean(axis=0)  # km3/yr annual mean
    finally:
        nc.close()
    return out


def write_output(out_path: Path, var: str, members: list[tuple[str, str]],
                 data: np.ndarray, area_m2: np.ndarray, provenance: dict) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds = Dataset(str(out_path), "w")
    try:
        ds.createDimension("member", len(members))
        ds.createDimension("scenario", len(SCENARIOS))
        ds.createDimension("decade", len(DECADES))
        ds.createDimension("value_type", N_VT)
        ds.createDimension("lat", 360)
        ds.createDimension("lon", 720)
        lat = ds.createVariable("lat", "f8", ("lat",))
        lat[:] = np.arange(89.75, -90, -0.5)
        lon = ds.createVariable("lon", "f8", ("lon",))
        lon[:] = np.arange(-179.75, 180, 0.5)
        dv = ds.createVariable("decade", "i4", ("decade",))
        dv[:] = DECADES
        mv = ds.createVariable("member", str, ("member",))
        for i, (model, gcm) in enumerate(members):
            mv[i] = f"{model}/{gcm}"
        sv = ds.createVariable("scenario", str, ("scenario",))
        for i, s in enumerate(SCENARIOS):
            sv[i] = s
        va = ds.createVariable("cell_area", "f8", ("lat", "lon"),
                               zlib=True, complevel=4)
        va[:] = area_m2
        va.units, va.long_name = "m2", "grid cell area (spherical band)"
        vv = ds.createVariable(f"{var}_km3", "f4",
                               ("member", "scenario", "decade", "value_type",
                                "lat", "lon"),
                               zlib=True, complevel=4,
                               fill_value=np.float32(np.nan))
        vv[:] = data
        vv.long_name = f"{var} demand volume"
        vv.units = "km3/month (vt 0-11), km3/yr (vt 12)"
        vv.value_type_names = json.dumps(VALUE_TYPE_NAMES)
        for k, val in provenance.items():
            ds.setncattr(k, val if isinstance(val, str) else json.dumps(val))
    finally:
        ds.close()


def validate_output(out_path: Path, var: str,
                    members: list[tuple[str, str]]) -> None:
    """Stage 4 (fatal): per-panel finiteness, nonnegativity, vt12 identity.

    Nonnegativity contract (extended DELIBERATELY 2026-09-03, never
    silently): potential demand (p*) can never be negative — fatal.
    ACTUAL consumption (a*use) is computed by the models as a residual
    (abstraction − return flows) and is legitimately negative in
    return-dominated cells — measured on WaterGAP atotuse: 1.54% of
    finite values, clustered in canal-irrigation regions (Punjab ~30N
    73-74E), −5.8 vs +86 km³ in the probe month. For a* variables the
    negatives are ALLOWED and a per-member census is printed and stamped
    into provenance by main(); wholesale negativity (>10% of finite
    values) still refuses.
    """
    # exactly the documented a*use class — a future a*ww must NOT inherit
    # the exception silently (lane-review H3)
    allow_negative = var.startswith("a") and var.endswith("use")
    neg_census = {}
    ds = Dataset(str(out_path))
    try:
        vv = ds.variables[f"{var}_km3"]
        for mi, (model, gcm) in enumerate(members):
            for si, scen in enumerate(SCENARIOS):
                panel = vv[mi, si]
                panel = panel.filled(np.nan) if hasattr(panel, "filled") \
                    else np.asarray(panel)
                for di, dec in enumerate(DECADES):
                    vt12 = panel[di, 12]
                    finite = np.isfinite(vt12)
                    if int(finite.sum()) == 0:
                        print(f"VALIDATE FAIL: {var} {model}/{gcm} {scen} "
                              f"{dec}s: zero finite cells")
                        sys.exit(1)
                    # Review 2026-09-02 M7: a finite annual cell with missing
                    # monthly members must not pass — assert the finite mask
                    # is IDENTICAL across all 13 value types.
                    for vt in range(12):
                        if not np.array_equal(np.isfinite(panel[di, vt]),
                                              finite):
                            print(f"VALIDATE FAIL: {var} {model}/{gcm} "
                                  f"{scen} {dec}s: finite mask of vt{vt} "
                                  f"differs from vt12")
                            sys.exit(1)
                    n_neg = int((panel[di] < 0).sum())
                    if n_neg and not allow_negative:
                        print(f"VALIDATE FAIL: {var} {model}/{gcm} {scen} "
                              f"{dec}s: negative volume")
                        sys.exit(1)
                    if n_neg:
                        # the wholesale gate runs PER VALUE TYPE — pooling
                        # all 13 lets an all-negative annual layer hide at
                        # 1/13 = 7.7% (lane-review H3)
                        for vt in range(13):
                            layer = panel[di, vt]
                            nv = int((layer < 0).sum())
                            nf = int(np.isfinite(layer).sum())
                            if nf and nv > 0.10 * nf:
                                print(f"VALIDATE FAIL: {var} {model}/{gcm} "
                                      f"{scen} {dec}s vt{vt}: {nv}/{nf} "
                                      f"negative (>10% — wholesale, not "
                                      f"return-credit)")
                                sys.exit(1)
                        n_fin = int(np.isfinite(panel[di]).sum())
                        k = f"{model}/{gcm}"
                        c = neg_census.setdefault(
                            k, {"n_negative_values": 0,
                                "min_km3": 0.0,
                                "negative_annual_km3_worst_decade": 0.0})
                        c["n_negative_values"] += n_neg
                        c["min_km3"] = min(c["min_km3"],
                                           float(np.nanmin(panel[di])))
                        neg_ann = float(np.nansum(
                            np.where(vt12 < 0, vt12, 0.0)))
                        c["negative_annual_km3_worst_decade"] = min(
                            c["negative_annual_km3_worst_decade"], neg_ann)
                    monthly_sum = np.nansum(panel[di, 0:12], axis=0)
                    # tolerance scales with the SUM OF ABSOLUTE months, not
                    # the net: in return-credit cells (a*use) positive and
                    # negative months cancel, so a net-relative tolerance
                    # explodes near zero while the arithmetic is fine. For
                    # all-positive variables sum|months| equals the monthly
                    # sum (not the stored vt12), so this bound is EQUIVALENT
                    # up to that scale basis — marginally more permissive
                    # only when vt12 < sum(months) by more than the
                    # tolerance itself (lane-review L2).
                    month_abs = np.nansum(np.abs(panel[di, 0:12]), axis=0)
                    bad = (np.abs(monthly_sum[finite] - vt12[finite])
                           > 1e-4 * month_abs[finite] + 1e-9)
                    if bool(bad.any()):
                        print(f"VALIDATE FAIL: {var} {model}/{gcm} {scen} "
                              f"{dec}s: vt12 != sum(vt0-11) "
                              f"({int(bad.sum())} cells)")
                        sys.exit(1)
            first = vv[mi, 0, 0, 12]
            first = first.filled(np.nan) if hasattr(first, "filled") \
                else np.asarray(first)
            print(f"    {model}/{gcm}: valid cells {int(np.isfinite(first).sum())}, "
                  f"global {np.nansum(first):,.1f} km3/yr "
                  f"(2010s, {SCENARIOS[0]})")
    finally:
        ds.close()
    if neg_census:
        print(f"  negative-value census ({var} — return-credit cells, "
              f"allowed for actuals): "
              + "; ".join(f"{k}: n={v['n_negative_values']:,}, worst-decade "
                          f"negative sum {v['negative_annual_km3_worst_decade']:.1f} km3/yr"
                          for k, v in sorted(neg_census.items())))
    print(f"  validation PASSED: {var}")
    return neg_census


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--out-dir", type=Path,
                    default=BASE / "data" / "processed" / "ws2_role1")
    ap.add_argument("--variables", type=str, default=None,
                    help="comma-separated subset; default = all in manifest")
    ap.add_argument("--scenarios", type=str, default="ssp126,ssp370,ssp585",
                    help="comma-separated scenario scope (default: the full "
                         "family). The rev-6 balance lane acquires ssp370 "
                         "only (baseline-only, ACTUALUSE_20260903) — the "
                         "ensemble assertion and output dims follow this "
                         "scope; the scope is stamped into provenance")
    ap.add_argument("--skip-hash", action="store_true")
    args = ap.parse_args()

    global SCENARIOS
    SCENARIOS = tuple(s.strip() for s in args.scenarios.split(",") if s.strip())
    if not SCENARIOS:
        print("REFUSED: empty --scenarios")
        sys.exit(1)

    rows = list(csv.DictReader(open(args.manifest)))
    off_scope = sorted({r["scenario"] for r in rows} - set(SCENARIOS))
    if off_scope:
        print(f"REFUSED: manifest contains scenarios outside --scenarios: "
              f"{off_scope}")
        sys.exit(1)
    results_by_dest = {r["dest"]: r for r in csv.DictReader(open(args.results))}
    unresolved = [r for r in csv.DictReader(open(args.results))
                  if r["status"] not in ("downloaded", "verified-existing")]
    if unresolved:
        print(f"REFUSED: {len(unresolved)} manifest rows unresolved in the "
              f"download results (first: {unresolved[0]['dest']})")
        sys.exit(1)

    variables = (args.variables.split(",") if args.variables
                 else sorted({r["variable"] for r in rows}))
    area = cell_area_m2(np.arange(89.75, -90, -0.5), 720)

    for var in variables:
        vrows = [r for r in rows if r["variable"] == var]
        if not vrows:
            print(f"REFUSED: no manifest rows for '{var}'")
            sys.exit(1)
        print(f"=== {var} ({len(vrows)} raw files) ===", flush=True)
        audit_path = args.out_dir / f"{var}_raw_audit.csv"
        args.out_dir.mkdir(parents=True, exist_ok=True)
        calendars = audit_raw(var, vrows, results_by_dest,
                              args.skip_hash, audit_path)
        members = assert_ensemble(vrows)

        data = np.full((len(members), len(SCENARIOS), len(DECADES), N_VT,
                        360, 720), np.nan, dtype=np.float32)
        soc_by_model, cal_by_member = {}, {}
        for row in vrows:
            mi = members.index((row["model"], row["gcm"]))
            si = SCENARIOS.index(row["scenario"])
            dest = BASE / row["dest"]
            cal = calendars[row["dest"]]
            data[mi, si] = process_member_file(dest, var, area, cal)
            soc_by_model[row["model"]] = row["soc"]
            cal_by_member[f"{row['model']}/{row['gcm']}"] = cal
            print(f"    processed {row['model']}/{row['gcm']}/{row['scenario']}",
                  flush=True)

        out_path = args.out_dir / f"demand_{var}_2015soc.nc"
        provenance = {
            "title": f"WRI Role-1 demand product: {var}, fixed-2015soc",
            "product_contract": (
                "WS2 Role 1 / Gate G10-E (decided 2026-08-28): 13 value "
                "types (monthly climatology + annual mean per decade), NO "
                "demand quantiles. NOT a water-index family file "
                "(GUARDRAILS §6/§7 govern waterIndexUnderlyingData_* only); "
                "ws2_role3 pelec precedent."),
            "monthly_means_assumption": (
                "Raw files carry no time bounds and no cell_methods; monthly "
                "values are ASSUMED monthly-mean fluxes (GUARDRAILS §18.3) "
                "and converted with each file's own calendar month lengths."),
            "members_resolved": [f"{m}/{g}" for m, g in members],
            "soc_by_model": soc_by_model,
            "calendar_by_member": cal_by_member,
            "scenario_scope": list(SCENARIOS),
            "source_manifest": str(args.manifest),
            "source_results": str(args.results),
            "manifest_sha256": sha256_file(args.manifest),
            "results_sha256": sha256_file(args.results),
            "raw_audit_csv": str(audit_path),
            "raw_audit_sha256": sha256_file(audit_path),
            # Review 2026-09-02 H1: a commit hash identifies the code only
            # if the code is tracked and clean at that commit — record the
            # script's own sha256 unconditionally, and say when the commit
            # does not contain it.
            "processing_code_commit": git_commit(),
            "processor_file_sha256": sha256_file(Path(__file__)),
            # tracked (ls-files) AND unmodified vs HEAD — `git diff --quiet`
            # alone exits 0 for an untracked file, the exact H1 failure mode
            "processor_tracked_clean": str(
                subprocess.call(
                    ["git", "-C", str(BASE), "ls-files", "--error-unmatch",
                     str(Path(__file__).relative_to(BASE))],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
                and subprocess.call(
                    ["git", "-C", str(BASE), "diff", "--quiet", "HEAD", "--",
                     str(Path(__file__).relative_to(BASE))],
                    stderr=subprocess.DEVNULL) == 0),
            "created_utc": datetime.now(timezone.utc).isoformat(
                timespec="seconds"),
        }
        write_output(out_path, var, members, data, area, provenance)
        neg_census = validate_output(out_path, var, members)
        if neg_census:
            # attributes-only amendment (the established in-place pattern):
            # the return-credit census becomes part of the product's record
            with Dataset(str(out_path), "a") as dsa:
                dsa.negative_value_census = json.dumps(neg_census)
                dsa.negative_value_policy = (
                    "actual consumption (a*use) is a model residual "
                    "(abstraction - return flows); locally negative values "
                    "are return-credit cells, retained as-is (contract "
                    "extended deliberately 2026-09-03; wholesale negativity "
                    ">10% refuses)")
        print(f"=== {var} complete: {out_path} ===", flush=True)

    print("\nALL VARIABLES COMPLETE")


if __name__ == "__main__":
    main()
