"""B1 precipitation pipeline: streamed daily ISIMIP3b forcing -> monthly-mean flux.

WS1 Stage B (waterRiskIndex_beta/v2/PLAN_ws1_ssp_reingestion.md, Gate B approved
2026-08-21): downloads the daily bias-adjusted `pr` InputData chunk by chunk with
the hardened fetch (fresh connection, stall guard, length+sha256+NetCDF checks),
aggregates each chunk to MONTHLY-MEAN FLUX (kg m-2 s-1 — the measured legacy
value-type contract: vt0-11 are monthly means of the rate, so mean-of-days is the
correct reduction), writes the monthly intermediate in the model-output layout
(`{interim}/w5e5/{gcm}_{scenario}/pr_monthly_*.nc`, variable `pr`, monthly time
axis) so `process_water_variable.py` computes the 20 value types with the SAME
pooling code as every other variable, records a per-chunk receipt (url, raw
sha256/bytes, output sha256, code commit), and only then deletes the raw chunk
(retention exception approved at Gate B, delete-after-receipt).

Resumable: a chunk whose receipt exists and whose intermediate verifies is
skipped. Burst pacing between chunk downloads.

Usage:
    # offline unit test of the aggregation + engine conversion arithmetic (B0):
    python scripts/download_and_aggregate_precip.py --self-test

    # build the chunk manifest from the A1 inventory:
    python scripts/download_and_aggregate_precip.py --build-manifest A1_files.csv --manifest OUT.csv

    # run the pipeline:
    python scripts/download_and_aggregate_precip.py --manifest OUT.csv [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from download_water_manifest import fetch  # hardened: fresh conn, stall guard, .part

BASE_DIR = SCRIPTS.parent
RAW_DIR = BASE_DIR / "data" / "raw" / "water_precip_daily_staging"
INTERIM_DIR = BASE_DIR / "data" / "interim" / "water_precip_monthly"
RECEIPTS = INTERIM_DIR / "chunk_receipts.jsonl"
INPUTDATA_BASE = "https://files.isimip.org/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily"
BURST, COOLDOWN = 4, 90  # chunks are ~1.9 GB, so bursts are smaller than for model output


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def code_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(BASE_DIR), "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def aggregate_daily_to_monthly(ds: xr.Dataset, var: str = "pr") -> xr.DataArray:
    """Daily flux -> monthly-MEAN flux (kg m-2 s-1), calendar-aware.

    Mirrors the measured legacy contract: monthly value types are monthly means
    of the rate. `resample(time="MS").mean()` averages exactly the days of each
    calendar month under any CF calendar decoded by xarray/cftime.
    """
    return ds[var].resample(time="MS").mean()


def self_test() -> None:
    """B0: synthetic-field verification of aggregation + engine conversion."""
    import pandas as pd

    sec_per_year = 60 * 60 * 24 * 365.25
    time_idx = xr.cftime_range("2015-01-01", "2016-12-31", freq="D", calendar="proleptic_gregorian")
    lat, lon = [0.25, 0.75], [10.25]

    # Case 1: constant field c — every monthly mean must equal c exactly
    c = 3.0e-5
    ds = xr.Dataset(
        {"pr": (("time", "lat", "lon"), np.full((len(time_idx), 2, 1), c))},
        coords={"time": time_idx, "lat": lat, "lon": lon},
    )
    monthly = aggregate_daily_to_monthly(ds)
    assert monthly.shape[0] == 24, f"expected 24 months, got {monthly.shape[0]}"
    np.testing.assert_allclose(monthly.values, c, rtol=1e-12)

    # Engine conversion arithmetic (v2 wri_engine): raw * sec_per_year, monthlies / 12
    mm_per_month = c * sec_per_year / 12
    np.testing.assert_allclose(float(monthly.isel(time=0, lat=0, lon=0)) * sec_per_year / 12,
                               mm_per_month, rtol=1e-12)

    # Case 2: month-constant steps — each monthly mean must recover its month's value
    # exactly (a mean over identical days), proving no cross-month leakage.
    month_vals = {m: 1e-5 * m for m in range(1, 13)}
    vals = np.array([month_vals[t.month] for t in time_idx])[:, None, None] * np.ones((1, 2, 1))
    ds2 = ds.copy(); ds2["pr"].values = vals
    monthly2 = aggregate_daily_to_monthly(ds2)
    got = {pd.Timestamp(str(t)).month: float(monthly2.sel(time=t).isel(lat=0, lon=0))
           for t in monthly2.time.values[:12]}
    for m in range(1, 13):
        np.testing.assert_allclose(got[m], month_vals[m], rtol=1e-12,
                                   err_msg=f"month {m} mean wrong")

    # Case 3: NaN cell stays NaN, finite cell unaffected
    ds3 = ds.copy(); ds3["pr"].values[:, 1, 0] = np.nan
    monthly3 = aggregate_daily_to_monthly(ds3)
    assert np.isnan(monthly3.isel(lat=1, lon=0)).all()
    np.testing.assert_allclose(monthly3.isel(lat=0, lon=0).values, c, rtol=1e-12)

    print("B0 self-test PASSED: constant-field exactness, per-month recovery, "
          "NaN propagation, engine mm/month arithmetic")


def build_manifest(inventory_csv: Path, out_csv: Path) -> None:
    rows = []
    with open(inventory_csv) as f:
        for r in csv.DictReader(f):
            if r["tree"] == "inputdata_daily" and r["scenario"] in ("ssp126", "ssp370", "ssp585"):
                rows.append({
                    "gcm": r["dir_gcm"], "scenario": r["scenario"], "file": r["file"],
                    "url": f"{INPUTDATA_BASE}/{r['scenario']}/{r['dir_gcm']}/{r['file']}",
                    "listed_size": r["listed_size"],
                })
    rows.sort(key=lambda r: (r["gcm"], r["scenario"], r["file"]))
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"manifest: {len(rows)} daily chunks -> {out_csv}")


def already_done(row: dict, receipts_seen: set) -> bool:
    key = row["file"]
    if key not in receipts_seen:
        return False
    out_path = INTERIM_DIR / "w5e5" / f"{row['gcm'].lower()}_{row['scenario']}" / monthly_name(row["file"])
    if not out_path.exists():
        return False
    try:
        xr.open_dataset(out_path).close()
        return True
    except Exception:
        return False


def monthly_name(daily_fname: str) -> str:
    return daily_fname.replace("_daily_", "_monthly-from-daily_")


def run_pipeline(manifest_csv: Path, limit: int | None) -> None:
    rows = list(csv.DictReader(open(manifest_csv)))
    if limit:
        rows = rows[:limit]
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    receipts_seen = set()
    if RECEIPTS.exists():
        for line in open(RECEIPTS):
            try:
                receipts_seen.add(json.loads(line)["file"])
            except Exception:
                pass

    commit = code_commit()
    failed, done_dl = 0, 0
    t0 = time.time()
    for i, row in enumerate(rows, 1):
        if already_done(row, receipts_seen):
            print(f"[{i}/{len(rows)}] receipt+intermediate verified: {row['file']}", flush=True)
            continue
        raw = RAW_DIR / row["file"]
        out_dir = INTERIM_DIR / "w5e5" / f"{row['gcm'].lower()}_{row['scenario']}"
        out_path = out_dir / monthly_name(row["file"])
        try:
            nbytes, digest = fetch(row["url"], raw)
            ds = xr.open_dataset(raw)
            monthly = aggregate_daily_to_monthly(ds)
            n_months = int(monthly.shape[0])
            out_dir.mkdir(parents=True, exist_ok=True)
            monthly.to_dataset(name="pr").to_netcdf(
                out_path, encoding={"pr": {"dtype": "float32", "zlib": True, "complevel": 4}})
            ds.close()
            with open(RECEIPTS, "a") as f:
                f.write(json.dumps({
                    "file": row["file"], "url": row["url"], "raw_bytes": nbytes,
                    "raw_sha256": digest, "n_months": n_months,
                    "monthly_out": str(out_path.relative_to(BASE_DIR)),
                    "monthly_sha256": sha256_file(out_path), "code_commit": commit,
                    "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                }) + "\n")
            raw.unlink()  # delete-after-receipt (Gate B retention policy)
            done_dl += 1
            print(f"[{i}/{len(rows)}] aggregated {row['gcm']}/{row['scenario']} "
                  f"{row['file'].split('_')[-2]}-{row['file'].split('_')[-1][:4]} "
                  f"({nbytes/2**20:.0f} MB raw -> {n_months} months, t+{(time.time()-t0)/60:.1f}m)",
                  flush=True)
            if done_dl % BURST == 0:
                print(f"--- burst of {BURST} chunks done, cooling {COOLDOWN}s ---", flush=True)
                time.sleep(COOLDOWN)
        except Exception as e:
            failed += 1
            raw.unlink(missing_ok=True)
            print(f"[{i}/{len(rows)}] FAILED {row['file']}: {e}", flush=True)
    print(f"\n=== precip pipeline: {len(rows)-failed}/{len(rows)} chunks resolved, {failed} FAILED ===")
    sys.exit(1 if failed else 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--build-manifest", type=Path)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    if args.self_test:
        self_test()
    elif args.build_manifest:
        build_manifest(args.build_manifest, args.manifest)
    elif args.manifest:
        run_pipeline(args.manifest, args.limit)
    else:
        ap.error("one of --self-test / --build-manifest / --manifest required")


if __name__ == "__main__":
    main()
