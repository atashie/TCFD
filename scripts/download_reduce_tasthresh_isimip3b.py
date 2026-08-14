"""Stage 1: stream ISIMIP3b daily tasmax/tasmin and reduce to annual threshold counts.

THE RAW DATA IS NEVER RETAINED. 12 GCMs x 3 scenarios x 2 variables x 9 chunks is ~1.34 TB
and this machine has ~600 GB of usable temporary space, so retention is not an option and
streaming is not an optimisation. Each chunk is downloaded, verified against the publisher
sidecar's sha512, reduced to per-year counts, and DELETED. Peak disk is one chunk per
worker plus the counts.

The counts are a SUFFICIENT STATISTIC for every rung of the ladder, which is why the whole
ladder is computed in one pass: once the day is read, testing it against nine thresholds
instead of one costs nothing, and a rung we do not publish today needs no second 1.34 TB.

WHAT IS RECORDED INSTEAD OF THE FILE. Provenance (url, size, sha512, UTC timestamp) is
appended to download_provenance.csv BEFORE the file is deleted, so the ingest remains
auditable against a file we no longer hold. That is a deliberate deviation from the
data/raw/{layer_id}/ retention convention and is declared in the interim file's attributes.

WHY A DAY COUNT NEEDS THE CALENDAR MEASURED (GUARDRAILS 9). A 360_day member gets 360
chances a year to cross a threshold and a proleptic_gregorian one gets 365.25. That is a
~1.5% low bias which is CONSTANT PER MEMBER, so it does not average out -- it moves with
ensemble composition and lands in the trend. This script records days-per-year per member
so the pooling decision rests on a measurement. It also detects Kelvin from the VALUES
rather than the declared unit, because the five group-I GCMs publish sidecars with no
netcdf_header at all (verified 2026-08-14) and their declared unit is unreadable without
opening the file.

Run:
    .venv/bin/python3 scripts/download_reduce_tasthresh_isimip3b.py --plan
    .venv/bin/python3 scripts/download_reduce_tasthresh_isimip3b.py --run [--workers 4]
    .venv/bin/python3 scripts/download_reduce_tasthresh_isimip3b.py --run --gcm GFDL-ESM4
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_tas_thresholds_nature import (  # noqa: E402
    BASE, COLD_MAX, COLD_MIN, EXTENDED, EXTENDED_GCMS, GROUP_I, GROUP_I_GCMS,
    HOT, TROPICAL_NIGHT,
)

LAYER_ID = "heatdays-isimip3b_tasthresh_annual"
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
INTERIM = Path("data/interim/tasthresh")
STAGING = Path("data/raw/_tasthresh_staging")
CKPT = INTERIM / "_ckpt"
PROV = INTERIM / "download_provenance.csv"

#: rung -> (variable, comparison, threshold degC). Order is the output variable order.
LADDER: list[tuple[str, str, str, float]] = (
    [(f"hd{t:g}", "tasmax", "gt", t) for t in HOT]
    + [("ID", "tasmax", "lt", COLD_MAX[0])]
    + [("FD", "tasmin", "lt", COLD_MIN[0]),
       ("FDm10", "tasmin", "lt", COLD_MIN[1])]
    + [(f"TR{t:g}", "tasmin", "gt", t) for t in TROPICAL_NIGHT]
)


def member_root(gcm: str):
    if gcm in GROUP_I_GCMS:
        return GROUP_I, GROUP_I_GCMS[gcm], "group-I"
    return EXTENDED, EXTENDED_GCMS[gcm], "extended"


ALL_GCMS = sorted(GROUP_I_GCMS) + sorted(EXTENDED_GCMS)


def http(url: str, timeout=120):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "curl/8.4"}), timeout=timeout)


def list_chunks(gcm: str, scenario: str, var: str) -> list[str]:
    """Discover chunks from the server. GUARDRAILS 3 -- never hardcode a span list."""
    root, (tok, ens), _ = member_root(gcm)
    url = f"{BASE}/{root}/{scenario}/{gcm}/"
    for attempt in range(4):
        try:
            body = http(url, 120).read().decode("utf-8", "replace")
            if len(body) < 200:
                raise OSError("short body -- rate limited, NOT an empty directory")
            names = sorted(set(re.findall(r'href="([^"]+\.nc)"', body)))
            pat = re.compile(rf"_{re.escape(var)}_global_daily_(\d{{4}})_(\d{{4}})\.nc$")
            spans = sorted({m.group(0)[-12:-3] for n in names if (m := pat.search(n))})
            if not spans:
                raise OSError(f"no {var} files matched at {url} -- wrong depth or variable")
            return spans
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"    listing retry {attempt + 1} ({str(e)[:60]})")
            time.sleep(5 * (attempt + 1))
    return []


def sha512_of(path: Path) -> str:
    h = hashlib.sha512()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def record(row: dict):
    PROV.parent.mkdir(parents=True, exist_ok=True)
    new = not PROV.exists()
    with open(PROV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)


def stream_chunk(gcm: str, scenario: str, var: str, span: str, tag: str) -> Path:
    """Download + verify one chunk into staging. Returns the path; caller deletes it."""
    root, (tok, ens), _ = member_root(gcm)
    stem = f"{tok}_{ens}_w5e5_{scenario}_{var}_global_daily_{span}"
    url = f"{BASE}/{root}/{scenario}/{gcm}/{stem}.nc"
    STAGING.mkdir(parents=True, exist_ok=True)
    dest = STAGING / f"{stem}.nc"

    side = {}
    try:
        side = json.loads(http(f"{BASE}/{root}/{scenario}/{gcm}/{stem}.json", 60)
                          .read().decode())
    except Exception:  # noqa: BLE001
        pass

    if not (dest.exists() and side.get("size") and dest.stat().st_size == side["size"]):
        t0 = time.time()
        with http(url, 3600) as r, open(dest, "wb") as fh:
            while True:
                b = r.read(1 << 22)
                if not b:
                    break
                fh.write(b)
        mbps = dest.stat().st_size / 1e6 / max(1e-9, time.time() - t0)
        print(f"    [{tag}] {stem}.nc {dest.stat().st_size / 1e9:.2f} GB @ {mbps:.1f} MB/s")

    digest = sha512_of(dest)
    if side.get("checksum") and digest != side["checksum"]:
        dest.unlink()
        raise SystemExit(f"sha512 MISMATCH for {stem}.nc -- refusing to use it")
    record({"utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gcm": gcm, "scenario": scenario, "variable": var, "span": span,
            "url": url, "size": dest.stat().st_size, "sha512": digest,
            "sidecar_checksum": "match" if side.get("checksum") else "NO SIDECAR",
            "retained": "no -- deleted after reduction"})
    return dest


def reduce_member(gcm: str, scenario: str, force=False, max_chunks: int = 0) -> Path:
    """Produce {gcm}_{scenario}_counts.nc. Idempotent: an existing complete file is kept.

    max_chunks > 0 truncates the span list -- a SMOKE TEST ONLY. The resulting file is
    written with partial_smoke_test set, so it can never be mistaken for a real member.
    """
    out = INTERIM / (f"SMOKE_{gcm}_{scenario}_counts.nc" if max_chunks
                     else f"{gcm}_{scenario}_counts.nc")
    if out.exists() and not force:
        print(f"  {out.name} exists -- skipping (use --force to rebuild)")
        return out
    INTERIM.mkdir(parents=True, exist_ok=True)
    tag = f"{gcm}/{scenario}"

    acc: dict[str, dict[int, np.ndarray]] = {r[0]: {} for r in LADDER}
    meta: dict[str, object] = {}
    lat = lon = None

    for var in ("tasmax", "tasmin"):
        rungs = [r for r in LADDER if r[1] == var]
        spans = list_chunks(gcm, scenario, var)
        if max_chunks:
            spans = spans[:max_chunks]
        for span in spans:
            # CHUNK-LEVEL RESUME. Member-level resume alone is too coarse for a ~41 h job:
            # a member that dies on its ninth chunk would re-download the eight before it,
            # ~30 GB and two hours. The first run WAS killed (process-group teardown at
            # 2026-08-14 13:45), which is why this exists. Each chunk's per-year counts are
            # a few MB compressed, so checkpointing them is far cheaper than re-fetching.
            ck = CKPT / f"{gcm}_{scenario}_{var}_{span}.npz"
            if ck.exists() and not force:
                try:
                    z = np.load(ck, allow_pickle=False)
                    for name, _, _, _ in rungs:
                        for k, y in enumerate(z["years"]):
                            acc[name][int(y)] = z[name][k]
                    if lat is None:
                        lat, lon = z["lat"], z["lon"]
                    print(f"    [{tag}] {var} {span} from checkpoint")
                    continue
                except Exception as e:  # noqa: BLE001
                    print(f"    [{tag}] checkpoint {ck.name} unreadable ({str(e)[:50]}) "
                          "-- refetching")
                    ck.unlink(missing_ok=True)

            path = stream_chunk(gcm, scenario, var, span, tag)
            try:
                ds = Dataset(path)
                v = ds.variables[var]
                t = ds.variables["time"]
                cal = getattr(t, "calendar", None) or "standard"
                tu = getattr(t, "units")
                years = np.array([d.year for d in num2date(t[:], tu, calendar=cal)])
                if lat is None:
                    lat = np.asarray(ds.variables["lat"][:])
                    lon = np.asarray(ds.variables["lon"][:])
                # Record EVERY distinct value seen across chunks, not the first. A unit or
                # calendar that changes mid-member is a documented ISIMIP failure mode
                # (burntarea classic is fraction-scaled for one GCM and percent for its
                # sibling), and a setdefault would silently keep only chunk 1's answer.
                for key, val in ((f"{var}_units", getattr(v, "units", "UNDECLARED")),
                                 (f"{var}_calendar", cal),
                                 (f"{var}_days_per_year", f"{len(years) / len(set(years)):.2f}")):
                    meta.setdefault(key, set()).add(str(val))

                probe = np.asarray(v[0], dtype="float64")
                pf = probe[np.isfinite(probe)]
                kelvin = bool(pf.size and pf.min() > 150.0)
                meta.setdefault(f"{var}_scale", set()).add(
                    "K" if kelvin else "degC (MEASURED from values, not declared)")
                off = 273.15 if kelvin else 0.0

                chunk_years = sorted(set(int(y) for y in years))
                local = {n: {y: None for y in chunk_years} for n, _, _, _ in rungs}
                step = 120
                for i in range(0, len(years), step):
                    a = np.asarray(v[i:i + step], dtype="float32") - off
                    yr = years[i:i + step]
                    for y in np.unique(yr):
                        sel = a[yr == y]
                        for name, _, cmp_, thr in rungs:
                            hit = (sel > thr) if cmp_ == "gt" else (sel < thr)
                            c = hit.sum(axis=0).astype("int16")
                            prev = local[name][int(y)]
                            local[name][int(y)] = c if prev is None else prev + c
                    del a
                ds.close()

                CKPT.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(
                    ck, years=np.array(chunk_years, "i4"), lat=lat, lon=lon,
                    **{n: np.stack([local[n][y] for y in chunk_years]) for n, _, _, _ in rungs})
                for name, _, _, _ in rungs:
                    for y in chunk_years:
                        acc[name][y] = local[name][y]
            finally:
                path.unlink(missing_ok=True)

    if not meta:
        meta = {"note": {"every chunk came from a checkpoint; declared units/calendar were "
                         "recorded on the run that produced them"}}
    all_years = sorted(set().union(*[set(d) for d in acc.values()]))
    with Dataset(out, "w") as ds:
        ds.createDimension("year", len(all_years))
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        ds.createVariable("year", "i4", ("year",))[:] = all_years
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for name, var, cmp_, thr in LADDER:
            arr = np.stack([acc[name][y] for y in all_years]).astype("int16")
            nc = ds.createVariable(name, "i2", ("year", "lat", "lon"),
                                   zlib=True, complevel=4)
            nc[:] = arr
            nc.units = "days per year"
            nc.long_name = (f"Annual count of days with {var} "
                            f"{'>' if cmp_ == 'gt' else '<'} {thr:g} degC")
            nc.source_variable = var
            nc.threshold_degC = thr
        ds.gcm, ds.scenario, ds.layer_id = gcm, scenario, LAYER_ID
        if max_chunks:
            ds.partial_smoke_test = (f"ONLY the first {max_chunks} chunk(s) were read. "
                                     "NOT a complete member. Do not process this file.")
        ds.root = member_root(gcm)[2]
        ds.raw_retention = ("NOT RETAINED -- each daily chunk was verified against the "
                            "publisher sidecar sha512, reduced, and deleted. See "
                            "download_provenance.csv.")
        ds.ladder = "; ".join(f"{n}={v}{'>' if c == 'gt' else '<'}{t:g}C"
                              for n, v, c, t in LADDER)
        for k, val in meta.items():
            joined = " | ".join(sorted(val)) if isinstance(val, set) else str(val)
            setattr(ds, k, joined)
            if isinstance(val, set) and len(val) > 1:
                setattr(ds, k + "_HETEROGENEOUS",
                        "chunks of this member disagree -- do not pool without deciding")
        ds.created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"  wrote {out} ({out.stat().st_size / 1e6:.0f} MB, {len(all_years)} years)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--gcm", action="append", default=None)
    ap.add_argument("--scenario", action="append", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--max-chunks", type=int, default=0,
                    help="SMOKE TEST: read only the first N chunks; output is prefixed SMOKE_")
    a = ap.parse_args()

    gcms = a.gcm or ALL_GCMS
    scens = a.scenario or SCENARIOS
    jobs = [(g, s) for g in gcms for s in scens]

    print(f"layer_id      {LAYER_ID}")
    print(f"members       {len(gcms)} GCMs "
          f"({sum(1 for g in gcms if g in GROUP_I_GCMS)} group-I + "
          f"{sum(1 for g in gcms if g in EXTENDED_GCMS)} extended)")
    print(f"scenarios     {scens}")
    print(f"jobs          {len(jobs)} (GCM x scenario), 2 variables x ~9 chunks each")
    print(f"ladder        {', '.join(r[0] for r in LADDER)}")
    est = len(jobs) * 2 * 9 * 2.07
    print(f"bytes moved   ~{est:.0f} GB  (~{est * 1000 / 9 / 3600:.1f} h at the measured "
          f"9 MB/s aggregate)")
    print(f"peak disk     ~{a.workers * 2.07:.1f} GB staging + "
          f"~{len(jobs) * 0.401:.1f} GB counts")
    print(f"raw retained  NO -- verified, reduced, deleted")

    if not a.run:
        print("\nPLANNING ONLY. Re-run with --run to execute.")
        for g, s in jobs:
            print(f"  {member_root(g)[2]:9s} {g:14s} {s}")
        return 0

    done, failed = [], []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(reduce_member, g, s, a.force, a.max_chunks): (g, s)
                for g, s in jobs}
        for f in as_completed(futs):
            g, s = futs[f]
            try:
                f.result()
                done.append((g, s))
            except Exception as e:  # noqa: BLE001
                failed.append((g, s, str(e)[:160]))
                print(f"  FAILED {g} {s}: {str(e)[:160]}")
    print(f"\n{len(done)} complete, {len(failed)} failed")
    for g, s, e in failed:
        print(f"  {g} {s}: {e}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
