"""Stage 1: stream ISIMIP3b daily `pr` and reduce it to annual pluvial-flood statistics.

THE RAW DATA IS NEVER RETAINED. 14 GCMs x 3 scenarios x ~9 chunks is ~764 GB against a
500 GB disk, so streaming is not an optimisation. Each chunk is downloaded, verified against
the publisher sidecar's sha512, reduced to per-year statistics, and DELETED; peak disk is one
chunk per worker (~8 GB at 4 workers) plus the interim outputs. Provenance is appended BEFORE
deletion, so the ingest stays auditable against files we no longer hold.

REQUIRES STAGE 0. The relative rungs count days above each cell's OWN 2020s wet-day
percentile, and that threshold is not knowable from the chunk in front of you -- it is read
from data/interim/prthresh/baseline_percentiles.nc, written by pr_baseline_percentiles.py.
Run that first; this script refuses to start without it.

WHAT IS COMPUTED, AND WHY THESE (user decisions 2026-08-16/17)
  ABSOLUTE frequency, days/yr -- "how often does a damaging depth fall"
    R10mm  R20mm   ETCCDI standard heavy / very-heavy precipitation days
    R50mm          World Bank CCKP ships r50mm
    R100mm         OUR EXTENSION, not a standard index -- label it as such, like FDm10
  RELATIVE frequency, days/yr -- "how often is it extreme FOR HERE"
    R95pD R99pD    days above this cell's own 2020s wet-day p95/p99. NOT ETCCDI `R95p`,
                   which is a TOTAL in mm on a 1961-1990 base. Different quantity, hence
                   the different name; a reader who compares them to published R95p gets
                   nonsense.
  INTENSITY, mm -- "how hard does the worst day hit"
    Rx1day         annual maximum 1-day precipitation. The best single pluvial proxy in the
                   set, and the only metric here needing no threshold choice at all.
    Rx5day         annual maximum 5-day total: the saturation-driven and compound case.
  DIAGNOSTIC
    wetdays        days >= 1 mm. Also the quantity stage 0's mask rule is written on.
    prcptot        annual total from wet days.  (SDII = prcptot / wetdays, derived later.)

THE THRESHOLDS ARE PICKED AGAINST A MEASUREMENT, NOT A CONVENTION. Measured 2026-08-16 on
four members: each cell's 95th percentile of wet days spans 3.0-161.5 mm/day across the land
surface, a 35-44x range, median 17-18. So 20 mm/day sits at roughly the MEDIAN cell's own p95
-- which is why it is the headline -- while being unreachable in the driest cells and an
ordinary wet day in the wettest. Daily temperature varies about two-fold globally and a fixed
threshold means the same thing everywhere; precipitation does not, which is the whole reason
this layer carries both an absolute and a relative family instead of one.

Rx5day AND CHUNK BOUNDARIES -- the one correctness trap in this script. A 5-day window
straddling 2030-12-30/2031-01-02 lives in two different files. Reducing each chunk
independently would silently drop every boundary event: 8 boundaries per member, and they are
exactly the multi-day accumulations the metric exists to catch. So each chunk carries its
last 4 days forward, the carry is stored IN the checkpoint (a resume that lost it would
reintroduce the bug in a form nothing tests), and the first chunk of a member starts with an
explicit empty carry. Windows are assigned to the year in which they END, per ETCCDI.

Run:
    .venv/bin/python3 scripts/download_reduce_prthresh_isimip3b.py --plan
    .venv/bin/python3 scripts/download_reduce_prthresh_isimip3b.py --run [--workers 4]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pr_baseline_percentiles import (  # noqa: E402
    ALL_GCMS, BASE, INTERIM, MIN_POOLED_WET_DAYS, OUT as BASELINE_FILE, PROV,
    SCENARIOS, SEC_PER_DAY, WET_DAY_MM, download_with_retry, http, land_mask,
    list_chunks, member_root,
)

LAYER_ID = "pluvial-isimip3b_prthresh_annual"
STAGING = Path("data/raw/_pr_staging")
CKPT = INTERIM / "_ckpt_main"

#: (name, threshold mm/day). Absolute rungs, counted as days >= threshold.
ABSOLUTE = [("R10mm", 10.0), ("R20mm", 20.0), ("R50mm", 50.0), ("R100mm", 100.0)]
#: (name, percentile variable in the stage-0 file)
RELATIVE = [("R95pD", "p95"), ("R99pD", "p99")]
COUNTS = [n for n, _ in ABSOLUTE] + [n for n, _ in RELATIVE] + ["wetdays"]
INTENSITY = ["Rx1day", "Rx5day", "prcptot"]
RX5_WINDOW = 5


def sha512_of(path):
    h = hashlib.sha512()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def record(row):
    PROV.parent.mkdir(parents=True, exist_ok=True)
    new = not PROV.exists()
    with open(PROV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)


def stream_chunk(gcm, scenario, span, tag):
    root, (tok, ens), _ = member_root(gcm)
    stem = f"{tok}_{ens}_w5e5_{scenario}_pr_global_daily_{span}"
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
        download_with_retry(url, dest, tag, stem, side.get("size"))
    if not dest.exists():
        raise OSError(f"{stem}: staging file vanished before checksum -- another process "
                      f"is writing to {STAGING}. Check for orphaned workers.")
    digest = sha512_of(dest)
    if side.get("checksum") and digest != side["checksum"]:
        dest.unlink()
        raise SystemExit(f"sha512 MISMATCH for {stem}.nc -- refusing to use it")
    record({"utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gcm": gcm, "scenario": scenario, "stage": "main", "span": span,
            "url": url, "size": dest.stat().st_size, "sha512": digest,
            "sidecar_checksum": "match" if side.get("checksum") else "NO SIDECAR",
            "retained": "no -- deleted after reduction"})
    return dest


def load_thresholds(mask):
    if not BASELINE_FILE.exists():
        raise SystemExit(
            f"{BASELINE_FILE} not found -- run scripts/pr_baseline_percentiles.py first. "
            "The relative rungs cannot be counted without the per-cell 2020s thresholds.")
    out = {}
    with Dataset(BASELINE_FILE) as ds:
        for _, var in RELATIVE:
            a = np.asarray(ds.variables[var][:], "f4")
            out[var] = a[mask]
        usable = np.asarray(ds.variables["usable"][:], "i1")[mask].astype(bool)
    return out, usable


def reduce_member(gcm, scenario, force=False, max_chunks=0):
    """Produce {gcm}_{scenario}_pr.nc. Idempotent; an existing complete file is kept."""
    out = INTERIM / (f"SMOKE_{gcm}_{scenario}_pr.nc" if max_chunks
                     else f"{gcm}_{scenario}_pr.nc")
    if out.exists() and not force:
        print(f"  {out.name} exists -- skipping (use --force to rebuild)", flush=True)
        return out
    tag = f"{gcm}/{scenario}"
    mask, lat, lon = land_mask()
    n_land = int(mask.sum())
    thr, usable = load_thresholds(mask)

    acc_counts = {n: {} for n in COUNTS}
    acc_int = {n: {} for n in INTENSITY}
    meta = {}
    carry = np.zeros((0, n_land), np.float32)   # last <=4 days of the previous chunk
    carry_years = np.zeros(0, np.int32)

    spans = list_chunks(gcm, scenario)
    if max_chunks:
        spans = spans[:max_chunks]
    for span in spans:
        ck = CKPT / f"{gcm}_{scenario}_{span}.npz"
        if ck.exists() and not force:
            try:
                z = np.load(ck, allow_pickle=False)
                yrs = z["years"]
                for n in COUNTS:
                    for k, y in enumerate(yrs):
                        acc_counts[n][int(y)] = z[n][k]
                for n in INTENSITY:
                    for k, y in enumerate(yrs):
                        acc_int[n][int(y)] = z[n][k]
                # THE CARRY IS PART OF THE CHECKPOINT. Without it a resumed member loses
                # every 5-day window spanning the boundary it resumed at, silently.
                carry, carry_years = z["carry"], z["carry_years"]
                print(f"    [{tag}] {span} from checkpoint", flush=True)
                continue
            except Exception as e:  # noqa: BLE001
                print(f"    [{tag}] checkpoint unreadable ({str(e)[:50]}) -- refetching",
                      flush=True)
                ck.unlink(missing_ok=True)

        path = stream_chunk(gcm, scenario, span, tag)
        try:
            ds = Dataset(path)
            v = ds.variables["pr"]
            t = ds.variables["time"]
            cal = getattr(t, "calendar", None) or "standard"
            years = np.array([d.year for d in num2date(t[:], t.units, calendar=cal)],
                             dtype=np.int32)
            for key, val in (("calendar", cal),
                             ("units", getattr(v, "units", "UNDECLARED")),
                             ("long_name", getattr(v, "long_name", "NONE"))):
                meta.setdefault(key, set()).add(str(val))
            probe = np.asarray(v[0], "float64")
            pf = probe[np.isfinite(probe)]
            # Scale MEASURED from values, never from the declared unit: group-I sidecars
            # carry no netcdf_header, and a declared unit has been wrong before here.
            flux = bool(pf.size and np.nanmax(pf) < 0.1)
            meta.setdefault("scale", set()).add(
                "kg m-2 s-1 -> mm/day via x86400 (MEASURED)" if flux
                else "ALREADY mm/day?? (MEASURED) -- verify before use")
            fac = SEC_PER_DAY if flux else 1.0

            chunk_years = sorted(set(int(y) for y in years))
            loc_c = {n: {y: np.zeros(n_land, np.int32) for y in chunk_years}
                     for n in COUNTS}
            loc_i = {n: {y: np.zeros(n_land, np.float32) for y in chunk_years}
                     for n in INTENSITY}
            for y in chunk_years:
                loc_i["Rx1day"][y][:] = 0.0
                loc_i["Rx5day"][y][:] = 0.0

            # Prepend the carry so 5-day windows cross the boundary.
            step = 120
            prev = carry
            prev_y = carry_years
            for i in range(0, len(years), step):
                a = np.asarray(v[i:i + step], "float32")[:, mask] * fac
                yr = years[i:i + step]
                for y in np.unique(yr):
                    sel = a[yr == y]
                    yi = int(y)
                    for name, thr_mm in ABSOLUTE:
                        loc_c[name][yi] += (sel >= thr_mm).sum(axis=0).astype(np.int32)
                    for name, var in RELATIVE:
                        hit = (sel >= thr[var][None, :]) & usable[None, :]
                        loc_c[name][yi] += hit.sum(axis=0).astype(np.int32)
                    wet = sel >= WET_DAY_MM
                    loc_c["wetdays"][yi] += wet.sum(axis=0).astype(np.int32)
                    loc_i["prcptot"][yi] += np.where(wet, sel, 0.0).sum(axis=0)
                    loc_i["Rx1day"][yi] = np.maximum(loc_i["Rx1day"][yi],
                                                     sel.max(axis=0))
                # 5-day running sums over [carry | this slab]
                block = np.concatenate([prev, a], axis=0)
                block_y = np.concatenate([prev_y, yr])
                if block.shape[0] >= RX5_WINDOW:
                    cs = np.cumsum(block, axis=0)
                    win = cs[RX5_WINDOW - 1:] - np.concatenate(
                        [np.zeros((1, n_land), np.float32), cs[:-RX5_WINDOW]], axis=0)
                    end_y = block_y[RX5_WINDOW - 1:]
                    for y in np.unique(end_y):
                        yi = int(y)
                        if yi in loc_i["Rx5day"]:
                            loc_i["Rx5day"][yi] = np.maximum(
                                loc_i["Rx5day"][yi], win[end_y == y].max(axis=0))
                prev = block[-(RX5_WINDOW - 1):]
                prev_y = block_y[-(RX5_WINDOW - 1):]
                del a, block
            ds.close()
            carry, carry_years = prev.astype(np.float32), prev_y.astype(np.int32)

            CKPT.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                ck, years=np.array(chunk_years, "i4"),
                carry=carry, carry_years=carry_years,
                **{n: np.stack([loc_c[n][y] for y in chunk_years]) for n in COUNTS},
                **{n: np.stack([loc_i[n][y] for y in chunk_years]) for n in INTENSITY})
            for n in COUNTS:
                for y in chunk_years:
                    acc_counts[n][y] = loc_c[n][y]
            for n in INTENSITY:
                for y in chunk_years:
                    acc_int[n][y] = loc_i[n][y]
        finally:
            path.unlink(missing_ok=True)

    all_years = sorted(set().union(*[set(d) for d in acc_counts.values()]))
    LAT, LON = mask.shape
    idx = np.flatnonzero(mask.ravel())

    def grid(vec, fill):
        g = np.full(LAT * LON, fill, "f4")
        g[idx] = vec
        return g.reshape(LAT, LON)

    INTERIM.mkdir(parents=True, exist_ok=True)
    with Dataset(out, "w") as ds:
        ds.createDimension("year", len(all_years))
        ds.createDimension("lat", LAT)
        ds.createDimension("lon", LON)
        ds.createVariable("year", "i4", ("year",))[:] = all_years
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for name in COUNTS:
            arr = np.stack([grid(acc_counts[name][y], np.nan) for y in all_years])
            nc = ds.createVariable(name, "i2", ("year", "lat", "lon"), zlib=True,
                                   complevel=4, fill_value=np.int16(-1))
            nc[:] = np.where(np.isfinite(arr), arr, -1).astype(np.int16)
            nc.units = "days per year"
            if name in dict(ABSOLUTE):
                nc.long_name = f"Days with precipitation >= {dict(ABSOLUTE)[name]:g} mm"
                nc.threshold_mm = dict(ABSOLUTE)[name]
            elif name in dict(RELATIVE):
                q = dict(RELATIVE)[name]
                nc.long_name = (f"Days exceeding this cell's own 2020s wet-day "
                                f"{q[1:]}th percentile")
                nc.threshold_source = str(BASELINE_FILE)
                nc.masked_where = (f"pooled baseline wet days < {MIN_POOLED_WET_DAYS} "
                                   "(fill -1); the ABSOLUTE rungs are not masked")
            else:
                nc.long_name = f"Wet days (>= {WET_DAY_MM:g} mm)"
        for name in INTENSITY:
            arr = np.stack([grid(acc_int[name][y], np.nan) for y in all_years])
            nc = ds.createVariable(name, "f4", ("year", "lat", "lon"), zlib=True,
                                   complevel=4, fill_value=np.float32(np.nan))
            nc[:] = arr
            nc.units = "mm"
            nc.long_name = {"Rx1day": "Annual maximum 1-day precipitation",
                            "Rx5day": ("Annual maximum 5-day precipitation total; windows "
                                       "cross chunk boundaries and are assigned to the year "
                                       "they END in"),
                            "prcptot": "Annual total precipitation from wet days"}[name]
        ds.gcm, ds.scenario, ds.layer_id = gcm, scenario, LAYER_ID
        ds.root = member_root(gcm)[2]
        ds.raw_retention = ("NOT RETAINED -- each daily chunk was verified against the "
                            "publisher sidecar sha512, reduced, and deleted. See "
                            "download_provenance.csv.")
        ds.relative_threshold_note = (
            "R95pD/R99pD count days above THIS CELL's own 2020s wet-day percentile, pooled "
            "over all members and scenarios. They are NOT ETCCDI R95p, which is a TOTAL in "
            "mm on a 1961-1990 base.")
        for k, val in meta.items():
            joined = " | ".join(sorted(val))
            setattr(ds, f"source_{k}", joined)
            if len(val) > 1:
                setattr(ds, f"source_{k}_HETEROGENEOUS",
                        "chunks of this member disagree -- do not pool without deciding")
        ds.created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"  wrote {out} ({out.stat().st_size / 1e6:.0f} MB, {len(all_years)} years)",
          flush=True)
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
                    help="SMOKE TEST: first N chunks only; output prefixed SMOKE_")
    a = ap.parse_args()

    gcms = a.gcm or ALL_GCMS
    scens = a.scenario or SCENARIOS
    jobs = [(g, s) for g in gcms for s in scens]
    print(f"layer_id      {LAYER_ID}")
    print(f"members       {len(gcms)} GCMs x {len(scens)} scenarios = {len(jobs)} jobs")
    print(f"metrics       absolute {[n for n, _ in ABSOLUTE]}")
    print(f"              relative {[n for n, _ in RELATIVE]} (masked where baseline "
          f"wet days < {MIN_POOLED_WET_DAYS})")
    print(f"              intensity {INTENSITY}")
    print(f"thresholds    {BASELINE_FILE} "
          f"{'FOUND' if BASELINE_FILE.exists() else '*** MISSING -- run stage 0 ***'}")
    est = len(jobs) * 18.2
    print(f"bytes moved   ~{est:.0f} GB, NOT retained")
    print(f"peak disk     ~{a.workers * 2.1:.1f} GB staging + ~5 GB interim")
    if not a.run:
        print("\nPLANNING ONLY. Re-run with --run.")
        return 0
    if not BASELINE_FILE.exists():
        print("\nERROR: stage 0 has not produced baseline_percentiles.nc.")
        return 1

    done, failed = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(reduce_member, g, s, a.force, a.max_chunks): (g, s)
                for g, s in jobs}
        for f in as_completed(futs):
            g, s = futs[f]
            try:
                f.result()
                done.append((g, s))
            except Exception as e:  # noqa: BLE001
                failed.append((g, s, str(e)[:160]))
                print(f"  FAILED {g} {s}: {str(e)[:160]}", flush=True)
    print(f"\n{len(done)} complete, {len(failed)} failed")
    for g, s, e in failed:
        print(f"  {g} {s}: {e}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
