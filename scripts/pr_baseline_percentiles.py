"""Stage 0: derive each land cell's 2020s wet-day precipitation percentiles.

The relative rungs (R95pD / R99pD -- days exceeding a cell's OWN historical extreme) cannot
be counted in the same streaming pass that reads the data, because the threshold is not known
until the baseline has been seen. So the 2020s chunks are read once here, reduced to a
per-cell distribution, and the thresholds written out for stage 1 to count against. Only the
two baseline chunks per member are moved (~147 GB of the ~911 GB total).

WHY THE PERCENTILE IS POOLED ACROSS MEMBERS, NOT PER MEMBER (user decision 2026-08-17).
Arik's sizing argument was framed on the pooled sample -- 14 GCMs x ~10 wet days a year is
~1,200 wet days a decade -- and pooling is right for two further reasons. A PER-MEMBER
threshold would pin every member at exactly 5% exceedance across the baseline decade by
construction, collapsing the inter-member spread that the confidence band is supposed to
carry. And the product's decadal statistics already pool (year x member) into one sample, so
a pooled threshold is the consistent choice. One threshold per cell, from every member and
every scenario in the 2020s window.

WHY WET DAYS AND NOT ALL DAYS. In an arid cell 95% of days are dry, so the 95th percentile of
ALL days is ~0 mm and "exceeding it" degenerates into "it rained at all". ETCCDI defines its
precipitation percentile indices on wet days (>= 1 mm) for exactly this reason, and so does
this script. WET_DAY_MM is the definition, not a tuning knob.

WHY A HISTOGRAM AND NOT THE VALUES. Pooled over 42 member-scenarios x 10 years, a cell holds
~153,000 daily values of which ~15,000 are wet; retaining them for 65,797 cells is ~2 TB.
A log-spaced histogram over [1, 1000] mm/day gives ~2.7% relative bin width -- far finer than
the decision needs, since a threshold of 17.5 vs 17.9 mm/day changes nothing -- for 67 MB.
Percentiles are read back by linear interpolation within the containing bin.

WHY A MINIMUM WET-DAY RULE. A p99 estimated from a handful of wet days is noise dressed as a
threshold. Arik set the bar deliberately low (2026-08-17): "we can be pretty liberal in our
openness to low-rainfall-frequency data, and accept as few as (say) 2 days per year per model
on average". Over 420 model-years that is MIN_POOLED_WET_DAYS pooled wet days. Cells below it
are masked FOR THE RELATIVE RUNGS ONLY -- the absolute rungs and the intensity metrics are
unaffected and stay global.

Run:
    .venv/bin/python3 scripts/pr_baseline_percentiles.py --plan
    .venv/bin/python3 scripts/pr_baseline_percentiles.py --run [--workers 4]
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BASE = "https://files.isimip.org/ISIMIP3b"
GROUP_I = "InputData/climate/atmosphere/bias-adjusted/global/daily"
EXTENDED = "SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily"

#: 14 GCMs publish daily `pr` -- TWO MORE than publish tasmax/tasmin. CESM2-WACCM and
#: IITM-ESM have no tasmax, so this layer does NOT share a member set with the heat ladder;
#: that is a declared, accepted difference (user decision 2026-08-17: prefer 14 over 12).
GROUP_I_GCMS = {
    "GFDL-ESM4": ("gfdl-esm4", "r1i1p1f1"),
    "IPSL-CM6A-LR": ("ipsl-cm6a-lr", "r1i1p1f1"),
    "MPI-ESM1-2-HR": ("mpi-esm1-2-hr", "r1i1p1f1"),
    "MRI-ESM2-0": ("mri-esm2-0", "r1i1p1f1"),
    "UKESM1-0-LL": ("ukesm1-0-ll", "r1i1p1f2"),
}
EXTENDED_GCMS = {
    "CESM2-WACCM": ("cesm2-waccm", "r1i1p1f1"),
    "CNRM-CM6-1": ("cnrm-cm6-1", "r1i1p1f2"),
    "CNRM-ESM2-1": ("cnrm-esm2-1", "r1i1p1f2"),
    "CanESM5": ("canesm5", "r1i1p1f1"),
    "EC-Earth3": ("ec-earth3", "r1i1p1f1"),
    "IITM-ESM": ("iitm-esm", "r1i1p1f1"),
    "KACE-1-0-G": ("kace-1-0-g", "r1i1p1f1"),
    "MIROC6": ("miroc6", "r1i1p1f1"),
    "TaiESM1": ("taiesm1", "r1i1p1f1"),
}
ALL_GCMS = sorted(GROUP_I_GCMS) + sorted(EXTENDED_GCMS)
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

BASELINE = (2020, 2029)
WET_DAY_MM = 1.0
PERCENTILES = (95.0, 99.0)
#: 2 wet days per model-year x 14 GCMs x 3 scenarios x 10 years = 840.
MIN_POOLED_WET_DAYS = 2 * 14 * 3 * 10

#: Log-spaced edges over [1, 1000] mm/day. 256 bins -> ~2.7% relative width.
N_BINS = 256
EDGES = np.geomspace(WET_DAY_MM, 1000.0, N_BINS + 1)

INTERIM = Path("data/interim/prthresh")
STAGING = Path("data/raw/_pr_staging")
CKPT = INTERIM / "_ckpt_baseline"
PROV = INTERIM / "download_provenance.csv"
MASK = Path("data/masks/ISIMIP3b_landseamask_no-ant.nc")
OUT = INTERIM / "baseline_percentiles.nc"

SEC_PER_DAY = 86400.0


def member_root(gcm):
    if gcm in GROUP_I_GCMS:
        return GROUP_I, GROUP_I_GCMS[gcm], "group-I"
    return EXTENDED, EXTENDED_GCMS[gcm], "extended"


#: Seconds a single socket read may block. NOT the download budget -- urllib's `timeout` is
#: PER SOCKET OPERATION, so passing 3600 here does not mean "an hour to fetch a 2 GB file",
#: it means "wait an hour for one read() to return". Measured 2026-08-17: after ~6 h of
#: sustained concurrent transfers the file server left connections half-open, and eight
#: workers sat at 0.0% CPU with byte counts frozen, each due to block for the full hour.
#: A stalled read must fail in minutes so the chunk can be retried on a fresh connection.
READ_TIMEOUT_S = 120
#: Per-chunk download attempts before the member is failed to the driver, which retries it.
CHUNK_ATTEMPTS = 4


def http(url, timeout=READ_TIMEOUT_S):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "curl/8.4"}), timeout=timeout)


def download_with_retry(url, dest, tag, stem, expect_size=None):
    """Fetch to `dest`, restarting the whole chunk on a stall or truncation.

    Deliberately restarts rather than resuming with a Range header: a stalled transfer that
    silently resumes at the wrong offset produces a file that passes a size check and fails
    sha512, and diagnosing that costs more than re-fetching. The publisher sidecar's size is
    checked here so a truncated body is caught immediately rather than at the checksum.
    """
    last = None
    for attempt in range(CHUNK_ATTEMPTS):
        t0 = time.time()
        try:
            with http(url, READ_TIMEOUT_S) as r, open(dest, "wb") as fh:
                while True:
                    b = r.read(1 << 22)
                    if not b:
                        break
                    fh.write(b)
            got = dest.stat().st_size
            if expect_size and got != expect_size:
                raise OSError(f"truncated: {got} of {expect_size} bytes")
            mbps = got / 1e6 / max(1e-9, time.time() - t0)
            print(f"    [{tag}] {stem}.nc {got / 1e9:.2f} GB @ {mbps:.1f} MB/s", flush=True)
            return
        except Exception as e:  # noqa: BLE001
            last = e
            dest.unlink(missing_ok=True)
            wait = 30 * (attempt + 1)
            print(f"    [{tag}] {stem} attempt {attempt + 1}/{CHUNK_ATTEMPTS} failed "
                  f"({type(e).__name__}: {str(e)[:60]}) -- retrying in {wait}s", flush=True)
            if attempt < CHUNK_ATTEMPTS - 1:
                time.sleep(wait)
    raise OSError(f"{stem}: {CHUNK_ATTEMPTS} download attempts failed; last: {last}")


def land_mask():
    with Dataset(MASK) as ds:
        raw = ds.variables["mask"][:]
        a = np.asarray(raw.filled(np.nan) if np.ma.isMaskedArray(raw) else raw,
                       "f8").squeeze()
        lat = np.asarray(ds.variables["lat"][:])
        lon = np.asarray(ds.variables["lon"][:])
    return (np.isfinite(a) & (a > 0.5)), lat, lon


def list_chunks(gcm, scenario):
    """Discover spans from the server (GUARDRAILS 3). An empty listing is a FAILURE."""
    root, _, _ = member_root(gcm)
    url = f"{BASE}/{root}/{scenario}/{gcm}/"
    for attempt in range(4):
        try:
            body = http(url, 120).read().decode("utf-8", "replace")
            if len(body) < 200:
                raise OSError("short body -- rate limited, NOT an empty directory")
            pat = re.compile(r"_pr_global_daily_(\d{4})_(\d{4})\.nc$")
            spans = sorted({f"{m.group(1)}_{m.group(2)}"
                            for n in set(re.findall(r'href="([^"]+\.nc)"', body))
                            if (m := pat.search(n))})
            if not spans:
                raise OSError(f"no pr files matched at {url}")
            return spans
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"    listing retry {attempt + 1} ({str(e)[:60]})", flush=True)
            time.sleep(5 * (attempt + 1))
    return []


def baseline_spans(gcm, scenario):
    """Only the chunks overlapping the 2020s window."""
    out = []
    for s in list_chunks(gcm, scenario):
        y0, y1 = (int(x) for x in s.split("_"))
        if y1 >= BASELINE[0] and y0 <= BASELINE[1]:
            out.append(s)
    return out


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
            "gcm": gcm, "scenario": scenario, "stage": "baseline", "span": span,
            "url": url, "size": dest.stat().st_size, "sha512": digest,
            "sidecar_checksum": "match" if side.get("checksum") else "NO SIDECAR",
            "retained": "no -- deleted after reduction"})
    return dest


def histogram_member(args):
    """One (GCM, scenario): return (hist, wet_count, meta). Raw is deleted."""
    gcm, scenario = args
    tag = f"{gcm}/{scenario}"
    mask, _, _ = land_mask()
    n_land = int(mask.sum())
    hist = np.zeros((n_land, N_BINS), np.int32)
    wet = np.zeros(n_land, np.int64)
    meta = {}
    for span in baseline_spans(gcm, scenario):
        ck = CKPT / f"{gcm}_{scenario}_{span}.npz"
        if ck.exists():
            z = np.load(ck)
            hist += z["hist"]
            wet += z["wet"]
            print(f"    [{tag}] {span} from checkpoint", flush=True)
            continue
        path = stream_chunk(gcm, scenario, span, tag)
        try:
            ds = Dataset(path)
            v = ds.variables["pr"]
            t = ds.variables["time"]
            cal = getattr(t, "calendar", None) or "standard"
            years = np.array([d.year for d in num2date(t[:], t.units, calendar=cal)])
            meta.setdefault("calendar", set()).add(cal)
            meta.setdefault("units", set()).add(getattr(v, "units", "UNDECLARED"))
            h = np.zeros((n_land, N_BINS), np.int32)
            w = np.zeros(n_land, np.int64)
            sel = np.flatnonzero((years >= BASELINE[0]) & (years <= BASELINE[1]))
            for i in range(0, sel.size, 120):
                idx = sel[i:i + 120]
                a = np.asarray(v[idx[0]:idx[-1] + 1], "float32")[
                    idx - idx[0]][:, mask] * SEC_PER_DAY
                wetmask = a >= WET_DAY_MM
                w += wetmask.sum(axis=0)
                # Flatten (cell, bin) to ONE index and use bincount. np.add.at is a
                # documented slow path in numpy (it takes an unbuffered ufunc route), and
                # this loop runs 84 times over ~3,650 days each; bincount does the same
                # scatter-add in compiled code.
                if wetmask.any():
                    cells = np.broadcast_to(np.arange(n_land, dtype=np.int64),
                                            a.shape)[wetmask]
                    bins = np.clip(np.digitize(a[wetmask], EDGES) - 1, 0, N_BINS - 1)
                    h += np.bincount(cells * N_BINS + bins.astype(np.int64),
                                     minlength=n_land * N_BINS
                                     ).reshape(n_land, N_BINS).astype(np.int32)
                del a, wetmask
            ds.close()
            CKPT.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(ck, hist=h, wet=w)
            hist += h
            wet += w
        finally:
            path.unlink(missing_ok=True)
    return hist, wet, {k: sorted(v) for k, v in meta.items()}


def percentile_from_hist(hist, q):
    """Linear interpolation inside the containing log bin. hist is (cell, bin)."""
    tot = hist.sum(axis=1)
    out = np.full(hist.shape[0], np.nan, np.float32)
    ok = tot > 0
    if not ok.any():
        return out
    cum = np.cumsum(hist[ok], axis=1).astype(np.float64)
    target = tot[ok] * (q / 100.0)
    idx = np.array([np.searchsorted(c, tg) for c, tg in zip(cum, target)])
    idx = np.clip(idx, 0, N_BINS - 1)
    lo = EDGES[idx]
    hi = EDGES[idx + 1]
    below = np.where(idx > 0, cum[np.arange(cum.shape[0]), np.maximum(idx - 1, 0)], 0.0)
    below = np.where(idx > 0, below, 0.0)
    inbin = hist[ok][np.arange(idx.size), idx].astype(np.float64)
    frac = np.where(inbin > 0, (target - below) / np.maximum(inbin, 1e-9), 0.0)
    out[ok] = (lo + np.clip(frac, 0.0, 1.0) * (hi - lo)).astype(np.float32)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    mask, lat, lon = land_mask()
    jobs = [(g, s) for g in ALL_GCMS for s in SCENARIOS]
    print(f"stage         0 -- 2020s wet-day percentiles for the relative rungs")
    print(f"members       {len(ALL_GCMS)} GCMs x {len(SCENARIOS)} scenarios = {len(jobs)}")
    print(f"land cells    {int(mask.sum()):,} (landseamask_no-ant.nc, Antarctica EXCLUDED)")
    print(f"baseline      {BASELINE[0]}-{BASELINE[1]}, wet day >= {WET_DAY_MM:g} mm")
    print(f"percentiles   {PERCENTILES} of WET days, POOLED over all members and scenarios")
    print(f"min wet days  {MIN_POOLED_WET_DAYS} pooled (= 2 per model-year x "
          f"{len(ALL_GCMS)} GCMs x {len(SCENARIOS)} scen x 10 yr)")
    print(f"histogram     {N_BINS} log bins over [{WET_DAY_MM:g}, 1000] mm/day "
          f"(~{100 * (EDGES[1] / EDGES[0] - 1):.1f}% relative width)")
    print(f"bytes moved   ~{len(jobs) * 3.5:.0f} GB (2 baseline chunks/member), NOT retained")
    if not a.run:
        print("\nPLANNING ONLY. Re-run with --run.")
        return 0

    INTERIM.mkdir(parents=True, exist_ok=True)
    n_land = int(mask.sum())
    hist = np.zeros((n_land, N_BINS), np.int64)
    wet = np.zeros(n_land, np.int64)
    meta_all, done, failed = {}, 0, []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(histogram_member, j): j for j in jobs}
        for f in as_completed(futs):
            g, s = futs[f]
            try:
                h, w, m = f.result()
                hist += h
                wet += w
                for k, v in m.items():
                    meta_all.setdefault(k, set()).update(v)
                done += 1
                print(f"  [{done}/{len(jobs)}] {g} {s}: "
                      f"{w.sum():,} wet cell-days", flush=True)
            except Exception as e:  # noqa: BLE001
                failed.append((g, s, str(e)[:160]))
                print(f"  FAILED {g} {s}: {str(e)[:160]}", flush=True)

    if failed:
        print(f"\n{len(failed)} member(s) failed -- NOT writing thresholds. Re-run to resume.")
        for g, s, e in failed:
            print(f"  {g} {s}: {e}")
        return 1

    thresholds = {q: percentile_from_hist(hist, q) for q in PERCENTILES}
    usable = wet >= MIN_POOLED_WET_DAYS
    print(f"\nPOOLED wet days per cell: min {wet.min():,} median {int(np.median(wet)):,} "
          f"max {wet.max():,}")
    print(f"cells meeting the >= {MIN_POOLED_WET_DAYS} rule: {int(usable.sum()):,} of "
          f"{n_land:,} ({usable.mean():.2%}) -- {int((~usable).sum()):,} masked for the "
          f"RELATIVE rungs only")
    for q in PERCENTILES:
        v = thresholds[q][usable]
        print(f"  p{q:g} threshold (mm/day): min {np.nanmin(v):.1f} "
              f"median {np.nanmedian(v):.1f} max {np.nanmax(v):.1f}")

    LAT, LON = mask.shape
    idx = np.flatnonzero(mask.ravel())

    def scatter(vec, dtype="f4"):
        g = np.full(LAT * LON, np.nan, "f4")
        g[idx] = vec
        return g.reshape(LAT, LON).astype(dtype)

    with Dataset(OUT, "w") as ds:
        ds.createDimension("lat", LAT)
        ds.createDimension("lon", LON)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for q in PERCENTILES:
            name = f"p{q:g}"
            var = ds.createVariable(name, "f4", ("lat", "lon"), zlib=True, complevel=4,
                                    fill_value=np.float32(np.nan))
            var[:] = scatter(np.where(usable, thresholds[q], np.nan))
            var.units = "mm day-1"
            var.long_name = (f"{q:g}th percentile of WET-day (>= {WET_DAY_MM:g} mm) daily "
                             f"precipitation, {BASELINE[0]}-{BASELINE[1]}")
        wv = ds.createVariable("pooled_wet_days", "i4", ("lat", "lon"), zlib=True,
                               complevel=4, fill_value=-1)
        wv[:] = np.nan_to_num(scatter(wet.astype("f4")), nan=-1).astype("i4")
        wv.long_name = "pooled wet-day count in the baseline window across all members"
        uv = ds.createVariable("usable", "i1", ("lat", "lon"), zlib=True, complevel=4)
        uv[:] = np.nan_to_num(scatter(usable.astype("f4")), nan=0).astype("i1")
        uv.long_name = (f"1 where pooled wet days >= {MIN_POOLED_WET_DAYS}; the relative "
                        "rungs are masked elsewhere. The ABSOLUTE rungs are not.")
        ds.baseline_window = f"{BASELINE[0]}-{BASELINE[1]}"
        ds.wet_day_definition = f">= {WET_DAY_MM:g} mm/day"
        ds.pooling = ("ONE threshold per cell from every member and scenario pooled. A "
                      "per-member threshold would pin each member at exactly (100-q)% "
                      "exceedance in the baseline and collapse the inter-member spread the "
                      "confidence band carries.")
        ds.min_wet_days_rule = (
            f"{MIN_POOLED_WET_DAYS} pooled wet days = 2 per model-year across "
            f"{len(ALL_GCMS)} GCMs x {len(SCENARIOS)} scenarios x 10 years (user decision "
            f"2026-08-17). Deliberately liberal.")
        ds.histogram = (f"{N_BINS} log-spaced bins over [{WET_DAY_MM:g}, 1000] mm/day; "
                        "percentiles by linear interpolation within the containing bin.")
        ds.gcms = ",".join(ALL_GCMS)
        for k, v in meta_all.items():
            setattr(ds, f"source_{k}", " | ".join(sorted(v)))
        ds.created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"\nwrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
