"""Download the Lange 2020 `let` ensemble for the TCFD tropical-cyclone layer.

`let` = fraction of each cell's LAND AREA exposed to tropical cyclones in a given year
(Lange 2020, ISIMIP2b DerivedOutputData, impact model ke-tg-meanfield -- an Emanuel
kinetic-energy wind model driven by each GCM's storm field).

Ensemble: 1 impact model x 4 CMIP5 GCMs = 4 members per scenario x {rcp26, rcp60}
= 8 files, ~58 MB. Composition is IDENTICAL across the two scenarios, so the shared
2020s baseline is valid (OUTPUT-SPEC "Shared 2020s baseline").

Coverage facts, enumerated 2026-08-11 rather than assumed:

  * There is NO rcp85 for this family -- rcp26/rcp60 only, the same gap `lew` has.
  * There is NO ISIMIP3b/SSP re-issue of TC exposure. The Lange 2020 exposure concept
    WAS carried into 3b, but split across Heinicke2026 (driedarea, floodedarea) and
    Zantout2025 (heatwave, wildfire, cropfailure); ke-tg-meanfield appears in no 3a or
    3b publication directory. ISIMIP3b's newer TC product is a different thing entirely
    -- per-storm wind footprints under InputData/climate/tropical_cyclones/MIT, noted
    for a future delivery.
  * soc/sens are uniformly `nosoc` / `co2`, so no harmonization compromise is needed
    (contrast the burntarea ensemble's mixed soc).

`historical` (1861-2005) and `pre-industrial` are NOT fetched: the future files start in
2006, which already covers the 2010s and the 2020s baseline decade.

The repository API (data.isimip.org) is behind Anubis anti-bot and the CLI downloader
needs `httpx`, which is absent from this venv -- so URLs are constructed deterministically
and fetched with the stdlib.

Unlike the ISIMIP3b `biomes` files, ISIMIP2b publishes a `.json` sidecar carrying `size`
AND `sha512` for every file, so downloads here are verified by CHECKSUM, not merely by
byte count. Downloads are resumable via Range.

Usage:
    python scripts/download_let_cyclone.py [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP2b/DerivedOutputData/Lange2020/KE-TG-meanfield"
LAYER_ID = "cyclone_let_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "let"
MODEL = "ke-tg-meanfield"
GCMS = ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
SCENARIOS = ["rcp26", "rcp60"]
YEARS = "2006_2099"

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP2b let ingest)"
CHUNK = 1 << 20


def build_items():
    """Deterministic URL per (gcm, scenario). Every file spans 2006_2099."""
    items = []
    for gcm in GCMS:
        for scen in SCENARIOS:
            stem = (f"lange2020_{MODEL}_{gcm}_ewembi_{scen}_nosoc_co2_"
                    f"{VAR}_global_annual_{YEARS}")
            items.append(dict(
                fname=f"{stem}.nc4",
                url=f"{BASE}/{gcm}/future/{stem}.nc4",
                sidecar=f"{BASE}/{gcm}/future/{stem}.json",
                gcm=gcm, scenario=scen, member=f"{MODEL}_{gcm}",
            ))
    return items


def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def sha512_of(path):
    h = hashlib.sha512()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(CHUNK), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(url, dest, expected_size):
    """Download `url` to `dest`, resuming a partial file. Returns bytes written."""
    have = dest.stat().st_size if dest.exists() else 0
    if have == expected_size:
        return 0
    if have > expected_size:          # stale/corrupt -- start over
        dest.unlink()
        have = 0

    headers = {"User-Agent": USER_AGENT}
    mode = "wb"
    if have:
        headers["Range"] = f"bytes={have}-"
        mode = "ab"

    req = urllib.request.Request(url, headers=headers)
    written = 0
    with urllib.request.urlopen(req, timeout=600) as r, open(dest, mode) as fh:
        while True:
            buf = r.read(CHUNK)
            if not buf:
                break
            fh.write(buf)
            written += len(buf)
    return written


def one(item):
    dest = OUT_DIR / item["fname"]
    try:
        meta = get_json(item["sidecar"])
    except urllib.error.HTTPError as e:
        return dict(item, ok=False, msg=f"sidecar HTTP {e.code}")
    except Exception as e:  # noqa: BLE001 - report, never abort the batch
        return dict(item, ok=False, msg=f"sidecar {type(e).__name__}: {e}")

    size = int(meta["size"])
    want = meta.get("checksum", "")
    if meta.get("checksum_type", "sha512") != "sha512":
        return dict(item, ok=False, msg=f"unexpected checksum_type {meta.get('checksum_type')}")

    try:
        wrote = fetch(item["url"], dest, size)
    except Exception as e:  # noqa: BLE001
        return dict(item, ok=False, msg=f"GET {type(e).__name__}: {e}")

    got = dest.stat().st_size if dest.exists() else 0
    if got != size:
        return dict(item, ok=False, msg=f"size {got} != {size}")

    digest = sha512_of(dest)
    if want and digest != want:
        return dict(item, ok=False, msg="sha512 MISMATCH")

    return dict(item, ok=True, bytes=size, sha512=digest,
                wrote=wrote, skipped=(wrote == 0), msg="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    items = build_items()
    print(f"{len(items)} files = {len(GCMS)} members x {len(SCENARIOS)} scenarios")
    if args.dry_run:
        for it in items:
            print(it["url"])
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Serial: 8 small files, and parallel requests to files.isimip.org get rate-limited
    # into empty/parital responses (see the isimip-search-download skill).
    results = []
    for i, it in enumerate(items, 1):
        r = one(it)
        results.append(r)
        tag = "skip" if r.get("skipped") else ("ok" if r["ok"] else "FAIL")
        print(f"[{i}/{len(items)}] {tag:4s} {r['fname']}"
              + ("" if r["ok"] else f"  <-- {r['msg']}"), flush=True)

    good = [r for r in results if r["ok"]]
    bad = [r for r in results if not r["ok"]]
    total = sum(r["bytes"] for r in good)
    print(f"\n{len(good)}/{len(items)} verified by sha512, {total/2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "member", "gcm", "scenario"])
            w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                w.writerow({k: r[k] for k in w.fieldnames})
        print(f"provenance -> {OUT_DIR / 'download_provenance.csv'}")

    if bad:
        print(f"\n{len(bad)} FAILED (re-run to resume):")
        for r in bad:
            print(f"  {r['fname']}: {r['msg']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
