"""Download the Heinicke2026 `driedarea` ensemble for the ISIMIP3b/SSP drought layer.

`driedarea` is the ISIMIP3b re-issue of the Lange 2020 drought-exposure concept, renamed by
hazard word rather than `le*` code. It is a DIFFERENT product from ISIMIP2b `led`, not a
newer version of the same files: different models, different GCM generation, different
scenarios.

Ensemble: 3 impact models x 5 CMIP6 GCMs x 3 SSPs = 45 files, ~168 MB.

Coverage ENUMERATED 2026-08-11 by listing all 15 `{MODEL}/{gcm}/future/` directories
serially (15 x HTTP 200, zero failures, zero gaps):

  * COMPLETE matrix -- every (model, gcm, scenario) combination exists. No holes to work
    around, so the shared 2020s baseline is valid on a uniform 15-member ensemble.
  * soc/sens are UNIFORM `2015soc` / `default` across all 45 files -- no harmonization
    compromise (contrast ISIMIP2b `led`, where soc is mixed per model, and the water_global
    sector, which offers SIX soc variants per scenario).
  * Bias adjustment `w5e5`; annual; 2015-2100.
  * `.json` sidecars WITH sha512 are published here, so downloads are verified by CHECKSUM
    (3b `biomes` largely lacks these -- do not assume from the round).

FILENAME GRAMMAR WARNING: unlike ISIMIP2b's Lange2020 files, Heinicke2026 filenames carry
NO leading publication token --

    h08_gfdl-esm4_w5e5_ssp126_2015soc_default_driedarea_global_annual_2015_2100.nc
    lange2020_lpjml_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099.nc4

so DerivedOutputData grammar is per-PUBLICATION, not per-product. Parse from the END.

DATA NATURE IS UNVERIFIED AT DOWNLOAD TIME and must be value-checked before processing
(GUARDRAILS 9). Do NOT assume it inherits `led`'s binary nature, and do NOT assume it
matches its own sibling: `floodedarea` -- same publication, same models -- is binary {0,1}
despite declaring units="1", and is non-NaN over 94.7% of the globe INCLUDING OCEAN, while
`driedarea` is reported land-masked. Check both nature and mask.

Usage:
    python scripts/download_driedarea_isimip3b.py [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/Heinicke2026"
LAYER_ID = "drought-isimip3b_driedarea_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "driedarea"
#: (server directory, filename model token)
MODELS = [("H08", "h08"), ("JULES-W2", "jules-w2"), ("WaterGAP2-2e", "watergap2-2e")]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
FORCING = "w5e5"
SOC = "2015soc"
SENS = "default"
YEARS = "2015_2100"

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b driedarea ingest)"
CHUNK = 1 << 20


def build_items():
    """One item per (model, gcm, scenario). The enumerated matrix is complete: no gaps."""
    items = []
    for mdir, model in MODELS:
        for gcm in GCMS:
            for scen in SCENARIOS:
                stem = (f"{model}_{gcm}_{FORCING}_{scen}_{SOC}_{SENS}_"
                        f"{VAR}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.nc",
                    url=f"{BASE}/{mdir}/{gcm}/future/{stem}.nc",
                    sidecar=f"{BASE}/{mdir}/{gcm}/future/{stem}.json",
                    model=model, gcm=gcm, scenario=scen, soc=SOC,
                    member=f"{model}_{gcm}",
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
        return dict(item, ok=False,
                    msg=f"unexpected checksum_type {meta.get('checksum_type')}")

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
    print(f"{len(items)} files = {len(MODELS)} models x {len(GCMS)} GCMs "
          f"x {len(SCENARIOS)} scenarios")
    if args.dry_run:
        for it in items:
            print(it["url"])
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Serial: parallel requests to files.isimip.org get rate-limited into empty/partial
    # responses (see the isimip-search-download skill).
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
                "fname", "url", "bytes", "sha512", "member", "model", "gcm",
                "scenario", "soc"])
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
