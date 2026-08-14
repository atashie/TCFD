"""Download the Zantout2025 `cropfailure` ensemble for the ISIMIP3b/SSP crop-failure layer.

`cropfailure` is the ISIMIP3b re-issue of the Lange 2020 crop-failure-exposure concept
(`lec`), renamed by hazard word rather than `le*` code. It is a DIFFERENT product from
ISIMIP2b `lec`, not a newer version of the same files: 8 crop models vs 2, CMIP6 vs CMIP5,
SSPs vs RCPs.

Ensemble: 8 impact models x 5 CMIP6 GCMs x 3 SSPs = 120 files, ~705 MB.

Coverage ENUMERATED 2026-08-13 by listing all 8 `{MODEL}/` and 40 `{MODEL}/{gcm}/future/`
directories serially (48 x HTTP 200, zero empty listings, zero gaps). The variable field
was projected at `$(NF-4)` rather than grepped for a believed-in token (GUARDRAILS 8):
160/160 files in `future/` are `cropfailure`, with no crop suffix and no PFT split.

  * COMPLETE matrix -- every (model, gcm, scenario) combination exists. 40 members per
    scenario, and the composition is IDENTICAL across ssp126/ssp370/ssp585, so the shared
    2020s baseline is valid on a uniform 40-member ensemble.
  * soc/sens are UNIFORM `2015soc` / `default` across all 160 files -- no harmonization
    compromise, no model dropped to make the ensemble uniform (contrast ISIMIP2b `led`,
    where soc is mixed per model and harmonizing would have cost 2 of 8 models).
  * Bias adjustment `w5e5`; annual; 2015-2100, so the 2020s baseline decade is complete.
  * GUARDRAILS 1 does not apply: `annual` is the ONLY cadence published here (all 160
    files are `..._global_annual_2015_2100.nc`), so there is no resolution to choose.

NO `.json` SIDECARS. Unlike Heinicke2026 `driedarea` -- same round, same product type --
Zantout2025 publishes no sidecars at all (`{stem}.nc.json` returns 404, verified
2026-08-13). Two consequences, both real downgrades to state plainly rather than paper
over:

  * There is no upstream sha512 to verify against. Integrity is checked against the
    `Content-Length` returned by a HEAD on the same URL, which catches truncation and
    resumption errors but NOT silent corruption in transit.
  * The sha512 recorded in `download_provenance.csv` is computed HERE, from the bytes on
    disk. It is a self-consistency receipt for later re-verification (did this file change
    since ingest?), NOT a confirmation that we received what the publisher intended.

FILENAME GRAMMAR WARNING: Zantout2025 filenames DO carry a leading publication token,
where its 3b sibling Heinicke2026 does not --

    zantout2025_crover_gfdl-esm4_w5e5_ssp126_2015soc_default_cropfailure_global_annual_2015_2100.nc
    h08_gfdl-esm4_w5e5_ssp126_2015soc_default_driedarea_global_annual_2015_2100.nc

so DerivedOutputData grammar is per-PUBLICATION, not per-round and not per-product. Parse
from the END regardless.

DATA NATURE IS UNVERIFIED AT DOWNLOAD TIME and must be value-checked before processing
(GUARDRAILS 9). Do NOT assume it inherits `lec`'s nature, and do NOT assume it matches a
sibling in its own publication: within Heinicke2026, `floodedarea` is binary {0,1} while
`driedarea` is not, and the two differ on masking as well. Run
`scripts/check_cropfailure_nature.py` after this script and before any processor.

FILE SIZES VARY 2.8x ACROSS MODELS (promet 3.13 MB -> isam 8.74 MB at the same GCM and
scenario), which points at heterogeneous land/cropland masks rather than a storage
difference. Measure per-member cell counts at processing time; it drives the
minimum-model publication mask.

Usage:
    python scripts/download_cropfailure_isimip3b.py [--dry-run]
"""

import argparse
import csv
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/Zantout2025"
LAYER_ID = "cropfailure-isimip3b_cropfailure_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "cropfailure"
PUBLICATION = "zantout2025"
#: (server directory, filename model token) -- both enumerated 2026-08-13, not derived by
#: lowercasing: the server dir is CYGMA1P74 while the filename token is cygma1p74, and
#: LPJML/lpjml likewise. Do not compute one from the other.
MODELS = [
    ("CROVER", "crover"),
    ("CYGMA1P74", "cygma1p74"),
    ("EPIC-IIASA", "epic-iiasa"),
    ("ISAM", "isam"),
    ("LDNDC", "ldndc"),
    ("LPJML", "lpjml"),
    ("PEPIC", "pepic"),
    ("PROMET", "promet"),
]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
FORCING = "w5e5"
SOC = "2015soc"
SENS = "default"
YEARS = "2015_2100"

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b cropfailure ingest)"
CHUNK = 1 << 20


def build_items():
    """One item per (model, gcm, scenario). The enumerated matrix is complete: no gaps."""
    items = []
    for mdir, model in MODELS:
        for gcm in GCMS:
            for scen in SCENARIOS:
                stem = (f"{PUBLICATION}_{model}_{gcm}_{FORCING}_{scen}_{SOC}_{SENS}_"
                        f"{VAR}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.nc",
                    url=f"{BASE}/{mdir}/{gcm}/future/{stem}.nc",
                    model=model, gcm=gcm, scenario=scen, soc=SOC,
                    member=f"{model}_{gcm}",
                ))
    return items


def head_size(url):
    """Expected byte count from a HEAD. This is the ONLY upstream integrity signal
    Zantout2025 offers -- there are no sidecars and therefore no published checksum."""
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        length = r.headers.get("Content-Length")
        if length is None:
            raise ValueError("no Content-Length on HEAD")
        return int(length)


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
        size = head_size(item["url"])
    except urllib.error.HTTPError as e:
        return dict(item, ok=False, msg=f"HEAD HTTP {e.code}")
    except Exception as e:  # noqa: BLE001 - report, never abort the batch
        return dict(item, ok=False, msg=f"HEAD {type(e).__name__}: {e}")

    try:
        wrote = fetch(item["url"], dest, size)
    except Exception as e:  # noqa: BLE001
        return dict(item, ok=False, msg=f"GET {type(e).__name__}: {e}")

    got = dest.stat().st_size if dest.exists() else 0
    if got != size:
        return dict(item, ok=False, msg=f"size {got} != {size}")

    return dict(item, ok=True, bytes=size, sha512=sha512_of(dest),
                wrote=wrote, skipped=(wrote == 0), msg="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    items = build_items()
    print(f"{len(items)} files = {len(MODELS)} models x {len(GCMS)} GCMs "
          f"x {len(SCENARIOS)} scenarios")
    print("integrity: Content-Length only -- Zantout2025 publishes no .json sidecars, "
          "so there is no upstream sha512 to verify against")
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
    print(f"\n{len(good)}/{len(items)} size-verified, {total/2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "sha512_source", "member", "model",
                "gcm", "scenario", "soc"])
            w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                row = {k: r[k] for k in w.fieldnames if k != "sha512_source"}
                # Say where the digest came from. Heinicke2026 rows in a sibling layer's
                # provenance mean "matched the publisher"; these mean "computed locally".
                row["sha512_source"] = "computed-locally (no upstream sidecar)"
                w.writerow(row)
        print(f"provenance -> {OUT_DIR / 'download_provenance.csv'}")

    if bad:
        print(f"\n{len(bad)} FAILED (re-run to resume):")
        for r in bad:
            print(f"  {r['fname']}: {r['msg']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
