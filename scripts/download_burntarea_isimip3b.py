"""Download the ISIMIP3b burntarea-total ensemble for the TCFD wildfire (SSP) layer.

Ensemble (user decision 2026-08-08 -- maximize members across BOTH GCMs and LSMs, and
prefer the recent SSP experiments over the RCP round):

    LSM                soc                   sens     cadence   GCMs
    classic            2015soc               default  monthly   gfdl-esm4, ukesm1-0-ll
    elm-eca            2015soc-from-histsoc  default  monthly   all 5
    lpjml5-7-10-fire   2015soc               default  monthly   all 5
    visit              2015soc               default  monthly   all 5
    mc2-usfs-r87g5c1   nat                   default  ANNUAL    all 5

  = 22 members per scenario x {ssp126, ssp370, ssp585} = 66 files, ~6.0 GB.
    Member composition is IDENTICAL across the three scenarios, so the shared 2020s
    baseline is valid (isimip-process-visualize: pool per-scenario only when composition
    differs).

Three selection decisions, each recorded rather than inherited:

  1. SECTOR: ISIMIP3b publishes burntarea under BOTH `fire` and `biomes`, and they are the
     SAME files -- 2001/2031 basenames shared, spot-checked pair identical on
     Content-Length AND ETag. We take `fire` for the four monthly models and `biomes` only
     for mc2-usfs (which appears nowhere else). Ingesting both would double-weight models.

  2. MIXED soc, deliberately. No single soc token spans the ensemble: elm-eca publishes
     ONLY `2015soc-from-histsoc`, visit ONLY `2015soc`, mc2-usfs ONLY `nat`. A uniform
     filter would drop 5-15 members. We keep all 22 and record `soc_by_model` in the
     processed output (GUARDRAILS 9).
     For `classic` -- which offers both -- we deliberately pick `2015soc` and AVOID
     `2015soc-from-histsoc`, which is documented MIS-SCALED ACROSS GCMs within that one
     model (gfdl-esm4 fraction-scaled, ukesm1-0-ll percent, with identical `units` and
     `long_name`). See isimip-process-visualize SKILL.md.

  3. sens=`default` (transient CO2) is available for every member, so CO2 treatment is
     uniform -- unlike the csoil ensemble, no mixed-CO2 compromise is needed here.

`jules-es-vn6p3` is absent by fact, not by filter: it publishes NO burntarea at all (0
files) despite appearing in the 3b `fire` sector listing.

The repository API (data.isimip.org) is behind Anubis anti-bot and the CLI downloader needs
`httpx`, which is not installed in this venv -- so URLs are constructed deterministically
(one file per member per scenario, all spanning 2015_2100) and fetched with the stdlib.

Downloads are resumable and VERIFIED BY SIZE against the server's Content-Length. The CLI
downloader's `should_skip` only tests existence, which silently accepts a truncated file;
this script re-checks size on every run and resumes with a Range request. Provenance
(source_url, bytes, ETag) is written to download_provenance.csv, because the anti-bot means
an unrecorded source is not casually re-downloadable.

Usage:
    python scripts/download_burntarea_isimip3b.py [--dry-run] [--jobs N]
"""

import argparse
import csv
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/OutputData"
LAYER_ID = "wildfire-isimip3b_burntarea-total_annual"
OUT_DIR = Path("data/raw") / LAYER_ID
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
VAR = "burntarea-total"

ALL_GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]

# (sector, dir_name, model_token, soc, cadence, gcms)
MEMBERS = [
    ("fire", "CLASSIC", "classic", "2015soc", "monthly", ["gfdl-esm4", "ukesm1-0-ll"]),
    ("fire", "ELM-ECA", "elm-eca", "2015soc-from-histsoc", "monthly", ALL_GCMS),
    ("fire", "LPJmL5-7-10-fire", "lpjml5-7-10-fire", "2015soc", "monthly", ALL_GCMS),
    ("fire", "VISIT", "visit", "2015soc", "monthly", ALL_GCMS),
    ("biomes", "MC2-USFS-r87g5c1", "mc2-usfs-r87g5c1", "nat", "annual", ALL_GCMS),
]

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b burntarea ingest)"
CHUNK = 1 << 20


def build_urls():
    """Deterministic URL per (model, gcm, scenario); every file spans 2015_2100."""
    out = []
    for sector, dirname, model, soc, cadence, gcms in MEMBERS:
        for gcm in gcms:
            for scen in SCENARIOS:
                fname = (
                    f"{model}_{gcm}_w5e5_{scen}_{soc}_default_"
                    f"{VAR}_global_{cadence}_2015_2100.nc"
                )
                url = f"{BASE}/{sector}/{dirname}/{gcm}/future/{fname}"
                out.append((url, fname, model, gcm, scen, cadence, soc))
    return out


def remote_meta(url):
    """HEAD the URL -> (size, etag). Raises on HTTP error."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return int(r.headers["Content-Length"]), r.headers.get("ETag", "")


def fetch(url, dest, expected):
    """Download `url` to `dest`, resuming a partial file. Returns bytes written."""
    have = dest.stat().st_size if dest.exists() else 0
    if have == expected:
        return 0
    if have > expected:  # corrupt/stale -- start over
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
    url, fname, model, gcm, scen, cadence, soc = item
    dest = OUT_DIR / fname
    try:
        size, etag = remote_meta(url)
    except urllib.error.HTTPError as e:
        return dict(ok=False, fname=fname, msg=f"HEAD {e.code}")
    except Exception as e:  # noqa: BLE001 - report, never abort the batch
        return dict(ok=False, fname=fname, msg=f"HEAD {type(e).__name__}: {e}")

    try:
        wrote = fetch(url, dest, size)
    except Exception as e:  # noqa: BLE001
        return dict(ok=False, fname=fname, msg=f"GET {type(e).__name__}: {e}")

    got = dest.stat().st_size if dest.exists() else 0
    if got != size:
        return dict(ok=False, fname=fname, msg=f"size {got} != {size}")

    return dict(
        ok=True, fname=fname, url=url, bytes=size, etag=etag, wrote=wrote,
        model=model, gcm=gcm, scenario=scen, cadence=cadence, soc=soc,
        skipped=(wrote == 0),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    items = build_urls()
    print(f"{len(items)} files, {len(items)//len(SCENARIOS)} members x {len(SCENARIOS)} scenarios")
    if args.dry_run:
        for it in items:
            print(it[0])
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(one, it): it for it in items}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            tag = "skip" if r.get("skipped") else ("ok" if r["ok"] else "FAIL")
            print(f"[{i}/{len(items)}] {tag:4s} {r['fname']}"
                  + ("" if r["ok"] else f"  <-- {r['msg']}"), flush=True)

    good = [r for r in results if r["ok"]]
    bad = [r for r in results if not r["ok"]]
    total = sum(r["bytes"] for r in good)
    print(f"\n{len(good)}/{len(items)} verified, {total/2**30:.2f} GiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "etag", "model", "gcm", "scenario", "cadence", "soc"])
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
