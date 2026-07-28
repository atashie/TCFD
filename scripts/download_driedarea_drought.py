"""Ingest ISIMIP3b Heinicke2026 `driedarea` (drought exposure) into S3 raw.

SOURCE: ISIMIP3b/DerivedOutputData/Heinicke2026/{MODEL}/{gcm}/future/ --
    {model}_{gcm}_w5e5_{ssp}_2015soc_default_driedarea_global_annual_2015_2100.nc

45 files = 3 GHMs x 5 CMIP6 GCMs x {ssp126, ssp370, ssp585}, ~3.8 MB each (~170 MB).

WHY THIS SOURCE
---------------
`driedarea` is the ISIMIP3b/SSP sibling of the ISIMIP2b Lange2020 `led` exposure
member: same binary per-cell annual drought flag, same 0.5 deg grid, but on CMIP6
GCMs and the SSP scenario family rather than CMIP5/RCP. It is published under a
DIFFERENT publication directory (Heinicke2026, not Lange2020), which is why a
search for the variable name `led` across rounds finds nothing in 3b -- the
publication changes, not just the variable. See GUARDRAILS §8 and
config/isimip_search_catalog.yaml search_results.drought.

`historical` (1850-2014) is deliberately NOT ingested: the shared 2020s baseline
lives inside the future files, as for the csoil and burntarea layers.

`floodedarea` sits in the same directories and is handled by a separate layer --
do not fold it in here.

VALUE-CHECKED BEFORE INGEST (2026-07-28, all 45 files)
------------------------------------------------------
Every member is BINARY {0,1} -- n_unique == 2, range [0,1], zero intermediate
values -- on calendar `proleptic_gregorian` with time units `days since 1601-01-01`.
The three GHMs do NOT share a land mask:

    h08           53,713 land cells
    jules-w2      57,523
    watergap2-2e  56,410
    union         63,455      intersection  46,024

so 17,431 cells (27% of the union) are missing at least one model. Those cells are
KEPT, not masked; the processor emits per-cell `n_members` / `n_models` so thin
coverage is auditable downstream (csoil precedent).

Usage:
    python scripts/download_driedarea_drought.py [--dry-run] [--force]
"""

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402

LAYER_ID = "drought_driedarea_annual"
BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/Heinicke2026"

# {directory name: filename token} -- the file server capitalises model directories
MODELS = {"H08": "h08", "JULES-W2": "jules-w2", "WaterGAP2-2e": "watergap2-2e"}
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
EXPECTED = len(MODELS) * len(GCMS) * len(SCENARIOS)  # 45


def log(msg):
    print(msg, flush=True)


def members():
    """Enumerate every (url, filename) pair in the expected 3 x 5 x 3 matrix."""
    for directory, model in MODELS.items():
        for gcm in GCMS:
            for scen in SCENARIOS:
                name = (f"{model}_{gcm}_w5e5_{scen}_2015soc_default"
                        f"_driedarea_global_annual_2015_2100.nc")
                yield {
                    "model": model,
                    "gcm": gcm,
                    "scenario": scen,
                    "name": name,
                    "url": f"{BASE}/{directory}/{gcm}/future/{name}",
                }


def fetch(url, dest, tries=3):
    """Download to `dest`, retrying transient failures."""
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            return
        except Exception as exc:  # noqa: BLE001 - retry any transport error
            if attempt == tries:
                raise
            log(f"    retry {attempt}/{tries - 1} after {type(exc).__name__}: {exc}")
            time.sleep(3 * attempt)


def sidecar(url, tries=3):
    """Fetch the ISIMIP `.json` sidecar carrying `size` and `sha512`."""
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url[:-3] + ".json", timeout=120) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            if attempt == tries:
                raise
            time.sleep(3 * attempt)


def sha512_of(path):
    h = hashlib.sha512()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="re-download and re-upload even if already present")
    args = ap.parse_args()

    todo = list(members())
    assert len(todo) == EXPECTED, f"expected {EXPECTED} members, built {len(todo)}"

    if args.dry_run:
        for m in todo[:5]:
            log(f"  {m['url']}")
        log(f"  ... ({len(todo)} total)")
        return

    cache = storage.cache_root() / LAYER_ID
    cache.mkdir(parents=True, exist_ok=True)

    fs = storage.s3_filesystem()
    prefix = storage._p(f"{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}")
    present = ({Path(k).name for k in fs.glob(f"{prefix}/*.nc")}
               if fs.exists(prefix) else set())
    log(f"{len(present)} of {EXPECTED} already in S3 raw")

    paths, urls, failures = [], {}, []
    for i, m in enumerate(todo, 1):
        dest = cache / m["name"]
        tag = f"{m['model']:<13} {m['gcm']:<14} {m['scenario']}"

        if m["name"] in present and not args.force:
            log(f"[{i}/{EXPECTED}] skip  {tag}  (already in S3)")
            continue

        try:
            meta = sidecar(m["url"])
            want_size = meta["size"]
            want_sha = meta["checksum"]
            assert meta["checksum_type"] == "sha512", meta["checksum_type"]

            # Resumable: a cached copy that already matches the sidecar is reused,
            # which is what lets the pre-seeded scratch downloads count.
            if not (dest.exists() and dest.stat().st_size == want_size) or args.force:
                fetch(m["url"], dest)

            got_size = dest.stat().st_size
            if got_size != want_size:
                raise ValueError(f"size {got_size} != sidecar {want_size}")
            got_sha = sha512_of(dest)
            if got_sha != want_sha:
                raise ValueError(f"sha512 mismatch\n  got  {got_sha}\n  want {want_sha}")

            paths.append(dest)
            urls[m["name"]] = m["url"]
            log(f"[{i}/{EXPECTED}] ok    {tag}  {got_size / 1e6:.1f} MB  sha512 verified")
        except Exception as exc:  # noqa: BLE001 - report and continue, fail at the end
            failures.append((tag, str(exc)))
            log(f"[{i}/{EXPECTED}] ERROR {tag}  {exc}")

    if paths:
        log(f"\nUploading {len(paths)} verified members to "
            f"s3://{storage.BUCKET}/{storage.raw_prefix(LAYER_ID)}/ ...")
        entries = storage.ingest_raw(paths, LAYER_ID, source_urls=urls)
        log(f"  ingested {len(entries)} objects")

    final = ({Path(k).name for k in fs.glob(f"{prefix}/*.nc")}
             if fs.exists(prefix) else set())
    log(f"\nraw/isimip/{LAYER_ID}/ now holds {len(final)} of {EXPECTED} members")

    if failures:
        log(f"\n{len(failures)} FAILED:")
        for tag, err in failures:
            log(f"  {tag}: {err}")
        raise SystemExit(1)
    if len(final) != EXPECTED:
        raise SystemExit(f"incomplete ingest: {len(final)}/{EXPECTED}")
    log("\nDone.")


if __name__ == "__main__":
    main()
