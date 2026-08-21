"""Manifest-driven, verification-hardened downloader for Water Index raw data.

WS1 Stage A2 (waterRiskIndex_beta/v2/PLAN_ws1_ssp_reingestion.md): replaces the
per-variable download_water_*.py skip-on-exists behavior. Every file comes from
an enumeration-derived manifest (url, dest, expected identity); each download
streams to `.part`, is verified against the response Content-Length, hashed
(sha256) while streaming, test-opened as NetCDF, then atomically renamed.
Pre-existing files are re-verified against a HEAD Content-Length, never trusted
on existence alone. A results CSV records what was actually fetched; exit is
nonzero if any manifest row is unresolved.

Usage:
    python scripts/download_water_manifest.py MANIFEST.csv [--base-dir TCFD_ROOT] [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import httpx
from netCDF4 import Dataset

RETRIES = 3
CHUNK = 1024 * 1024


def verify_open(path: Path) -> bool:
    try:
        nc = Dataset(str(path))
        ok = len(nc.variables) > 0
        nc.close()
        return ok
    except Exception:
        return False


def head_length(client: httpx.Client, url: str) -> int | None:
    try:
        r = client.head(url)
        r.raise_for_status()
        return int(r.headers.get("content-length", -1))
    except Exception:
        return None


def fetch(client: httpx.Client, url: str, dest: Path) -> tuple[int, str]:
    """Stream url -> dest.part with sha256; verify length; rename. Returns (bytes, sha256)."""
    part = dest.with_suffix(dest.suffix + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    with client.stream("GET", url) as r:
        r.raise_for_status()
        expected = int(r.headers.get("content-length", -1))
        n = 0
        with open(part, "wb") as f:
            for chunk in r.iter_bytes(CHUNK):
                f.write(chunk)
                h.update(chunk)
                n += len(chunk)
    if expected >= 0 and n != expected:
        part.unlink(missing_ok=True)
        raise IOError(f"length mismatch: got {n}, expected {expected}")
    if not verify_open(part):
        part.unlink(missing_ok=True)
        raise IOError("downloaded file failed NetCDF open check")
    part.rename(dest)
    return n, h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--base-dir", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    base = Path(args.base_dir)
    rows = list(csv.DictReader(open(args.manifest)))
    if args.limit:
        rows = rows[: args.limit]

    results_path = Path(args.manifest).with_name(
        Path(args.manifest).stem + "_results.csv"
    )
    results = []
    failed = 0
    t0 = time.time()

    with httpx.Client(timeout=120, follow_redirects=True) as client:
        for i, row in enumerate(rows, 1):
            url, dest = row["url"], base / row["dest"]
            status, nbytes, digest = "", 0, ""
            for attempt in range(1, RETRIES + 1):
                try:
                    if dest.exists():
                        expected = head_length(client, url)
                        actual = dest.stat().st_size
                        if expected is not None and actual == expected and verify_open(dest):
                            status, nbytes = "verified-existing", actual
                            break
                        print(f"[{i}/{len(rows)}] stale/short existing file, re-fetching: {dest.name}")
                        dest.unlink()
                    nbytes, digest = fetch(client, url, dest)
                    status = "downloaded"
                    break
                except Exception as e:
                    print(f"[{i}/{len(rows)}] attempt {attempt}/{RETRIES} failed: {dest.name}: {e}")
                    time.sleep(2 * attempt)
            else:
                status = "FAILED"
                failed += 1
            elapsed = time.time() - t0
            print(f"[{i}/{len(rows)}] {status:18s} {row['variable']:9s} "
                  f"{row['model']}/{row['gcm']}/{row['scenario']} "
                  f"({nbytes/2**20:.0f} MB, t+{elapsed/60:.1f}m)", flush=True)
            results.append({**row, "status": status, "bytes": nbytes,
                            "sha256": digest,
                            "utc": datetime.now(timezone.utc).isoformat(timespec="seconds")})

    with open(results_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    print(f"\n=== reconciliation: {len(rows) - failed}/{len(rows)} resolved, "
          f"{failed} FAILED; results -> {results_path} ===")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
