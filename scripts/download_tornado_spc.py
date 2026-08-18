#!/usr/bin/env python3
"""Ingest the NOAA SPC severe weather database (tornado tracks) for the CONUS
historical tornado hazard layer.

PRODUCT: TCFD/CDP. This is NOT an ISIMIP ingest -- there is no tornado variable
anywhere in ISIMIP2a/2b/3a/3b, and no tornado-environment proxy can be built from
ISIMIP forcing either (the bias-adjusted atmosphere publishes 11 SURFACE variables
and nothing aloft, so the vertical shear term is simply absent). Enumeration
receipt and the full options review: docs/tornado-data-options-2026-08-18.md.

WHAT THIS DOWNLOADS
    1. SPC "actual tornadoes" database -- one row per tornado, single tracks with
       state segments already merged, so there is no double counting. 1950-2025.
       https://www.spc.noaa.gov/wcm/ -- US Government work, public domain.
    2. Natural Earth 1:50m admin-0 country boundaries (GeoJSON), used to build the
       CONUS coverage mask. scripts/utils/natural_earth.py cannot be used: it
       imports geopandas at module level and geopandas is NOT installed in this
       venv (shapely is).

WHY A COVERAGE MASK IS NOT OPTIONAL
    SPC observes the United States and nothing else. A grid cell in Sonora or
    Manitoba is zero because nobody reported to SPC there, not because tornadoes
    do not occur -- and an unmasked percentile would rank those cells as least
    hazardous. That is the same defect as ranking Bangladesh at the 1st percentile
    on a global grid, just at the border rather than the ocean. The mask makes the
    scored domain equal to the observed domain.

PROVENANCE
    Every file is written with a .json sidecar carrying source_url, sha256, byte
    size and retrieval date. data/ is gitignored and ephemeral, so the sidecar is
    the only record of what was actually ingested (GUARDRAILS 11 -- a claim needs
    its receipt).

USAGE
    python3 scripts/download_tornado_spc.py            # resumable; skips verified files
    python3 scripts/download_tornado_spc.py --force    # re-fetch regardless
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

RAW_DIR = Path("data/raw/tornado-spc")

SOURCES = {
    "spc_tornadoes_1950_2025.csv": {
        "url": "https://www.spc.noaa.gov/wcm/data/1950-2025_actual_tornadoes.csv",
        "what": "SPC severe weather database, tornado tracks, one row per tornado",
        "licence": "US Government work -- public domain",
        "caveat": (
            "SPC: 'these data are used by the NWS for verification purposes and may "
            "not accurately reflect all storm events. Monetary loss information is "
            "highly suspect and should be used with caution, if at all.'"
        ),
    },
    "ne_50m_admin_0_countries.geojson": {
        "url": (
            "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
            "master/geojson/ne_50m_admin_0_countries.geojson"
        ),
        "what": "Natural Earth 1:50m country boundaries, for the CONUS coverage mask",
        "licence": "Natural Earth -- public domain",
        "caveat": "1:50m generalisation; cell-centre containment at 0.25 deg.",
    },
}

USER_AGENT = "TCFD-climate-pipeline/1.0 (research ingest; contact repo owner)"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def already_good(path: Path, sidecar: Path) -> bool:
    """True if the file is present and matches the checksum we recorded for it."""
    if not (path.exists() and sidecar.exists()):
        return False
    try:
        meta = json.loads(sidecar.read_text())
    except json.JSONDecodeError:
        return False
    if meta.get("sha256") != sha256_of(path):
        print(f"  ! {path.name}: checksum mismatch against sidecar -- re-fetching")
        return False
    if meta.get("size_bytes") != path.stat().st_size:
        return False
    return True


def fetch(name: str, spec: dict, force: bool) -> bool:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    dest = RAW_DIR / name
    sidecar = RAW_DIR / f"{name}.json"

    if not force and already_good(dest, sidecar):
        print(f"  = {name}: present and verified ({dest.stat().st_size:,} bytes)")
        return True

    print(f"  > {name}: fetching {spec['url']}")
    req = urllib.request.Request(spec["url"], headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            payload = resp.read()
    except Exception as exc:  # noqa: BLE001 -- report and continue to next source
        print(f"  ! {name}: FAILED -- {exc}")
        return False

    if not payload:
        print(f"  ! {name}: FAILED -- empty response (treat as failure, not as absence)")
        return False

    dest.write_bytes(payload)
    digest = sha256_of(dest)
    sidecar.write_text(
        json.dumps(
            {
                "file": name,
                "source_url": spec["url"],
                "what": spec["what"],
                "licence": spec["licence"],
                "caveat": spec["caveat"],
                "sha256": digest,
                "size_bytes": dest.stat().st_size,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"  + {name}: {dest.stat().st_size:,} bytes, sha256 {digest[:16]}...")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="re-fetch even if present and verified")
    args = ap.parse_args()

    print(f"Ingesting CONUS tornado occurrence sources into {RAW_DIR}/")
    # NOT `all(genexp)` -- that short-circuits on the first failure and silently skips
    # the remaining sources, contradicting the "report and continue" contract in fetch().
    # A partial ingest that reports only the first problem is how a second, unrelated
    # outage gets discovered one run later.
    results = [fetch(name, spec, args.force) for name, spec in SOURCES.items()]
    ok = all(results)

    if not ok:
        print("\nFAILED -- at least one source did not download. Nothing was inferred "
              "from the failure; re-run before drawing any conclusion from missing data.")
        return 1

    print("\nAll sources present and verified. Next: scripts/process_tornado_spc.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
