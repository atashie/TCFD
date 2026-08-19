#!/usr/bin/env python3
"""Ingest the World Bank / GFDRR "Global landslide hazard map" (Arup, 2021).

PRODUCT: TCFD/CDP. This is NOT an ISIMIP ingest. ISIMIP publishes no landslide,
mass-movement, slope or geotechnical output in any round -- enumeration receipt in
`config/isimip_search_catalog.yaml` under `negative_results.landslide`, full options
review in docs/landslide-data-options-2026-08-19.md.

WHY THIS SOURCE AND NOT ANOTHER
    Of the global landslide products enumerated 2026-08-18/19, this is the ONLY one
    that publishes a RATE with physical units (annual frequency of significant
    landslides per km^2). NASA's global susceptibility map is a 5-class rank; GIRI is
    a 5-class scenario hazard behind registration; the Duan et al. CMIP6 projections
    publish no grids at all. A rank cannot be aggregated to a coarser cell without
    inventing a meaning for the average of an ordinal; a rate can, because the areal
    mean of a per-km^2 frequency is arithmetic rather than a modelling choice.

    See the options doc for what this costs: this layer is HISTORICAL ONLY (1980-2018)
    and carries no scenario axis, so it cannot answer a forward-looking question.

LICENCE -- CLEARED 2026-08-19, WITH THE SOURCE AMBIGUITY RECORDED
    Cleared by the user for our limited commercial use cases on 2026-08-19. ATTRIBUTION
    to World Bank / GFDRR and Arup is required wherever a value from this layer is
    published. The underlying ambiguity is retained below because it is a real defect in
    the source's own records, and anyone re-deriving the provenance will hit it:

    Two publisher-side catalogue records for the same dataset disagree:

        datacatalog.worldbank.org  (dataset 0037584)  -> "Creative Commons
                                                          Attribution-Non Commercial 4.0"
        energydata.info CKAN mirror                   -> "CC-BY-4.0"

    The Arup project report itself carries no licence statement in its front or back
    matter (checked 2026-08-19, 113 pages). NC and BY are materially different for a
    commercial climate-risk product. The user's 2026-08-19 determination resolves it for
    our use; the sidecars and the processed file record both the determination and the
    discrepancy. Recorded here and in `source_licence_status` rather than in a caveat
    attribute on purpose -- caveat attributes are promoted verbatim into customer
    reports, and our licensing position is not the customer's concern.

WHAT THIS DOWNLOADS
    1. LS_RF_Median_1980-2018_COG.tif -- rainfall-triggered landslide hazard, median
       over 1980-2018, cloud-optimised GeoTIFF, EPSG:4326, 30 arcsec (~1 km),
       43201 x 15779, float32, 124 MB.
    2. global-landslide-hazard-map-report.pdf -- the method report (121 MB).

    The EARTHQUAKE-triggered layer (ls_eq_tiled.tif) is deliberately NOT ingested.
    Seismic triggering is stationary under every emissions pathway, and mixing it into
    a climate-risk layer would attribute a tectonic hazard to a climate driver. It is
    listed in the options doc as a separate product for a separate question.

PROVENANCE
    Every file gets a .json sidecar carrying source_url, sha256, byte size and
    retrieval date. data/ is gitignored and ephemeral, so the sidecar is the only
    record of what was actually ingested (GUARDRAILS 11 -- a claim needs its receipt).

USAGE
    python3 scripts/download_landslide_arup.py            # resumable; skips verified files
    python3 scripts/download_landslide_arup.py --force    # re-fetch regardless
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

RAW_DIR = Path("data/raw/landslide-arup")

BASE = "https://datacatalogfiles.worldbank.org/ddh-published/0037584"

FILES = {
    "LS_RF_Median_1980-2018_COG.tif": f"{BASE}/DR0045419/LS_RF_Median_1980-2018_COG.tif",
    "global-landslide-hazard-map-report.pdf": f"{BASE}/DR0045411/global-landslide-hazard-map-report.pdf",
}

#: Not ingested, recorded so the decision is visible rather than looking like an oversight.
DELIBERATELY_NOT_INGESTED = {
    "ls_eq_tiled.tif": "earthquake-triggered -- stationary under emissions pathways, separate question",
    "LS_RF_Mean_1980-2018_COG.tif": "mean over 1980-2018; the median is the published headline layer",
    "LS_TH_COG.tif": "'TH ranks' -- an ordinal, not a rate",
}

DATASET_PAGE = "https://datacatalog.worldbank.org/search/dataset/0037584/global-landslide-hazard-map"


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fetch(name: str, url: str, force: bool) -> bool:
    """Download one file and write its sidecar. Returns True if it is present and verified."""
    dest = RAW_DIR / name
    side = dest.with_suffix(dest.suffix + ".json")

    if dest.exists() and side.exists() and not force:
        meta = json.loads(side.read_text())
        if dest.stat().st_size == meta.get("bytes"):
            print(f"  SKIP  {name} ({dest.stat().st_size:,} B, sidecar matches)")
            return True
        print(f"  STALE {name} -- size differs from sidecar, re-fetching")

    print(f"  GET   {name}")
    req = urllib.request.Request(url, headers={"User-Agent": "TCFD-pipeline/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=900) as resp:
            declared = resp.headers.get("Content-Length")
            with dest.open("wb") as out:
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
    except Exception as exc:  # noqa: BLE001
        print(f"  FAIL  {name}: {exc}")
        return False

    size = dest.stat().st_size
    if declared is not None and int(declared) != size:
        # A truncated download is the failure mode that looks like success downstream.
        print(f"  FAIL  {name}: got {size:,} B, server declared {declared} B")
        dest.unlink(missing_ok=True)
        return False

    side.write_text(json.dumps({
        "source_url": url,
        "dataset_page": DATASET_PAGE,
        "bytes": size,
        "sha256": sha256_of(dest),
        "retrieved_on": dt.date.today().isoformat(),
        "licence_status": (
            "CLEARED for our limited commercial use -- user determination 2026-08-19. "
            "Attribution to World Bank / GFDRR and Arup required. Source records disagree "
            "(WB DDH: CC BY-NC 4.0; energydata.info mirror: CC-BY-4.0; report: neither)."
        ),
    }, indent=2) + "\n")
    print(f"  OK    {name} ({size:,} B)")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="re-fetch even if present and verified")
    args = ap.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Ingesting the Arup global landslide hazard map into {RAW_DIR}")
    print("LICENCE: cleared for our limited commercial use 2026-08-19; "
          "ATTRIBUTION to World Bank / GFDRR and Arup is required.\n")

    ok = all(fetch(name, url, args.force) for name, url in FILES.items())

    print("\nDeliberately not ingested:")
    for name, why in DELIBERATELY_NOT_INGESTED.items():
        print(f"  {name:34s} {why}")

    if not ok:
        print("\nFAILED -- at least one file did not verify.")
        return 1
    print("\nAll files present and verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
