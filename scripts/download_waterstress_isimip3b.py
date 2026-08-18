"""Download the ISIMIP3b water-stress inputs -- WaterGAP2-2e only.

MODEL SCOPE IS A DELIBERATE, TEMPORARY NARROWING (user decision 2026-08-17).
Four ISIMIP3b models publish an all-sector total withdrawal; three of them
(H08, LPJmL5-7-10-fire, WaterGAP2-2e) also publish the `1850soc` naturalised
supply this layer needs.  We measured all three against each other on one GCM
and found they do NOT differ by a scale factor -- the driest-to-wettest ordering
REVERSES in 68.4% of cells, and the spread reaches 29.7x on the Colorado and
5.3x on the Nile while staying inside 1.1x on the Amazon.  See
docs/water-stress-status-2026-08-17.md section 6 for the measurements.

WaterGAP2-2e was selected because it is the only one of the three that
simulates lake, wetland and floodplain evaporation, and it is correspondingly
the only one whose naturalised Nile (3,053 m3/s vs ~2,700-2,800 observed) and
Murray (232 vs ~380) are near reality; H08 and LPJmL5 overshoot the Nile by
3x and 6x.  REVISIT THIS.  A single-model layer has no structural uncertainty
in its confidence interval at all, and we have now measured that the missing
component is large.

What this pulls, per scenario (ssp126, ssp370, ssp585) x 5 GCMs:
    ptotww  @ 2015soc   total potential withdrawal, all sectors   [kg m-2 s-1]
    qtot    @ 1850soc   total runoff, naturalised                 [kg m-2 s-1]
    dis     @ 1850soc   discharge, naturalised                    [m3 s-1]
plus the fixed cellarea/contfrac grids.

The numerator is 2015soc and the denominator 1850soc BY DESIGN: dividing demand
by post-abstraction discharge would put demand on top and also subtract it from
the bottom.  Pairing is member-wise -- same model, same GCM, same scenario.

Usage:
    python scripts/download_waterstress_isimip3b.py            # plan only
    python scripts/download_waterstress_isimip3b.py --run
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/OutputData/water_global"
MODEL_DIR = "WaterGAP2-2e"
MODEL = "watergap2-2e"
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

# (variable, soc) -- the soc split is the point of the layer, see module docstring.
MONTHLY = [("ptotww", "2015soc"), ("qtot", "1850soc"), ("dis", "1850soc")]
FIXED = ["cellarea", "contfrac"]

OUT = Path("data/raw/water_stress")
RETRIES = 3
BACKOFF = 5


def url_for(gcm: str, scenario: str, soc: str, var: str, ext: str = "nc") -> str:
    stem = f"{MODEL}_{gcm}_w5e5_{scenario}_{soc}_default_{var}_global_monthly_2015_2100"
    return f"{BASE}/{MODEL_DIR}/{gcm}/future/{stem}.{ext}"


def fixed_url(gcm: str, var: str, ext: str = "nc") -> str:
    stem = f"{MODEL}_{gcm}_w5e5_ssp370_2015soc_default_{var}_global_fixed_0000_0000"
    return f"{BASE}/{MODEL_DIR}/{gcm}/future/{stem}.{ext}"


def plan() -> list[tuple[str, Path]]:
    jobs: list[tuple[str, Path]] = []
    for scenario in SCENARIOS:
        for gcm in GCMS:
            for var, soc in MONTHLY:
                u = url_for(gcm, scenario, soc, var)
                jobs.append((u, OUT / Path(u).name))
    for var in FIXED:
        u = fixed_url(GCMS[0], var)
        jobs.append((u, OUT / Path(u).name))
    return jobs


def remote_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except urllib.error.HTTPError:
        return None


def sidecar_checksum(url: str) -> str | None:
    """sha512 from the published .json sidecar, when one exists.

    NOTE: sidecars are not universally reliable across this repository -- CWatM's
    carry a superseded `soc` token and 404 against the real filenames (measured
    2026-08-16).  WaterGAP2-2e's were spot-checked and resolve, but a miss is
    treated as "no checksum available", never as a failure.
    """
    try:
        with urllib.request.urlopen(url.rsplit(".", 1)[0] + ".json", timeout=60) as r:
            return json.load(r).get("checksum")
    except Exception:
        return None


VERIFIED = OUT / "verification.log"


def _record(dest: Path, how: str) -> None:
    """Persist how each file was verified, so provenance is auditable after the fact."""
    VERIFIED.parent.mkdir(parents=True, exist_ok=True)
    with VERIFIED.open("a") as f:
        f.write(f"{dest.name}\t{dest.stat().st_size}\t{how}\n")


def sha512(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha512()
    with path.open("rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> tuple[bool, str]:
    size = remote_size(url)
    if size is None:
        return False, "HEAD failed -- file may not exist"
    if dest.exists() and dest.stat().st_size == size:
        # A size match is NOT integrity. Verify the publisher checksum even for a file we
        # did not fetch this run -- a seeded or resumed file would otherwise never be
        # checked (review round 3). Result is persisted to VERIFIED for provenance.
        want = sidecar_checksum(url)
        if want is None:
            _record(dest, "size-only (no sidecar checksum published)")
            return True, "already complete, size-only"
        if sha512(dest) != want:
            return False, "PRESENT BUT sha512 MISMATCH -- delete and re-fetch"
        _record(dest, f"sha512 {want[:16]}...")
        return True, "already complete, sha512 verified"
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=1800) as r, tmp.open("wb") as f:
                while True:
                    blk = r.read(1 << 22)
                    if not blk:
                        break
                    f.write(blk)
            if tmp.stat().st_size != size:
                raise IOError(f"size {tmp.stat().st_size} != {size}")
            want = sidecar_checksum(url)
            if want and sha512(tmp) != want:
                raise IOError("sha512 mismatch against sidecar")
            tmp.replace(dest)
            _record(dest, f"sha512 {want[:16]}..." if want else "size-only (no sidecar)")
            return True, "ok" + ("" if want else " (no sidecar checksum)")
        except Exception as exc:  # noqa: BLE001 -- retry any transport error
            tmp.unlink(missing_ok=True)
            if attempt == RETRIES:
                return False, f"failed after {RETRIES}: {exc}"
            time.sleep(BACKOFF * attempt)
    return False, "unreachable"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="store_true",
                    help="actually download; default is to print the plan only")
    args = ap.parse_args()

    jobs = plan()
    print(f"WaterGAP2-2e water-stress ingest -- {len(jobs)} files")
    print(f"  {len(SCENARIOS)} scenarios x {len(GCMS)} GCMs x "
          f"{len(MONTHLY)} monthly variables, plus {len(FIXED)} fixed grids")
    print(f"  numerator ptotww@2015soc / denominator qtot,dis@1850soc")
    print(f"  destination {OUT}")
    if not args.run:
        print("\nplan only -- pass --run to download")
        for u, d in jobs[:4]:
            print(f"    {Path(u).name}")
        print(f"    ... and {len(jobs)-4} more")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    ok = bad = 0
    for i, (u, d) in enumerate(jobs, 1):
        good, msg = fetch(u, d)
        ok, bad = ok + good, bad + (not good)
        mb = d.stat().st_size / 1048576 if d.exists() else 0
        print(f"[{i:>2}/{len(jobs)}] {d.name}  {mb:8.1f} MB  {msg}", flush=True)
    print(f"\ndone: {ok} ok, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
