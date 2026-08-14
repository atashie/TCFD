#!/usr/bin/env python3
"""Download the ISIMIP3b TipESM2025 `fldfrcmax` ensemble (CaMa-Flood inundation).

    python scripts/download_fldfrcmax_isimip3b.py --dry-run      # plan only (default)
    python scripts/download_fldfrcmax_isimip3b.py --run
    python scripts/download_fldfrcmax_isimip3b.py --run --protection none flopros

WHAT THIS IS
------------
Annual MAXIMUM flooded fraction of each cell, from CaMa-Flood v4.0.0 routing of seven
ISIMIP3b GHMs onto a 15 arcmin (0.25 deg) floodplain. Unlike `floodedarea` -- the 0.5 deg
binary exposure flag that reads 0.000 across the Amazon because it scores departure from a
preindustrial reference -- this is an absolute inundation field. Measured on the probe
member, the Amazon main stem box runs 0.0743 unprotected.

It is NOT the same variable as ISIMIP2b Zimmer2023 `fldfrc`: that is a mean flooded
fraction at 150 arcsec under RCP forcing, this is an annual maximum at 15 arcmin under SSP.
The two cannot be pooled.

DIRECTORY SHAPE -- FIVE LEVELS, SCENARIO ABOVE MODEL
    TipESM2025/MIT/{scenario}/{GHM}/{gcm}/
not the {MODEL}/{gcm}/{period}/ that every other publication uses. A three-level walker
returns this publication as EMPTY, which is a wrong-depth signal, not a negative
(GUARDRAILS 8). This script enumerates the leaves at run time rather than constructing
paths, so a change in composition surfaces as a diff instead of a silent 404.

ENSEMBLE COMPOSITION -- ENUMERATED 2026-08-14, 96 leaves, 0 empty, 915 files
    protection  none                    -> 7 GHMs, 32 members per scenario
    protection  2yr / 40yr / flopros / hanze -> 6 GHMs, 27 members (CWATM publishes `none` ONLY)
So the unprotected field and the protected fields have DIFFERENT ensembles. A `none` minus
`flopros` difference is not a clean defence effect unless it is restricted to the common 27.

    soc   HETEROGENEOUS, and no single value covers every model:
            CLASSIC {2015soc, 2015soc-from-histsoc}   H08 {1850soc, 2015soc}
            JULES-W2 {2015soc, 2015soc-from-histsoc}  MIROC-INTEG-LAND {1850soc, 2015soc-from-histsoc}
            WATERGAP2-2E {all three}                  WEB-DHM-SG {2015soc}
            CWATM {2015soc-from-histsoc}
          Rule applied here: prefer `2015soc`, fall back to `2015soc-from-histsoc`, NEVER
          `1850soc` (a preindustrial-society counterfactual, not a projection). That keeps
          all seven models and puts two of them on the fallback -- a declared heterogeneity,
          the same compromise `led` records. Harmonising strictly on 2015soc would drop
          MIROC-INTEG-LAND and CWATM.

    sens  `default` only. `2015co2` exists for ssp585 ONLY (75 files, 3 GHMs), so it cannot
          enter a scenario-symmetric ensemble. It is a CO2-sensitivity experiment; ingest it
          separately if that question is ever asked.

TWO DEFECTS THE PROBE FOUND -- BOTH SILENT
  1. LATITUDE ORIENTATION DIFFERS BETWEEN PROTECTION VARIANTS OF THE SAME MEMBER, and it
     is systematic: verified across all 258 ingested files, `none` is DESCENDING
     (89.875 -> -89.875) in all 96, `40yr`/`flopros` ASCENDING in all 162.
     What differs is the COORDINATE order (and the order the file declares its dimensions,
     visible in `ds.sizes`). The variable's own dim order is ('time','lat','lon')
     everywhere -- the arrays are NOT transposed, so do not go looking for one.
     A `.sel(lat=slice(hi, lo))` silently returns ZERO cells on whichever half does not
     match, and code that indexes by position rather than coordinate builds an upside-down
     layer that passes every contract check. Always `.sortby("lat")` on load.
  2. PERMANENT WATER READS 1.000. Cells fully flooded in every year: unprotected median
     511 (494-531), protected median ~230. The Caspian and other inland water bodies. Not
     flood risk, and they will sit at the top of any percentile ranking. Needs a mask
     decision before processing.

Neither is visible in the metadata, which is otherwise the cleanest of any layer here:
`units='1'`, `long_name='flooded fraction'`, and a full provenance block naming
protection_level, resolution, soc and sens.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/TipESM2025/MIT"
LAYER_ID = "flood-isimip3b_fldfrcmax_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "fldfrcmax"
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
GHMS = ["CLASSIC", "CWATM", "H08", "JULES-W2", "MIROC-INTEG-LAND", "WATERGAP2-2E", "WEB-DHM-SG"]
#: In preference order. `1850soc` is deliberately absent -- see the docstring.
SOC_PREFERENCE = ["2015soc", "2015soc-from-histsoc"]
SENS = "default"
DEFAULT_PROTECTION = ["none", "40yr", "flopros"]
ALL_PROTECTION = ["none", "2yr", "40yr", "flopros", "hanze"]

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b fldfrcmax ingest)"
CHUNK = 1 << 20
RETRIES = 3


def _get(url: str, binary: bool = False):
    last = None
    for _ in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read() if binary else resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = exc
        except Exception as exc:  # noqa: BLE001 - network flake, retry
            last = exc
    raise RuntimeError(f"{url}: {last}")


def list_files(url: str) -> list[str]:
    """Directory listing, asserted non-empty by the caller.

    An empty body here is a FAILURE SIGNAL, not a negative -- parallel or rate-limited
    requests return empty listings indistinguishable from absent data.
    """
    body = _get(url)
    if body is None:
        return []
    return re.findall(r'href="([^"]+\.nc4?)"', body)


def parse(name: str) -> dict:
    """Parse from the END. This publication has no bias-adjustment token, so its field
    count differs from Zimmer2023's -- never carry an offset between publications."""
    f = name.split("_")
    return {
        "scenario": f[-10], "soc": f[-9], "sens": f[-8], "variable": f[-7],
        "protection": f[-6], "resolution": f[-5], "start": f[-2], "end": f[-1].split(".")[0],
    }


def enumerate_members(protections: list[str]) -> tuple[list[dict], list[str]]:
    """Walk every leaf and select ONE soc per (scenario, GHM, gcm), held constant across
    every protection level that model publishes.

    Choosing soc per (member, protection) independently looks harmless and is not:
    MIROC-INTEG-LAND publishes `2015soc` for `none` but only `1850soc` /
    `2015soc-from-histsoc` for the protected variants, so a per-protection preference gives
    that model a DIFFERENT SIMULATION unprotected than protected -- and every
    none-minus-flopros statement about it would then conflate defences with the
    socioeconomic run. We intersect the soc sets across the protection levels the model
    actually publishes, then apply the preference to the intersection.
    """
    candidates: dict[tuple, dict[str, dict[str, str]]] = {}
    warnings: list[str] = []
    for scenario in SCENARIOS:
        for ghm in GHMS:
            gcm_body = _get(f"{BASE}/{scenario}/{ghm}/")
            if gcm_body is None:
                warnings.append(f"{scenario}/{ghm}: 404")
                continue
            gcms = [d.rstrip("/") for d in re.findall(r'href="([^"]+/)"', gcm_body)
                    if not d.startswith(".")]
            for gcm in gcms:
                url = f"{BASE}/{scenario}/{ghm}/{gcm}/"
                names = list_files(url)
                if not names:
                    warnings.append(f"EMPTY LISTING (failure signal, not a negative): {url}")
                    continue
                for name in names:
                    p = parse(name)
                    if (p["variable"] != VAR or p["sens"] != SENS
                            or p["protection"] not in protections
                            or p["soc"] not in SOC_PREFERENCE):
                        continue
                    candidates.setdefault((scenario, ghm, gcm), {}) \
                              .setdefault(p["protection"], {})[p["soc"]] = url + name

    members: list[dict] = []
    for (scenario, ghm, gcm), by_prot in candidates.items():
        common = set.intersection(*(set(s) for s in by_prot.values()))
        if common:
            soc = next(s for s in SOC_PREFERENCE if s in common)
            socs = {prot: soc for prot in by_prot}
        else:
            socs = {prot: next(s for s in SOC_PREFERENCE if s in by_prot[prot])
                    for prot in by_prot}
            warnings.append(
                f"{scenario}/{ghm}/{gcm}: no soc common to {sorted(by_prot)} -- using "
                f"{socs}. Protection differences for this member are NOT clean.")
        for prot, soc in socs.items():
            members.append({"url": by_prot[prot][soc], "name": by_prot[prot][soc].rsplit("/", 1)[-1],
                            "scenario": scenario, "ghm": ghm, "gcm": gcm,
                            "soc": soc, "protection": prot})
    return sorted(members, key=lambda m: (m["protection"], m["scenario"], m["ghm"], m["gcm"])), warnings


def check_scenario_symmetry(members: list[dict]) -> list[str]:
    """A shared 2020s baseline is only valid if every scenario has the same members."""
    problems = []
    by_prot: dict[str, dict[str, set]] = {}
    for m in members:
        by_prot.setdefault(m["protection"], {}).setdefault(m["scenario"], set()).add(
            (m["ghm"], m["gcm"], m["soc"]))
    for prot, per_scen in sorted(by_prot.items()):
        sets = list(per_scen.values())
        if not all(s == sets[0] for s in sets):
            for scen, s in per_scen.items():
                missing = sets[0] - s
                if missing:
                    problems.append(f"{prot}/{scen} missing {sorted(missing)}")
        counts = {s: len(v) for s, v in per_scen.items()}
        print(f"  {prot:<9} members per scenario: {counts}")
    return problems


def sidecar(url: str) -> dict | None:
    body = _get(re.sub(r"\.nc4?$", ".json", url))
    return json.loads(body) if body else None


def download(members: list[dict]) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows, ok, skipped, failed = [], 0, 0, 0
    for i, m in enumerate(members, 1):
        dest = OUT_DIR / m["name"]
        side = sidecar(m["url"])
        want = side.get("checksum") if side else None
        size = side.get("size") if side else None
        if dest.exists() and (size is None or dest.stat().st_size == size):
            if want is None or hashlib.sha512(dest.read_bytes()).hexdigest() == want:
                skipped += 1
                rows.append({**{k: m[k] for k in ("scenario", "ghm", "gcm", "soc", "protection")},
                             "file": m["name"], "url": m["url"], "bytes": dest.stat().st_size,
                             "sha512": want or "", "sha512_source": "sidecar" if want else "none"})
                continue
        blob = _get(m["url"], binary=True)
        if blob is None:
            print(f"  [{i}/{len(members)}] 404 {m['name']}", file=sys.stderr)
            failed += 1
            continue
        got = hashlib.sha512(blob).hexdigest()
        if want and got != want:
            print(f"  [{i}/{len(members)}] CHECKSUM MISMATCH {m['name']}", file=sys.stderr)
            failed += 1
            continue
        dest.write_bytes(blob)
        ok += 1
        print(f"  [{i}/{len(members)}] {m['name']}  {len(blob)/1e6:.1f} MB", flush=True)
        rows.append({**{k: m[k] for k in ("scenario", "ghm", "gcm", "soc", "protection")},
                     "file": m["name"], "url": m["url"], "bytes": len(blob),
                     "sha512": got, "sha512_source": "sidecar" if want else "computed-locally"})
    with (OUT_DIR / "provenance.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["scenario", "ghm", "gcm", "soc", "protection",
                                           "file", "url", "bytes", "sha512", "sha512_source"])
        w.writeheader()
        w.writerows(rows)
    print(f"\ndownloaded={ok} skipped={skipped} failed={failed} -> {OUT_DIR}")
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true", help="actually download (default is a plan)")
    ap.add_argument("--dry-run", action="store_true", help="explicit no-op, the default")
    ap.add_argument("--protection", nargs="+", choices=ALL_PROTECTION, default=DEFAULT_PROTECTION)
    args = ap.parse_args()

    print(f"Enumerating {len(SCENARIOS)} scenarios x {len(GHMS)} GHMs (96 leaves)...")
    members, warnings = enumerate_members(args.protection)
    for w in warnings:
        print(f"  WARNING {w}", file=sys.stderr)

    print(f"\nSelected {len(members)} files "
          f"(protection={' '.join(args.protection)}, sens={SENS}, soc preference={SOC_PREFERENCE}):")
    problems = check_scenario_symmetry(members)
    for p in problems:
        print(f"  ASYMMETRY {p}", file=sys.stderr)
    soc_used: dict[str, set] = {}
    for m in members:
        soc_used.setdefault(m["ghm"], set()).add(m["soc"])
    print("\n  soc actually selected per model:")
    for ghm in sorted(soc_used):
        print(f"    {ghm:<20}{', '.join(sorted(soc_used[ghm]))}")

    if not args.run:
        print("\nPLAN ONLY. Re-run with --run to download.")
        return 0
    if problems:
        print("\nRefusing to download: scenario composition is asymmetric, so a shared "
              "2020s baseline would not be valid. Resolve the asymmetry first.", file=sys.stderr)
        return 2
    return download(members)


if __name__ == "__main__":
    raise SystemExit(main())
