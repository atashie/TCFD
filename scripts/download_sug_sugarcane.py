"""Download the ISIMIP2b LPJmL sugarcane yield ensemble for the TCFD sugarcane layer.

`yield-sug-{noirr,firr}` = simulated sugarcane yield, ISIMIP2b `agriculture` sector,
impact model LPJmL. `sug` is the OUTPUT code for sugarcane; `sgc` is the ISIMIP3b
crop-calendar InputData code and no model publishes output under it (see GUARDRAILS §8
and the 2026-08-11 WORKFLOW-ISSUES entry).

Ensemble: 1 impact model x 4 CMIP5 GCMs = 4 members per scenario per irrigation regime.
LPJmL is the ONLY sugarcane source with future scenarios in the entire repository --
enumerated 2026-08-11 by projecting the variable field over all 150 agriculture
directories in ISIMIP3b/2b/3a/2a (0 empty listings):

  * ISIMIP3b publishes NO sugarcane at all (11 crops: bea cas mai mil pot ri1 ri2 sor
    soy swh wwh). There is therefore no SSP version of this layer.
  * ISIMIP2b: LPJmL only. CLM45/GEPIC/PEPIC carry no sugarcane in future/ or historical/.
  * ISIMIP2a: LPJmL + CLM-Crop, historical only (no scenarios).

Scenario/experiment coverage for `sug`, measured not assumed:
  * rcp26 and rcp60 only. NO rcp85 for this crop -- note that rcp85 DOES exist in 2b
    agriculture, but only from CLM45 and only for mai + soy, so the negative is scoped
    to sugarcane rather than to the round.
  * soc is uniformly `2005soc`; sens `co2` (transient) is available for both scenarios
    while `2005co2` exists for rcp26 only. We take `co2` for both so the two scenarios
    are experimentally symmetric (the skill's "prefer default/transient CO2").

`historical` (1861-2005) is NOT fetched: the future files start in 2006, which already
covers the 2020s baseline decade.

Irrigation regimes are DISTINCT FIELDS, not ensemble members -- `noirr` (rainfed) and
`firr` (fully irrigated) go to separate raw dirs and become separate processed layers.

Units and data nature are NOT knowable here: the sidecars carry `specifiers` but no
`netcdf_header`, so there is no units/long_name to read. Value-check after download per
GUARDRAILS §9 before choosing the decadal statistic.

Downloads are verified by sha512 from the `.json` sidecar (note the sidecar name drops
the `.nc4`: `{stem}.json`, not `{stem}.nc4.json`) and are resumable via Range.

Usage:
    python scripts/download_sug_sugarcane.py [--irrigation noirr|firr|both] [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP2b/OutputData/agriculture/LPJmL"

MODEL = "lpjml"
CROP = "sug"
GCMS = ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
SCENARIOS = ["rcp26", "rcp60"]
SOC = "2005soc"
SENS = "co2"
YEARS = "2006_2099"

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP2b yield-sug ingest)"
CHUNK = 1 << 20


def layer_id(irrigation):
    """Raw dir name mirrors the processed layer name (CLAUDE.md 'Downloading')."""
    return f"sugarcane_yield-{CROP}-{irrigation}_annual"


def build_items(irrigations):
    items = []
    for irr in irrigations:
        for gcm in GCMS:
            for scen in SCENARIOS:
                stem = (f"{MODEL}_{gcm}_ewembi_{scen}_{SOC}_{SENS}_"
                        f"yield-{CROP}-{irr}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.nc4",
                    url=f"{BASE}/{gcm}/future/{stem}.nc4",
                    sidecar=f"{BASE}/{gcm}/future/{stem}.json",
                    gcm=gcm, scenario=scen, irrigation=irr,
                    member=f"{MODEL}_{gcm}",
                    out_dir=Path("data/raw") / layer_id(irr),
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
    item["out_dir"].mkdir(parents=True, exist_ok=True)
    dest = item["out_dir"] / item["fname"]
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

    # The sidecar's parsed specifiers are the authority on what this file IS -- assert
    # rather than trust the filename we constructed.
    spec = meta.get("specifiers", {})
    for key, expect in (("crop", CROP), ("irrigation", item["irrigation"]),
                        ("climate_scenario", item["scenario"]), ("variable", "yield")):
        if spec.get(key) != expect:
            return dict(item, ok=False, msg=f"specifier {key}={spec.get(key)!r} != {expect!r}")

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
    ap.add_argument("--irrigation", choices=["noirr", "firr", "both"], default="both")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    irrigations = ["noirr", "firr"] if args.irrigation == "both" else [args.irrigation]
    items = build_items(irrigations)
    print(f"{len(items)} files = {len(GCMS)} members x {len(SCENARIOS)} scenarios "
          f"x {len(irrigations)} irrigation regime(s)")
    if args.dry_run:
        for it in items:
            print(it["url"])
        return 0

    # Serial: parallel requests to files.isimip.org get rate-limited into empty or
    # partial responses (see the isimip-search-download skill).
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

    for irr in irrigations:
        rows = [r for r in good if r["irrigation"] == irr]
        if not rows:
            continue
        out = Path("data/raw") / layer_id(irr) / "download_provenance.csv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "member", "gcm", "scenario", "irrigation"])
            w.writeheader()
            for r in sorted(rows, key=lambda x: x["fname"]):
                w.writerow({k: r[k] for k in w.fieldnames})
        print(f"provenance -> {out}")

    if bad:
        print(f"\n{len(bad)} FAILED (re-run to resume):")
        for r in bad:
            print(f"  {r['fname']}: {r['msg']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
