"""Download the Lange2020 `leh` ensemble for the ISIMIP2b/RCP heatwave layer.

`leh` is the ISIMIP2b heatwave member of the Lange 2020 extreme-event exposure family
(`le{d,r,w,c,h,t}` land-area-exposed, each paired with a `pe*` population-exposed twin).

WHY THIS OLDER DATASET IS BEING INGESTED ALONGSIDE THE NEWER ISIMIP3b `heatwave` LAYER:
the two are not versions of one measurement, they are DIFFERENT INDICES, and this one is
the only ISIMIP heatwave product with an ABSOLUTE, health-anchored criterion.

    ISIMIP2b  Lange2020 / `HWMId-humidex` -> HWMId (relative) AND Humidex >= 45 (absolute)
    ISIMIP3b  Zantout2025 / `HWMID-NONE`  -> HWMId ONLY; the "NONE" is the humidity term

Lange et al. 2020 (Earth's Future, doi:10.1029/2020EF001616) combines a RELATIVE indicator
-- the Heat Wave Magnitude Index daily, which depends only on daily maximum temperature --
with an ABSOLUTE one, the Humidex (daily max and mean temperature plus relative humidity),
so that the events counted "would also adversely affect human health". The reference
implementation states the combination explicitly: land area is exposed "if BOTH a relative
indicator based on temperature and an absolute indicator based on temperature and relative
humidity exceed their respective threshold values". The Humidex threshold is 45, inside
Environment Canada's 40-45 "great discomfort; avoid exertion" band.

READ THAT AS AN INTERSECTION, NOT AS AN ABSOLUTE INDEX. The Humidex gate makes the counted
events health-relevant, but the relative HWMId gate still has to open first, so `leh` ADDS
a health filter on top of a per-cell-distribution framing rather than replacing it.

THE EXPECTED FAILURE MODE IS THE OPPOSITE OF THE 3b LAYER'S, AND IT MUST BE MEASURED FIRST.
The ISIMIP3b `heatwave` layer SATURATES (45.9% of ssp585 2090s cells pinned at frequency
1.0). This one is expected to be SPARSE: Lange et al. report that below 2 degC of global
warming the Humidex "hardly ever exceeds their threshold value of 45 in the temperate
regions", and ISIMIP2b publishes `leh` for rcp26/rcp60 ONLY -- precisely the forcing range
where that gate mostly stays shut. So a near-zero field across Europe, Canada and the
northern US is the outcome to check for BEFORE building anything on this layer, not after.
That is a claim from the paper; the files are 35 MB and it is measurable directly. Run
`scripts/check_leh_nature.py` after this script.

Ensemble: 1 index model x 4 CMIP5 GCMs x 2 RCPs = 8 files, ~35 MB (`leh` alone).

Coverage ENUMERATED 2026-08-13 by listing `Lange2020/`, `HWMId-humidex/`, all 4
`HWMId-humidex/{gcm}/` and all 4 `{gcm}/future/` directories serially (zero empty listings,
zero retries). The variable field was projected at `$(NF-4)` rather than grepped for a
believed-in token (GUARDRAILS 8): 24 files in `future/`, exactly `leh` and `peh` across
picontrol/rcp26/rcp60.

  * COMPLETE matrix -- every (gcm, scenario) combination exists, 4 members per scenario,
    composition IDENTICAL across rcp26/rcp60, so the shared baseline is valid.
  * soc/sens are UNIFORM `nosoc` / `co2` across all files -- no harmonization compromise.
  * Bias adjustment `ewembi`; annual; 2006-2099, so unlike the ISIMIP3b layers this one
    carries a FULL 2010s panel and the baseline decade is index 1, not index 0.
  * NO rcp85. The 2b climate FORCING publishes rcp85 (4 GCMs), but this derived product
    does not, so there is no high-forcing member for this hazard.
  * `peh`, the population-exposed twin, is published 1:1 and is NOT downloaded by default
    -- it answers a different question (people exposed, not land exposed) and would be a
    separate layer. `--include-peh` fetches it.

FILENAME GRAMMAR: ISIMIP2b Lange2020 files carry a LEADING publication token AND use the
`.nc4` extension, where the 3b publications use `.nc`:

    lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099.nc4
    zantout2025_hwmid-none_gfdl-esm4_w5e5_ssp126_2015soc_default_heatwave_global_annual_2015_2100.nc

A `*.nc` filter silently drops the ENTIRE 2b round. Parse from the END regardless.

SERVER DIRECTORY CASING DIFFERS BETWEEN ROUNDS AND IS NOT DERIVABLE: the 2b directory is
`HWMId-humidex` (lower-case d in HWMId) while the 3b one is `HWMID-NONE`. Copy them; do not
compute one from the other or retype them.

DATA NATURE IS UNVERIFIED AT DOWNLOAD TIME and must be value-checked before processing
(GUARDRAILS 9). The 2b Lange2020 `.json` sidecars carry `size` and sha512 but NO
`netcdf_header` block, so units and `long_name` cannot be read without downloading. Do NOT
assume it matches a sibling: within this very publication `led` is binary {0,1} and `let` is
a continuous fraction in [0, 0.952), measured.

Usage:
    python scripts/download_leh_isimip2b.py [--dry-run] [--include-peh]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP2b/DerivedOutputData/Lange2020/HWMId-humidex"
LAYER_ID = "heatwave-isimip2b_leh_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

PUBLICATION = "lange2020"
#: Server directory casing enumerated 2026-08-13, not derived: `HWMId-humidex` on the
#: server, `hwmid-humidex` in the filename. The 3b sibling is `HWMID-NONE`/`hwmid-none`.
MODEL_TOKEN = "hwmid-humidex"
GCMS = ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
SCENARIOS = ["rcp26", "rcp60"]          # no rcp85 for this product
FORCING = "ewembi"
SOC = "nosoc"
SENS = "co2"
YEARS = "2006_2099"
EXT = "nc4"                              # ISIMIP2b publishes .nc4, 3a/3b publish .nc

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP2b leh ingest)"
CHUNK = 1 << 20


def build_items(variables):
    """One item per (variable, gcm, scenario). The enumerated matrix is complete."""
    items = []
    for var in variables:
        for gcm in GCMS:
            for scen in SCENARIOS:
                stem = (f"{PUBLICATION}_{MODEL_TOKEN}_{gcm}_{FORCING}_{scen}_{SOC}_{SENS}_"
                        f"{var}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.{EXT}",
                    url=f"{BASE}/{gcm}/future/{stem}.{EXT}",
                    sidecar=f"{BASE}/{gcm}/future/{stem}.json",
                    variable=var, model=MODEL_TOKEN, gcm=gcm, scenario=scen, soc=SOC,
                    member=f"{MODEL_TOKEN}_{gcm}",
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
    ap.add_argument("--include-peh", action="store_true",
                    help="also fetch the population-exposed twin (a DIFFERENT question and "
                         "a separate layer; not part of the land-exposure ingest)")
    args = ap.parse_args()

    variables = ["leh"] + (["peh"] if args.include_peh else [])
    items = build_items(variables)
    print(f"{len(items)} files = {len(variables)} variable(s) x {len(GCMS)} GCMs "
          f"x {len(SCENARIOS)} scenarios  [{','.join(variables)}]")
    print("index: HWMId (relative) AND Humidex >= 45 (absolute, health-anchored) -- an "
          "INTERSECTION, not an absolute index")
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
    print(f"\n{len(good)}/{len(items)} verified by sha512, {total / 2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "sha512_source", "variable", "member",
                "model", "gcm", "scenario", "soc"])
            w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                row = {k: r[k] for k in w.fieldnames if k != "sha512_source"}
                row["sha512_source"] = "publisher sidecar ({stem}.json), matched"
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
