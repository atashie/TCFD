"""Download the Zantout2025 `heatwave` ensemble for the ISIMIP3b/SSP heatwave layer.

`heatwave` is the ISIMIP3b re-issue of the Lange 2020 heatwave-exposure concept (`leh`),
renamed by hazard word rather than `le*` code. It is a DIFFERENT product from ISIMIP2b
`leh`, and in this case not only a different ensemble but a different INDEX -- see
"THE INDEX CHANGED BETWEEN ROUNDS" below, which is the single most important fact about
this layer.

Ensemble: 1 impact model x 5 CMIP6 GCMs x 3 SSPs = 15 files, ~63 MB.

Coverage ENUMERATED 2026-08-13 by listing `Zantout2025/`, `HWMID-NONE/`, all 5
`HWMID-NONE/{gcm}/` and all 5 `HWMID-NONE/{gcm}/future/` directories serially (zero empty
listings, zero retries). The variable field was projected at `$(NF-4)` rather than grepped
for a believed-in token (GUARDRAILS 8): 20/20 files in `future/` are `heatwave`.

  * COMPLETE matrix -- every (gcm, scenario) combination exists. 5 members per scenario and
    the composition is IDENTICAL across ssp126/ssp370/ssp585, so the shared 2020s baseline
    is valid on a uniform 5-member ensemble.
  * soc/sens are UNIFORM `2015soc` / `default` across all 15 files -- no harmonization
    compromise, no member dropped to make the ensemble uniform.
  * Bias adjustment `w5e5`; annual; 2015-2100, so the 2020s baseline decade is complete
    (there is no complete 2010s decade -- 3b future runs start in 2015).
  * ONE impact model. n_models = 1 for every member, so the ">=2 impact models"
    publication mask used by `led` HAS NO MEANING HERE. A mask, if any, must be defined on
    GCM count and declared as a different rule -- do not inherit one.

THE INDEX CHANGED BETWEEN ROUNDS -- this is not a version bump.

    ISIMIP2b  Lange2020 / `HWMId-humidex`  -> HWMId (relative) AND Humidex (absolute)
    ISIMIP3b  Zantout2025 / `HWMID-NONE`   -> HWMId ONLY; the "NONE" is the humidity term

Lange et al. 2020 (Earth's Future, doi:10.1029/2020EF001616) combines a RELATIVE indicator
-- the Heat Wave Magnitude Index daily (HWMId), which depends only on daily maximum
temperature -- with an ABSOLUTE indicator, the Humidex, which additionally uses daily mean
temperature and relative humidity, so that the events counted "would also adversely affect
human health". Zantout et al. 2025 (Nat. Commun., doi:10.1038/s41467-025-65600-7, Methods)
drops the absolute criterion: a grid cell is exposed "if the HWMId of that year exceeds the
97.5th percentile of the picontrol distribution of that grid cell". No humidity term.

Two consequences for anything customer-facing:

  * This layer measures UNPRECEDENTED heat relative to a per-cell preindustrial control
    distribution, NOT dangerous-to-humans heat. A hot-and-dry site that has always been hot
    can read low; a temperate site warming past its own preindustrial spread reads high.
    That is the same relative-baseline framing as the drought layers and must be disclosed
    the same way.
  * It is expected to flag MORE area than the 2b `leh` product, especially in arid and
    high-latitude cells where the Humidex gate would never have opened. Do not describe a
    3b-vs-2b difference as a trend or an improvement in the hazard.

DATA NATURE IS UNVERIFIED AT DOWNLOAD TIME and must be value-checked before processing
(GUARDRAILS 9). The published header declares `long_name = "heatwave area share"` and
`units = "1"`, which reads continuous -- but the Zantout 2025 Methods say "The exposed grid
cell area fraction is one if the percentile threshold is exceeded and zero else", i.e.
BINARY. A paper is not a measurement either: `floodedarea` (3b, sibling publication) also
declares units="1", is binary, and is non-NaN over 94.7% of the globe INCLUDING OCEAN.
Run `scripts/check_heatwave_nature.py` after this script and before any processor, and
check the land mask as well as the value set.

SIDECARS EXIST -- at `{stem}.json`, NOT `{stem}.nc.json`. Verified 2026-08-13: the `.json`
form returns HTTP 200 with `size`, sha512 `checksum` and a full `netcdf_header`, while the
`.nc.json` form returns 404 for the same file. Downloads here are CHECKSUM-verified against
the publisher. (This trap is worth knowing: `.nc.json` 404s look exactly like "this
publication has no sidecars".)

Usage:
    python scripts/download_heatwave_isimip3b.py [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/DerivedOutputData/Zantout2025"
LAYER_ID = "heatwave-isimip3b_heatwave_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "heatwave"
PUBLICATION = "zantout2025"
#: (server directory, filename model token) -- both enumerated 2026-08-13, not derived by
#: lowercasing. The server directory is `HWMID-NONE` while ISIMIP2b's sibling directory is
#: `HWMId-humidex`; the server is case-inconsistent between rounds. Copy, do not retype.
MODELS = [("HWMID-NONE", "hwmid-none")]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
FORCING = "w5e5"
SOC = "2015soc"
SENS = "default"
YEARS = "2015_2100"

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b heatwave ingest)"
CHUNK = 1 << 20


def build_items():
    """One item per (model, gcm, scenario). The enumerated matrix is complete: no gaps."""
    items = []
    for mdir, model in MODELS:
        for gcm in GCMS:
            for scen in SCENARIOS:
                # Zantout2025 filenames DO carry a leading publication token, where its 3b
                # sibling Heinicke2026 does not. DerivedOutputData grammar is
                # per-PUBLICATION, not per-round. Parse from the END regardless.
                stem = (f"{PUBLICATION}_{model}_{gcm}_{FORCING}_{scen}_{SOC}_{SENS}_"
                        f"{VAR}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.nc",
                    url=f"{BASE}/{mdir}/{gcm}/future/{stem}.nc",
                    sidecar=f"{BASE}/{mdir}/{gcm}/future/{stem}.json",
                    model=model, gcm=gcm, scenario=scen, soc=SOC,
                    member=f"{model}_{gcm}",
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
    args = ap.parse_args()

    items = build_items()
    print(f"{len(items)} files = {len(MODELS)} model x {len(GCMS)} GCMs "
          f"x {len(SCENARIOS)} scenarios")
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
    print(f"\n{len(good)}/{len(items)} verified by sha512, {total/2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "sha512_source", "member", "model",
                "gcm", "scenario", "soc"])
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
