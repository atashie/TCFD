"""Download temperate needleleaf evergreen NPP (ISIMIP2b biomes) for the conifer layer.

The PFT class differs per model -- there is no shared code, so each model's own class is
named explicitly rather than pattern-matched:

    CLM45     npp-needleleaf-evergreen-tree-temperate
    ORCHIDEE  npp-tendev
    LPJmL     npp-temperate-needleleaved-evergreen-tree

All three are genuinely TEMPERATE classes (unlike CLASSIC `evgndltr` and JULES `ndlevg`,
which merge boreal with temperate, and unlike MC2-USFS's needleleaf biome names, which
publish `pft` presence flags only). Enumerated 2026-08-12; see
`config/isimip_search_catalog.yaml` -> pft_equivalences.temperate_needleleaf_evergreen.

EXPERIMENT: `2005soc` + `co2` (transient CO2) for every member, so the ensemble is
experimentally uniform and no CO2-treatment mixing is needed (contrast `csoil`, where
JULES publishes only its fixed-CO2 run).

COVERAGE IS NOT UNIFORM ACROSS MODELS -- measured, not assumed:

    model     rcp26            rcp60            rcp85
    ORCHIDEE  4 GCMs           4 GCMs           4 GCMs
    LPJmL     4 GCMs           4 GCMs           4 GCMs
    CLM45     hadgem2, miroc5  ipsl only        gfdl, hadgem2

CLM45 publishes only 7 files for this PFT across ALL future scenarios and **no GCM of its
appears in all three RCPs**, so including it makes ensemble composition scenario-dependent
and the shared 2020s baseline can no longer be bit-identical across scenarios (OUTPUT-SPEC
allows this, but it must be declared via `members_by_scenario`). `--models` selects.

`pft-` cover fractions are fetched alongside for ORCHIDEE and LPJmL because the two models
report NPP on different denominators -- ORCHIDEE on the PFT's own tile area, LPJmL
cover-scaled (measured previously on `timber_*-tempnle`: raw spread 10.5x/177x vs 2.35x/
1.83x harmonized). CLM45 publishes NO `pft-` fraction at all, which is why its
denominator cannot be harmonized from within the model.

Raw goes to data/raw/conifer-temperate_npp-tempnle_annual/.

Usage:
    python scripts/download_tempnle_npp.py [--models orchidee,lpjml,clm45] [--no-pft] [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP2b/OutputData/biomes"
LAYER_ID = "conifer-temperate_npp-tempnle_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

GCMS = ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
SCENARIOS = ["rcp26", "rcp60", "rcp85"]
SOC = "2005soc"
SENS = "co2"
YEARS = "2006_2099"

#: model -> (directory name, filename prefix, its own temperate-NLE class code)
MODELS = {
    "orchidee": ("ORCHIDEE", "orchidee", "tendev"),
    "lpjml": ("LPJmL", "lpjml", "temperate-needleleaved-evergreen-tree"),
    "clm45": ("CLM45", "clm45", "needleleaf-evergreen-tree-temperate"),
}

#: Measured 2026-08-12 -- the ONLY (scenario, gcm) pairs CLM45 publishes at 2005soc/co2.
CLM45_AVAILABLE = {
    ("rcp26", "hadgem2-es"), ("rcp26", "miroc5"),
    ("rcp60", "ipsl-cm5a-lr"),
    ("rcp85", "gfdl-esm2m"), ("rcp85", "hadgem2-es"),
}

#: CLM45 publishes no pft- cover fraction for any PFT.
NO_PFT_FRACTION = {"clm45"}

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP2b tempnle npp ingest)"
CHUNK = 1 << 20


def build_items(models, want_pft):
    items = []
    for key in models:
        dirname, prefix, pft = MODELS[key]
        metrics = ["npp"] + ([] if (not want_pft or key in NO_PFT_FRACTION) else ["pft"])
        for gcm in GCMS:
            for scen in SCENARIOS:
                if key == "clm45" and (scen, gcm) not in CLM45_AVAILABLE:
                    continue
                for metric in metrics:
                    stem = (f"{prefix}_{gcm}_ewembi_{scen}_{SOC}_{SENS}_"
                            f"{metric}-{pft}_global_annual_{YEARS}")
                    items.append(dict(
                        fname=f"{stem}.nc4",
                        url=f"{BASE}/{dirname}/{gcm}/future/{stem}.nc4",
                        sidecar=f"{BASE}/{dirname}/{gcm}/future/{stem}.json",
                        model=prefix, gcm=gcm, scenario=scen, metric=metric, pft=pft,
                        member=f"{prefix}_{gcm}",
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
    have = dest.stat().st_size if dest.exists() else 0
    if have == expected_size:
        return 0
    if have > expected_size:
        dest.unlink()
        have = 0
    headers = {"User-Agent": USER_AGENT}
    mode = "wb"
    if have:
        headers["Range"] = f"bytes={have}-"
        mode = "ab"
    req = urllib.request.Request(url, headers=headers)
    written = 0
    with urllib.request.urlopen(req, timeout=900) as r, open(dest, mode) as fh:
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
    except Exception as e:  # noqa: BLE001
        return dict(item, ok=False, msg=f"sidecar {type(e).__name__}: {e}")

    size = int(meta["size"])
    want = meta.get("checksum", "")

    # The sidecar's parsed specifiers are the authority on what the file IS.
    spec = meta.get("specifiers", {})
    for key, expect in (("climate_scenario", item["scenario"]),
                        ("soc_scenario", SOC), ("sens_scenario", SENS)):
        if spec.get(key) not in (None, expect):
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
    return dict(item, ok=True, bytes=size, sha512=digest, wrote=wrote,
                skipped=(wrote == 0), msg="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="orchidee,lpjml",
                    help="comma-separated subset of orchidee,lpjml,clm45")
    ap.add_argument("--no-pft", action="store_true",
                    help="skip the pft- cover fractions (harmonization inputs)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    models = [m.strip().lower() for m in args.models.split(",") if m.strip()]
    bad = [m for m in models if m not in MODELS]
    if bad:
        print(f"unknown model(s): {bad}; choose from {list(MODELS)}")
        return 2

    items = build_items(models, want_pft=not args.no_pft)
    npp = sum(1 for i in items if i["metric"] == "npp")
    print(f"{len(items)} files ({npp} npp + {len(items) - npp} pft) for models={models}")
    if args.dry_run:
        for it in items:
            print(it["url"])
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    for i, it in enumerate(items, 1):  # serial: parallel requests get rate-limited
        r = one(it)
        results.append(r)
        tag = "skip" if r.get("skipped") else ("ok" if r["ok"] else "FAIL")
        print(f"[{i}/{len(items)}] {tag:4s} {r['fname'][:88]}"
              + ("" if r["ok"] else f"  <-- {r['msg']}"), flush=True)

    good = [r for r in results if r["ok"]]
    bad_r = [r for r in results if not r["ok"]]
    print(f"\n{len(good)}/{len(items)} verified by sha512, "
          f"{sum(r['bytes'] for r in good) / 2**20:.1f} MiB")

    if good:
        out = OUT_DIR / "download_provenance.csv"
        exists = out.exists()
        with open(out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "member", "model", "gcm",
                "scenario", "metric", "pft"])
            if not exists:
                w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                w.writerow({k: r[k] for k in w.fieldnames})
        print(f"provenance -> {out}")

    if bad_r:
        print(f"\n{len(bad_r)} FAILED (re-run to resume):")
        for r in bad_r:
            print(f"  {r['fname']}: {r['msg']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
