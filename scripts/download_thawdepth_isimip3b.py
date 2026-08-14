"""Download the ISIMIP3b `thawdepth` ensemble for a permafrost-thaw layer.

`thawdepth` is the ONLY variable in the ISIMIP repository that names this hazard. There is
no permafrost-extent, ground-ice, subsidence or solifluction variable in any round, and no
derived exposure index: the Lange2020 six-hazard family (drought, river flood, wildfire,
crop failure, heatwave, tropical cyclone) has no permafrost member and neither Heinicke2026
nor Zantout2025 re-issues one. So this layer starts from a PHYSICAL FIELD, not from a
ready-made "area exposed" index -- the opposite starting point from every exposure layer
shipped so far, and the reason the framing decisions below are open rather than inherited.

Ensemble: 3 impact models x {5,5,2} CMIP6 GCMs x 3 SSPs = 36 files, ~190 MB.

Coverage ENUMERATED 2026-08-14 (full matrices in config/isimip_search_catalog.yaml
search_results.permafrost_thaw). Serial listings, empty listing treated as failure, variable
field projected from the END at $(NF-4) rather than grepped for a believed-in token
(GUARDRAILS 8).

  * COMPLETE matrix -- every (model, gcm, scenario) combination exists, and the member
    composition is IDENTICAL across ssp126/ssp370/ssp585, so the shared 2020s baseline rests
    on a uniform 12-member ensemble.
  * annual, w5e5, 2015-2100. The 2020s baseline decade is complete (3b future runs start in
    2015, so there is no complete 2010s decade).
  * soc/sens are uniform `2015soc-from-histsoc` / `default` ACROSS THE FILES DOWNLOADED HERE,
    but that uniformity is a CHOICE and it is not the repository's usual preference. See
    "WHY 2015soc-from-histsoc" below.

THE SECTOR TRAP -- why this ensemble is 3 models and not 1.

`thawdepth` is not confined to the `permafrost` sector. ISIMIP3b `permafrost` holds two
models and only ONE of them (LPJmL5-7-10-fire) publishes thawdepth at all; JULES-ES-VN6P3
and CLASSIC publish it under `biomes` and are absent from the permafrost sector entirely.
A sector-only walk answers "1 model, 5 members" to a question whose true answer is "3
models, 12 members" (ELM-ECA, the permafrost sector's other model, publishes `tsl` and `snd`
but no thawdepth).

The files are republished BYTE-IDENTICALLY across `permafrost`, `biomes` and `water_global`
-- verified 2026-08-14, 3b LPJmL5 ssp370 matched on Content-Length (6,697,198) AND ETag
("6899d4e4-6630ee") between permafrost and biomes. Each URL below is taken from the ONE
directory that was actually listed for that model, and every member appears exactly once;
ingesting the same model from two sectors would double-weight it.

WHY 2015soc-from-histsoc, and what it costs.

The three models do not offer the same soc tokens:

    LPJmL5-7-10-fire   1850soc | 2015soc | 2015soc-from-histsoc
    CLASSIC                      2015soc | 2015soc-from-histsoc
    JULES-ES-VN6P3                         2015soc-from-histsoc

`2015soc-from-histsoc` is the ONLY combination all three share, and it is what makes this a
3-model ensemble. Choosing the repository's usual `2015soc` preference instead would drop
JULES-ES-VN6P3 entirely: 7 members, 2 models. That trade is recorded here rather than
buried -- if the user prefers the stricter socioeconomic isolation, re-run with SOC set to
`2015soc` and expect a thinner, differently-composed ensemble.

NO SIDECARS. `{stem}.json` 404s for this product in both the permafrost and biomes sectors
(verified 2026-08-14 on two LPJmL5 members), unlike Zantout2025 where the sidecar exists
under exactly that name. Downloads here are therefore verified against the server's
Content-Length only, and the sha512 recorded in the provenance CSV is COMPUTED LOCALLY --
it proves the file did not change after we wrote it, NOT that it matches what the publisher
intended. Do not describe it as publisher-verified.

DATA NATURE IS UNVERIFIED AT DOWNLOAD TIME and must be value-checked before any processor is
written (GUARDRAILS 9). The specific unknowns, none of which a listing can answer:

  * WHAT A NON-PERMAFROST CELL CONTAINS -- NaN, the full soil-column depth, or 0. This
    single fact decides the mask, the statistic branch, and whether the field saturates.
    Outside the permafrost domain a thaw depth is either meaningless or degenerate, and no
    permafrost mask ships with the product.
  * WHETHER THE ANNUAL VALUE IS A MAXIMUM OR A MEAN. An annual thaw depth is presumably the
    maximum, i.e. active-layer thickness, but the filename does not say so and mean-vs-max
    changes what the number means for ground stability.
  * units and long_name, which cannot be read remotely (NetCDF4 is HDF5 and places attribute
    headers at unpredictable offsets; a 4 MB prefix read of a sibling product surfaced no
    units at all).

Run `scripts/check_thawdepth_nature.py` after this script and before writing a processor.

Usage:
    python scripts/download_thawdepth_isimip3b.py [--dry-run]
"""

import argparse
import csv
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/OutputData"
LAYER_ID = "permafrost-isimip3b_thawdepth_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "thawdepth"
FORCING = "w5e5"
SOC = "2015soc-from-histsoc"
SENS = "default"
YEARS = "2015_2100"
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

#: (sector, model directory, filename model token, GCMs). The sector is the directory that
#: was actually LISTED for that model on 2026-08-14 -- not a guess, and not interchangeable
#: with the byte-identical republications elsewhere. CLASSIC runs only 2 GCMs in ISIMIP3b;
#: that is the model's coverage, not a missing download.
MEMBERS = [
    ("permafrost", "LPJmL5-7-10-fire", "lpjml5-7-10-fire",
     ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]),
    ("biomes", "JULES-ES-VN6P3", "jules-es-vn6p3",
     ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]),
    ("biomes", "CLASSIC", "classic",
     ["gfdl-esm4", "ukesm1-0-ll"]),
]

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b thawdepth ingest)"
CHUNK = 1 << 20


def build_items():
    items = []
    for sector, mdir, model, gcms in MEMBERS:
        for gcm in gcms:
            for scen in SCENARIOS:
                stem = (f"{model}_{gcm}_{FORCING}_{scen}_{SOC}_{SENS}_"
                        f"{VAR}_global_annual_{YEARS}")
                items.append(dict(
                    fname=f"{stem}.nc",
                    url=f"{BASE}/{sector}/{mdir}/{gcm}/future/{stem}.nc",
                    sector=sector, model=mdir, gcm=gcm, scenario=scen, soc=SOC,
                    member=f"{mdir}_{gcm}",
                ))
    return items


def head_size(url):
    """Content-Length from a HEAD. The only integrity reference this product offers."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        n = r.headers.get("Content-Length")
        if n is None:
            raise RuntimeError("no Content-Length on HEAD")
        return int(n)


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
        size = head_size(item["url"])
    except urllib.error.HTTPError as e:
        return dict(item, ok=False, msg=f"HEAD HTTP {e.code}")
    except Exception as e:  # noqa: BLE001 - report, never abort the batch
        return dict(item, ok=False, msg=f"HEAD {type(e).__name__}: {e}")

    try:
        wrote = fetch(item["url"], dest, size)
    except Exception as e:  # noqa: BLE001
        return dict(item, ok=False, msg=f"GET {type(e).__name__}: {e}")

    got = dest.stat().st_size if dest.exists() else 0
    if got != size:
        return dict(item, ok=False, msg=f"size {got} != {size}")

    return dict(item, ok=True, bytes=size, sha512=sha512_of(dest),
                wrote=wrote, skipped=(wrote == 0), msg="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    items = build_items()
    n_members = sum(len(g) for _, _, _, g in MEMBERS)
    print(f"{len(items)} files = {n_members} members x {len(SCENARIOS)} scenarios "
          f"({len(MEMBERS)} models, soc={SOC}, sens={SENS})")
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
    print(f"\n{len(good)}/{len(items)} verified by Content-Length, {total/2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "sha512_source", "member", "model",
                "gcm", "scenario", "soc", "sector"])
            w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                row = {k: r[k] for k in w.fieldnames if k != "sha512_source"}
                # NOT a publisher checksum. This product ships no sidecar, so the digest
                # proves only that the bytes on disk are the bytes we wrote.
                row["sha512_source"] = "computed locally (no publisher sidecar exists)"
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
