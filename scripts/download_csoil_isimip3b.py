"""Download the ISIMIP3b `csoil-total` ensemble for a soil-organic-carbon layer.

Soil organic carbon stock (kg C m-2) is the direct subsurface storage pool -- a STOCK, not
the net-sink FLUX `nbp`, and distinct from the vegetation pools `cveg`/`croot`/`cvegbg`. It
is the entry-level indicator for the `soil-degradation` hazard family: more soil carbon =
healthier soil, less = degrading soil. It addresses exactly ONE of the ten degradation
processes enumerated in recital 4 of Directive (EU) 2025/2360 -- loss of soil organic carbon
-- and the layer must be named for that, not for the family. See
config/hazard_taxonomy.yaml families.soil-degradation.

Ensemble: 4 impact models x {2,5,5,5} CMIP6 GCMs x 3 SSPs = 51 files, ~932 MB, 17 members
per scenario. Coverage ENUMERATED 2026-08-14 (matrices in
config/isimip_search_catalog.yaml search_results.soil_carbon).

THIS SUPERSEDES A 12-MEMBER, 3-MODEL COUNT THAT WAS WRONG BY OMISSION.

The 2026-07-25 enumeration recorded `csoil_annual: [classic, jules-es-vn6p3, mc2-usfs]`.
LPJmL5-7-10-fire publishes `csoil-total` annually for all three SSPs and was missed --
+5 members, +1 structurally independent model. The failure mode is the one this repository
keeps logging (GUARDRAILS 11): an UNDERSTATED POSITIVE in our own catalog reads as a
measured fact and is never re-checked, where a recorded negative at least looks like
something to test. The skill's own notes already warned that the catalog "later still missed
CLM45 and VEGAS for csoil" in ISIMIP2b; the same under-listing had happened again in 3b.

Method, so the count can be trusted this time: the variable field was projected from the END
at $(NF-4) over a full listing of one GCM per model -- never grepped for a believed-in token
(GUARDRAILS 8) -- and then EVERY (model, gcm, scenario) URL below was confirmed by HTTP
status. Each model's token was proved valid by requiring its gfdl-esm4 file to return 200
before any 404 elsewhere was believed, so CLASSIC's three 404s are the model's real coverage
(it runs 2 GCMs in ISIMIP3b) and not a filename error.

SOC/SENS ARE NOT UNIFORM, AND CANNOT BE MADE SO.

    model             soc                   sens      why not something else
    classic           2015soc-from-histsoc  default   also offers 2015soc/default
    jules-es-vn6p3    2015soc-from-histsoc  2015co2   publishes NO transient-CO2 csoil run
    lpjml5-7-10-fire  2015soc-from-histsoc  default   also offers 1850soc, 2015soc, nat
    mc2-usfs-r87g5c1  nat                   default   publishes ONLY nat

No combination is shared by all four. Demanding uniform `sens=default` drops JULES entirely
(12 members, 3 models); demanding uniform `soc` drops MC2 as well (12 members, 3 models).
The heterogeneity is therefore RETAINED and DECLARED, consistent with the 2026-07-25 user
decision on this same layer and with the standing preference for ensemble abundance.

What the heterogeneity costs, stated so it reaches the customer rather than dying here:

  * CO2: 12 of 17 members run transient CO2; JULES's 5 run fixed 2015 CO2. Adding LPJmL
    improved this ratio from 7/12 to 12/17 -- it takes `default`, the same tokens as
    CLASSIC, so it introduces NO new heterogeneity dimension and is purely additive.
    DO NOT CALL THE FIXED-CO2 TREND "MUTED". Measured 2026-08-15 on ssp585, those members
    show the LARGEST relative soil-carbon loss of the four models (-4.37%, against -2.75%
    lpjml, -0.05% mc2-usfs, +0.79% classic), because removing fertilisation removes litter
    input -- so the mixed treatment makes the ensemble decline somewhat STRONGER than a
    uniformly transient one would, not weaker. The direction a CO2 treatment biases a trend
    is a measurement, never a deduction from the mechanism's name (GUARDRAILS 9).
  * LAND USE: `nat` (MC2) is a natural-vegetation run with no land-use forcing at all, while
    `2015soc-from-histsoc` holds land use fixed at 2015 after a historical transient. NEITHER
    lets the layer see management-driven soil-carbon loss -- which is the point for a
    CLIMATE-hazard framing (it isolates the climate signal) and a hard limit on any claim
    about soil degradation in general, most of which is management.

The processor MUST write this into `interpretation_caveat`, which is on the
`LAYER_ATTRS_EXPORTED` allowlist in scripts/utils/delivery.py. The previous run put it in
`co2_treatment`, which is NOT on that closed allowlist and would have been silently dropped
before reaching layers.csv, the caveat generator or either report.

ONE SECTOR, WALKED EVERYWHERE. `csoil-total` is republished byte-identically across sectors
(`elm-eca` in biomes+permafrost, same sha512). Every member here is taken from `biomes`, the
one directory actually listed, so no model is double-weighted. The other sectors were probed
for ROSTER, not ingested: 3b `permafrost` holds only ELM-ECA and LPJmL5-7-10-fire and adds
no model biomes lacks, and `fire` is a subset of biomes. That probe is mandatory rather than
optional -- on `thawdepth` a sector-only walk answered "1 model, 5 members" against a true
"3 models, 12 members".

MODELS DELIBERATELY EXCLUDED, with the reason:

  * ELM-ECA  -- publishes csoil-total MONTHLY only (verified: its biomes vocabulary is
    entirely monthly). Annualising a smooth stock is cheap and would add 5 members, but it
    is a GUARDRAILS 1-2 aggregation decision that has not been taken. Also effectively
    ~4x5 degrees despite a 0.5-degree header, per the search skill.
  * VISIT    -- publishes csoil-total MONTHLY only. Same open decision.

Both are recorded rather than dropped silently: an unshipped member with a stated reason is
re-checkable, an absent one is invisible.

SIDECARS EXIST, at `{stem}.json` -- NOT `{stem}.nc.json`, which 404s and reads convincingly
as "this publication has no sidecars". The cropfailure ingest made exactly that error and ran
unverified until it was caught. They carry `checksum` (sha512), `size`, parsed `specifiers`
and a full `netcdf_header`, so every file here is verified against a PUBLISHER checksum, not
merely against Content-Length.

DATA NATURE IS DECLARED, NOT MEASURED, at download time (GUARDRAILS 9). The sidecar header
declares units `kg m-2` and long_name "Carbon Mass in Soil Pool" for LPJmL. What a listing
cannot answer, and what must be value-checked before the processor is trusted:

  * WHAT A NON-SOIL CELL CONTAINS -- NaN, 0, or a fill value. LPJmL's own global attributes
    say "pools are zero for cell fractions covered by inland waterbodies" and give
    reference_area as "continental area (including inland water bodies)", so a zero may be a
    water fraction rather than carbon-free soil. This decides the mask.
  * WHETHER THE FOUR MODELS' MAGNITUDES REMAIN COMPARABLE with LPJmL added. The 2026-07-25
    check found medians 5.76 / 10.28 / 7.67 across the three original models and concluded
    no normalization was needed; LPJmL has not been checked and could break that.
  * WHETHER THE POOLED SAMPLE IS UNIMODAL. With four models at different levels the
    `pooled_mean_multimodel` branch (OUTPUT-SPEC.md) must be tested against, not assumed away.

Run `scripts/check_csoil_nature.py` after this script and before writing the processor.

Usage:
    python scripts/download_csoil_isimip3b.py [--dry-run]
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://files.isimip.org/ISIMIP3b/OutputData"
SECTOR = "biomes"
LAYER_ID = "soilcarbon_csoil_annual"
OUT_DIR = Path("data/raw") / LAYER_ID

VAR = "csoil-total"
FORCING = "w5e5"
YEARS = "2015_2100"
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

FIVE_GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]

#: (model directory, filename token, soc, sens, GCMs). soc/sens are PER MODEL -- see the
#: module docstring. Every combination below returned HTTP 200 on 2026-08-14; CLASSIC's
#: 2-GCM roster is its real ISIMIP3b coverage, confirmed by genuine 404s on the other three
#: after its token was proved by a 200 on gfdl-esm4.
MEMBERS = [
    ("CLASSIC", "classic", "2015soc-from-histsoc", "default",
     ["gfdl-esm4", "ukesm1-0-ll"]),
    ("JULES-ES-VN6P3", "jules-es-vn6p3", "2015soc-from-histsoc", "2015co2", FIVE_GCMS),
    ("LPJmL5-7-10-fire", "lpjml5-7-10-fire", "2015soc-from-histsoc", "default", FIVE_GCMS),
    ("MC2-USFS-r87g5c1", "mc2-usfs-r87g5c1", "nat", "default", FIVE_GCMS),
]

USER_AGENT = "TCFD-pipeline/1.0 (ISIMIP3b csoil-total ingest)"
CHUNK = 1 << 20


def build_items():
    items = []
    for mdir, token, soc, sens, gcms in MEMBERS:
        for gcm in gcms:
            for scen in SCENARIOS:
                stem = (f"{token}_{gcm}_{FORCING}_{scen}_{soc}_{sens}_"
                        f"{VAR}_global_annual_{YEARS}")
                url = f"{BASE}/{SECTOR}/{mdir}/{gcm}/future/{stem}"
                items.append(dict(
                    fname=f"{stem}.nc",
                    url=f"{url}.nc",
                    sidecar_url=f"{url}.json",
                    sector=SECTOR, model=mdir, gcm=gcm, scenario=scen,
                    soc=soc, sens=sens, member=f"{mdir}_{gcm}",
                ))
    return items


def sidecar(url):
    """Publisher size + sha512. `{stem}.json`, NOT `{stem}.nc.json` (that 404s)."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        d = json.load(r)
    if d.get("checksum_type") != "sha512":
        raise RuntimeError(f"unexpected checksum_type {d.get('checksum_type')!r}")
    return int(d["size"]), d["checksum"]


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
        size, want_sha = sidecar(item["sidecar_url"])
    except urllib.error.HTTPError as e:
        return dict(item, ok=False, msg=f"sidecar HTTP {e.code}")
    except Exception as e:  # noqa: BLE001 - report, never abort the batch
        return dict(item, ok=False, msg=f"sidecar {type(e).__name__}: {e}")

    # A file already on disk at the right size is re-hashed, not trusted: resume can leave a
    # correct length with wrong bytes if the server served a different revision mid-range.
    try:
        wrote = fetch(item["url"], dest, size)
    except Exception as e:  # noqa: BLE001
        return dict(item, ok=False, msg=f"GET {type(e).__name__}: {e}")

    got = dest.stat().st_size if dest.exists() else 0
    if got != size:
        return dict(item, ok=False, msg=f"size {got} != {size}")

    have_sha = sha512_of(dest)
    if have_sha != want_sha:
        return dict(item, ok=False, msg="sha512 MISMATCH vs publisher sidecar")

    return dict(item, ok=True, bytes=size, sha512=have_sha,
                wrote=wrote, skipped=(wrote == 0), msg="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    items = build_items()
    n_members = sum(len(g) for *_, g in MEMBERS)
    print(f"{len(items)} files = {n_members} members x {len(SCENARIOS)} scenarios "
          f"({len(MEMBERS)} models, biomes sector)")
    for mdir, _, soc, sens, gcms in MEMBERS:
        print(f"  {mdir:<18} {len(gcms)} GCMs  soc={soc:<21} sens={sens}")
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
    print(f"\n{len(good)}/{len(items)} verified against publisher sha512, "
          f"{total/2**20:.1f} MiB")

    if good:
        with open(OUT_DIR / "download_provenance.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "fname", "url", "bytes", "sha512", "sha512_source", "member", "model",
                "gcm", "scenario", "soc", "sens", "sector"])
            w.writeheader()
            for r in sorted(good, key=lambda x: x["fname"]):
                row = {k: r[k] for k in w.fieldnames if k != "sha512_source"}
                # Publisher-verified, unlike thawdepth: this product DOES ship sidecars.
                row["sha512_source"] = "publisher sidecar {stem}.json (checksum_type=sha512)"
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
