"""Ingest ISIMIP2b Zimmer2023 `fldfrc` (CaMa-Flood annual flooded fraction) into S3 raw.

SOURCE: ISIMIP2b/DerivedOutputData/Zimmer2023/{MODEL}/{gcm}/future/ --
    cama-flood_{ghm}_{gcm}_ewembi_{rcp}_2005soc_co2_fldfrc_{protection}_150arcsec_global_annual_2006_2100.nc4

216 files = 6 GHMs x 4 GCMs x rcp26/60/85 x {none, 100yr, flopros}. Three PROTECTION
LEVELS are ingested as three SEPARATE layers, because the protection level is not a
detail -- it changes the sign of the story (see the module notes in
process_fldfrc_flood.py and config/isimip_search_catalog.yaml search_results.flooding).

WHY THIS SCRIPT COARSENS AT INGEST TIME (deliberate deviation, do not "fix")
---------------------------------------------------------------------------
The source is 150 arcsec (4320 x 8640) -- 36x finer than the 0.5 deg TCFD grid -- and the
full 216-file set is ~54 GB, against 19 GB of local scratch. So each member is streamed:
downloaded, checksum-verified against its ISIMIP sidecar, coarsened 12x12 to 0.5 deg,
and the 150 arcsec original deleted. What lands in `TCFD/raw/isimip/{layer_id}/` is the
0.5 deg annual field, ~36x smaller.

That is a departure from "raw is byte-for-byte what ISIMIP served", so it is made
auditable rather than silent: every ingested file records, in its own global attrs AND in
the layer.json manifest, the `source_url` and the **sha512 of the 150 arcsec original**,
plus the exact transform. Raw staging is transient by contract anyway (STORAGE.md:206 --
`cleanup_raw` deletes it once source_url + checksum are recorded), so uploading 54 GB
only to delete it later would be pure churn. files.isimip.org is NOT behind the Anubis
anti-bot that guards the data.isimip.org API, so re-fetching an original is routine.

THE COARSENING IS AREA-PRESERVING AND EXACTLY ALIGNED
-----------------------------------------------------
4320/12 = 360 and 8640/12 = 720 exactly, and each 12-cell block's centre coincides with
the ISIMIP 0.5 deg centre to 1.4e-14 deg, so no interpolation is involved -- it is a
block sum. Sub-cells are weighted by true spherical area (they differ across a block's
12 latitude rows), and the denominator is the FULL 0.5 deg cell area:

    value[cell] = sum(fldfrc_i * area_i) / sum(area_i over ALL 144 sub-cells)
                = flooded area / total cell area

Using only the *valid* sub-cells as the denominator was rejected: it would report
"fraction of the modelled floodplain that flooded" -- inflating cells that are only
partly in the CaMa-Flood domain, and making the field non-aggregatable. With the full-cell
denominator, sum(value * cell_area) over any region is exactly the flooded area in km2,
which is what "area impacted per year" has to mean. Verified area-conserving to 6.3e-10
relative error against the 150 arcsec field.

`floodplain_fraction` (the area share of each 0.5 deg cell that lies inside the
CaMa-Flood domain) is emitted alongside as the coverage diagnostic, so a partly-covered
cell is auditable instead of merely looking low.

COVERAGE IS A NORMAL ISIMIP LAND MASK, NOT A SPARSE SUBSET
----------------------------------------------------------
62,066 cells carry data = 128.8 million km2 = 95.7% of global land excluding Antarctica,
against 67,420 land cells for Lange2020 `led` (61,546 shared). 79.6% of them are >=99%
inside the domain and the median cell is 100%. An earlier note in this project described
the footprint as "~76% of the grid is NaN, i.e. only ~24% of land" -- that compared a
land variable against the whole globe INCLUDING OCEAN and was wrong. Cells outside the
domain stay NaN (no model data); nothing is zero-filled.

Usage:
    python scripts/download_fldfrc_flood.py [--protection none,100yr,flopros]
                                            [--workers 5] [--dry-run] [--limit N]

Resumable: a member already present in S3 raw is skipped unless --force.
"""

import argparse
import hashlib
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402

BASE = "https://files.isimip.org/ISIMIP2b/DerivedOutputData/Zimmer2023"
VAR = "fldfrc"
# 6 GHMs publish FUTURE runs. clm50, mpi-hm and pcr-globwb are historical-only
# (their future/ directories return a genuine 404) -- verified 2026-07-28.
MODELS = {  # directory name -> filename token
    "CLM45": "clm45",
    "CWatM": "cwatm",
    "H08": "h08",
    "LPJmL": "lpjml",
    "MATSIRO": "matsiro",
    "WaterGAP2": "watergap2",
}
GCMS = ["gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5"]
SCENARIOS = ["rcp26", "rcp60", "rcp85"]
PROTECTIONS = ["none", "100yr", "flopros"]
YEARS = (2020, 2099)          # decades 2020s..2090s; source spans 2006-2100
SRC_START_YEAR = 2006
FINE = (4320, 8640)
COARSE = (360, 720)
BLK = 12
R_EARTH = 6371.0


def layer_id(protection):
    return f"riverflood_fldfrc-{protection}_annual"


def src_name(model_tok, gcm, scenario, protection):
    return (f"cama-flood_{model_tok}_{gcm}_ewembi_{scenario}_2005soc_co2_"
            f"{VAR}_{protection}_150arcsec_global_annual_2006_2100.nc4")


def out_name(model_tok, gcm, scenario, protection):
    """0.5 deg product name. Keeps ISIMIP filename grammar so parse-from-the-end still
    works, but says `halfdeg` where the source says `150arcsec` -- the file is NOT the
    original and must not look like it."""
    return (f"cama-flood_{model_tok}_{gcm}_ewembi_{scenario}_2005soc_co2_"
            f"{VAR}_{protection}_halfdeg_global_annual_{YEARS[0]}_{YEARS[1]}.nc")


def src_url(model_dir, gcm, name):
    return f"{BASE}/{model_dir}/{gcm}/future/{name}"


def members():
    for protection in PROTECTIONS:
        for model_dir, tok in MODELS.items():
            for gcm in GCMS:
                for scen in SCENARIOS:
                    yield dict(
                        protection=protection, model_dir=model_dir, model=tok, gcm=gcm,
                        scenario=scen,
                        src=src_name(tok, gcm, scen, protection),
                        out=out_name(tok, gcm, scen, protection),
                        url=src_url(model_dir, gcm, src_name(tok, gcm, scen, protection)),
                    )


def sub_cell_area():
    """True spherical area (km2) of each 150 arcsec sub-cell, by latitude row."""
    d = 1.0 / 24.0
    lat = 90.0 - d / 2 - d * np.arange(FINE[0])
    half = np.deg2rad(d) / 2
    la = np.deg2rad(lat)
    return (R_EARTH ** 2) * np.deg2rad(d) * (np.sin(la + half) - np.sin(la - half))


def _fetch(url, dest, retries=4):
    """Download with retries. Returns bytes written."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as f:
                shutil.copyfileobj(r, f, length=1 << 22)
            return dest.stat().st_size
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
            print(f"    retry {attempt+1} after {type(e).__name__}: {e}", flush=True)
    raise RuntimeError("unreachable")


def sha512_file(path, chunk=1 << 22):
    h = hashlib.sha512()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def process_one(m, cache_dir, force=False):
    """Download -> verify -> coarsen -> write 0.5 deg NetCDF. Returns a result dict.

    Runs in a worker process. Deliberately does NOT touch S3: the parent uploads, so
    there is one credential lifecycle instead of N (STORAGE.md credentials note).
    """
    import cftime
    import netCDF4 as nc

    out_path = Path(cache_dir) / m["out"]
    if out_path.exists() and not force:
        return dict(m, status="cached", out_path=str(out_path))

    tmp = Path(cache_dir) / (m["src"] + ".part")
    src_path = Path(cache_dir) / m["src"]
    try:
        # --- sidecar first: gives us size + sha512 to verify against ------------
        side = {}
        try:
            import json
            with urllib.request.urlopen(m["url"][:-4] + ".json", timeout=120) as r:
                side = json.load(r)
        except Exception:                                          # noqa: BLE001
            pass

        n = _fetch(m["url"], tmp)
        if side.get("size") and n != side["size"]:
            raise ValueError(f"size mismatch: got {n}, sidecar says {side['size']}")
        digest = sha512_file(tmp)
        if side.get("checksum") and side.get("checksum_type") == "sha512" \
                and digest != side["checksum"]:
            raise ValueError("sha512 mismatch against ISIMIP sidecar")
        tmp.rename(src_path)

        # --- coarsen -----------------------------------------------------------
        a_sub = sub_cell_area()                                    # (4320,)
        A_cell = a_sub.reshape(COARSE[0], BLK).sum(1) * BLK        # (360,) km2
        w = np.broadcast_to(a_sub[:, None], FINE)

        ds = nc.Dataset(src_path)
        v = ds.variables[VAR]
        if tuple(v.shape[1:]) != FINE:
            raise ValueError(f"unexpected grid {v.shape[1:]}, expected {FINE}")
        t = ds.variables["time"]
        # Capture the time metadata NOW: `ds` is closed before the output is written,
        # and touching a variable attribute afterwards raises
        # "AttributeError: NetCDF: Attribute not found" from the C library.
        time_units = str(t.units)
        time_cal = str(getattr(t, "calendar", "standard"))
        years = np.array([d.year for d in np.atleast_1d(
            cftime.num2date(t[:], time_units, calendar=time_cal))])
        keep = np.where((years >= YEARS[0]) & (years <= YEARS[1]))[0]
        if len(keep) != YEARS[1] - YEARS[0] + 1:
            raise ValueError(f"expected {YEARS[1]-YEARS[0]+1} years in range, got {len(keep)}")

        vals = np.full((len(keep), *COARSE), np.nan, np.float32)
        fpf = None
        tot_fine = tot_coarse = 0.0
        for k, i in enumerate(keep):
            sl = np.ma.filled(v[i].astype("f4"), np.nan)
            valid = np.isfinite(sl)
            flooded = (np.where(valid, sl, 0.0) * w).reshape(
                COARSE[0], BLK, COARSE[1], BLK).sum(axis=(1, 3))          # km2
            vals[k] = (flooded / A_cell[:, None]).astype(np.float32)
            tot_fine += float((np.where(valid, sl, 0.0) * w).sum())
            tot_coarse += float(flooded.sum())
            if fpf is None:   # domain is static in time; compute once
                fpf = ((np.where(valid, 1.0, 0.0) * w).reshape(
                    COARSE[0], BLK, COARSE[1], BLK).sum(axis=(1, 3))
                    / A_cell[:, None]).astype(np.float32)
        ds.close()

        # Off-domain cells are NaN (no model data), never 0.
        vals[:, fpf <= 0] = np.nan
        rel_err = abs(tot_coarse / tot_fine - 1) if tot_fine > 0 else 0.0
        if rel_err > 1e-6:
            raise ValueError(f"coarsening not area-conserving: rel_err={rel_err:.3e}")

        # --- write 0.5 deg NetCDF ---------------------------------------------
        lat = 90.0 - 0.25 - 0.5 * np.arange(COARSE[0])
        lon = -180.0 + 0.25 + 0.5 * np.arange(COARSE[1])
        yrs = np.arange(YEARS[0], YEARS[1] + 1)
        tmp_out = out_path.with_suffix(".tmp.nc")
        with nc.Dataset(tmp_out, "w") as o:
            o.createDimension("year", len(yrs))
            o.createDimension("lat", COARSE[0])
            o.createDimension("lon", COARSE[1])
            ov = o.createVariable("year", "i4", ("year",)); ov[:] = yrs
            ov = o.createVariable("lat", "f8", ("lat",)); ov[:] = lat
            ov.units, ov.standard_name = "degrees_north", "latitude"
            ov = o.createVariable("lon", "f8", ("lon",)); ov[:] = lon
            ov.units, ov.standard_name = "degrees_east", "longitude"
            ov = o.createVariable(VAR, "f4", ("year", "lat", "lon"),
                                  zlib=True, complevel=5, fill_value=np.float32(np.nan))
            ov[:] = vals
            ov.units = "1"
            ov.long_name = "Annual flooded area fraction of grid cell"
            ov.definition = ("Area-weighted 12x12 block mean of the 150 arcsec CaMa-Flood "
                             "flooded fraction, divided by the FULL 0.5 deg cell area "
                             "(sub-cells outside the CaMa-Flood domain contribute 0 "
                             "flooded area but ARE counted in the denominator).")
            ov = o.createVariable("floodplain_fraction", "f4", ("lat", "lon"),
                                  zlib=True, complevel=5, fill_value=np.float32(np.nan))
            ov[:] = np.where(fpf > 0, fpf, np.nan)
            ov.units = "1"
            ov.long_name = "Area share of grid cell inside the CaMa-Flood domain"
            o.setncatts({
                "title": "ISIMIP2b Zimmer2023 fldfrc coarsened to the 0.5 deg ISIMIP grid",
                "source_url": m["url"],
                "source_sha512": digest,
                "source_size_bytes": int(n),
                "source_grid": "150arcsec (4320 x 8640)",
                "source_time_units": time_units,
                "source_time_calendar": time_cal,
                "transform": ("area-weighted 12x12 block aggregation to 0.5 deg; "
                              "denominator = full cell area; area-conserving to "
                              f"rel_err={rel_err:.2e}; NOT the original ISIMIP file"),
                "years_retained": f"{YEARS[0]}-{YEARS[1]} of 2006-2100",
                "model": m["model"], "gcm": m["gcm"],
                "climate_scenario": m["scenario"], "protection": m["protection"],
                "soc_scenario": "2005soc", "sens_scenario": "co2",
                "created_by": "scripts/download_fldfrc_flood.py",
            })
        tmp_out.rename(out_path)
        return dict(m, status="ok", out_path=str(out_path), source_sha512=digest,
                    source_bytes=int(n), rel_err=rel_err,
                    domain_cells=int((fpf > 0).sum()))
    except Exception as e:                                          # noqa: BLE001
        return dict(m, status="error", error=f"{type(e).__name__}: {e}")
    finally:
        for p in (tmp, src_path):
            try:
                if p.exists():
                    p.unlink()          # the 150 arcsec original never persists
            except OSError:
                pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protection", default=",".join(PROTECTIONS),
                    help="comma-separated subset of none,100yr,flopros")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0, help="process only the first N (smoke test)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    want = [p.strip() for p in args.protection.split(",") if p.strip()]
    bad = set(want) - set(PROTECTIONS)
    if bad:
        raise SystemExit(f"unknown protection level(s): {sorted(bad)}")

    todo = [m for m in members() if m["protection"] in want]
    if args.limit:
        todo = todo[: args.limit]

    print(f"{len(todo)} member-files across protections {want}")
    print(f"target layers: {', '.join(layer_id(p) for p in want)}")
    if args.dry_run:
        for m in todo[:8]:
            print("  ", m["url"])
        print(f"  ... ({len(todo)} total)")
        return

    fs = storage.s3_filesystem()
    cache_dir = storage.cache_root() / "fldfrc_coarsen"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # --- skip members already ingested (resumable) ---------------------------
    present = {}
    for p in want:
        prefix = storage._p(f"{storage.BUCKET}/{storage.raw_prefix(layer_id(p))}")
        present[p] = {os.path.basename(k) for k in fs.glob(f"{prefix}/*.nc")} \
            if fs.exists(prefix) else set()
        print(f"  {layer_id(p)}: {len(present[p])} already in S3 raw")
    if not args.force:
        todo = [m for m in todo if m["out"] not in present[m["protection"]]]
    print(f"{len(todo)} to fetch\n")
    if not todo:
        print("nothing to do")
        return

    done, errors, pending_upload = [], [], {p: [] for p in want}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_one, m, str(cache_dir), args.force): m for m in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            tag = f"{r['model']:<9} {r['gcm']:<13} {r['scenario']} {r['protection']:<7}"
            if r["status"] == "error":
                errors.append(r)
                print(f"[{i}/{len(todo)}] ERROR {tag}  {r['error']}", flush=True)
                continue
            done.append(r)
            pending_upload[r["protection"]].append(r)
            el = time.time() - t0
            print(f"[{i}/{len(todo)}] ok    {tag}  cells={r.get('domain_cells','-')} "
                  f"relerr={r.get('rel_err',0):.1e}  "
                  f"({el/60:.1f} min elapsed, ~{el/i*(len(todo)-i)/60:.0f} min left)",
                  flush=True)

            # upload in batches of 12 so a long run does not hold everything on disk
            batch = pending_upload[r["protection"]]
            if len(batch) >= 12:
                _upload(batch, r["protection"])
                pending_upload[r["protection"]] = []

    for p, batch in pending_upload.items():
        if batch:
            _upload(batch, p)

    print(f"\n{len(done)} ingested, {len(errors)} errors in {(time.time()-t0)/60:.1f} min")
    for e in errors:
        print(f"  FAILED {e['model']} {e['gcm']} {e['scenario']} {e['protection']}: {e['error']}")
    if errors:
        raise SystemExit(1)


def _upload(batch, protection):
    lid = layer_id(protection)
    paths = [Path(r["out_path"]) for r in batch]
    urls = {Path(r["out_path"]).name: r["url"] for r in batch}
    storage.ingest_raw(paths, lid, source_urls=urls)
    print(f"    -> uploaded {len(paths)} files to raw/isimip/{lid}/", flush=True)
    for p in paths:
        try:
            p.unlink()
        except OSError:
            pass


if __name__ == "__main__":
    main()
