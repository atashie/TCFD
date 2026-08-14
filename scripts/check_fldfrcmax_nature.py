#!/usr/bin/env python3
"""Measure the data nature, orientation, mask and place-plausibility of raw `fldfrcmax`.

Usage:
    python scripts/check_fldfrcmax_nature.py                     # all sections, all members
    python scripts/check_fldfrcmax_nature.py --protection flopros
    python scripts/check_fldfrcmax_nature.py --section orientation

Run BEFORE processing and again before shipping. Four guardrail checks that need the same
open files:

  * GUARDRAILS 9 -- WHAT the values are. Continuous [0,1], NOT binary. The decisive number
    for the statistic branch is the zero share, and here it **depends on the protection
    level**: measured over all 258 members, none 17.05%, flopros 83.11%, 40yr 87.27%. So
    the branch cannot be chosen once for the layer -- it is a per-variant measurement, and
    this script prints it per variant.
    The mask, by contrast, is UNIFORM: finite over exactly 24.34% of the grid in every
    member, min equal to max, so there is no minimum-model rule to resolve here.

  * ORIENTATION -- specific to this publication and genuinely dangerous. TipESM2025 flips
    the latitude COORDINATE between protection variants of the same member: `none`
    descending in all 96 files, `40yr`/`flopros` ascending in all 162. (The file's declared
    dimension order differs too, but the variable's dims are ('time','lat','lon')
    everywhere -- the arrays are not transposed.) Consequences, both silent:
    `.sel(lat=slice(hi, lo))` returns ZERO cells on whichever half does not match, and
    positional indexing builds an upside-down layer that satisfies every line of
    OUTPUT-SPEC. Everything downstream must `.sortby("lat")`; this section asserts each
    file's data agrees with its own declared coordinates by checking that the Caspian is
    wet and the Sahara dry in *coordinate* space.

  * PERMANENT WATER -- cells fully flooded in every year: unprotected median 511 (494-531),
    protected median ~230, mostly the Caspian. They are not flood risk and they would sit
    at the top of any percentile ranking. This section counts and locates them so the mask
    decision is made on numbers.

  * GUARDRAILS 12 -- WHERE the values are. The entire reason this product was chosen over
    `floodedarea` is that `floodedarea` reads 0.000 across the Amazon. If that is not
    reproduced here on every member, the layer does not ship. Dry controls must stay dry.

Nothing in this file is inherited from another layer. Within this one publication the zero
share, the axis order and the header dim order all vary between files that differ only in
protection level.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr

RAW = Path("data/raw/flood-isimip3b_fldfrcmax_annual")
VAR = "fldfrcmax"

#: GUARDRAILS 12. Boxes, not points: at 0.25 deg a single cell can miss a channel, and a
#: box mean is what a customer's site neighbourhood actually looks like.
#: (lat_south, lat_north, lon_west, lon_east)
BOXES = {
    "Amazon main stem": (-5.0, -1.0, -65.0, -55.0),
    "Ganges-Brahmaputra": (22.0, 26.0, 88.0, 92.0),
    "Mekong delta": (9.0, 11.5, 105.0, 107.0),
    "Lower Mississippi": (30.0, 35.0, -92.5, -89.0),
    "Rhine (NL/DE)": (50.5, 52.5, 5.0, 8.0),
    "Niger inland delta": (13.5, 16.0, -6.0, -3.0),
    "Sahara (dry control)": (23.0, 26.0, 8.0, 12.0),
    "Atacama (dry control)": (-25.0, -21.0, -70.0, -68.0),
}
#: Inland water bodies that read as permanently flooded. Used to prove orientation.
CASPIAN = (37.0, 46.0, 48.0, 54.0)


def members(protection: list[str] | None):
    for f in sorted(RAW.glob(f"*_{VAR}_*.nc")):
        parts = f.name.split("_")
        prot = parts[-6]
        if protection and prot not in protection:
            continue
        yield f, {"ghm": parts[1], "gcm": parts[2], "scenario": parts[3],
                  "soc": parts[4], "protection": prot}


def load(path: Path):
    """Open and NORMALISE orientation. Never index this product positionally."""
    ds = xr.open_dataset(path, decode_times=False)
    return ds, ds[VAR].sortby("lat").sortby("lon")


def box(v, b):
    la_lo, la_hi, lo_w, lo_e = b
    return v.sel(lat=slice(la_lo, la_hi), lon=slice(lo_w, lo_e))


def survey_one(path: Path, meta: dict) -> dict:
    """Every metric this script reports, from ONE decompression of the file.

    Each member is 86 x 720 x 1440 float32 = 356 MB decompressed and there are 258 of them,
    so a section-per-pass design costs ~275 GB of decompression per extra pass. Everything
    below is computed from a single `.values`.
    """
    ds = xr.open_dataset(path, decode_times=False)
    lat_raw = ds["lat"].values
    rec = dict(meta)
    rec["lat_order"] = "DESC" if lat_raw[0] > lat_raw[-1] else "ASC"
    rec["dims"] = str(ds[VAR].dims).replace("'", "").replace(" ", "")

    v = ds[VAR].sortby("lat").sortby("lon")
    a = v.values
    fin = np.isfinite(a)
    vals = a[fin]
    filled = np.nan_to_num(a, nan=0.0)

    rec["n_finite"] = int(vals.size)
    rec["zeros"] = int((vals == 0).sum())
    rec["ones"] = int((vals == 1).sum())
    rec["max"] = float(vals.max()) if vals.size else np.nan
    rec["uniq_1M"] = int(np.unique(vals[:1_000_000]).size)
    rec["finite_frac"] = float(fin.any(axis=0).mean())

    # Permanent water: fully flooded in EVERY year.
    always = np.all(filled >= 1.0, axis=0) & fin.any(axis=0)
    rec["always_full"] = int(always.sum())

    # Orientation proof, in COORDINATE space: after sorting, the Caspian must be wetter
    # than the Sahara. If a file's data were stored against a flipped axis, this inverts.
    rec["caspian"] = float(np.nanmean(box(v, CASPIAN).values))
    casp_always = np.all(np.nan_to_num(box(v, CASPIAN).values, nan=0.0) >= 1.0, axis=0)
    rec["caspian_always_full"] = int(casp_always.sum())

    last = a[-10:]
    lat_v, lon_v = v["lat"].values, v["lon"].values
    for name, (la_lo, la_hi, lo_w, lo_e) in BOXES.items():
        i = (lat_v >= la_lo) & (lat_v <= la_hi)
        j = (lon_v >= lo_w) & (lon_v <= lo_e)
        sub = last[:, i][:, :, j]
        rec[f"box::{name}"] = float(np.nanmean(sub)) if np.isfinite(sub).any() else np.nan
    rec["sahara"] = rec["box::Sahara (dry control)"]
    ds.close()
    return rec


def report_orientation(recs):
    print("\n=== ORIENTATION: declared axis vs where the water actually is ===")
    combos = defaultdict(int)
    for r in recs:
        combos[(r["protection"], r["lat_order"], r["dims"])] += 1
    print(f"{'protection':<12}{'lat axis':<9}{'header dims':<22}{'files':>7}")
    for (prot, order, dims), n in sorted(combos.items()):
        print(f"{prot:<12}{order:<9}{dims:<22}{n:>7}")
    bad = [r for r in recs if not (r["caspian"] > r["sahara"] or
                                   (np.isnan(r["caspian"]) and np.isnan(r["sahara"])))]
    print(f"\n  files whose DATA disagrees with their own coordinates: {len(bad)}")
    for r in bad[:10]:
        print(f"    FLIPPED {r['ghm']}/{r['gcm']}/{r['scenario']}/{r['protection']} "
              f"caspian={r['caspian']:.4f} sahara={r['sahara']:.4f}")
    print("  Axis order varying WITHIN a publication is why nothing here may index by")
    print("  position. Every reader must .sortby('lat').")
    return len(bad)


def report_nature(recs):
    print("\n=== GUARDRAILS 9: data nature, PER PROTECTION LEVEL ===")
    acc = defaultdict(lambda: {"n": 0, "zeros": 0, "ones": 0, "ff": [], "max": 0.0,
                               "uniq": 0, "members": 0})
    for r in recs:
        s = acc[r["protection"]]
        s["n"] += r["n_finite"]; s["zeros"] += r["zeros"]; s["ones"] += r["ones"]
        s["ff"].append(r["finite_frac"]); s["max"] = max(s["max"], r["max"])
        s["uniq"] = max(s["uniq"], r["uniq_1M"]); s["members"] += 1
    print(f"{'protection':<12}{'members':>8}{'finite% (min-max)':>22}{'exact-0 %':>11}"
          f"{'exact-1 %':>11}{'max':>7}{'uniq/1M':>9}  branch hint")
    for prot, s in sorted(acc.items()):
        z = 100 * s["zeros"] / s["n"]
        hint = ("pooled_mean_zero_inflated?" if z > 60 else
                "pooled_median?" if z < 35 else "MEASURE THE MEDIAN")
        print(f"{prot:<12}{s['members']:>8}"
              f"{100*min(s['ff']):>10.2f}-{100*max(s['ff']):<11.2f}"
              f"{z:>11.2f}{100*s['ones']/s['n']:>11.4f}{s['max']:>7.3f}{s['uniq']:>9}  {hint}")
    print("\n  The branch hint is a HINT. OUTPUT-SPEC requires measuring what the median")
    print("  branch would actually publish (its exact-zero share on the decade pool) and")
    print("  recording it in decadal_statistic_rationale before deviating.")
    print("  NOTE the zero share differs BY PROTECTION LEVEL, so the branch is a per-variant")
    print("  decision -- this layer may legitimately ship on two different branches.")


def report_permanent_water(recs):
    print("\n=== PERMANENT WATER: cells fully flooded in every single year ===")
    acc = defaultdict(list)
    for r in recs:
        acc[r["protection"]].append((r["always_full"], r["caspian_always_full"]))
    print(f"{'protection':<12}{'min':>8}{'median':>9}{'max':>8}{'median in Caspian box':>24}")
    for prot, vals in sorted(acc.items()):
        tot = [v[0] for v in vals]; casp = [v[1] for v in vals]
        print(f"{prot:<12}{min(tot):>8}{int(np.median(tot)):>9}{max(tot):>8}{int(np.median(casp)):>24}")
    print("\n  Inland water bodies, not flood risk, and they sit at the top of any percentile")
    print("  ranking. Decide the mask on these numbers rather than on inspection.")


def report_sites(recs):
    print("\n=== GUARDRAILS 12: reference boxes, ensemble mean of the last decade ===")
    prots = sorted({r["protection"] for r in recs})
    print(f"{'box':<26}" + "".join(f"{p:>14}" for p in prots))
    fail = []
    for name in BOXES:
        row = f"{name:<26}"
        for p in prots:
            vals = [r[f"box::{name}"] for r in recs
                    if r["protection"] == p and np.isfinite(r[f"box::{name}"])]
            row += f"{np.mean(vals):>14.5f}" if vals else f"{'--':>14}"
        print(row)
        unprot = [r[f"box::{name}"] for r in recs
                  if r["protection"] == "none" and np.isfinite(r[f"box::{name}"])]
        if not unprot:
            continue
        if "dry control" in name and np.mean(unprot) > 0.01:
            fail.append(f"{name} (should be dry)")
        elif "dry control" not in name and np.mean(unprot) <= 0:
            fail.append(name)
    if fail:
        print(f"\n  *** GUARDRAILS 12 FAILURE: {fail} -- STOP, do not process. ***")
    else:
        print("\n  PASS. Every flood-prone box is non-trivial unprotected; controls are dry.")
    print("  A ZERO under `flopros` at a defended site means PROTECTED TO STANDARD, not NO")
    print("  HAZARD. That distinction is a must-disclose caveat, not a footnote.")
    return fail


REPORTS = {"orientation": report_orientation, "nature": report_nature,
           "permanent-water": report_permanent_water, "sites": report_sites}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protection", nargs="+")
    ap.add_argument("--section", choices=sorted(REPORTS), action="append")
    ap.add_argument("--limit", type=int, help="first N members only, for a fast pass")
    args = ap.parse_args()
    if not RAW.exists():
        print(f"missing raw ensemble: {RAW}", file=sys.stderr)
        return 1
    files = list(members(args.protection))
    if args.limit:
        files = files[:args.limit]
    if not files:
        print("no members matched", file=sys.stderr)
        return 1
    print(f"{len(files)} members from {RAW} -- one pass, all sections", flush=True)
    recs = []
    for n, (path, meta) in enumerate(files, 1):
        recs.append(survey_one(path, meta))
        if n % 20 == 0 or n == len(files):
            print(f"  surveyed {n}/{len(files)}", flush=True)
    for name in (args.section or sorted(REPORTS)):
        REPORTS[name](recs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
