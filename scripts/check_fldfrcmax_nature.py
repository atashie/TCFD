#!/usr/bin/env python3
"""Measure the data nature, orientation, mask and place-plausibility of raw `fldfrcmax`.

Usage:
    python scripts/check_fldfrcmax_nature.py                     # all sections, all members
    python scripts/check_fldfrcmax_nature.py --protection flopros
    python scripts/check_fldfrcmax_nature.py --section orientation

Run BEFORE processing and again before shipping. Four guardrail checks that need the same
open files:

  * GUARDRAILS 9 -- WHAT the values are. Measured on the probe: continuous [0,1], 283,749
    unique values, NOT binary. The decisive number for the statistic branch is the zero
    share, and here it **depends on the protection level**: 27.0% exact zeros unprotected
    against 87.5% under flopros, on the same member. So the branch cannot be chosen once
    for the layer -- it is a per-variant measurement, and this script prints it per variant.

  * ORIENTATION -- specific to this publication and genuinely dangerous. TipESM2025 flips
    the latitude axis BETWEEN PROTECTION VARIANTS OF THE SAME MEMBER: `none` is descending
    (89.875 -> -89.875) with header dims (time, lon, lat), `2yr`/`flopros` ascending with
    (time, lat, lon). Consequences, both silent: `.sel(lat=slice(hi, lo))` returns ZERO
    cells on half the files, and positional indexing builds an upside-down layer that
    satisfies every line of OUTPUT-SPEC. Everything downstream must `.sortby("lat")`; this
    section asserts that each file's data agrees with its own declared coordinates by
    checking that the Sahara is dry and the Caspian is wet in *coordinate* space.

  * PERMANENT WATER -- 515 cells read exactly 1.000 in the unprotected 2090s mean: the
    Caspian and other inland water bodies, permanently "flooded". They are not flood risk
    and they would dominate a percentile ranking. This section counts and locates them so
    the mask decision is made on numbers.

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


def section_orientation(files):
    print("\n=== ORIENTATION: declared axis vs where the water actually is ===")
    print(f"{'ghm':<18}{'gcm':<15}{'scen':<9}{'prot':<9}{'lat axis':<7}{'dims':<20}"
          f"{'Caspian':>9}{'Sahara':>9}{'verdict':>9}")
    bad = 0
    for path, m in files:
        ds = xr.open_dataset(path, decode_times=False)
        lat = ds["lat"].values
        order = "DESC" if lat[0] > lat[-1] else "ASC"
        dims = str(ds[VAR].dims).replace("'", "").replace(" ", "")
        v = ds[VAR].sortby("lat").sortby("lon")
        casp = float(np.nanmean(box(v, CASPIAN).values))
        sah = float(np.nanmean(np.nan_to_num(box(v, BOXES["Sahara (dry control)"]).values, nan=0.0)))
        # After sorting, the Caspian must be wetter than the Sahara. If a file's data were
        # stored against a flipped axis, this inverts -- the Sahara would sit where the
        # Caspian is and the check fires.
        good = casp > sah
        bad += 0 if good else 1
        print(f"{m['ghm']:<18}{m['gcm']:<15}{m['scenario']:<9}{m['protection']:<9}{order:<7}"
              f"{dims:<20}{casp:>9.4f}{sah:>9.4f}{'ok' if good else 'FLIPPED':>9}")
        ds.close()
    print(f"\n  files whose data disagrees with their own coordinates: {bad}")
    return bad


def section_nature(files):
    print("\n=== GUARDRAILS 9: data nature, PER PROTECTION LEVEL ===")
    acc = defaultdict(lambda: {"n": 0, "zeros": 0, "ones": 0, "finite_frac": [],
                               "max": 0.0, "uniq": 0, "members": 0})
    for path, m in files:
        ds, v = load(path)
        a = v.values
        fin = np.isfinite(a)
        vals = a[fin]
        s = acc[m["protection"]]
        s["n"] += vals.size
        s["zeros"] += int((vals == 0).sum())
        s["ones"] += int((vals == 1).sum())
        s["finite_frac"].append(float(fin.any(axis=0).mean()))
        s["max"] = max(s["max"], float(vals.max()))
        s["uniq"] = max(s["uniq"], int(np.unique(vals[:1_000_000]).size))
        s["members"] += 1
        ds.close()
    print(f"{'protection':<12}{'members':>8}{'finite% (min-max)':>22}{'exact-0 %':>11}"
          f"{'exact-1 %':>11}{'max':>7}{'uniq/1M':>9}{'branch hint':>26}")
    for prot, s in sorted(acc.items()):
        z = 100 * s["zeros"] / s["n"]
        hint = ("pooled_mean_zero_inflated?" if z > 60 else
                "pooled_median?" if z < 35 else "MEASURE THE MEDIAN")
        print(f"{prot:<12}{s['members']:>8}"
              f"{100*min(s['finite_frac']):>10.2f}-{100*max(s['finite_frac']):<11.2f}"
              f"{z:>11.2f}{100*s['ones']/s['n']:>11.4f}{s['max']:>7.3f}{s['uniq']:>9}{hint:>26}")
    print("\n  The branch hint is a HINT. OUTPUT-SPEC requires measuring what the median")
    print("  branch would actually publish (its exact-zero share on the decade pool) and")
    print("  recording it in decadal_statistic_rationale before deviating.")


def section_permanent_water(files):
    print("\n=== PERMANENT WATER: cells that are always fully flooded ===")
    print(f"{'ghm':<18}{'scen':<9}{'prot':<9}{'cells==1.0 all yrs':>20}{'in Caspian box':>16}")
    for path, m in files:
        ds, v = load(path)
        a = v.values
        always = np.all(np.nan_to_num(a, nan=0.0) >= 1.0, axis=0) & np.isfinite(a).any(axis=0)
        casp = box(v, CASPIAN)
        casp_always = np.all(np.nan_to_num(casp.values, nan=0.0) >= 1.0, axis=0).sum()
        print(f"{m['ghm']:<18}{m['scenario']:<9}{m['protection']:<9}{int(always.sum()):>20}{int(casp_always):>16}")
        ds.close()
    print("\n  These are inland water bodies, not flood risk. They will sit at the top of")
    print("  any percentile ranking. Decide the mask on these numbers, not on inspection.")


def section_sites(files):
    print("\n=== GUARDRAILS 12: reference boxes, mean of the last decade ===")
    per = defaultdict(lambda: defaultdict(list))
    for path, m in files:
        ds, v = load(path)
        last = v.isel(time=slice(-10, None))
        for name, b in BOXES.items():
            sub = box(last, b).values
            fin = np.isfinite(sub)
            per[name][m["protection"]].append(float(np.nanmean(sub)) if fin.any() else np.nan)
        ds.close()
    prots = sorted({p for d in per.values() for p in d})
    print(f"{'box':<26}" + "".join(f"{p:>14}" for p in prots))
    fail = []
    for name, d in per.items():
        row = f"{name:<26}"
        for p in prots:
            vals = [x for x in d.get(p, []) if np.isfinite(x)]
            row += f"{np.mean(vals):>14.5f}" if vals else f"{'--':>14}"
        print(row)
        if "dry control" not in name:
            unprot = [x for x in d.get("none", []) if np.isfinite(x)]
            if unprot and np.mean(unprot) <= 0:
                fail.append(name)
        else:
            unprot = [x for x in d.get("none", []) if np.isfinite(x)]
            if unprot and np.mean(unprot) > 0.01:
                fail.append(f"{name} (should be dry)")
    if fail:
        print(f"\n  *** GUARDRAILS 12 FAILURE: {fail} -- STOP, do not process. ***")
    else:
        print("\n  PASS. Every flood-prone box is non-trivial unprotected; controls are dry.")
    print("  Reading a ZERO under `flopros` at a defended site means PROTECTED TO STANDARD,")
    print("  not NO HAZARD. That distinction is a must-disclose caveat, not a footnote.")
    return fail


SECTIONS = {"orientation": section_orientation, "nature": section_nature,
            "permanent-water": section_permanent_water, "sites": section_sites}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protection", nargs="+")
    ap.add_argument("--section", choices=sorted(SECTIONS), action="append")
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
    print(f"{len(files)} members from {RAW}")
    for name in (args.section or sorted(SECTIONS)):
        SECTIONS[name](files)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
