"""Measure the data nature, mask and COVERAGE SPARSITY of raw `leh` before processing.

Same three guardrail checks as the ISIMIP3b heatwave check (GUARDRAILS 9 nature, the mask,
GUARDRAILS 12 reference sites) plus a fourth that is specific to this layer and is the
reason it was ingested at all.

THE SPARSITY TEST IS THE POINT OF THIS SCRIPT.

`leh` counts a cell as exposed only when BOTH criteria fire: the relative HWMId AND the
absolute Humidex >= 45. Lange et al. 2020 report that below 2 degC of global warming the
Humidex "hardly ever exceeds their threshold value of 45 in the temperate regions", and
ISIMIP2b publishes `leh` for rcp26/rcp60 ONLY -- exactly the forcing range where that gate
mostly stays shut. If that is right, the field is near-zero across Europe, Canada and the
northern US, and a layer built on it would report "no heat hazard" for a large part of any
Northern-Hemisphere portfolio.

THAT IS A CLAIM FROM A PAPER. This script turns it into a measurement, because the
difference between "sparse by design, and here is where" and "sparse everywhere, unusable"
decides whether this can be a primary layer. A sparse ABSOLUTE index is not a defect -- it
is what an absolute threshold does -- but shipping one without knowing WHERE it is silent
would repeat the sugarcane failure exactly: a layer that passes every contract check while
reading zero across the region a reader cares about.

The failure mode to look for is therefore NOT "lots of zeros". It is a reference site that
is zero in every member and every year, and a temperate belt with no signal in any decade.

Usage:
    python scripts/check_leh_nature.py [--limit N]
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr

VAR = "leh"
LAYER_ID = "heatwave-isimip2b_leh_annual"

#: Same sites as the ISIMIP3b heatwave check, so the two indices can be read against each
#: other cell for cell. The hot/temperate contrast is the whole question here.
REFERENCE_SITES = [
    ("Phoenix, AZ (hot arid)",       33.4, -112.1),
    ("Kuwait City (hot arid)",       29.4,   47.9),
    ("Jacobabad, PK (hot humid)",    28.3,   68.4),
    ("Delhi, IN",                    28.6,   77.2),
    ("Cordoba, ES",                  37.9,   -4.8),
    ("Paris, FR (temperate)",        48.9,    2.4),
    ("Chicago, US (temperate)",      41.9,  -87.6),
    ("Frankfurt, DE (data centre)",  50.1,    8.7),
    ("Singapore (tropical)",          1.4,  103.8),
    ("Lagos, NG",                     6.5,    3.4),
    ("Sao Paulo, BR",               -23.5,  -46.6),
    ("Sydney, AU",                  -33.9,  151.2),
    ("Yakutsk, RU (boreal)",         62.0,  129.7),
]

#: Latitude belts for the sparsity survey. "Temperate" is the band the Humidex gate is
#: reported to leave shut below 2 degC, and is where most disclosed assets sit.
BELTS = [(-90, -60, "60-90S"), (-60, -35, "35-60S"), (-35, -23.5, "23.5-35S"),
         (-23.5, 23.5, "tropics"), (23.5, 35, "23.5-35N"), (35, 60, "35-60N (temperate)"),
         (60, 90, "60-90N (boreal)")]


def parse_name(fpath):
    """(model, gcm, scenario, member) read from the END of the filename.

    ISIMIP2b Lange2020 files carry a LEADING publication token AND the `.nc4` extension:
    lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099.nc4
    Reading from the end is index-stable across the round's grammar differences.
    """
    p = os.path.basename(fpath).split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def decode_years(ds, fpath):
    """Years per record, without trusting the time axis to declare itself.

    The time convention is per-PUBLICATION in ISIMIP: Lange2020 (2b), Heinicke2026 and
    Zantout2025 (both 3b) each use a different epoch and step. The span is taken from the
    filename -- which IS part of the ISIMIP grammar -- the record count checked against it,
    and the step convention classified and reported rather than assumed.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    y0, y1 = int(p[-2]), int(p[-1])
    t = np.asarray(ds["time"].values, dtype="float64")
    n = y1 - y0 + 1
    if t.size != n:
        raise ValueError(f"{os.path.basename(fpath)}: {t.size} records but the filename "
                         f"declares {y0}-{y1} ({n} years)")
    d = np.diff(t)
    if d.size and (d <= 0).any():
        raise ValueError(f"{os.path.basename(fpath)}: time axis is not increasing")
    step = float(np.median(d)) if d.size else float("nan")
    if d.size and np.allclose(d, 1.0):
        conv = "integer years (step 1)"
    elif d.size and 359.0 <= step <= 367.0:
        conv = f"days (median step {step:.2f})"
    else:
        conv = f"UNRECOGNISED (median step {step:.4g})"
    return np.arange(y0, y1 + 1, dtype=int), conv


def site_indices(lats, lons, lat, lon):
    lon_q = lon if lons.min() < 0 else (lon % 360.0)
    return int(np.abs(lats - lat).argmin()), int(np.abs(lons - lon_q).argmin())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc4")))
    if not files:
        print(f"ERROR: no {VAR} files in {raw_dir}")
        return 2
    if args.limit:
        files = files[:args.limit]

    print("=" * 78)
    print(f"leh (ISIMIP2b heatwave exposure, HWMId AND Humidex>=45) -- {len(files)} members")
    print("=" * 78)

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    print(f"grid: lat {lats.min():.2f}..{lats.max():.2f} ({len(lats)}), "
          f"lon {lons.min():.2f}..{lons.max():.2f} ({len(lons)})")
    sites = [(name, *site_indices(lats, lons, la, lo), la, lo)
             for name, la, lo in REFERENCE_SITES]

    rows, attr_sets = [], defaultdict(set)
    finite_mask, nonzero_mask = {}, {}
    site_series, site_late = defaultdict(dict), defaultdict(dict)
    belt_acc = defaultdict(list)
    global_min, global_max = np.inf, -np.inf
    any_non_binary = False
    ever_nonzero = None

    for i, f in enumerate(files, 1):
        info = parse_name(f)
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[VAR]
            attr_sets["units"].add(str(da.attrs.get("units", "<absent>")))
            attr_sets["long_name"].add(str(da.attrs.get("long_name", "<absent>")))
            attr_sets["time_units"].add(str(ds.time.attrs.get("units", "<absent>")))
            attr_sets["calendar"].add(str(ds.time.attrs.get("calendar", "<absent>")))
            yrs, conv = decode_years(ds, f)
            attr_sets["time_step_convention"].add(conv)
            vals = da.values.astype("float32")
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))

        if fill is not None:
            vals = np.where(np.isclose(vals, np.float32(fill), rtol=1e-6), np.nan, vals)
        vals[~np.isfinite(vals)] = np.nan

        fin = vals[np.isfinite(vals)]
        if fin.size == 0:
            print(f"  [{i}/{len(files)}] {info['member']:<28} ALL-NaN MEMBER -- STOP")
            return 3

        n_uniq = np.unique(fin).size
        exact0 = float((fin == 0.0).mean())
        exact1 = float((fin == 1.0).mean())
        binary = bool(n_uniq <= 2 and np.all((fin == 0.0) | (fin == 1.0)))
        if not binary:
            any_non_binary = True
        global_min = min(global_min, float(fin.min()))
        global_max = max(global_max, float(fin.max()))

        valid = np.isfinite(vals)
        mask2d = valid.any(axis=0)
        time_invariant = bool(np.array_equal(valid, np.broadcast_to(mask2d, valid.shape)))
        nz2d = (np.nan_to_num(vals) > 0).any(axis=0)
        key = f"{info['member']}|{info['scenario']}"
        finite_mask[key] = mask2d
        nonzero_mask[key] = nz2d
        ever_nonzero = nz2d.copy() if ever_nonzero is None else (ever_nonzero | nz2d)

        # Per-belt mean over ALL years, and over the final decade, this member.
        late = yrs >= 2090
        for lo, hi, name in BELTS:
            sel = (lats >= lo) & (lats < hi)
            band_all = vals[:, sel, :]
            band_late = vals[late][:, sel, :]
            fa = band_all[np.isfinite(band_all)]
            fl = band_late[np.isfinite(band_late)]
            if fa.size:
                belt_acc[name].append((float(fa.mean()),
                                       float(fl.mean()) if fl.size else np.nan,
                                       float((fa == 0).mean())))

        for name, ilat, ilon, _, _ in sites:
            v = vals[:, ilat, ilon]
            site_series[name][key] = (float(np.nanmean(v)) if np.isfinite(v).any()
                                      else float("nan"))
            vl = v[late]
            site_late[name][key] = (float(np.nanmean(vl)) if np.isfinite(vl).any()
                                    else float("nan"))

        rows.append(dict(member=info["member"], gcm=info["gcm"], scenario=info["scenario"],
                         n_unique=int(n_uniq), vmin=float(fin.min()), vmax=float(fin.max()),
                         exact0=exact0, exact1=exact1, binary=binary,
                         finite_cells=int(mask2d.sum()), nonzero_cells=int(nz2d.sum()),
                         time_invariant_mask=time_invariant,
                         n_years=int(yrs.size), y0=int(yrs.min()), y1=int(yrs.max())))
        print(f"  [{i}/{len(files)}] {info['gcm']:<14} {info['scenario']:<7} "
              f"uniq={n_uniq:<5} range=[{fin.min():.4g},{fin.max():.4g}] "
              f"zero={exact0:.2%}  nonzero-ever cells={int(nz2d.sum()):,}", flush=True)

    print("\n" + "-" * 78)
    print("1. METADATA (per-member, not assumed uniform)")
    print("-" * 78)
    for k, v in attr_sets.items():
        flag = "" if len(v) == 1 else "   <-- DIVERGES ACROSS MEMBERS"
        print(f"  {k:<22} {sorted(v)}{flag}")

    print("\n" + "-" * 78)
    print("2. DATA NATURE (GUARDRAILS 9 -- from the values)")
    print("-" * 78)
    n_bin = sum(r["binary"] for r in rows)
    z = np.array([r["exact0"] for r in rows])
    print(f"  binary members: {n_bin}/{len(rows)}   "
          f"distinct-value counts: {sorted({r['n_unique'] for r in rows})}")
    print(f"  global value range: [{global_min:.6g}, {global_max:.6g}]")
    print(f"  exact-0 share: min {z.min():.2%}  median {np.median(z):.2%}  max {z.max():.2%}")
    if n_bin == len(rows):
        verdict = ("BOOLEAN {0,1} in every member -> pooled_mean_boolean branch. Note this "
                   "was NOT inherited: the sibling `let` in this same publication is a "
                   "CONTINUOUS fraction.")
    elif any_non_binary and n_bin > 0:
        verdict = "MIXED -- STOP: the contract takes ONE branch per layer."
    else:
        verdict = ("CONTINUOUS -> median/IQR branch unless the decade pool is degenerate "
                   "at zero; check the exact-zero share against let's 97.84%")
    print(f"  VERDICT: {verdict}")

    print("\n" + "-" * 78)
    print("3. MASK")
    print("-" * 78)
    bad = [r for r in rows if not r["time_invariant_mask"]]
    print(f"  time-invariant per member: {len(rows) - len(bad)}/{len(rows)}")
    fin_cells = sorted({r["finite_cells"] for r in rows})
    print(f"  finite cells per member: {fin_cells}"
          f"  ({100.0 * max(fin_cells) / (len(lats) * len(lons)):.1f}% of the grid)")
    land = None
    try:
        sys.path.insert(0, str(root / "scripts"))
        from utils.land_mask import get_isimip_landmask  # noqa: E402
        lm = xr.open_dataset(get_isimip_landmask("2b"))
        # The mask VARIABLE NAME differs by round: ISIMIP2b publishes `LSM` (with a
        # singleton time axis), ISIMIP3b publishes `mask`. Resolve rather than assume --
        # assuming `mask` here silently disabled the land denominator on the first run and
        # the sparsity figures came out quoted against the whole globe, ocean included.
        name = "LSM" if "LSM" in lm.variables else "mask"
        arr = np.squeeze(lm[name].values)
        land = np.nan_to_num(arr) > 0
        lm.close()
    except Exception as e:  # noqa: BLE001
        print(f"  (land mask unavailable: {type(e).__name__}: {e})")
    if land is not None:
        fu = np.zeros(land.shape, bool)
        for a in finite_mask.values():
            fu |= a
        print(f"  finite union {int(fu.sum()):,} | on land {int((fu & land).sum()):,} | "
              f"off land {int((fu & ~land).sum()):,}")
        print(f"  ISIMIP2b land cells: {int(land.sum()):,}")
        if fu.mean() > 0.9:
            print("  -> THE PRODUCT ZERO-FILLS THE WHOLE GRID: `isfinite` is NOT a mask "
                  "here, so every share below is quoted against LAND, not against finite "
                  "cells. Same publisher behaviour as `cropfailure`.")

    print("\n" + "-" * 78)
    print("4. SPARSITY -- WHERE IS THE ABSOLUTE HUMIDEX GATE SILENT? (this layer's question)")
    print("-" * 78)
    fu = np.zeros(lats.shape + lons.shape, bool)
    for a in finite_mask.values():
        fu |= a
    # Quote sparsity against LAND. The product zero-fills ocean, so a finite-cell
    # denominator would dilute every share with 66% ocean and make the layer look far more
    # active than it is -- the same dilution that once turned a 66-76% sen-zero share into
    # a reported 20.7% by counting ocean.
    domain = (fu & land) if land is not None else fu
    dom_label = "land" if land is not None else "finite (NO LAND MASK AVAILABLE)"
    silent = domain & ~ever_nonzero
    print(f"  {dom_label} cells NEVER non-zero in ANY member or year: "
          f"{int(silent.sum()):,} of {int(domain.sum()):,} "
          f"({100.0 * silent.sum() / max(domain.sum(), 1):.1f}%)")
    if silent.any():
        si, _ = np.where(silent)
        print(f"    their latitude range: {lats[si].min():.1f} .. {lats[si].max():.1f}; "
              f"share poleward of 35 deg: {(np.abs(lats[si]) > 35).mean():.1%}")
    print(f"\n  {'belt':<20} {'mean all yr':>12} {'mean 2090s':>11} {'exact-0':>9}"
          f" {'silent / ' + dom_label:>20}")
    for lo, hi, name in BELTS:
        acc = belt_acc.get(name)
        if not acc:
            continue
        a = np.array(acc)
        sel = (lats >= lo) & (lats < hi)
        sil = int((silent[sel, :]).sum())
        tot = int((domain[sel, :]).sum())
        pct = f"{100.0 * sil / tot:.0f}%" if tot else "  -"
        print(f"  {name:<20} {a[:, 0].mean():>12.4f} {np.nanmean(a[:, 1]):>11.4f} "
              f"{a[:, 2].mean():>8.1%} {sil:>8,}/{tot:<7,} {pct:>5}")
    print("\n  Reading: a HIGH silent-cell count in the temperate belt is the EXPECTED "
          "consequence of the absolute Humidex>=45 gate under rcp26/rcp60, not a defect -- "
          "but it decides whether this layer can carry a Northern-Hemisphere portfolio. "
          "A belt that is silent everywhere cannot be reported as 'low heat risk'; it has "
          "to be reported as 'below the health threshold this index measures'.")

    print("\n" + "-" * 78)
    print("5. REFERENCE SITES (GUARDRAILS 12)")
    print("-" * 78)
    print(f"  {'site':<28} {'lat':>6} {'lon':>7} {'mem':>4} {'NaN':>4} "
          f"{'mean all yr':>12} {'mean 2090s':>11}")
    stops = []
    for name, ilat, ilon, la, lo in sites:
        vals = np.array(list(site_series[name].values()), dtype=float)
        late = np.array(list(site_late[name].values()), dtype=float)
        n_nan = int(np.isnan(vals).sum())
        f_ = vals[np.isfinite(vals)]
        l_ = late[np.isfinite(late)]
        mean = float(f_.mean()) if f_.size else float("nan")
        lmean = float(l_.mean()) if l_.size else float("nan")
        flag = ""
        if f_.size == 0:
            flag = "  <-- NaN in EVERY member"
            stops.append(name)
        elif np.allclose(f_, 0.0) and np.allclose(l_, 0.0):
            flag = "  <-- ZERO in every member, every year"
            stops.append(name)
        print(f"  {name:<28} {la:>6.1f} {lo:>7.1f} {len(vals):>4} {n_nan:>4} "
              f"{mean:>12.4f} {lmean:>11.4f}{flag}")

    out = raw_dir / "value_check.json"
    with open(out, "w") as fh:
        json.dump(dict(n_members=len(rows), binary_members=n_bin,
                       value_range=[global_min, global_max],
                       exact0_median=float(np.median(z)),
                       silent_cells=int(silent.sum()), finite_union=int(fu.sum()),
                       attrs={k: sorted(v) for k, v in attr_sets.items()}, rows=rows,
                       sites_all_years={n: site_series[n] for n in site_series},
                       sites_2090s={n: site_late[n] for n in site_late}), fh, indent=2)
    print(f"\nwrote {out}")

    if stops:
        print(f"\nATTENTION -- {len(stops)} reference site(s) are NaN or identically zero: "
              f"{stops}.\nOn an ABSOLUTE-threshold index a zero can be a true 'never crosses "
              "the health threshold' rather than a defect (contrast the sugarcane case, "
              "where the model simply did not run). Decide which it is from the belt table "
              "above BEFORE building a layer, and carry the answer as a caveat either way.")
        return 4
    print("\nNo reference site is silent. Data nature above governs the processor's branch.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
