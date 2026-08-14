"""Measure the data nature, mask and site coverage of raw `heatwave` before processing.

Three INDEPENDENT guardrail checks, run together because they need the same open files but
answer different questions:

  * GUARDRAILS 9 -- WHAT the values are. Data nature is taken from the values, never from
    the variable name, the CF `long_name`, a sibling variable, another layer, or a paper.
    This one has all four temptations pointing in conflicting directions: the header
    declares `long_name = "heatwave area share"` and `units = "1"` (reads continuous), the
    Zantout 2025 Methods say "The exposed grid cell area fraction is one if the percentile
    threshold is exceeded and zero else" (reads binary), the 2b sibling `leh` is unverified,
    and within Heinicke2026 `floodedarea` is binary while `driedarea` is not. Only the
    values decide, and what they decide is the decadal-statistic branch, the CI definition
    and the percentile tiering.

  * The MASK. `floodedarea` -- same round, adjacent publication -- is non-NaN over 94.7% of
    the globe INCLUDING open ocean, which no contract check would catch. `cropfailure`
    (same publication as this layer) zero-fills the globe instead, so `isfinite` is not a
    footprint there either. Both failure modes are checked here against the official
    ISIMIP3b land-sea mask: finite extent, non-zero-ever extent, and how much of each
    falls off land.

  * GUARDRAILS 12 -- WHERE the values are. A layer about extreme heat must be non-trivial
    where extreme heat demonstrably occurs. The withdrawn sugarcane layers passed every
    contract check while reading exactly 0 across the entire cane belt. Counting zeros does
    not locate them.

READ THE SITE TABLE CORRECTLY -- this index is RELATIVE, not absolute. A cell is exposed
when its HWMId exceeds the 97.5th percentile of ITS OWN picontrol distribution, so the
field measures departure from a per-cell preindustrial baseline, not absolute heat. A
permanently hot site (Kuwait, Jacobabad) reading LOWER than a temperate one (Paris,
Chicago) is the index behaving as designed, NOT a defect to fix. The failure mode is a NaN
in every member, or an exact zero in every member across every year and scenario.

Usage:
    python scripts/check_heatwave_nature.py [--limit N]
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

VAR = "heatwave"
LAYER_ID = "heatwave-isimip3b_heatwave_annual"

#: GUARDRAILS 12 reference sites. Chosen to span the two things this layer must get right:
#: places where extreme heat is a first-order risk to people and equipment (the customer
#: framing), and a deliberate hot/temperate/cold contrast that exposes the relative-baseline
#: behaviour so it is visible in evidence rather than discovered by a reader.
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


def parse_name(fpath):
    """(model, gcm, scenario, member) read from the END of the filename.

    Zantout2025 carries a LEADING publication token (`zantout2025_`) where its 3b sibling
    Heinicke2026 does not, so a forward index reads the wrong column. From the end it is
    invariant under either convention.
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
    """Years for each record, WITHOUT trusting the time axis to declare itself.

    The time convention is NOT uniform within the Zantout2025 publication: the sidecar for
    `heatwave` declares `units = "days since 2015-01-01"`, while the same publication's
    `cropfailure` files carry no `units` at all and store a bare integer year sequence. So
    the span is taken from the filename (which IS part of the ISIMIP grammar), the record
    count is checked against it, and the step convention is CLASSIFIED and reported rather
    than assumed -- a member that disagrees with the rest is then visible instead of
    silently mis-binned by a decade.
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
    """Nearest grid cell, handling either a -180..180 or a 0..360 longitude axis."""
    lon_q = lon if lons.min() < 0 else (lon % 360.0)
    return int(np.abs(lats - lat).argmin()), int(np.abs(lons - lon_q).argmin())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="check only the first N members (iteration only -- the shipped "
                         "check must be exhaustive, a sampled scan can step over a rare "
                         "continuous tail and silently switch the published statistic)")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc")))
    if not files:
        print(f"ERROR: no {VAR} files in {raw_dir}")
        return 2
    if args.limit:
        files = files[:args.limit]

    print("=" * 78)
    print(f"heatwave raw value check -- {len(files)} members")
    print("=" * 78)

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    print(f"grid: lat {lats.min():.2f}..{lats.max():.2f} ({len(lats)}), "
          f"lon {lons.min():.2f}..{lons.max():.2f} ({len(lons)})")
    sites = [(name, *site_indices(lats, lons, la, lo), la, lo)
             for name, la, lo in REFERENCE_SITES]

    # ---- per-member scan ------------------------------------------------------ #
    rows = []
    attr_sets = defaultdict(set)
    finite_mask = {}
    nonzero_mask = {}
    site_series = defaultdict(dict)     # site -> member|scenario -> mean over all years
    site_late = defaultdict(dict)       # site -> member|scenario -> mean over 2090s
    global_min, global_max = np.inf, -np.inf
    any_non_binary = False

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
            print(f"  [{i}/{len(files)}] {info['member']:<28} {info['scenario']}  "
                  "ALL-NaN MEMBER -- STOP")
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
        time_invariant = bool(np.array_equal(
            valid, np.broadcast_to(mask2d, valid.shape)))
        nz2d = (np.nan_to_num(vals) > 0).any(axis=0)
        key = info["member"]
        finite_mask.setdefault(key, np.zeros(mask2d.shape, bool))
        finite_mask[key] |= mask2d
        nonzero_mask.setdefault(key, np.zeros(nz2d.shape, bool))
        nonzero_mask[key] |= nz2d

        late = (yrs >= 2090)
        for name, ilat, ilon, _, _ in sites:
            v = vals[:, ilat, ilon]
            tag = f"{info['member']}|{info['scenario']}"
            site_series[name][tag] = (
                float(np.nanmean(v)) if np.isfinite(v).any() else float("nan"))
            vl = v[late]
            site_late[name][tag] = (
                float(np.nanmean(vl)) if np.isfinite(vl).any() else float("nan"))

        rows.append(dict(
            member=info["member"], model=info["model"], gcm=info["gcm"],
            scenario=info["scenario"],
            n_unique=int(n_uniq), vmin=float(fin.min()), vmax=float(fin.max()),
            exact0=exact0, exact1=exact1, binary=binary,
            finite_cells=int(mask2d.sum()),
            finite_frac=float(mask2d.mean()),
            nonzero_cells=int(nz2d.sum()),
            time_invariant_mask=time_invariant,
            n_years=int(yrs.size), y0=int(yrs.min()), y1=int(yrs.max()),
        ))
        print(f"  [{i}/{len(files)}] {info['member']:<26} {info['scenario']:<7} "
              f"uniq={n_uniq:<6} range=[{fin.min():.4g},{fin.max():.4g}] "
              f"zero={exact0:.2%}", flush=True)

    # ---- 1. metadata consistency --------------------------------------------- #
    print("\n" + "-" * 78)
    print("1. METADATA (per-member, not assumed uniform)")
    print("-" * 78)
    for k, v in attr_sets.items():
        flag = "" if len(v) == 1 else "   <-- DIVERGES ACROSS MEMBERS"
        print(f"  {k:<22} {sorted(v)}{flag}")

    # ---- 2. data nature ------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("2. DATA NATURE (GUARDRAILS 9 -- from the values, not the header or the paper)")
    print("-" * 78)
    n_bin = sum(r["binary"] for r in rows)
    uniq_counts = sorted({r["n_unique"] for r in rows})
    print(f"  binary members: {n_bin}/{len(rows)}")
    print(f"  distinct-value counts across members: {uniq_counts[:8]}"
          + (" ..." if len(uniq_counts) > 8 else ""))
    print(f"  global value range: [{global_min:.6g}, {global_max:.6g}]")
    z = np.array([r["exact0"] for r in rows])
    o = np.array([r["exact1"] for r in rows])
    print(f"  exact-0 share of finite cell-years: min {z.min():.2%}  "
          f"median {np.median(z):.2%}  max {z.max():.2%}")
    print(f"  exact-1 share: min {o.min():.2%}  median {np.median(o):.2%}  "
          f"max {o.max():.2%}")
    print(f"  declared long_name says 'area share'; measured verdict follows the values:")
    if n_bin == len(rows):
        verdict = ("BOOLEAN {0,1} in every member -> pooled_mean_boolean branch "
                   "(mean +/- 1 SD). The declared 'area share' is a MISNOMER: the field is "
                   "an occurrence flag, so its decadal mean is a FREQUENCY, not an area.")
    elif any_non_binary and n_bin > 0:
        verdict = ("MIXED -- some members binary, some continuous. This is a STOP: the "
                   "contract takes ONE branch per layer and a mixed ensemble means the "
                   "members are not the same quantity.")
    else:
        verdict = ("CONTINUOUS -> median/IQR branch UNLESS the decade pool is degenerate "
                   "at zero; see the zero-inflation section below")
    print(f"  VERDICT: {verdict}")

    # ---- 3. zero-inflation regime -------------------------------------------- #
    print("\n" + "-" * 78)
    print("3. ZERO-INFLATION (decides whether the third branch is even in scope)")
    print("-" * 78)
    print(f"  annual exact-zero share, ensemble median: {np.median(z):.2%}")
    print("  reference points: let 97.84% (took pooled_mean_zero_inflated),")
    print("                    burntarea 29.2% (did NOT qualify, median branch)")
    print("  NOTE: this is the ANNUAL share. The third branch is justified only by the "
          "exact-zero share of the MEDIAN of the decade POOL, which the processor must "
          "measure on the 2020s panel. Do not decide from this number alone.")

    # ---- 4. mask survey ------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("4. MASK -- finite extent vs non-zero extent vs land")
    print("-" * 78)
    bad_mask = [r for r in rows if not r["time_invariant_mask"]]
    print(f"  time-invariant per member: {len(rows) - len(bad_mask)}/{len(rows)}")
    if bad_mask:
        print(f"  <-- {len(bad_mask)} member(s) have a TIME-VARYING mask; "
              "load_member() must not assume otherwise")
    fin_cells = sorted({r["finite_cells"] for r in rows})
    print(f"  finite cells per member: {fin_cells}"
          f"  ({100.0 * max(fin_cells) / (len(lats) * len(lons)):.1f}% of the grid at most)")

    land = None
    try:
        sys.path.insert(0, str(root / "scripts"))
        from utils.land_mask import get_isimip_landmask  # noqa: E402
        lm = xr.open_dataset(get_isimip_landmask("3b"))["mask"].values
        land = np.nan_to_num(lm) > 0
    except Exception as e:  # noqa: BLE001
        print(f"  (land mask unavailable: {type(e).__name__}: {e})")

    if land is not None:
        fin_union = np.zeros(land.shape, bool)
        for a in finite_mask.values():
            fin_union |= a
        nz_union = np.zeros(land.shape, bool)
        for a in nonzero_mask.values():
            nz_union |= a
        print(f"\n  ISIMIP3b land cells: {int(land.sum()):,}")
        print(f"  finite-ever union   : {int(fin_union.sum()):,}  "
              f"on land {int((fin_union & land).sum()):,}  "
              f"OFF LAND {int((fin_union & ~land).sum()):,}")
        print(f"  non-zero-ever union : {int(nz_union.sum()):,}  "
              f"on land {int((nz_union & land).sum()):,}  "
              f"OFF LAND {int((nz_union & ~land).sum()):,}")
        # Cells off the ISIMIP3b land-sea mask are NOT automatically an ocean leak. That
        # mask is coarse: it drops small islands the derived products keep, and the two
        # disagree along every coastline. What distinguishes a leak is SCALE and PLACEMENT
        # -- `floodedarea` is non-NaN over 94.7% of the GRID, values sitting on open ocean.
        # Measured reference points on this machine: `driedarea` (3b, shipped) carries
        # 1,236 finite cells off this same mask and `led` (2b, shipped) carries 46.
        off = int((nz_union & ~land).sum())
        grid_frac = float(fin_union.mean())
        pad = np.pad(land, 1, mode="edge")
        neigh = np.zeros(land.shape, dtype=int)
        for dy in (0, 1, 2):
            for dx in (0, 1, 2):
                neigh += pad[dy:dy + land.shape[0], dx:dx + land.shape[1]].astype(int)
        neigh -= land.astype(int)
        oi, oj = np.where(nz_union & ~land)
        isolated = int((neigh[oi, oj] == 0).sum()) if oi.size else 0
        print(f"    of which coastal fringe (>=1 land neighbour): {off - isolated:,}")
        print(f"    of which isolated (no land neighbour)       : {isolated:,}")
        print(f"  finite extent as a share of the GRID: {grid_frac:.1%}  "
              f"(floodedarea's leak was 94.7%)")
        if grid_frac > 0.5:
            print("    -> OCEAN LEAK: the field covers most of the globe. This is the "
                  "`floodedarea` failure mode; STOP and mask before processing.")
        else:
            print("    -> NOT an ocean leak: the extent is a land footprint. The off-mask "
                  "cells are a mask-CONVENTION difference (islands and coastline), the "
                  "same disagreement the shipped `driedarea` and `led` layers carry.")
        allzero_land = int((land & ~nz_union).sum())
        az_i, _ = np.where(land & ~nz_union)
        polar = float((lats[az_i] < -60).mean()) if az_i.size else 0.0
        print(f"  land cells all-zero in EVERY member: {allzero_land:,} "
              f"({100.0 * allzero_land / land.sum():.1f}% of land), of which "
              f"{polar:.1%} are south of -60 deg (Antarctica)")

        # Per-member spread decides the publication mask. n_models = 1 here, so `led`'s
        # ">=2 impact models" rule cannot apply -- any mask must be a GCM-count rule and
        # has to be measured, not inherited.
        stack = np.array([finite_mask[k] for k in sorted(finite_mask)])
        nmem = stack.sum(axis=0)
        print(f"\n  finite coverage tiers (cells inside exactly K of "
              f"{len(finite_mask)} GCM members):")
        for k in range(1, len(finite_mask) + 1):
            c = int((nmem == k).sum())
            if c:
                print(f"    {k:>2} member(s): {c:>8,}")

    # ---- 5. reference sites (GUARDRAILS 12) ----------------------------------- #
    print("\n" + "-" * 78)
    print("5. REFERENCE SITES (GUARDRAILS 12 -- non-trivial where extreme heat occurs?)")
    print("-" * 78)
    print(f"  {'site':<28} {'lat':>6} {'lon':>7} {'mem':>4} {'NaN':>4} "
          f"{'mean all yr':>12} {'mean 2090s':>11}")
    stops = []
    for name, ilat, ilon, la, lo in sites:
        vals = np.array(list(site_series[name].values()), dtype=float)
        late = np.array(list(site_late[name].values()), dtype=float)
        n_nan = int(np.isnan(vals).sum())
        finite = vals[np.isfinite(vals)]
        lfin = late[np.isfinite(late)]
        mean = float(finite.mean()) if finite.size else float("nan")
        lmean = float(lfin.mean()) if lfin.size else float("nan")
        flag = ""
        if finite.size == 0:
            flag = "  <-- STOP: NaN in EVERY member"
            stops.append(name)
        elif np.allclose(finite, 0.0) and np.allclose(lfin, 0.0):
            flag = "  <-- STOP: exactly 0 in every member and every year"
            stops.append(name)
        elif n_nan:
            flag = f"  <-- masked in {n_nan} member(s)"
        print(f"  {name:<28} {la:>6.1f} {lo:>7.1f} {len(vals):>4} {n_nan:>4} "
              f"{mean:>12.4f} {lmean:>11.4f}{flag}")

    print("\n  Reading: the 'mean all yr' column is 2015-2100, the '2090s' column is "
          "2090-2100. A rising pair is the expected signature of an unprecedented-event "
          "index. A permanently hot site reading BELOW a temperate one is the "
          "relative-baseline definition working as designed -- exposure is measured "
          "against each cell's OWN preindustrial control, not against an absolute "
          "temperature. Do not 'fix' it; disclose it.")

    # ---- summary -------------------------------------------------------------- #
    out = raw_dir / "value_check.json"
    with open(out, "w") as fh:
        json.dump(dict(
            n_members=len(rows), binary_members=n_bin,
            value_range=[global_min, global_max],
            exact0_median=float(np.median(z)), exact1_median=float(np.median(o)),
            attrs={k: sorted(v) for k, v in attr_sets.items()},
            rows=rows,
            sites_all_years={n: site_series[n] for n in site_series},
            sites_2090s={n: site_late[n] for n in site_late},
        ), fh, indent=2)
    print(f"\nwrote {out}")

    if stops:
        print(f"\nSTOP -- {len(stops)} reference site(s) are NaN or identically zero: "
              f"{stops}. Investigate before writing any output (GUARDRAILS 12).")
        return 4
    print("\nNo reference-site STOP condition. Data nature recorded above governs the "
          "processor's branch -- copy the measured numbers into its docstring.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
