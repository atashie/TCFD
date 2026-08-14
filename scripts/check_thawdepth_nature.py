"""Measure the data nature, domain and site behaviour of raw `thawdepth` before processing.

This layer is unlike every exposure layer shipped so far: `thawdepth` is a PHYSICAL FIELD in
metres, not a 0/1 exposure flag produced by a publisher who already made the framing
decisions. Nothing about the statistic, the mask or the direction is inherited -- all of it
is decided by what is measured here.

Five checks, run together because they need the same open files:

  * GUARDRAILS 9 -- WHAT the values are, PER MODEL. Three independent impact models
    (LPJmL5-7-10-fire, JULES-ES-VN6P3, CLASSIC) each implement their own soil column, and
    nothing guarantees they publish the same quantity on the same scale. A unit divergence
    (m vs mm) or a column-depth divergence between models does not break any contract check
    and silently makes the ensemble a mixture of incompatible fields, so per-model
    distributions are reported side by side and never pooled before they are compared.

  * THE ENCODING OF A NON-PERMAFROST CELL -- the decisive unknown, and one no listing can
    answer. A cell with no permafrost can be published as NaN, as 0, or as the full soil
    column depth. Each choice implies a different mask, a different statistic and a
    different meaning for a high value. The check looks for the "pinned at the column
    bottom" signature: a large share of cells sitting at exactly the per-model maximum.

  * SATURATION / CENSORING. If thaw reaches the bottom of the modelled soil column the value
    stops rising, so a cell that has completely lost its permafrost has slope ~0 -- reading
    as "no trend" while meaning "no permafrost left". That is the `heatwave` failure mode
    arriving from the opposite end, and it is measured here as the share of cells at the
    per-model maximum by decade and scenario, not assumed.

  * THE DOMAIN. No permafrost mask ships with this product. The diagnostic used here is the
    only one available from the field itself: a cell is treated as PERMAFROST-BEARING in a
    decade when its thaw depth stays strictly below that model's maximum. The count of such
    cells, and how it falls between the 2020s and the 2090s, is both the sanity check and
    the candidate hazard signal.

  * GUARDRAILS 12 -- WHERE the values are. Reference sites span continuous permafrost,
    discontinuous permafrost, alpine permafrost and permafrost-free controls. The controls
    are not decoration: what a permafrost-free cell reads is exactly the encoding question
    above, and Nairobi answers it more directly than any global statistic.

READ THE SITE TABLE CORRECTLY. A DEEPER thaw is not automatically a worse outcome in the
sense the product's other layers use -- it is the loss of the frozen layer that destabilises
foundations, and beyond the point of total loss the number stops moving. Do not read this
column as a hazard score; it is evidence for choosing one.

Usage:
    python scripts/check_thawdepth_nature.py [--limit N]
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

VAR = "thawdepth"
LAYER_ID = "permafrost-isimip3b_thawdepth_annual"

#: GUARDRAILS 12 reference sites. Three permafrost classes plus deliberate controls -- the
#: controls are how the non-permafrost encoding gets measured rather than inferred.
REFERENCE_SITES = [
    # continuous permafrost
    ("Yakutsk, RU (continuous)",        62.0,  129.7),
    ("Tiksi, RU (continuous)",          71.6,  128.9),
    ("Norilsk, RU (continuous)",        69.3,   88.2),
    ("Deadhorse, AK (continuous)",      70.2, -148.5),
    ("Inuvik, CA (continuous)",         68.4, -133.7),
    ("Resolute, CA (continuous)",       74.7,  -94.8),
    ("Longyearbyen, SJ (continuous)",   78.2,   15.6),
    # discontinuous / sporadic
    ("Fairbanks, AK (discontinuous)",   64.8, -147.7),
    ("Yellowknife, CA (discontinuous)", 62.5, -114.4),
    ("Churchill, CA (discontinuous)",   58.8,  -94.2),
    ("Salekhard, RU (discontinuous)",   66.5,   66.6),
    # alpine
    ("Golmud, CN (Tibetan plateau)",    36.4,   94.9),
    ("Nagqu, CN (Tibetan plateau)",     31.5,   92.1),
    # permafrost-free controls -- these measure the encoding, not the hazard
    ("Stockholm, SE (seasonal frost)",  59.3,   18.1),
    ("Paris, FR (no permafrost)",       48.9,    2.4),
    ("Nairobi, KE (no permafrost)",     -1.3,   36.8),
]


def parse_name(fpath):
    """(model, gcm, scenario, member) read from the END of the filename.

    OutputData filenames carry no leading publication token, but parsing from the end costs
    nothing and is the only form that survives a publication that adds one.
    """
    p = os.path.basename(fpath).rsplit(".", 1)[0].split("_")
    info = dict(model=p[-11], gcm=p[-10], forcing=p[-9], scenario=p[-8], soc=p[-7],
                sens=p[-6], variable=p[-5], cadence=p[-3])
    info["member"] = f"{info['model']}_{info['gcm']}"
    if info["variable"] != VAR:
        raise ValueError(f"{os.path.basename(fpath)}: parsed variable "
                         f"{info['variable']!r} != {VAR!r} -- filename grammar changed")
    return info


def decode_years(ds, fpath):
    """Years per record, taken from the filename span and CHECKED, never trusted blind."""
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
                    help="check only the first N files (iteration only -- a sampled scan "
                         "can miss the one model whose units differ)")
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
    print(f"thawdepth raw value check -- {len(files)} files")
    print("=" * 78)

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    print(f"grid: lat {lats.min():.2f}..{lats.max():.2f} ({len(lats)}), "
          f"lon {lons.min():.2f}..{lons.max():.2f} ({len(lons)})")
    sites = [(name, *site_indices(lats, lons, la, lo), la, lo)
             for name, la, lo in REFERENCE_SITES]

    rows = []
    attr_sets = defaultdict(set)
    attr_by_model = defaultdict(lambda: defaultdict(set))
    finite_mask = {}
    site_base = defaultdict(dict)      # site -> tag -> 2020s mean
    site_late = defaultdict(dict)      # site -> tag -> 2090s mean
    per_model_vals = defaultdict(list)  # model -> subsampled finite values

    for i, f in enumerate(files, 1):
        info = parse_name(f)
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[VAR]
            for k in ("units", "long_name", "standard_name", "cell_methods", "comment"):
                v = str(da.attrs.get(k, "<absent>"))
                attr_sets[k].add(v)
                attr_by_model[info["model"]][k].add(v)
            attr_sets["time_units"].add(str(ds.time.attrs.get("units", "<absent>")))
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

        vmax = float(fin.max())
        # "Pinned at the column bottom": the share of finite cell-years sitting at the
        # member's own maximum. A physical field should not concentrate there; a
        # column-limited one does, and that concentration IS the censoring.
        #
        # The tolerance is RELATIVE and deliberately loose. An absolute 1e-4 reported
        # lpjml's at-max share as 0.0% while its MEDIAN sat at the ceiling, because that
        # model's maximum is 13.001 and the pinned mass is at exactly 13.000 -- a
        # 1-millimetre gap that inverted the finding. Any ceiling test on a float field
        # needs slack proportional to the ceiling.
        tol = max(1e-3, 1e-4 * abs(vmax))
        ceiling = vmax - tol
        at_max = float((fin >= ceiling).mean())
        exact0 = float((fin == 0.0).mean())

        valid = np.isfinite(vals)
        mask2d = valid.any(axis=0)
        time_invariant = bool(np.array_equal(valid,
                                             np.broadcast_to(mask2d, valid.shape)))
        finite_mask.setdefault(info["member"], np.zeros(mask2d.shape, bool))
        finite_mask[info["member"]] |= mask2d

        base = (yrs >= 2020) & (yrs <= 2029)
        late = (yrs >= 2090)
        # Domain proxy per decade: cells whose thaw stays strictly below the member ceiling.
        dom_base = int((np.nanmax(vals[base], axis=0) < ceiling).sum())
        dom_late = int((np.nanmax(vals[late], axis=0) < ceiling).sum())

        tag = f"{info['member']}|{info['scenario']}"
        for name, ilat, ilon, _, _ in sites:
            v = vals[:, ilat, ilon]
            vb, vl = v[base], v[late]
            site_base[name][tag] = (float(np.nanmean(vb))
                                    if np.isfinite(vb).any() else float("nan"))
            site_late[name][tag] = (float(np.nanmean(vl))
                                    if np.isfinite(vl).any() else float("nan"))

        step = max(1, fin.size // 200_000)
        per_model_vals[info["model"]].append(fin[::step])

        rows.append(dict(
            member=info["member"], model=info["model"], gcm=info["gcm"],
            scenario=info["scenario"],
            n_unique=int(np.unique(fin[::step]).size), vmin=float(fin.min()), vmax=vmax,
            ceiling=ceiling,
            median=float(np.median(fin[::step])), at_max_frac=at_max, exact0=exact0,
            finite_cells=int(mask2d.sum()), finite_frac=float(mask2d.mean()),
            time_invariant_mask=time_invariant,
            domain_cells_2020s=dom_base, domain_cells_2090s=dom_late,
            mean_2020s=float(np.nanmean(vals[base])),
            mean_2090s=float(np.nanmean(vals[late])),
            n_years=int(yrs.size), y0=int(yrs.min()), y1=int(yrs.max()),
        ))
        print(f"  [{i}/{len(files)}] {info['member']:<26} {info['scenario']:<7} "
              f"range=[{fin.min():.4g},{vmax:.4g}] med={np.median(fin[::step]):.4g} "
              f"at_max={at_max:.1%} zero={exact0:.1%}", flush=True)

    # ---- 1. metadata ---------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("1. METADATA (per-member, not assumed uniform)")
    print("-" * 78)
    for k, v in attr_sets.items():
        flag = "" if len(v) == 1 else "   <-- DIVERGES ACROSS MEMBERS"
        print(f"  {k:<22} {sorted(v)}{flag}")
    print("\n  per model (a divergence here is the ensemble being a mixture):")
    for model in sorted(attr_by_model):
        a = attr_by_model[model]
        print(f"    {model:<18} units={sorted(a['units'])} "
              f"long_name={sorted(a['long_name'])}")
        print(f"    {'':<18} cell_methods={sorted(a['cell_methods'])} "
              f"standard_name={sorted(a['standard_name'])}")
    print("\n  NOTE: `cell_methods` is the field that says whether the annual value is a "
          "MAXIMUM (active-layer thickness) or a MEAN. If it is absent, the aggregation is "
          "UNDOCUMENTED and must be recorded as such -- do not assume maximum.")

    # ---- 2. data nature, per model -------------------------------------------- #
    print("\n" + "-" * 78)
    print("2. DATA NATURE PER MODEL (GUARDRAILS 9 -- never pooled before compared)")
    print("-" * 78)
    print(f"  {'model':<18} {'min':>8} {'p50':>8} {'p90':>8} {'p99':>8} {'max':>8} "
          f"{'zero%':>7} {'at-max%':>8}")
    model_max = {}
    for model in sorted(per_model_vals):
        v = np.concatenate(per_model_vals[model])
        mx = float(v.max())
        model_max[model] = mx
        print(f"  {model:<18} {v.min():>8.3f} {np.percentile(v, 50):>8.3f} "
              f"{np.percentile(v, 90):>8.3f} {np.percentile(v, 99):>8.3f} {mx:>8.3f} "
              f"{(v == 0).mean():>6.1%} {np.isclose(v, mx, atol=1e-4).mean():>7.1%}")
    mx_vals = np.array(list(model_max.values()))
    if mx_vals.size > 1 and mx_vals.max() / max(mx_vals.min(), 1e-9) > 5:
        print("\n  <-- STOP: per-model maxima differ by more than 5x. Check units before "
              "pooling: a metres/millimetres mixture reads as a huge inter-model spread "
              "and would be published as structural uncertainty.")

    # ---- 3. saturation / censoring -------------------------------------------- #
    print("\n" + "-" * 78)
    print("3. SATURATION (does the value stop moving once the column is thawed?)")
    print("-" * 78)
    print(f"  {'model':<18} {'scenario':<8} {'mean 2020s':>11} {'mean 2090s':>11} "
          f"{'at-max share':>13}")
    for model in sorted({r["model"] for r in rows}):
        for scen in sorted({r["scenario"] for r in rows}):
            sel = [r for r in rows if r["model"] == model and r["scenario"] == scen]
            if not sel:
                continue
            print(f"  {model:<18} {scen:<8} "
                  f"{np.mean([r['mean_2020s'] for r in sel]):>11.4f} "
                  f"{np.mean([r['mean_2090s'] for r in sel]):>11.4f} "
                  f"{np.mean([r['at_max_frac'] for r in sel]):>12.1%}")
    print("\n  Reading: a rising pair is thaw deepening. A large at-max share means cells "
          "have hit the bottom of the modelled soil column, where BOTH slopes go to ~0 and "
          "mean 'no permafrost left', not 'no trend' -- the trend-censoring problem, and "
          "the reason the two-slope disagreement rule cannot rescue this layer either.")

    # ---- 4. domain ------------------------------------------------------------ #
    print("\n" + "-" * 78)
    print("4. DOMAIN (no permafrost mask ships -- this is the only field-derived proxy)")
    print("-" * 78)
    print(f"  {'model':<18} {'scenario':<8} {'cells<max 2020s':>16} "
          f"{'cells<max 2090s':>16} {'change':>9}")
    for model in sorted({r["model"] for r in rows}):
        for scen in sorted({r["scenario"] for r in rows}):
            sel = [r for r in rows if r["model"] == model and r["scenario"] == scen]
            if not sel:
                continue
            b = np.mean([r["domain_cells_2020s"] for r in sel])
            l = np.mean([r["domain_cells_2090s"] for r in sel])
            print(f"  {model:<18} {scen:<8} {b:>16,.0f} {l:>16,.0f} "
                  f"{(l - b) / max(b, 1):>8.1%}")
    print("\n  'cells<max' = cells whose thaw stays strictly below that member's maximum, "
          "i.e. something frozen remains. If this count is close to the finite extent the "
          "proxy is NOT separating permafrost from seasonally frozen ground and a real mask "
          "must come from elsewhere.")

    # ---- 5. mask -------------------------------------------------------------- #
    print("\n" + "-" * 78)
    print("5. MASK -- finite extent vs land, and member agreement")
    print("-" * 78)
    bad_mask = [r for r in rows if not r["time_invariant_mask"]]
    print(f"  time-invariant per file: {len(rows) - len(bad_mask)}/{len(rows)}")
    for r in bad_mask:
        print(f"  <-- TIME-VARYING mask: {r['member']} {r['scenario']}")
    fin_cells = sorted({r["finite_cells"] for r in rows})
    print(f"  finite cells per file: {fin_cells[:6]}"
          + (" ..." if len(fin_cells) > 6 else "")
          + f"  (max {100.0 * max(fin_cells) / (len(lats) * len(lons)):.1f}% of the grid)")

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
        print(f"  ISIMIP3b land cells: {int(land.sum()):,}")
        print(f"  finite-ever union  : {int(fin_union.sum()):,}  "
              f"on land {int((fin_union & land).sum()):,}  "
              f"OFF LAND {int((fin_union & ~land).sum()):,}")
        print(f"  finite extent as a share of the GRID: {fin_union.mean():.1%}")
        stack = np.array([finite_mask[k] for k in sorted(finite_mask)])
        nmem = stack.sum(axis=0)
        print(f"\n  coverage tiers (cells finite in exactly K of {len(finite_mask)} "
              f"members) -- this is what a publication mask would be defined on:")
        for k in range(1, len(finite_mask) + 1):
            c = int((nmem == k).sum())
            if c:
                print(f"    {k:>2} member(s): {c:>8,}")

    # ---- 6. reference sites --------------------------------------------------- #
    print("\n" + "-" * 78)
    print("6. REFERENCE SITES (GUARDRAILS 12)")
    print("-" * 78)
    # PER MODEL, never pooled. Each model pins its permafrost-free cells at its OWN soil
    # column depth, so an across-model mean at a site is an average of three different
    # ceilings and means nothing -- the first version of this table printed 16.900 m for
    # Nairobi, which is (61.4x6 + 13x15 + 3x15)/36 and not a depth anything reaches.
    models = sorted({r["model"] for r in rows})
    stops = []
    for model in models:
        ceil = np.mean([r["ceiling"] for r in rows if r["model"] == model])
        print(f"\n  {model}  (soil column ceiling ~{ceil:.3f} m)")
        print(f"  {'site':<32} {'lat':>6} {'lon':>7} {'NaN':>4} "
              f"{'2020s':>9} {'2090s':>9} {'delta':>8}  at ceiling?")
        for name, ilat, ilon, la, lo in sites:
            keys = [k for k in site_base[name] if k.rsplit("_", 1)[0] == model]
            b = np.array([site_base[name][k] for k in keys], dtype=float)
            l = np.array([site_late[name][k] for k in keys], dtype=float)
            n_nan = int(np.isnan(b).sum())
            bf, lf = b[np.isfinite(b)], l[np.isfinite(l)]
            bm = float(bf.mean()) if bf.size else float("nan")
            lm_ = float(lf.mean()) if lf.size else float("nan")
            pin = ""
            if bf.size == 0:
                pin = "NaN in every member"
                if "no permafrost" not in name:
                    stops.append(f"{model}/{name}")
            else:
                at_b, at_l = bm >= ceil, lm_ >= ceil
                pin = ("pinned throughout" if at_b and at_l else
                       "REACHES ceiling by 2090s" if at_l else
                       "below ceiling" if not at_b else "")
            print(f"  {name:<32} {la:>6.1f} {lo:>7.1f} {n_nan:>4} "
                  f"{bm:>9.3f} {lm_:>9.3f} {lm_ - bm:>8.3f}  {pin}")
    print("\n  Reading: a continuous-permafrost site should sit WELL BELOW its model's "
          "ceiling and deepen. A site pinned at the ceiling in the 2020s has no permafrost "
          "in that model at all, and a site that REACHES the ceiling has lost it within the "
          "record -- which is the hazard, and is exactly where the raw depth stops moving.")

    # ---- 7. domain in AREA terms, and member agreement ------------------------ #
    #
    # A cell count is not comparable to anything published: 0.5 deg cells shrink with
    # latitude and this hazard lives at 60-80 deg N, where a cell is a third the area of an
    # equatorial one. Area is the only form in which a modelled permafrost domain can be
    # checked against the outside world at all.
    print("\n" + "-" * 78)
    print("7. DOMAIN AREA and ENSEMBLE AGREEMENT (2020s baseline decade)")
    print("-" * 78)
    R = 6371.0088
    dlat = float(abs(lats[1] - lats[0]))
    dlon = float(abs(lons[1] - lons[0]))
    cell_km2 = (np.deg2rad(dlat) * R) * (np.deg2rad(dlon) * R
                                         * np.cos(np.deg2rad(lats)))[:, None]
    cell_km2 = np.broadcast_to(cell_km2, (len(lats), len(lons)))

    present_stack, present_keys = [], []
    for f in files:
        info = parse_name(f)
        if info["scenario"] != "ssp370":      # one scenario: the 2020s baseline is shared
            continue
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[VAR]
            yrs, _ = decode_years(ds, f)
            v = da.values.astype("float32")
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))
        if fill is not None:
            v = np.where(np.isclose(v, np.float32(fill), rtol=1e-6), np.nan, v)
        v[~np.isfinite(v)] = np.nan
        vmax = float(np.nanmax(v))
        ceil = vmax - max(1e-3, 1e-4 * abs(vmax))
        base = (yrs >= 2020) & (yrs <= 2029)
        present = np.nanmax(v[base], axis=0) < ceil
        present_stack.append(present)
        present_keys.append(info["member"])

    print(f"  {'model':<18} {'members':>8} {'permafrost area 2020s (M km2)':>32}")
    for model in models:
        idx = [i for i, k in enumerate(present_keys) if k.rsplit("_", 1)[0] == model]
        areas = [float((present_stack[i] * cell_km2).sum()) / 1e6 for i in idx]
        print(f"  {model:<18} {len(idx):>8} "
              f"{np.mean(areas):>21.2f}  (range {min(areas):.2f}-{max(areas):.2f})")
    stack = np.array(present_stack)
    agree = stack.sum(axis=0)
    print(f"\n  ensemble agreement -- cells called permafrost-bearing by exactly K of "
          f"{len(present_stack)} members:")
    for k in range(1, len(present_stack) + 1):
        c = int((agree == k).sum())
        if c:
            a = float(((agree == k) * cell_km2).sum()) / 1e6
            print(f"    {k:>2}: {c:>7,} cells  {a:>7.2f} M km2")
    for k in (1, 6, len(present_stack)):
        a = float(((agree >= k) * cell_km2).sum()) / 1e6
        print(f"  union at >={k:>2} members: {int((agree >= k).sum()):>7,} cells  "
              f"{a:>7.2f} M km2")
    print("\n  Compare against the observed Northern Hemisphere permafrost area. Any model "
          "whose 2020s domain is far from it is not a candidate for a site-level product, "
          "whatever its ensemble contribution.")

    out = raw_dir / "value_check.json"
    with open(out, "w") as fh:
        json.dump(dict(
            n_files=len(rows),
            attrs={k: sorted(v) for k, v in attr_sets.items()},
            per_model_max=model_max,
            rows=rows,
            sites_2020s={n: site_base[n] for n in site_base},
            sites_2090s={n: site_late[n] for n in site_late},
        ), fh, indent=2)
    print(f"\nwrote {out}")

    if stops:
        print(f"\nSTOP -- {len(stops)} permafrost reference site(s) are NaN in every "
              f"member: {stops}. Investigate before writing any output (GUARDRAILS 12).")
        return 4
    print("\nNo reference-site STOP condition. The measured numbers above -- not this "
          "script's prose -- govern the processor's statistic, mask and direction.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
