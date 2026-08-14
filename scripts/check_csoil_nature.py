"""Measure the data nature, domain and site behaviour of raw `csoil-total` before processing.

Soil organic carbon is a STOCK in kg C m-2 -- a physical field, not a 0/1 exposure flag with
the framing already decided by the publisher. Everything downstream (which statistic branch,
which mask, whether to normalise, which slope to read) follows from what is measured here,
and NONE of it is inherited from the 2026-07-25 run: that run had three models, this one has
four, and a fourth model can break a comparability finding that held for three.

Six checks, run together because they need the same open files:

  * GUARDRAILS 9 -- WHAT the values are, PER MODEL, never pooled before they are compared.
    Four independent land models each implement their own soil carbon scheme and vertical
    extent. A unit divergence, or a magnitude divergence large enough to make one model
    dominate the pool, breaks no contract check and silently turns the ensemble into a
    mixture of incompatible fields. The 2026-07-25 finding -- medians 5.76 / 10.28 / 7.67,
    "comparable, so no normalization" -- is treated here as a HYPOTHESIS TO RE-TEST with
    LPJmL added, not as a settled fact.

  * THE ENCODING OF A NON-SOIL CELL -- NaN, exact 0, or a fill value. This decides the mask,
    and it is not guessable: LPJmL's own global attributes say "pools are zero for cell
    fractions covered by inland waterbodies" and give reference_area as "continental area
    (including inland water bodies)", so a zero in that model may mean water rather than
    carbon-free soil. `cropfailure` is the standing warning here -- its publisher zero-fills
    the entire globe, so `isfinite` carried no footprint at all and using it as a mask would
    have put ocean into the percentile baseline.

  * MULTIMODALITY. OUTPUT-SPEC's fourth branch (`pooled_mean_multimodel`) exists because the
    median SELECTS the larger cluster when members separate, and jumps when the balance tips.
    With four models at potentially different levels this must be tested, not assumed away.
    Reported as per-model 2020s distributions side by side plus the gap between adjacent
    model medians relative to the within-model spread. `permafrost-3b` needed this branch at
    a separation of 0.035 / 0.046 / 0.951; a comparable separation here would too.

  * ZERO INFLATION. The third branch (`pooled_mean_zero_inflated`) applies when the decade
    pool is degenerate at zero. Soil carbon should not be -- deserts are low, not absent --
    but the share of exact zeros is measured per model rather than assumed, because that is
    exactly the check that would have caught `let` and `cropfailure` earlier.

  * BOUNDS / CENSORING. A field pinned at a bound makes BOTH slopes go to ~0 and AGREE, so
    the dual-slope disagreement rule gives no warning (OUTPUT-SPEC, and `heatwave-3b`). Soil
    carbon has a floor at 0 and no physical ceiling, so the floor is what is measured: the
    share of cells at exactly 0 by model. If it is large the layer needs a `sparsity_caveat`.

  * GUARDRAILS 12 -- WHERE the values are. A layer can pass every contract check and be
    meaningless; the sugarcane withdrawal is the precedent. Soil carbon has an unusually
    strong prior: peatlands and chernozems are the highest-carbon soils on the planet and
    hot deserts are near the bottom, so a model that does not reproduce that ordering is not
    reporting soil carbon. Controls (ocean, ice sheet) measure the encoding question above
    more directly than any global statistic.

READ THE SITE TABLE WITH THE DIRECTION IN MIND. High soil carbon is GOOD -- this layer is
`higher_is_better` and its percentile is inverted downstream. A high number at Ukraine or
the Flow Country is the layer working, not a hazard.

Usage:
    python scripts/check_csoil_nature.py [--limit N]
"""

import argparse
import glob
import os
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr

VAR = "csoil-total"
LAYER_ID = "soilcarbon_csoil_annual"
BASELINE = (2020, 2029)
LATE = (2090, 2099)

#: GUARDRAILS 12 reference sites. The prior is strong and well documented in soil science,
#: which is what makes this check able to fail: organic soils and chernozems must sit far
#: above hot deserts, and the controls must show what a non-soil cell contains.
REFERENCE_SITES = [
    # -- organic / peat soils: the global maxima -------------------------------------
    ("W Siberian Lowland peat, RU",   61.3,   73.4, "high"),
    ("Hudson Bay Lowlands, CA",       51.3,  -80.6, "high"),
    ("Cuvette Centrale peat, CD",     -1.0,   18.5, "high"),
    ("Flow Country blanket bog, GB",  58.4,   -3.9, "high"),
    ("Sodankyla peat, FI",            67.4,   26.6, "high"),
    ("North Slope tundra, AK",        70.2, -148.5, "high"),
    # -- chernozem / mollisol: the agricultural high-carbon soils --------------------
    ("Ukraine chernozem",             49.0,   32.0, "high"),
    ("Iowa mollisol, US",             42.0,  -93.6, "high"),
    ("Pampas, AR",                   -34.0,  -61.0, "mid"),
    # -- humid tropical: high productivity, moderate soil carbon --------------------
    ("Amazon terra firme, BR",        -3.1,  -60.0, "mid"),
    # -- hot desert: the global minima ---------------------------------------------
    ("Tamanrasset, Sahara, DZ",       22.8,    5.5, "low"),
    ("Atacama, CL",                  -24.5,  -69.3, "low"),
    ("Rub al Khali, SA",              20.0,   50.0, "low"),
    ("Taklamakan, CN",                39.0,   83.0, "low"),
    # -- controls: these measure the ENCODING, not the hazard -----------------------
    ("Greenland ice sheet",           72.0,  -40.0, "control"),
    ("Pacific ocean (5N,150W)",        5.0, -150.0, "control"),
]


def parse_name(fpath):
    """(model, gcm, scenario, member) read from the END of the filename.

    Model and GCM tokens contain hyphens, never underscores, so a from-the-end index is
    stable; a fixed forward index breaks on `lpjml5-7-10-fire` and `mc2-usfs-r87g5c1`.
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
                    help="check only the first N files (ITERATION ONLY -- a sampled scan "
                         "can step over the one model whose units or encoding differ)")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    raw_dir = root / "data" / "raw" / LAYER_ID
    files = sorted(glob.glob(str(raw_dir / f"*_{VAR}_global_annual_*.nc")))
    if not files:
        print(f"ERROR: no {VAR} files in {raw_dir}")
        return 2
    if args.limit:
        files = files[:args.limit]

    print("=" * 84)
    print(f"csoil-total raw value check -- {len(files)} files")
    print("=" * 84)

    with xr.open_dataset(files[0], decode_times=False) as ds0:
        lats, lons = ds0["lat"].values, ds0["lon"].values
    print(f"grid: lat {lats.min():.2f}..{lats.max():.2f} ({len(lats)}), "
          f"lon {lons.min():.2f}..{lons.max():.2f} ({len(lons)})")
    sites = [(name, *site_indices(lats, lons, la, lo), la, lo, klass)
             for name, la, lo, klass in REFERENCE_SITES]

    attr_by_model = defaultdict(lambda: defaultdict(set))
    stats = defaultdict(list)                 # model -> per-file summary dicts
    base_vals = defaultdict(list)             # model -> 2020s spatial means (flat, land)
    site_base = defaultdict(lambda: defaultdict(list))   # site -> model -> 2020s values
    site_late = defaultdict(lambda: defaultdict(list))   # site -> model -> 2090s values
    finite_counts = {}

    for i, f in enumerate(files, 1):
        info = parse_name(f)
        m = info["model"]
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[VAR]
            for k in ("units", "long_name", "standard_name", "cell_methods", "comment"):
                attr_by_model[m][k].add(str(da.attrs.get(k, "<absent>")))
            attr_by_model[m]["time_units"].add(str(ds.time.attrs.get("units", "<absent>")))
            attr_by_model[m]["calendar"].add(str(ds.time.attrs.get("calendar", "<absent>")))
            yrs, conv = decode_years(ds, f)
            attr_by_model[m]["time_step_convention"].add(conv)
            vals = da.values.astype("float32")
            fill = da.attrs.get("_FillValue", da.attrs.get("missing_value", None))

        raw_nan = int(np.sum(~np.isfinite(vals)))
        if fill is not None:
            hit = np.isclose(vals, np.float32(fill), rtol=1e-6)
            vals = np.where(hit, np.nan, vals)
        vals[~np.isfinite(vals)] = np.nan

        b0, b1 = BASELINE
        l0, l1 = LATE
        bsel = (yrs >= b0) & (yrs <= b1)
        lsel = (yrs >= l0) & (yrs <= l1)
        with warnings.catch_warnings():
            # Ocean cells are all-NaN through the window by construction; nanmean says so
            # once per cell and drowns the report. The NaN it returns is the right answer.
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            bmean = np.nanmean(vals[bsel], axis=0)      # (lat, lon) 2020s mean
            lmean = np.nanmean(vals[lsel], axis=0)

        fin = bmean[np.isfinite(bmean)]
        if fin.size == 0:
            print(f"  [{i}/{len(files)}] {info['member']:<32} ALL-NaN MEMBER -- STOP")
            return 3

        stats[m].append(dict(
            member=info["member"], scenario=info["scenario"],
            n_finite=int(fin.size),
            frac_finite=float(fin.size) / bmean.size,
            frac_zero=float(np.mean(fin == 0.0)),
            p50=float(np.median(fin)), p05=float(np.percentile(fin, 5)),
            p95=float(np.percentile(fin, 95)), vmax=float(fin.max()),
            raw_nan=raw_nan,
        ))
        if info["scenario"] == "ssp126":
            base_vals[m].append(fin[::7])           # subsample: distribution shape only
        finite_counts[info["member"]] = int(fin.size)

        for name, ii, jj, _, _, klass in sites:
            site_base[name][m].append(float(bmean[ii, jj]))
            site_late[name][m].append(float(lmean[ii, jj]))

        if i % 10 == 0 or i == len(files):
            print(f"  ... {i}/{len(files)} files read")

    models = sorted(stats)

    # -- GUARDRAILS 9: declared metadata, per model ---------------------------------
    print("\n" + "=" * 84)
    print("1. DECLARED METADATA, PER MODEL (GUARDRAILS 9 -- divergence here is fatal)")
    print("=" * 84)
    for m in models:
        print(f"\n{m}")
        for k in ("units", "long_name", "cell_methods", "comment", "calendar",
                  "time_units", "time_step_convention"):
            v = sorted(attr_by_model[m][k])
            print(f"   {k:<22} {' | '.join(x[:70] for x in v)}")
    all_units = {u for m in models for u in attr_by_model[m]["units"]}
    print(f"\n   UNITS ACROSS ALL MODELS: {sorted(all_units)}"
          + ("   <-- OK, single unit" if len(all_units) == 1
             else "   <-- DIVERGENT, DO NOT POOL"))

    # -- magnitude comparability ----------------------------------------------------
    print("\n" + "=" * 84)
    print("2. MAGNITUDE PER MODEL, 2020s (re-testing the 3-model 'comparable' finding)")
    print("=" * 84)
    print(f"{'model':<20} {'members':>7} {'finite%':>8} {'zero%':>7} "
          f"{'p05':>8} {'p50':>8} {'p95':>8} {'max':>8}")
    med_by_model = {}
    for m in models:
        rows = stats[m]
        n_mem = len({r['member'] for r in rows})
        med = float(np.median([r["p50"] for r in rows]))
        med_by_model[m] = med
        print(f"{m:<20} {n_mem:>7} {np.mean([r['frac_finite'] for r in rows])*100:>7.2f}% "
              f"{np.mean([r['frac_zero'] for r in rows])*100:>6.2f}% "
              f"{np.median([r['p05'] for r in rows]):>8.3f} {med:>8.3f} "
              f"{np.median([r['p95'] for r in rows]):>8.3f} "
              f"{max(r['vmax'] for r in rows):>8.1f}")
    lo, hi = min(med_by_model.values()), max(med_by_model.values())
    print(f"\n   median spread across models: {lo:.3f} .. {hi:.3f}  = {hi/max(lo,1e-9):.2f}x")
    print("   A large ratio does NOT automatically mean normalise -- it means the pooled"
          "\n   sample may be multimodal. See check 4.")

    # -- encoding of a non-soil cell ------------------------------------------------
    print("\n" + "=" * 84)
    print("3. WHAT A NON-SOIL CELL CONTAINS (decides the mask)")
    print("=" * 84)
    ncell = len(lats) * len(lons)
    for m in models:
        rows = stats[m]
        ff = np.mean([r["frac_finite"] for r in rows])
        fz = np.mean([r["frac_zero"] for r in rows])
        print(f"{m:<20} finite {ff*100:6.2f}% of {ncell} cells "
              f"({np.mean([r['n_finite'] for r in rows]):.0f} cells), "
              f"exact-zero {fz*100:5.2f}% of finite")
    print("\n   Land is ~27-29% of a 360x720 grid. A finite share MUCH larger than that"
          "\n   means the publisher fills non-land, and `isfinite` is NOT a footprint"
          "\n   (the `cropfailure` failure mode). Compare against the site controls below.")

    # -- multimodality --------------------------------------------------------------
    print("\n" + "=" * 84)
    print("4. IS THE POOLED 2020s SAMPLE UNIMODAL? (decides branch 4)")
    print("=" * 84)
    order = sorted(models, key=lambda m: med_by_model[m])
    print(f"{'model':<20} {'p50':>8} {'IQR':>16}   gap to next model / own IQR")
    for k, m in enumerate(order):
        v = np.concatenate(base_vals[m]) if base_vals[m] else np.array([np.nan])
        q1, q3 = np.percentile(v, [25, 75])
        iqr = q3 - q1
        if k + 1 < len(order):
            gap = med_by_model[order[k + 1]] - med_by_model[m]
            rel = gap / iqr if iqr > 0 else np.inf
            note = f"gap {gap:+.3f}  = {rel:.2f} x IQR"
            if rel >= 1.0:
                note += "   <-- SEPARATED"
        else:
            note = ""
        print(f"{m:<20} {med_by_model[m]:>8.3f} {q1:>7.3f}-{q3:<8.3f} {note}")
    print("\n   Clusters separated by >~1 IQR are the `permafrost-3b` regime, where the"
          "\n   pooled MEDIAN selects a cluster instead of summarising. If nothing is"
          "\n   separated, `pooled_median` stands and branch 4 is declined ON MEASUREMENT.")

    # -- GUARDRAILS 12 --------------------------------------------------------------
    print("\n" + "=" * 84)
    print("5. GUARDRAILS 12 -- REFERENCE SITES, 2020s mean kg C m-2 (ensemble mean/model)")
    print("=" * 84)
    hdr = f"{'site':<30} {'class':<8}" + "".join(f"{m[:11]:>12}" for m in models)
    print(hdr)
    print("-" * len(hdr))
    for name, _, _, la, lo, klass in sites:
        cells = "".join(
            f"{np.nanmean(site_base[name][m]):>12.2f}" if site_base[name][m]
            and np.isfinite(np.nanmean(site_base[name][m])) else f"{'NaN':>12}"
            for m in models)
        print(f"{name:<30} {klass:<8}{cells}")
    print("\n   EXPECTED: high >> mid >> low, and controls NaN. If a hot desert reads at or"
          "\n   above a peatland in some model, that model is not reporting soil carbon and"
          "\n   the layer must not ship on it (the sugarcane precedent).")

    # -- direction of change --------------------------------------------------------
    print("\n" + "=" * 84)
    print("6. 2090s MINUS 2020s AT THE REFERENCE SITES (ensemble mean, all scenarios)")
    print("=" * 84)
    hdr = f"{'site':<30} {'class':<8}" + "".join(f"{m[:11]:>12}" for m in models)
    print(hdr)
    print("-" * len(hdr))
    for name, _, _, la, lo, klass in sites:
        cells = ""
        for m in models:
            b = np.nanmean(site_base[name][m]) if site_base[name][m] else np.nan
            l = np.nanmean(site_late[name][m]) if site_late[name][m] else np.nan
            cells += f"{l-b:>+12.2f}" if np.isfinite(b) and np.isfinite(l) else f"{'NaN':>12}"
        print(f"{name:<30} {klass:<8}{cells}")
    print("\n   Sign matters more than magnitude here: LOSS is the risk, and this layer is"
          "\n   higher_is_better with an INVERTED percentile downstream.")

    print("\n" + "=" * 84)
    print("Record these numbers in the processor docstring and in "
          "config/isimip_search_catalog.yaml")
    print("search_results.soil_carbon before writing or trusting the processor.")
    print("=" * 84)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
