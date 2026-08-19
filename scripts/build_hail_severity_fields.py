#!/usr/bin/env python3
"""Gridded hail SEVERITY response from the Nature 2026 deposit: intensity and mean size.

READ `scripts/measure_hail_nature2026.py` FIRST. Its docstring carries the semantics of the
deposit -- the event index, the `w_perc` common-support rule, what the zeros mean -- and the
evidence gate that established them. This script assumes all of it and adds nothing to it.

WHAT THIS PRODUCES, AND WHY THESE QUANTITIES
    Maximum stone diameter is not the hazard. What damages an asset is how hard the stones
    hit and how big they typically are, so the fields here are:

      mean_ke        mean kinetic energy of the stones that actually reach the ground (J).
                     Pure intensity: conditioned on a stone existing.
      mean_diameter  mean diameter of the stones that reach the ground (mm). Pure size.
      ake            mean kinetic energy per SEEDED embryo (J) -- the authors' own
                     `severity_calc_n`, which divides by the embryo count and therefore
                     counts melted embryos as zero. It BLENDS production and intensity and
                     is the quantity behind the paper's "damage potential" headline.
      production     share of seeded embryos that reach the ground as hail. The production
                     term that `ake` folds in, published separately so the two are legible.

WEIGHTING -- IT CHANGES THE SIGN, SO BOTH ARE PUBLISHED
    A cell holds many storms and each storm drops a different number of stones. Averaging
    over STONES and averaging over STORMS are different questions and on this data they
    disagree. Measured, US Plains ssp585: pooled over landed stones the median diameter rises
    37.3 -> 42.6 mm, while the storm-weighted mean of each storm's own mean FALLS. Warming
    lets more storms produce hail at all, and those newly-producing storms are weak ones, so
    they drag the per-storm average down while the stones themselves are not smaller.

      mean_diameter / mean_ke / production / ake   POOLED over stones (and embryos). Primary.
      storm_mean_diameter / storm_mean_ke          each storm counts once. Secondary.

    Pooled is the right default for "mean size of stones at this location"; storm-weighted
    answers "what does a typical storm here drop". Quoting either without naming it is how a
    sign error reaches a customer.

    `mean_ke` x `production` == `ake` EXACTLY, per event and per cell, for the pooled
    fields -- all three are ratios of the same summed numerators and denominators, and the
    build asserts it. Publishing all three is deliberate: they move by different amounts and
    a single blended number hides which mechanism moved.

FREQUENCY IS NOT HERE AND CANNOT BE
    Storm frequency is not derivable from this deposit at any resolution. The future runs
    re-simulate THE SAME satellite-detected storms under perturbed environments; there is no
    occurrence model, no time-area exposure and no detection probability, so the event count
    in a cell measures GPM sampling, not hail climatology. `production` is a per-storm
    productivity term and must never be relabelled as frequency. A frequency field needs a
    different source entirely (an ERA5/AR-CHaMo style occurrence model, or SPC reports over
    CONUS).

ONE MODEL, MAPPED; THREE MODELS, GLOBAL ONLY
    The spatial fields are built from the 20 EC-Earth3 realizations, because only those
    carry a recoverable event index (per-member sounding profiles). The MPI-ESM1-2-LR and
    NorESM2-LM runs are deposited as ensemble MEANS whose `sample` coordinate is a contiguous
    renumbering with no profile file and a different length per model, so they cannot be
    placed in space at all. Their structural spread is real and large -- mean kinetic energy
    rises +28.4% (MPI) to +49.1% (NorESM2/EC-Earth3) globally -- and is reported by
    `measure_hail_nature2026.py --reproduction-only`, but it CANNOT be mapped. Every
    magnitude in these fields is one model's.

AGGREGATION
    Per member: pair events against the historical run on common support, compute the four
    per-event quantities, bin by the event's own lat/lon, and take the unweighted mean over
    events in the cell. Historical and future are computed on the SAME event set within each
    member, so every change is paired. Across members the published value is the median and
    the range is [min, max] -- an initial-condition spread, not a confidence interval.

    Cells are reported with their event count. `--min-events` masks thin cells rather than
    smoothing them; at 5 degrees only 136 cells hold 30+ events while at 10 degrees 117 do,
    holding 92% of all events, which is why 10 degrees is the default.

USAGE
    python3 scripts/build_hail_severity_fields.py                      # 10 deg, n>=30
    python3 scripts/build_hail_severity_fields.py --cell 5 --min-events 10
    python3 scripts/build_hail_severity_fields.py --members 1 2 3      # quick pass
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from measure_hail_nature2026 import (  # noqa: E402
    COMMON_W_PERC,
    MEMBER_IDS,
    N_EMBRYO_COMMON,
    TYPES,
    fut_event_index,
    fut_radii,
    hist_event_index,
    hist_radii,
    kinetic_energy_j,
)

REGIONS = {
    "US Plains":        (30, 50, -105, -90),
    "US Southeast":     (25, 37, -95, -75),
    "Argentina/Uruguay": (-40, -20, -70, -53),
    "Southern Brazil":  (-33, -20, -58, -45),
    "Europe":           (40, 55, -5, 30),
    "South Asia":       (15, 35, 68, 92),
    "East China":       (25, 45, 105, 125),
    "Southern Africa":  (-35, -20, 18, 35),
    "East Africa":      (-10, 15, 30, 50),
    "Australia":        (-40, -20, 138, 155),
    "Maritime SE Asia": (-10, 20, 95, 130),
}


def event_quantities(radii_mm: np.ndarray) -> dict[str, np.ndarray]:
    """Per-event severity terms from an (event, r_int, h_int) block of landing sizes.

    Returns SUMS and COUNTS, not just means, because the two ways of averaging to a cell give
    different answers and both are wanted -- see `WEIGHTING` in the module docstring.
    """
    flat = radii_mm.reshape(radii_mm.shape[0], -1)
    produced = flat > 0
    n_landed = produced.sum(axis=1).astype(float)
    ke = kinetic_energy_j(flat)
    sum_d = np.nansum(np.where(produced, flat, 0.0), axis=1)
    sum_ke = np.nansum(np.where(produced, ke, 0.0), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_d = np.where(n_landed > 0, sum_d / np.maximum(n_landed, 1), np.nan)
        mean_ke = np.where(n_landed > 0, sum_ke / np.maximum(n_landed, 1), np.nan)
    return {
        "sum_diameter": sum_d,
        "sum_ke": sum_ke,
        "n_landed": n_landed,
        "n_embryo": np.full(flat.shape[0], float(flat.shape[1])),
        "storm_mean_diameter": mean_d,
        "storm_mean_ke": mean_ke,
    }


def member_fields(raw: Path, work: Path, scenario: str, member: int,
                  hidx: dict[int, np.ndarray], cat: xr.Dataset,
                  names: dict[str, int]) -> dict[str, np.ndarray]:
    """Paired per-event quantities plus each event's location, for one member."""
    fidx = fut_event_index(raw, work, scenario, member)
    keys = ("sum_diameter", "sum_ke", "n_landed", "n_embryo", "storm_mean_diameter", "storm_mean_ke")
    out: dict[str, list[np.ndarray]] = {k: [] for k in
                                        ("lat", "lon", *[f"{s}_{k}" for s in "hf" for k in keys])}
    for t in TYPES:
        hr = hist_radii(raw, work, t)
        fr = fut_radii(raw, work, scenario, member, t)
        keep_h = np.isin(hidx[t], fidx[t])
        keep_f = np.isin(fidx[t], hidx[t])
        events = hidx[t][keep_h]
        hq = event_quantities(hr[keep_h])
        fq = event_quantities(fr[keep_f])
        params = cat[f"para33_Idx{t}"].values
        out["lat"].append(params[names["lat"]][events])
        lon = params[names["lon"]][events]
        out["lon"].append(((lon + 180.0) % 360.0) - 180.0)  # catalogue is 0..360, grids are -180..180
        for k, v in hq.items():
            out[f"h_{k}"].append(v)
        for k, v in fq.items():
            out[f"f_{k}"].append(v)
    return {k: np.concatenate(v) for k, v in out.items()}


def bin_sum(values: np.ndarray, iy: np.ndarray, ix: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Sum of `values` per cell, ignoring NaN."""
    ok = np.isfinite(values)
    total = np.zeros(shape)
    np.add.at(total, (iy[ok], ix[ok]), values[ok])
    return total


def bin_mean(values: np.ndarray, iy: np.ndarray, ix: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Unweighted mean of `values` per cell, ignoring NaN."""
    ok = np.isfinite(values)
    total = np.zeros(shape)
    count = np.zeros(shape)
    np.add.at(total, (iy[ok], ix[ok]), values[ok])
    np.add.at(count, (iy[ok], ix[ok]), 1.0)
    with np.errstate(invalid="ignore"):
        return np.where(count > 0, total / count, np.nan)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=Path("data/raw/hail-nature2026"))
    ap.add_argument("--scenarios", nargs="+", default=["ssp245", "ssp370", "ssp585"])
    ap.add_argument("--members", nargs="+", type=int, default=list(range(1, 21)))
    ap.add_argument("--cell", type=float, default=10.0, help="grid cell size in degrees")
    ap.add_argument("--min-events", type=int, default=30, help="cells with fewer paired events are masked")
    ap.add_argument("--out-dir", type=Path, default=Path("reports/maps/hail-severity"))
    args = ap.parse_args()

    raw = args.raw_dir
    if not raw.is_dir():
        print(f"raw dir not found: {raw}", file=sys.stderr)
        return 2
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cat = xr.open_dataset(raw / "data_SOM5_global_after200_rightMUCAPE_mdf.nc")
    names = {str(n): i for i, n in enumerate(cat["dim2"].values)}

    cell = args.cell
    # int(0.5) == 0, so every sub-degree cell size would write to the same "0deg" name.
    tag = f"{cell:g}".replace(".", "p")
    ny, nx = int(round(180 / cell)), int(round(360 / cell))
    lats = -90 + cell * (np.arange(ny) + 0.5)
    lons = -180 + cell * (np.arange(nx) + 0.5)

    # POOLED fields are ratios of summed numerators and denominators, so the decomposition
    # ake == production * mean_ke holds EXACTLY at cell level as well as per event. STORM
    # fields average each storm's own mean, so a storm counts once regardless of how many
    # stones it drops. They are not interchangeable and on this data they can differ in SIGN.
    POOLED = ("mean_diameter", "mean_ke", "production", "ake")
    STORM = ("storm_mean_diameter", "storm_mean_ke")

    # Empty cells are the normal case on a sparse event population, not a defect.
    warnings.filterwarnings("ignore", message="All-NaN slice encountered")
    warnings.filterwarnings("ignore", message="Mean of empty slice")

    work = Path(tempfile.mkdtemp(prefix="hail_fields_"))
    panels: dict[str, dict[str, list[np.ndarray]]] = {}
    counts: dict[str, np.ndarray] = {}
    try:
        hidx = hist_event_index(raw, work)
        for scenario in args.scenarios:
            if not (raw / f"{scenario}-simulation_results.7z").exists():
                print(f"[skip] {scenario}: archive not present")
                continue
            stack = {f"{side}_{q}": [] for side in ("h", "f") for q in (*POOLED, *STORM)}
            n_cells = None
            for member in args.members:
                d = member_fields(raw, work, scenario, member, hidx, cat, names)
                iy = np.clip(((d["lat"] + 90) / cell).astype(int), 0, ny - 1)
                ix = np.clip(((d["lon"] + 180) / cell).astype(int), 0, nx - 1)
                for side in ("h", "f"):
                    s_d = bin_sum(d[f"{side}_sum_diameter"], iy, ix, (ny, nx))
                    s_ke = bin_sum(d[f"{side}_sum_ke"], iy, ix, (ny, nx))
                    n_l = bin_sum(d[f"{side}_n_landed"], iy, ix, (ny, nx))
                    n_e = bin_sum(d[f"{side}_n_embryo"], iy, ix, (ny, nx))
                    with np.errstate(invalid="ignore", divide="ignore"):
                        stack[f"{side}_mean_diameter"].append(np.where(n_l > 0, s_d / n_l, np.nan))
                        stack[f"{side}_mean_ke"].append(np.where(n_l > 0, s_ke / n_l, np.nan))
                        stack[f"{side}_production"].append(np.where(n_e > 0, n_l / n_e, np.nan))
                        stack[f"{side}_ake"].append(np.where(n_e > 0, s_ke / n_e, np.nan))
                    # Exact by construction for THIS member: all three are ratios of the same
                    # summed numerators and denominators. Assert here, not after the median --
                    # a median is not multiplicative, so median(ake) != median(production) *
                    # median(mean_ke) and asserting there tests arithmetic that was never true.
                    a = stack[f"{side}_ake"][-1]
                    b = stack[f"{side}_production"][-1] * stack[f"{side}_mean_ke"][-1]
                    fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
                    if fin.any():
                        resid = float(np.max(np.abs(a[fin] - b[fin]) / a[fin]))
                        assert resid < 1e-9, f"pooled decomposition broken per member: {resid:.2e}"
                    for q in STORM:
                        stack[f"{side}_{q}"].append(bin_mean(d[f"{side}_{q}"], iy, ix, (ny, nx)))
                c = np.zeros((ny, nx))
                np.add.at(c, (iy, ix), 1.0)
                n_cells = c if n_cells is None else n_cells + c
                print(f"  {scenario} m{member:<2} events={d['lat'].size:,}", end="\r", flush=True)
            counts[scenario] = n_cells / len(args.members)
            panels[scenario] = stack
            print(f"  {scenario}: {len(args.members)} members, "
                  f"{int(np.isfinite(stack['h_mean_ke'][0]).sum())} occupied cells, "
                  f"{int((counts[scenario] >= args.min_events).sum())} at n>={args.min_events}"
                  + " " * 20)

        if not panels:
            print("nothing processed", file=sys.stderr)
            return 1

        quantities = (*POOLED, *STORM)
        ds = xr.Dataset(coords={"lat": ("lat", lats), "lon": ("lon", lons),
                                "scenario": ("scenario", list(panels))})
        shape3 = (len(panels), ny, nx)
        for q in quantities:
            for side, label in (("h", "historical"), ("f", "future")):
                med = np.full(shape3, np.nan)
                for i, sc in enumerate(panels):
                    med[i] = np.nanmedian(np.stack(panels[sc][f"{side}_{q}"]), axis=0)
                ds[f"{q}_{label}"] = (("scenario", "lat", "lon"), med)
            chg = np.full(shape3, np.nan)
            lo_ = np.full(shape3, np.nan)
            hi_ = np.full(shape3, np.nan)
            for i, sc in enumerate(panels):
                per_member = np.stack([
                    100.0 * (f / np.where(h > 0, h, np.nan) - 1.0)
                    for h, f in zip(panels[sc][f"h_{q}"], panels[sc][f"f_{q}"])
                ])
                chg[i] = np.nanmedian(per_member, axis=0)
                lo_[i] = np.nanmin(per_member, axis=0)
                hi_[i] = np.nanmax(per_member, axis=0)
            ds[f"{q}_change_pct"] = (("scenario", "lat", "lon"), chg)
            ds[f"{q}_change_pct_member_min"] = (("scenario", "lat", "lon"), lo_)
            ds[f"{q}_change_pct_member_max"] = (("scenario", "lat", "lon"), hi_)

        ds["n_events"] = (("scenario", "lat", "lon"), np.stack([counts[sc] for sc in panels]))
        thin = ds["n_events"] < args.min_events
        for v in list(ds.data_vars):
            if v != "n_events":
                ds[v] = ds[v].where(~thin)

        ds["mean_ke_historical"].attrs = {"units": "J", "long_name":
            "mean kinetic energy per landed stone, pooled over stones in the cell (historical)"}
        ds["mean_diameter_historical"].attrs = {"units": "mm", "long_name":
            "mean diameter per landed stone, pooled over stones in the cell (historical)"}
        ds["ake_historical"].attrs = {"units": "J", "long_name":
            "mean kinetic energy per seeded embryo (authors' severity_calc_n; melted embryos count as zero)"}
        ds["production_historical"].attrs = {"units": "1", "long_name":
            "share of seeded embryos reaching the ground as hail"}
        ds["storm_mean_diameter_historical"].attrs = {"units": "mm", "long_name":
            "mean over storms of each storm's own mean landed diameter (storm-weighted)"}
        ds["storm_mean_ke_historical"].attrs = {"units": "J", "long_name":
            "mean over storms of each storm's own mean landed kinetic energy (storm-weighted)"}

        nc = args.out_dir / f"hail_severity_{tag}deg.nc"
        enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in ds.data_vars}
        ds.to_netcdf(nc, encoding=enc)
        print(f"\nwrote {nc}")

        # Reported, not asserted. The identity is exact per member (asserted above); after
        # taking the across-member MEDIAN it no longer holds, because the median of a product
        # is not the product of the medians. The residual measures that non-linearity only.
        for label in ("historical", "future"):
            lhs = ds[f"ake_{label}"].values
            rhs = (ds[f"production_{label}"] * ds[f"mean_ke_{label}"]).values
            ok = np.isfinite(lhs) & np.isfinite(rhs) & (lhs > 0)
            if ok.any():
                rel = np.abs(lhs[ok] - rhs[ok]) / lhs[ok]
                print(f"  median non-linearity ({label}): max {np.max(rel):.2e}, "
                      f"median {np.median(rel):.2e} over {ok.sum()} cells")

        rows = []
        for sc in panels:
            for name, (la0, la1, lo0, lo1) in REGIONS.items():
                sel = ds.sel(scenario=sc).where(
                    (ds.lat >= la0) & (ds.lat <= la1) & (ds.lon >= lo0) & (ds.lon <= lo1))
                n = float(np.nansum(sel["n_events"].values))
                if n < args.min_events:
                    continue
                row = {"scenario": sc, "region": name, "n_events": round(n)}
                for q in quantities:
                    row[f"{q}_hist"] = float(np.nanmean(sel[f"{q}_historical"].values))
                    row[f"{q}_fut"] = float(np.nanmean(sel[f"{q}_future"].values))
                    row[f"{q}_change_pct"] = float(np.nanmean(sel[f"{q}_change_pct"].values))
                rows.append(row)
        csv_path = args.out_dir / f"hail_severity_regions_{tag}deg.csv"
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {csv_path}")

        print("\nREGIONAL SUMMARY -- change vs historical, median over members")
        print(f"{'region':<20}{'scen':<9}{'n':>6}  {'mean KE (J)':>22}  {'mean size (mm)':>22}  {'production':>20}")
        for r in rows:
            print(f"{r['region']:<20}{r['scenario']:<9}{r['n_events']:>6}  "
                  f"{r['mean_ke_hist']:7.3f}->{r['mean_ke_fut']:6.3f} {r['mean_ke_change_pct']:+6.1f}%  "
                  f"{r['mean_diameter_hist']:7.1f}->{r['mean_diameter_fut']:6.1f} {r['mean_diameter_change_pct']:+6.1f}%  "
                  f"{r['production_hist']:5.3f}->{r['production_fut']:5.3f} {r['production_change_pct']:+6.1f}%")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
