"""Is any sea-level member defective, and where?

MIROC5's Mediterranean 2020s->2090s delta is -1.635 m mean / -4.076 m min while the other
three members give +0.154, +0.368 and +0.183. A sea-level FALL of that size in a basin
ringed with coastal assets is not a fingerprint, it is a model failing to exchange water
through a strait it cannot resolve. Its GLOBAL mean (+0.311 m) is unremarkable, which is
exactly why this would have shipped: every aggregate summary looks fine.

With four members, one member reading -1.6 m drags the Mediterranean ensemble mean to
about -0.23 m -- the layer would tell a customer in Venice, Alexandria or the Po delta
that their sea level is FALLING.

This script does not assume the problem is MIROC5 or the Mediterranean. For every member
it compares that member's delta against the median of the OTHER members (a leave-one-out
consensus), so each member is judged by the same rule, and reports where the deviations
concentrate.

Run:  .venv/bin/python3 scripts/check_slr_member_outliers.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

RAW = Path("data/raw/sealevelrise_2b")
GCMS = ["GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5"]
SCENARIOS = ["rcp26", "rcp60"]
PCT_INDEX = 3

LAT = np.arange(-89.75, 90, 0.5)
LON = np.arange(-179.75, 180, 0.5)

# Semi-enclosed basins -- the places a coarse ocean model is most likely to get wrong,
# named so the finding is checkable rather than a blob of coordinates.
BASINS = {
    "Mediterranean":  (30.0, 46.0, -6.0, 36.0),
    "Black Sea":      (40.5, 47.5, 27.0, 42.0),
    "Baltic":         (53.5, 66.0, 10.0, 30.0),
    "Red Sea":        (12.0, 30.0, 32.0, 44.0),
    "Persian Gulf":   (24.0, 30.5, 48.0, 57.0),
    "Hudson Bay":     (51.0, 65.0, -95.0, -75.0),
    "Kara/Barents":   (68.0, 78.0, 40.0, 80.0),
    "Sea of Japan":   (34.0, 46.0, 128.0, 142.0),
    "Open Atlantic":  (20.0, 50.0, -60.0, -20.0),
    "Open Pacific":   (-20.0, 20.0, -160.0, -120.0),
}


def delta_2090s(gcm: str, scen: str) -> np.ndarray:
    ds = xr.open_dataset(RAW / f"total_year_{gcm}_{scen}_2006_2099.nc4", decode_times=False)
    yrs = 1661 + ds["time"].values.astype(int)
    v = ds["total"].isel(percentile=PCT_INDEX).values
    out = v[yrs >= 2090].mean(0) - v[(yrs >= 2020) & (yrs < 2030)].mean(0)
    ds.close()
    return out


def main():
    for scen in SCENARIOS:
        print("=" * 78)
        print(f"{scen}: leave-one-out consensus deviation, 2020s -> 2090s delta")
        print("=" * 78)
        stack = np.stack([delta_2090s(g, scen) for g in GCMS])       # (member, lat, lon)
        valid = np.isfinite(stack).all(0)
        print(f"cells finite in all members: {valid.sum():,}\n")

        print(f"{'member':16s} {'median dev':>11s} {'|dev|>0.5m':>11s} {'|dev|>1m':>9s} "
              f"{'|dev|>2m':>9s}")
        devs = {}
        for k, g in enumerate(GCMS):
            others = np.median(np.delete(stack, k, axis=0), axis=0)
            dev = np.where(valid, stack[k] - others, np.nan)
            devs[g] = dev
            n = np.isfinite(dev).sum()
            print(f"{g:16s} {np.nanmedian(dev):+11.3f} "
                  f"{100*np.nansum(np.abs(dev) > 0.5)/n:10.2f}% "
                  f"{100*np.nansum(np.abs(dev) > 1.0)/n:8.2f}% "
                  f"{100*np.nansum(np.abs(dev) > 2.0)/n:8.2f}%")

        print(f"\n{'basin':16s} " + " ".join(f"{g[:12]:>13s}" for g in GCMS)
              + "   <- mean delta (m), and [dev from the other three]")
        for name, (la0, la1, lo0, lo1) in BASINS.items():
            sj = (LAT >= la0) & (LAT < la1)
            si = (LON >= lo0) & (LON < lo1)
            cells = []
            for k, g in enumerate(GCMS):
                sub = stack[k][np.ix_(sj, si)]
                dsub = devs[g][np.ix_(sj, si)]
                cells.append(f"{np.nanmean(sub):+6.3f}[{np.nanmean(dsub):+5.2f}]")
            print(f"{name:16s} " + " ".join(f"{c:>13s}" for c in cells))

        # What the ensemble would publish, with and without the worst member.
        print("\n--- what the 4-member ensemble mean would publish vs a 3-member one ---")
        worst = max(GCMS, key=lambda g: np.nansum(np.abs(devs[g]) > 1.0))
        keep = [i for i, g in enumerate(GCMS) if g != worst]
        print(f"worst member by |dev|>1m count: {worst}")
        print(f"{'basin':16s} {'all 4':>10s} {'without':>10s}  {'shift':>8s}")
        for name, (la0, la1, lo0, lo1) in BASINS.items():
            sj = (LAT >= la0) & (LAT < la1)
            si = (LON >= lo0) & (LON < lo1)
            a4 = np.nanmean(stack[:, sj][:, :, si])
            a3 = np.nanmean(stack[keep][:, sj][:, :, si])
            flag = "  <-- SIGN FLIP" if a4 * a3 < 0 else ""
            print(f"{name:16s} {a4:+10.3f} {a3:+10.3f}  {a3-a4:+8.3f}{flag}")
        print()


if __name__ == "__main__":
    main()
