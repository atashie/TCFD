#!/usr/bin/env python3
"""Reproduce the ISIMIP3b `floodedarea` characterisation the flood decision rests on.

    python scripts/measure_floodedarea_characteristics.py            # all sections
    python scripts/measure_floodedarea_characteristics.py --section mask

WHY THIS EXISTS
---------------
For eighteen days the repository carried one sentence about `floodedarea`, repeated in six
files: "binary {0,1} despite units='1', non-NaN over 94.7% of the globe INCLUDING OCEAN --
a per-variable mask defect". It was a true spot observation made while building `driedarea`,
and it was the ONLY thing anyone knew about the dataset. It had no retained receipt: no
script, no member list, no record of whether it came from one file or forty-five. It was
then quoted in `config/hazard_taxonomy.yaml` as the blocker on the hazard the taxonomy calls
the most consequential gap in the pipeline -- so a decision not to ingest river flood rested
on a number nobody could reproduce.

Re-measured 2026-08-13 across ALL 45 members, the observation holds and its FRAMING does not.
Ocean is published as exact zero, not as garbage: of 165,774 finite ocean cells only 258
(0.16%) are ever flagged, and 78.3% of those touch a land cell -- the model grid and the
generic ISIMIP mask disagreeing at fjords and small islands. That is a masking CONVENTION
with a one-line fix, not a corrupt field. The mask is still mandatory: left in place, 64% of
the grid enters the percentile baseline as zeros.

What the same pass established that the spot check could not:

  * inter-model spread is 1.18-1.25x, far tighter than `driedarea` (2.69x) or `led` (7.8x),
    while the GCM spread inside one model is 1.76x -- uncertainty here is forcing-dominated
  * the publication-mask question matches NEITHER shipped precedent (see below)
  * the Amazon floodplain reads ~0.000 in every member, which is the family's
    departure-from-preindustrial definition showing through, not a hole in the data

That last one is why this script prints reference sites: a layer that reads zero across the
world's largest floodplain will be asked about by the first person who looks at the map, and
"measured, and here is why" is the only good answer.

Requires the 45-member raw ensemble in data/raw/flood-isimip3b_floodedarea_annual/.
Nothing here is layer-registry aware -- this runs BEFORE any decision to process.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.land_mask import get_isimip_landmask, load_landmask  # noqa: E402

RAW = Path("data/raw/flood-isimip3b_floodedarea_annual")
VAR = "floodedarea"
MODELS = ["h08", "jules-w2", "watergap2-2e"]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

#: GUARDRAILS 12 -- places river flooding demonstrably occurs, plus two dry controls that
#: SHOULD come back unpublished. A layer can satisfy every line of OUTPUT-SPEC and still be
#: about nothing; this is the check that asks whether it is about rivers.
SITES = {
    "Ganges-Brahmaputra delta": (23.75, 90.25),
    "Mekong delta": (10.25, 105.75),
    "Lower Mississippi": (32.25, -91.25),
    "Rhine (NL/DE)": (51.75, 6.25),
    "Amazon floodplain": (-3.25, -60.25),
    "Niger inland delta": (14.75, -4.25),
    "Sahara (dry control)": (24.75, 10.25),
    "Atacama (dry control)": (-23.25, -69.25),
}


def member_path(model: str, gcm: str, scenario: str) -> Path:
    return RAW / f"{model}_{gcm}_w5e5_{scenario}_2015soc_default_{VAR}_global_annual_2015_2100.nc"


def load(model: str, gcm: str, scenario: str):
    """Return (values[time, lat, lon], years, attrs). Time is NOT decoded -- see below."""
    path = member_path(model, gcm, scenario)
    if not path.exists():
        raise SystemExit(f"missing member: {path}\nDownload the ensemble first.")
    # decode_times=False deliberately: the axis is `days since 1601-01-01` on
    # proleptic_gregorian, and this script only needs the year INDEX, not calendar dates.
    # A processor must decode it properly with cftime -- days/365 drifts ~0.4 yr over a
    # 400-year offset, enough to push a January record into the wrong decade bin.
    ds = xr.open_dataset(path, decode_times=False)
    var = ds[VAR]
    values = var.values
    attrs = {
        "units": var.attrs.get("units", "?"),
        "long_name": var.attrs.get("long_name", "?"),
        "time_units": ds["time"].attrs.get("units", "?"),
        "calendar": ds["time"].attrs.get("calendar", "?"),
    }
    ds.close()
    return values, 2015 + np.arange(values.shape[0]), attrs


def landmask():
    mask = load_landmask(get_isimip_landmask("3b"))
    land = np.nan_to_num(np.asarray(mask.values), nan=0.0) > 0.5
    return land, np.asarray(mask["lat"].values)


def section_nature(land):
    """Data nature, metadata, and the mask convention -- the GUARDRAILS 9 check."""
    print("\n=== nature, metadata and mask, all 45 members ===")
    seen_attrs, rows = {}, []
    for model in MODELS:
        for gcm in GCMS:
            for scenario in SCENARIOS:
                a, _, attrs = load(model, gcm, scenario)
                key = tuple(attrs.values())
                seen_attrs[key] = seen_attrs.get(key, 0) + 1
                finite = np.isfinite(a)
                vals = a[finite]
                uniq = np.unique(vals[:2_000_000])
                cell_any = finite.any(axis=0)
                ocean = cell_any & ~land
                land_c = cell_any & land
                oc = a[:, ocean][np.isfinite(a[:, ocean])]
                ln = a[:, land_c][np.isfinite(a[:, land_c])]
                rows.append((model, gcm, scenario, cell_any.mean(), uniq.size,
                             bool(np.all(np.isin(uniq, (0.0, 1.0)))), float(vals.max()),
                             int(ocean.sum()), float((oc > 0).mean()),
                             int(land_c.sum()), float(ln.mean()), float((ln == 0).mean()),
                             bool((cell_any == finite.all(axis=0)).all())))
    for key, n in seen_attrs.items():
        print(f"  units={key[0]!r} long_name={key[1]!r} time={key[2]!r} cal={key[3]!r}  x{n}")
    print(f"\n{'model':<14}{'gcm':<15}{'scen':<8}{'finite%':>8}{'uniq':>5}{'bin':>6}"
          f"{'max':>6}{'ocn cells':>11}{'ocn>0%':>8}{'land cells':>11}{'land mean':>10}{'land 0%':>9}")
    for r in rows:
        print(f"{r[0]:<14}{r[1]:<15}{r[2]:<8}{100*r[3]:>8.1f}{r[4]:>5}{str(r[5]):>6}"
              f"{r[6]:>6.2f}{r[7]:>11}{100*r[8]:>8.2f}{r[9]:>11}{r[10]:>10.4f}{100*r[11]:>9.2f}")
    print(f"\n  finite mask time-invariant in every member: {all(r[12] for r in rows)}")
    print("  READ: binary {0,1} everywhere -> a decadal mean is a FREQUENCY of flagged years,")
    print("        not the area share that units='1' and long_name='Exposed Area Share' imply.")


def section_mask(land, lat):
    """What the non-NaN ocean cells actually contain."""
    from scipy import ndimage

    print("\n=== the 'ocean' cells: convention or corruption? (jules-w2/gfdl-esm4/ssp585) ===")
    a, _, _ = load("jules-w2", "gfdl-esm4", "ssp585")
    finite = np.isfinite(a).any(axis=0)
    ever = (np.nan_to_num(a, nan=0.0) > 0).any(axis=0)
    hit = ever & ~land & finite
    near_land = ndimage.binary_dilation(land, iterations=1)
    print(f"  finite ocean cells                    : {(finite & ~land).sum()}")
    print(f"  ...ever flagged flooded               : {hit.sum()} ({100*hit.sum()/max((finite&~land).sum(),1):.2f}%)")
    print(f"  ...of those, adjacent to a land cell  : {100*(hit & near_land).sum()/max(hit.sum(),1):.1f}%")
    ii, jj = np.where(hit)
    freq = np.nan_to_num(a, nan=0.0)[:, ii, jj].mean(axis=0)
    lon = np.asarray(xr.open_dataset(member_path("jules-w2", "gfdl-esm4", "ssp585"),
                                     decode_times=False)["lon"].values)
    for k in np.argsort(-freq)[:6]:
        print(f"     lat {lat[ii[k]]:>7.2f} lon {lon[jj[k]]:>8.2f}  mean flag {freq[k]:.3f}")
    print(f"\n  ISIMIP3b land cells                   : {land.sum()}")
    print(f"  published land cells                  : {(finite & land).sum()} "
          f"({100*(finite&land).sum()/land.sum():.1f}% of land)")
    unpub = land & ~finite
    print(f"  land cells published by this member   : missing {unpub.sum()}, "
          f"median |lat| {np.median(np.abs(lat[np.where(unpub)[0]])):.1f} (the arid belt: no river network)")


def section_signal(land):
    """Decadal signal, inter-model spread, and the GCM-vs-model uncertainty split."""
    print("\n=== global land-mean exposure frequency by decade ===")
    print(f"{'scen':<9}{'model':<14}{'2020s':>9}{'2050s':>9}{'2090s':>9}{'90s/20s':>10}")
    store = {}
    for scenario in SCENARIOS:
        for model in MODELS:
            dec = {}
            for d0 in (2020, 2050, 2090):
                per_gcm = []
                for gcm in GCMS:
                    a, years, _ = load(model, gcm, scenario)
                    sel = (years >= d0) & (years < d0 + 10)
                    per_gcm.append(np.nanmean(a[sel][:, land]))
                dec[d0] = float(np.mean(per_gcm))
            store[(scenario, model)] = dec
            print(f"{scenario:<9}{model:<14}{dec[2020]:>9.4f}{dec[2050]:>9.4f}"
                  f"{dec[2090]:>9.4f}{dec[2090]/max(dec[2020],1e-9):>9.2f}x")
    print("\n=== inter-model spread at the 2090s ===")
    for scenario in SCENARIOS:
        v = [store[(scenario, m)][2090] for m in MODELS]
        print(f"  {scenario}: {' / '.join(f'{m}={x:.4f}' for m, x in zip(MODELS, v))}"
              f"   spread={max(v)/min(v):.2f}x")
    print("\n=== GCM spread WITHIN one model (h08, ssp585, 2090s) -- the larger term ===")
    for gcm in GCMS:
        a, _, _ = load("h08", gcm, "ssp585")
        print(f"  {gcm:<15}{np.nanmean(a[-11:-1][:, land]):.4f}")


def section_publication_mask(land, lat):
    """The >=N-model publication question. Measured per layer, NEVER inherited."""
    print("\n=== publication mask evidence (ssp585, gfdl-esm4) ===")
    finite = {m: np.isfinite(load(m, "gfdl-esm4", "ssp585")[0]).any(axis=0) for m in MODELS}
    count = sum(f.astype(int) for f in finite.values())
    for k in range(4):
        n = int(((count == k) & land).sum())
        print(f"  {k} model(s): {n:>7} land cells ({100*n/land.sum():>5.1f}% of land)")
    means = {}
    for m in MODELS:
        a, _, _ = load(m, "gfdl-esm4", "ssp585")
        means[m] = np.nanmean(np.where(np.isfinite(a), a, np.nan), axis=0)
    ens = np.nanmean(np.stack([means[m] for m in MODELS]), axis=0)
    for k in (1, 2, 3):
        sel = (count == k) & land & np.isfinite(ens)
        print(f"  ensemble-mean exposure on {k}-model cells: {np.nanmean(ens[sel]):.4f} (n={sel.sum()})")
    print("\n  who publishes the solo cells:")
    for m in MODELS:
        solo = (count == 1) & finite[m] & land
        print(f"    {m:<14}{solo.sum():>6} cells, median |lat| "
              f"{np.median(np.abs(lat[np.where(solo)[0]])):.1f}")
    print("  READ: solo cells run ~2.2x the 3-model level while inter-model spread is only")
    print("        ~1.2x, and they are split across all three models -- so this is NOT led's")
    print("        'one hot model undiluted' mechanism, and it is NOT driedarea's null result.")
    print("        The rule has to be decided on THIS evidence.")


def section_sites():
    """GUARDRAILS 12 -- is this layer about rivers?"""
    print("\n=== reference sites: mean annual flag, first 10 yr -> last 10 yr ===")
    print(f"{'site':<26}{'model':<14}{'ssp126':>16}{'ssp370':>16}{'ssp585':>16}")
    for site, (la, lo) in SITES.items():
        for i, model in enumerate(MODELS):
            cells = []
            for scenario in SCENARIOS:
                acc_first, acc_last = [], []
                for gcm in GCMS:
                    path = member_path(model, gcm, scenario)
                    ds = xr.open_dataset(path, decode_times=False)
                    s = ds[VAR].sel(lat=la, lon=lo, method="nearest").values
                    ds.close()
                    acc_first.append(np.nanmean(s[:10]))
                    acc_last.append(np.nanmean(s[-10:]))
                with np.errstate(invalid="ignore"):
                    cells.append(f"{np.nanmean(acc_first):>7.3f}->{np.nanmean(acc_last):<8.3f}")
            print(f"{site if i == 0 else '':<26}{model:<14}{''.join(cells)}")
    print("\n  READ: the Amazon floodplain reading ~0 is the family's departure-from-")
    print("        preindustrial definition, not a hole. This is NOT a flood-hazard map.")


SECTIONS = {
    "nature": lambda land, lat: section_nature(land),
    "mask": section_mask,
    "signal": lambda land, lat: section_signal(land),
    "publication-mask": section_publication_mask,
    "sites": lambda land, lat: section_sites(),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--section", choices=sorted(SECTIONS), action="append",
                    help="run only these sections (default: all)")
    args = ap.parse_args()
    if not RAW.exists():
        print(f"missing raw ensemble: {RAW}", file=sys.stderr)
        return 1
    land, lat = landmask()
    print(f"ISIMIP3b land-sea mask: {land.sum()} land cells of {land.size}")
    for name in (args.section or sorted(SECTIONS)):
        SECTIONS[name](land, lat)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
