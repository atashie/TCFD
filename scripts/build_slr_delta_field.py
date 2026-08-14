"""Stage 1 of the coastal layer: the sea-level DELTA field, QC'd and filled to land.

Turns the eight ISIMIP2b `sealevelrise` files into an annual, per-member field of
sea-level change relative to the shared 2020s baseline, defined on LAND cells -- which is
where the DEM needs it and where the published field is not.

WHY A FILL STAGE EXISTS AT ALL
------------------------------
Measured 2026-08-14 (scripts/measure_sealevelrise_support.py): the published field is an
OCEAN field. It is finite on 83.4% of ocean cells and on only 2.9% of land cells, and of
9,493 coastal land cells 81% carry no value. The elevation calculation needs a sea level
at the LAND cell, so those have to borrow one from the nearest cell where the ocean model
actually solved. That borrow distance is a real quantity a customer can be misled by, so
it is computed exactly per member, stored, and carried into the layer as a QA field.

THE MIROC5 DEFECT -- WHY THERE IS A CONSENSUS MASK
--------------------------------------------------
Measured 2026-08-14 (scripts/check_slr_member_outliers.py). MIROC5's 2020s->2090s delta
in the Mediterranean is -1.635 m mean / -4.076 m min while GFDL-ESM2M, HadGEM2-ES and
IPSL-CM5A-LR give +0.154, +0.368 and +0.183. Against a leave-one-out consensus, MIROC5
deviates by more than 0.5 m on 10.31% of cells at rcp26 and 12.32% at rcp60; no other
member exceeds 1.28% anywhere.

Pooled with all four members the ensemble publishes a FALLING sea level:

    Mediterranean   4-member -0.753 m   vs 3-member +0.238 m   SIGN FLIP
    Black Sea       4-member -1.053 m   vs 3-member +0.217 m   SIGN FLIP

which would tell Venice, Alexandria, the Po and Nile deltas, Istanbul and Odesa that
their sea level is falling. MIROC5's GLOBAL mean delta is +0.311 m, squarely mid-pack --
so no global summary catches this. It is a coarse ocean model failing to exchange water
through straits it cannot resolve, and it is sound in the open ocean (deviation +0.06 to
+0.09 m).

The member is therefore kept where it agrees and dropped where it does not (user decision
2026-08-14), rather than discarded wholesale. `n_members` becomes 3 or 4 per cell and the
rule is declared in the file.

TWO CHOSEN CONSTANTS, BOTH RECORDED IN THE OUTPUT
-------------------------------------------------
CONSENSUS_TOL_M -- how far a member may sit from the median of ALL members before it is
dropped at that cell. 0.5 m separates the populations cleanly: no sound member reaches
even 1.3% of cells at this threshold while MIROC5 reaches 10-12%.

MAX_BORROW_KM -- how far a land cell may reach for a sea level. Without it the mask
backfires: MIROC5 masked across the Mediterranean would simply borrow from the Atlantic,
reinstating a wrong value by another route. 250 km passes the legitimate coastal borrows
(measured p90 = 144 km on a single member's own support) and rejects cross-basin
substitution.

THE MASK IS COMPUTED ACROSS BOTH SCENARIOS AT ONCE
--------------------------------------------------
A member masked at a cell in rcp60 but not rcp26 would give the two scenarios different
ensemble composition, and OUTPUT-SPEC's shared 2020s baseline is only bit-identical when
composition is uniform. A member failing in EITHER scenario is therefore dropped at that
cell in BOTH, and the deviation is judged on every decade rather than on the endpoint
alone -- a member can agree in the 2090s having been wild in the 2050s.

PERCENTILE LEVEL
----------------
The file carries a 7-level `percentile` axis [5, 16.6, 33, 50, 66, 83.3, 95] whose
long_name is the unhelpful string "generic". p50 is taken, and that choice is close to
free for a DELTA design: across the whole axis the 2020s->2090s global-mean delta spans
+0.191 to +0.223 m at rcp26, a 3 cm spread on a 21 cm signal. The absolute levels differ
far more (2090s mean +0.357 at p5 vs +0.448 at p95), so this reasoning would NOT hold for
an absolute-height design.

Run:  .venv/bin/python3 scripts/build_slr_delta_field.py
Out:  data/interim/slr_delta/slr_delta_{scenario}.nc
"""

from pathlib import Path

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.spatial import cKDTree

RAW = Path("data/raw/sealevelrise_2b")
OUT = Path("data/interim/slr_delta")
MASK = RAW / "landseamask_generic.nc4"

GCMS = ["GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5"]
SCENARIOS = ["rcp26", "rcp60"]
BASELINE_DECADE = 2020
DECADES = list(range(2010, 2100, 10))
PCT_INDEX = 3               # p50 of [5, 16.6, 33, 50, 66, 83.3, 95]
CONSENSUS_TOL_M = 0.5
MAX_BORROW_KM = 250.0
BODY_FAIL_FRAC = 0.25
EARTH_R_KM = 6371.0


def _path(gcm: str, scen: str) -> Path:
    return RAW / f"total_year_{gcm}_{scen}_2006_2099.nc4"


def load_stack():
    """(scenario -> (member, year, lat, lon)) at p50, plus the shared coords."""
    stacks, years_ref, lat, lon = {}, None, None, None
    for scen in SCENARIOS:
        members = []
        for gcm in GCMS:
            ds = xr.open_dataset(_path(gcm, scen), decode_times=False)
            years = 1661 + ds["time"].values.astype(int)
            if years_ref is None:
                years_ref, lat, lon = years, ds["lat"].values, ds["lon"].values
            elif not np.array_equal(years, years_ref):
                raise ValueError(f"{gcm}/{scen} year axis differs from the first member")
            members.append(ds["total"].isel(percentile=PCT_INDEX).values.astype("float32"))
            ds.close()
        stacks[scen] = np.stack(members)
    return stacks, years_ref, lat, lon


def water_bodies(ocean: np.ndarray) -> np.ndarray:
    """Label connected water bodies, joined across the antimeridian.

    The basins where a coarse ocean model fails are, conveniently, the ones a 0.5 deg grid
    cannot connect to the open ocean: Gibraltar is ~14 km wide and the Bosphorus ~1 km, so
    the Mediterranean and Black Sea come out as their own components. That is what makes a
    basin-level verdict possible at all.
    """
    lab, n = ndimage.label(ocean, structure=np.ones((3, 3), int))
    # Union labels that touch across the lon seam -- ndimage does not wrap.
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for r in range(ocean.shape[0]):
        for dr in (-1, 0, 1):
            rr = r + dr
            if not (0 <= rr < ocean.shape[0]):
                continue
            if lab[r, 0] and lab[rr, -1]:
                a, b = find(lab[r, 0]), find(lab[rr, -1])
                if a != b:
                    parent[b] = a
    remap = np.array([find(i) for i in range(n + 1)])
    remap[0] = 0
    return remap[lab]


def consensus_mask(deltas_by_scen, years, ocean):
    """True where a member may be used, judged per cell and then per water body.

    JUDGED AGAINST THE MEDIAN OF ALL MEMBERS, NOT A LEAVE-ONE-OUT MEDIAN. With four
    members, leave-one-out leaves three, and where coverage is ragged it can leave two --
    at which point the median sits exactly between a sound member and a broken one and
    flags BOTH. Measured in the Black Sea, where HadGEM2-ES publishes nothing: against the
    other two, GFDL-ESM2M deviates 1.14 m from a median dragged down by MIROC5 and would
    have been thrown out for being right. The median of ALL finite members is resistant to
    one outlier as soon as three are present, and picks out only the member that is wrong.

    ESCALATED TO THE WHOLE WATER BODY. A per-cell verdict alone does not work here, and
    the failure is not subtle: MIROC5 resolves 862 Mediterranean cells where the other
    members resolve ~660, so ~180 near-shore cells have too few members to judge. Left as
    "keep by default" they survive -- and because they are the cells closest to land, they
    win the nearest-neighbour fill and put MIROC5's -2.35 m straight back onto the
    Mediterranean coast. Trusting a member exactly where nothing can check it is the wrong
    default, so a member failing on more than BODY_FAIL_FRAC of the cells it WAS judged on
    within a water body is dropped from that body entirely, unjudged cells included.
    """
    shape = deltas_by_scen[SCENARIOS[0]].shape[2:]
    keep = np.ones((len(GCMS),) + shape, dtype=bool)
    judged = np.zeros_like(keep)
    failed = np.zeros((len(GCMS),) + shape, dtype=bool)

    for scen in SCENARIOS:
        d = deltas_by_scen[scen]
        for dec in DECADES:
            sel = (years >= dec) & (years < dec + 10)
            if not sel.any():
                continue
            with np.errstate(invalid="ignore"):
                panel = np.nanmean(d[:, sel], axis=1)            # (member, lat, lon)
                med = np.nanmedian(panel, axis=0)                # median of ALL members
            n_fin = np.isfinite(panel).sum(axis=0)
            for k in range(len(GCMS)):
                ok = (n_fin >= 3) & np.isfinite(panel[k]) & np.isfinite(med)
                judged[k] |= ok
                bad = ok & (np.abs(panel[k] - med) > CONSENSUS_TOL_M)
                keep[k] &= ~bad
                failed[k] |= bad

    bodies = water_bodies(ocean)
    for k in range(len(GCMS)):
        for lbl in np.unique(bodies[bodies > 0]):
            sel = bodies == lbl
            nj = int((judged[k] & sel).sum())
            if nj == 0:
                continue
            frac = (failed[k] & sel).sum() / nj
            if frac > BODY_FAIL_FRAC:
                keep[k] &= ~sel
    return keep, judged, bodies


def build_lookup(valid_2d, land, lat, lon):
    """Nearest valid cell for every land cell -- exact great-circle, via a KD-tree.

    A grid-index distance transform is wrong at both ends of this problem: it does not
    wrap at the antimeridian, and at 60N a longitude step is half a latitude step, so it
    biases the pick east-west. Unit vectors handle both for free.
    """
    LON, LAT = np.meshgrid(np.radians(lon), np.radians(lat))
    xyz = np.stack([np.cos(LAT) * np.cos(LON),
                    np.cos(LAT) * np.sin(LON),
                    np.sin(LAT)], axis=-1).reshape(-1, 3)
    src = np.flatnonzero(valid_2d.ravel())
    if src.size == 0:
        raise ValueError("no valid source cells")
    tgt = np.flatnonzero(land.ravel())
    chord, k = cKDTree(xyz[src]).query(xyz[tgt], k=1)
    km = 2 * np.arcsin(np.clip(chord / 2, 0, 1)) * EARTH_R_KM
    return tgt, src[k], km


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stacks, years, lat, lon = load_stack()
    print(f"loaded {len(SCENARIOS)} scenarios x {len(GCMS)} members, "
          f"years {years.min()}-{years.max()}")

    # The land mask is stored NORTH-DOWN while the sea-level files are SOUTH-UP. Align on
    # the coordinate, never on array order -- this put Antarctica on the Arctic once.
    lsm = xr.open_dataset(MASK)["LSM"].squeeze().sortby("lat")
    assert np.allclose(lsm["lat"].values, lat) and np.allclose(lsm["lon"].values, lon), \
        "land mask does not align with the sea-level grid"
    land = np.isfinite(np.asarray(lsm.values))
    ocean = ~land

    # --- shared 2020s reference, averaged across scenarios FIRST (OUTPUT-SPEC) ---------
    base = (years >= BASELINE_DECADE) & (years < BASELINE_DECADE + 10)
    ref = np.mean([stacks[s][:, base].mean(axis=1) for s in SCENARIOS], axis=0)
    deltas = {s: stacks[s] - ref[:, None] for s in SCENARIOS}

    keep, judged, bodies = consensus_mask(deltas, years, ocean)
    print(f"\nconsensus mask (|dev from the median of all members| <= {CONSENSUS_TOL_M} m, "
          f"every decade, both scenarios; >=3 members needed to judge):")
    for k, g in enumerate(GCMS):
        nj = int(judged[k].sum())
        print(f"  {g:14s} judged {nj:7,} cells   dropped {int((judged[k] & ~keep[k]).sum()):7,} "
              f"({100 * (judged[k] & ~keep[k]).sum() / max(nj, 1):5.2f}%)")

    big = np.bincount(bodies.ravel())[1:].argmax() + 1 if bodies.max() else 0
    print(f"  water bodies: {int(bodies.max())}, largest (open ocean) = {int((bodies==big).sum()):,} cells")

    coastal = land & (np.stack([np.roll(np.roll(ocean, dj, 0), di, 1)
                                for dj in (-1, 0, 1) for di in (-1, 0, 1)
                                if (dj or di)]).sum(0) > 0)
    print(f"\ncoastal land cells: {coastal.sum():,}")

    nlat, nlon = land.shape
    land_flat = np.flatnonzero(land.ravel())
    coastal_flat = coastal.ravel()[land_flat]

    for scen in SCENARIOS:
        d = deltas[scen].reshape(len(GCMS), len(years), -1)
        filled = np.full((len(GCMS), len(years), land_flat.size), np.nan, "float32")
        borrow = np.full((len(GCMS), land_flat.size), np.nan, "float32")

        for k, g in enumerate(GCMS):
            # A member's usable sources: finite for every year, on ocean, consensus-passing.
            # A member is a source ONLY where it was actually judged and passed.
            # `judged & keep` rather than `keep` alone: `keep` defaults to True at cells
            # too thinly covered to check, and those are exactly the shallow near-shore
            # cells only the finer model resolves -- the ones that then win the
            # nearest-neighbour fill. Escalating to the whole water body does not catch
            # them either: at 0.5 deg with diagonal connectivity the Mediterranean and
            # Black Sea link to the open ocean through Gibraltar and the Bosphorus, so
            # 98.9% of ocean cells label as ONE body and no basin-level verdict can fire.
            valid = np.isfinite(deltas[scen][k]).all(axis=0) & ocean & keep[k] & judged[k]
            tgt, src, km = build_lookup(valid, land, lat, lon)
            within = km <= MAX_BORROW_KM
            filled[k][:, within] = d[k][:, src[within]]
            borrow[k][within] = km[within]
            cb = km[coastal_flat & within]
            print(f"  [{scen}] {g:14s} sources {valid.sum():7,}  "
                  f"coastal cells served {int((coastal_flat & within).sum()):5,}/"
                  f"{int(coastal_flat.sum()):5,}  borrow p50 {np.median(cb):5.1f} "
                  f"p90 {np.percentile(cb, 90):6.1f} km")

        n_members = np.isfinite(filled).all(axis=1).sum(axis=0).astype("int16")
        ds = xr.Dataset(
            {
                "slr_delta": (("member", "year", "land_cell"), filled,
                              {"units": "m",
                               "long_name": "sea-level change relative to the shared 2020s "
                                            "baseline",
                               "description": "Taken from the nearest ocean cell where this "
                                              "member solved AND passed the consensus check; "
                                              "NaN where none lies within max_borrow_km."}),
                "borrow_km": (("member", "land_cell"), borrow,
                              {"units": "km",
                               "long_name": "great-circle distance to the sea-level cell "
                                            "this member's value was taken from"}),
                "n_members": (("land_cell",), n_members,
                              {"long_name": "members contributing at this land cell"}),
                "is_coastal": (("land_cell",), coastal_flat.astype("int8"),
                               {"long_name": "land cell adjacent to an ocean cell"}),
                "land_flat_index": (("land_cell",), land_flat.astype("int32"),
                                    {"long_name": f"flat index into the {nlat}x{nlon} grid"}),
            },
            coords={"member": GCMS, "year": years},
            attrs={
                "source": "ISIMIP2b InputData/sealevelrise, variable `total`, p50 of the "
                          "7-level percentile axis",
                "method": "Mengel et al. PNAS (2016) and Bamber and Riva (2010)",
                "scenario": scen,
                "baseline_decade": f"{BASELINE_DECADE}s",
                "baseline_rule": "each member's 2020s window averaged across rcp26 and rcp60 "
                                 "before subtraction, so the baseline panel is identical "
                                 "across scenarios (OUTPUT-SPEC)",
                "ocean_definition": "the sea-level field's own finite support, NOT the "
                                    "inverted land-sea mask -- the mask counts the Caspian "
                                    "and large lakes as water",
                "consensus_tol_m": CONSENSUS_TOL_M,
                "consensus_rule": "a member is used only where >=3 members are finite AND "
                                  "it sits within consensus_tol_m of the median of all "
                                  "members in EVERY decade of BOTH scenarios; cells too "
                                  "thinly covered to judge are not used as sources",
                "consensus_rationale": "MIROC5 gives -1.635 m mean Mediterranean delta vs "
                                       "+0.15/+0.37/+0.18 for the other members; the "
                                       "4-member mean flips sign in the Mediterranean "
                                       "(-0.753 -> +0.238) and Black Sea (-1.053 -> +0.217). "
                                       "See scripts/check_slr_member_outliers.py",
                "max_borrow_km": MAX_BORROW_KM,
                "max_borrow_rationale": "without a cutoff a member masked across a basin "
                                        "borrows from the next one, reinstating the wrong "
                                        "value by another route",
                "grid": f"{nlat}x{nlon} at 0.5 deg, south-up",
            },
        )
        p = OUT / f"slr_delta_{scen}.nc"
        ds.to_netcdf(p)
        c = coastal_flat
        print(f"  [{scen}] wrote {p} ({p.stat().st_size / 1e6:.0f} MB); coastal n_members "
              f"mean {n_members[c].mean():.2f}, "
              f"=4 on {100 * (n_members[c] == 4).mean():.1f}% of coastal cells\n")
        ds.close()


if __name__ == "__main__":
    main()
