#!/usr/bin/env python3
"""Evidence gate for the Nature 2026 global hailstone deposit (Zhang et al.).

WHAT THIS IS. `Rising global hail damage potential in a warming world`, Nature 653,
1069-1077 (2026), doi 10.1038/s41586-026-10543-2, deposits its processed model output
under CC-BY-4.0: figshare 10.6084/m9.figshare.30103471 (data) and
10.6084/m9.figshare.30103474 (code). This script measures that deposit against the
published claims BEFORE any decision to build a layer from it. It writes nothing to
data/processed and satisfies no part of OUTPUT-SPEC.md -- see WHY IT IS NOT A LAYER.

WHAT THE DATA IS. A semi-three-dimensional hail trajectory model was run on 12,412
satellite-detected severe hailstorm events (Zhou et al. 2021 SOM classification, five
environment types) spanning 2014-04-01 to 2021-03-31. Historical runs are driven by ERA5
soundings; future runs re-drive THE SAME EVENTS with pseudo-global-warming deltas, so the
design is a paired perturbation experiment with no time axis. The ensemble is 20
EC-Earth3 realizations per scenario (n_models = 1); MPI-ESM1-2-LR and NorESM2-LM appear
only as ensemble-mean sensitivity runs. Scenarios: ssp245, ssp370, ssp585.

THE JOIN, WHICH IS NOT THE OBVIOUS ONE. `hailstone_growth_radii_*` files renumber events
0..n-1 and carry no coordinates, so a positional HIST-row to FUT-row join is WRONG
(measured r = -0.0002, i.e. noise). The real index lives in the sounding profile files as
`date{type}`, which are integer row indices into the catalogue's `para33_Idx{type}`; the
growth-duration files copy the same array into their `sample` coordinate. The authors pair
exactly this way -- see `hailstone_trajectory_composite_for_HMU-LA_and_HMU-HA...py`:

    condisf = np.isin(proff['date'+str(i+1)].values, profh['date'+str(i+1)].values)
    condish = np.isin(profh['date'+str(i+1)].values, proff['date'+str(i+1)].values)

Confirmed four ways: per-type counts match the radii arrays exactly; paired
Spearman(HIST, FUT) = 0.248 against a permutation null of 0.001; hail diameter correlates
with the CATALOGUE's MUCAPE for the mapped event (rho = 0.36 historical, 0.45 future,
0.45 vs SHR6); and the authors' own code performs this join.

THE SUPPORT ASYMMETRY, WHICH IS THE CENTRAL LIMITATION. The paper seeds 330 embryos per
event = 6 `w_perc` (horizontal launch position, 0.2..1.0, larger is closer to the
convective core) x 5 `r_int` (embryo radius) x 11 `h_int` (launch height). The deposited
FUTURE radii carry all 330. The deposited HISTORICAL radii carry 55 -- `w_perc = 0.5`
only -- although the historical growth-duration files carry all six positions and
`severity_calc_n` hard-codes `ndia = 6*5*11`, so the full historical run exists and was
not archived. Therefore every comparison here is made on the COMMON SUPPORT `w_perc=0.5`,
which is a reduced experiment, not the published configuration. Reducing over all six
future positions while the historical side has one inflates the median diameter shift from
+30% to +97%: that comparison is not an estimate of anything. Assertions below refuse to
run if the common-support rule is violated.

UNITS ARE EXTERNAL, NOT SELF-DESCRIBED. The NetCDFs carry no `units`, `long_name` or
global attributes at all. Metres and RADIUS come from the figshare description and are
corroborated by the authors' code (`severity_calc_n(nc['r_hail']*2000.)`, diameter in mm).
`r_hail` is the radius at the surface and `r_grow` the radius at the melting level (their
code comment); note `r_hail > r_grow` in 2,506 historical embryo runs, so treat that
ordering as data, not an invariant.

ZEROS ARE PHYSICS, NOT MISSINGNESS. `r_hail == 0` means the embryo melted completely
before landing (~50% of historical embryo runs). NaN is separate and rare (1e-4). Never
conflate them, and never turn an absent event into a zero.

THREE DIFFERENT ESTIMANDS GET CONFLATED, SO ALL THREE ARE REPORTED SEPARATELY:
  * embryo share      -- fraction of ALL seeded embryos landing >= 30 mm. This is the
                         paper's headline "frequency of >=30 mm hailstones".
  * produced share    -- the same fraction among embryos that produced a stone at all.
  * event share       -- fraction of EVENTS whose maximum stone reaches 30 mm.
They move very differently (+46.0%, +21.3%, +19.5% on ssp245 member 1), so a number quoted
without its denominator is unusable.

WHY IT IS NOT A LAYER (and why this script writes no NetCDF). The deposit has no time
axis, no decades, no trend and no 2020s baseline -- one ERA5 control against one 2071-2100
perturbation per scenario. It is event-conditioned: the events are satellite detections,
so spatial density is a sampling property of GPM, not a hail climatology, and simulation
eligibility is itself climate-dependent (the model requires a peak updraft above a
threshold, so events enter and leave the population under warming; the transition counts
are reported below). It therefore cannot express hail FREQUENCY, only conditional
severity, and at 12,412 events a regular global grid averages 4.8 events per cell at 5
degrees. Any product built from this belongs outside the OUTPUT-SPEC contract.

USAGE
    python3 scripts/measure_hail_nature2026.py --raw-dir data/raw/hail-nature2026
    python3 scripts/measure_hail_nature2026.py --members 1 2 3 --scenarios ssp245

Requires `py7zr` (the deposit ships nested .7z archives and this machine has no 7z binary).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import py7zr
import xarray as xr

TYPES = (1, 2, 3, 4, 5)
COMMON_W_PERC = 0.5           # the only launch position present on both sides
N_EMBRYO_COMMON = 5 * 11      # embryos per event on the common support
N_EMBRYO_PAPER = 6 * 5 * 11   # what the publication used
MM_PER_RADIUS_M = 2000.0      # radius in metres -> diameter in millimetres
LARGE_MM = 30.0               # the publication's large-hail threshold

HIST_PROFILE = "Profile_of_meterological_variables_in_hailstorm_environments_new_w_params_adiab_dir.nc"

# EC-Earth3 realization order as declared in the authors' post-processing script; member N
# in the archive layout is MEMBER_IDS[N-1].
MEMBER_IDS = [
    "r140i1p1f1", "r119i1p1f1", "r126i1p1f1", "r114i1p1f1", "r107i1p1f1",
    "r143i1p1f1", "r102i1p1f1", "r108i1p1f1", "r113i1p1f1", "r133i1p1f1",
    "r132i1p1f1", "r103i1p1f1", "r128i1p1f1", "r145i1p1f1", "r142i1p1f1",
    "r111i1p1f1", "r147i1p1f1", "r124i1p1f1", "r118i1p1f1", "r141i1p1f1",
]


def kinetic_energy_j(diameter_mm: np.ndarray) -> np.ndarray:
    """Per-stone kinetic energy, verbatim from the authors' `severity_calc_n`.

    Piecewise, and the break at 22.1 mm matters: quoting only the 3.5-power branch
    misstates every stone below it.
    """
    d = diameter_mm
    return np.where(
        np.isnan(d),
        np.nan,
        np.where(d < 22.1, 0.0107 * (d / 10.0) ** 4.59, 0.0243 * (d / 10.0) ** 3.5),
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _extract(archive: Path, targets: list[str], dest: Path) -> None:
    with py7zr.SevenZipFile(archive, "r") as z:
        z.extract(path=str(dest), targets=targets)


def hist_event_index(raw_dir: Path, work: Path) -> dict[int, np.ndarray]:
    """Catalogue row indices of the events simulated under the historical (ERA5) run."""
    target = work / HIST_PROFILE
    if not target.exists():
        _extract(raw_dir / "Soundings_to_drive_model_simulations.7z", [HIST_PROFILE], work)
    with xr.open_dataset(target) as ds:
        return {t: ds[f"date{t}"].values.astype(int) for t in TYPES}


def fut_event_index(raw_dir: Path, work: Path, scenario: str, member: int) -> dict[int, np.ndarray]:
    mid = MEMBER_IDS[member - 1]
    name = f"{scenario}/Profile_of_meterological_variables_CMIP6_EC-Earth3_{scenario}_{mid}_PGW_adiab.nc"
    target = work / name
    if not target.exists():
        _extract(raw_dir / "Soundings_to_drive_model_simulations.7z", [name], work)
    with xr.open_dataset(target) as ds:
        return {t: ds[f"date{t}"].values.astype(int) for t in TYPES}


def hist_radii(raw_dir: Path, work: Path, t: int) -> np.ndarray:
    name = f"hailstone_growth_radii_of_type{t}_hist.nc"
    target = work / "hist" / name
    if not target.exists():
        _extract(raw_dir / "HIST-simulation_results.7z", [name], work / "hist")
    with xr.open_dataset(target) as ds:
        w = np.atleast_1d(ds["w_perc"].values)
        if w.size != 1 or not np.isclose(w[0], COMMON_W_PERC):
            raise SystemExit(
                f"historical type {t} carries w_perc={w}; the common-support rule assumes "
                f"the deposited historical file is the {COMMON_W_PERC} slice. If the full "
                "six-position historical run has been obtained, this script must be "
                "updated to select w_perc on BOTH sides rather than assume it."
            )
        arr = ds["r_hail"].values * MM_PER_RADIUS_M
    if arr.ndim != 3:
        raise SystemExit(f"historical type {t} r_hail has dims {arr.shape}, expected (sample, r_int, h_int)")
    return arr


def fut_radii(raw_dir: Path, work: Path, scenario: str, member: int, t: int) -> np.ndarray:
    mid = MEMBER_IDS[member - 1]
    outer = work / scenario
    inner = outer / f"{scenario}_member{member}.7z"
    if not inner.exists():
        _extract(raw_dir / f"{scenario}-simulation_results.7z", [f"{scenario}_member{member}.7z"], outer)
    name = f"member{member}/hailstone_growth_radii_of_type{t}_EC-Earth3_{scenario}_{mid}_FUT.nc"
    target = outer / name
    if not target.exists():
        _extract(inner, [name], outer)
    with xr.open_dataset(target) as ds:
        if "w_perc" not in ds["r_hail"].dims:
            raise SystemExit(f"{scenario} member {member} type {t}: r_hail has no w_perc dimension")
        arr = ds["r_hail"].sel(w_perc=COMMON_W_PERC).values * MM_PER_RADIUS_M
    return arr


def ensemble_mean_radii(raw_dir: Path, work: Path, gcm: str, scenario: str, t: int) -> np.ndarray:
    """Ensemble-mean future radii. EC-Earth3, MPI-ESM1-2-LR and NorESM2-LM are all present
    here, which is the only place structural spread across driving models can be measured."""
    name = (f"simulation results/{gcm}/{scenario.upper()}/"
            f"hailstone_growth_radii_of_type{t}_{gcm}_{scenario}_ensemble_mean_FUT.nc")
    target = work / name
    if not target.exists():
        _extract(raw_dir / "ensemble_mean_simulations.7z", [name], work)
    with xr.open_dataset(target) as ds:
        return ds["r_hail"].sel(w_perc=COMMON_W_PERC).values * MM_PER_RADIUS_M


def reproduction_table(raw_dir: Path, work: Path) -> list[dict]:
    """Reproduce the publication's headline numbers.

    UNPAIRED and over ALL simulated events on each side, because that is what the paper
    reports: its quoted ranges span GCM x scenario, not scenario alone. Pairing on common
    support (the `run_member` path) is the right basis for a PRODUCT claim but is a
    different, slightly higher-signal population.
    """
    hist = [hist_radii(raw_dir, work, t) for t in TYPES]
    h_flat = np.concatenate([h.ravel() for h in hist])
    ake_h = float(np.nanmean(np.concatenate(
        [np.nansum(kinetic_energy_j(h), axis=(1, 2)) / N_EMBRYO_COMMON for h in hist])))
    rows = []
    for gcm in ("EC-Earth3", "MPI-ESM1-2-LR", "NorESM2-LM"):
        for scenario in ("ssp245", "ssp370", "ssp585"):
            try:
                fut = [ensemble_mean_radii(raw_dir, work, gcm, scenario, t) for t in TYPES]
            except Exception as exc:  # noqa: BLE001 - report and continue
                print(f"  [skip] {gcm} {scenario}: {type(exc).__name__}")
                continue
            f_flat = np.concatenate([f.ravel() for f in fut])
            ake_f = float(np.nanmean(np.concatenate(
                [np.nansum(kinetic_energy_j(f), axis=(1, 2)) / N_EMBRYO_COMMON for f in fut])))
            rows.append({
                "gcm": gcm, "scenario": scenario,
                "ge30_hist": float(np.mean(h_flat >= LARGE_MM)),
                "ge30_fut": float(np.mean(f_flat >= LARGE_MM)),
                "lt30_hist": float(np.mean(h_flat < LARGE_MM)),
                "lt30_fut": float(np.mean(f_flat < LARGE_MM)),
                "ake_hist_j": ake_h, "ake_fut_j": ake_f,
            })
    return rows


def estimands(hist_mm: np.ndarray, fut_mm: np.ndarray) -> dict[str, float]:
    """All estimands on paired events, common support. Denominators are named."""
    h_evt = np.nanmax(hist_mm, axis=(1, 2))
    f_evt = np.nanmax(fut_mm, axis=(1, 2))
    h_flat, f_flat = hist_mm.ravel(), fut_mm.ravel()
    h_prod, f_prod = h_flat[h_flat > 0], f_flat[f_flat > 0]
    ake_h = np.nanmean(np.nansum(kinetic_energy_j(hist_mm), axis=(1, 2)) / N_EMBRYO_COMMON)
    ake_f = np.nanmean(np.nansum(kinetic_energy_j(fut_mm), axis=(1, 2)) / N_EMBRYO_COMMON)
    return {
        "n_paired_events": float(hist_mm.shape[0]),
        "embryo_melted_share_hist": float(np.mean(h_flat == 0)),
        "embryo_melted_share_fut": float(np.mean(f_flat == 0)),
        "embryo_share_ge30_hist": float(np.mean(h_flat >= LARGE_MM)),
        "embryo_share_ge30_fut": float(np.mean(f_flat >= LARGE_MM)),
        "embryo_share_lt30_hist": float(np.mean(h_flat < LARGE_MM)),
        "embryo_share_lt30_fut": float(np.mean(f_flat < LARGE_MM)),
        "produced_share_ge30_hist": float(np.mean(h_prod >= LARGE_MM)),
        "produced_share_ge30_fut": float(np.mean(f_prod >= LARGE_MM)),
        "event_share_ge30_hist": float(np.mean(h_evt >= LARGE_MM)),
        "event_share_ge30_fut": float(np.mean(f_evt >= LARGE_MM)),
        "event_share_nohail_hist": float(np.mean(h_evt == 0)),
        "event_share_nohail_fut": float(np.mean(f_evt == 0)),
        "event_mhd_median_hist": float(np.median(h_evt)),
        "event_mhd_median_fut": float(np.median(f_evt)),
        "paired_mhd_delta_median": float(np.median(f_evt - h_evt)),
        "paired_mhd_delta_mean": float(np.mean(f_evt - h_evt)),
        "paired_share_increasing": float(np.mean(f_evt > h_evt)),
        "ake_hist_j": float(ake_h),
        "ake_fut_j": float(ake_f),
    }


def run_member(raw_dir: Path, work: Path, scenario: str, member: int,
               hidx: dict[int, np.ndarray]) -> dict[str, float]:
    fidx = fut_event_index(raw_dir, work, scenario, member)
    H, F = [], []
    n_hist = n_fut = 0
    for t in TYPES:
        hr = hist_radii(raw_dir, work, t)
        fr = fut_radii(raw_dir, work, scenario, member, t)
        if hr.shape[0] != hidx[t].size:
            raise SystemExit(f"type {t}: historical radii {hr.shape[0]} rows vs profile index {hidx[t].size}")
        if fr.shape[0] != fidx[t].size:
            raise SystemExit(f"{scenario} m{member} type {t}: future radii {fr.shape[0]} rows vs profile index {fidx[t].size}")
        keep_h = np.isin(hidx[t], fidx[t])
        keep_f = np.isin(fidx[t], hidx[t])
        # np.isin preserves each side's own order; both sides are sorted ascending, so the
        # kept subsequences correspond element-for-element. Refuse to guess if they do not.
        if not (np.all(np.diff(hidx[t]) > 0) and np.all(np.diff(fidx[t]) > 0)):
            raise SystemExit(f"type {t}: event indices are not strictly increasing; the pairing assumption fails")
        H.append(hr[keep_h])
        F.append(fr[keep_f])
        n_hist += hidx[t].size
        n_fut += fidx[t].size
    out = estimands(np.concatenate(H), np.concatenate(F))
    out.update({
        "scenario_n_hist_simulated": float(n_hist),
        "scenario_n_fut_simulated": float(n_fut),
        "eligibility_hist_only": float(n_hist - out["n_paired_events"]),
        "eligibility_fut_only": float(n_fut - out["n_paired_events"]),
    })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=Path("data/raw/hail-nature2026"))
    ap.add_argument("--scenarios", nargs="+", default=["ssp245", "ssp370", "ssp585"])
    ap.add_argument("--members", nargs="+", type=int, default=list(range(1, 21)))
    ap.add_argument("--out", type=Path, default=Path("reports/maps/hail-severity/gate.csv"))
    ap.add_argument("--keep-work", action="store_true", help="do not delete the extraction scratch dir")
    ap.add_argument("--reproduction-only", action="store_true",
                    help="run only the GCM x scenario reproduction of the published headline numbers")
    args = ap.parse_args()

    raw = args.raw_dir
    if not raw.is_dir():
        print(f"raw dir not found: {raw}", file=sys.stderr)
        return 2

    print("provenance (sha256, first 16):")
    for p in sorted(raw.iterdir()):
        if p.is_file() and p.suffix in {".7z", ".nc", ".zip", ".docx"}:
            print(f"  {sha256(p)[:16]}  {p.stat().st_size:>13,}  {p.name}")

    work = Path(tempfile.mkdtemp(prefix="hail_gate_"))
    try:
        if (raw / "ensemble_mean_simulations.7z").exists():
            print("\nREPRODUCTION OF THE PUBLISHED HEADLINES (ensemble-mean runs, unpaired,")
            print(f"w_perc={COMMON_W_PERC} common support -- {N_EMBRYO_COMMON} of the publication's {N_EMBRYO_PAPER} embryos):")
            print("  %-16s %-8s %-22s %-22s %-22s" % ("GCM", "scenario", "embryo >=30 mm", "embryo <30 mm", "mean kinetic energy"))
            repro = reproduction_table(raw, work)
            for r in repro:
                print("  %-16s %-8s %.3f->%.3f %+6.1f%%   %.3f->%.3f %+6.1f%%   %.3f->%.3f %+6.1f%%" % (
                    r["gcm"], r["scenario"],
                    r["ge30_hist"], r["ge30_fut"], 100 * (r["ge30_fut"] / r["ge30_hist"] - 1),
                    r["lt30_hist"], r["lt30_fut"], 100 * (r["lt30_fut"] / r["lt30_hist"] - 1),
                    r["ake_hist_j"], r["ake_fut_j"], 100 * (r["ake_fut_j"] / r["ake_hist_j"] - 1)))
            if repro:
                g = [100 * (r["ge30_fut"] / r["ge30_hist"] - 1) for r in repro]
                l = [100 * (r["lt30_fut"] / r["lt30_hist"] - 1) for r in repro]
                k = [100 * (r["ake_fut_j"] / r["ake_hist_j"] - 1) for r in repro]
                print(f"  RANGE over GCM x scenario:  >=30 mm [{min(g):+.1f}, {max(g):+.1f}]   "
                      f"<30 mm [{min(l):+.1f}, {max(l):+.1f}]   KE [{min(k):+.1f}, {max(k):+.1f}]")
                print( "  PUBLISHED:                  >=30 mm [+37.9, +51.8]   <30 mm [-4.2, -12.3]   "
                       "damage potential [+36.5, +42.1]")
                print("  The >=30 mm range matches. <30 mm and kinetic energy run more negative / higher")
                print("  than published, which is the expected direction if the missing five launch")
                print("  positions matter -- see the module docstring and reports/maps/hail-severity/author_query.md.")
            if args.reproduction_only:
                return 0

        hidx = hist_event_index(raw, work)
        print(f"\nhistorical events simulated: {sum(v.size for v in hidx.values()):,} of 12,412 catalogue events")

        rows = []
        for scenario in args.scenarios:
            if not (raw / f"{scenario}-simulation_results.7z").exists():
                print(f"  [skip] {scenario}: archive not present")
                continue
            for member in args.members:
                try:
                    rec = run_member(raw, work, scenario, member, hidx)
                except SystemExit as exc:
                    print(f"  [fail] {scenario} member {member}: {exc}")
                    continue
                rec.update({"scenario": scenario, "member": member, "member_id": MEMBER_IDS[member - 1]})
                rows.append(rec)
                print(f"  {scenario} m{member:<2} "
                      f"paired={rec['n_paired_events']:>6.0f} "
                      f"embryo>=30 {rec['embryo_share_ge30_hist']:.3f}->{rec['embryo_share_ge30_fut']:.3f} "
                      f"({100*(rec['embryo_share_ge30_fut']/rec['embryo_share_ge30_hist']-1):+.1f}%)  "
                      f"AKE {rec['ake_hist_j']:.3f}->{rec['ake_fut_j']:.3f} J "
                      f"({100*(rec['ake_fut_j']/rec['ake_hist_j']-1):+.1f}%)")

        if not rows:
            print("no members processed", file=sys.stderr)
            return 1

        args.out.parent.mkdir(parents=True, exist_ok=True)
        cols = ["scenario", "member", "member_id"] + [k for k in rows[0] if k not in {"scenario", "member", "member_id"}]
        with args.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.out} ({len(rows)} member-runs)")

        print("\nENSEMBLE RANGE across members -- this is initial-condition spread of ONE model,")
        print("not a confidence interval, and it excludes structural uncertainty entirely.")
        for scenario in sorted({r["scenario"] for r in rows}):
            sub = [r for r in rows if r["scenario"] == scenario]
            for label, num, den in (
                ("embryo share >=30 mm", "embryo_share_ge30_fut", "embryo_share_ge30_hist"),
                ("embryo share <30 mm", "embryo_share_lt30_fut", "embryo_share_lt30_hist"),
                ("mean kinetic energy", "ake_fut_j", "ake_hist_j"),
            ):
                pct = np.array([100 * (r[num] / r[den] - 1) for r in sub])
                print(f"  {scenario}  {label:<22} {np.median(pct):+6.1f}%  "
                      f"[{pct.min():+.1f}, {pct.max():+.1f}] over n={len(sub)} members")
        print("\nPublished for comparison: >=30 mm +37.9 to +51.8%, <30 mm -4.2 to -12.3%,")
        print(f"damage potential +36.5 to +42.1% -- all computed on {N_EMBRYO_PAPER} embryos per event,")
        print(f"whereas the deposited historical file supports only {N_EMBRYO_COMMON}. Agreement here is")
        print("corroboration on reduced support, NOT a reproduction of the published figures.")
        return 0
    finally:
        if args.keep_work:
            print(f"work dir kept: {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
