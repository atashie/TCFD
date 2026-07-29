"""Generate a QA/QC report for a published TCFD layer and publish it to qa/.

Runs the invariant checks that the 6-value-class contract depends on, plus
per-scenario summary statistics, and writes both a machine-readable
``qa_report.json`` and a human-readable ``qa_report.html`` into the layer
version's ``qa/`` prefix.

The checks are driven by each file's OWN declared attributes (``trend_definition``,
``percentile_direction``, ``baseline_source``) rather than hardcoded per variable,
so the same report serves every annualized layer (led, let, burntarea, csoil, ...).

``qa/`` is deliberately OUTSIDE the ``_COMPLETE.json`` gate (STORAGE.md): the report
is regenerable evidence pinned to its data, so re-running this never invalidates a
published layer.

Usage:
    python scripts/generate_qa_report.py soilcarbon_csoil_annual
    python scripts/generate_qa_report.py drought_led_annual --version v2026-07-24_abc1234
    python scripts/generate_qa_report.py cyclone_let_annual --local-only
"""

import argparse
import html
import json
import re
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.stats import rankdata   # average-rank ties, for the Spearman direction check

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "isimip-pipeline" / "src"))
from isimip_pipeline import storage  # noqa: E402

VALUE_CLASSES = ["median", "percentile", "trend", "lower_ci", "upper_ci"]
TOL = 1e-4          # absolute tolerance for the trend<->change identity
CI_TOL = 1e-6       # float32 slack for CI ordering comparisons


def log(msg):
    print(msg, flush=True)


class Check:
    """One pass/fail assertion plus the evidence behind it."""

    def __init__(self, results):
        self.results = results

    def __call__(self, scope, name, passed, detail="", severity="error"):
        self.results.append({
            "scope": scope, "check": name, "passed": bool(passed),
            "detail": str(detail), "severity": severity,
        })
        mark = "PASS" if passed else ("WARN" if severity == "warning" else "FAIL")
        log(f"  {mark:4s}  [{scope}] {name}" + (f"  -> {detail}" if detail else ""))
        return bool(passed)


def _composition_of(attrs, scenario):
    """Return the declared ensemble composition signature for one scenario, or "".

    Reads the layer's own `members_by_scenario` attribute, written as
    ``"rcp26=[member,member,...]; rcp60=[...]"`` -- member IDENTITY, not a count. Two
    scenarios can hold the same NUMBER of members without holding the same ones (clm45
    contributes different GCM pairs to different RCPs), and a count-based signature would
    then demand bit-identity between baseline panels that legitimately differ.

    Layers that do not declare the attribute (every layer with a uniform ensemble) get ""
    for every scenario, so all scenarios fall into a single group and the strict
    cross-scenario baseline identity check behaves exactly as before.
    """
    spec = str(attrs.get("members_by_scenario", ""))
    for part in spec.split(";"):
        name, _, sig = part.strip().partition("=")
        if name.strip() == scenario and sig:
            return sig.strip()
    return ""


def _finite(a):
    return a[np.isfinite(a)]


def seam_spacing(a, axis=1, min_share=0.4):
    """Most common spacing between gradient 'seams', in cells, or None if smooth.

    Finds columns (or rows) whose mean |gradient| far exceeds its local median, then
    reports the modal gap between them. Unlike grid_periodicity() this assumes neither a
    candidate block width nor alignment to the array origin, and it works on blocks that
    are merely SMOOTH inside rather than exactly constant -- the case that defeated both
    the exact-tie and fixed-modulo tests when a ~4x5 deg member was hiding in a 0.5 deg
    ensemble.

    Returns ``(spacing_cells, n_gaps_at_that_spacing, n_gaps_total)`` when one spacing
    accounts for at least ``min_share`` of the gaps, else None.
    """
    import collections
    g = np.abs(np.diff(a, axis=axis))
    with np.errstate(invalid="ignore"):
        prof = np.nanmean(g, axis=1 - axis)
    if not np.isfinite(prof).any() or np.nanmean(prof) == 0:
        return None
    prof = prof / np.nanmean(prof)
    w = 9
    pad = np.pad(prof, (w // 2, w // 2), mode="edge")
    local = np.nanmedian(np.lib.stride_tricks.sliding_window_view(pad, w), axis=1)
    peaks = np.flatnonzero(np.isfinite(prof) & np.isfinite(local) & (prof > 2.0 * local))
    if peaks.size < 5:
        return None
    gaps = np.diff(peaks)
    gaps = gaps[gaps > 1]
    if gaps.size < 4:
        return None
    spacing, count = collections.Counter(gaps.tolist()).most_common(1)[0]
    if count / gaps.size < min_share:
        return None
    return int(spacing), int(count), int(gaps.size)


def grid_periodicity(a, k):
    """Ratio of mean |gradient| on columns/rows INSIDE a k-block vs on its seams.

    A model that runs natively coarser than 0.5 deg and was replicated onto the ISIMIP
    grid keeps the declared 360x720 dims, so a dimension check cannot see it — but its
    values are constant within each k x k block, leaving zero gradient inside a block
    and all the gradient on the seams. Ratio ~1 means genuinely 0.5 deg; ratio well
    below 1 means a coarse component at k*0.5 deg is present.

    Returns ``(lon_ratio, lat_ratio)``. Detects a coarse component even when it is
    diluted by finer members in an ensemble mean, which is why it is worth running on
    the published field rather than trusting the raw files' declared grid.
    """
    out = []
    for axis in (1, 0):
        g = np.abs(np.diff(a, axis=axis))
        with np.errstate(invalid="ignore"):
            prof = np.nanmean(g, axis=1 - axis)
        idx = np.arange(prof.size)
        inside, seam = prof[idx % k != 0], prof[idx % k == 0]
        inside, seam = inside[np.isfinite(inside)], seam[np.isfinite(seam)]
        if inside.size < 5 or seam.size < 5 or np.mean(seam) == 0:
            out.append(float("nan"))
        else:
            out.append(float(np.mean(inside) / np.mean(seam)))
    return out[0], out[1]


def summarize(ds, var, decades):
    """Per-decade summary statistics for one value class."""
    out = {}
    a = ds[var].values
    for i, d in enumerate(decades):
        f = _finite(a[i])
        if f.size == 0:
            out[str(d)] = None
            continue
        out[str(d)] = {
            "n_finite": int(f.size),
            "min": round(float(f.min()), 6),
            "p05": round(float(np.percentile(f, 5)), 6),
            "median": round(float(np.median(f)), 6),
            "mean": round(float(f.mean()), 6),
            "p95": round(float(np.percentile(f, 95)), 6),
            "max": round(float(f.max()), 6),
        }
    return out


def check_layer(ds_by_scen, layer_id, version):
    results = []
    chk = Check(results)
    scenarios = sorted(ds_by_scen)
    first = ds_by_scen[scenarios[0]]
    attrs = dict(first.attrs)
    decades = [int(d) for d in first.decade.values]
    declared_baseline = attrs.get("baseline_decade")
    baseline = int(declared_baseline) if declared_baseline is not None else decades[0]
    # If a file DECLARES a baseline that is not on its own decade axis, silently
    # falling back to index 0 would test a different decade while every label still
    # said the declared one -- a passing report about the wrong slice.
    baseline_ok = baseline in decades
    b_idx = decades.index(baseline) if baseline_ok else 0
    # Match hyphen OR underscore: processors write the prose form
    # "baseline-anchored rate: ...", the manifest uses baseline_anchored. Testing
    # only one spelling silently SKIPS the identity check below, which is worse
    # than failing it.
    trend_def = str(attrs.get("trend_definition", "")).replace("-", "_")
    anchored = "baseline_anchored" in trend_def
    hib = attrs.get("percentile_direction", "") == "higher_is_better"
    shared_baseline = "shared" in str(attrs.get("baseline_source", ""))

    log(f"\nLayer {layer_id} @ {version}")
    log(f"  scenarios={scenarios}  decades={decades}  baseline={baseline}s")
    log(f"  percentile_direction={attrs.get('percentile_direction','?')}  "
        f"trend={'baseline-anchored' if anchored else attrs.get('trend_definition','?')}")

    log("\n--- structure ---")
    chk("layer", "declared baseline_decade exists on the decade axis", baseline_ok,
        f"declared={declared_baseline!r}, decades={decades}"
        + ("" if baseline_ok else f" -- FALLING BACK to {decades[0]}s; every "
                                  f"baseline check below targets that instead"))
    for s in scenarios:
        ds = ds_by_scen[s]
        missing = [v for v in VALUE_CLASSES if v not in ds.data_vars]
        chk(s, "all value classes present", not missing, f"missing={missing}")
        chk(s, "dims are (decade, lat, lon)",
            ds["median"].dims == ("decade", "lat", "lon"), str(ds["median"].dims))
        chk(s, "decade axis matches",
            [int(d) for d in ds.decade.values] == decades,
            str([int(d) for d in ds.decade.values]))

    log("\n--- value-class invariants ---")
    for s in scenarios:
        ds = ds_by_scen[s]
        m, lo, hi = (ds[k].values for k in ("median", "lower_ci", "upper_ci"))
        fin = np.isfinite(m) & np.isfinite(lo) & np.isfinite(hi)
        viol = int(np.sum((lo[fin] > m[fin] + CI_TOL) | (m[fin] > hi[fin] + CI_TOL)))
        chk(s, "CI ordering lower <= median <= upper", viol == 0, f"violations={viol}")

        # A zero-width CI is legitimate where every contributing member reports
        # exactly 0 (barren desert/ice), and suspect anywhere else -- report the
        # split rather than one ambiguous count.
        collapsed = np.isclose(lo[fin], hi[fin])
        at_zero = collapsed & np.isclose(m[fin], 0.0)
        n_all, n_zero = int(collapsed.sum()), int(at_zero.sum())
        extra = ""
        if "n_models" in ds.data_vars and n_all > n_zero:
            nmo_f = ds["n_models"].values[fin]
            resid = collapsed & ~at_zero
            single = int(np.sum(nmo_f[resid] <= 1))
            extra = (f"; of the {n_all - n_zero:,} non-zero ones, {single:,} are "
                     f"single-model cells (a lone model's GCMs agreeing exactly)")
        chk(s, "CI zero-width only on all-zero cells", n_all == n_zero,
            f"{n_all:,} zero-width of {int(fin.sum()):,} finite "
            f"({100*n_all/max(int(fin.sum()),1):.3f}%); {n_zero:,} are all-zero "
            f"(expected), {n_all - n_zero:,} are not{extra}", severity="warning")

        p = ds["percentile"].values
        pf = _finite(p)
        chk(s, "percentile within [1, 100]",
            bool(pf.min() >= 1 - CI_TOL and pf.max() <= 100 + CI_TOL),
            f"[{pf.min():.3f}, {pf.max():.3f}]")

        # Declared direction must match the realized sign of corr(value, risk pct).
        # RANK correlation, not Pearson: percentile is a percentile-of-score, i.e. a
        # monotone rank transform of the value, so the relationship is monotone but
        # deliberately NON-LINEAR. On a heavy-tailed, zero-inflated variable Pearson
        # therefore reads far below 1 even when the mapping is perfectly ordered --
        # burntarea (45% exact zeros, tail past 100%) scored +0.53 and FAILed a
        # correct layer. Spearman is ~1 by construction and still catches a genuinely
        # inverted or scrambled percentile, which is what this check is for.
        b, bp = m[b_idx], p[b_idx]
        g = np.isfinite(b) & np.isfinite(bp)
        r = float("nan")
        if int(g.sum()) > 2:
            bv, pv = b[g], bp[g]
            if np.ptp(bv) == 0 or np.ptp(pv) == 0:
                r = float("nan")            # degenerate: no ordering to verify
            else:
                r = float(np.corrcoef(rankdata(bv), rankdata(pv))[0, 1])
        if hib:
            chk(s, "percentile INVERTED (higher_is_better)", r < -0.9,
                f"spearman(median, percentile)={r:+.4f}, expected strongly negative")
        else:
            chk(s, "percentile aligned (higher_is_worse)", r > 0.9,
                f"spearman(median, percentile)={r:+.4f}, expected strongly positive")

        t = ds["trend"].values
        chk(s, f"trend == 0 in baseline decade ({baseline}s)",
            bool(np.all(_finite(t[b_idx]) == 0)),
            f"max|trend[{baseline}]|={np.nanmax(np.abs(t[b_idx])):.3e}")

        if anchored:
            worst, n_compared = 0.0, 0
            for i in range(len(decades)):
                span = i - b_idx
                if span == 0:
                    continue
                dev = np.abs(t[i] * span - (m[i] - m[b_idx]))
                df = _finite(dev)
                if df.size:
                    worst = max(worst, float(df.max()))
                    n_compared += int(df.size)
            # worst stays 0.0 if NOTHING was comparable, which would report PASS
            # having compared nothing -- the same vacuous-pass class as the
            # hyphen/underscore bug above. Require actual comparisons.
            chk(s, "trend x elapsed_decades == change map",
                n_compared > 0 and worst < TOL,
                f"max abs deviation={worst:.3e} over {n_compared:,} cells"
                + ("" if n_compared else "  -- NOTHING COMPARABLE: trend and median "
                                         "finite masks may be disjoint"))
        else:
            # Never let an unrecognized trend definition quietly drop the check.
            chk(s, "trend<->change identity CHECKED", False,
                f"skipped: trend_definition not recognized as baseline-anchored "
                f"({trend_def[:60]!r})", severity="warning")

        if "n_members" in ds.data_vars:
            nm = ds["n_members"].values
            nmf = _finite(nm)
            declared = attrs.get("n_members")
            chk(s, "n_members never exceeds declared ensemble size",
                bool(nmf.max() <= float(declared)) if declared else True,
                f"max={nmf.max():.0f} declared={declared}")
            chk(s, "coverage counts align with finite data",
                int(np.sum(np.isfinite(nm[b_idx]))) == int(np.sum(np.isfinite(m[b_idx]))),
                f"counts={int(np.sum(np.isfinite(nm[b_idx]))):,} "
                f"data={int(np.sum(np.isfinite(m[b_idx]))):,}")

    # --- zonal profile -------------------------------------------------------
    # A global or area-weighted statistic cannot see a defect confined to one latitude
    # zone: polar cells carry almost no area, so a member projecting absurd values there
    # passes every aggregate check. The wildfire layer's `visit` members do exactly that
    # (~100%/yr burnt area on Arctic islands). Always REPORT the profile so a reviewer can
    # see the shape, and warn when a layer that declares itself low-latitude-dominated
    # has a polar band exceeding its tropical band.
    log("\n--- zonal profile (land-mean of median by latitude band) ---")
    zonal = {}
    for s in scenarios:
        ds = ds_by_scen[s]
        lat = ds["lat"].values
        a = ds["median"].values[-1]              # last decade = strongest signal
        bands = [(-90, -60), (-60, -23), (-23, 23), (23, 45), (45, 60),
                 (60, 70), (70, 75), (75, 90)]
        prof = {}
        for lo, hi in bands:
            sel = (lat >= lo) & (lat < hi)
            f = _finite(a[sel])
            prof[f"{lo}..{hi}"] = round(float(f.mean()), 4) if f.size else None
        zonal[s] = prof
        log(f"  [{s}] {decades[-1]}s: " +
            "  ".join(f"{k}={'--' if v is None else f'{v:.3f}'}" for k, v in prof.items()))

        if str(attrs.get("zonal_expectation", "")).startswith("low_latitude"):
            trop = prof.get("-23..23")
            polar = max((prof.get("70..75") or 0.0), (prof.get("75..90") or 0.0))
            chk(s, "polar band does not exceed tropical band",
                not (trop and polar > trop),
                f"tropics(-23..23)={trop}, worst polar band={polar:.3f} "
                f"-- layer declares zonal_expectation={attrs.get('zonal_expectation')!r}",
                severity="warning")

        # --- latitude-seam detection -----------------------------------------
        # A band profile is too coarse to see a discontinuity that falls INSIDE a band,
        # and a smooth geophysical field should not jump sharply between two adjacent
        # 0.5 deg rows. A sharp step usually means the INPUT changed source there --
        # `fldfrc` halves across exactly 60.0N (11.06% -> 5.93%, 1.87x) because
        # CaMa-Flood's floodplain topography switches DEM at the SRTM/HydroSHEDS
        # coverage limit. That is a property of the source, not of the processing, but
        # it silently biases every zonal or regional aggregate that straddles it.
        # Scored as an outlier against the field's own row-to-row variability so the
        # check is scale-free and needs no per-layer threshold.
        # Only compare rows whose zonal mean rests on enough cells to be stable. Without
        # this the check is dominated by the poles and the far south, where a row holds a
        # dozen cells and its mean swings wildly for reasons that are not seams: it
        # false-fired at -50.75 (12-19 cells), 83.75 N (19-32) and 80.25 N (97-236) on
        # three existing layers. The genuine `fldfrc` seam sits on rows of 451-475 cells.
        MIN_ROW_CELLS = 150
        counts = np.array([int(_finite(a[i]).size) for i in range(a.shape[0])])
        zrow = np.array([
            (lambda f: float(f.mean()) if f.size else np.nan)(_finite(a[i]))
            for i in range(a.shape[0])])
        well_sampled = counts >= MIN_ROW_CELLS
        ok_rows = (np.isfinite(zrow[:-1]) & np.isfinite(zrow[1:])
                   & well_sampled[:-1] & well_sampled[1:])
        jumps = np.abs(np.diff(zrow))
        if ok_rows.sum() >= 20:
            worst_i = int(np.nanargmax(np.where(ok_rows, jumps, np.nan)))
            # Scored against the LOCAL median row-to-row change (a +/-10-row window,
            # excluding the candidate), not the global one. This is what separates a SEAM
            # from a steep GRADIENT, and the distinction is not academic: on `fldfrc` the
            # global-median form flagged 70 N and the equator as loudly as the real 60 N
            # seam. Checked against the native 150 arcsec data, only 60 N is a
            # discontinuity -- its sharpest native step lands on exactly
            # 60.0208 -> 59.9792 in all three protection levels -- whereas near 70 N and
            # the equator the sharpest native step sits at 70.27 and 0.23, i.e. the
            # profile is simply declining/rising steeply over many rows (fewer Arctic
            # rivers; the Congo-Amazon equatorial belt). A true seam is QUIET on both
            # sides, so a local denominator keeps its ratio huge and collapses a
            # gradient's. A MAD z-score was tried first and rejected: it moved enough
            # between scenarios to put the SAME seam on both sides of a fixed cut
            # (z = 12.2 / 11.7 / 11.7 for rcp26/60/85).
            lo_w, hi_w = max(worst_i - 10, 0), min(worst_i + 11, len(jumps))
            neigh = np.array([jumps[k] for k in range(lo_w, hi_w)
                              if k != worst_i and ok_rows[k]])
            med = float(np.median(neigh)) if neigh.size >= 5 else float(
                np.median(jumps[ok_rows]))
            worst_ratio = jumps[worst_i] / med if med > 0 else float("inf")
            lat_hi, lat_lo = float(lat[worst_i]), float(lat[worst_i + 1])
            step = (zrow[worst_i + 1] / zrow[worst_i]
                    if zrow[worst_i] else float("nan"))
            # Require the step to be MATERIAL as well as anomalous. On a field that is
            # near-zero in the band concerned, the median row-to-row change is tiny and
            # the ratio alone inflates: `driedarea` scored 8-11x on a level change of only
            # 1.25-1.33x, which is not a seam. A real source discontinuity moves the level
            # by a large factor -- `fldfrc` halves (0.52x = 1.92-fold).
            fold = max(step, 1.0 / step) if step and np.isfinite(step) and step > 0 \
                else float("inf")
            declared = str(attrs.get("known_latitude_seams", ""))
            expected = any(
                abs(lat_hi - float(tok)) < 1.0
                for tok in re.findall(r"-?\d+(?:\.\d+)?", declared))
            chk(s, "no undeclared sharp latitude seam in the zonal profile",
                not (worst_ratio >= 9.0 and fold >= 1.5) or expected,
                f"largest single-row jump {lat_hi:.2f}N->{lat_lo:.2f}N: "
                f"{zrow[worst_i]:.4f}->{zrow[worst_i+1]:.4f} ({step:.2f}x level), "
                f"{worst_ratio:.1f}x the LOCAL median row-to-row change ({med:.4f}), "
                f"{fold:.2f}-fold level step"
                + (f"; DECLARED in known_latitude_seams={declared!r}" if expected
                   else "; not declared -- investigate whether the INPUT changes "
                        "source at this latitude, then record it in "
                        "known_latitude_seams"),
                severity="warning")

    if shared_baseline and len(scenarios) > 1:
        # The baseline panel must be identical across scenarios that share the same
        # ENSEMBLE COMPOSITION -- which is every scenario whenever the ensemble is uniform,
        # the case for every layer before timber_*-tempnle. Where a member is missing from
        # some scenario (orchidee-dgvm publishes no rcp85), the baseline must be pooled over
        # that scenario's own members, or differencing baseline against decades manufactures
        # trend (WORKFLOW-ISSUES 2026-07-28; GUARDRAILS S13). Requiring bit-identity there
        # would demand the very defect it is meant to prevent, so scenarios are grouped by
        # the composition the layer declares and identity is asserted WITHIN each group.
        groups = {}
        for s in scenarios:
            groups.setdefault(_composition_of(attrs, s), []).append(s)
        if len(groups) > 1:
            log(f"\n--- shared baseline, grouped by ensemble composition: "
                + "; ".join(f"{sig or 'unspecified'} -> {grp}" for sig, grp in groups.items()))
            log("    (composition varies by scenario, so the baseline panel is pooled per "
                "scenario; identity is required only within a group)")
        else:
            log("\n--- shared baseline (must be identical across scenarios) ---")
        if all(len(grp) < 2 for grp in groups.values()):
            # Every scenario has a distinct member set, so this invariant cannot be tested
            # at all. Say so as a WARNING rather than letting the section pass in silence --
            # a check that quietly tests nothing is worse than one that fails.
            chk("cross-scenario",
                f"{baseline}s baseline identity is testable across scenarios", False,
                f"NOT TESTED: every scenario has a distinct ensemble composition "
                f"({'; '.join(f'{s}->{_composition_of(attrs, s)}' for s in scenarios)}), so "
                f"no two baseline panels are required to match. The per-member baseline is "
                f"still scenario-independent by construction, but that is not visible in the "
                f"published artifact.", severity="warning")
        for sig, grp in groups.items():
            if len(grp) < 2:
                log(f"  [cross-scenario] {grp[0]} is the only scenario with composition "
                    f"{sig or 'unspecified'}; nothing to compare it against")
                continue
            label = "" if len(groups) == 1 else f" within {sig}"
            for k in VALUE_CLASSES:
                ref = ds_by_scen[grp[0]][k].values[b_idx]
                same = all(np.allclose(ref, ds_by_scen[s][k].values[b_idx], equal_nan=True)
                           for s in grp[1:])
                chk("cross-scenario", f"{baseline}s {k} identical across scenarios{label}",
                    same, f"scenarios compared: {grp}")

    log("\n--- land coverage ---")
    for s in scenarios:
        m = ds_by_scen[s]["median"].values
        n = int(np.sum(np.isfinite(m[b_idx])))
        chk(s, "land cells non-empty", n > 0,
            f"{n:,} cells ({100*n/m[b_idx].size:.2f}% of grid)")

    log("\n--- effective spatial resolution ---")
    # Declared dims cannot reveal a natively coarse member replicated onto the 0.5 deg
    # grid; only the values can. Warning, not failure: a coarse member may be a
    # deliberate, documented trade-off for ensemble depth.
    #
    # Test a WIDE range of block widths. A first version of this check stopped at k=4
    # (2 deg) and reported the layer as essentially clean, while the dominant artifact
    # was a member on a ~4x5 deg grid (k=8..10) -- an order of magnitude coarser than
    # anything it looked for. Also report the measured seam SPACING, which needs no
    # candidate list and no assumption that blocks align to the array origin.
    s0 = scenarios[0]
    base = ds_by_scen[s0]["median"].values[b_idx]
    for k in (2, 3, 4, 6, 8, 10, 12):
        rlon, rlat = grid_periodicity(base, k)
        clean = all(np.isnan(r) or r > 0.90 for r in (rlon, rlat))
        chk("grid", f"no {k*0.5:g}-deg block structure in the ensemble", clean,
            f"inside/seam gradient ratio lon={rlon:.3f} lat={rlat:.3f} "
            f"(1.0 = native 0.5 deg; <0.90 = a coarser member is present)",
            severity="warning")
    for axis, name in ((1, "lon"), (0, "lat")):
        sp = seam_spacing(base, axis)
        chk("grid", f"no dominant {name} seam spacing", sp is None,
            "none detected (field is smooth at 0.5 deg)" if sp is None else
            f"seams every {sp[0]} cells = {sp[0]*0.5:g} deg "
            f"({sp[1]} of {sp[2]} gaps) -> a member is effectively that coarse",
            severity="warning")
    if "n_models" in ds_by_scen[s0].data_vars:
        nmo = ds_by_scen[s0]["n_models"].values[b_idx]
        one = int(np.sum(np.isfinite(base) & (nmo == 1)))
        n = int(np.sum(np.isfinite(base)))
        chk("grid", "no cells rendered by a single model alone", one == 0,
            f"{one:,} of {n:,} land cells ({100*one/max(n,1):.2f}%) come from one model, "
            f"so they inherit that model's native resolution unchanged",
            severity="warning")

    return results, dict(scenarios=scenarios, decades=decades, baseline=baseline,
                         attrs=attrs, anchored=anchored, higher_is_better=hib,
                         zonal_profile=zonal)


def build_report(layer_id, version, ds_by_scen, results, meta):
    scenarios = meta["scenarios"]
    decades = meta["decades"]
    stats = {}
    for s in scenarios:
        ds = ds_by_scen[s]
        present = [v for v in VALUE_CLASSES + ["n_members", "n_models"] if v in ds.data_vars]
        stats[s] = {v: summarize(ds, v, decades) for v in present}
        m = ds["median"].values
        g = [float(np.nanmean(m[i])) for i in range(len(decades))]
        stats[s]["_global_mean_by_decade"] = {str(d): round(v, 6) for d, v in zip(decades, g)}
        stats[s]["_global_mean_change_pct"] = (
            round(100 * (g[-1] / g[0] - 1), 4) if g[0] else None)

    n_fail = sum(1 for r in results if not r["passed"] and r["severity"] == "error")
    n_warn = sum(1 for r in results if not r["passed"] and r["severity"] == "warning")
    return {
        "layer_id": layer_id,
        "version": version,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generated_by": "scripts/generate_qa_report.py",
        "summary": {
            "checks_run": len(results),
            "passed": sum(1 for r in results if r["passed"]),
            "failed": n_fail,
            "warnings": n_warn,
            "verdict": "PASS" if n_fail == 0 else "FAIL",
        },
        "layer_attributes": {k: str(v) for k, v in meta["attrs"].items()},
        "checks": results,
        "statistics": stats,
        # Land-mean by latitude band. An aggregate statistic is blind to a defect
        # confined to one zone (polar cells carry almost no area), so the profile is
        # recorded for every layer whether or not it triggered the warning.
        "zonal_profile": meta.get("zonal_profile", {}),
    }


def render_html(rep):
    s = rep["summary"]
    colour = "#2e7d32" if s["verdict"] == "PASS" else "#c62828"
    rows = []
    for r in rep["checks"]:
        if r["passed"]:
            badge, bg = "PASS", "#e8f5e9"
        elif r["severity"] == "warning":
            badge, bg = "WARN", "#fff8e1"
        else:
            badge, bg = "FAIL", "#ffebee"
        rows.append(
            f'<tr style="background:{bg}"><td><b>{badge}</b></td>'
            f'<td>{html.escape(r["scope"])}</td><td>{html.escape(r["check"])}</td>'
            f'<td><code>{html.escape(r["detail"])}</code></td></tr>')

    stat_blocks = []
    for scen, sv in rep["statistics"].items():
        gm = sv.get("_global_mean_by_decade", {})
        chg = sv.get("_global_mean_change_pct")
        heads = "".join(f"<th>{d}s</th>" for d in gm)
        cells = "".join(f"<td>{v:.4f}</td>" for v in gm.values())
        stat_blocks.append(f"<h3>{html.escape(scen)}</h3>")
        if chg is not None:
            stat_blocks.append(
                f"<p>Global-mean change, first to last decade: <b>{chg:+.2f}%</b></p>")
        stat_blocks.append(
            f'<table><tr><th>global mean</th>{heads}</tr>'
            f'<tr><td>median</td>{cells}</tr></table>')

    # Zonal profile: the one view that exposes a defect confined to a latitude zone,
    # which every area-weighted or global statistic is blind to.
    zonal = rep.get("zonal_profile") or {}
    zonal_html = ""
    if zonal:
        bands = list(next(iter(zonal.values())).keys())
        head = "".join(f"<th>{html.escape(b)}&deg;</th>" for b in bands)
        body = []
        for scen, prof in zonal.items():
            vals = [prof.get(b) for b in bands]
            mx = max((v for v in vals if v is not None), default=None)
            tds = "".join(
                "<td>&mdash;</td>" if v is None else
                (f'<td style="background:#ffe0b2"><b>{v:.3f}</b></td>'
                 if mx is not None and v == mx else f"<td>{v:.3f}</td>")
                for v in vals)
            body.append(f"<tr><td>{html.escape(scen)}</td>{tds}</tr>")
        zonal_html = (
            "<h2>Zonal profile</h2>"
            "<p>Land-mean of <code>median</code> by latitude band, final decade. "
            "Highest band is shaded. A global or area-weighted mean cannot see a defect "
            "confined to one zone &mdash; polar cells carry almost no area &mdash; so "
            "check the <i>shape</i> here, not just the totals.</p>"
            f'<div class="wrap"><table><tr><th>scenario</th>{head}</tr>'
            f'{"".join(body)}</table></div>')

    attrs = "".join(
        f"<tr><td><code>{html.escape(k)}</code></td><td>{html.escape(v[:400])}"
        f"{'&hellip;' if len(v) > 400 else ''}</td></tr>"
        for k, v in sorted(rep["layer_attributes"].items()))

    return f"""<title>QA report — {html.escape(rep['layer_id'])}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:2rem auto;max-width:1100px;
       padding:0 1rem;line-height:1.5;color:#222}}
 table{{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.9rem}}
 th,td{{border:1px solid #ddd;padding:.4rem .6rem;text-align:left;vertical-align:top}}
 th{{background:#f5f5f5}} code{{font-size:.85em;word-break:break-word}}
 .verdict{{display:inline-block;padding:.3rem .9rem;border-radius:4px;color:#fff;
           background:{colour};font-weight:700}}
 .wrap{{overflow-x:auto}}
</style>
<h1>QA/QC report</h1>
<p><b>{html.escape(rep['layer_id'])}</b> &middot; version
<code>{html.escape(rep['version'])}</code><br>
generated {html.escape(rep['generated_utc'])} by
<code>{html.escape(rep['generated_by'])}</code></p>
<p><span class="verdict">{s['verdict']}</span>
&nbsp; {s['passed']}/{s['checks_run']} checks passed &middot;
{s['failed']} failed &middot; {s['warnings']} warnings</p>
<h2>Checks</h2>
<div class="wrap"><table>
<tr><th>result</th><th>scope</th><th>check</th><th>evidence</th></tr>
{''.join(rows)}
</table></div>
<h2>Statistics</h2>
<div class="wrap">{''.join(stat_blocks)}</div>
{zonal_html}
<h2>Layer attributes</h2>
<div class="wrap"><table><tr><th>attribute</th><th>value</th></tr>{attrs}</table></div>
"""


def run(layer_id, variable=None, version=None, local_only=False):
    """Build and publish the QA report. Returns the report dict.

    Importable so processors can emit QA evidence as part of a normal run
    (scripts/utils/finalize.py), not only via the CLI.
    """
    warnings.filterwarnings("ignore")
    version = version or storage.resolve_current(layer_id)
    vprefix = storage.version_prefix(layer_id, version)
    storage.verify_complete(vprefix, require=["layer.json"])   # never QA ungated data

    data_dir = storage.pull_prefix(f"{vprefix}/data")
    stems = {p.stem.replace("_processed", "").rsplit("_", 1)[0]
             for p in data_dir.glob("*_processed.nc")}
    stem = variable or (stems.pop() if len(stems) == 1 else None)
    if not stem:
        raise ValueError(f"cannot infer variable from {sorted(stems)}; pass one explicitly")

    ds_by_scen = {}
    for p in sorted(data_dir.glob(f"{stem}_*_processed.nc")):
        scen = p.stem.replace("_processed", "").rsplit("_", 1)[-1]
        ds_by_scen[scen] = xr.open_dataset(p)
    if not ds_by_scen:
        raise FileNotFoundError(f"no {stem}_*_processed.nc under {data_dir}")

    log("=" * 66)
    log(f"QA/QC report: {layer_id} @ {version}")
    log("=" * 66)

    results, meta = check_layer(ds_by_scen, layer_id, version)
    rep = build_report(layer_id, version, ds_by_scen, results, meta)

    out = storage.cache_root() / "_qa" / layer_id / version
    out.mkdir(parents=True, exist_ok=True)
    (out / "qa_report.json").write_text(json.dumps(rep, indent=1))
    (out / "qa_report.html").write_text(render_html(rep))

    s = rep["summary"]
    log(f"\n{'='*66}")
    log(f"{s['verdict']}: {s['passed']}/{s['checks_run']} passed, "
        f"{s['failed']} failed, {s['warnings']} warnings")

    if local_only:
        log(f"Report written locally: {out}")
        rep["published_to"] = str(out)
    else:
        prefix = storage.publish_derived(layer_id, "qa", out, version=version)
        log(f"Published: s3://{storage.BUCKET}/{prefix}")
        rep["published_to"] = f"s3://{storage.BUCKET}/{prefix}"
    log("=" * 66)
    for ds in ds_by_scen.values():
        ds.close()
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("layer_id")
    ap.add_argument("variable", nargs="?", default=None,
                    help="File stem; inferred from the data directory if omitted")
    ap.add_argument("--version", default=None, help="Default: the layer's current version")
    ap.add_argument("--local-only", action="store_true",
                    help="Write the report to the local cache without uploading")
    args = ap.parse_args()
    try:
        rep = run(args.layer_id, args.variable, args.version, args.local_only)
    except (ValueError, FileNotFoundError) as e:
        ap.error(str(e))
    return 1 if rep["summary"]["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
