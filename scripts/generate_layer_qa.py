"""QA/QC report for a processed TCFD layer, read from the files themselves.

WHY THIS EXISTS. The pipeline's documented QA command, `isimip-pipeline report`, predates
OUTPUT-SPEC.md: it derives a variable name from the filename and calls `ds[that_name]`, so on
a contract file -- which carries `median`, `lower_ci`, ... and no variable of its own name --
it raises `KeyError: 'thawfrac_ssp126'`. Verified 2026-08-14; it fails the same way on every
layer built to the current contract. (`generate_qa_report.py`, hardcoded for the retired
`qg` layer, was removed 2026-08-21.)

What this writes is a MARKDOWN QA/QC record, not a picture. The three QA artifacts are
complementary and none replaces another:

    scripts/test_shared_baseline.py   the CONTRACT -- is the file shaped right
    scripts/generate_maps.py          the EYES -- six-tab dashboard, per-member contact sheet
    this script                       the RECORD -- what the layer declares about itself,
                                      the per-decade numbers, and the checks with verdicts

Every number here is READ FROM THE FILES at run time. Nothing is typed in, because a QA
figure quoted from memory is indistinguishable from an invented one, and the framing
decisions live in the NetCDF global attributes precisely so they cannot drift from what was
actually run.

A PASS here means the file is internally consistent and its declarations are present. It does
NOT mean the layer is about what its name says -- both withdrawn sugarcane layers passed every
algebraic check (GUARDRAILS 12). View the maps.

Usage:
    python scripts/generate_layer_qa.py data/processed/{layer_folder} [-o reports/qa/x.md]
"""

import argparse
import glob
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")

#: Attributes that record a load-bearing CHOICE. Absent ones are reported as gaps, because a
#: layer that does not declare how it was framed cannot be audited from the file alone.
DECLARED = [
    "decadal_statistic", "decadal_statistic_rationale", "field_nature", "normalization",
    "spatial_smoothing", "percentile_direction", "percentile_baseline", "baseline_decade",
    "baseline_source", "slope_units", "mask_rule", "value_note", "ci_definition",
    "slope_definition", "interpretation_caveat", "source_dataset",
]
VALUE_CLASSES = ["median", "lower_ci", "upper_ci", "percentile",
                 "ols_slope", "sen_slope", "n_members", "n_models"]


def fmt(x, n=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:,.{n}f}" if abs(x) < 1e4 else f"{x:.3e}"


def land_stats(a):
    f = a[np.isfinite(a)]
    if f.size == 0:
        return {}
    return dict(n=int(f.size), mean=float(f.mean()), p50=float(np.median(f)),
                p90=float(np.percentile(f, 90)), vmin=float(f.min()), vmax=float(f.max()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_dir")
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    pdir = Path(args.processed_dir)
    files = sorted(glob.glob(str(pdir / "*_processed.nc")))
    if not files:
        print(f"ERROR: no *_processed.nc in {pdir}")
        return 2
    prefix = Path(files[0]).name.rsplit("_", 2)[0]
    out = Path(args.output or f"reports/qa/{pdir.name}_qa.md")
    out.parent.mkdir(parents=True, exist_ok=True)

    L = []
    L.append(f"# QA/QC — `{pdir.name}`\n")
    L.append(f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} by "
             f"`scripts/generate_layer_qa.py` from the processed files. Every number below is "
             f"read from the NetCDF at run time.\n")

    ds0 = xr.open_dataset(files[0])
    attrs = ds0.attrs

    # ---- 1. identity + declarations ---------------------------------------- #
    L.append("## 1. What the layer declares about itself\n")
    L.append("| attribute | value |\n|---|---|")
    for k in ("variable", "source_variable", "long_name", "units", "n_members",
              "impact_models", "gcms", "soc_treatment", "co2_treatment",
              "ensemble_uniform_across_scenarios", "window_years", "minimum_models"):
        if k in attrs:
            L.append(f"| `{k}` | {str(attrs[k])[:400]} |")
    L.append("")
    missing = [k for k in DECLARED if k not in attrs]
    L.append(f"**Declared choices present**: {len(DECLARED) - len(missing)} of {len(DECLARED)}."
             + ("" if not missing else
                f" **Missing — cannot be audited from the file**: `{'`, `'.join(missing)}`."))
    L.append("")
    for k in DECLARED:
        if k in attrs and len(str(attrs[k])) > 120:
            L.append(f"**`{k}`** — {attrs[k]}\n")

    # ---- 2. structure ------------------------------------------------------- #
    L.append("## 2. Structure\n")
    decades = [int(d) for d in ds0["decade"].values]
    L.append(f"- {len(files)} scenario file(s): "
             + ", ".join(f"`{Path(f).name}`" for f in files))
    L.append(f"- dims `{dict(ds0.sizes)}`; decades {decades}")
    present = [v for v in VALUE_CLASSES if v in ds0]
    extra = [v for v in ds0.data_vars if v not in VALUE_CLASSES]
    L.append(f"- value classes present: {len(present)}/8"
             + (f"; **unexpected extras**: {extra}" if extra else ""))
    ds0.close()

    # ---- 3. per-decade table per scenario ----------------------------------- #
    L.append("\n## 3. Per-decade values\n")
    base_dec = int(attrs.get("baseline_decade", decades[0]))
    anomaly = {}
    for f in files:
        ds = xr.open_dataset(f)
        scen = ds.attrs.get("scenario", Path(f).stem)
        med = ds["median"].values
        L.append(f"\n### {scen}\n")
        L.append("| decade | cells | mean | p50 | p90 | max | CI width p50 | "
                 "ols_slope mean | sen_slope mean | Δ vs baseline |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        b = med[decades.index(base_dec)]
        for i, d in enumerate(decades):
            s = land_stats(med[i])
            if not s:
                continue
            w = ds["upper_ci"].values[i] - ds["lower_ci"].values[i]
            delta = med[i] - b
            dstat = land_stats(delta)
            anomaly.setdefault(scen, {})[d] = dstat.get("mean")
            L.append(
                f"| {d}s | {s['n']:,} | {fmt(s['mean'])} | {fmt(s['p50'])} | "
                f"{fmt(s['p90'])} | {fmt(s['vmax'])} | {fmt(np.nanmedian(w))} | "
                f"{fmt(np.nanmean(ds['ols_slope'].values[i]))} | "
                f"{fmt(np.nanmean(ds['sen_slope'].values[i]))} | {fmt(dstat.get('mean'))} |")
        ds.close()

    # ---- 4. checks ---------------------------------------------------------- #
    L.append("\n## 4. Checks\n")
    checks = []
    base_panels = []
    for f in files:
        ds = xr.open_dataset(f)
        i0 = decades.index(base_dec)
        base_panels.append(np.nan_to_num(ds["median"].values[i0], nan=-9e9))
        med = ds["median"].values
        for i, d in enumerate(decades):
            fin = np.isfinite(med[i])
            lo, hi = ds["lower_ci"].values[i], ds["upper_ci"].values[i]
            checks.append((f"{ds.attrs.get('scenario', '')} {d}s CI ordering",
                           bool(np.all(lo[fin] <= med[i][fin] + 1e-5)
                                and np.all(hi[fin] >= med[i][fin] - 1e-5)), ""))
            for k in ("ols_slope", "sen_slope"):
                sl = ds[k].values[i]
                leak = int((np.isfinite(sl) & ~fin).sum())
                if d <= base_dec:
                    checks.append((f"{ds.attrs.get('scenario', '')} {d}s {k} is NaN "
                                   "(baseline)", bool(np.all(np.isnan(sl))), ""))
                else:
                    checks.append((f"{ds.attrs.get('scenario', '')} {d}s {k} mask agrees "
                                   "with median", leak == 0,
                                   f"{leak} cells finite where median is NaN"))
            pc = ds["percentile"].values[i]
            pf = pc[np.isfinite(pc)]
            checks.append((f"{ds.attrs.get('scenario', '')} {d}s percentile in [1,100]",
                           bool(pf.size == 0 or (pf.min() >= 1 and pf.max() <= 100)), ""))
        ds.close()
    if len(base_panels) > 1:
        same = all(np.array_equal(base_panels[0], b) for b in base_panels[1:])
        checks.append(("shared baseline bit-identical across scenarios", same, ""))
    else:
        checks.append(("shared baseline across scenarios", None, "only one scenario — NOT TESTED"))

    failed = [c for c in checks if c[1] is False]
    skipped = [c for c in checks if c[1] is None]
    L.append(f"{len(checks) - len(failed) - len(skipped)} passed, **{len(failed)} failed**, "
             f"{len(skipped)} not tested.\n")
    for name, ok, note in checks:
        if ok is False:
            L.append(f"- **[FAIL]** {name} — {note}")
        elif ok is None:
            L.append(f"- [SKIP] {name} — {note}")
    if not failed and not skipped:
        L.append("- all checks passed")

    # ---- 5. contract test verbatim ------------------------------------------ #
    L.append("\n## 5. `test_shared_baseline.py` (the contract)\n")
    try:
        r = subprocess.run([sys.executable, "scripts/test_shared_baseline.py", str(pdir)],
                           capture_output=True, text=True, timeout=1800)
        L.append("```\n" + r.stdout.strip() + "\n```")
        L.append(f"\nexit code {r.returncode}"
                 + (" — **CONTRACT VIOLATION**" if r.returncode else " — passed"))
    except Exception as e:  # noqa: BLE001
        L.append(f"could not run: {type(e).__name__}: {e}")

    # ---- 6. members --------------------------------------------------------- #
    mem_path = pdir / f"{prefix}_members.nc"
    L.append("\n## 6. Members\n")
    if mem_path.exists():
        m = xr.open_dataset(mem_path)
        L.append(f"`{mem_path.name}` — {m.sizes.get('member', 0)} members. "
                 f"{m.attrs.get('member_field', '')}\n")
        L.append("| member | cells | mean | p50 | p90 | max |")
        L.append("|---|---|---|---|---|---|")
        for i, nm in enumerate(m["member"].values):
            s = land_stats(m["value"].values[i])
            if s:
                L.append(f"| {nm} | {s['n']:,} | {fmt(s['mean'])} | {fmt(s['p50'])} | "
                         f"{fmt(s['p90'])} | {fmt(s['vmax'])} |")
        if m.attrs.get("note"):
            L.append(f"\n{m.attrs['note']}")
        m.close()
    else:
        L.append(f"**No `{mem_path.name}`** — the dashboard's Members tab will be skipped "
                 "and a single-member defect cannot be seen. Every processor should emit it.")

    # ---- 7. what this does not establish ------------------------------------ #
    L.append("\n## 7. What this report does NOT establish\n")
    L.append("These are algebraic self-consistency checks. A layer can satisfy every one of "
             "them and still be geophysically wrong: both withdrawn sugarcane layers passed "
             "the full contract while being identically zero across the entire sugarcane "
             "belt. Before this layer is called verified, a human must view the maps "
             "(`scripts/generate_maps.py` output) and confirm the geography is plausible, "
             "then set `qa_reviewed_on` in `config/layer_registry.yaml`. If you have not "
             "looked at an image, the layer is *unreviewed*, not *verified*.")

    out.write_text("\n".join(L) + "\n")
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")
    print(f"{len(checks) - len(failed) - len(skipped)} passed, {len(failed)} failed, "
          f"{len(skipped)} not tested")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
