"""Generate QA report for processed qg data."""

import xarray as xr
import numpy as np
from pathlib import Path
import json

def generate_report():
    """Generate summary statistics and save as JSON report."""
    # Get project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    processed_dir = project_root / "data" / "processed"
    output_path = project_root / "reports" / "qg_qa_report.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "variable": "qg",
        "description": "Groundwater runoff (baseflow)",
        "processing_parameters": {},
        "scenarios": {}
    }

    for nc_file in sorted(processed_dir.glob("qg_*.nc")):
        print(f"Analyzing {nc_file.name}...")

        ds = xr.open_dataset(nc_file)
        scenario = ds.attrs.get("scenario", nc_file.stem.split("_")[1])

        # Store processing parameters (from first file)
        if not report["processing_parameters"]:
            report["processing_parameters"] = {
                "window_years": int(ds.attrs.get("window_years", 20)),
                "baseline_decade": int(ds.attrs.get("baseline_decade", 2020)),
                "decades": [int(d) for d in ds.decade.values],
            }

        scenario_stats = {
            "file": nc_file.name,
            "dimensions": {
                "decades": len(ds.decade),
                "lat": len(ds.lat),
                "lon": len(ds.lon),
            },
            "decades": {}
        }

        for decade in ds.decade.values:
            decade_data = ds.sel(decade=decade)
            median_vals = decade_data["median"].values
            pct_vals = decade_data["percentile"].values

            # Calculate summary stats (excluding NaN)
            valid_median = median_vals[~np.isnan(median_vals)]
            valid_pct = pct_vals[~np.isnan(pct_vals)]

            decade_stats = {
                "median": {
                    "mean": float(np.mean(valid_median)) if len(valid_median) > 0 else None,
                    "min": float(np.min(valid_median)) if len(valid_median) > 0 else None,
                    "max": float(np.max(valid_median)) if len(valid_median) > 0 else None,
                    "std": float(np.std(valid_median)) if len(valid_median) > 0 else None,
                },
                "percentile": {
                    "mean": float(np.mean(valid_pct)) if len(valid_pct) > 0 else None,
                    "p25": float(np.percentile(valid_pct, 25)) if len(valid_pct) > 0 else None,
                    "p50": float(np.percentile(valid_pct, 50)) if len(valid_pct) > 0 else None,
                    "p75": float(np.percentile(valid_pct, 75)) if len(valid_pct) > 0 else None,
                },
                "valid_cells": int(len(valid_median)),
                "total_cells": int(median_vals.size),
                "coverage_pct": float(len(valid_median) / median_vals.size * 100),
            }
            scenario_stats["decades"][str(int(decade))] = decade_stats

        report["scenarios"][scenario] = scenario_stats
        ds.close()

    # Save report
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nReport saved to: {output_path}")

    # Print summary
    print("\n" + "=" * 60)
    print("QA REPORT SUMMARY")
    print("=" * 60)
    print(f"\nVariable: {report['variable']}")
    print(f"Description: {report['description']}")
    print(f"Window size: {report['processing_parameters']['window_years']} years")
    print(f"Baseline decade: {report['processing_parameters']['baseline_decade']}s")
    print(f"Decades processed: {report['processing_parameters']['decades']}")

    for scenario, stats in report["scenarios"].items():
        print(f"\n--- {scenario.upper()} ---")
        print(f"Grid: {stats['dimensions']['lat']} x {stats['dimensions']['lon']}")

        # Show 2020s and 2080s comparison
        d2020 = stats["decades"].get("2020", {})
        d2080 = stats["decades"].get("2080", {})

        if d2020 and d2080:
            pct_2020 = d2020.get("percentile", {}).get("mean")
            pct_2080 = d2080.get("percentile", {}).get("mean")

            print(f"  2020s: percentile mean = {pct_2020:.1f}%" if pct_2020 else "  2020s: no data")
            print(f"  2080s: percentile mean = {pct_2080:.1f}%" if pct_2080 else "  2080s: no data")

            if pct_2020 and pct_2080:
                change = pct_2080 - pct_2020
                print(f"  Change: {change:+.1f} percentile points")

if __name__ == "__main__":
    generate_report()
