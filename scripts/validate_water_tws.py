"""Statistical QA/QC validation for processed water index files.

Runs dimensional, range, consistency, and spot-check validations.
Supports any water variable via --variable flag.

Usage:
    python scripts/validate_water_tws.py                          # Validate TWS (default)
    python scripts/validate_water_tws.py --variable rootmoist     # Validate rootmoist
    python scripts/validate_water_tws.py --file /path/to/file.nc  # Custom file path
    python scripts/validate_water_tws.py --json-report report.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import xarray as xr

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_water_variables import get_variable_config, SHARED_CONFIG, VALUE_TYPE_NAMES
from utils.water_index_compare import (
    summary_statistics, flag_suspicious_patterns, compare_seasonal_cycle,
)

# Rich is optional — fall back to plain print
try:
    from rich.console import Console
    from rich.table import Table
    console = Console()
    HAS_RICH = True
except ImportError:
    HAS_RICH = False
    class Console:
        def print(self, *a, **kw):
            print(*a)
    console = Console()


# Reference locations for spot checks (name, lat, lon)
REFERENCE_LOCATIONS = [
    ("Amazon Basin", -3.0, -60.0),
    ("Sahara Desert", 23.0, 5.0),
    ("Central US", 40.0, -95.0),
    ("Ganges Delta", 23.0, 89.0),
    ("Greenland", 72.0, -40.0),
]


def log(msg: str, style: str = ""):
    if HAS_RICH and style:
        console.print(msg, style=style)
    else:
        print(msg, flush=True)


def check_dimensions(ds: xr.Dataset, variable: str, scenarios: list,
                     decades: list, n_value_types: int, n_lat: int, n_lon: int) -> List[Dict]:
    """Check that dimensions match expected shape."""
    issues = []
    expected = {
        "lat": n_lat,
        "lon": n_lon,
        "scenario": len(scenarios),
        "value_type": n_value_types,
        "decade": len(decades),
    }

    for dim, expected_size in expected.items():
        if dim not in ds.dims:
            issues.append({"check": "dimensions", "severity": "ERROR",
                           "message": f"Missing dimension: {dim}"})
        elif ds.dims[dim] != expected_size:
            issues.append({"check": "dimensions", "severity": "WARNING",
                           "message": f"{dim}: got {ds.dims[dim]}, expected {expected_size}"})

    if variable not in ds.data_vars:
        issues.append({"check": "dimensions", "severity": "ERROR",
                       "message": f"Missing variable: {variable}"})

    return issues


def check_coordinates(ds: xr.Dataset, scenarios: list, decades: list) -> List[Dict]:
    """Check coordinate values."""
    issues = []

    # Lat should go from ~89.75 to ~-89.75
    if "lat" in ds.coords:
        lats = ds.lat.values
        if lats[0] < lats[-1]:
            issues.append({"check": "coordinates", "severity": "WARNING",
                           "message": "Latitudes are ascending (expected descending N->S)"})
        if abs(lats[0] - 89.75) > 0.1:
            issues.append({"check": "coordinates", "severity": "WARNING",
                           "message": f"First lat is {lats[0]}, expected ~89.75"})

    # Lon should go from ~-179.75 to ~179.75
    if "lon" in ds.coords:
        lons = ds.lon.values
        if abs(lons[0] - (-179.75)) > 0.1:
            issues.append({"check": "coordinates", "severity": "WARNING",
                           "message": f"First lon is {lons[0]}, expected ~-179.75"})

    # Scenarios
    if "scenario" in ds.coords:
        actual_scenarios = list(ds.scenario.values)
        if actual_scenarios != scenarios:
            issues.append({"check": "coordinates", "severity": "INFO",
                           "message": f"Scenarios: {actual_scenarios} (expected {scenarios})"})

    # Decades
    if "decade" in ds.coords:
        actual_decades = list(ds.decade.values)
        if actual_decades != decades:
            issues.append({"check": "coordinates", "severity": "INFO",
                           "message": f"Decades: {actual_decades} (expected {decades})"})

    return issues


def check_nan_patterns(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Check NaN patterns for ocean mask consistency."""
    issues = []

    # Use vt 12 (Annual_Mean) as reference
    ref_slice = data.isel(scenario=0, value_type=12, decade=-1).values
    ref_mask = np.isnan(ref_slice)
    total_land = (~ref_mask).sum()
    total_cells = ref_mask.size

    issues.append({"check": "nan_patterns", "severity": "INFO",
                   "message": f"Land cells: {total_land:,}/{total_cells:,} "
                              f"({total_land/total_cells*100:.1f}%)"})

    # Compare NaN mask across scenarios for same decade
    for d_idx in range(len(decades)):
        masks = []
        for s_idx in range(len(scenarios)):
            m = np.isnan(data.isel(scenario=s_idx, value_type=12, decade=d_idx).values)
            masks.append(m)
        for s_idx in range(1, len(scenarios)):
            diff = (masks[0] != masks[s_idx]).sum()
            if diff > 0:
                issues.append({"check": "nan_patterns", "severity": "WARNING",
                               "message": f"Decade {decades[d_idx]}: NaN mask differs by "
                                          f"{diff} cells between {scenarios[0]} and {scenarios[s_idx]}"})

    return issues


def check_value_ranges(data: xr.DataArray, n_value_types: int) -> List[Dict]:
    """Check that each value_type has plausible ranges."""
    issues = []

    for vt in range(n_value_types):
        name = VALUE_TYPE_NAMES.get(vt, f"vt{vt}")
        vals = data.isel(value_type=vt).values
        valid = vals[~np.isnan(vals)]
        if len(valid) == 0:
            issues.append({"check": "value_ranges", "severity": "WARNING",
                           "message": f"vt {vt} ({name}): ALL NaN"})
            continue

        actual_min, actual_max = float(np.min(valid)), float(np.max(valid))
        issues.append({"check": "value_ranges", "severity": "OK",
                       "message": f"vt {vt} ({name}): [{actual_min:.4g}, {actual_max:.4g}]"})

    return issues


def check_annual_mean_consistency(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Check that vt12 (Annual_Mean) equals mean of vt0-11 (monthly means)."""
    issues = []

    for s_idx, scenario in enumerate(scenarios):
        for d_idx in [0, len(decades) // 2, -1]:
            monthly = data.isel(scenario=s_idx, decade=d_idx, value_type=slice(0, 12)).values
            annual_mean = data.isel(scenario=s_idx, decade=d_idx, value_type=12).values
            computed_mean = np.nanmean(monthly, axis=2)  # axis 2 is value_type

            mask = ~np.isnan(annual_mean) & ~np.isnan(computed_mean)
            if mask.sum() == 0:
                continue

            max_diff = float(np.max(np.abs(annual_mean[mask] - computed_mean[mask])))
            decade = decades[d_idx]
            severity = "OK" if max_diff < 1e-3 else "WARNING"
            issues.append({
                "check": "annual_mean", "severity": severity,
                "message": f"{scenario}/{decade}s: vt12 vs mean(vt0-11) max_diff={max_diff:.6g}"
            })

    return issues


def check_quantile_ordering(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Check that quantiles are monotonically non-decreasing (Q05 <= Q15 <= ... <= Q95)."""
    issues = []
    quantile_vts = list(range(13, 20))  # vt13=Q05 through vt19=Q95

    for s_idx, scenario in enumerate(scenarios):
        for d_idx in [0, len(decades) // 2, -1]:
            decade = decades[d_idx]
            violations = 0
            total_cells = 0

            for q_i in range(len(quantile_vts) - 1):
                lower = data.isel(scenario=s_idx, decade=d_idx, value_type=quantile_vts[q_i]).values
                upper = data.isel(scenario=s_idx, decade=d_idx, value_type=quantile_vts[q_i + 1]).values
                mask = ~np.isnan(lower) & ~np.isnan(upper)
                total_cells += mask.sum()
                violations += (lower[mask] > upper[mask] + 1e-6).sum()

            if total_cells > 0:
                pct = violations / total_cells * 100
                severity = "OK" if pct < 0.01 else ("WARNING" if pct < 1 else "ERROR")
                issues.append({
                    "check": "quantile_ordering", "severity": severity,
                    "message": f"{scenario}/{decade}s: {violations}/{total_cells} "
                               f"quantile ordering violations ({pct:.4f}%)"
                })

    return issues


def check_cross_scenario_consistency(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Check that 2010s-2020s values are similar across scenarios."""
    issues = []

    for d_idx in [0, 1]:
        if d_idx >= len(decades):
            continue
        decade = decades[d_idx]

        # Compare vt 12 (Annual_Mean) across scenarios
        slices = []
        for s_idx in range(len(scenarios)):
            slices.append(data.isel(scenario=s_idx, value_type=12, decade=d_idx).values)

        for i in range(1, len(scenarios)):
            mask = ~np.isnan(slices[0]) & ~np.isnan(slices[i])
            if mask.sum() < 100:
                continue
            corr = np.corrcoef(slices[0][mask], slices[i][mask])[0, 1]
            diff_pct = np.mean(np.abs(slices[0][mask] - slices[i][mask])) / (np.mean(np.abs(slices[0][mask])) + 1e-10) * 100

            severity = "OK" if corr > 0.99 else ("WARNING" if corr > 0.95 else "ERROR")
            issues.append({
                "check": "cross_scenario", "severity": severity,
                "message": f"Decade {decade}: {scenarios[0]} vs {scenarios[i]} "
                           f"correlation={corr:.4f}, mean_diff={diff_pct:.2f}%"
            })

    return issues


def check_seasonal_sanity(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Check that monthly values (vt 0-11) show realistic seasonal patterns."""
    issues = []

    for name, lat, lon in REFERENCE_LOCATIONS:
        for s_idx, scenario in enumerate(scenarios[:1]):  # Just first scenario
            lat_idx = int(np.argmin(np.abs(data.lat.values - lat)))
            lon_idx = int(np.argmin(np.abs(data.lon.values - lon)))
            monthly = np.array([
                float(data.isel(lat=lat_idx, lon=lon_idx).sel(
                    scenario=scenario, value_type=m, decade=decades[-1]).values)
                for m in range(12)
            ])

            if np.all(np.isnan(monthly)):
                issues.append({"check": "seasonal", "severity": "INFO",
                               "message": f"{name} ({scenario}): all NaN (likely ocean)"})
                continue

            valid_months = (~np.isnan(monthly)).sum()
            if valid_months < 6:
                issues.append({"check": "seasonal", "severity": "WARNING",
                               "message": f"{name} ({scenario}): only {valid_months}/12 months have data"})
                continue

            cv = np.nanstd(monthly) / (np.abs(np.nanmean(monthly)) + 1e-10)
            if cv < 0.001:
                issues.append({"check": "seasonal", "severity": "WARNING",
                               "message": f"{name} ({scenario}): no seasonality (CV={cv:.4f})"})
            else:
                peak_month = np.nanargmax(monthly) + 1
                trough_month = np.nanargmin(monthly) + 1
                issues.append({"check": "seasonal", "severity": "OK",
                               "message": f"{name} ({scenario}): peak month={peak_month}, "
                                          f"trough={trough_month}, CV={cv:.3f}"})

    return issues


def check_spot_locations(data: xr.DataArray, scenarios: list, decades: list) -> List[Dict]:
    """Spot-check reference locations for physically plausible values."""
    issues = []

    for name, lat, lon in REFERENCE_LOCATIONS:
        for s_idx, scenario in enumerate(scenarios[:1]):  # Just first scenario
            lat_idx = int(np.argmin(np.abs(data.lat.values - lat)))
            lon_idx = int(np.argmin(np.abs(data.lon.values - lon)))
            cell = data.isel(lat=lat_idx, lon=lon_idx).sel(scenario=scenario, decade=decades[-1])
            vt12 = float(cell.sel(value_type=12).values)
            vt16 = float(cell.sel(value_type=16).values)  # Q50 (median)

            if np.isnan(vt12):
                issues.append({"check": "spot_check", "severity": "INFO",
                               "message": f"{name}: NaN (likely ocean/ice)"})
            else:
                issues.append({"check": "spot_check", "severity": "OK",
                               "message": f"{name}: Annual_Mean={vt12:.4g}, Q50={vt16:.4g}"})

    return issues


def run_all_checks(file_path: Path, variable: str, scenarios: list,
                   decades: list, n_value_types: int, n_lat: int, n_lon: int) -> Dict:
    """Run all validation checks and return results."""
    log(f"\nValidating: {file_path}")
    log(f"Variable: {variable}")
    log("=" * 70)

    ds = xr.open_dataset(file_path)
    data = ds[variable] if variable in ds.data_vars else ds[list(ds.data_vars)[0]]

    all_issues = []
    checks = [
        ("Dimensions", lambda arg: check_dimensions(arg, variable, scenarios, decades, n_value_types, n_lat, n_lon), ds),
        ("Coordinates", lambda arg: check_coordinates(arg, scenarios, decades), ds),
        ("NaN Patterns", lambda arg: check_nan_patterns(arg, scenarios, decades), data),
        ("Value Ranges", lambda arg: check_value_ranges(arg, n_value_types), data),
        ("Annual Mean Consistency", lambda arg: check_annual_mean_consistency(arg, scenarios, decades), data),
        ("Quantile Ordering", lambda arg: check_quantile_ordering(arg, scenarios, decades), data),
        ("Cross-Scenario Consistency", lambda arg: check_cross_scenario_consistency(arg, scenarios, decades), data),
        ("Seasonal Sanity", lambda arg: check_seasonal_sanity(arg, scenarios, decades), data),
        ("Spot Checks", lambda arg: check_spot_locations(arg, scenarios, decades), data),
    ]

    for check_name, check_fn, check_arg in checks:
        log(f"\n--- {check_name} ---")
        issues = check_fn(check_arg)
        all_issues.extend(issues)

        for issue in issues:
            sev = issue["severity"]
            msg = issue["message"]
            if sev == "ERROR":
                log(f"  [ERROR] {msg}", "bold red")
            elif sev == "WARNING":
                log(f"  [WARN]  {msg}", "yellow")
            elif sev == "OK":
                log(f"  [OK]    {msg}", "green")
            else:
                log(f"  [INFO]  {msg}")

    ds.close()

    # Summary
    errors = sum(1 for i in all_issues if i["severity"] == "ERROR")
    warnings = sum(1 for i in all_issues if i["severity"] == "WARNING")
    oks = sum(1 for i in all_issues if i["severity"] == "OK")

    log(f"\n{'='*70}")
    log(f"Summary: {errors} errors, {warnings} warnings, {oks} passed")

    if errors == 0 and warnings == 0:
        log("RESULT: ALL CHECKS PASSED", "bold green")
    elif errors == 0:
        log(f"RESULT: PASSED WITH {warnings} WARNING(S)", "bold yellow")
    else:
        log(f"RESULT: FAILED WITH {errors} ERROR(S)", "bold red")

    return {
        "file": str(file_path),
        "variable": variable,
        "summary": {"errors": errors, "warnings": warnings, "ok": oks},
        "issues": all_issues,
    }


def main():
    parser = argparse.ArgumentParser(description="Validate water index NetCDF file")
    parser.add_argument(
        "--variable", type=str, default="tws",
        help="Variable to validate (default: tws). Options: tws, rootmoist, dis, qr, potevap, precip",
    )
    parser.add_argument(
        "--file", type=Path, default=None,
        help="Path to NetCDF file to validate (auto-derived from variable if not specified)",
    )
    parser.add_argument(
        "--json-report", type=Path, default=None,
        help="Save validation report as JSON",
    )
    args = parser.parse_args()

    # Derive constants from config_water_variables
    var_config = get_variable_config(args.variable)
    variable = var_config.name
    scenarios = list(SHARED_CONFIG.scenarios)
    decades = list(SHARED_CONFIG.decades)
    n_value_types = SHARED_CONFIG.n_value_types
    n_lat = SHARED_CONFIG.n_lat
    n_lon = SHARED_CONFIG.n_lon

    file_path = args.file or Path(rf"C:\Cai_data\WaterIndex\{var_config.output_filename}")

    if not file_path.exists():
        log(f"File not found: {file_path}")
        sys.exit(1)

    report = run_all_checks(file_path, variable, scenarios, decades, n_value_types, n_lat, n_lon)

    if args.json_report:
        def _convert(obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj

        with open(args.json_report, "w") as f:
            json.dump(report, f, indent=2, default=_convert)
        log(f"\nJSON report saved: {args.json_report}")

    sys.exit(1 if report["summary"]["errors"] > 0 else 0)


if __name__ == "__main__":
    main()
