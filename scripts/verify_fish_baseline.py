"""Verify shared 2020s baseline in processed fish TCB data."""
import xarray as xr
import numpy as np
from pathlib import Path

output_dir = Path("data/processed/fish-tcb_tcb_monthly")

print("=" * 60)
print("Verifying Shared 2020s Baseline for Fish TCB Data")
print("=" * 60)

# Load both scenarios
files = sorted(output_dir.glob("tcb_*_processed.nc"))
print(f"\nFound {len(files)} processed files:")
for f in files:
    print(f"  {f.name}")

datasets = {}
for f in files:
    scenario = f.stem.split("_")[1]  # tcb_ssp126_processed -> ssp126
    datasets[scenario] = xr.open_dataset(f)
    print(f"\nLoaded {scenario}:")
    print(f"  Decades: {list(datasets[scenario].decade.values)}")
    print(f"  baseline_source: {datasets[scenario].attrs.get('baseline_source', 'N/A')}")
    print(f"  n_models: {datasets[scenario].attrs.get('n_models', 'N/A')}")
    print(f"  models: {datasets[scenario].attrs.get('models', 'N/A')}")

# Test 1: Check baseline_source attribute
print("\n" + "=" * 60)
print("TEST 1: baseline_source attribute")
print("=" * 60)
for scenario, ds in datasets.items():
    attr = ds.attrs.get('baseline_source', 'MISSING')
    status = "PASS" if attr == "shared_across_all_scenarios" else "FAIL"
    print(f"  {scenario}: {attr} [{status}]")

# Test 2: Check 2020s values are identical
print("\n" + "=" * 60)
print("TEST 2: 2020s values identical across scenarios")
print("=" * 60)

scenarios = list(datasets.keys())
if len(scenarios) >= 2:
    ds1 = datasets[scenarios[0]]
    ds2 = datasets[scenarios[1]]

    # Get 2020s median values
    median_2020s_1 = ds1['median'].sel(decade=2020).values
    median_2020s_2 = ds2['median'].sel(decade=2020).values

    # Compare (ignoring NaN)
    valid_mask = ~np.isnan(median_2020s_1) & ~np.isnan(median_2020s_2)
    diff = np.abs(median_2020s_1 - median_2020s_2)
    max_diff = np.nanmax(diff)
    mean_diff = np.nanmean(diff)

    if max_diff < 1e-6:
        print(f"  PASS: 2020s values IDENTICAL between {scenarios[0]} and {scenarios[1]}")
        print(f"        Max difference: {max_diff:.2e} g/m²")
    else:
        print(f"  FAIL: 2020s values DIFFER between scenarios")
        print(f"        Max difference: {max_diff:.2e} g/m²")
        print(f"        Mean difference: {mean_diff:.2e} g/m²")

    # Sample values
    print(f"\n  Sample 2020s values at (0, 180) [equator, prime meridian]:")
    val1 = median_2020s_1[90, 180]
    val2 = median_2020s_2[90, 180]
    print(f"    {scenarios[0]}: {val1:.3f} g/m²")
    print(f"    {scenarios[1]}: {val2:.3f} g/m²")

# Test 3: Check 2090s values differ (they should diverge)
print("\n" + "=" * 60)
print("TEST 3: 2090s values DIFFER across scenarios (expected)")
print("=" * 60)

if len(scenarios) >= 2:
    median_2090s_1 = ds1['median'].sel(decade=2090).values
    median_2090s_2 = ds2['median'].sel(decade=2090).values

    diff_2090s = np.abs(median_2090s_1 - median_2090s_2)
    max_diff = np.nanmax(diff_2090s)
    mean_diff = np.nanmean(diff_2090s[~np.isnan(diff_2090s)])

    if mean_diff > 0.01:  # Expect significant divergence
        print(f"  PASS: 2090s values DIFFER between {scenarios[0]} and {scenarios[1]}")
        print(f"        Max difference: {max_diff:.3f} g/m²")
        print(f"        Mean difference: {mean_diff:.3f} g/m²")
    else:
        print(f"  WARNING: 2090s values nearly identical (unexpected)")

    # Sample values
    print(f"\n  Sample 2090s values at (0, 180) [equator, prime meridian]:")
    val1 = median_2090s_1[90, 180]
    val2 = median_2090s_2[90, 180]
    print(f"    {scenarios[0]}: {val1:.3f} g/m²")
    print(f"    {scenarios[1]}: {val2:.3f} g/m²")
    print(f"    Change from 2020s:")
    print(f"      {scenarios[0]}: {((val1 - median_2020s_1[90, 180]) / median_2020s_1[90, 180] * 100):.1f}%")
    print(f"      {scenarios[1]}: {((val2 - median_2020s_2[90, 180]) / median_2020s_2[90, 180] * 100):.1f}%")

# Test 4: Summary statistics
print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)

for scenario, ds in datasets.items():
    print(f"\n{scenario}:")
    for decade in [2020, 2050, 2090]:
        median = ds['median'].sel(decade=decade).values
        valid = ~np.isnan(median)
        print(f"  {decade}s: mean={np.nanmean(median):.2f}, "
              f"min={np.nanmin(median):.2f}, max={np.nanmax(median):.2f} g/m², "
              f"valid cells={np.sum(valid)}")

# Close datasets
for ds in datasets.values():
    ds.close()

print("\n" + "=" * 60)
print("VERIFICATION COMPLETE")
print("=" * 60)
