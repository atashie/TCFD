"""Tests for shared 2020s baseline functionality.

Run after processing to verify:
1. 2020s values are identical across all scenarios
2. 2030s+ values differ across scenarios (as expected)
"""

import numpy as np
import xarray as xr
from pathlib import Path


def test_2020s_identical_across_scenarios():
    """Verify that 2020s values are identical across all processed scenarios."""
    # This test assumes processing has already been run
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load all scenario outputs
    scenarios = ["ssp126", "ssp370", "ssp585"]
    datasets = {}

    for scenario in scenarios:
        fpath = processed_dir / f"qg_{scenario}_processed.nc"
        if fpath.exists():
            datasets[scenario] = xr.open_dataset(fpath)

    if len(datasets) < 2:
        print("SKIP: Need at least 2 scenario outputs for comparison")
        print(f"  Found: {list(datasets.keys())}")
        return False

    # Get 2020s decade index (decade=2020)
    scenario_list = list(datasets.keys())
    first_scenario = scenario_list[0]

    all_pass = True

    # Compare 2020s values across all pairs
    for metric in ["median", "percentile", "lower_ci", "upper_ci"]:
        first_2020s = datasets[first_scenario][metric].sel(decade=2020).values

        for other_scenario in scenario_list[1:]:
            other_2020s = datasets[other_scenario][metric].sel(decade=2020).values

            # Check if values are identical (allowing for NaN)
            try:
                np.testing.assert_array_equal(
                    first_2020s,
                    other_2020s,
                    err_msg=f"{metric} differs between {first_scenario} and {other_scenario} for 2020s"
                )
                print(f"  PASS: {metric} identical for 2020s ({first_scenario} vs {other_scenario})")
            except AssertionError as e:
                print(f"  FAIL: {e}")
                all_pass = False

    # Close datasets
    for ds in datasets.values():
        ds.close()

    if all_pass:
        print("\nAll 2020s values are identical across scenarios!")
    else:
        print("\nSome 2020s values differ across scenarios (UNEXPECTED)")

    return all_pass


def test_2030s_differs_across_scenarios():
    """Verify that 2030s+ values differ across scenarios (as expected)."""
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    scenarios = ["ssp126", "ssp585"]  # Compare extreme scenarios
    datasets = {}

    for scenario in scenarios:
        fpath = processed_dir / f"qg_{scenario}_processed.nc"
        if fpath.exists():
            datasets[scenario] = xr.open_dataset(fpath)

    if len(datasets) < 2:
        print("SKIP: Need ssp126 and ssp585 outputs for comparison")
        return True  # Not a failure, just can't test

    all_pass = True

    # Compare 2080s values (should be different for extreme scenarios)
    for metric in ["median"]:
        ssp126_2080s = datasets["ssp126"][metric].sel(decade=2080).values
        ssp585_2080s = datasets["ssp585"][metric].sel(decade=2080).values

        # They should NOT be identical (assuming different scenario data)
        if np.allclose(ssp126_2080s, ssp585_2080s, equal_nan=True):
            print(f"  WARNING: {metric} is identical for 2080s (unexpected for extreme scenarios)")
            all_pass = False
        else:
            print(f"  PASS: {metric} differs for 2080s (as expected)")

    for ds in datasets.values():
        ds.close()

    return all_pass


def test_baseline_source_attribute():
    """Verify that processed files have the correct baseline_source attribute."""
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    scenarios = ["ssp126", "ssp370", "ssp585"]
    all_pass = True

    for scenario in scenarios:
        fpath = processed_dir / f"qg_{scenario}_processed.nc"
        if fpath.exists():
            ds = xr.open_dataset(fpath)
            baseline_source = ds.attrs.get("baseline_source", "NOT_SET")
            ds.close()

            if baseline_source == "shared_across_all_scenarios":
                print(f"  PASS: {scenario} has correct baseline_source attribute")
            else:
                print(f"  FAIL: {scenario} baseline_source = '{baseline_source}' (expected 'shared_across_all_scenarios')")
                all_pass = False
        else:
            print(f"  SKIP: {fpath.name} not found")

    return all_pass


if __name__ == "__main__":
    print("=" * 60)
    print("Testing shared 2020s baseline functionality")
    print("=" * 60)

    print("\n1. Testing baseline_source attribute:")
    attr_ok = test_baseline_source_attribute()

    print("\n2. Testing 2020s identical across scenarios:")
    identical_ok = test_2020s_identical_across_scenarios()

    print("\n3. Testing 2030s+ differs across scenarios:")
    differs_ok = test_2030s_differs_across_scenarios()

    print("\n" + "=" * 60)
    if attr_ok and identical_ok and differs_ok:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED - check output above")
    print("=" * 60)
