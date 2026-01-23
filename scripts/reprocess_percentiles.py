"""Reprocess percentile values in processed NetCDF files.

Inverts percentile values for "higher is better" variables so that:
- Percentile 1-20 = Very Low risk (safe/good)
- Percentile 81-100 = Very High risk (dangerous/bad)

For "higher is better" variables (biomass, carbon, productivity), high raw values
should map to LOW percentiles (good), but the original processing did not account
for this. This script fixes that by inverting: percentile = 100 - percentile.

Usage:
    python scripts/reprocess_percentiles.py
"""

from pathlib import Path

import xarray as xr


# Directories with "higher is better" variables (need percentile inversion)
HIGHER_IS_BETTER_DIRS = [
    "evgndltr_npp-cveg",
    "fish-b30cm_b30cm_monthly",
    "loblolly-pine-proxy_cwood-evgndltr_annual",
    "loblolly-temperate_cveg-needleleaf-evergreen-tree-temperate_annual",
    "loblolly-temperate_npp-needleleaf-evergreen-tree-temperate_annual",
    "oak-timber_cwood_annual",
    "tebrsu_cveg_annual",
    "tebrsu_npp_annual",
]

# Directories with "higher is worse" variables (keep percentile as-is)
HIGHER_IS_WORSE_DIRS = [
    "health-mortality_an-tot-heat_annual",
]


def reprocess_file(nc_file: Path, direction: str) -> None:
    """Reprocess a single NetCDF file.

    Args:
        nc_file: Path to the processed NetCDF file
        direction: "higher_is_better" or "higher_is_worse"
    """
    # Load dataset into memory and close file immediately
    with xr.open_dataset(nc_file) as ds:
        # Check if already processed
        if ds.attrs.get("percentile_direction") == direction:
            print(f"  Skipping (already processed): {nc_file.name}")
            return

        # Load all data into memory
        ds_memory = ds.load()

    # Invert percentile for "higher is better" variables
    if direction == "higher_is_better" and "percentile" in ds_memory:
        old_min = float(ds_memory["percentile"].min())
        old_max = float(ds_memory["percentile"].max())

        ds_memory["percentile"] = 100 - ds_memory["percentile"]

        new_min = float(ds_memory["percentile"].min())
        new_max = float(ds_memory["percentile"].max())

        print(f"  Inverted percentile: [{old_min:.1f}-{old_max:.1f}] -> [{new_min:.1f}-{new_max:.1f}]")

    # Add percentile_direction attribute
    ds_memory.attrs["percentile_direction"] = direction

    # Save back to the same file (file is now closed)
    ds_memory.to_netcdf(nc_file)

    print(f"  Saved: {nc_file.name}")


def reprocess_directory(dir_path: Path, direction: str) -> int:
    """Reprocess all NetCDF files in a directory.

    Args:
        dir_path: Path to the processed data directory
        direction: "higher_is_better" or "higher_is_worse"

    Returns:
        Number of files processed
    """
    if not dir_path.exists():
        print(f"  Directory not found: {dir_path}")
        return 0

    nc_files = list(dir_path.glob("*_processed.nc"))
    if not nc_files:
        print(f"  No processed files found in: {dir_path}")
        return 0

    print(f"\nProcessing {dir_path.name} ({direction}):")
    for nc_file in sorted(nc_files):
        reprocess_file(nc_file, direction)

    return len(nc_files)


def main():
    """Main entry point."""
    processed_dir = Path("data/processed")

    if not processed_dir.exists():
        print(f"Error: {processed_dir} does not exist")
        return

    total_files = 0

    # Process "higher is better" directories (invert percentile)
    print("=" * 60)
    print("Processing 'higher is better' datasets (inverting percentile)")
    print("=" * 60)

    for dir_name in HIGHER_IS_BETTER_DIRS:
        dir_path = processed_dir / dir_name
        total_files += reprocess_directory(dir_path, "higher_is_better")

    # Process "higher is worse" directories (keep percentile, add attribute)
    print("\n" + "=" * 60)
    print("Processing 'higher is worse' datasets (adding attribute only)")
    print("=" * 60)

    for dir_name in HIGHER_IS_WORSE_DIRS:
        dir_path = processed_dir / dir_name
        total_files += reprocess_directory(dir_path, "higher_is_worse")

    print("\n" + "=" * 60)
    print(f"Complete! Processed {total_files} files.")
    print("=" * 60)


if __name__ == "__main__":
    main()
