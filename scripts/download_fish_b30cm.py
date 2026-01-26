"""Download ISIMIP2b marine fisheries b30cm data for large fish analysis.

Downloads biomass density of consumers >30cm (b30cm) from ISIMIP2b marine-fishery_global sector.
This variable captures large fish species like tuna and salmon based on asymptotic length.
"""
import os
import sys
from pathlib import Path
from isimip_client.client import ISIMIPClient
import urllib.request

# Configuration
OUTPUT_DIR = Path("data/raw/fish-b30cm")
VARIABLE = "b30cm"  # Biomass density of consumers with L_infinity > 30cm
SIMULATION_ROUND = "ISIMIP2b"
SECTOR = "marine-fishery_global"

# Only download RCP projection scenarios (not historical/picontrol)
TARGET_SCENARIOS = ["rcp26", "rcp45", "rcp60", "rcp85"]

# Only include no-fishing scenarios (natural conditions)
FISHING_FILTER = "no-fishing"

def get_all_datasets(client, **kwargs):
    """Paginate through all API results."""
    all_datasets = []
    response = client.datasets(**kwargs)
    all_datasets.extend(response.get('results', []))

    while response.get('next'):
        print(f"  Fetching next page... ({len(all_datasets)} datasets so far)")
        import re
        page_match = re.search(r'page=(\d+)', response['next'])
        if page_match:
            page_num = int(page_match.group(1))
            response = client.datasets(**kwargs, page=page_num)
            all_datasets.extend(response.get('results', []))
        else:
            break

    return all_datasets

def main():
    client = ISIMIPClient()

    print(f"Searching {SIMULATION_ROUND} {SECTOR} for variable: {VARIABLE}")
    print("=" * 60)

    # Get all datasets
    all_datasets = get_all_datasets(
        client,
        simulation_round=SIMULATION_ROUND,
        sector=SECTOR,
        variable=VARIABLE
    )

    print(f"\nFound {len(all_datasets)} total {VARIABLE} datasets")

    # Filter for RCP scenarios only (projection data) and no-fishing
    rcp_datasets = [
        ds for ds in all_datasets
        if ds.get('specifiers', {}).get('climate_scenario', '') in TARGET_SCENARIOS
        and ds.get('specifiers', {}).get('fishing_type', '') == FISHING_FILTER
    ]

    print(f"Filtered to {len(rcp_datasets)} RCP projection datasets (no-fishing only)")

    # Show what we'll download
    print("\n=== Datasets to Download ===")
    total_size = 0
    files_to_download = []

    for ds in rcp_datasets:
        spec = ds.get('specifiers', {})
        model = spec.get('model', 'unknown')
        gcm = spec.get('climate_forcing', 'unknown')
        scenario = spec.get('climate_scenario', 'unknown')

        for f in ds.get('files', []):
            fname = f.get('name', '')
            fsize = f.get('size', 0)
            furl = f.get('file_url', '')

            total_size += fsize
            files_to_download.append({
                'name': fname,
                'url': furl,
                'size': fsize,
                'model': model,
                'gcm': gcm,
                'scenario': scenario
            })
            print(f"  {model} / {gcm} / {scenario}")
            print(f"    {fname} ({fsize / 1e6:.1f} MB)")

    print(f"\nTotal: {len(files_to_download)} files, {total_size / 1e9:.2f} GB")

    # Confirm download
    if "--yes" not in sys.argv:
        response = input("\nProceed with download? [y/N]: ")
        if response.lower() != 'y':
            print("Download cancelled.")
            return

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Download files
    print("\n=== Downloading ===")
    for i, f in enumerate(files_to_download, 1):
        dest = OUTPUT_DIR / f['name']

        if dest.exists():
            print(f"[{i}/{len(files_to_download)}] Skipping (exists): {f['name']}")
            continue

        print(f"[{i}/{len(files_to_download)}] Downloading: {f['name']} ({f['size'] / 1e6:.1f} MB)")

        try:
            urllib.request.urlretrieve(f['url'], dest)
            print(f"    Saved to: {dest}")
        except Exception as e:
            print(f"    ERROR: {e}")

    print("\n=== Download Complete ===")
    print(f"Files saved to: {OUTPUT_DIR.absolute()}")

if __name__ == "__main__":
    main()
