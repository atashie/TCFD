"""Download qr (groundwater recharge) monthly data from ISIMIP3b for water index processing.

Queries the ISIMIP API for all qr monthly files across impact models, GCMs,
and SSP scenarios. Reports availability matrix, downloads with resume support,
and validates each file.

Usage:
    python scripts/download_water_qr.py                  # Check availability only
    python scripts/download_water_qr.py --download        # Download all files
    python scripts/download_water_qr.py --download --dry-run  # Show what would download
    python scripts/download_water_qr.py --validate        # Check existing files
"""

import argparse
import asyncio
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "isimip-pipeline" / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_water_variables import get_variable_config, SHARED_CONFIG
from isimip_pipeline.search.isimip_query import ISIMIPQuery, SearchFilters, DatasetInfo
from isimip_pipeline.download.downloader import Downloader

# Load variable config
cfg = get_variable_config("qr")
VARIABLE = cfg.name
SECTOR = cfg.sector
SIMULATION_ROUND = cfg.simulation_round
TIMESTEP = cfg.timestep
SOCIAL_FORCING = cfg.social_forcing
IMPACT_MODELS = cfg.models
GCMS = list(SHARED_CONFIG.gcms)
SCENARIOS = list(SHARED_CONFIG.scenarios)
DOWNLOAD_BASE_DIR = cfg.download_subdir


def get_download_subdir(model: str, gcm: str, scenario: str) -> str:
    """Return subdirectory path for organizing downloads."""
    return f"{model}/{gcm}_{scenario}"


def log(msg: str):
    print(msg, flush=True)


def filter_by_social_forcing(datasets: List[DatasetInfo]) -> List[DatasetInfo]:
    """Filter datasets to prefer 2015soc, fall back to 2015soc-from-histsoc.

    For each model/GCM/scenario combination:
    - Use '2015soc' if available (standard constant social forcing)
    - Fall back to '2015soc-from-histsoc' (e.g., CWatM only has this variant)
    - Skip '1850soc' (pre-industrial counterfactual, not relevant)
    """
    # Group by model+GCM+scenario key
    groups: Dict[str, List[DatasetInfo]] = {}
    for ds in datasets:
        # Extract key from filename: model_gcm_*_scenario
        name = ds.name
        key = f"{ds.model}_{ds.climate_scenario}"
        # Also include GCM in key
        for gcm in GCMS:
            if gcm in name:
                key = f"{ds.model}_{gcm}_{ds.climate_scenario}"
                break
        if key not in groups:
            groups[key] = []
        groups[key].append(ds)

    filtered = []
    for key, group in groups.items():
        # Prefer 2015soc, fall back to 2015soc-from-histsoc
        soc_2015 = [d for d in group if "_2015soc_" in d.name]
        soc_from_hist = [d for d in group if "2015soc-from-histsoc" in d.name]

        if soc_2015:
            filtered.extend(soc_2015)
        elif soc_from_hist:
            filtered.extend(soc_from_hist)
        else:
            # No 2015soc variant found — skip (1850soc only)
            log(f"  WARNING: {key} has no 2015soc variant, skipping")

    return filtered


def search_all_files(query: ISIMIPQuery) -> List[DatasetInfo]:
    """Search ISIMIP API for all qr monthly files across all combinations."""
    all_datasets = []

    for model in IMPACT_MODELS:
        for scenario in SCENARIOS:
            filters = SearchFilters(
                simulation_round=SIMULATION_ROUND,
                variable=VARIABLE,
                model=model,
                climate_scenario=scenario,
                timestep=TIMESTEP,
            )
            results = query.search(filters)

            for ds in results:
                all_datasets.append(ds)

            if results:
                log(f"  {model} / {scenario}: {len(results)} datasets (unfiltered)")
            else:
                log(f"  {model} / {scenario}: not found")

    # Filter to one social forcing variant per combination
    log(f"\n  Total before social forcing filter: {len(all_datasets)}")
    filtered = filter_by_social_forcing(all_datasets)
    log(f"  Total after filter (2015soc preferred, 2015soc-from-histsoc fallback): {len(filtered)}")

    return filtered


def search_exploratory(query: ISIMIPQuery):
    """Run exploratory searches to discover all available models."""
    log("\n--- Exploratory Search: ISIMIP3b qr (no model filter) ---")
    filters = SearchFilters(
        simulation_round="ISIMIP3b",
        variable=VARIABLE,
        timestep=TIMESTEP,
    )
    results = query.search(filters)
    if results:
        models_found = set()
        for ds in results:
            models_found.add(ds.model or "unknown")
        log(f"  Found {len(results)} datasets across models: {sorted(models_found)}")
    else:
        log("  No ISIMIP3b qr data found")

    log("\n--- Exploratory Search: ISIMIP2b qr (old RCP ensemble) ---")
    filters_2b = SearchFilters(
        simulation_round="ISIMIP2b",
        variable=VARIABLE,
        timestep="monthly",
    )
    results_2b = query.search(filters_2b)
    if results_2b:
        models_2b = set()
        for ds in results_2b:
            models_2b.add(ds.model or "unknown")
        log(f"  Found {len(results_2b)} datasets across models: {sorted(models_2b)}")
    else:
        log("  No ISIMIP2b qr data found")


def build_availability_matrix(
    datasets: List[DatasetInfo],
) -> Dict[str, Dict[str, Dict[str, List[DatasetInfo]]]]:
    """Organize datasets into model -> gcm -> scenario structure."""
    matrix: Dict[str, Dict[str, Dict[str, List[DatasetInfo]]]] = {}

    for ds in datasets:
        model = ds.model or "unknown"
        scenario = ds.climate_scenario or "unknown"

        gcm = "unknown"
        for known_gcm in GCMS:
            if known_gcm in ds.name:
                gcm = known_gcm
                break

        if model not in matrix:
            matrix[model] = {}
        if gcm not in matrix[model]:
            matrix[model][gcm] = {}
        if scenario not in matrix[model][gcm]:
            matrix[model][gcm][scenario] = []
        matrix[model][gcm][scenario].append(ds)

    return matrix


def print_availability_report(
    matrix: Dict[str, Dict[str, Dict[str, List[DatasetInfo]]]],
):
    """Print a formatted availability table."""
    log("\n" + "=" * 80)
    log("ISIMIP3b Groundwater Recharge (qr) Monthly Data Availability")
    log("=" * 80)

    total_files = 0
    total_size_mb = 0

    for model in sorted(matrix.keys()):
        log(f"\n  Model: {model}")
        log(f"  {'GCM':<20} {'ssp126':>8} {'ssp370':>8} {'ssp585':>8}")
        log(f"  {'-'*20} {'-'*8} {'-'*8} {'-'*8}")

        for gcm in sorted(matrix[model].keys()):
            counts = []
            for scenario in SCENARIOS:
                files = matrix[model].get(gcm, {}).get(scenario, [])
                n = len(files)
                total_files += n
                for f in files:
                    if f.size:
                        total_size_mb += f.size / (1024 * 1024)
                counts.append(f"{n:>8}")
            log(f"  {gcm:<20} {''.join(counts)}")

    log(f"\n  Total files: {total_files}")
    log(f"  Estimated total size: {total_size_mb:.0f} MB ({total_size_mb/1024:.1f} GB)")
    log("=" * 80)


def collect_download_urls(
    matrix: Dict[str, Dict[str, Dict[str, List[DatasetInfo]]]],
) -> List[Tuple[str, str, DatasetInfo]]:
    """Collect all (output_subdir, url, dataset) tuples for download."""
    downloads = []
    for model in matrix:
        for gcm in matrix[model]:
            for scenario in matrix[model][gcm]:
                for ds in matrix[model][gcm][scenario]:
                    subdir = get_download_subdir(model, gcm, scenario)
                    downloads.append((subdir, ds.url, ds))
    return downloads


async def download_with_organization(
    downloads: List[Tuple[str, str, DatasetInfo]],
    base_dir: Path,
    max_concurrent: int = 4,
    dry_run: bool = False,
) -> Dict[str, int]:
    """Download files organized into subdirectories.

    Returns dict with counts: downloaded, skipped, failed.
    """
    import httpx

    stats = {"downloaded": 0, "skipped": 0, "failed": 0}

    if dry_run:
        for subdir, url, ds in downloads:
            dest = base_dir / subdir / ds.name
            exists = dest.exists()
            size_str = f"{ds.size / 1024 / 1024:.1f} MB" if ds.size else "? MB"
            status = "SKIP (exists)" if exists else "DOWNLOAD"
            log(f"  [{status}] {dest} ({size_str})")
            if exists:
                stats["skipped"] += 1
            else:
                stats["downloaded"] += 1
        return stats

    semaphore = asyncio.Semaphore(max_concurrent)

    async def download_one(
        client: httpx.AsyncClient,
        subdir: str,
        url: str,
        ds: DatasetInfo,
    ):
        dest_dir = base_dir / subdir
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / ds.name

        if dest.exists():
            log(f"  SKIP (exists): {ds.name}")
            stats["skipped"] += 1
            return

        async with semaphore:
            try:
                log(f"  Downloading: {ds.name}...")
                async with client.stream("GET", url) as response:
                    response.raise_for_status()
                    with open(dest, "wb") as f:
                        async for chunk in response.aiter_bytes(8192):
                            f.write(chunk)
                size_mb = dest.stat().st_size / (1024 * 1024)
                log(f"  OK: {ds.name} ({size_mb:.1f} MB)")
                stats["downloaded"] += 1
            except Exception as e:
                log(f"  FAILED: {ds.name} - {e}")
                stats["failed"] += 1
                if dest.exists():
                    dest.unlink()

    async with httpx.AsyncClient(timeout=600) as client:
        tasks = [
            download_one(client, subdir, url, ds)
            for subdir, url, ds in downloads
        ]
        await asyncio.gather(*tasks)

    return stats


def validate_downloads(base_dir: Path) -> List[str]:
    """Validate downloaded NetCDF files: check dimensions, time range, variable."""
    import xarray as xr

    issues = []
    nc_files = sorted(base_dir.rglob("*.nc"))

    for fpath in nc_files:
        try:
            with xr.open_dataset(fpath) as ds:
                if VARIABLE not in ds.data_vars:
                    issues.append(f"{fpath.name}: variable '{VARIABLE}' not found (has: {list(ds.data_vars)})")
                    continue

                if "time" not in ds.dims:
                    issues.append(f"{fpath.name}: missing 'time' dimension")
                if "lat" not in ds.dims:
                    issues.append(f"{fpath.name}: missing 'lat' dimension")
                if "lon" not in ds.dims:
                    issues.append(f"{fpath.name}: missing 'lon' dimension")

                if "lat" in ds.dims and len(ds.lat) != 360:
                    issues.append(f"{fpath.name}: unexpected lat size {len(ds.lat)} (expected 360)")
                if "lon" in ds.dims and len(ds.lon) != 720:
                    issues.append(f"{fpath.name}: unexpected lon size {len(ds.lon)} (expected 720)")

        except Exception as e:
            issues.append(f"{fpath.name}: failed to open - {e}")

    return issues


def save_manifest(
    matrix: Dict[str, Dict[str, Dict[str, List[DatasetInfo]]]],
    output_path: Path,
):
    """Save download manifest as JSON."""
    manifest = {
        "variable": VARIABLE,
        "simulation_round": SIMULATION_ROUND,
        "timestep": TIMESTEP,
        "models": {},
    }

    for model in matrix:
        manifest["models"][model] = {}
        for gcm in matrix[model]:
            manifest["models"][model][gcm] = {}
            for scenario in matrix[model][gcm]:
                files = [
                    {
                        "name": ds.name,
                        "url": ds.url,
                        "size": ds.size,
                    }
                    for ds in matrix[model][gcm][scenario]
                ]
                manifest["models"][model][gcm][scenario] = files

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(manifest, f, indent=2)
    log(f"\nManifest saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Download qr (groundwater recharge) monthly data from ISIMIP3b")
    parser.add_argument("--download", action="store_true", help="Actually download files")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be downloaded")
    parser.add_argument("--validate", action="store_true", help="Validate existing downloads")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help=f"Output directory (default: {DOWNLOAD_BASE_DIR})")
    args = parser.parse_args()

    base_dir = args.output_dir or (PROJECT_ROOT / DOWNLOAD_BASE_DIR)

    log("=" * 60)
    log("Groundwater Recharge (qr) Monthly Data Download (ISIMIP3b)")
    log("=" * 60)
    log(f"Variable: {VARIABLE}")
    log(f"Models: {', '.join(IMPACT_MODELS)}")
    log(f"GCMs: {', '.join(GCMS)}")
    log(f"Scenarios: {', '.join(SCENARIOS)}")
    log(f"Timestep: {TIMESTEP}")

    if args.validate:
        log(f"\nValidating downloads in {base_dir}...")
        issues = validate_downloads(base_dir)
        if issues:
            log(f"\n{len(issues)} issue(s) found:")
            for issue in issues:
                log(f"  - {issue}")
        else:
            log("All files valid.")
        return

    # Exploratory search first (discover all available models)
    log("\nRunning exploratory search...")
    query = ISIMIPQuery()
    search_exploratory(query)

    # Search ISIMIP API for configured models
    log("\nSearching ISIMIP repository for configured models...")
    datasets = search_all_files(query)

    if not datasets:
        log("\nNo datasets found! Check search parameters.")
        return

    # Build and display availability matrix
    matrix = build_availability_matrix(datasets)
    print_availability_report(matrix)

    # Save manifest
    manifest_path = base_dir / "download_manifest.json"
    save_manifest(matrix, manifest_path)

    if args.download or args.dry_run:
        downloads = collect_download_urls(matrix)
        log(f"\n{'DRY RUN: ' if args.dry_run else ''}Downloading {len(downloads)} files to {base_dir}...")

        stats = asyncio.run(
            download_with_organization(downloads, base_dir, dry_run=args.dry_run)
        )

        log(f"\nResults: {stats['downloaded']} downloaded, {stats['skipped']} skipped, {stats['failed']} failed")

        if not args.dry_run and stats["downloaded"] > 0:
            log("\nValidating downloads...")
            issues = validate_downloads(base_dir)
            if issues:
                log(f"\n{len(issues)} validation issue(s):")
                for issue in issues:
                    log(f"  - {issue}")
            else:
                log("All downloads validated successfully.")
    else:
        log("\nTo download, run with --download flag.")
        log("To preview, run with --download --dry-run.")


if __name__ == "__main__":
    main()
