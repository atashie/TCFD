"""Generic water variable processor for ISIMIP3b monthly data.

Processes any of the 6 water variables (tws, dis, potevap, qr, rootmoist, precip)
into the water index format with 20 value_types.

Produces a NetCDF file with dimensions (lat, lon, scenario, value_type, decade) where:
  - value_types 0-11: Per-month ensemble means (Jan-Dec) within each decade
  - value_type 12: Annual mean (= mean of monthly means)
  - value_types 13-19: Annual quantiles Q05, Q15, Q25, Q50, Q75, Q85, Q95

This is the generalized version of process_water_tws.py, parameterized via
WaterVariableConfig from config_water_variables.py.

Usage:
    python scripts/process_water_variable.py tws
    python scripts/process_water_variable.py dis --data-dir /path/to/data
    python scripts/process_water_variable.py potevap --output /path/to/output.nc
    python scripts/process_water_variable.py --list   # Show available variables
"""

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import xarray as xr

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_water_variables import (
    WaterVariableConfig, SharedConfig, SHARED_CONFIG,
    WATER_VARIABLES, VALUE_TYPE_NAMES,
    get_variable_config, list_variables,
)

# Quantile levels for vt13-19
ANNUAL_QUANTILES = [5, 15, 25, 50, 75, 85, 95]

# Normalization defaults (same as TWS)
NORM_TARGET_MEAN = 1000.0
NORM_TARGET_SD = 200.0
NORM_REF_YEARS = (2015, 2024)


def log(msg: str):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Model normalization
# ---------------------------------------------------------------------------

def compute_normalization_stats(
    inventory: Dict[str, Dict[str, Dict[str, List[Path]]]],
    variable_name: str,
    ref_years: tuple = NORM_REF_YEARS,
    cache_path: Optional[Path] = None,
) -> Dict[str, Dict[str, float]]:
    """Compute per-model median and IQR of annual-mean values from reference period.

    Streams one year at a time per file and spatially subsamples (every 4th cell)
    to keep memory manageable.

    Returns: {model: {"median": float, "iqr": float, "p25": float, "p75": float, "n": int}}
    """
    # Check cache first
    if cache_path and cache_path.exists():
        log(f"Loading normalization stats from cache: {cache_path}")
        with open(cache_path) as f:
            stats = json.load(f)
        for model, s in stats.items():
            log(f"  {model}: median={s['median']:.4g}, IQR={s['iqr']:.4g}")
        return stats

    ref_start, ref_end = ref_years
    stride = 4  # spatial subsampling for efficiency

    stats = {}

    for model in sorted(inventory.keys()):
        log(f"  Computing stats for {model}...")
        all_annual_values = []

        for gcm in sorted(inventory[model].keys()):
            for scenario in sorted(inventory[model][gcm].keys()):
                files = inventory[model][gcm][scenario]
                for fpath in files:
                    try:
                        ds = xr.open_dataset(fpath)
                    except Exception:
                        continue

                    var_name = variable_name if variable_name in ds.data_vars else list(ds.data_vars)[0]
                    da = ds[var_name]

                    try:
                        time_years = da.time.dt.year
                    except Exception:
                        ds.close()
                        continue

                    # Process one year at a time
                    for yr in range(ref_start, ref_end + 1):
                        yr_mask = time_years == yr
                        if yr_mask.values.sum() == 0:
                            continue
                        monthly = da.isel(time=yr_mask).values  # (months, lat, lon)
                        # Subsample spatially
                        monthly_sub = monthly[:, ::stride, ::stride]
                        # Annual mean
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            ann_mean = np.nanmean(monthly_sub, axis=0)  # (lat_sub, lon_sub)
                        valid = ann_mean[~np.isnan(ann_mean)]
                        if len(valid) > 0:
                            all_annual_values.append(valid)

                    ds.close()

        if all_annual_values:
            pooled = np.concatenate(all_annual_values)
            p25 = float(np.percentile(pooled, 25))
            p75 = float(np.percentile(pooled, 75))
            median = float(np.median(pooled))
            iqr = p75 - p25

            stats[model] = {
                "median": median,
                "iqr": iqr,
                "p25": p25,
                "p75": p75,
                "n": int(len(pooled)),
            }
            log(f"    {model}: median={median:.4g}, IQR={iqr:.4g} "
                f"(P25={p25:.4g}, P75={p75:.4g}, n={len(pooled):,})")
        else:
            log(f"    WARNING: no data for {model}")

    # Save cache
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(stats, f, indent=2)
        log(f"  Normalization stats cached: {cache_path}")

    return stats


def apply_normalization(
    data: np.ndarray,
    model_median: float,
    model_iqr: float,
    target_mean: float = NORM_TARGET_MEAN,
    target_sd: float = NORM_TARGET_SD,
) -> np.ndarray:
    """Apply robust z-score normalization: target_mean + (value - median) / IQR * target_sd.

    When normalizing a single outlier model to the ensemble, set target_mean and
    target_sd to the reference models' pooled median and IQR respectively.
    """
    if model_iqr <= 0:
        return data
    return target_mean + (data - model_median) / model_iqr * target_sd


# ---------------------------------------------------------------------------
# File filtering
# ---------------------------------------------------------------------------

def _filter_standard_run(nc_files: List[Path]) -> List[Path]:
    """Filter NetCDF files to keep only the standard social forcing variant.

    Preference order:
      1. Files with '2015soc_default' (standard run)
      2. Files with '2015soc-from-histsoc_default' (transition variant)
    Excludes: '1850soc' (pre-industrial) and '2015co2' (fixed CO2 sensitivity).
    """
    standard = [f for f in nc_files if "_2015soc_default_" in f.name]
    if standard:
        return standard

    histsoc = [f for f in nc_files if "_2015soc-from-histsoc_default_" in f.name]
    if histsoc:
        return histsoc

    fallback = [f for f in nc_files
                if "2015soc" in f.name and "1850soc" not in f.name and "2015co2" not in f.name]
    return fallback if fallback else nc_files


# ---------------------------------------------------------------------------
# WaterVariableProcessor
# ---------------------------------------------------------------------------

class WaterVariableProcessor:
    """Generic processor for water index variables.

    Parameterized by WaterVariableConfig (variable-specific) and
    SharedConfig (grid, decades, scenarios).
    """

    def __init__(
        self,
        var_config: WaterVariableConfig,
        shared: SharedConfig = SHARED_CONFIG,
        norm_stats: Optional[Dict[str, Dict[str, float]]] = None,
    ):
        self.var = var_config
        self.shared = shared
        self.norm_stats = norm_stats

        # Precompute coordinate arrays
        self.lats = np.arange(89.75, -90, -0.5, dtype=np.float32)[:shared.n_lat]
        self.lons = np.arange(-179.75, 180, 0.5, dtype=np.float32)[:shared.n_lon]

    # --- File discovery ---

    def discover_files(self, data_dir: Path) -> Dict[str, Dict[str, Dict[str, List[Path]]]]:
        """Discover NetCDF files organized as model/gcm_scenario/*.nc

        Filters to keep only the standard social forcing variant.
        Returns: {model: {gcm: {scenario: [paths]}}}
        """
        inventory: Dict[str, Dict[str, Dict[str, List[Path]]]] = {}

        if not data_dir.exists():
            return inventory

        # Only include models from the variable config (if specified)
        allowed_models = set(self.var.models) if self.var.models else None

        # Try organized subdirectory structure
        for model_dir in sorted(data_dir.iterdir()):
            if not model_dir.is_dir():
                continue
            model = model_dir.name
            if allowed_models and model not in allowed_models:
                continue
            inventory[model] = {}

            for sub in sorted(model_dir.iterdir()):
                if not sub.is_dir():
                    continue
                parts = sub.name.split("_", 1)
                if len(parts) != 2:
                    continue
                gcm, scenario = parts
                if gcm not in inventory[model]:
                    inventory[model][gcm] = {}
                nc_files = sorted(sub.glob("*.nc"))
                nc_files = _filter_standard_run(nc_files)
                if nc_files:
                    inventory[model][gcm][scenario] = nc_files

        # Fallback: flat directory, parse filenames
        if not inventory:
            known_models = self.var.models or []
            nc_files = sorted(data_dir.rglob("*.nc"))
            for fpath in nc_files:
                name = fpath.stem.lower()
                model = gcm = scenario = None

                for m in known_models:
                    if m in name:
                        model = m
                        break
                for g in self.shared.gcms:
                    if g.replace("-", "") in name.replace("-", ""):
                        gcm = g
                        break
                for s in self.shared.scenarios:
                    if s in name:
                        scenario = s
                        break

                if model and gcm and scenario:
                    inventory.setdefault(model, {}).setdefault(gcm, {}).setdefault(scenario, []).append(fpath)

            # Filter each group
            for model in inventory:
                for gcm in inventory[model]:
                    for scenario in inventory[model][gcm]:
                        inventory[model][gcm][scenario] = _filter_standard_run(inventory[model][gcm][scenario])

        return inventory

    # --- Data loading ---

    def load_monthly_data(
        self, files: List[Path], year_start: int, year_end: int,
    ) -> Optional[xr.DataArray]:
        """Load monthly data from files, filtering to year range.

        Applies unit conversion if configured.
        """
        arrays = []
        for fpath in files:
            try:
                ds = xr.open_dataset(fpath, use_cftime=True)
            except Exception:
                try:
                    ds = xr.open_dataset(fpath)
                except Exception as e:
                    log(f"    WARNING: Cannot open {fpath.name}: {e}")
                    continue

            var_name = self.var.name if self.var.name in ds.data_vars else list(ds.data_vars)[0]
            da = ds[var_name]

            # Apply unit conversion
            if self.var.unit_conversion_factor != 1.0:
                da = da * self.var.unit_conversion_factor

            # Filter by year range
            try:
                years = da.time.dt.year
                mask = (years >= year_start) & (years <= year_end)
                da = da.isel(time=mask)
            except Exception:
                pass

            if len(da.time) > 0:
                arrays.append(da)
            ds.close()

        if not arrays:
            return None

        combined = xr.concat(arrays, dim="time")
        combined = combined.sortby("time")
        return combined

    # --- Core chunk processing ---

    def process_chunk(
        self,
        inventory: Dict[str, Dict[str, Dict[str, List[Path]]]],
        lat_start: int,
        lat_end: int,
    ) -> np.ndarray:
        """Process a latitude chunk for all scenarios, decades, and value_types.

        Returns array of shape (lat_chunk, n_lon, n_scenarios, n_value_types, n_decades).
        """
        n_lat_chunk = lat_end - lat_start
        n_lon = self.shared.n_lon
        scenarios = list(self.shared.scenarios)
        decades = list(self.shared.decades)
        n_scenarios = len(scenarios)
        n_decades = len(decades)

        result = np.full(
            (n_lat_chunk, n_lon, n_scenarios, self.shared.n_value_types, n_decades),
            np.nan, dtype=np.float32,
        )

        for s_idx, scenario in enumerate(scenarios):
            log(f"    Scenario {scenario}...")

            # Collect ensemble members
            ensemble_monthly = []
            for model in inventory:
                for gcm in inventory[model]:
                    if scenario not in inventory[model][gcm]:
                        continue
                    files = inventory[model][gcm][scenario]
                    da = self.load_monthly_data(files, self.shared.min_year, self.shared.max_year)
                    if da is None:
                        continue

                    chunk_lats = self.lats[lat_start:lat_end]
                    try:
                        da_chunk = da.sel(lat=chunk_lats, method="nearest")
                    except Exception:
                        da_chunk = da.isel(lat=slice(lat_start, lat_end))

                    # Apply per-model normalization if stats available
                    if self.norm_stats and model in self.norm_stats:
                        ms = self.norm_stats[model]
                        da_chunk = da_chunk.copy(data=apply_normalization(
                            da_chunk.values, ms["median"], ms["iqr"],
                            target_mean=ms.get("target_mean", NORM_TARGET_MEAN),
                            target_sd=ms.get("target_sd", NORM_TARGET_SD)))

                    ensemble_monthly.append(da_chunk)

            if not ensemble_monthly:
                log(f"      No data for {scenario}")
                continue

            n_members = len(ensemble_monthly)
            log(f"      {n_members} ensemble members loaded")

            for d_idx, decade in enumerate(decades):
                year_start, year_end = self.shared.get_decade_years(decade)

                # --- Value types 0-11: Monthly ensemble means ---
                for month in range(12):
                    month_values = []
                    for da in ensemble_monthly:
                        try:
                            years_mask = (da.time.dt.year >= year_start) & (da.time.dt.year <= year_end)
                            month_mask = da.time.dt.month == (month + 1)
                            subset = da.isel(time=years_mask & month_mask)
                            if len(subset.time) > 0:
                                month_values.append(subset.values)
                        except Exception:
                            continue

                    if month_values:
                        stacked = np.concatenate(month_values, axis=0)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            result[:, :, s_idx, month, d_idx] = np.nanmean(stacked, axis=0)

                # --- Value type 12: Annual mean (= mean of vt0-11) ---
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result[:, :, s_idx, 12, d_idx] = np.nanmean(
                        result[:, :, s_idx, 0:12, d_idx], axis=2)

                # --- Value types 13-19: Annual quantiles Q05..Q95 ---
                # Pool raw annual values across all ensemble members and years
                member_annual_values = []
                for da in ensemble_monthly:
                    try:
                        years_mask = (da.time.dt.year >= year_start) & (da.time.dt.year <= year_end)
                        subset = da.isel(time=years_mask)
                        if len(subset.time) == 0:
                            continue

                        # Always use mean to keep quantiles in the same
                        # units as vt0-12 (rate units for flux variables).
                        annual = subset.groupby("time.year").mean(dim="time")

                        # annual has shape (n_years_this_member, lat_chunk, lon)
                        member_annual_values.append(annual.values)
                    except Exception:
                        continue

                if member_annual_values:
                    # Pool: (n_members * n_years, lat_chunk, lon)
                    pooled = np.concatenate(member_annual_values, axis=0)

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        for q_idx, q_level in enumerate(ANNUAL_QUANTILES):
                            result[:, :, s_idx, 13 + q_idx, d_idx] = np.nanpercentile(
                                pooled, q_level, axis=0)

        return result

    # --- Output writer ---

    def write_output(self, result: np.ndarray, output_path: Path):
        """Write processed result to NetCDF matching target format."""
        scenarios = list(self.shared.scenarios)
        decades = list(self.shared.decades)

        attrs = {
            "Conventions": "CF-1.8",
            "title": f"Water Index Underlying Data - {self.var.long_name}",
            "variable": self.var.name,
            "variable_long_name": self.var.long_name,
            "units": self.var.units_output,
            "units_raw": self.var.units_raw,
            "source": self.var.simulation_round,
            "scenarios": ", ".join(scenarios),
            "data_source": "SSP projections only (no historical model runs)",
            "aggregation_method": self.var.aggregation,
            "unit_conversion_factor": str(self.var.unit_conversion_factor),
            "impact_models": ", ".join(self.var.models),
            "gcms": ", ".join(self.shared.gcms),
            "value_type_names": str(VALUE_TYPE_NAMES),
        }

        if self.norm_stats:
            attrs["normalization_method"] = "robust_zscore"
            attrs["normalization_formula"] = (
                "target_mean + (value - model_median) / model_IQR * target_sd"
            )
            attrs["normalization_ref_period"] = f"{NORM_REF_YEARS[0]}-{NORM_REF_YEARS[1]}"
            attrs["normalization_models_normalized"] = ", ".join(sorted(self.norm_stats.keys()))
            for model, s in self.norm_stats.items():
                attrs[f"normalization_{model}_median"] = s["median"]
                attrs[f"normalization_{model}_iqr"] = s["iqr"]
                attrs[f"normalization_{model}_target_mean"] = s.get("target_mean", NORM_TARGET_MEAN)
                attrs[f"normalization_{model}_target_sd"] = s.get("target_sd", NORM_TARGET_SD)

        ds = xr.Dataset(
            {
                self.var.name: (
                    ["lat", "lon", "scenario", "value_type", "decade"],
                    result,
                ),
            },
            coords={
                "lat": ("lat", self.lats, {"units": "degrees_north", "long_name": "latitude"}),
                "lon": ("lon", self.lons, {"units": "degrees_east", "long_name": "longitude"}),
                "scenario": ("scenario", scenarios),
                "value_type": ("value_type", list(range(self.shared.n_value_types))),
                "decade": ("decade", decades),
            },
            attrs=attrs,
        )

        encoding = {
            self.var.name: {
                "dtype": "float32",
                "zlib": True,
                "complevel": 4,
                "_FillValue": np.float32(np.nan),
            },
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(output_path, encoding=encoding)
        log(f"Output saved: {output_path}")
        log(f"File size: {output_path.stat().st_size / (1024 * 1024):.1f} MB")

    # --- Main orchestrator ---

    def run(self, data_dir: Path, output_path: Path, chunk_size: int = None):
        """Run the full processing pipeline."""
        chunk_size = chunk_size or self.shared.lat_chunk_size

        log("=" * 70)
        log(f"Processing {self.var.long_name} ({self.var.name})")
        log(f"  Aggregation: {self.var.aggregation}")
        log(f"  Unit conversion: {self.var.units_raw} -> {self.var.units_output} "
            f"(x{self.var.unit_conversion_factor})")
        log(f"  Models: {', '.join(self.var.models) if self.var.models else 'climate forcing'}")
        log("=" * 70)

        # Discover files
        log(f"\nDiscovering files in {data_dir}...")
        inventory = self.discover_files(data_dir)

        if not inventory:
            log("ERROR: No data files found!")
            log(f"Expected: {data_dir}/{{model}}/{{gcm}}_{{scenario}}/*.nc")
            sys.exit(1)

        # Report inventory
        n_members = 0
        for model in sorted(inventory):
            for gcm in sorted(inventory[model]):
                n_members += 1
                scenarios = sorted(inventory[model][gcm].keys())
                n_files = sum(len(inventory[model][gcm][s]) for s in scenarios)
                log(f"  {model} / {gcm}: {', '.join(scenarios)} ({n_files} files)")
        log(f"\nTotal ensemble members: {n_members}")

        # Initialize output
        full_result = np.full(
            (self.shared.n_lat, self.shared.n_lon, len(self.shared.scenarios),
             self.shared.n_value_types, len(self.shared.decades)),
            np.nan, dtype=np.float32,
        )

        # Process in latitude chunks
        n_chunks = (self.shared.n_lat + chunk_size - 1) // chunk_size
        log(f"\nProcessing {self.shared.n_lat} lat rows in {n_chunks} chunks of {chunk_size}...")

        for chunk_idx in range(n_chunks):
            lat_start = chunk_idx * chunk_size
            lat_end = min(lat_start + chunk_size, self.shared.n_lat)

            log(f"\n{'='*50}")
            log(f"Chunk {chunk_idx + 1}/{n_chunks}: lats [{lat_start}:{lat_end}] "
                f"({self.lats[lat_start]:.2f} to {self.lats[lat_end - 1]:.2f})")
            log(f"{'='*50}")

            chunk_result = self.process_chunk(inventory, lat_start, lat_end)
            full_result[lat_start:lat_end, :, :, :, :] = chunk_result

        # Write output
        log(f"\nWriting output to {output_path}...")
        self.write_output(full_result, output_path)

        # Summary
        log("\nSummary statistics:")
        for s_idx, scenario in enumerate(self.shared.scenarios):
            data = full_result[:, :, s_idx, 12, :]
            valid = ~np.isnan(data)
            if valid.sum() > 0:
                log(f"  {scenario}: {valid.sum():,} valid cells, "
                    f"range [{np.nanmin(data):.4g}, {np.nanmax(data):.4g}]")
            else:
                log(f"  {scenario}: no valid cells")

        log("\n" + "=" * 70)
        log("Processing complete!")
        log("=" * 70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Process water variable monthly ISIMIP3b data into water index format",
    )
    parser.add_argument(
        "variable", nargs="?", default=None,
        help=f"Variable to process ({', '.join(list_variables())})",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available variables and exit",
    )
    parser.add_argument(
        "--data-dir", type=Path, default=None,
        help="Directory with downloaded raw data (default: data/raw/water_{variable})",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help=r"Output file (default: C:\Cai_data\WaterIndex\waterIndexUnderlyingData_{var}_ssp.nc)",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=None,
        help=f"Latitude rows per chunk (default: {SHARED_CONFIG.lat_chunk_size})",
    )
    parser.add_argument(
        "--normalize", action="store_true",
        help="Enable per-model robust z-score normalization before ensemble averaging",
    )
    parser.add_argument(
        "--normalize-models", type=str, default=None,
        help=("Comma-separated list of models to normalize (others become the reference). "
              "Target median/IQR are derived from the reference models. "
              "E.g., --normalize-models h08 normalizes only h08 to match the other models."),
    )
    parser.add_argument(
        "--norm-cache", type=Path, default=None,
        help="Path to cache normalization statistics JSON (default: data/normalization_stats_{var}.json)",
    )
    args = parser.parse_args()

    if args.list:
        log("Available water variables:\n")
        log(f"  {'Name':<12} {'Aggregation':<8} {'Unit Conversion':<30} {'Models'}")
        log(f"  {'-'*12} {'-'*8} {'-'*30} {'-'*40}")
        for name in list_variables():
            cfg = WATER_VARIABLES[name]
            conv = f"{cfg.units_raw} -> {cfg.units_output}" if cfg.unit_conversion_factor != 1.0 else "none"
            models = ", ".join(cfg.models) if cfg.models else "(climate forcing)"
            log(f"  {name:<12} {cfg.aggregation:<8} {conv:<30} {models}")
        return

    if args.variable is None:
        parser.error("Please specify a variable (e.g., 'tws') or use --list")

    var_config = get_variable_config(args.variable)
    data_dir = args.data_dir or (PROJECT_ROOT / var_config.download_subdir)
    output_path = args.output or Path(rf"C:\Cai_data\WaterIndex\{var_config.output_filename}")

    # Normalization pre-pass (if requested)
    norm_stats = None
    if args.normalize or args.normalize_models:
        # Need to discover files first for normalization stats
        temp_processor = WaterVariableProcessor(var_config)
        inventory = temp_processor.discover_files(data_dir)
        if not inventory:
            log("ERROR: No data files found for normalization!")
            sys.exit(1)

        norm_cache = args.norm_cache or (PROJECT_ROOT / "data" / f"normalization_stats_{var_config.name}.json")
        all_stats = compute_normalization_stats(
            inventory, var_config.name, cache_path=norm_cache)
        if not all_stats:
            log("WARNING: No normalization stats computed, proceeding without normalization")
        elif args.normalize_models:
            # Selective normalization: normalize only specified models,
            # using the other models' pooled stats as target
            models_to_normalize = [m.strip() for m in args.normalize_models.split(",")]
            reference_models = [m for m in all_stats if m not in models_to_normalize]

            if not reference_models:
                log("ERROR: --normalize-models includes all models, no reference models left!")
                sys.exit(1)

            # Pool reference models' values to get target median/IQR
            ref_medians = [all_stats[m]["median"] for m in reference_models]
            ref_iqrs = [all_stats[m]["iqr"] for m in reference_models]
            target_median = float(np.median(ref_medians))
            target_iqr = float(np.median(ref_iqrs))

            log(f"\n--- Selective Model Normalization ---")
            log(f"Models to normalize: {', '.join(models_to_normalize)}")
            log(f"Reference models: {', '.join(reference_models)}")
            log(f"Target (from reference): median={target_median:.4g}, IQR={target_iqr:.4g}")

            norm_stats = {}
            for model in models_to_normalize:
                if model not in all_stats:
                    log(f"WARNING: {model} not found in stats, skipping")
                    continue
                norm_stats[model] = {
                    **all_stats[model],
                    "target_mean": target_median,
                    "target_sd": target_iqr,
                }
                ms = all_stats[model]
                log(f"  {model}: median={ms['median']:.4g}, IQR={ms['iqr']:.4g} "
                    f"-> target median={target_median:.4g}, IQR={target_iqr:.4g}")
        else:
            # Full normalization (all models to synthetic target)
            log(f"\n--- Model Normalization ---")
            log(f"Reference period: {NORM_REF_YEARS[0]}-{NORM_REF_YEARS[1]}")
            log(f"Target: mean={NORM_TARGET_MEAN}, SD={NORM_TARGET_SD}")
            norm_stats = all_stats

    processor = WaterVariableProcessor(var_config, norm_stats=norm_stats)
    processor.run(data_dir, output_path, args.chunk_size)


if __name__ == "__main__":
    main()
