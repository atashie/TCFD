"""NetCDF output generation for processed ISIMIP data.

Creates standardized NetCDF4 files with dimensions:
(lon, lat, decade, scenario, value_class)

Value classes (6 types):
1. Absolute smoothed median value
2. Percentile rank (1-100)
3. Decadal trend (slope)
4. Trend significance (p-value)
5. Lower confidence bound
6. Upper confidence bound
"""

from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

import numpy as np
import xarray as xr


# Value class definitions
VALUE_CLASSES = [
    "smoothed_median",
    "percentile",
    "trend",
    "significance",
    "lower_bound",
    "upper_bound",
]


def create_output_dataset(
    lat: np.ndarray,
    lon: np.ndarray,
    decades: List[int],
    scenarios: List[str],
    variable: str,
    data: Optional[np.ndarray] = None,
) -> xr.Dataset:
    """Create standardized xarray Dataset for output.

    Args:
        lat: Array of latitude values.
        lon: Array of longitude values.
        decades: List of decade numbers (10, 20, ..., 90).
        scenarios: List of scenario names (e.g., ["ssp126", "ssp585"]).
        variable: Variable name for the dataset.
        data: Optional data array with shape
              (lon, lat, decades, scenarios, value_classes).

    Returns:
        xarray Dataset with CF-compliant metadata.
    """
    # Create coordinate arrays
    coords = {
        "lon": ("lon", lon),
        "lat": ("lat", lat),
        "decade": ("decade", decades),
        "scenario": ("scenario", scenarios),
        "value_class": ("value_class", VALUE_CLASSES),
    }

    # Initialize data array if not provided
    if data is None:
        data = np.full(
            (len(lon), len(lat), len(decades), len(scenarios), len(VALUE_CLASSES)),
            np.nan,
            dtype=np.float32,
        )

    # Create data variable
    data_vars = {
        variable: (
            ["lon", "lat", "decade", "scenario", "value_class"],
            data,
            {
                "long_name": f"{variable} climate features",
                "units": "varies by value_class",
                "_FillValue": np.float32(np.nan),
            },
        )
    }

    # Create dataset
    ds = xr.Dataset(data_vars, coords=coords)

    # Add CF-compliant attributes to coordinates
    ds.lon.attrs.update(
        {
            "long_name": "longitude",
            "units": "degrees_east",
            "standard_name": "longitude",
            "axis": "X",
        }
    )

    ds.lat.attrs.update(
        {
            "long_name": "latitude",
            "units": "degrees_north",
            "standard_name": "latitude",
            "axis": "Y",
        }
    )

    ds.decade.attrs.update(
        {
            "long_name": "decade",
            "units": "1",
            "description": "Decade index (10=2010s, 20=2020s, ...)",
        }
    )

    ds.scenario.attrs.update(
        {
            "long_name": "climate scenario",
            "description": "CMIP6 SSP scenario or CMIP5 RCP scenario",
        }
    )

    ds.value_class.attrs.update(
        {
            "long_name": "feature type",
            "description": (
                "1=smoothed_median, 2=percentile, 3=trend, "
                "4=significance, 5=lower_bound, 6=upper_bound"
            ),
        }
    )

    # Add global attributes
    ds.attrs.update(
        {
            "Conventions": "CF-1.8",
            "title": f"Processed ISIMIP {variable} data",
            "institution": "ISIMIP Pipeline",
            "source": "ISIMIP repository",
            "history": f"Created {datetime.now().isoformat()}",
            "references": "https://www.isimip.org/",
        }
    )

    return ds


def write_netcdf(
    ds: xr.Dataset,
    output_path: Path,
    compression: bool = True,
    complevel: int = 4,
) -> None:
    """Write xarray Dataset to NetCDF file.

    Args:
        ds: xarray Dataset to write.
        output_path: Path to output file.
        compression: Whether to enable zlib compression.
        complevel: Compression level (1-9).
    """
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Set up encoding for compression
    encoding = {}
    if compression:
        for var in ds.data_vars:
            encoding[var] = {
                "zlib": True,
                "complevel": complevel,
                "dtype": "float32",
            }

    # Write to NetCDF4
    ds.to_netcdf(
        output_path,
        format="NETCDF4",
        encoding=encoding if compression else None,
    )


class OutputWriter:
    """Manages writing processed data to NetCDF files.

    Provides a high-level interface for creating standardized
    output files from processed ISIMIP data.
    """

    def __init__(
        self,
        output_dir: Path,
        compression: bool = True,
        complevel: int = 4,
    ):
        """Initialize output writer.

        Args:
            output_dir: Directory for output files.
            compression: Whether to enable compression.
            complevel: Compression level (1-9).
        """
        self.output_dir = Path(output_dir)
        self.compression = compression
        self.complevel = complevel

    def write(
        self,
        ds: xr.Dataset,
        filename: str,
    ) -> Path:
        """Write dataset to file.

        Args:
            ds: xarray Dataset to write.
            filename: Output filename (without path).

        Returns:
            Path to written file.
        """
        output_path = self.output_dir / filename
        write_netcdf(
            ds,
            output_path,
            compression=self.compression,
            complevel=self.complevel,
        )
        return output_path

    def create_and_write(
        self,
        lat: np.ndarray,
        lon: np.ndarray,
        decades: List[int],
        scenarios: List[str],
        variable: str,
        data: np.ndarray,
        filename: str,
    ) -> Path:
        """Create dataset and write to file.

        Args:
            lat: Latitude array.
            lon: Longitude array.
            decades: List of decades.
            scenarios: List of scenarios.
            variable: Variable name.
            data: Feature data array.
            filename: Output filename.

        Returns:
            Path to written file.
        """
        ds = create_output_dataset(
            lat=lat,
            lon=lon,
            decades=decades,
            scenarios=scenarios,
            variable=variable,
            data=data,
        )
        return self.write(ds, filename)
