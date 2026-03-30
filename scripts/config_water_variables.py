"""Configuration for all 6 water variables in the water index pipeline.

Defines per-variable parameters: aggregation method, unit conversion,
impact models, and data source. Shared constants (decades, grid, value_types)
are defined once and reused across all variables.

Usage:
    from config_water_variables import WATER_VARIABLES, SHARED_CONFIG
    cfg = WATER_VARIABLES["tws"]
    print(cfg.long_name, cfg.aggregation, cfg.models)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# Per-variable configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WaterVariableConfig:
    """Configuration for a single water variable."""

    name: str                        # ISIMIP variable code
    long_name: str                   # Human-readable name
    units_raw: str                   # Units in raw ISIMIP files
    units_output: str                # Units after conversion (if any)
    aggregation: str                 # "mean" for stocks, "sum" for fluxes
    unit_conversion_factor: float    # Multiply raw values by this; 1.0 = no conversion
    models: List[str]                # Impact models with monthly data in ISIMIP3b
    sector: str = "water_global"     # ISIMIP sector
    product: str = "OutputData"      # "OutputData" for model outputs, "InputData" for climate forcing
    simulation_round: str = "ISIMIP3b"
    timestep: str = "monthly"
    social_forcing: str = "2015soc"
    notes: str = ""

    @property
    def download_subdir(self) -> str:
        return f"data/raw/water_{self.name}"

    @property
    def output_filename(self) -> str:
        return f"waterIndexUnderlyingData_{self.name}_ssp.nc"


# Conversion constant: kg/m2/s → mm/month
# 1 kg/m2/s = 1 mm/s (water density ~1000 kg/m3, 1 mm = 1 kg/m2)
# seconds per average month: 60 * 60 * 24 * 30.4375 = 2,629,800
KG_M2_S_TO_MM_MONTH = 60 * 60 * 24 * 30.4375  # 2,629,800


WATER_VARIABLES: Dict[str, WaterVariableConfig] = {
    "tws": WaterVariableConfig(
        name="tws",
        long_name="Total Water Storage",
        units_raw="kg m-2",
        units_output="kg m-2",
        aggregation="mean",
        unit_conversion_factor=1.0,
        models=["cwatm", "h08", "jules-w2", "miroc-integ-land"],
        notes="5 GCMs, best coverage among water variables",
    ),
    "dis": WaterVariableConfig(
        name="dis",
        long_name="Streamflow (Discharge)",
        units_raw="m3 s-1",
        units_output="m3 s-1",
        aggregation="mean",
        unit_conversion_factor=1.0,
        models=["cwatm", "h08", "jules-w2", "miroc-integ-land", "watergap2-2e"],
        notes="6 models, largest ensemble",
    ),
    "potevap": WaterVariableConfig(
        name="potevap",
        long_name="Potential Evapotranspiration",
        units_raw="kg m-2 s-1",
        units_output="mm month-1",
        aggregation="sum",
        unit_conversion_factor=KG_M2_S_TO_MM_MONTH,
        models=["cwatm", "h08", "miroc-integ-land", "watergap2-2e"],
        notes="Flux variable, requires unit conversion to mm/month",
    ),
    "qr": WaterVariableConfig(
        name="qr",
        long_name="Groundwater Runoff (Recharge)",
        units_raw="kg m-2 s-1",
        units_output="mm month-1",
        aggregation="sum",
        unit_conversion_factor=KG_M2_S_TO_MM_MONTH,
        models=["cwatm", "h08", "miroc-integ-land", "watergap2-2e"],
        notes="Flux variable, requires unit conversion to mm/month",
    ),
    "rootmoist": WaterVariableConfig(
        name="rootmoist",
        long_name="Root Zone Soil Moisture",
        units_raw="kg m-2",
        units_output="% max capacity",
        aggregation="mean",
        unit_conversion_factor=100.0 / 1187.2939,
        models=["web-dhm-sg"],
        notes=(
            "WEB-DHM-SG only (1 model). miroc-integ-land excluded: rootmoist includes "
            "glacier/ice mass (max ~45,000 kg/m², exceeds ISIMIP valid_max 10,000 by 4.5x), "
            "84%% Arctic zeros, and undefined root zone depth. WEB-DHM-SG root zone ~1.5m, "
            "max 1,187.29 kg/m². Values normalized to %% of max capacity (1187.29 kg/m² = 100%%)."
        ),
    ),
    "precip": WaterVariableConfig(
        name="pr",
        long_name="Precipitation",
        units_raw="kg m-2 s-1",
        units_output="mm month-1",
        aggregation="sum",
        unit_conversion_factor=KG_M2_S_TO_MM_MONTH,
        models=[],  # Climate forcing, not model output
        sector="",
        product="InputData",
        notes="Climate forcing (InputData), not model output. Different download path.",
    ),
}


# ---------------------------------------------------------------------------
# Shared constants (identical across all variables)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SharedConfig:
    """Constants shared across all water variables."""

    # GCMs
    gcms: tuple = (
        "gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr",
        "mri-esm2-0", "ukesm1-0-ll",
    )

    # Scenarios (projection-only)
    scenarios: tuple = ("ssp126", "ssp370", "ssp585")

    # Temporal
    decades: tuple = (2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090)
    min_year: int = 2015
    max_year: int = 2099

    # Spatial grid
    n_lat: int = 360
    n_lon: int = 720
    lat_resolution: float = 0.5
    lon_resolution: float = 0.5

    # Processing
    lat_chunk_size: int = 36

    # Value types (20 total)
    n_value_types: int = 20
    n_monthly: int = 12
    n_annual_stats: int = 8

    # Scenario mapping for QA
    ssp_to_rcp: dict = field(default_factory=lambda: {
        "ssp126": "rcp26", "ssp370": "rcp60", "ssp585": "rcp85",
    })

    def get_decade_years(self, decade: int) -> tuple:
        """Return (start_year, end_year) for a decade."""
        if decade == 2010:
            return (2015, 2019)
        return (decade, decade + 9)


# Singleton shared config
SHARED_CONFIG = SharedConfig()

# Value type names (shared across all variables)
VALUE_TYPE_NAMES = {
    0: "Jan_Ensemble_Mean",
    1: "Feb_Ensemble_Mean",
    2: "Mar_Ensemble_Mean",
    3: "Apr_Ensemble_Mean",
    4: "May_Ensemble_Mean",
    5: "Jun_Ensemble_Mean",
    6: "Jul_Ensemble_Mean",
    7: "Aug_Ensemble_Mean",
    8: "Sep_Ensemble_Mean",
    9: "Oct_Ensemble_Mean",
    10: "Nov_Ensemble_Mean",
    11: "Dec_Ensemble_Mean",
    12: "Annual_Mean",
    13: "Annual_Q05",
    14: "Annual_Q15",
    15: "Annual_Q25",
    16: "Annual_Q50",
    17: "Annual_Q75",
    18: "Annual_Q85",
    19: "Annual_Q95",
}


def get_variable_config(name: str) -> WaterVariableConfig:
    """Get config for a variable by name. Raises KeyError if not found."""
    if name not in WATER_VARIABLES:
        available = ", ".join(sorted(WATER_VARIABLES.keys()))
        raise KeyError(f"Unknown variable '{name}'. Available: {available}")
    return WATER_VARIABLES[name]


def list_variables() -> List[str]:
    """List all configured variable names."""
    return sorted(WATER_VARIABLES.keys())
