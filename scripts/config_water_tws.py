"""Configuration for TWS (Total Water Storage) monthly processing.

Defines ISIMIP3b search parameters, value_type mappings, aggregation method,
decade windows, and baseline period for the water index pipeline.
"""

# --- ISIMIP3b Data Parameters ---

VARIABLE = "tws"
VARIABLE_LONG_NAME = "Total Water Storage"
VARIABLE_UNITS = "kg m-2"
SECTOR = "water_global"
SIMULATION_ROUND = "ISIMIP3b"
TIMESTEP = "monthly"
SOCIAL_FORCING = "2015soc"

# Aggregation: mean for stock variables (tws, rootmoist), sum for fluxes (dis, qr, potevap, precip)
AGGREGATION_METHOD = "mean"

# Impact models with monthly tws in ISIMIP3b
IMPACT_MODELS = [
    "cwatm",
    "h08",
    "jules-w2",
    "miroc-integ-land",
    # jules-es-vn6p3 listed in plan but not confirmed in API search — include if found
]

# GCMs (ISIMIP3b standard set)
GCMS = [
    "gfdl-esm4",
    "ipsl-cm6a-lr",
    "mpi-esm1-2-hr",
    "mri-esm2-0",
    "ukesm1-0-ll",
]

# SSP scenarios only (no historical model runs — projection-only approach)
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

# --- Temporal Configuration ---

# Projection-only: all data from SSP files (2015-2100)
# 2010s decade uses 2015-2019 (incomplete but acceptable)
DECADES = [2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
MIN_YEAR = 2015  # Start of SSP projections
MAX_YEAR = 2099

def get_decade_years(decade: int) -> tuple:
    """Return (start_year, end_year) for a decade.

    Special case: 2010s uses 2015-2019 (projection start).
    All others: standard 10-year windows.
    """
    if decade == 2010:
        return (2015, 2019)
    return (decade, decade + 9)

# --- Output Value Types (20 total) ---

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
    12: "Annual_Mean",               # Mean of monthly means (= mean of vt0-11)
    13: "Annual_Q05",                 # 5th percentile of annual values
    14: "Annual_Q15",                 # 15th percentile of annual values
    15: "Annual_Q25",                 # 25th percentile of annual values
    16: "Annual_Q50",                 # 50th percentile (median) of annual values
    17: "Annual_Q75",                 # 75th percentile of annual values
    18: "Annual_Q85",                 # 85th percentile of annual values
    19: "Annual_Q95",                 # 95th percentile of annual values
}

N_VALUE_TYPES = 20
N_MONTHLY = 12  # value_types 0-11
N_ANNUAL_STATS = 8  # value_types 12-19

# --- Spatial Grid ---

N_LAT = 360
N_LON = 720
LAT_RESOLUTION = 0.5
LON_RESOLUTION = 0.5

# --- Processing Parameters ---

LAT_CHUNK_SIZE = 36    # Rows per chunk (~2 GB memory target)

# --- Model Normalization ---
# Robust z-score normalization to align impact models before ensemble averaging.
# Each model's values are transformed: target_mean + (value - model_median) / model_IQR * target_sd
# Stats computed from annual-mean TWS during NORM_REF_YEARS.
NORMALIZE_MODELS = True
NORM_TARGET_MEAN = 1000.0
NORM_TARGET_SD = 200.0
NORM_REF_YEARS = (2015, 2024)  # Reference period for computing per-model median/IQR

# --- Scenario Mapping (for QA comparison with old RCP files) ---

SSP_TO_RCP_MAPPING = {
    "ssp126": "rcp26",
    "ssp370": "rcp60",
    "ssp585": "rcp85",
}

# --- Download Configuration ---

DOWNLOAD_BASE_DIR = "data/raw/water_tws"

def get_download_subdir(model: str, gcm: str, scenario: str) -> str:
    """Return subdirectory path for organizing downloads."""
    return f"{model}/{gcm}_{scenario}"
