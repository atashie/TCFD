"""Keyword-based natural language query parsing for ISIMIP searches."""

import re
from dataclasses import dataclass
from typing import Optional, Dict, Any

from isimip_pipeline.search.isimip_query import SearchFilters


# Known ISIMIP variable mappings for keyword matching
VARIABLE_KEYWORDS = {
    # Exposure metrics (ISIMIP3b SecondaryOutputData)
    "drought": "led",
    "drought exposure": "led",
    "heatwave": "leh",
    "heat wave": "leh",
    "heat exposure": "leh",
    "wildfire exposure": "lew",
    "flood": "ler",
    "river flood": "ler",
    "flood exposure": "ler",
    "crop failure": "lec",
    "crop exposure": "lec",
    # Fire/burnt area
    "fire": "burntarea",
    "burnt area": "burntarea",
    "burned area": "burntarea",
    "wildfire": "burntarea",
    "fire emissions": "ffire",
    # Hydrology
    "evapotranspiration": "potevap",
    "pet": "potevap",
    "potential evapotranspiration": "potevap",
    "evaporation": "evap",
    "discharge": "dis",
    "river discharge": "dis",
    "runoff": "qtot",
    "total runoff": "qtot",
    "surface runoff": "qs",
    "groundwater": "groundwstor",
    "snow": "swe",
    "snow water": "swe",
    # Temperature
    "temperature": "tas",
    "air temperature": "tas",
    "max temperature": "tasmax",
    "min temperature": "tasmin",
    # Precipitation
    "precipitation": "pr",
    "rainfall": "pr",
    "snowfall": "prsn",
    # Agriculture - common crops
    "maize": "yield-mai",
    "maize yield": "yield-mai",
    "corn": "yield-mai",
    "wheat": "yield-swh",
    "wheat yield": "yield-swh",
    "spring wheat": "yield-swh",
    "winter wheat": "yield-wwh",
    "rice": "yield-ric",
    "rice yield": "yield-ric",
    "soybean": "yield-soy",
    "soy": "yield-soy",
    "barley": "yield-bar",
    "barley yield": "yield-bar",
    "sorghum": "yield-sor",
    "millet": "yield-mil",
    "potato": "yield-pot",
    "cotton": "yield-cot",
    "crop yield": "yield",
    # Marine fishery
    "fish": "tcb",
    "fishery": "tcb",
    "fish biomass": "tcb",
    "large fish": "b30cm",
    "large fish biomass": "b30cm",
    "total consumer biomass": "tcb",
    # Vegetation/Carbon
    "gpp": "gpp",
    "gross primary": "gpp",
    "npp": "npp",
    "net primary": "npp",
    "carbon": "cveg",
    "vegetation carbon": "cveg",
    "wood carbon": "cwood",
    "soil carbon": "csoil",
}

# Known scenario patterns
SCENARIO_PATTERNS = {
    r"ssp\s*5[\-\.]?8\.?5": "ssp585",
    r"ssp\s*3[\-\.]?7\.?0": "ssp370",
    r"ssp\s*2[\-\.]?4\.?5": "ssp245",
    r"ssp\s*1[\-\.]?2\.?6": "ssp126",
    r"ssp585": "ssp585",
    r"ssp370": "ssp370",
    r"ssp245": "ssp245",
    r"ssp126": "ssp126",
    r"rcp\s*8\.?5": "rcp85",
    r"rcp\s*6\.?0": "rcp60",
    r"rcp\s*4\.?5": "rcp45",
    r"rcp\s*2\.?6": "rcp26",
    r"rcp85": "rcp85",
    r"rcp60": "rcp60",
    r"rcp45": "rcp45",
    r"rcp26": "rcp26",
    r"historical": "historical",
}

# Known simulation rounds
SIMULATION_ROUNDS = ["ISIMIP2a", "ISIMIP2b", "ISIMIP3a", "ISIMIP3b"]

# Crops only available in specific rounds
ISIMIP2A_ONLY_CROPS = {"bar", "rye", "pot", "cas", "cot", "ben"}


@dataclass
class ParsedQuery:
    """Result of parsing a natural language query.

    Attributes:
        filters: SearchFilters extracted from the query.
        explanation: Human-readable explanation of what was found.
        raw_response: Raw response data (for compatibility).
    """

    filters: SearchFilters
    explanation: str = ""
    raw_response: Optional[Dict[str, Any]] = None


def parse_query(query: str) -> ParsedQuery:
    """Parse natural language query using keyword matching.

    Args:
        query: Natural language query.

    Returns:
        ParsedQuery with extracted filters.
    """
    query_lower = query.lower()

    # Extract variable
    variable = None
    matched_keyword = None
    # Sort by length (descending) to match longer phrases first
    for keyword in sorted(VARIABLE_KEYWORDS.keys(), key=len, reverse=True):
        if keyword in query_lower:
            variable = VARIABLE_KEYWORDS[keyword]
            matched_keyword = keyword
            break

    # Extract scenario
    climate_scenario = None
    for pattern, scenario in SCENARIO_PATTERNS.items():
        if re.search(pattern, query_lower):
            climate_scenario = scenario
            break

    # Extract simulation round
    simulation_round = None
    for round_name in SIMULATION_ROUNDS:
        if round_name.lower() in query_lower:
            simulation_round = round_name
            break

    # Auto-detect round for ISIMIP2a-only crops
    if variable and simulation_round is None:
        crop_code = variable.split("-")[1] if "-" in variable else None
        if crop_code in ISIMIP2A_ONLY_CROPS:
            simulation_round = "ISIMIP2a"

    filters = SearchFilters(
        variable=variable,
        climate_scenario=climate_scenario,
        simulation_round=simulation_round,
    )

    # Build explanation
    explanation = "Extracted from keywords: "
    parts = []
    if variable:
        parts.append(f"variable={variable}")
    if climate_scenario:
        parts.append(f"scenario={climate_scenario}")
    if simulation_round:
        parts.append(f"round={simulation_round}")

    explanation += ", ".join(parts) if parts else "no specific filters found"

    return ParsedQuery(
        filters=filters,
        explanation=explanation,
    )


# Backwards compatibility alias
def parse_natural_query(
    query: str,
    api_key: str = "",
    agent_id: str = "",
) -> ParsedQuery:
    """Parse natural language query (compatibility wrapper).

    Args:
        query: Natural language query.
        api_key: Ignored (kept for backwards compatibility).
        agent_id: Ignored (kept for backwards compatibility).

    Returns:
        ParsedQuery with extracted filters.
    """
    return parse_query(query)


# For backwards compatibility with existing imports
class LLMParser:
    """Keyword-based query parser (legacy compatibility class)."""

    def __init__(self, api_key: str = "", agent_id: str = ""):
        """Initialize parser (api_key and agent_id ignored)."""
        pass

    def parse_with_fallback(self, query: str) -> ParsedQuery:
        """Parse query using keyword matching."""
        return parse_query(query)

    def keyword_fallback(self, query: str) -> ParsedQuery:
        """Parse query using keyword matching."""
        return parse_query(query)
