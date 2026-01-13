"""LLM-based natural language query parsing for ISIMIP searches."""

import json
import re
from dataclasses import dataclass, field
from typing import Optional, Dict, Any

import httpx

from isimip_pipeline.search.isimip_query import SearchFilters


# Known ISIMIP variable mappings for keyword fallback
# These are checked before calling the LLM API
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
    # Agriculture
    "maize": "yield-mai-noirr",
    "maize yield": "yield-mai-noirr",
    "wheat": "yield-whe-noirr",
    "wheat yield": "yield-whe-noirr",
    "rice": "yield-ric-noirr",
    "rice yield": "yield-ric-noirr",
    "soybean": "yield-soy-noirr",
    "crop yield": "yield-mai-noirr",
    # Vegetation/Carbon
    "gpp": "gpp",
    "gross primary": "gpp",
    "npp": "npp",
    "net primary": "npp",
    "carbon": "cveg",
    "vegetation carbon": "cveg",
    "soil carbon": "csoil",
}

# Known scenario patterns
SCENARIO_PATTERNS = {
    r"ssp\d{3}": lambda m: m.group(0),  # ssp126, ssp370, ssp585
    r"rcp\d{2}": lambda m: m.group(0),  # rcp26, rcp60, rcp85
    r"historical": lambda m: "historical",
}

# Known simulation rounds
SIMULATION_ROUNDS = ["ISIMIP2a", "ISIMIP2b", "ISIMIP3a", "ISIMIP3b"]


@dataclass
class ParsedQuery:
    """Result of parsing a natural language query.

    Attributes:
        filters: SearchFilters extracted from the query.
        explanation: Human-readable explanation of what was found.
        raw_response: Raw response from LLM (if used).
    """

    filters: SearchFilters
    explanation: str = ""
    raw_response: Optional[Dict[str, Any]] = None


class LLMParser:
    """Parses natural language queries using LLM API.

    Uses you.com Agent API to convert natural language queries
    into structured ISIMIP search filters.
    """

    API_URL = "https://api.you.com/v1/agents/runs"

    SYSTEM_PROMPT = """You are an ISIMIP dataset search assistant. Convert natural language queries into structured ISIMIP API parameters.

Available parameters:
- simulation_round: ISIMIP2a, ISIMIP2b, ISIMIP3a, ISIMIP3b
- climate_scenario: historical, picontrol, rcp26, rcp60, rcp85, ssp126, ssp370, ssp585
- variable: Common codes include:
  - led: Land area fraction exposed to drought
  - leh: Land area fraction exposed to heatwave
  - lew: Land area fraction exposed to wildfire
  - ler: Land area fraction exposed to river flood
  - lec: Land area fraction exposed to crop failure
  - burntarea: Fire burnt area fraction
  - potevap: Potential evapotranspiration
- climate_forcing: gfdl-esm2m, hadgem2-es, ipsl-cm5a-lr, miroc5, gfdl-esm4, ukesm1-0-ll
- timestep: daily, monthly, annual
- product: InputData, OutputData, DerivedOutputData

Respond with ONLY valid JSON in this format:
{
  "filters": {
    "variable": "...",
    "simulation_round": "...",
    "climate_scenario": "..."
  },
  "explanation": "Brief description of what these datasets contain"
}

Only include filters that are clearly specified or implied in the query."""

    def __init__(self, api_key: str, agent_id: str):
        """Initialize the parser.

        Args:
            api_key: you.com API key.
            agent_id: you.com Agent ID.
        """
        self.api_key = api_key
        self.agent_id = agent_id

    def build_prompt(self, query: str) -> str:
        """Build the prompt for the LLM.

        Args:
            query: Natural language query.

        Returns:
            Full prompt including system instructions.
        """
        return f"{self.SYSTEM_PROMPT}\n\nUser query: {query}\n\nRespond with JSON:"

    async def parse_async(self, query: str) -> ParsedQuery:
        """Parse query using LLM API asynchronously.

        Args:
            query: Natural language query.

        Returns:
            ParsedQuery with extracted filters.
        """
        if not self.api_key or not self.agent_id:
            return self.keyword_fallback(query)

        prompt = self.build_prompt(query)

        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(
                    self.API_URL,
                    headers={
                        "Authorization": f"Bearer {self.api_key}",
                        "Content-Type": "application/json",
                    },
                    json={
                        "agent": self.agent_id,
                        "input": prompt,
                        "stream": False,
                    },
                    timeout=120.0,
                )

                if response.status_code == 200:
                    data = response.json()
                    response_text = data.get("output", "{}")
                    return self.parse_response(response_text)

            except (httpx.RequestError, json.JSONDecodeError):
                pass

        return self.keyword_fallback(query)

    def parse_sync(self, query: str) -> ParsedQuery:
        """Parse query using LLM API synchronously.

        Args:
            query: Natural language query.

        Returns:
            ParsedQuery with extracted filters.
        """
        import asyncio
        return asyncio.run(self.parse_async(query))

    def parse_with_fallback(self, query: str) -> ParsedQuery:
        """Parse query with keyword fallback on failure.

        Args:
            query: Natural language query.

        Returns:
            ParsedQuery with extracted filters.
        """
        try:
            return self.parse_sync(query)
        except Exception:
            return self.keyword_fallback(query)

    def parse_response(self, response_text: str) -> ParsedQuery:
        """Parse LLM response text into ParsedQuery.

        Args:
            response_text: JSON response from LLM.

        Returns:
            ParsedQuery with extracted filters.
        """
        try:
            # Try to extract JSON from response
            data = json.loads(response_text)
        except json.JSONDecodeError:
            # Try to find JSON in the response
            json_match = re.search(r'\{[^{}]*\}', response_text)
            if json_match:
                try:
                    data = json.loads(json_match.group(0))
                except json.JSONDecodeError:
                    data = {}
            else:
                data = {}

        filters_data = data.get("filters", {})
        explanation = data.get("explanation", "")

        filters = SearchFilters(
            simulation_round=filters_data.get("simulation_round"),
            climate_scenario=filters_data.get("climate_scenario"),
            variable=filters_data.get("variable"),
            climate_forcing=filters_data.get("climate_forcing"),
            model=filters_data.get("model"),
            timestep=filters_data.get("timestep"),
            product=filters_data.get("product"),
        )

        return ParsedQuery(
            filters=filters,
            explanation=explanation,
            raw_response=data,
        )

    def keyword_fallback(self, query: str) -> ParsedQuery:
        """Extract filters using keyword matching.

        Args:
            query: Natural language query.

        Returns:
            ParsedQuery with extracted filters.
        """
        query_lower = query.lower()

        # Extract variable
        variable = None
        for keyword, var_code in VARIABLE_KEYWORDS.items():
            if keyword in query_lower:
                variable = var_code
                break

        # Extract scenario
        climate_scenario = None
        for pattern, extractor in SCENARIO_PATTERNS.items():
            match = re.search(pattern, query_lower)
            if match:
                climate_scenario = extractor(match)
                break

        # Extract simulation round
        simulation_round = None
        for round_name in SIMULATION_ROUNDS:
            if round_name.lower() in query_lower:
                simulation_round = round_name
                break

        filters = SearchFilters(
            variable=variable,
            climate_scenario=climate_scenario,
            simulation_round=simulation_round,
        )

        explanation = f"Extracted from keywords: "
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


def parse_natural_query(
    query: str,
    api_key: str = "",
    agent_id: str = "",
) -> ParsedQuery:
    """Convenience function to parse natural language query.

    Args:
        query: Natural language query.
        api_key: you.com API key.
        agent_id: you.com Agent ID.

    Returns:
        ParsedQuery with extracted filters.
    """
    parser = LLMParser(api_key=api_key, agent_id=agent_id)
    return parser.parse_with_fallback(query)
