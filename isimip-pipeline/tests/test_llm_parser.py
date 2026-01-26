"""Tests for keyword-based query parsing."""

import pytest

from isimip_pipeline.search.llm_parser import (
    LLMParser,
    ParsedQuery,
    parse_natural_query,
    parse_query,
    VARIABLE_KEYWORDS,
    ISIMIP2A_ONLY_CROPS,
)
from isimip_pipeline.search.isimip_query import SearchFilters


class TestParsedQuery:
    """Test ParsedQuery data structure."""

    def test_parsed_query_with_filters(self):
        """ParsedQuery should hold filters and explanation."""
        query = ParsedQuery(
            filters=SearchFilters(variable="led", simulation_round="ISIMIP3b"),
            explanation="Searching for drought exposure data.",
            raw_response={"filters": {"variable": "led"}},
        )
        assert query.filters.variable == "led"
        assert "drought" in query.explanation.lower()

    def test_parsed_query_to_search_filters(self):
        """ParsedQuery should convert to SearchFilters."""
        query = ParsedQuery(
            filters=SearchFilters(
                variable="burntarea",
                climate_scenario="ssp585",
            ),
            explanation="Fire burnt area under high emissions.",
        )
        filters = query.filters
        assert filters.variable == "burntarea"
        assert filters.climate_scenario == "ssp585"


class TestLLMParser:
    """Test LLMParser compatibility class."""

    def test_parser_initializes(self):
        """LLMParser should initialize (API credentials ignored)."""
        parser = LLMParser(api_key="test-key", agent_id="test-agent")
        assert parser is not None

    def test_parser_initializes_without_credentials(self):
        """LLMParser should work without credentials."""
        parser = LLMParser()
        assert parser is not None

    def test_parser_parse_with_fallback(self):
        """Parser should use keyword matching."""
        parser = LLMParser(api_key="", agent_id="")
        result = parser.parse_with_fallback("drought exposure ssp585")
        assert result.filters.variable == "led"
        assert result.filters.climate_scenario == "ssp585"

    def test_parser_keyword_fallback(self):
        """Parser keyword_fallback should work."""
        parser = LLMParser()
        result = parser.keyword_fallback("wildfire burned area ssp370")
        assert result.filters.variable == "burntarea"
        assert result.filters.climate_scenario == "ssp370"


class TestParseQuery:
    """Test the main parse_query function."""

    def test_extract_variable_drought(self):
        """Should extract drought variable."""
        result = parse_query("drought exposure data")
        assert result.filters.variable == "led"

    def test_extract_variable_fire(self):
        """Should extract fire/burntarea variable."""
        result = parse_query("wildfire burnt area")
        assert result.filters.variable == "burntarea"

    def test_extract_variable_temperature(self):
        """Should extract temperature variable."""
        result = parse_query("air temperature data")
        assert result.filters.variable == "tas"

    def test_extract_variable_agriculture(self):
        """Should extract agriculture variables."""
        result = parse_query("maize yield")
        assert result.filters.variable == "yield-mai"

    def test_extract_variable_barley(self):
        """Should extract barley and auto-detect ISIMIP2a."""
        result = parse_query("barley yield")
        assert result.filters.variable == "yield-bar"
        assert result.filters.simulation_round == "ISIMIP2a"

    def test_extract_variable_fish(self):
        """Should extract marine fishery variable."""
        result = parse_query("fish biomass")
        assert result.filters.variable == "tcb"

    def test_extract_variable_large_fish(self):
        """Should extract large fish size class."""
        result = parse_query("large fish biomass")
        assert result.filters.variable == "b30cm"


class TestScenarioExtraction:
    """Test scenario extraction from queries."""

    def test_extract_ssp585(self):
        """Should extract SSP5-8.5 scenario."""
        result = parse_query("ssp585 future projections")
        assert result.filters.climate_scenario == "ssp585"

    def test_extract_ssp585_with_dash(self):
        """Should extract SSP5-8.5 with various formats."""
        result = parse_query("ssp 5-8.5 data")
        assert result.filters.climate_scenario == "ssp585"

    def test_extract_ssp370(self):
        """Should extract SSP3-7.0 scenario."""
        result = parse_query("middle of the road ssp370")
        assert result.filters.climate_scenario == "ssp370"

    def test_extract_rcp85(self):
        """Should extract RCP8.5 scenario."""
        result = parse_query("rcp85 projections")
        assert result.filters.climate_scenario == "rcp85"

    def test_extract_rcp_with_dot(self):
        """Should extract RCP with decimal format."""
        result = parse_query("rcp 8.5 data")
        assert result.filters.climate_scenario == "rcp85"

    def test_extract_historical(self):
        """Should extract historical scenario."""
        result = parse_query("historical baseline data")
        assert result.filters.climate_scenario == "historical"


class TestSimulationRoundExtraction:
    """Test simulation round extraction from queries."""

    def test_extract_isimip3b(self):
        """Should extract ISIMIP3b round."""
        result = parse_query("ISIMIP3b data")
        assert result.filters.simulation_round == "ISIMIP3b"

    def test_extract_isimip2a(self):
        """Should extract ISIMIP2a round."""
        result = parse_query("isimip2a historical")
        assert result.filters.simulation_round == "ISIMIP2a"

    def test_auto_detect_isimip2a_for_barley(self):
        """Should auto-detect ISIMIP2a for barley."""
        result = parse_query("barley yield projections")
        assert result.filters.variable == "yield-bar"
        assert result.filters.simulation_round == "ISIMIP2a"

    def test_auto_detect_isimip2a_for_cotton(self):
        """Should auto-detect ISIMIP2a for cotton."""
        result = parse_query("cotton yield data")
        assert result.filters.variable == "yield-cot"
        assert result.filters.simulation_round == "ISIMIP2a"

    def test_no_auto_detect_for_common_crops(self):
        """Should not auto-detect round for common crops."""
        result = parse_query("maize yield data")
        assert result.filters.variable == "yield-mai"
        assert result.filters.simulation_round is None


class TestCombinedExtraction:
    """Test extraction of multiple fields."""

    def test_extract_variable_and_scenario(self):
        """Should extract both variable and scenario."""
        result = parse_query("drought exposure ssp585")
        assert result.filters.variable == "led"
        assert result.filters.climate_scenario == "ssp585"

    def test_extract_all_fields(self):
        """Should extract variable, scenario, and round."""
        result = parse_query("ISIMIP3b wildfire burnt area ssp370")
        assert result.filters.variable == "burntarea"
        assert result.filters.climate_scenario == "ssp370"
        assert result.filters.simulation_round == "ISIMIP3b"

    def test_explanation_contains_extracted_info(self):
        """Explanation should list extracted filters."""
        result = parse_query("drought ssp585")
        assert "variable=led" in result.explanation
        assert "scenario=ssp585" in result.explanation


class TestEdgeCases:
    """Test edge cases and special scenarios."""

    def test_empty_query(self):
        """Should handle empty query."""
        result = parse_query("")
        assert result.filters.variable is None
        assert "no specific filters" in result.explanation

    def test_no_matches(self):
        """Should handle query with no matches."""
        result = parse_query("random unrelated text")
        assert result.filters.variable is None
        assert "no specific filters" in result.explanation

    def test_longer_phrase_priority(self):
        """Should match longer phrases before shorter ones."""
        # "drought exposure" should match before "drought"
        result = parse_query("drought exposure")
        assert result.filters.variable == "led"

    def test_case_insensitivity(self):
        """Should be case insensitive."""
        result = parse_query("DROUGHT EXPOSURE SSP585")
        assert result.filters.variable == "led"
        assert result.filters.climate_scenario == "ssp585"


class TestParseNaturalQueryFunction:
    """Test the convenience function."""

    def test_parse_natural_query_works(self):
        """parse_natural_query should work (API params ignored)."""
        result = parse_natural_query(
            "drought data ssp585",
            api_key="ignored",
            agent_id="ignored",
        )
        assert result.filters.variable == "led"
        assert result.filters.climate_scenario == "ssp585"

    def test_parse_natural_query_without_params(self):
        """parse_natural_query should work without optional params."""
        result = parse_natural_query("wildfire burnt area")
        assert result.filters.variable == "burntarea"


class TestVariableKeywordsCoverage:
    """Test that important variables are covered."""

    def test_exposure_metrics_covered(self):
        """Exposure metrics should be in keywords."""
        assert "drought" in VARIABLE_KEYWORDS
        assert "heatwave" in VARIABLE_KEYWORDS
        assert "flood" in VARIABLE_KEYWORDS

    def test_agriculture_covered(self):
        """Agriculture crops should be in keywords."""
        assert "maize" in VARIABLE_KEYWORDS
        assert "wheat" in VARIABLE_KEYWORDS
        assert "barley" in VARIABLE_KEYWORDS
        assert "rice" in VARIABLE_KEYWORDS

    def test_marine_covered(self):
        """Marine fishery should be in keywords."""
        assert "fish" in VARIABLE_KEYWORDS
        assert "large fish" in VARIABLE_KEYWORDS

    def test_isimip2a_only_crops(self):
        """ISIMIP2a-only crops should be tracked."""
        assert "bar" in ISIMIP2A_ONLY_CROPS
        assert "cot" in ISIMIP2A_ONLY_CROPS
        assert "pot" in ISIMIP2A_ONLY_CROPS
