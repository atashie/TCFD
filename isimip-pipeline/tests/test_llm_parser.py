"""Tests for LLM-based query parsing."""

import pytest
from unittest.mock import Mock, patch, AsyncMock
import json

from isimip_pipeline.search.llm_parser import (
    LLMParser,
    ParsedQuery,
    parse_natural_query,
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
    """Test LLMParser class."""

    def test_parser_initializes(self):
        """LLMParser should initialize with API credentials."""
        parser = LLMParser(api_key="test-key", agent_id="test-agent")
        assert parser is not None
        assert parser.api_key == "test-key"
        assert parser.agent_id == "test-agent"

    def test_parser_builds_prompt(self):
        """Parser should build prompt with query."""
        parser = LLMParser(api_key="test-key", agent_id="test-agent")
        prompt = parser.build_prompt("drought exposure metrics")
        assert "drought" in prompt
        assert "ISIMIP" in prompt or "variable" in prompt.lower()

    @pytest.mark.asyncio
    @patch("httpx.AsyncClient.post")
    async def test_parser_calls_api(self, mock_post):
        """Parser should call you.com API."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "output": json.dumps({
                "filters": {"variable": "led", "simulation_round": "ISIMIP3b"},
                "explanation": "Drought exposure data",
            })
        }
        mock_post.return_value = mock_response

        parser = LLMParser(api_key="test-key", agent_id="test-agent")
        result = await parser.parse_async("drought metrics")

        assert result.filters.variable == "led"

    def test_parser_handles_fallback(self):
        """Parser should fall back to keyword extraction on API failure."""
        parser = LLMParser(api_key="", agent_id="")  # No valid credentials

        # Should attempt keyword fallback
        result = parser.parse_with_fallback("wildfire burned area ssp370")

        # Should extract some keywords even without API
        assert result is not None


class TestLLMResponseParsing:
    """Test parsing of LLM responses."""

    def test_parse_json_response(self):
        """Parser should extract filters from JSON response."""
        parser = LLMParser(api_key="test-key", agent_id="test-agent")

        response_text = json.dumps({
            "filters": {
                "variable": "burntarea",
                "simulation_round": "ISIMIP2b",
                "climate_scenario": "rcp85",
            },
            "explanation": "Wildfire burnt area data under high emissions scenario.",
        })

        parsed = parser.parse_response(response_text)

        assert parsed.filters.variable == "burntarea"
        assert parsed.filters.simulation_round == "ISIMIP2b"
        assert parsed.filters.climate_scenario == "rcp85"

    def test_parse_handles_missing_filters(self):
        """Parser should handle response with missing filters."""
        parser = LLMParser(api_key="test-key", agent_id="test-agent")

        response_text = json.dumps({
            "explanation": "Could not determine filters.",
        })

        parsed = parser.parse_response(response_text)
        assert parsed.filters is not None  # Should return empty filters


class TestKeywordFallback:
    """Test keyword-based fallback parsing."""

    def test_extract_variable_from_keywords(self):
        """Should extract known variables from query."""
        parser = LLMParser(api_key="", agent_id="")

        result = parser.keyword_fallback("drought exposure data")
        assert result.filters.variable == "led"

    def test_extract_scenario_from_keywords(self):
        """Should extract scenario from query."""
        parser = LLMParser(api_key="", agent_id="")

        result = parser.keyword_fallback("ssp585 future projections")
        assert result.filters.climate_scenario == "ssp585"

    def test_extract_simulation_round_from_keywords(self):
        """Should extract simulation round from query."""
        parser = LLMParser(api_key="", agent_id="")

        result = parser.keyword_fallback("ISIMIP3b data")
        assert result.filters.simulation_round == "ISIMIP3b"


class TestParseNaturalQueryFunction:
    """Test the convenience function."""

    @patch("isimip_pipeline.search.llm_parser.LLMParser")
    def test_parse_natural_query_creates_parser(self, mock_parser_class):
        """parse_natural_query should create LLMParser."""
        mock_parser = Mock()
        mock_parser.parse_with_fallback.return_value = ParsedQuery(
            filters=SearchFilters(),
            explanation="Test",
        )
        mock_parser_class.return_value = mock_parser

        result = parse_natural_query(
            "drought data",
            api_key="test-key",
            agent_id="test-agent",
        )

        mock_parser_class.assert_called_once_with(
            api_key="test-key",
            agent_id="test-agent",
        )
