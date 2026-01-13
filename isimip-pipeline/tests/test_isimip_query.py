"""Tests for ISIMIP API query functionality."""

import pytest
from unittest.mock import Mock, patch

from isimip_pipeline.search.isimip_query import (
    ISIMIPQuery,
    DatasetInfo,
    SearchFilters,
    search_datasets,
)


class TestSearchFilters:
    """Test SearchFilters data structure."""

    def test_create_empty_filters(self):
        """Should be able to create filters with no constraints."""
        filters = SearchFilters()
        assert filters.simulation_round is None
        assert filters.climate_scenario is None
        assert filters.variable is None

    def test_create_filters_with_values(self):
        """Should be able to create filters with specific values."""
        filters = SearchFilters(
            simulation_round="ISIMIP3b",
            climate_scenario="ssp370",
            variable="led",
        )
        assert filters.simulation_round == "ISIMIP3b"
        assert filters.climate_scenario == "ssp370"
        assert filters.variable == "led"

    def test_filters_to_dict(self):
        """Filters should convert to dict for API calls."""
        filters = SearchFilters(
            simulation_round="ISIMIP3b",
            variable="burntarea",
        )
        d = filters.to_dict()
        assert d["simulation_round"] == "ISIMIP3b"
        assert d["variable"] == "burntarea"
        # None values should be excluded
        assert "climate_scenario" not in d


class TestDatasetInfo:
    """Test DatasetInfo data structure."""

    def test_dataset_info_has_required_fields(self):
        """DatasetInfo should have essential metadata."""
        info = DatasetInfo(
            id="dataset-123",
            name="test_dataset.nc",
            url="https://files.isimip.org/test.nc",
            simulation_round="ISIMIP3b",
            climate_scenario="ssp370",
            variable="led",
            model="gfdl-esm4",
            timestep="annual",
        )
        assert info.id == "dataset-123"
        assert info.name == "test_dataset.nc"
        assert info.url == "https://files.isimip.org/test.nc"
        assert info.simulation_round == "ISIMIP3b"
        assert info.climate_scenario == "ssp370"

    def test_from_api_response(self):
        """DatasetInfo should parse API response correctly."""
        api_data = {
            "id": "abc123",
            "name": "test.nc",
            "path": "ISIMIP3b/Output/test.nc",
            "specifiers": {
                "simulation_round": "ISIMIP3b",
                "climate_scenario": "ssp370",
                "variable": "led",
                "model": "gfdl-esm4",
                "timestep": "annual",
            },
            "size": 1024,
        }
        info = DatasetInfo.from_api_response(api_data)

        assert info.id == "abc123"
        assert info.name == "test.nc"
        assert info.url == "https://files.isimip.org/ISIMIP3b/Output/test.nc"
        assert info.simulation_round == "ISIMIP3b"
        assert info.climate_scenario == "ssp370"
        assert info.variable == "led"
        assert info.size == 1024


class TestISIMIPQuery:
    """Test ISIMIPQuery client wrapper."""

    def test_query_client_initializes(self):
        """Query client should initialize without errors."""
        query = ISIMIPQuery()
        assert query is not None
        assert hasattr(query, "_client")

    def test_query_has_search_method(self):
        """Query should have search method."""
        query = ISIMIPQuery()
        assert hasattr(query, "search")
        assert callable(query.search)

    def test_query_has_search_by_query_method(self):
        """Query should have search_by_query method for text search."""
        query = ISIMIPQuery()
        assert hasattr(query, "search_by_query")
        assert callable(query.search_by_query)

    def test_search_with_mocked_client(self):
        """Search should call client and return DatasetInfo list."""
        query = ISIMIPQuery()

        # Mock the internal client
        mock_client = Mock()
        mock_client.datasets.return_value = {
            "results": [
                {
                    "id": "abc123",
                    "name": "test.nc",
                    "path": "ISIMIP3b/Output/test.nc",
                    "specifiers": {
                        "simulation_round": "ISIMIP3b",
                        "climate_scenario": "ssp370",
                        "variable": "led",
                        "model": "gfdl-esm4",
                        "timestep": "annual",
                    },
                }
            ],
            "count": 1,
        }
        query._client = mock_client

        filters = SearchFilters(simulation_round="ISIMIP3b", variable="led")
        results = query.search(filters)

        mock_client.datasets.assert_called_once_with(
            simulation_round="ISIMIP3b", variable="led"
        )
        assert len(results) == 1
        assert isinstance(results[0], DatasetInfo)
        assert results[0].id == "abc123"


class TestSearchDatasetsFunction:
    """Test the convenience search_datasets function."""

    @patch("isimip_pipeline.search.isimip_query.ISIMIPQuery")
    def test_search_datasets_with_kwargs(self, mock_query_class):
        """search_datasets should accept keyword arguments."""
        mock_query = Mock()
        mock_query.search.return_value = []
        mock_query_class.return_value = mock_query

        results = search_datasets(
            simulation_round="ISIMIP3b",
            variable="burntarea",
        )

        assert results == []
        mock_query.search.assert_called_once()
