"""ISIMIP API query functionality."""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any

from isimip_client.client import ISIMIPClient


@dataclass
class SearchFilters:
    """Filters for ISIMIP dataset search.

    Attributes:
        simulation_round: ISIMIP simulation round (e.g., "ISIMIP3b", "ISIMIP2b")
        climate_scenario: Climate scenario (e.g., "ssp370", "rcp85", "historical")
        variable: Variable code (e.g., "led", "burntarea", "potevap")
        climate_forcing: Climate model (e.g., "gfdl-esm4", "ukesm1-0-ll")
        model: Impact model
        timestep: Time resolution (e.g., "annual", "monthly", "daily")
        product: Data product type (e.g., "OutputData", "DerivedOutputData")
    """

    simulation_round: Optional[str] = None
    climate_scenario: Optional[str] = None
    variable: Optional[str] = None
    climate_forcing: Optional[str] = None
    model: Optional[str] = None
    timestep: Optional[str] = None
    product: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert filters to dict, excluding None values."""
        return {k: v for k, v in self.__dict__.items() if v is not None}


@dataclass
class DatasetInfo:
    """Information about an ISIMIP dataset.

    Attributes:
        id: Unique dataset identifier
        name: Dataset filename
        url: Download URL
        simulation_round: ISIMIP simulation round
        climate_scenario: Climate scenario
        variable: Variable code
        model: Impact model name
        timestep: Time resolution
        size: File size in bytes (optional)
    """

    id: str
    name: str
    url: str
    simulation_round: str
    climate_scenario: Optional[str] = None
    variable: Optional[str] = None
    model: Optional[str] = None
    timestep: Optional[str] = None
    size: Optional[int] = None

    @classmethod
    def from_api_response(cls, data: Dict[str, Any]) -> "DatasetInfo":
        """Create DatasetInfo from ISIMIP API response.

        Args:
            data: Dataset dict from API response.

        Returns:
            DatasetInfo instance.
        """
        specifiers = data.get("specifiers", {})

        # Get file URL from the files array (actual downloadable file)
        files = data.get("files", [])
        if files:
            # Use the first file's URL
            file_info = files[0]
            url = file_info.get("file_url", "")
            file_name = file_info.get("name", data.get("name", ""))
            file_size = file_info.get("size", data.get("size"))
            file_id = file_info.get("id", data.get("id", ""))
        else:
            # Fallback to dataset path (less reliable)
            path = data.get("path", "")
            base_url = "https://files.isimip.org"
            url = f"{base_url}/{path}" if path else ""
            file_name = data.get("name", "")
            file_size = data.get("size")
            file_id = data.get("id", "")

        return cls(
            id=file_id,
            name=file_name,
            url=url,
            simulation_round=specifiers.get("simulation_round", ""),
            climate_scenario=specifiers.get("climate_scenario"),
            variable=specifiers.get("variable"),
            model=specifiers.get("model"),
            timestep=specifiers.get("time_step"),  # Note: API uses time_step not timestep
            size=file_size,
        )

    @classmethod
    def from_api_response_all_files(cls, data: Dict[str, Any]) -> List["DatasetInfo"]:
        """Create DatasetInfo for each file in the API response.

        Use this when a dataset contains multiple files and you want all of them.

        Args:
            data: Dataset dict from API response.

        Returns:
            List of DatasetInfo instances, one per file.
        """
        specifiers = data.get("specifiers", {})
        files = data.get("files", [])

        if not files:
            return [cls.from_api_response(data)]

        results = []
        for file_info in files:
            results.append(cls(
                id=file_info.get("id", ""),
                name=file_info.get("name", ""),
                url=file_info.get("file_url", ""),
                simulation_round=specifiers.get("simulation_round", ""),
                climate_scenario=specifiers.get("climate_scenario"),
                variable=specifiers.get("variable"),
                model=specifiers.get("model"),
                timestep=specifiers.get("time_step"),
                size=file_info.get("size"),
            ))

        return results


class ISIMIPQuery:
    """Client for querying the ISIMIP data repository.

    Wraps the isimip-client library to provide a simplified interface
    for searching and retrieving dataset information.
    """

    def __init__(self, timeout: int = 30):
        """Initialize the query client.

        Args:
            timeout: Request timeout in seconds.
        """
        self._client = ISIMIPClient()
        self._timeout = timeout

    def search(self, filters: SearchFilters) -> List[DatasetInfo]:
        """Search for datasets matching the given filters.

        Args:
            filters: SearchFilters with query constraints.

        Returns:
            List of DatasetInfo objects for matching datasets.
        """
        filter_dict = filters.to_dict()
        response = self._client.datasets(**filter_dict)

        results = response.get("results", [])
        return [DatasetInfo.from_api_response(item) for item in results]

    def search_by_query(self, query: str) -> List[DatasetInfo]:
        """Search using a text query string.

        Args:
            query: Free-text search query.

        Returns:
            List of DatasetInfo objects for matching datasets.
        """
        response = self._client.datasets(query=query)
        results = response.get("results", [])
        return [DatasetInfo.from_api_response(item) for item in results]

    def get_files(self, dataset_id: str) -> List[Dict[str, Any]]:
        """Get files for a specific dataset.

        Args:
            dataset_id: Dataset ID to get files for.

        Returns:
            List of file information dicts.
        """
        response = self._client.files(dataset=dataset_id)
        return response.get("results", [])


def search_datasets(**kwargs) -> List[DatasetInfo]:
    """Convenience function to search ISIMIP datasets.

    Args:
        **kwargs: Filter parameters (simulation_round, variable, etc.)

    Returns:
        List of DatasetInfo objects.
    """
    query = ISIMIPQuery()
    filters = SearchFilters(**kwargs)
    return query.search(filters)
