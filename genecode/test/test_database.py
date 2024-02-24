"""
Test Database
~~~~~~~~~~~~~
"""

from pathlib import Path
import pytest
import tealogger

from genecode import common
from genecode.database import Database


class TestDatabase:
    """Test Database"""

    # def __init__(self) -> None:
    #     """Constructor"""
    #     # Initialize Logger
    #     self._test_database_logger = None

    # def setup_method(self) -> None:
    #     """Setup Method"""
    #     # Configure Logger
    #     self._test_database_logger = tealogger.Logger('test-database-logger')
    #     self._test_database_logger.setLevel(tealogger.DEBUG)

    @pytest.mark.asyncio
    async def test_search(
        self,
        email: str,  # Fixture
        api_key: str  # Fixture
    ):
        """Test Search

        Test the `search` functionality of the `Database` class with
        default argument.

        :param email: The `email` address fixture
        :type email: str
        :param api_key: The `api_key` fixture
        :type api_key: str
        """

        # Get the search `term` (query)
        query_path = Path(__file__).parents[1] / 'query.txt'
        term = common.load_text_file(path=query_path)[0]

        # Create a `Database` instance
        database = Database(
            email=email,
            api_key=api_key,
        )

        # Execute
        result_list = await database.search(
            database='Nucleotide',
            term=term,
        )

        assert result_list and isinstance(result_list, list)

    @pytest.mark.asyncio
    async def test_search_retmax(
        self,
        email: str,  # Fixture
        api_key: str  # Fixture
    ):
        """Test Search Retmax

        Test the `search` functionality of the `Database` class with
        the `retmax` argument.

        :param email: The `email` address fixture
        :type email: str
        :param api_key: The `api_key` fixture
        :type api_key: str
        """
