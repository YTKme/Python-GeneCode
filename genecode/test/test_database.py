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

    def setup_method(self) -> None:
        """Setup Method"""
        # Configure Logger
        self._test_database_logger = tealogger.Logger('test-database-logger')
        self._test_database_logger.setLevel(tealogger.DEBUG)

    @pytest.mark.asyncio
    async def test_search(self, email: str):
        """Test Search

        :param email: The `email` address fixture
        :type email: str
        """

        query_path = Path(__file__).parents[1] / 'query.txt'
        term = common.load_text_file(path=query_path)[0]

        database = Database(
            email=email,
        )

        await database.search(
            database='Nucleotide',
            term=term,
        )
