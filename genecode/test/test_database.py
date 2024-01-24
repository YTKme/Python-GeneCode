"""
Test Database
~~~~~~~~~~~~~
"""

from genecode.database import Database


class TestDatabase:
    """Test Database"""

    def test_search(self):
        """Test Search"""

        database = Database(
            email='ytkme@outlook.com',
        )

        database.search(
            database='Nucleotide',
            term='"escherichia coli"[Organism]',
        )
