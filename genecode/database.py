"""
Database
~~~~~~~~

The `database` module provide functionality for the database.
"""

from Bio import Entrez
from Bio import SeqIO
import tealogger


class Database:
    """Database"""

    def __new__(
        cls,
        *args,
        **kwargs
    ) -> None:
        """New"""
        return super().__new__(cls)

    def __init__(
        self,
        email: str,
    ) -> None:
        """Constructor

        :param email: the email address
        :type email: str
        """

        # Configure Logger
        self._logger = tealogger.Logger('database-logger')
        self._logger.setLevel(tealogger.DEBUG)

        self._email = email

    def search(
        self,
        database: str,
        term: str,
        **kwargs,
    ) -> None:
        """Search

        :param database: the database to search
        :type database: str
        :param term: the query to search
        :type term: str
        """

        # Configure the Entrez email
        Entrez.email = self._email

        # Execute the search and get a `handle` for the result(s)
        handle = Entrez.esearch(
            db=database,
            term=term,
            **kwargs,
        )

        # Execute an initial `read` to parse the result(s)
        record = Entrez.read(handle)
        # Get a `Count` of the result(s)
        record_count = record['Count']
        self._logger.debug(f'Record Count: {record_count}')


# query = '"escherichia coli"[Organism]'

# Entrez.email = 'ytkme@outlook.com'
# handle = Entrez.esearch(db='Nucleotide', term=query, retmax=200, retstart=201)
# record = Entrez.read(handle)
# handle.close()

# print(f'Record Count: {record["Count"]}')

# print(f'Data: {record["IdList"]}')
