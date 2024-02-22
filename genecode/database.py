"""
Database
~~~~~~~~

The `database` module provide functionality for the database.
"""

from Bio import Entrez
from Bio import SeqIO
import tealogger

from genecode.constant import (
    RETRIEVE_MAXIMUM,
)


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
        api_key: str = None,
    ) -> None:
        """Constructor

        :param email: The `email` address of the account
        :type email: str
        :param api_key: The `api_key` of the account, defaults to `None`
        :type api_key: str, optional
        """

        # Configure Logger
        self._logger = tealogger.Logger('database-logger')
        self._logger.setLevel(tealogger.DEBUG)

        self._email = email
        self._api_key = api_key

    async def search(
        self,
        database: str,
        term: str,
        **kwargs,
    ) -> None:
        """Search

        :param database: The `database` to search
        :type database: str
        :param term: The `term` (query) to search
        :type term: str
        """

        # Configure the Entrez email
        Entrez.email = self._email
        Entrez.api_key = self._api_key

        # Execute the search and get a `handle` for the result(s)
        handle = Entrez.esearch(
            db=database,
            term=term,
            **kwargs,
        )

        # Execute an initial `read` to parse the result(s)
        record = Entrez.read(handle)
        # Get a `Count` of the result(s), convert to `int`
        record_count = int(record['Count'])
        self._logger.debug('Record Count: %s', record_count)
        self._logger.debug(f'Record: {record}')

        # if (kwargs.get('retmax') and kwargs.get('retmax') > record_count) or (constant.RETRIEVE_MAX > record_count):
        #     ...

    async def _query(
        self,
        database: str,
        term: str,
        retrieve_start: int = 0,
        retrieve_max: int = RETRIEVE_MAXIMUM,
    ) -> list[str]:
        """Query

        :param database: The `database` to query
        :type database: str
        :param term: The `term` (query) to search
        :type term: str
        :param retrieve_start: The sequential index of the first UID in
            the retrieved set in output, defaults to `0`
        :type retrieve_start: int, optional
        :param retrieve_max: The total number of UIDs from the retrieved
            set in output, defaults to `RETRIEVE_MAX`
        :type retrieve_max: int, optional
        :return: The list of UID(s)
        :rtype: list[str]
        """
        # Configure the Entrez email
        Entrez.email = self._email

        # Execute the search and get a `handle` for the result(s)
        handle = Entrez.esearch(
            db=database,
            term=term,
            retstart=retrieve_start,
            retmax=retrieve_max,
        )

        # Execute an initial `read` to parse the result(s)
        record = Entrez.read(handle)

        return record['IdList']


# query = '"escherichia coli"[Organism]'

# Entrez.email = 'ytkme@outlook.com'
# handle = Entrez.esearch(db='Nucleotide', term=query, retmax=200, retstart=201)
# record = Entrez.read(handle)
# handle.close()

# print(f'Record Count: {record["Count"]}')

# print(f'Data: {record["IdList"]}')
