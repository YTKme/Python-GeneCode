"""
Database
~~~~~~~~

The `database` module provide functionality for the database.
"""

from Bio import Entrez
from Bio import SeqIO
import tealogger

from genecode.constant import (
    DEFAULT_RETMAX,
)


class Database:
    """Database"""

    # def __new__(
    #     cls,
    #     *args,
    #     **kwargs
    # ) -> None:
    #     """New"""
    #     return super().__new__(cls)

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
    ) -> list:
        """Search

        :param database: The `database` to search
        :type database: str
        :param term: The `term` (query) to search
        :type term: str
        :return: The list of `IdList` from the search result
        :rtype: list
        """

        # Configure the Entrez email
        Entrez.email = self._email
        Entrez.api_key = self._api_key

        try:
            # Execute the search and get a `handle` for the result(s)
            handle = Entrez.esearch(
                db=database,
                term=term,
                **kwargs,
            )

            # Execute an initial `read` to parse the result(s)
            record = Entrez.read(handle)
            # Get the `Count` of the result(s), convert to `int`
            record_count = int(record['Count'])
            self._logger.debug('Record Count: %s', record_count)
            # Get the `RetMax` of the result(s), convert to `int`
            record_retmax = int(record['RetMax'])
            self._logger.debug('Record RetMax: %s', record_retmax)

            # self._logger.debug(f'Record: {record}')

            # Parse the `IdList` as current `result_list`
            result_list = record['IdList']

            # Compare record `Count` (total) to `RetMax` (retrieved)
            if record_retmax < record_count:
                # If the `RetMax` is less than the `Count`
                record_remain = record_count - record_retmax
                self._logger.debug('Record Remain: %s', record_remain)

                # Update the `retmax` with the remaining record(s)
                # - Use as `kwargs` for the next retrieve
                # - Account for previous `kwargs` function argument
                kwargs['retmax'] = record_remain

                # Need to retrieve more record(s)
                handle = Entrez.esearch(
                    db=database,
                    term=term,
                    **kwargs,
                )
                record = Entrez.read(handle)

                # Combine the new `IdList` with current `result_list`
                result_list += record['IdList']

            # Return the list of result(s)
            self._logger.debug('Result List: %s', result_list)
            return result_list
        except Exception as e:
            # TODO:
            pass
        finally:
            # Clean up by flush and close the `handle` stream
            handle.close()


    async def download(
        self,
    ) -> None:
        """Download"""
        pass
