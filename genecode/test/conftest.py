"""
Configure Test
~~~~~~~~~~~~~~
"""

from pathlib import Path
import pytest
import tealogger

from genecode import common


# Configure Logger
_conftest_logger = tealogger.Logger('conftest-logger')
_conftest_logger.setLevel(tealogger.DEBUG)

CURRENT_MODULE_PATH = Path(__file__).parent.expanduser().resolve()


@pytest.fixture(scope='session')
def email() -> str:
    """Email

    Get and return the `email` address.

    :return: The `email` address
    :rtype: str
    """

    data = common.load_json_file(CURRENT_MODULE_PATH.parent / 'metadata.json')

    return data['database']['email']


@pytest.fixture(scope='session')
def api_key() -> str:
    """API Key

    Get and return the `api_key`.

    :return: The `api_key`
    :rtype: str
    """

    data = common.load_json_file(CURRENT_MODULE_PATH.parent / 'metadata.json')

    return data['database']['api_key']


@pytest.fixture(scope='session')
def setup_teardown_database():
    """Setup and Teardown Database

    Prepare and finalize test for `database` functionality.
    """

    # Execute before test

    yield

    # Execute after test
