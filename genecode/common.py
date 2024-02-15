"""
Common
~~~~~~

The `common` module provide common functionality.
"""

import csv
import json
from os import PathLike
import random
import time
import yaml

from tealogger import logger


# Configure Logger
common_logger = logger.Logger('common-logger')
common_logger.setLevel(logger.DEBUG)


def fuzz(is_fuzz: bool = False) -> None:
    """Random Fuzz

    Fuzz(ing) is an automated software testing technique that involves
    providing invalid, unexpected, or random data as inputs to a
    computer program.

    Fuzz(ing) is also a technique for amplifying race condition errors
    to make them more visible.

    :param is_fuzz: Whether or not to enable fuzzing, defaults to
        `False`
    :type is_fuzz: bool, optional
    """
    if is_fuzz:
        time.sleep(random.random())


def load_text_file(path: PathLike) -> list[str]:
    """Load a text file.

    :param path: The path to the text file
    :type path: PathLike
    :returns: The text data as a list of string
    :rtype: list[str]
    """
    common_logger.info('Loading text file %s', path)

    with open(path, 'r', encoding='utf-8') as file:
        data = file.readlines()

    common_logger.info('Done loading text file %s', path)
    common_logger.debug('Text data: %s', data)

    return data


def dump_text_file(path: PathLike, data: list[str]) -> None:
    """Dump a text file.

    :param path: The path to write the text file
    :type path: PathLike
    :param data: The text data as a list of string
    :type data: list[str]
    """
    ...


def load_csv_file(path: PathLike) -> list[dict]:
    """Load a CSV file.

    :param path: The path to the CSV file
    :type path: PathLike
    :returns: The CSV data as a list of dictionary
    :rtype: list[dict]
    """

    common_logger.info('Loading CSV file %s', path)

    with open(path, 'r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        data = list(csv_reader)

    common_logger.info('Done loading CSV file %s', path)
    common_logger.debug('CSV data: %s', data)

    return data


def dump_csv_file(path: PathLike, data: list[dict]) -> None:
    """Dump a CSV file.

    :param path: The path to write the CSV file
    :type path: PathLike
    :param data: The CSV data as a list of dictionary
    :type data: list[dict]
    """

    common_logger.info(f'Dumping CSV file {path!r}')

    # Assume consistent `keys` for `data`
    header = data[0].keys()

    with open(path, 'w', encoding='utf-8') as file:
        csv_writer = csv.DictWriter(file, fieldnames=header)
        csv_writer.writeheader()
        csv_writer.writerows(data)

    common_logger.info(f'Done dumping CSV file {path!r}')


def load_json_file(path: PathLike) -> dict:
    """Load a JSON file.

    :param path: The path to the JSON file
    :type path: PathLike
    :returns: The JSON data as a dictionary
    :rtype: dict
    """

    common_logger.info(f'Loading JSON file {path!r}')

    with open(path, 'r', encoding='utf-8') as file:
        data = json.load(file)

    common_logger.info(f'Done loading JSON file {path!r}')
    common_logger.debug(f'JSON data: {data!r}')

    return data


def dump_json_file(path: PathLike, data: dict) -> None:
    """Dump a JSON file.

    :param path: The path to write the JSON file
    :type path: PathLike
    :param data: The JSON data as a dictionary
    :type data: dict
    """

    common_logger.info(f'Dumping JSON file {path!r}')

    with open(path, 'w', encoding='utf-8') as file:
        json.dump(data, file)

    common_logger.info(f'Done dumping JSON file {path!r}')


def load_yaml_file(path: PathLike) -> dict:
    """Load a YAML file.

    :param path: The path to the YAML file
    :type path: PathLike
    :returns: The YAML data as a dictionary
    :rtype: dict
    """

    common_logger.info(f'Loading YAML file {path!r}')

    with open(path, 'r', encoding='utf-8') as file:
        data = yaml.safe_load(file)

    common_logger.info(f'Done loading YAML file {path!r}')
    common_logger.debug(f'YAML data: {data!r}')

    return data


def dump_yaml_file(path: PathLike, data: dict) -> None:
    """Dump a YAML file.

    :param path: The path to write the YAML file
    :type path: PathLike
    :param data: The YAML data as a dictionary
    :type data: dict
    """

    common_logger.info(f'Dumping YAML file {path!r}')

    with open(path, 'w', encoding='utf-8') as file:
        yaml.dump(data, file)

    common_logger.info(f'Done dumping YAML file {path!r}')
