"""
CLI (Command Line Interface)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CLI (Command Line Interface) functionality.
"""

import argparse


def parse_argument():
    """Parse Argument

    :returns: the namespace with the parsed arguments
    :rtype: namespace
    """

    # Parent Parser
    parent_parser = argparse.ArgumentParser(
        description='Gene Code is a Python package and application for genetic development.',
        add_help=False,
    )
    parent_parser.add_argument(
        '--database',
        action='store_true',
        required=False,
        help='Execute `database` functionality.',
    )
    parent_parser.add_argument(
        '--alignment',
        action='store_true',
        required=False,
        help='Execute `alignment` functionality.',
    )

    parent_argument, _ = parent_parser.parse_known_args()

    # Child Parser
    child_parser = argparse.ArgumentParser(
        parents=[parent_parser],
    )

    # Database Argument
    child_parser.add_argument(
        '-e',
        '--email',
        type=str,
        action='store',  # Store the value
        dest='email',  # Destination to store
        required=parent_argument.database,
        help='The `email` of the user.'
    )

    # Alignment Argument
    child_parser.add_argument(
        '-q',
        '--query',
        type=str,
        action='store',  # Store the value
        dest='query',  # Destination to store
        required=parent_argument.alignment,
        help='The `query` to execute.'
    )

    return child_parser.parse_args()


def main():
    """Main"""
    argument = parse_argument()
    print(argument)


if __name__ == '__main__':
    main()
