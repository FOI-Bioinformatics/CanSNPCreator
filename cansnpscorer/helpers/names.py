"""Cleaning of sequence and strain names.

Names reach the tools from FASTA headers, group files and annotation files, and
the same genome is often written differently in each. They are reduced to their
alphanumeric characters so that the three can be matched, which is also why a
cluster name must be compared before it is cleaned: a cluster called A.1 must
stay A.1 to match the hierarchy.
"""

import logging
import re

NOT_ALPHANUMERIC = re.compile("[^A-Za-z0-9]+")


def clean_name(name: str) -> str:
    """
    Reduce a name to its alphanumeric characters.

    Args:
        name: The name as the file gives it.

    Returns:
        The name with every other character removed.

    Example:
        >>> clean_name("GCA_000007125.1 ASM712v1")
        'GCA0000071251ASM712v1'
    """
    cleaned = NOT_ALPHANUMERIC.sub("", name)
    if cleaned != name:
        logging.debug(f"Cleaned the name {name!r} to {cleaned!r}")
    return cleaned
