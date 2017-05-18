# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This docstring will be used to print a description of the main program.
"""
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The line above will help with 2to3 support.


def main():
    """Entry-point for command-line scripts.

    Returns
    -------
    :class:`int`
        Exit status that will be passed to :func:`sys.exit`.
    """
    from sys import argv
    from os.path import basename
    from argparse import ArgumentParser
    #
    # Parse arguments
    #
    executable = basename(argv[0])
    parser = ArgumentParser(description=__doc__, prog=executable)
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Print extra information.')
    options = parser.parse_args()
    #
    #
    #
    print('Hello World!')
    if options.verbose:
        print('Verbose selected!')
    return 0
