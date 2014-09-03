# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB main module

This module runs MESS.DB. The primary responibility of this code is to load a
mess tool and execute it.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import inspect
import logging
import os
import sys
from datetime import date
from distutils.version import LooseVersion

MESS_DIR = os.path.dirname(inspect.getfile(inspect.currentframe()))
PARENT_DIR = os.path.realpath(os.path.join(MESS_DIR, '..'))
if PARENT_DIR not in sys.path:
    sys.path.insert(1, PARENT_DIR)

from mess._tool import ToolManager
from mess.decorators import decorate, GetLoggerDecorator, UnicodeDecorator
from mess.utils import (CustomArgparseFormatter,
                        CustomLoggingFormatter,
                        get_mem_usage)


def check_dependencies():
    """Make sure pybel is availiable, as basically every tool uses it."""
    emsg = ('MESS.DB requires Open Babel (and its python module, pybel) '
            'version >=2.3.0.')
    try:
        import pybel
        decorate(pybel, UnicodeDecorator)
        version = pybel.ob.OBReleaseVersion()
        if LooseVersion(version) < LooseVersion('2.3.0'):
            sys.exit('%s Current version is %s.' % (emsg, version))
    except (AttributeError, ImportError):
        sys.exit(emsg)


def setup_logger(verbose):
    """Setup logging environemnt, including stream and file handlers."""
    logging_formatter = CustomLoggingFormatter(('[%(asctime)s] '
                                                '%(modifiedname)s '
                                                '%(levelname)s: '
                                                '%(message)s'),
                                               '%Y-%m-%d %H:%M:%S')
    logging._defaultFormatter = logging_formatter
    logging.getLogger = GetLoggerDecorator(logging.getLogger)
    # setup main logger and stream handler
    logger = logging.getLogger('mess')
    stream_handler = logging.StreamHandler()
    if verbose:
        logger.setLevel(logging.DEBUG)
        stream_handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    # setup file handlers
    central_log_name = '%s.log' % date.today().strftime('%Y-%m-%d')
    central_log_path = os.path.join(MESS_DIR, '../logs/%s' % central_log_name)
    central_log_handler = logging.FileHandler(central_log_path)
    central_log_handler.setLevel(logging.INFO)
    logger.addHandler(central_log_handler)
    molecule_log_handler = logging.FileHandler('/dev/null')
    molecule_log_handler.setLevel(logging.INFO)
    logger.addHandler(molecule_log_handler)
    return logger


def main():
    """Parse args, load the user-specified tool and execute it."""
    check_dependencies()
    toolmanager = ToolManager()
    parser = argparse.ArgumentParser(
        description='A collection of tools for interacting with MESS.DB',
        formatter_class=CustomArgparseFormatter)
    parser.add_argument('-v', '--verbose', action='store_true',
                        help=('increase output verbosity to include '
                              'debugging messages'))
    toolmanager.populate_parser(parser)
    args = parser.parse_args()
    logger = setup_logger(args.verbose)
    logger.debug(get_mem_usage())
    tool = toolmanager.load_tool(args.subparser_name)
    tool.execute(args)
    logger.debug(get_mem_usage())

if __name__ == '__main__':
    try:
        main()
    except IOError, err:
        if err.errno == 32:
            pass  # broken pipe, fail silently
        else:
            raise
