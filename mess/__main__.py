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
import os
import sys
from distutils.version import LooseVersion

MESS_DIR = os.path.dirname(inspect.getfile(inspect.currentframe()))
PARENT_DIR = os.path.realpath(os.path.join(MESS_DIR, '..'))
if PARENT_DIR not in sys.path:
    sys.path.insert(1, PARENT_DIR)

from mess._tool import ToolManager
from mess.decorators import decorate, UnicodeDecorator
from mess.log import Log
from mess.utils import CustomArgparseFormatter, get_mem_usage


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
    log = Log('console')
    log.setup(MESS_DIR, args.verbose)
    log.debug(get_mem_usage())
    tool = toolmanager.load_tool(args.subparser_name)
    tool.execute(args)
    log.debug(get_mem_usage())

if __name__ == '__main__':
    try:
        main()
    except IOError, err:
        if err.errno == 32:
            pass  # broken pipe, fail silently
        else:
            raise
