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

mess_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
parent_dir = os.path.realpath(os.path.join(mess_dir, '..'))
if parent_dir not in sys.path:
    sys.path.insert(1, parent_dir)

from mess._tool import ToolManager
from mess.utils import get_mem_usage

DEBUG_MEM = 0  # set to 1 to print memory info


def check_dependencies():
    """Make sure pybel is availiable, as basically every tool uses it."""
    emsg = ('MESS.DB requires Open Babel (and its python module, pybel) '
            'version >=2.3.0.')
    try:
        import pybel
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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    toolmanager.populate_parser(parser)
    args = parser.parse_args()
    tool = toolmanager.load_tool(args.subparser_name)
    tool.execute(args)

if __name__ == '__main__':
    try:
        if DEBUG_MEM:
            print(get_mem_usage('Memory Start'), file=sys.stderr)
        main()
        if DEBUG_MEM:
            print(get_mem_usage('Memeory Finish'), file=sys.stderr)
    except IOError, err:
        if err.errno == 32:
            pass  # broken pipe, fail silently
        else:
            raise
