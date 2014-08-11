# -*- coding: utf-8 -*-
"""MESS.DB main module

This module runs MESS.DB. The primary responibility of this code is to load a
user tool and execute it.

Victor Amin 2013-2014
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
from distutils.version import LooseVersion
from signal import signal, SIGPIPE, SIG_DFL

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from _tool import ToolManager
from utils import get_mem_usage

DEBUG_MEM = 0  # set to 1 to print memory info


def check_dependencies():
    """Make sure pybel is availiable, as basically every tool uses it."""
    emsg = ('MESS.DB requires Open Babel (and its python module, pybel) '
            'version >=2.3.0.')
    try:
        import pybel
        version = pybel.ob.OBReleaseVersion()
        if LooseVersion(version) < LooseVersion('2.3.0'):
            sys.exit(emsg + ' Current version is %s.' % version)
    except AttributeError, ImportError:
        sys.exit(emsg)


def main():
    """Parse args and load the user-specified tool and execute it."""
    toolmanager = ToolManager()
    parser = argparse.ArgumentParser(
        description='A collection of tools for interacting with MESS.DB',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    toolmanager.populate_parser(parser)
    args = parser.parse_args()
    check_dependencies()
    tool = toolmanager.load_tool(args.subparser_name)
    tool.execute(args)

if __name__ == '__main__':
    signal(SIGPIPE, SIG_DFL)  # ignore SIG_PIPE, don't throw exception
    if DEBUG_MEM:
        print(get_mem_usage('Memory Start'), file=sys.stderr)
    main()
    if DEBUG_MEM:
        print(get_mem_usage('Memeory Finish'), file=sys.stderr)
