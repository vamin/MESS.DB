#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..' )) # add parent
                                                                  # dir to path
from _tools import ToolsManager

def main():
    """Parse args and load the specified tool.
    
    """
    toolsmanager = ToolsManager()
    parser = argparse.ArgumentParser(
        description='A collection of tools for interacting with MESS.DB',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    toolsmanager.populate_parser(parser)
    args = parser.parse_args()
    tool = toolsmanager.load_tool(args.subparser_name)
    if tool.check_dependencies():
        tool.execute(args)

if __name__ == '__main__':
    signal(SIGPIPE, SIG_DFL) # ignore SIG_PIPE, don't throw exception about it
    main()