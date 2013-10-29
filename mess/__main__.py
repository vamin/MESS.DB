#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import os
import sys

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..' )) # add parent
                                                                  # dir to path
from _tools import ToolsManager

def main():
    toolsmanager = ToolsManager()
    parser = argparse.ArgumentParser(
        description="A collection of tools for interacting with the mess db.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)#,
	    #add_help=False)
    toolsmanager.populate_parser(parser)
    args = parser.parse_args()
    tool = toolsmanager.load_tool(args.subparser_name)
    if tool.check_dependencies():
        tool.execute(args)

if __name__ == '__main__':
    main()