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
    toolsmanager.populate_parser(parser) # load tool-specific arguments
    # run tool with args
    args = parser.parse_args()
    tool = toolsmanager.load_tool(args.subparser_name)
    tool.execute(args)

def cursor():
    # connect to db
    try:
        mess_db_conn = sqlite3.connect(os.path.join(os.path.dirname(__file__), '../db/mess.db'))
    except IOError:
        sys.exit('could not find mess.db')
    mess_db_conn.row_factory = sqlite3.Row
    return mess_db_conn.cursor()

if __name__ == '__main__':
    main()