#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import csv
import os
import re
import sqlite3
import sys

from _db import MessDB
from _tools import AbstractTool

class Select(AbstractTool):
    def __init__(self):
        self.description = 'returns a list of molecules based on SQL query on mess.db'
    
    def subparse(self, subparser):
        subparser.add_argument("sql", help='an SQL statement or file that returns inchikeys in first column')
        subparser.add_argument("-s", "--subset", type=str, help="subset the output by first letter(s) of inchikey")
        subparser.add_argument("-r", "--regex-subset", type=str, help="subset the output by regex on inchikey (superceded by subset)")
        subparser.add_argument("-d", "--delimiter", type=str, default='\t', help="choose a delimiter for output files, tab is default")
        subparser.add_argument("-hd", "--headers", action='store_true', help="include headers in output, not recommended if piping to 'mess calculate'")

    def check_dependencies(self):
        return True
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        try:
            c.execute(open(args.sql).read())
        except sqlite3.OperationalError:
            sys.exit('"' + args.sql + '" does not contain valid sql.')
        except IOError:
            try:
                c.execute(args.sql)
            except sqlite3.OperationalError:
                sys.exit('"' + args.sql + '" is neither valid sql nor a path to a file containing valid sql.')
        # check that sql returns inchikey in first column
        if not (c.description[0][0].lower() == 'inchikey'):
            sys.exit('SQL must return inchikey in first column.')
        # print inchikeys
        writer = csv.writer(sys.stdout, delimiter = args.delimiter)
        if (args.headers):
            writer.writerow(list(h[0] for h in c.description))
        for row in c:
            if (not args.subset and not args.regex_subset) \
                or (args.subset and row[0].startswith(args.subset)) \
                or (args.regex_subset and re.match(args.regex_subset, row[0], re.IGNORECASE)):
                writer.writerow(list(self.xstr(v).decode('utf-8') for v in row))

    def xstr(self, s):
        if s is None:
            return ''
        return str(s)


def load():
    # loads the current plugin
    return Select()
