from __future__ import print_function
from __future__ import unicode_literals

import codecs
import csv
import os
import re
import sqlite3
import string
import sys

from _db import MessDB
from _tools import AbstractTool
from _utils import xstr

class Select(AbstractTool):
    def __init__(self):
        self.description = ('Select a list of molecules '
                            'based on SQL query')
        self.epilog = ''
    
    def subparse(self, subparser):
        subparser.add_argument('sql', nargs='?', 
                               default='SELECT inchikey FROM molecule',
                               help=('an SQL statement or file that returns '
                                     'inchikeys in first column'))
        subparser.add_argument('-p', '--part', type=int, 
                               help='subset, --part n --of N subsets')
        subparser.add_argument('-f', '--of', type=int, 
                               help='subset, --part n --of N subsets')
        subparser.add_argument('-r', '--regex-subset', type=str, 
                               help=('subset the output by regex on inchikey '
                                     '(cumulative with -p/-f)'))
        subparser.add_argument('-d', '--delimiter', type=str, default='\t', 
                               help=('choose a delimiter for output files, '
                                     'tab is default'))
        subparser.add_argument('-hd', '--headers', action='store_true', 
                               help=('include headers in output, not '
                                     'recommended if piping to '
                                     "'mess calculate'"))

    def check_dependencies(self):
        return True
    
    def execute(self, args):
        if (args.part or args.of) and not (args.part and args.of):
            sys.exit(('If you specify a --part n, you must also specify --of '
                      'N (e.g. something like --part 1 --of 5).'))
        if args.part and args.of:
            if args.part > args.of:
                sys.exit('--part must be smaller than --of.')
            if args.part < 1:
                sys.exit('--part must be >=1.')
            alpha = string.ascii_uppercase
            alpha3 = [''.join([a,b,c]) for a in alpha 
                                       for b in alpha 
                                       for c in alpha] # AAA to ZZZ
            if args.of > len(alpha3):
                sys.exit(('MESS.DB does not support subsetting into more than '
                          '%i parts.' % len(alpha3)))
            subsets = [ alpha3[i::args.of] for i in xrange(args.of) ]
            subset = subsets[args.part-1]
        db = MessDB(isolation_level='DEFERRED')
        c = db.cursor()
        try:
            c.execute(codecs.open(args.sql, encoding='utf-8').read())
        except sqlite3.OperationalError:
            sys.exit("'%s' does not contain valid sql." % args.sql)
        except IOError:
            try:
                c.execute(args.sql)
            except sqlite3.OperationalError:
                sys.exit(("'%s' is neither valid sql nor a path "
                          'to a file containing valid sql.') % args.sql)
        # check that sql returns inchikey in first column
        if not c.description[0][0].lower() == 'inchikey':
            sys.exit('Query must return inchikey in first column.')
        # print table
        writer = csv.writer(sys.stdout, delimiter=args.delimiter)
        if args.headers:
            writer.writerow(list(h[0] for h in c.description))
        for r in c:
            if args.regex_subset and not re.match(args.regex_subset, r[0], 
                                                  re.IGNORECASE):
                continue
            if args.part and args.of:
                if not any(r[0].startswith(a) for a in subset):
                    continue
            writer.writerow(list(xstr(v).decode('utf-8') for v in r))
        db.close() # must be closed manually to prevent db locking during pipe


def load():
    # loads the current plugin
    return Select()