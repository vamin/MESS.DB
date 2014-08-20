from __future__ import print_function
from __future__ import unicode_literals

import codecs
import csv
import re
import string
import sys

import pybel

from mess._db import MessDB
from mess._tool import AbstractTool
from mess.utils import xstr
from mess.tools.match import Match


class Select(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Select a list of molecules based on properties'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('-q', '--query', type=str,
                               help=('An SQL query string or file that '
                                     'returns inchikeys in first column'))
        subparser.add_argument('-s', '--smarts', type=str,
                               help=('Subset using SMARTS pattern or file '
                                     'listing a series of SMARTS patterns'))
        subparser.add_argument('-n', '--property-name', type=str,
                               help='Name of numeric propery to filter by')
        subparser.add_argument('-o', '--property-operator', type=str,
                               choices=['=', '!=', '>', '<', '>=', '<='],
                               help='How to compare property')
        subparser.add_argument('-v', '--property-value', type=float,
                               help='What value of property to compare to')
        subparser.add_argument('-p', '--path', type=int, default=None,
                               help=('Specify a path id to restrict to a '
                                     'particular calculation.'))
        subparser.add_argument('-t', '--part', type=int,
                               help='Subset, --part n --of N subsets')
        subparser.add_argument('-f', '--of', type=int,
                               help='Subset, --part n --of N subsets')
        subparser.add_argument('-r', '--regex-subset', type=str,
                               help=('Subset the output by regex on inchikey '
                                     '(cumulative with -p/-f)'))
        subparser.add_argument('-d', '--delimiter', type=str, default='\t',
                               help=('Choose a delimiter for output files, '
                                     'tab is default'))
        subparser.add_argument('-hd', '--headers', action='store_true',
                               help=('Include headers in output, not '
                                     'recommended if piping to '
                                     "'mess calculate'"))
    
    def execute(self, args):
        """Run select query, output table."""
        if args.query and (args.property_name
                           or args.property_operator
                           or args.property_value):
            sys.exit(('Custom SQL queries are mutually exclusive with '
                      'property filtering.'))
        if (args.part or args.of) and not (args.part and args.of):
            sys.exit(('If you specify a --part n, you must also specify --of '
                      'N (e.g. something like --part 1 --of 5).'))
        if args.part and args.of:
            if args.part > args.of:
                sys.exit('--part must be smaller than --of.')
            if args.part < 1:
                sys.exit('--part must be >=1.')
            alpha = string.ascii_uppercase
            alpha3 = [''.join([a, b, c]) for a in alpha
                                         for b in alpha
                                         for c in alpha]  # AAA to ZZZ
            if args.of > len(alpha3):
                sys.exit(('MESS.DB does not support subsetting into more than '
                          '%i parts.' % len(alpha3)))
            subsets = [alpha3[i::args.of] for i in xrange(args.of)]
            subset = subsets[args.part - 1]
        db = MessDB()
        cur = db.cursor()
        if args.query:
            try:
                cur.execute(codecs.open(args.query,
                                                 encoding='utf-8').read())
            except sqlite3.OperationalError:
                sys.exit("'%s' does not contain valid sql." % args.query)
            except IOError:
                try:
                    cur.execute(args.query)
                except sqlite3.OperationalError:
                    sys.exit(("'%s' is neither valid sql nor a path "
                              'to a file containing valid sql.') % args.query)
        elif (args.property_name and args.property_operator and
              args.property_value is not None):
            query, values = self.property_query(args.property_name,
                                                args.property_operator,
                                                args.property_value,
                                                args.path)
            cur.execute(query, values)
        else:
            cur.execute('SELECT inchikey FROM molecule')
        # check that sql returns inchikey in first column
        if not cur.description[0][0].lower() == 'inchikey':
            sys.exit('Query must return inchikey in first column.')
        # print table
        writer = csv.writer(sys.stdout, delimiter=args.delimiter)
        if args.headers:
            writer.writerow(list(h[0] for h in cur.description))
        for result in cur:
            if args.regex_subset and not re.match(args.regex_subset, result[0],
                                                  re.IGNORECASE):
                continue
            if args.part and args.of:
                if not any(result[0].startswith(a) for a in subset):
                    continue
            if args.smarts:
                matches = 0
                query = 'SELECT inchi FROM molecule WHERE inchikey = ?'
                inchi = db.execute(query, (result[0],)).fetchone()[0]
                mol = pybel.readstring('inchi', 'InChI=%s' % inchi)
                for (smarts_obj,
                     smarts_str) in Match.smarts_generator(args.smarts):
                    matches += len(smarts_obj.findall(mol))
                if not matches:
                    continue
            writer.writerow(list(xstr(v).decode('utf-8') for v in result))
        db.close()  # must be closed manually to prevent db locking during pipe
    
    @classmethod
    def property_query(cls, prop, operator, value, path=None):
        query = ('SELECT inchikey, name, result, units '
                 'FROM molecule_method_property_denorm '
                 'WHERE name = ? AND result %s ?') % operator
        if path is not None:
            query += ' AND method_path_id = ?'
            return query, (prop, value, path)
        else:
            return query, (prop, value)


def load():
    """Load Select()."""
    return Select()
