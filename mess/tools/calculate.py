from __future__ import print_function
from __future__ import unicode_literals

import argparse
import imp
import sys
from distutils.version import LooseVersion

import pybel
from mess.mincemeatpy import mincemeat

from _db import MessDB
from _path import Path
from _tool import AbstractTool
from decorators import decorate, UnicodeDecorator
from utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)

import hashlib
import json
import time

class Calculate(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Calculate properties for a set of molecules'
        self.epilog = ''
        self.mince = 0
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
        subparser.add_argument('-m', '--method', type=str, required=True,
                               help='A method name')
        subparser.add_argument('-p', '--parent-path', type=int, default=0,
                               help=('Specify a parent path id. If not set, '
                                     'start from InChI'))
        subparser.add_argument('-mr', '--mapreduce', action='store_true',
                               help='Start a mapreduce server')
    
    def execute(self, args):
        """Run calculations."""
        if (args.mapreduce):
            self.execute_mapreduce(args)
            return 0
        m = load_method(args.method, MessDB())
        p = Path(MessDB())
        p.setup(m.method_id, args.parent_path)
        pybel.ob.obErrorLog.StopLogging()
        method_args = {}
        method_args['path'] = p
        for row in args.inchikeys:
            method_args['inchikey'] = row.split()[0].strip()
            if not is_inchikey(method_args['inchikey'], enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % 
                         method_args['inchikey'])
            status = m.execute(method_args)
            print('%s: %s' % (method_args['inchikey'], status))

    def execute_mapreduce(self, args):
        m = load_method(args.method, MessDB())
        p = Path(MessDB())
        p.setup(m.method_id, args.parent_path)
        s = mincemeat.Server()
        source = {}
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            source[inchikey] = [p.path_id, p.method_dir, p.parent_method_dir]
        s.datasource = source
        s.mapfn = Calculate.map
        s.reducefn = Calculate.reduce
        pwd = hashlib.sha1(json.dumps(time.time(), sort_keys=True)).hexdigest()
        results = s.run_server(password=pwd)
        print(results)

    @staticmethod
    def map(k, v):
        inchikey = k
        (path_id, method_dir, parent_method_dir) = v
        for query, values in m.execute(inchikey, path_id, method_dir, parent_method_dir):           yield query, values

    @staticmethod
    def reduce(k, vs):
        #query = k
        #db = MessDB()
        #db.cursor().executemany(query, vs)
        print(k)
        print(vs)
        return k, 1

def load():
    """Load Calculate()."""
    return Calculate()