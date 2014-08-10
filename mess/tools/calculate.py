from __future__ import print_function
from __future__ import unicode_literals

import argparse
import imp
import sys
from collections import OrderedDict
from distutils.version import LooseVersion

import pybel
#from mess.mincemeatpy import mincemeat

from _db import MessDB
from _path import Path
from _tool import AbstractTool
import mapreduce
from decorators import decorate, UnicodeDecorator
from utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)

import hashlib
import json
import time
import logging

class Calculate(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Calculate properties for a set of molecules'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('method', type=str, help='A method name')
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
        subparser.add_argument('-p', '--parent-path', type=int, default=0,
                               help=('Specify a parent path id. If not set, '
                                     'start from InChI'))
        subparser.add_argument('-s', '--mapreduce', action='store_true',
                               help='Start a mapreduce server')
        subparser.add_argument('-m', '--map', action='store_true',
                               help='Start a mapreduce mapper client')
        subparser.add_argument('-r', '--reduce', action='store_true',
                               help='Start a mapreduce reducer client')
    
    def execute(self, args):
        """Run calculations."""
        method = load_method(args.method, MessDB())
        path = Path(MessDB())
        print(args.map)
        print(args.reduce)
        if (args.map):
            self.map_client(method)
            return
        elif (args.reduce):
            self.reduce_client(method)
            return
        if args.inchikeys.name == '<stdin>' and args.inchikeys.isatty():
            sys.exit('No input specified.')
        method.setup()
        path.setup(method.method_id, args.parent_path)
        if (args.mapreduce):
            self.mapreduce_server(args.inchikeys, method, path)
        else:
            self.mapreduce_local(args.inchikeys, method, path)

    def mapreduce_local(self, inchikeys, method, path):
        count = 0
        for row in inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            queries = dict()
            for query, values in method.map(inchikey, 
                                            [path.path_id, 
                                             path.method_dir, 
                                             path.parent_method_dir]):
                try:
                    queries[query].append(values)
                except KeyError:
                    queries[query] = [values]
            print('mapreduce_local:')
            print(queries)
            for query, values in queries.iteritems():
                count += method.reduce(query, values)
        print(count)

    def mapreduce_server(self, inchikeys, method, path):
        datasource = {}
        for row in inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            datasource[inchikey] = [path.path_id, 
                                    path.method_dir, 
                                    path.parent_method_dir]
        server = mapreduce.Server()
        server.datasource = datasource
        server.password = method.__hash__()
        print(server.password)
        results = server.run(debug=1)
        print(results)
    
    def map_client(self, method):
        while True:
            try:
                client = mapreduce.Client()
                client.password = method.__hash__()
                print(client.password)
                client.mapfn = method.map
                #client.reducefn = method.reduce
                client.run('localhost', mapreduce.DEFAULT_PORT, debug=1)
            except KeyboardInterrupt:
                break
            except:
                exc_info = sys.exc_info()
                logging.exception('%s:%s', exc_info[0], exc_info[1])
                break
            finally:
                break
    
    def reduce_client(self, method):
        client = mapreduce.Client()
        client.password = method.__hash__()
        print(client.password)
        client.reducefn = method.reduce
        client.run('localhost', mapreduce.DEFAULT_PORT, debug=1)
    
    
def load():
    """Load Calculate()."""
    return Calculate()