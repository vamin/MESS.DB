from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
from socket import gethostname

import pybel

import mapreduce
from _db import MessDB
from _path import Path
from _tool import AbstractTool
from decorators import decorate, UnicodeDecorator
from utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)


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
        subparser.add_argument('-s', '--mapreduce-server', action='store_true',
                               help='Start a mapreduce server')
        subparser.add_argument('-m', '--map', action='store_true',
                               help='Start a mapreduce mapper client')
        subparser.add_argument('-r', '--reduce', action='store_true',
                               help='Start a mapreduce reducer client')
        subparser.add_argument('-n', '--hostname', type=str,
                               default='localhost',
                               help=('Set mapreduce server hostname for '
                                     'client, defaults to \'localhost\''))
    
    def execute(self, args):
        """Run calculations."""
        method = load_method(args.method, MessDB())
        if args.map:
            self.map_client(method, args.hostname)
            return
        elif args.reduce:
            self.reduce_client(method, args.hostname)
            return
        if args.inchikeys.name == '<stdin>' and args.inchikeys.isatty():
            sys.exit('No input specified.')
        method.setup()
        path = Path(MessDB())
        path.setup(method.method_id, args.parent_path)
        if args.mapreduce_server:
            self.mapreduce_server(args.inchikeys, method, path)
        else:
            self.mapreduce_local(args.inchikeys, method, path)

    def mapreduce_local(self, inchikeys, method, path):
        keys = {}
        for row in inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            for key, values in method.map(inchikey,
                                          [path.path_id,
                                           path.method_dir,
                                           path.parent_method_dir]):
                try:
                    keys[key].append(values)
                except:
                    keys[key] = [values] 
        for key, values in keys.iteritems():
            method.reduce(key, values)

    def mapreduce_server(self, inchikeys, method, path):
        print('Hostname is: %s' % gethostname(), file=sys.stderr)
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
        server.password = hash(method)
        results = server.run(debug=0)
        print('Done!', file=sys.stderr)
    
    def map_client(self, method, hostname):
        print('Connecting to mapreduce server at: %s' % hostname,
              file=sys.stderr)
        client = mapreduce.Client()
        client.password = hash(method)
        client.mapfn = method.map
        client.run(hostname, mapreduce.DEFAULT_PORT, debug=0)
        print('Done!', file=sys.stderr)
    
    def reduce_client(self, method, hostname):
        print('Connecting to mapreduce server at: %s' % hostname,
              file=sys.stderr)
        client = mapreduce.Client()
        client.password = hash(method)
        client.reducefn = method.reduce
        client.run(hostname, mapreduce.DEFAULT_PORT, debug=0)
        print('Done!', file=sys.stderr)
    
    
def load():
    """Load Calculate()."""
    return Calculate()
