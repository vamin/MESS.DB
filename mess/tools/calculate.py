from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
from socket import gethostname

import pybel

import mess.mapreduce as mapreduce
from mess._db import MessDB
from mess._path import MethodPath
from mess._tool import AbstractTool
from mess.decorators import decorate, UnicodeDecorator
from mess.utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)


class Calculate(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Calculate molecular properties with mapreduce'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('method', type=str, help='A method name')
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
        subparser.add_argument('-p', '--parent-path', type=int, default=None,
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
        path = MethodPath()
        path.setup_path(method.method_id, args.parent_path)
        map_args = {'path_id': path.get_path_id(),
                    'method_dir': path.get_path_directory(),
                    'parent_method_dir': path.get_parent_path_directory()}
        if args.mapreduce_server:
            self.mapreduce_server(args.inchikeys, method, map_args)
        else:
            self.mapreduce_local(args.inchikeys, method, map_args)

    def mapreduce_local(self, inchikeys, method, map_args):
        keys = {}
        for row in inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            for key, values in method.map(inchikey, map_args):
                try:
                    keys[key].append(values)
                except KeyError:
                    keys[key] = [values]
        for key, values in keys.iteritems():
            method.reduce(key, values)

    def mapreduce_server(self, inchikeys, method, map_args):
        print('Hostname is: %s' % gethostname(), file=sys.stderr)
        datasource = {}
        for row in inchikeys:
            inchikey = row.split()[0].strip()
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            datasource[inchikey] = map_args
        server = mapreduce.Server()
        server.datasource = datasource
        server.password = method.hash
        results = server.run(debug=0)
        print('Done!', file=sys.stderr)
    
    def map_client(self, method, hostname):
        print('Connecting to mapreduce server at: %s' % hostname,
              file=sys.stderr)
        client = mapreduce.Client()
        client.password = method.hash
        client.mapfn = method.map
        client.run(hostname, mapreduce.DEFAULT_PORT, debug=0)
        print('Done!', file=sys.stderr)
    
    def reduce_client(self, method, hostname):
        print('Connecting to mapreduce server at: %s' % hostname,
              file=sys.stderr)
        client = mapreduce.Client()
        client.password = method.hash
        client.reducefn = method.reduce
        client.run(hostname, mapreduce.DEFAULT_PORT, debug=0)
        print('Done!', file=sys.stderr)
    
    
def load():
    """Load Calculate()."""
    return Calculate()
