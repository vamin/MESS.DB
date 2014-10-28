# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB calculate module

This module contains the calculate tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
from socket import gethostname

import mess.mapreduce as mapreduce
from mess.tool import AbstractTool
from mess.utils import get_inchikey_dir, is_inchikey, load_method


class Calculate(AbstractTool):
    """This tool runs calculations on molecules in MESS.DB."""

    def __init__(self):
        """Set description of tool."""
        self.description = 'Calculate molecular properties with mapreduce'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('method', type=str, help='A method name')
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help='a list of inchikeys (default: STDIN)')
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
        if args.inchikeys.name == '<stdin>' and args.inchikeys.isatty():
            sys.exit('No input specified.')
        else:
            try:
                inchikeys = set(row.split()[0].strip()
                                for row in args.inchikeys)
            except IndexError:
                inchikeys = set([])
                return
        method = load_method(args.method)
        method.set_parent_path(args.parent_path)
        if args.map:
            args.inchikeys.close()
            self.map_client(method, args.hostname)
            return
        elif args.reduce:
            args.inchikeys.close()
            self.reduce_client(method, args.hostname)
            return
        method.setup()
        if args.mapreduce_server:
            self.mapreduce_server(inchikeys, method)
        else:
            self.mapreduce_local(inchikeys, method)

    def mapreduce_local(self, inchikeys, method):
        """Run a method's map and reduce functions locally."""
        keys = {}
        for inchikey in inchikeys:
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            for key, values in method.map(inchikey,
                                          get_inchikey_dir(inchikey)):
                try:
                    keys[key].append(values)
                except KeyError:
                    keys[key] = [values]
        for key, values in keys.iteritems():
            method.reduce(key, values)

    def mapreduce_server(self, inchikeys, method):
        """Start a mapreduce server."""
        self.log_console.info('hostname is %s' % gethostname())
        datasource = {}
        for inchikey in inchikeys:
            if not is_inchikey(inchikey, enforce_standard=True):
                sys.exit('%s is not a valid InChIKey.' % inchikey)
            datasource[inchikey] = get_inchikey_dir(inchikey)
        server = mapreduce.Server()
        server.datasource = datasource
        server.password = method.hash
        hostfile = os.path.join(os.path.dirname(__file__),
                                '../../temp/%s.host' % server.password)
        with open(hostfile, 'w') as f:
            f.write(gethostname())
        server.run()
        self.log_console.info('all mappers and reducers have finished')
    
    def map_client(self, method, hostname):
        """Start a map client."""
        self.log_console.info(('map client connecting to mapreduce server '
                               'at %s'), hostname)
        client = mapreduce.Client()
        client.password = method.hash
        client.mapfn = method.map
        hostfile = os.path.join(os.path.dirname(__file__),
                                '../../temp/%s.host' % client.password)
        if hostname == 'localhost' and os.path.isfile(hostfile):
            with open(hostfile, 'r') as f:
                hostname = f.readline()
        client.run(hostname, mapreduce.DEFAULT_PORT)
        self.log_console.info('map client done')
    
    def reduce_client(self, method, hostname):
        """Start a reduce client."""
        self.log_console.info(('reduce client connecting to mapreduce server '
                               'at %s'), hostname)
        client = mapreduce.Client()
        client.password = method.hash
        client.reducefn = method.reduce
        hostfile = os.path.join(os.path.dirname(__file__),
                                '../../temp/%s.host' % client.password)
        if hostname == 'localhost' and os.path.isfile(hostfile):
            with open(hostfile, 'r') as f:
                hostname = f.readline()
        client.run(hostname, mapreduce.DEFAULT_PORT)
        self.log_console.info('reduce client done')
    
    
def load():
    """Load Calculate()."""
    return Calculate()
