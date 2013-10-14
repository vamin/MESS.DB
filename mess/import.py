#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import imp
import os
import sqlite3
import sys
#from cStringIO import StringIO

import pybel

from _db import MessDB
from _paths import Path
from _tools import AbstractTool
from _sources import Source

class Import(AbstractTool):
    def __init__(self):
        self.description = 'imports molecule file, multi-molecule file, or dir of molecules into mess.db'
            
    def subparse(self, subparser):
        subparser.add_argument("source", help='A molecule source file or directory.')
    
    def execute(self, args):
        self.db = MessDB()
        self.c = self.db.cursor()
        # setup import method
        try:
            method = imp.load_source('method', os.path.join(os.path.dirname( __file__ ), '..', 'methods', 'import', 'method.py'))
        except IOError:
            sys.exit("Can't find 'import' method in the 'methods' directory.")
        m = method.Method(self.db)
        m.setup()
        # setup source
        s = Source(self.db)
        s.setup(args.source)
        # setup path
        p = Path(self.db)
        p.setup(m.get_method_id())
        # turn off pybel logging
        pybel.ob.obErrorLog.StopLogging() 
        for f in s.files():
            if (f.split('.')[-1] == 'sql'):
                continue
            for mol in pybel.readfile(f.split('.')[-1], os.path.join(args.source, f)):
                #ob_logs = []
                # generate inchikey
                pybel.ob.obErrorLog.StartLogging()
                inchikey = mol.write('inchikey').rstrip()
                pybel.ob.obErrorLog.StopLogging()
                method_args = {}
                method_args['inchikey'] = inchikey
                method_args['mol'] = mol
                method_args['source'] = s
                method_args['path'] = p
                status = m.execute(method_args)
                print method_args['inchikey'] + ': ' + status


def load():
    # loads the current plugin
    return Import()