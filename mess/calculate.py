from __future__ import print_function
from __future__ import unicode_literals

import argparse
import imp
import os
import sqlite3
import sys
from distutils.version import LooseVersion

import pybel

from _db import MessDB
from _decorators import decorate, UnicodeDecorator
from _paths import Path
from _tools import AbstractTool
from _utils import load_method

decorate(pybel, UnicodeDecorator)

class Calculate(AbstractTool):
    def __init__(self):
        self.description = 'Calculate properties for a set of molecules'
        self.epilog = ''
    
    def subparse(self, subparser):
        subparser.add_argument('inchikeys', nargs='?', 
                               type=argparse.FileType('r'), default=sys.stdin, 
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
        subparser.add_argument('-m', '--method', required=True, 
                               help='A method name.')
        subparser.add_argument('-p', '--parent-path', default='', 
                               help=('Specify a parent path id. If not set, '
                                     'start from InChI.'))
    
    def check_dependencies(self):
        try:
            if (LooseVersion(pybel.ob.OBReleaseVersion()) < 
                LooseVersion('2.3.0')):
                sys.exit(('This tool requires Open Babel (and its python '
                          'module, pybel) version >=2.3.0.'))
        except OSError:
            sys.exit(('This tool requires Open Babel (and its python module, '
                      'pybel) version >=2.3.0.'))
        return True
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        # setup import method
        m = load_method(args.method, db)
        # setup path
        p = Path(db)
        p.setup(m.method_id, args.parent_path)
        pybel.ob.obErrorLog.StopLogging() 
        method_args = {}
        method_args['path'] = p
        for row in args.inchikeys:
            method_args['inchikey'] = row.split()[0].strip()
            status = m.execute(method_args)
            print(method_args['inchikey'] + ': ' + status)
        

def load():
    # loads the current plugin
    return Calculate()