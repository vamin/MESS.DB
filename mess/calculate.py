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
from _utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)

class Calculate(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Calculate properties for a set of molecules'
        self.epilog = ''
    
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
    
    def check_dependencies(self):
        """Check for dependencies (Open Babel >=2.3).
        
        Returns:
            True if all dependencies are met.
        
        """
        try:
            if (LooseVersion(pybel.ob.OBReleaseVersion()) <
                LooseVersion('2.3.0')):
                sys.exit(('This tool requires Open Babel (and its python '
                          'module, pybel) version >=2.3.0.'))
        except AttributeError, OSError:
            sys.exit(('This tool requires Open Babel (and its python module, '
                      'pybel) version >=2.3.0.'))
        return True
    
    def execute(self, args):
        """Run calculations."""
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


def load():
    """Load Calculate()."""
    return Calculate()