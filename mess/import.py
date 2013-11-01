import argparse
import imp
import os
import sqlite3
import subprocess
import sys
from distutils.version import LooseVersion

import pybel

from _db import MessDB
from _paths import Path
from _tools import AbstractTool
from _sources import Source

class Import(AbstractTool):
    def __init__(self):
        self.description = ('Import molecule file, multi-molecule file, '
                            'or dir of molecules into MESS.DB')
        self.epilog = ''
            
    def subparse(self, subparser):
        subparser.add_argument('source', 
                               help='A molecule source file or directory.')

    def check_dependencies(self):
        try:
            babel = subprocess.Popen(['babel', '-V'], 
                                     stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
            babel_version = babel.stdout.read().split()[2]
            if (LooseVersion(babel_version) < LooseVersion('2.3.0')):
                sys.exit(('This tool requires Open Babel (and its python '
                         'module, pybel) version >=2.3.0.'))
        except OSError:
            sys.exit(('This tool requires Open Babel (and its python module, '
                      'pybel) version >=2.3.0.'))
        return True
    
    def execute(self, args):
        self.db = MessDB()
        self.c = self.db.cursor()
        # setup import method
        try:
            method = imp.load_source('method', 
                                     os.path.join(os.path.dirname( __file__ ), 
                                     '..', 'methods', 'import', 'method.py'))
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
            if (f.split('.')[-1] == 'sql' or 
                f.split('.')[-1] == 'txt' or
                f[-1] == '~'):
                continue
            for mol in pybel.readfile(f.split('.')[-1], 
                                      os.path.join(args.source, f)):
                pybel.ob.obErrorLog.StartLogging()
                inchikey = mol.write('inchikey').rstrip()
                pybel.ob.obErrorLog.StopLogging()
                cansmi = mol.write('can').split()[0]
                frag_count = cansmi.count('.') + 1
                for f in cansmi.split('.'):
                    method_args = {}
                    if (frag_count > 1):
                        frag = pybel.readstring('can', f)
                        # neutralize fragments
                        #if (frag.charge != 0):
                        #    for atom in frag.atoms:
                        #        atom.OBAtom.SetFormalCharge(0)
                        pybel.ob.obErrorLog.ClearLog()
                        pybel.ob.obErrorLog.StartLogging()
                        frag_inchikey = frag.write('inchikey').rstrip()
                        if not s.is_inchikey(frag_inchikey):
                            sys.stderr.write("'" + 
                                             f + 
                                             "' is not an importable " +
                                             'molecule.\n')
                            continue
                        pybel.ob.obErrorLog.StopLogging()
                        method_args['parent'] = 'from:' + mol.title
                        method_args['inchikey'] = frag_inchikey
                        method_args['mol'] = frag
                    else:
                        if not s.is_inchikey(inchikey):
                            sys.stderr.write("'" + 
                                             f + 
                                             "' is not an importable " +
                                             'molecule.\n')
                            continue
                        method_args['inchikey'] = inchikey
                        method_args['mol'] = mol
                    method_args['source'] = s
                    method_args['path'] = p
                    status = m.execute(method_args)
                    print method_args['inchikey'] + ': ' + status


def load():
    # loads the current plugin
    return Import()