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
from _sources import Source
from _tools import AbstractTool
from _utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)

class Import(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = ('Import molecule file, multi-molecule file, '
                            'or dir of molecules into MESS.DB')
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('source',
                               help='A molecule source file or directory')
        subparser.add_argument('-s', '--skip-cir', action='store_true',
                               help=('Do not use CIR web service to import '
                                     'IUPAC names and other synonyms'))
    
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
        """Run import method for every molecule in source."""
        m = load_method('import', MessDB())
        s = Source(MessDB())
        s.setup(args.source)
        p = Path(MessDB())
        p.setup(m.method_id)
        pybel.ob.obErrorLog.StopLogging()
        for f in s.files():
            if (f.split('.')[-1] == 'sql' or
                f.split('.')[-1] == 'txt' or
                f[-1] == '~'):
                continue
            for mol in pybel.readfile(f.split('.')[-1],
                                      os.path.join(s.source_dir, f)):
                decorate(mol, UnicodeDecorator)
                pybel.ob.obErrorLog.ClearLog()
                pybel.ob.obErrorLog.StartLogging()
                inchikey = mol.write('inchikey').rstrip()
                pybel.ob.obErrorLog.StopLogging()
                cansmi = mol.write('can').split()[0]
                frag_count = cansmi.count('.') + 1
                for f in cansmi.split('.'):
                    method_args = {}
                    if (frag_count > 1):
                        frag = pybel.readstring('can', f)
                        decorate(frag, UnicodeDecorator)
                        # neutralize fragments
                        #if (frag.charge != 0):
                        #    for atom in frag.atoms:
                        #        atom.OBAtom.SetFormalCharge(0)
                        pybel.ob.obErrorLog.ClearLog()
                        pybel.ob.obErrorLog.StartLogging()
                        frag_inchikey = frag.write('inchikey').rstrip()
                        pybel.ob.obErrorLog.StopLogging()
                        if not is_inchikey(frag_inchikey):
                            print("'%s' is not an importable molecule.\n" % f, 
                                  file=sys.stderr)
                            continue
                        method_args['parent'] = ('from: %s' % 
                                                 unicode(mol.title,
                                                         'utf-8', 'replace'))
                        method_args['inchikey'] = frag_inchikey
                        method_args['mol'] = frag
                    else:
                        if not is_inchikey(inchikey):
                            print("'%s' is not an importable molecule.\n" % f, 
                                  file=sys.stderr)
                            continue
                        method_args['inchikey'] = inchikey
                        method_args['mol'] = mol
                    method_args['source'] = s
                    method_args['path'] = p
                    method_args['skip_cir'] = args.skip_cir
                    status = m.execute(method_args)
                    print('%s: %s' % (method_args['inchikey'], status))


def load():
    """Load Import()."""
    return Import()