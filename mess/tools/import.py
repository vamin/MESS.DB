from __future__ import print_function
from __future__ import unicode_literals

import os
import sys

import pybel

from mess._db import MessDB
from mess._path import MethodPath
from mess._source import Source
from mess._tool import AbstractTool
from mess.decorators import decorate, UnicodeDecorator
from mess.utils import is_inchikey, load_method

decorate(pybel, UnicodeDecorator)


class Import(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Import molecules into MESS.DB'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('source',
                               help='A molecule source file or directory')
    
    def execute(self, args):
        """Run import method for every molecule in source."""
        imp = load_method('_import', MessDB())
        imp.setup()
        source = Source(MessDB())
        source.setup(args.source)
        path = MethodPath()
        path.setup_path(imp.method_id)
        path_id = path.get_path_id()
        pybel.ob.obErrorLog.StopLogging()
        keys = {}
        for file_ in source.files():
            for mol in pybel.readfile(file_.split('.')[-1],
                                      os.path.join(source.source_dir, file_)):
                decorate(mol, UnicodeDecorator)
                pybel.ob.obErrorLog.ClearLog()
                pybel.ob.obErrorLog.StartLogging()
                pybel.ob.obErrorLog.StopLogging()
                cansmi = mol.write('can').split()[0]
                for fragment in cansmi.split('.'):
                    method_args = {}
                    if cansmi.count('.') > 0:
                        frag = pybel.readstring('can', fragment)
                        decorate(frag, UnicodeDecorator)
                        pybel.ob.obErrorLog.ClearLog()
                        pybel.ob.obErrorLog.StartLogging()
                        inchikey = frag.write('inchikey').rstrip()
                        pybel.ob.obErrorLog.StopLogging()
                        if not is_inchikey(inchikey):
                            print(("'%s' is not "
                                   'an importable molecule.\n') % fragment,
                                  file=sys.stderr)
                            continue
                        method_args['parent'] = ('from: %s' %
                                                 unicode(mol.title,
                                                         'utf-8', 'replace'))
                        method_args['mol'] = frag
                    else:
                        inchikey = mol.write('inchikey').rstrip()
                        if not is_inchikey(inchikey):
                            print(("'%s' is not "
                                   'an importable molecule.\n') % fragment,
                                  file=sys.stderr)
                            continue
                        method_args['mol'] = mol
                    method_args['source'] = source
                    method_args['path_id'] = path_id
                    imp.execute(inchikey, method_args)


def load():
    """Load Import()."""
    return Import()
