from __future__ import print_function
from __future__ import unicode_literals

import argparse
import shutil
import os
import sys

from _db import MessDB
from _tools import AbstractTool

### TODO: handle source, method (level, program, method*, parameter,
### property) pruning

class Remove(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = ('Remove molecules from MESS.DB')
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
    
    def execute(self, args):
        """Remove specified elements."""
        db = MessDB()
        c = db.cursor()
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            try:
                inchikey_dir = self.get_inchikey_dir(inchikey)
                shutil.rmtree(inchikey_dir)
                print('%s dir removed\n' % inchikey, file=sys.stderr)
            except OSError:
                print('%s did not have a directory\n' % inchikey, 
                      file=sys.stderr)
            try:
                parent = os.path.relpath(os.path.join(inchikey_dir, '../'))
                os.removedirs(parent)
            except OSError:
                pass
            records = 0
            q = 'DELETE from molecule WHERE inchikey=?'
            c.execute(q, (inchikey,))
            records += c.rowcount
            q = 'DELETE from molecule_synonym WHERE inchikey=?'
            c.execute(q, (inchikey,))
            records += c.rowcount
            q = 'DELETE from molecule_source WHERE inchikey=?'
            c.execute(q, (inchikey,))
            records += c.rowcount
            q = 'DELETE from molecule_state_method_property WHERE inchikey=?'
            c.execute(q, (inchikey,))
            records += c.rowcount
            db.commit()
            print('%i %s records removed from db\n\n' % (records, inchikey),
                  file=sys.stderr)


def load():
    """Load Remove()."""
    return Remove()