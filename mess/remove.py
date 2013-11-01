import argparse
import shutil
import os
import sqlite3
import sys

from _db import MessDB
from _tools import AbstractTool

### TODO: handle source, method (level, program, method*, parameter, 
### property) pruning 

class Remove(AbstractTool):
    def __init__(self):
        self.description = ('Remove molecules from MESS.DB')
        self.epilog = ''
    
    def subparse(self, subparser):
        subparser.add_argument('inchikeys', nargs='?', 
                               type=argparse.FileType('r'), default=sys.stdin, 
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))

    def check_dependencies(self):
        return True
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            try:
                inchikey_dir = self.get_inchikey_dir(inchikey)
                shutil.rmtree(inchikey_dir)
                sys.stderr.write(inchikey + ' dir removed\n')
            except OSError:
                sys.stderr.write(inchikey + ' did not have a directory\n')
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
            #db.commit()
            sys.stderr.write(str(records) + ' ' + inchikey + 
                             ' records removed from db\n')
            sys.stderr.write('\n')



    # dirty! need to set up a separate class for helper code like this
    def get_inchikey_dir(self, inchikey):
        molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
        return os.path.relpath(os.path.join(molecules_dir, inchikey[:1], 
                               inchikey[1:3], inchikey[3:]))

def load():
    # loads the current plugin
    return Remove()