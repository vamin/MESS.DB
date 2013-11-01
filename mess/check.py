import argparse
import glob
import os
import sqlite3
import sys

from _db import MessDB
from _tools import AbstractTool

class Check(AbstractTool):
    def __init__(self):
        self.description = ('Check integrity of mess.db and concordance with '
                            'molecules dir')
        self.epilog = ''
    
    def subparse(self, subparser):
        pass

    def check_dependencies(self):
        return True
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        self.check_dir_structure()
        self.check_db_structure()
        self.check_db_dir_inchikey_concordance(c)
        self.check_db_dir_method_concordance(c)
    
    def check_db_dir_inchikey_concordance(self, c):
        # get list of inchikeys from db:
        db_inchikeys = set()
        c.execute('SELECT inchikey FROM molecule;')
        for row in c:
            db_inchikeys.add(row['inchikey'])
        # get list of inchikeys from molecules/ dir
        dir_inchikeys = set()
        inchis = glob.glob(os.path.join(os.path.dirname(__file__), 
                           '../molecules/*/*/*/', '*.inchi'))
        for i in inchis:
            s = i.split('/')[-1]
            dir_inchikeys.add(s.split('.')[0])
        # compare inchikeys from db vs dir
        in_db_not_dir = db_inchikeys - dir_inchikeys
        in_dir_not_db = dir_inchikeys - db_inchikeys
        print str(len(in_db_not_dir)) + (' InChIKeys in mess.db that are not '
                                         'in molecules dir:')
        print "\n".join(i for i in in_db_not_dir)
        print str(len(in_dir_not_db)) + (' InChIKeys in molecules dir that are '
                                         'not in mess.db:')
        print "\n".join(i for i in in_dir_not_db)
    
    def check_db_dir_method_concordance(self, c):
        pass
    
    def check_dir_structure(self):
        pass
    
    def check_db_structure(self):
        pass

def load():
    # loads the current plugin
    return Check()