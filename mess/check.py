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
        self.c = db.cursor()
        self.c.execute('SELECT inchikey FROM molecule')
        self.db_inchikeys = set()
        # check that inchikeys are all valid
        for row in self.c:
            if self.is_inchikey(row['inchikey']):
                self.db_inchikeys.add(row['inchikey'])
            else:
                print row['inchikey'] + ' is not a valid InChiKey!'
        self.check_dir_structure()
        self.check_db_structure()
        self.check_db_dir_inchikey_concordance()
    
    def check_db_dir_inchikey_concordance(self):
        # get list of inchikeys from molecules/ dir
        dir_inchikeys = set()
        inchis = glob.glob(os.path.join(os.path.dirname(__file__), 
                           '../molecules/*/*/*/', '*.inchi'))
        for i in inchis:
            s = i.split('/')[-1]
            dir_inchikeys.add(s.split('.')[0])
        # compare inchikeys from db vs dir
        in_db_not_dir = self.db_inchikeys - dir_inchikeys
        in_dir_not_db = dir_inchikeys - self.db_inchikeys
        print str(len(in_db_not_dir)) + (' InChIKeys in mess.db that are not '
                                         'in molecules dir:')
        print "\n".join(i for i in in_db_not_dir)
        print str(len(in_dir_not_db)) + (' InChIKeys in molecules dir that '
                                         'are not in mess.db:')
        print "\n".join(i for i in in_dir_not_db)
    
    def check_dir_structure(self):
        pass
    
    def check_db_structure(self):
        # check inchikey foreign keys
        self.c.execute('SELECT DISTINCT inchikey FROM molecule_synonym')
        self.c.execute('SELECT DISTINCT inchikey FROM molecule_source')
        self.c.execute('SELECT DISTINCT inchikey FROM molecule_method_property')
        # check that sources in db exist in sources dir
        self.c.execute('SELECT source_id, dirname FROM source')
        # check source foreign keys
        self.c.execute('SELECT DISTINCT source_id FROM molecule_source')
        # check level foreign keys
        self.c.execute('SELECT level_id FROM level')
        self.c.execute('SELECT DISTINCT level_id FROM method')
        # check program foreign keys
        self.c.execute('SELECT program_id FROM program')
        self.c.execute('SELECT DISTINCT program_id FROM method')
        # check parameter foreign keys
        self.c.execute('SELECT parameter_id FROM parameter')
        self.c.execute('SELECT DISTINCT parameter_id FROM method_parameter')
        # check property foreign keys
        self.c.execute('SELECT property_id FROM property')
        self.c.execute(('SELECT DISTINCT property_id '
                        'FROM molecule_method_property'))
        # check that methods in db exist in methods dir
        self.c.execute('SELECT method_id, name FROM method')
        # check method foreign keys
        self.c.execute('SELECT DISTINCT method_id FROM method_parameter')
        self.c.execute('SELECT DISTINCT parent_method_id FROM method_edge')
        self.c.execute('SELECT DISTINCT child_method_id FROM method_edge')
        self.c.execute('SELECT DISTINCT method_id FROM method_path_parent')
        # check edge foreign keys
        self.c.execute('SELECT method_edge_id FROM method_edge')
        self.c.execute('SELECT DISTINCT method_edge_id FROM method_path_edge')
        # check path foreign keys
        self.c.execute('SELECT method_path_id FROM method_path')
        self.c.execute('SELECT DISTINCT method_path_id FROM method_path_edge')
        self.c.execute('SELECT DISTINCT method_path_id FROM method_path_parent')
        self.c.execute(('SELECT DISTINCT parent_method_path_id '
                        'FROM method_path_parent'))
        # check edge closure
        # check path completeness, connectedness, and length concordance
    
    # TODO: move is_inchikey to some sort of helper class
    def is_inchikey(self, inchikey, enforce_standard=False):
        if ('=' in inchikey):
            inchikey = inchikey.split('=')[1]
        if (len(inchikey) == 27):
            s = inchikey.split('-')
            try:
                if (len(s[0]) == 14 and len(s[1]) == 10 and len(s[2]) == 1):
                    if (s[0].isalpha() and s[1].isalpha() and s[2].isalpha()):
                        if (not enforce_standard or s[1][-2] == 'S'):
                            return True
            except IndexError:
                pass
        return False
        

def load():
    # loads the current plugin
    return Check()