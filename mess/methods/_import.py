from __future__ import print_function
from __future__ import unicode_literals

import codecs
import math
import os
import re
import sys
import time
import urllib2

import pybel

from _method import AbstractMethod
from _path import Path
from _source import Source
from decorators import decorate, UnicodeDecorator
from utils import get_inchikey_dir, setup_dir, touch, write_to_log

decorate(pybel, UnicodeDecorator)

class Import(AbstractMethod):
    # method info
    description = 'Initial import'
    geop = 0
    # program info
    prog_name = 'Open Babel'
    prog_version = ''
    prog_url = 'http://openbabel.org/wiki/Main_Page'
    # parameters
    parameters = {}
    tags = []
    
    def check_dependencies(self):
        """Return True, no external dependencies to check."""
        return True
    
    def map(self):
        pass
    
    def execute(self, inchikey, args):
        """Import molecule into MESS.DB."""
        mol = args['mol']
        s = args['source']
        p = args['path']
        inchikey_dir = get_inchikey_dir(inchikey)
        inchikey_basename = os.path.join(inchikey_dir, inchikey)
        setup_dir(inchikey_dir)
        try:
            identifier = args['parent']
        except KeyError:
            identifier = unicode(mol.title, 'utf-8', 'replace')
        # import basic properties
        if not self.check(inchikey):
            mol.title = b''
            mol.write('inchi',
                      (inchikey_basename + '.inchi'),
                      overwrite=True)
            if not os.path.exists(inchikey_basename + '.png'):
                mol.write('_png2',
                          (inchikey_basename + '.png'))
            touch(inchikey_basename + '.log')
            touch(inchikey_basename + '.notes')
            touch(os.path.join(inchikey_dir, 'sources.tsv'))
            with codecs.open(os.path.join(inchikey_dir, 'charge.txt'), 
                             'w', 'utf-8') as c:
                c.write('%i' % mol.charge)
            self.update_molecule(inchikey, mol)
            self.import_properties(inchikey, p.path_id, mol)
            self.check(inchikey)
        else:
            self.update_molecule(inchikey, mol)
        s.update_molecule_source(inchikey, identifier)
        s.update_source_tsv(inchikey_dir, identifier)
        self.log(inchikey, self.status)
        print('%s: %s' % (inchikey, self.status))
    
    def check(self, inchikey):
        """Check that a valid molecule folder was created.
        
        Args:
            inchikey: The valid InChIKey for the molecule.
            inchikey_dir: The full path to the molecule's dir.
        
        Sets:
            self.status: Message describing whether import succeeded or failed.
        
        Returns:
            True if everything is fine, False otherwise.
        
        """
        inchikey_dir = get_inchikey_dir(inchikey)
        inchi = os.path.join(inchikey_dir, inchikey + '.inchi')
        log = os.path.join(inchikey_dir, inchikey + '.log')
        notes = os.path.join(inchikey_dir, inchikey + '.notes')
        png = os.path.join(inchikey_dir, inchikey + '.png')
        charge = os.path.join(inchikey_dir, 'charge.txt')
        sources = os.path.join(inchikey_dir, 'sources.tsv')
        try:
            with codecs.open(inchi, encoding='utf-8') as f:
                inchi_str = f.readline().split('=')[1].strip()
                q = 'SELECT inchikey FROM molecule WHERE inchi=?'
                row = self.cur.execute(q, (inchi_str,)).fetchone()
                try:
                    if (row.inchikey != inchikey):
                        self.status = 'import failed'
                        return False
                except AttributeError:
                    self.status = 'import failed'
                    return False
            with codecs.open(log, encoding='utf-8'):
                pass
            with codecs.open(notes, encoding='utf-8'):
                pass
            with codecs.open(png, encoding='utf-8'):
                pass
            with codecs.open(charge, encoding='utf-8'):
                pass
            with codecs.open(sources, encoding='utf-8'):
                pass
            self.status = 'import successful'
            return True
        except IOError:
            self.status = 'import failed'
            return False
    
    def log(self, inchikey, status):
        """Log messages to base log and method log.
        
        Args:
            inchikey:
            status:
        
        """
        base_log_path = os.path.join(get_inchikey_dir(inchikey),
                                     '%s.log' % inchikey)
        ob_logs_raw = []
        ob_logs = []
        for i in range(3):
            # (0-error, 1-warning, 2-info, 3-audit)
            ob_logs_raw.append(pybel.ob.obErrorLog.GetMessagesOfLevel(i))
        for logs in ob_logs_raw:
            for log in logs:
                for line in iter(log.splitlines()):
                    if not line.startswith('==='):
                        ob_logs.append(line)
        write_to_log(base_log_path, ob_logs)
        write_to_log(base_log_path, ['status: %s' % status])
        pybel.ob.obErrorLog.ClearLog()
    
    def import_properties(self, inchikey, method_path_id, mol):
        """Load properties available in Open Babel into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            method_path_id: Path id for import.
            mol: A pybel mol object for the molecule.
        
        """
        # insert Open Babel molecule attributes
        query, values = self.get_insert_property_query(
            inchikey, method_path_id,
            'charge', 'Open Babel molecule attribute',
            type(mol.charge).__name__, mol.charge, '')
        all_values = [values]
        query, values = self.get_insert_property_query(
            inchikey, method_path_id,
            'exactmass', 'Open Babel molecule attribute',
            type(mol.exactmass).__name__, mol.exactmass, 'g/mol')
        all_values.append(values)
        query, values = self.get_insert_property_query(
            inchikey, method_path_id,
            'molwt', 'Open Babel descriptor value', type(mol.molwt).__name__,
            mol.molwt, 'g/mol')
        all_values.append(values)
        query, values = self.get_insert_property_query(
            inchikey, method_path_id,
            'spin', 'Open Babel descriptor value', type(mol.spin).__name__,
            mol.spin, '')
        all_values.append(values)
        # insert Open Babel descriptors
        for property_name, property_value in mol.calcdesc().iteritems():
            if (math.isnan(property_value)):
                continue
            query, values = self.get_insert_property_query(
                inchikey, method_path_id,
                property_name, 'Open Babel descriptor value',
                type(property_value).__name__, property_value, '')
            all_values.append(values)
        self.reduce(query, all_values)
        
        
    
    def update_molecule(self, inchikey, mol):
        """Load basic molecule attributes into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            mol: A pybel mol object for the molecule.
        """
        # calculate identifiers with ob/cir, unless entry exists
        query = 'SELECT inchi FROM molecule WHERE inchikey=?'
        inchikey_check_row = self.cur.execute(query, (inchikey,)).fetchone()
        inchi = mol.write('inchi').rstrip().split('=')[1]
        if (inchikey_check_row is not None and
            inchikey_check_row.inchi == inchi):
            self.status = 'updated'
        else:
            smiles = mol.write('can').rstrip() # canonical smiles
            formula = mol.formula
            # insert molecule identifiers
            query = ('INSERT OR IGNORE INTO molecule '
                     '(inchikey, inchi, smiles, formula) '
                     'VALUES (?, ?, ?, ?)')
            self.reduce(query, [(inchikey, inchi, smiles, formula)])


def load(db):
    """Load Import(db)."""
    return Import(db)