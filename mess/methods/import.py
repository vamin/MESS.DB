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
from utils import get_inchikey_dir, setup_dir, touch

decorate(pybel, UnicodeDecorator)

class Import_(AbstractMethod):
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
    
    def execute(self, args):
        """Import molecule into MESS.DB."""
        inchikey = args['inchikey']
        mol = args['mol']
        s = args['source']
        p = args['path']
        skip_cir = args['skip_cir']
        inchikey_dir = get_inchikey_dir(inchikey)
        inchikey_basename = os.path.join(inchikey_dir, inchikey)
        setup_dir(inchikey_dir)
        try:
            identifier = args['parent']
        except KeyError:
            identifier = unicode(mol.title, 'utf-8', 'replace')
        # import basic properties
        if not self.check(inchikey, inchikey_dir):
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
            self.update_molecule(inchikey, mol, skip_cir)
            self.import_properties(inchikey, p.path_id, mol)
            self.check(inchikey, inchikey_dir)
        else:
            self.update_molecule(inchikey, mol, skip_cir)
            #self.import_properties(inchikey, p.path_id, mol)
        s.update_molecule_source(inchikey, identifier)
        s.update_source_tsv(inchikey_dir, identifier)
        #if not skip_cir:
        #    self.update_synonyms(inchikey)
        #    self.update_iupac(inchikey)
        self.log(args, inchikey_dir)
        return self.status
    
    def check(self, inchikey, inchikey_dir):
        """Check that a valid molecule folder was created.
        
        Args:
            inchikey: The valid InChIKey for the molecule.
            inchikey_dir: The full path to the molecule's dir.
        
        Sets:
            self.status: Message describing whether import succeeded or failed.
        
        Returns:
            True if everything is fine, False otherwise.
        
        """
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
                row = self.c.execute(q, (inchi_str,)).fetchone()
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
    
    def log(self, args, inchikey_dir):
        """Log messages to base log and method log.
        
        Args:
            args: The args parameter passed to the method (dict expected to
                  contain 'inchikey')
            inchikey_dir: Directory of molecule.
        
        """
        base_log_path = os.path.join(inchikey_dir, '%s.log' % args['inchikey'])
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
        self.add_messages_to_log(base_log_path, self.method_name,
                                 ob_logs)
        self.add_messages_to_log(base_log_path, self.method_name,
                                 ['status: %s' % self.status])
        pybel.ob.obErrorLog.ClearLog()
    
    def import_properties(self, inchikey, method_path_id, mol):
        """Load properties available in Open Babel into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            method_path_id: Path id for import.
            mol: A pybel mol object for the molecule.
        
        """
        # insert Open Babel descriptors
        for property_name, property_value in mol.calcdesc().iteritems():
            if (math.isnan(property_value)):
                continue
            self.insert_property(
                inchikey, method_path_id,
                property_name, 'Open Babel descriptor value',
                type(property_value).__name__, property_value, '')
        # insert Open Babel molecule attributes
        self.insert_property(
            inchikey, method_path_id,
            'charge', 'Open Babel molecule attribute',
            type(mol.charge).__name__, mol.charge, '')
        self.insert_property(
            inchikey, method_path_id,
            'exactmass', 'Open Babel molecule attribute',
            type(mol.exactmass).__name__, mol.exactmass, 'g/mol')
        self.insert_property(
            inchikey, method_path_id,
            'molwt', 'Open Babel descriptor value', type(mol.molwt).__name__,
            mol.molwt, 'g/mol')
        self.insert_property(
            inchikey, method_path_id,
            'spin', 'Open Babel descriptor value', type(mol.spin).__name__,
            mol.spin, '')
        self.db.commit()
    
    def update_molecule(self, inchikey, mol, skip_cir=False):
        """Load basic molecule attributes into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            mol: A pybel mol object for the molecule.
            skip_cir: Set to true to skip querying the CIR web service.
        
        """
        # calculate identifiers with ob/cir, unless entry exists
        q = 'SELECT inchi FROM molecule WHERE inchikey=?'
        inchikey_check_row = self.c.execute(q, (inchikey,)).fetchone()
        inchi = mol.write('inchi').rstrip().split('=')[1]
        if (inchikey_check_row is not None and
            inchikey_check_row.inchi == inchi):
            #if not skip_cir:
            #    self.update_synonyms(inchikey)
            self.status = 'updated'
            return 0 # this molecule is already correct in the db
        smiles = mol.write('can').rstrip() # canonical smiles
        formula = mol.formula
        # insert molecule identifiers
        q = ('INSERT OR IGNORE INTO molecule '
             '(inchikey, inchi, smiles, formula) '
             'VALUES (?, ?, ?, ?)')
        self.c.execute(q, (inchikey, inchi, smiles, formula))
        self.db.commit()
    
    def update_synonyms(self, inchikey):
        """Get synonyms from CIR and load them into mess.db."""
        synonyms = self.cir_request(inchikey, 'names')
        if (synonyms):
            q = ('INSERT OR IGNORE INTO molecule_synonym (inchikey, name) '
                 'VALUES (?, ?)')
            for synonym in (synonyms.split('\n')):
                self.c.execute(q, (inchikey, synonym))
        self.db.commit()
    
    def update_iupac(self, inchikey):
        iupacs = []
        iupac = ''
        try:
            iupacs = self.cir_request(inchikey,
                                      'iupac_name').splitlines(True)
            # if multiple iupacs, take the longest (most specific) one
            iupac = max(iupacs, key=len).rstrip()
        except AttributeError:
            pass
        if (len(iupacs) > 1): # if multiple, add others as synonym
            q = ('INSERT OR IGNORE INTO molecule_synonym '
                 '(inchikey, name) '
                 'VALUES (?, ?)')
            for i in iupacs:
                if i != max(iupacs, key=len):
                    self.c.execute(q, (inchikey, i.rstrip()))
        self.db.commit()
    
    def cir_request(self, inchikey, representation):
        """Make request to CIR (Chemical Information Resolver).
        
        Args:
            inchikey: A valid InChIKey.
            representation: The representation desired from CIR.
        
        Sets:
            self.cir: Sets to True/False based on whether CIR is available to
                      avoid waiting on timeouts.
        
        Returns:
            CIR's response, or None if there isn't one.
        
        """
        try:
            if not self.cir:
                return None
        except AttributeError:
            self.cir = True
        url = 'http://cactus.nci.nih.gov/chemical/structure/%s/%s' %\
              (inchikey, representation)
        time.sleep(0.2) # protect cactus from hammering
        try:
            r = urllib2.urlopen(url)
            if (r.getcode() == 200):
                return r.read()
        except urllib2.URLError as e:
            if hasattr(e, 'reason'):
                print(e.reason, file=sys.stderr)
                print(('CIR is down. Proceeding without importing IUPAC names '
                       'or other synonyms.'), file=sys.stderr)
                self.cir = False
        return None


def load(db):
    """Load Import_(db)."""
    return Import_(db)