# import molecule into mess.db
# Victor Amin 2013

from __future__ import print_function
from __future__ import unicode_literals

import codecs
import math
import os
import re
import sqlite3
import sys
import time
import urllib2
from datetime import datetime

import pybel

from _decorators import decorate, UnicodeDecorator
from _methods import AbstractMethod
from _paths import Path
from _sources import Source
from _utils import get_inchikey_dir, setup_dir, touch

decorate(pybel, UnicodeDecorator)

class Import_(AbstractMethod):
    # method info
    description = 'initial import'
    geop = 0
    # program info
    prog_name = 'Open Babel'
    prog_version = pybel.ob.OBReleaseVersion()
    prog_url = 'http://openbabel.org/wiki/Main_Page'
    # parameters
    parameters = {}
    tags = []

    def check_dependencies(self):
        return True

    def execute(self, args):
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
            # import
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
            # insert molecule to db
            self.update_molecule(inchikey, mol, skip_cir)
            self.import_properties(inchikey, p.path_id, mol)
            # update source catalog in db
            s.update_molecule_source(inchikey, identifier)
            # update source list tsv
            s.update_source_tsv(inchikey_dir, identifier)
            self.check(inchikey, inchikey_dir)
        else:
            # update source catalog in db
            s.update_molecule_source(inchikey, identifier)
            # update source list tsv
            s.update_source_tsv(inchikey_dir, identifier)
            self.status = 'updated'
        self.log(args, inchikey_dir)
        self.db.commit()
        return self.status

    def check(self, inchikey, inchikey_dir):
        inchi = os.path.join(inchikey_dir, inchikey + '.inchi')
        log = os.path.join(inchikey_dir, inchikey + '.log')
        notes = os.path.join(inchikey_dir, inchikey + '.notes')
        png = os.path.join(inchikey_dir, inchikey + '.png')
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
            with codecs.open(sources, encoding='utf-8'):
                pass
            self.status = 'import successful'
            return True
        except IOError:
            self.status = 'import failed'
            return False
    
    def log(self, args, inchikey_dir):
        base_log_path = os.path.join(inchikey_dir, args['inchikey'] + '.log')
        ob_logs_raw = []
        ob_logs = []
        for i in range(3):
            # (0-error, 1-warning, 2-info, 3-audit)
            ob_logs_raw.append(pybel.ob.obErrorLog.GetMessagesOfLevel(i))
        for ll in ob_logs_raw:
            for l in ll:
                ob_logs.append(l)
        self.add_messages_to_log(base_log_path, self.method_name, 
                                 ['status: ' + self.status])
        pybel.ob.obErrorLog.ClearLog()
    
    def import_properties(self, inchikey, method_path_id, mol):
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
        # calculate identifiers with ob/cir, unless entry exists
        q = 'SELECT inchi FROM molecule WHERE inchikey=?'
        inchikey_check_row = self.c.execute(q, (inchikey,)).fetchone()
        inchi = mol.write('inchi').rstrip().split('=')[1]
        if (inchikey_check_row is not None and 
            inchikey_check_row.inchi == inchi):
            if not skip_cir:
                self.update_synonyms(inchikey)
            self.status = 'updated'
            return 0 # this molecule is already correct in the db
        smiles = mol.write('can').rstrip() # canonical smiles
        formula = mol.formula
        iupacs = []
        iupac = ''
        # get identifiers from CIR
        if not skip_cir:
            self.update_synonyms(inchikey)
            try:
                iupacs = self.cir_request(inchikey, 
                                          'iupac_name').splitlines(True)
                # if multiple iupacs, take the longest (most specific) one
                iupac = max(iupacs, key=len).rstrip()
            except AttributeError:
                pass
        # insert molecule identifiers
        self.c.execute("INSERT OR IGNORE INTO molecule \
                (inchikey, inchi, smiles, formula, iupac) \
                VALUES \
                (?, ?, ?, ?, ?)", \
                (inchikey, inchi, smiles, formula, iupac))
        if (len(iupacs) > 1): # if multiple, add others as synonym
            q = ('INSERT OR IGNORE INTO molecule_synonym '
                 '(inchikey, name) '
                 'VALUES (?, ?)')
            for i in iupacs:
                if i != max(iupacs, key=len):
                    self.c.execute(q, (inchikey, i.rstrip()))
    
    def update_synonyms(self, inchikey):
        synonyms = self.cir_request(inchikey, 'names')
        if (synonyms):
            q = ('INSERT OR IGNORE INTO molecule_synonym (inchikey, name) '
                 'VALUES (?, ?)')
            for synonym in (synonyms.split("\n")):
                self.c.execute(q, (inchikey, synonym))
    
    def cir_request(self, inchikey, representation):
        try:
            if not self.cir:
                return None
        except AttributeError:
            self.cir = True
        url = ('http://cactus.nci.nih.gov/chemical/structure/' + inchikey + 
              '/' + representation)
        time.sleep(0.1) # protect cactus from hammering
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
    # loads the current plugin
    return Import_(db)