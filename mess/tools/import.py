# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB import module

This module contains the import tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import codecs
import math
import os

import pybel

from mess._db import MessDB
from mess._method import AbstractMethod
from mess._path import MethodPath
from mess._source import Source
from mess._tool import AbstractTool
from mess.decorators import decorate, UnicodeDecorator
from mess.utils import get_inchikey_dir, is_inchikey, setup_dir, touch


class ImportTool(AbstractTool):
    """This tool imports molecules into MESS.DB from a source directory."""
    
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
        imp = ImportMethod(MessDB())
        imp.setup()
        source = Source(MessDB())
        source.setup(args.source)
        path = MethodPath()
        path.setup_path(imp.method_id)
        path_id = path.get_path_id()
        pybel.ob.obErrorLog.StopLogging()
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
                            self.log_console.warning(
                                "'%s' is not an importable molecule.",
                                fragment)
                            continue
                        method_args['identifier'] = ('fragment of: %s' %
                                                     unicode(mol.title,
                                                             'utf-8',
                                                             'replace'))
                        method_args['mol'] = frag
                    else:
                        inchikey = mol.write('inchikey').rstrip()
                        if not is_inchikey(inchikey):
                            self.log_console.warning(
                                "'%s' is not an importable molecule.",
                                fragment)
                            continue
                        method_args['mol'] = mol
                    method_args['source'] = source
                    method_args['path_id'] = path_id
                    imp.execute(inchikey, method_args)


class ImportMethod(AbstractMethod):
    """This class adds an individual molecule to MESS.DB."""
    
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
        # setup local variables
        self.inchikey = inchikey
        mol = args['mol']
        source = args['source']
        path_id = args['path_id']
        inchikey_dir = get_inchikey_dir(self.inchikey)
        inchikey_basename = os.path.join(inchikey_dir, self.inchikey)
        try:
            identifier = args['identifier']
        except KeyError:
            identifier = unicode(mol.title, 'utf-8', 'replace')
        # setup directory
        setup_dir(inchikey_dir)
        if not self.check():
            mol.title = b''
            mol.write('inchi',
                      (inchikey_basename + '.inchi'),
                      overwrite=True)
            if not os.path.exists(inchikey_basename + '.png'):
                mol.write('_png2',
                          (inchikey_basename + '.png'))
            touch(inchikey_basename + '.log')
            touch(inchikey_basename + '.notes')
            touch(os.path.join(inchikey_dir,
                               '%s.sources.tsv' % inchikey_basename))
            self.log_all.info('%s molecule directory initialized',
                              self.inchikey)
        self.insert_molecule(self.inchikey, mol)
        self.insert_basic_properties(self.inchikey, path_id, mol)
        source.update_molecule_source(self.inchikey, identifier)
        source.update_source_tsv(self.inchikey, identifier)
        if self.check():
            self.log_all.info('import of %s successful', self.inchikey)
        else:
            self.log_console.warning('import failed for %s', self.inchikey)
        
    def check(self):
        """Check that a valid molecule folder was created and that there is
        a matching molecule in the database.
        
        Args:
            inchikey: The valid InChIKey for the molecule.
            inchikey_dir: The full path to the molecule's dir.
        
        Returns:
            True if everything is fine, False otherwise.
        """
        inchikey_dir = get_inchikey_dir(self.inchikey)
        inchi = os.path.join(inchikey_dir, '%s.inchi' % self.inchikey)
        log = os.path.join(inchikey_dir, '%s.log' % self.inchikey)
        notes = os.path.join(inchikey_dir, '%s.notes' % self.inchikey)
        png = os.path.join(inchikey_dir, '%s.png' % self.inchikey)
        sources = os.path.join(inchikey_dir, '%s.sources.tsv' % self.inchikey)
        try:
            with codecs.open(inchi, encoding='utf-8') as file_:
                inchi_str = file_.readline().split('=')[1].strip()
                query = 'SELECT inchikey FROM molecule WHERE inchi=?'
                row = self.db.execute(query, (inchi_str,)).fetchone()
                try:
                    if row.inchikey != self.inchikey:
                        return False
                except AttributeError:
                    return False
            with codecs.open(log, encoding='utf-8'):
                pass
            with codecs.open(notes, encoding='utf-8'):
                pass
            with codecs.open(png, encoding='utf-8'):
                pass
            with codecs.open(sources, encoding='utf-8'):
                pass
            return True
        except IOError:
            return False
    
    def insert_molecule(self, inchikey, mol):
        """Load basic molecule attributes into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            mol: A pybel mol object for the molecule.
        """
        query = 'SELECT inchi FROM molecule WHERE inchikey=?'
        inchikey_check_row = self.db.execute(query, (inchikey,)).fetchone()
        inchi = mol.write('inchi').rstrip().split('=')[1]
        if (inchikey_check_row is not None and
                inchikey_check_row.inchi == inchi):
            self.status = 'updated'
        else:
            smiles = mol.write('can').rstrip()  # canonical smiles
            formula = mol.formula
            # insert molecule identifiers
            query = ('INSERT OR IGNORE INTO molecule '
                     '(inchikey, inchi, smiles, formula) '
                     'VALUES (?, ?, ?, ?)')
            self.reduce(query, [(inchikey, inchi, smiles, formula)])
    
    def insert_basic_properties(self, inchikey, method_path_id, mol):
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
            if math.isnan(property_value):
                continue
            query, values = self.get_insert_property_query(
                inchikey, method_path_id,
                property_name, 'Open Babel descriptor value',
                type(property_value).__name__, property_value, '')
            all_values.append(values)
        self.reduce(query, all_values)


def load():
    """Load Import()."""
    return ImportTool()
