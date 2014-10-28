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

from mess.decorators import decorate, UnicodeDecorator
from mess.method import AbstractMethod
from mess.source import Source
from mess.tool import AbstractTool
from mess.utils import get_inchikey_dir, is_inchikey, setup_dir, touch


class Import(AbstractTool):
    """This tool imports molecules into MESS.DB from a source directory."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Import molecules into MESS.DB'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('source',
                               help='a molecule source file or directory')
        subparser.add_argument('-k', '--skip-fragments', action='store_true',
                               help=('do not attempt to separate and import '
                                     'non-covalently bound fragments'))
    
    def execute(self, args):
        """Run import method for every molecule in source."""
        source = Source()
        source.setup(args.source)
        self.log_console.info('reading molecules')
        molecules = {}
        threedee = False
        pybel.ob.obErrorLog.SetOutputLevel(-1)
        for source_file in source.files():
            for mol in pybel.readfile(source_file.split('.')[-1],
                                      os.path.join(source.source_dir,
                                                   source_file)):
                if not threedee and mol.dim == 3:
                    threedee = True
                try:
                    decorate(mol, UnicodeDecorator)
                except IndexError:
                    self.log_console.error('Unexpected error importing %s.',
                                           mol.title)
                    continue
                inchikey = mol.write('inchikey').rstrip()
                if not is_inchikey(inchikey):
                    self.log_console.warning(
                        ("'%s' is not an importable molecule."), mol.title)
                    continue
                molecules[inchikey] = (mol, source)
                if not args.skip_fragments:
                    cansmi = mol.write('can').split()[0]
                    if cansmi.count('.') > 0:
                        for fragment in cansmi.split('.'):
                            fragmol = pybel.readstring('can', fragment)
                            decorate(fragmol, UnicodeDecorator)
                            inchikey = fragmol.write('inchikey').rstrip()
                            if not is_inchikey(inchikey):
                                self.log_console.warning(
                                    ("'%s' fragment in %s "
                                     "is not an importable molecule."),
                                    fragment, mol.title)
                            else:
                                fragmol.title = mol.title
                                molecules[inchikey] = (fragmol, source)
        import0d = Import0D()
        import0d.setup()
        if threedee:
            import3d = Import3D()
            import3d.shortdesc = source.dirname
            import3d.setup()
        self.log_console.info('setting up molecule dirs')
        queries = {}
        for inchikey, (mol, source) in molecules.iteritems():
            for query, values in import0d.map(mol, source):
                try:
                    queries[query].append(values)
                except KeyError:
                    queries[query] = [values]
            if mol.dim == 3:
                import3d.map(mol, source)
        self.log_console.info('loading simple properties')
        for query, values in queries.iteritems():
            import0d.reduce(query, values)


class Import0D(AbstractMethod):
    """This class adds an individual 0D molecule to MESS.DB."""
    
    # method info
    description = 'import0d'
    geop = 0
    # program info
    prog_name = 'Open Babel'
    prog_version = pybel.ob.OBReleaseVersion()
    prog_url = 'http://openbabel.org/wiki/Main_Page'
    prog_citation = ('Noel M. O’Boyle, Michael Banck, Craig A. James, '
                     'Chris Morley, Tim Vandermeersch, Geoffrey R. Hutchison '
                     'Open Babel: An open chemical toolbox. J. Cheminf. '
                     '2011, 3, 33.')
    # parameters
    parameters = {}
    
    def check_dependencies(self):
        """Return True, no external dependencies to check."""
        return True

    def map(self, mol, source):
        """Import molecule into MESS.DB."""
        # setup local variables
        self.inchikey = mol.write('inchikey').rstrip()
        inchikey_dir = get_inchikey_dir(self.inchikey)
        inchikey_basename = os.path.join(inchikey_dir, self.inchikey)
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
        source.update_source_tsv(self.inchikey, identifier)
        yield source.update_molecule_source_query(self.inchikey, identifier)
        yield self.insert_molecule_query(self.inchikey, mol)
        for query, values in self.insert_property_queries(self.inchikey,
                                                          self.path_id,
                                                          mol):
            yield query, values
        
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
    
    def insert_molecule_query(self, inchikey, mol):
        """Load basic molecule attributes into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            mol: A pybel mol object for the molecule.
        """
        inchi = mol.write('inchi').rstrip().split('=')[1]
        smiles = mol.write('can').rstrip()  # canonical smiles
        formula = mol.formula
        # insert molecule identifiers
        query = ('INSERT OR IGNORE INTO molecule '
                 '(inchikey, inchi, smiles, formula) '
                 'VALUES (?, ?, ?, ?)')
        return (query, (inchikey, inchi, smiles, formula))
    
    def insert_property_queries(self, inchikey, method_path_id, mol):
        """Load properties available in Open Babel into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            method_path_id: Path id for import.
            mol: A pybel mol object for the molecule.
        
        """
        # insert Open Babel molecule attributes
        yield self.get_insert_property_query(
            inchikey, method_path_id,
            'charge', 'Open Babel molecule attribute',
            type(mol.charge).__name__, mol.charge, '')
        yield self.get_insert_property_query(
            inchikey, method_path_id,
            'exactmass', 'Open Babel molecule attribute',
            type(mol.exactmass).__name__, mol.exactmass, 'g/mol')
        yield self.get_insert_property_query(
            inchikey, method_path_id,
            'molwt', 'Open Babel descriptor value', type(mol.molwt).__name__,
            mol.molwt, 'g/mol')
        yield self.get_insert_property_query(
            inchikey, method_path_id,
            'spin', 'Open Babel descriptor value', type(mol.spin).__name__,
            mol.spin, '')
        # insert Open Babel descriptors
        for property_name, property_value in mol.calcdesc().iteritems():
            if math.isnan(property_value):
                continue
            yield self.get_insert_property_query(
                inchikey, method_path_id,
                property_name, 'Open Babel descriptor value',
                type(property_value).__name__, property_value, '')


class Import3D(AbstractMethod):
    """This class adds an individual 3D molecule to MESS.DB."""
    
    # method info
    description = 'import3d'
    geop = 1
    # program info
    prog_name = 'Open Babel'
    prog_version = pybel.ob.OBReleaseVersion()
    prog_url = 'http://openbabel.org/wiki/Main_Page'
    prog_citation = ('Noel M. O’Boyle, Michael Banck, Craig A. James, '
                     'Chris Morley, Tim Vandermeersch, Geoffrey R. Hutchison '
                     'Open Babel: An open chemical toolbox. J. Cheminf. '
                     '2011, 3, 33.')
    # parameters
    parameters = {}
    
    def check_dependencies(self):
        """Return True, no external dependencies to check."""
        return True

    def map(self, mol, source):
        """Import molecule into MESS.DB."""
        self.inchikey = mol.write('inchikey').rstrip()
        inchikey_dir = get_inchikey_dir(self.inchikey)
        setup_dir(os.path.join(inchikey_dir, self.method_dir))
        mol.write('xyz',
                  os.path.join(inchikey_dir,
                               self.method_dir,
                               '%s.xyz' % self.inchikey),
                  overwrite=True)
        self.log_all.info('%s 3D structure from %s added',
                          self.inchikey, source.dirname)
    
    def check(self):
        pass


def load():
    """Load Import()."""
    return Import()
