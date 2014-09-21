# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB annotate module

This module contains the annotate tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import json
import os
import sys
import time
import urllib2

import pybel

from mess._db import MessDB
from mess._path import MethodPath
from mess._tool import AbstractTool
from mess.utils import get_inchikey_dir
from mess.tools.match import Match


class Annotate(AbstractTool):
    """This tool annotates molecules with synonyms (common names, CAS, etc) and
    various fingerprints.
    """
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Annotate molecules with synonyms and fingerprints'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('a list of inchikeys (default: STDIN)'))
        subparser.add_argument('-c', '--cir', action='store_true',
                               help=('get IUPAC names and other synonyms from '
                                     'the Chemical Information Resolver web '
                                     'service'))
        subparser.add_argument('-f', '--fingerprint', type=str,
                               choices=[b'FP2', b'FP3', b'FP4',
                                        b'MACCS', b'MNA', b'MPD'],
                               help=('calculate fingerprint'))
        subparser.add_argument('-s', '--spectrophore', action='store_true',
                               help=('calculate Spectrophore '
                                     'descriptor/fingerprint; '
                                     'requires 3D geometry (i.e., you must '
                                     'set a path to a method that has '
                                     'generated xyz coordinates)'))
        subparser.add_argument('-p', '--path', type=int, default=None,
                               help=('specify a path id, only used for '
                                     'Spectrophore'))
        sp_group = subparser.add_argument_group(('Spectrophore optional '
                                                 'arguments'))
        sp_group.add_argument('-sn', '--spectrophore-normalization', type=str,
                              default='No', choices=[b'No',
                                                     b'ZeroMean',
                                                     b'UnitStd',
                                                     b'ZeroMeanAndUnitStd'],
                              help=('perform normalization of Spectrophore'))
        sp_group.add_argument('-sa', '--spectrophore-accuracy', type=int,
                              default=20, choices=[1, 2, 5, 10, 15,
                                                   20, 30, 36, 45, 60],
                              help=('Spectrophore accuracy expressed as '
                                    'angular stepsize; lower is more accurate'
                                    'but slower'))
        sp_group.add_argument('-ss', '--spectrophore-stereospecificity',
                              type=str, default='No',
                              choices=[b'No', b'Unique', b'Mirror', b'All'],
                              help=('cage type in terms of the underlying '
                                    'pointgroup: P1 or P-1'))
        sp_group.add_argument('-sr', '--spectrophore-resolution', type=float,
                              metavar='FLOAT', default=3.0,
                              help=('required Spectrophore resolution in '
                                    'Angstroms'))
    
    def execute(self, args):
        """Match molecules to SMARTS patterns."""
        if args.inchikeys.name == '<stdin>' and args.inchikeys.isatty():
            sys.exit('No input specified.')
        if not (args.cir or args.fingerprint or args.spectrophore):
            sys.exit('You did not request any annotations.')
        if args.spectrophore:
            if args.path is None:
                sys.exit(('Spectrophore calculation requires 3D geometry. '
                          'You must specify a 3D geometry with --path.'))
            else:
                path = MethodPath()
                path.set_path(args.path)
                method_dir = path.get_path_directory()
                sp_args = {'normalization': args.spectrophore_normalization,
                           'accuracy': args.spectrophore_accuracy,
                           'stereo': args.spectrophore_stereospecificity,
                           'resolution': args.spectrophore_resolution}
        self.db = MessDB()
        inchi_select_query = 'SELECT inchi FROM molecule WHERE inchikey = ?'
        fp_select_query = ('SELECT fingerprint FROM molecule_fingerprint '
                           'WHERE inchikey = ? '
                           'AND name = ? '
                           'AND settings = ? '
                           'AND method_path_id = ?')
        fp_insert_query = ('INSERT INTO molecule_fingerprint '
                           '(inchikey, name, settings, '
                           'fingerprint, method_path_id) '
                           'VALUES (?, ?, ?, ?, ?)')
        for row in args.inchikeys:
            self.inchikey = row.split()[0].strip()
            if args.cir:
                self.update_iupac(self.inchikey)
                self.update_synonyms(self.inchikey)
            if args.fingerprint:
                inchi = self.db.execute(inchi_select_query,
                                        (self.inchikey,)).fetchone()[0]
                mol = pybel.readstring('inchi', 'InChI=%s' % inchi)
                canonical = pybel.ob.OBOp.FindType(b'canonical')
                canonical.Do(mol.OBMol)
                fp = Match.calculate_fingerprint(mol, args.fingerprint)
                try:
                    db_fp = self.db.execute(fp_select_query,
                                            (self.inchikey,
                                             args.fingerprint,
                                             '',
                                             '')).fetchone()[0]
                    if not str(fp) == db_fp:
                        self.log_console.warning(('new %s fingerprint '
                                                  'for %s did not match '
                                                  'fingerprint in db, '
                                                  'db not updated'),
                                                 args.fingerprint,
                                                 self.inchikey)
                except TypeError:
                    self.db.execute(fp_insert_query, (self.inchikey,
                                                      args.fingerprint,
                                                      '',
                                                      str(fp),
                                                      ''))
                    self.log_all.info('%s fingerprint for %s added to db',
                                      args.fingerprint, self.inchikey)
            if args.spectrophore:
                xyz_file = os.path.join(get_inchikey_dir(self.inchikey),
                                        method_dir,
                                        '%s.xyz' % self.inchikey)
                mol = pybel.readfile('xyz', xyz_file).next()
                sp = Match.calculate_spectrophore(mol, sp_args)
                try:
                    db_sp = self.db.execute(fp_select_query,
                                            (self.inchikey,
                                             'Spectrophore',
                                             json.dumps(sp_args,
                                                        sort_keys=True),
                                             args.path)).fetchone()[0]
                    if not str(sp) == db_sp:
                        self.log_console.warning(('new Spectrophore '
                                                  'fingerprint for '
                                                  '%s did not match '
                                                  'fingerprint in db, '
                                                  'db not updated'),
                                                 self.inchikey)
                except TypeError:
                    json_sp_args = json.dumps(sp_args, sort_keys=True)
                    self.db.execute(fp_insert_query, (self.inchikey,
                                                      'Spectrophore',
                                                      json_sp_args,
                                                      str(sp),
                                                      args.path))
                    self.log_all.info(('Spectrophore fingerprint for %s '
                                       'with parameters %s and '
                                       'geometry from path %i '
                                       'added to db'),
                                      self.inchikey, json_sp_args, args.path)
    
    def update_synonyms(self, inchikey):
        """Get synonyms from CIR and load them into mess.db."""
        new_synonyms = 0
        synonyms = self.cir_request(inchikey, 'names')
        if synonyms:
            select_query = ('SELECT inchikey FROM molecule_synonym '
                            'WHERE inchikey = ? AND name = ?')
            insert_query = ('INSERT INTO molecule_synonym (inchikey, name) '
                            'VALUES (?, ?)')
            for synonym in synonyms.split('\n'):
                if self.db.execute(select_query,
                                   (inchikey, synonym)).fetchone() is None:
                    self.db.execute(insert_query, (inchikey, synonym))
                    new_synonyms += 1
            if new_synonyms > 0:
                if new_synonyms > 1:
                    plural = 's'
                else:
                    plural = ''
                self.log_all.info('%i new synonym%s for %s added',
                                  new_synonyms, plural, inchikey)
    
    def update_iupac(self, inchikey):
        """Get IUPAC name from CIR and load it into mess.db."""
        iupacs = []
        iupac = None
        new_synonyms = 0
        try:
            iupacs = self.cir_request(inchikey,
                                      'iupac_name').splitlines(True)
            # if multiple iupacs, take the longest (most specific) one
            iupac = max(iupacs, key=len).rstrip()
        except AttributeError:
            return
        if iupac is not None:
            iupac_select_query = ('SELECT iupac FROM molecule '
                                  'WHERE inchikey = ?')
            iupac_update_query = ('UPDATE molecule SET iupac = ? '
                                  'WHERE inchikey = ?')
            db_iupac = self.db.execute(iupac_select_query,
                                       (inchikey, )).fetchone()[0]
            if not db_iupac == iupac:
                self.db.execute(iupac_update_query, (iupac, inchikey))
                self.log_all.info('iupac name for %s updated', inchikey)
            if len(iupacs) > 1:  # if multiple, add others as synonym
                select_query = ('SELECT inchikey FROM molecule_synonym '
                                'WHERE inchikey = ? AND name = ?')
                insert_query = ('INSERT INTO molecule_synonym '
                                '(inchikey, name) VALUES (?, ?)')
                for i in iupacs:
                    if i != max(iupacs, key=len):  # ignore longest iupac
                        synonym = i.rstrip()
                        if self.db.execute(select_query,
                                           (inchikey,
                                            synonym)).fetchone() is None:
                            self.db.execute(insert_query, (inchikey,
                                                           synonym))
                            new_synonyms += 1
                if new_synonyms > 0:
                    if new_synonyms > 1:
                        plural = 's'
                    else:
                        plural = ''
                    self.log_all.info('%i new synonym%s for %s added',
                                      new_synonyms, plural, inchikey)
    
    def cir_request(self, inchikey, representation):
        """Make request to CIR (Chemical Information Resolver).
        
        Args:
            inchikey: A valid InChIKey.
            representation: The representation desired from CIR.
        
        Returns:
            CIR's response, or None if there isn't one.
        """
        url = 'http://cactus.nci.nih.gov/chemical/structure/%s/%s' %\
              (inchikey, representation)
        headers = {'User-Agent': 'MESS.DB'}
        request = urllib2.Request(url, None, headers)
        try:
            response = urllib2.urlopen(request)
            if response.getcode() == 200:
                time.sleep(0.2)  # protect cactus from hammering
                return response.read()
        except urllib2.URLError as err:
            if hasattr(err, 'reason'):
                reason = err.reason.lower()
                self.log_console.info('%s %s %s in cir',
                                      inchikey, representation, reason)
        return None


def load():
    """Load Annotate()."""
    return Annotate()
