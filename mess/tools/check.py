# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB check module

This module contains the check tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import glob
import os
import sys

from mess._db import MessDB
from mess._tool import AbstractTool
from mess.utils import is_inchikey


class Check(AbstractTool):
    """This tool verifies that data in the database and the molecules directory
    is self-consistent.
    """
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Check integrity of db/mess.db and molecules dir'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        pass  # no arguments
    
    def execute(self, args):
        """Run self checks."""
        self.db = MessDB()
        cur = self.db.execute('SELECT inchikey FROM molecule')
        self.db_inchikeys = set()
        # check that inchikeys are all valid
        for result in cur:
            if is_inchikey(result.inchikey, enforce_standard=True):
                self.db_inchikeys.add(result.inchikey)
            else:
                self.log_consoleonly.warning(('%s is not a valid standard '
                                              'InChiKey!'), result.inchikey)
        self.check_dir_structure()
        self.check_db_structure()
        self.check_db_dir_inchikey_concordance()
        self.summary()
    
    def check_db_dir_inchikey_concordance(self):
        """Check concordance between inchikey dirs and mess.db."""
        # get list of inchikeys from molecules/ dir
        dir_inchikeys = set()
        inchis = glob.glob(os.path.join(os.path.dirname(__file__),
                                        '../../molecules/*/*/*/', '*.inchi'))
        for i in inchis:
            s = i.split('/')[-1]
            dir_inchikeys.add(s.split('.')[0])
        # compare inchikeys from db vs dir
        in_db_not_dir = self.db_inchikeys - dir_inchikeys
        in_db_msg = '%d InChIKeys in mess.db that are not in molecules dir'
        if len(in_db_not_dir) > 0:
            self.log_consoleonly.warning(in_db_msg, len(in_db_not_dir))
            print('\n'.join(i for i in in_db_not_dir), file=sys.stderr)
        in_dir_not_db = dir_inchikeys - self.db_inchikeys
        in_dir_msg = '%d InChIKeys in molecules dir that are not in mess.db'
        if len(in_dir_not_db) > 0:
            self.log_consoleonly.warning(in_dir_msg, len(in_db_not_dir))
            print('\n'.join(i for i in in_dir_not_db), file=sys.stderr)
    
    def check_dir_structure(self):
        """Check that the structure of the molecules dir is consistent."""
        moldir = os.path.join(os.path.dirname(__file__), '../../molecules')
        for l in os.listdir(moldir):
            lp = os.path.join(moldir, l)
            if not os.path.isdir(lp):
                if 'README' not in l and not l.startswith('.'):
                    self.log_consoleonly.warning(('Unexpected file in '
                                                  'molecules dir: %s'), l)
                continue
            if not len(l) == 1:
                self.log_consoleonly.warning(('Unexpected dir in '
                                              'molecules dir: %s'), l)
                continue
            for ll in os.listdir(lp):
                llp = os.path.join(moldir, l, ll)
                if not os.path.isdir(llp):
                    self.log_consoleonly.warning(('Unexpected file '
                                                  'in molecules dir: '
                                                  '%s/%s'), l, ll)
                    continue
                if not (len(ll) == 2 and ll.isalpha()):
                    self.log_consoleonly.warning(('Unexpected dir '
                                                  'in molecules dir: '
                                                  '%s/%s'), l, ll)
                    continue
                for lll in os.listdir(llp):
                    lllp = os.path.join(moldir, l, ll, lll)
                    if not os.path.isdir(lllp):
                        self.log_consoleonly.warning(('Unexpected file '
                                                      'in molecules dir: '
                                                      '%s/%s/%s'), l, ll, lll)
                        continue
                    if not is_inchikey(l + ll + lll, enforce_standard=True):
                        self.log_consoleonly.warning(('Unexpected dir '
                                                      'in molecules dir: '
                                                      '%s/%s/%s'), l, ll, lll)
                        continue
                    self.check_molecule_dir(l + ll + lll, lllp)
    
    def check_molecule_dir(self, inchikey, d):
        """Check that molecule directory has the proper file in it."""
        l = os.listdir(d)
        if not inchikey + '.inchi' in l:
            self.log_consoleonly.warning('%s does not contain inchi',
                                         inchikey)
        if not inchikey + '.log' in l:
            self.log_consoleonly.warning('%s does not contain log file',
                                         inchikey)
        if not inchikey + '.notes' in l:
            self.log_consoleonly.warning('%s does not contain notes',
                                         inchikey)
        if not inchikey + '.png' in l:
            self.log_consoleonly.warning('%s does not contain png',
                                         inchikey)
        if 'sources.tsv' not in l:
            self.log_consoleonly.warning(('%s does not contain '
                                          'sources.tsv file'), inchikey)
        else:
            self.check_sources_tsv(os.path.join(d, 'sources.tsv'))
        for ll in l:
            if os.path.isdir(ll):
                self.check_method_dir(os.path.join(d, ll))
    
    def check_sources_tsv(self, s):
        """Check that sources.tsv files are consistent with mess.db."""
        pass
    
    def check_method_dir(self, d):
        """Check that method dirs are consistent with mess.db."""
        pass
    
    def check_db_structure(self):
        """Check that DB structure/data is self-consistent."""
        # check inchikey foreign keys
        molecule_synonym_inchikeys = set()
        q = 'SELECT DISTINCT inchikey FROM molecule_synonym'
        for r in self.db.execute(q).fetchall():
            molecule_synonym_inchikeys.add(r.inchikey)
        molecule_source_inchikeys = set()
        q = 'SELECT DISTINCT inchikey FROM molecule_source'
        for r in self.db.execute(q).fetchall():
            molecule_source_inchikeys.add(r.inchikey)
        q = 'SELECT DISTINCT inchikey FROM molecule_method_property'
        molecule_method_property_inchikeys = set()
        for r in self.db.execute(q).fetchall():
            molecule_method_property_inchikeys.add(r.inchikey)
        loose_keys_mol_syn = molecule_synonym_inchikeys - self.db_inchikeys
        loose_keys_mol_src = molecule_source_inchikeys - self.db_inchikeys
        loose_keys_mol_met_prp = (molecule_method_property_inchikeys -
                                  self.db_inchikeys)
        if len(loose_keys_mol_syn) > 0:
            self.log_consoleonly.warning(('%d loose InChIKeys in '
                                          'molecule_synonym table'),
                                         len(loose_keys_mol_syn))
            print('\n'.join(i for i in loose_keys_mol_syn), file=sys.stderr)
        if len(loose_keys_mol_src) > 0:
            self.log_consoleonly.warning(('%d loose InChIKeys in '
                                          'molecule_source table'),
                                         len(loose_keys_mol_src))
            print('\n'.join(i for i in loose_keys_mol_src), file=sys.stderr)
        if len(loose_keys_mol_met_prp) > 0:
            self.log_consoleonly.warning(('%d loose InChIKeys in '
                                          'molecule_method_property table'),
                                         len(loose_keys_mol_met_prp))
            print('\n'.join(i for i in loose_keys_mol_met_prp),
                  file=sys.stderr)
        # check that sources in db exist in sources dir
        q = 'SELECT source_id, dirname FROM source'
        source_ids = set()
        source_path = os.path.join(os.path.dirname(__file__), '../../sources')
        for r in self.db.execute(q).fetchall():
            source_ids.add(r.source_id)
            if not os.path.isdir(os.path.join(source_path, r.dirname)):
                self.log_consoleonly.warning('%s not in sources directory',
                                             r.dirname)
        # check source foreign keys
        molecule_source_ids = set()
        q = 'SELECT DISTINCT source_id FROM molecule_source'
        for r in self.db.execute(q).fetchall():
            molecule_source_ids.add(r.source_id)
        loose_source_ids = molecule_source_ids - source_ids
        if len(loose_source_ids) > 0:
            self.log_consoleonly.warning(('%d loose source_ids in '
                                          'molecule_source table'),
                                         len(loose_source_ids))
            print('\n'.join(i for i in loose_source_ids), file=sys.stderr)
        # check program foreign keys
        program_ids = set()
        q = 'SELECT program_id FROM program'
        for r in self.db.execute(q).fetchall():
            program_ids.add(r.program_id)
        method_program_ids = set()
        q = 'SELECT DISTINCT program_id FROM method'
        for r in self.db.execute(q).fetchall():
            method_program_ids.add(r.program_id)
        loose_program_ids = method_program_ids - program_ids
        if len(loose_program_ids) > 0:
            self.log_consoleonly.warning(('%d loose program_ids in '
                                          'method table'),
                                         len(loose_program_ids))
            print('\n'.join(i for i in loose_program_ids), file=sys.stderr)
        # check parameter foreign keys
        parameter_ids = set()
        q = 'SELECT parameter_id FROM parameter'
        for r in self.db.execute(q).fetchall():
            parameter_ids.add(r.parameter_id)
        method_parameter_ids = set()
        q = 'SELECT DISTINCT parameter_id FROM method_parameter'
        for r in self.db.execute(q).fetchall():
            method_parameter_ids.add(r.parameter_id)
        method_tag_parameter_ids = set()
        q = 'SELECT DISTINCT parameter_id FROM method_tag'
        for r in self.db.execute(q).fetchall():
            method_tag_parameter_ids.add(r.parameter_id)
        loose_m_pids = method_parameter_ids - parameter_ids
        if len(loose_m_pids) > 0:
            self.log_consoleonly.warning(('%d loose parameter_ids in '
                                          'method_parameter table'),
                                         len(loose_m_pids))
            print('\n'.join(i for i in loose_m_pids), file=sys.stderr)
        loose_mt_pids = method_tag_parameter_ids - parameter_ids
        if len(loose_mt_pids) > 0:
            self.log_consoleonly.warning(('%d loose parameter_ids in '
                                          'method_tag table'),
                                         len(loose_mt_pids))
            print('\n'.join(i for i in loose_mt_pids), file=sys.stderr)
        property_ids = set()
        q = 'SELECT property_id FROM property'
        for r in self.db.execute(q).fetchall():
            property_ids.add(r.property_id)
        mmp_property_ids = set()
        q = 'SELECT DISTINCT property_id FROM molecule_method_property'
        for r in self.db.execute(q).fetchall():
            mmp_property_ids.add(r.property_id)
        loose_mmp_property_ids = mmp_property_ids - property_ids
        if len(loose_mmp_property_ids) > 0:
            self.log_consoleonly.warning(('%d loose property_ids in '
                                          'method_tag table'),
                                         len(loose_mmp_property_ids))
            print('\n'.join(i for i in loose_mmp_property_ids),
                  file=sys.stderr)
        # check that methods in db exist in methods dir
        method_ids = set()
        q = 'SELECT method_id, name FROM method'
        for r in self.db.execute(q).fetchall():
            method_ids.add(r.method_id)
        # check method foreign keys
        method_parameter_mids = set()
        q = 'SELECT DISTINCT method_id FROM method_parameter'
        for r in self.db.execute(q).fetchall():
            method_parameter_mids.add(r.method_id)
        method_edge_mids = set()
        q = ('SELECT DISTINCT parent_method_id, child_method_id  '
             'FROM method_edge')
        for r in self.db.execute(q).fetchall():
            method_edge_mids.add(r.parent_method_id)
            method_edge_mids.add(r.child_method_id)
        loose_method_parameter_mids = method_parameter_mids - method_ids
        if len(loose_method_parameter_mids) > 0:
            self.log_consoleonly.warning(('%d loose method_ids in '
                                          'method_parameter table'),
                                         len(loose_method_parameter_mids))
            print('\n'.join(i for i in loose_method_parameter_mids),
                  file=sys.stderr)
        loose_method_edge_mids = method_edge_mids - method_ids
        if len(loose_method_edge_mids) > 0:
            self.log_consoleonly.warning(('%d loose method_ids in '
                                          'method_edge table'),
                                         len(loose_method_edge_mids))
            print('\n'.join(i for i in loose_method_edge_mids),
                  file=sys.stderr)
        # check edge foreign keys
        method_edge_ids = set()
        q = 'SELECT method_edge_id FROM method_edge'
        for r in self.db.execute(q).fetchall():
            method_edge_ids.add(r.method_edge_id)
        method_path_edge_ids = set()
        q = 'SELECT DISTINCT method_edge_id FROM method_path_edge'
        for r in self.db.execute(q).fetchall():
            method_path_edge_ids.add(r.method_edge_id)
        loose_method_path_edge_ids = method_path_edge_ids - method_edge_ids
        if len(loose_method_path_edge_ids) > 0:
            self.log_consoleonly.warning(('%d loose method_edge_ids in '
                                          'method_path_edge table'),
                                         len(loose_method_path_edge_ids))
            print('\n'.join(i for i in loose_method_path_edge_ids),
                  file=sys.stderr)
        # check path foreign keys
        method_path_ids = set()
        q = 'SELECT method_path_id FROM method_path'
        for r in self.db.execute(q).fetchall():
            method_path_ids.add(r.method_path_id)
        method_path_edge_pids = set()
        q = 'SELECT DISTINCT method_path_id FROM method_path_edge'
        for r in self.db.execute(q).fetchall():
            method_path_edge_pids.add(r.method_path_id)
        loose_method_path_edge_pids = method_path_edge_pids - method_path_ids
        if len(loose_method_path_edge_pids) > 0:
            self.log_consoleonly.warning(('%d loose method_path_ids in '
                                          'method_path_edge table'),
                                         len(loose_method_path_edge_pids))
            print('\n'.join(i for i in loose_method_path_edge_pids),
                  file=sys.stderr)
        # check edge closure
        # check path completeness, connectedness, and length concordance
    
    def summary(self):
        """Print summary statistics about molecules in MESS.DB."""
        self.log_consoleonly.info('%d molecules in MESS.DB',
                                  len(self.db_inchikeys))


def load():
    """Load Check()."""
    return Check()
