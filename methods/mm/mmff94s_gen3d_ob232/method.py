#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

# generate 3d structure with obabel

import math
import os
import sqlite3
import sys
from datetime import datetime

import pybel

from _methods import AbstractMethod

class Method(AbstractMethod):
    def execute(self, args):
        inchikey_dir = self.get_inchikey_dir(args['inchikey'])
        mol = pybel.readfile('inchi', os.path.join(inchikey_dir, args['inchikey'] + '.inchi')).next()
        mol.localopt('mmff94s', 1024)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, args['method_dir']))
        self.setup_dir(out_dir)
        xyz_out = os.path.join(out_dir, args['inchikey'] + '.xyz')
        if not self.check(xyz_out):
            mol.write('xyz', str(xyz_out))
            self.import_properties(args['inchikey'], args['path_id'], mol)
            self.check(xyz_out)
        else:
            self.status = 'skipped'
        self.log(args, inchikey_dir)
        return self.status

    def check(self, xyz_out):
        try:
            with open(xyz_out) as f:
                for i, l in enumerate(f):
                    if (i == 0):
                        atoms = l.strip()
            lines = i + 1
            if (int(atoms) == lines - 2):
                self.status = 'generating 3D coordinates succeeded'
                return True
            else:
                self.status = 'generating 3D coordinates failed'
                return False
        except IOError:
            self.status = 'generating 3D coordinates failed'
            return False
    
    def log(self, args, inchikey_dir):
        base_log_path = os.path.join(inchikey_dir, args['inchikey'] + '.log')
        method_log_path = os.path.join(inchikey_dir, args['method_dir'], args['inchikey'] + '.log')
        ob_logs_raw = []
        ob_logs = []
        for i in range(3):
            ob_logs_raw.append(pybel.ob.obErrorLog.GetMessagesOfLevel(i)) # (0-error, 1-warning, 2-info, 3-audit)
        for ll in ob_logs_raw:
            for l in ll:
                ob_logs.append(l)
        self.add_messages_to_log(base_log_path, args['method_name'], ['status: ' + self.status])
        self.add_messages_to_log(method_log_path, args['method_name'], ['status: ' + self.status] + ob_logs)
        pybel.ob.obErrorLog.ClearLog()
    
    def setup_properties(self):
        mol = pybel.readstring('inchi', 'InChI=1/H')
        # insert Open Babel descriptors
        for key, value in mol.calcdesc().iteritems():
            if (math.isnan(value)):
                format = ''
            else:
                format = type(value).__name__
            self.insert_property(key, 'Open Babel descriptor value', format)
        # insert Open Babel molecule attributes
        self.insert_property(
            'charge', 'Open Babel molecule attribute', type(mol.charge).__name__)
        self.insert_property(
            'exactmass', 'Open Babel molecule attribute', type(mol.exactmass).__name__)
        self.insert_property(
            'molwt', 'Open Babel molecule attribute', type(mol.molwt).__name__)
        self.insert_property(
            'spin', 'Open Babel molecule attribute', type(mol.spin).__name__)
        self.db.commit()
        
    def import_properties(self, inchikey, method_path_id, mol):
        # insert Open Babel descriptors
        for property_name, property_value in mol.calcdesc().iteritems():
            if (math.isnan(property_value)):
                continue
            self.insert_property_value(
                inchikey, '', method_path_id,
                property_name, 'Open Babel descriptor value', type(property_value).__name__,
                property_value, '')
        # insert Open Babel molecule attributes
        self.insert_property_value(
            inchikey, '', method_path_id,
            'charge', 'Open Babel molecule attribute', type(mol.charge).__name__,
            mol.charge, '')
        self.insert_property_value(
            inchikey, '', method_path_id,
            'exactmass', 'Open Babel molecule attribute', type(mol.exactmass).__name__,
            mol.exactmass, 'g/mol')
        self.insert_property_value(
            inchikey, '', method_path_id,
            'molwt', 'Open Babel descriptor value', type(mol.molwt).__name__,
            mol.molwt, 'g/mol')
        self.insert_property_value(
            inchikey, '', method_path_id,
            'spin', 'Open Babel descriptor value', type(mol.spin).__name__,
            mol.spin, '')
        self.db.commit()
