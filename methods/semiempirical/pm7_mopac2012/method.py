#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

# calculate mopac pm7 geometry and energies

import math
import os
import sqlite3
import subprocess
import sys

import pybel

from _methods import AbstractMethod

class Method(AbstractMethod):
    def execute(self, args):
        if args['parent_method_dir'] is None:
            sys.exit('This method requires a parent path with a valid xyz file (i.e., it cannot accept an InChI).')
        inchikey_dir = self.get_inchikey_dir(args['inchikey'])
        out_dir = os.path.realpath(os.path.join(inchikey_dir, args['method_dir']))
        self.setup_dir(out_dir)
        mop_file = os.path.join(out_dir, args['inchikey'] + '.mop')
        out_file = os.path.join(out_dir, args['inchikey'] + '.out')
        xyz_in = os.path.relpath(os.path.join(inchikey_dir, args['parent_method_dir'], args['inchikey'] + '.xyz').encode("ascii"))
        xyz_out = os.path.relpath(os.path.join(out_dir, args['inchikey'] + '.xyz').encode("ascii"))
        if not self.check(out_file, xyz_out):
            subprocess.call(['obabel', '-ixyz', xyz_in, '-omop', '-xkPM7 PRECISE LARGE ALLVEC BONDS LOCALIZE AUX GNORM=0.0 T=1W'], stdout=open(mop_file, 'w'))
            subprocess.call(['/home/victor/software/mopac/MOPAC2012.32.exe', mop_file])
            subprocess.call(['obabel', '-imoo', os.path.relpath(os.path.join(out_dir, args['inchikey'] + '.out').encode("ascii")), '-oxyz'], stdout=open(xyz_out, 'w'))
            if self.check(out_file, xyz_out):
                self.import_properties(out_file)
        else:
            self.status = 'calculation skipped'
            self.import_properties(out_file)
        self.log(args, inchikey_dir)
        return self.status

    def check(self, moo_out, xyz_out):
        xyz_check = False
        moo_check = False
        try: # xyz check
            with open(xyz_out) as f:
                for i, l in enumerate(f):
                    if (i == 0):
                        atoms = l.strip()
            lines = i + 1
            if (int(atoms) == lines - 2):
                xyz_check = True
        except IOError:
            pass
        # moo check
        try:
            proc = subprocess.Popen(['tail -n 2 ' + moo_out + ' | head -n 1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            if (stdout.strip() == '== MOPAC DONE =='):
                moo_check = True
        except OSError:
            pass
        if (moo_check and xyz_check):
            self.status = 'PM7 calculation completed successfully'
            return True
        else:
            self.status = 'PM7 calculation or xyz conversion failed'
            return False
    
    def log(self, args, inchikey_dir):
        base_log_path = os.path.join(inchikey_dir, args['inchikey'] + '.log')
        method_log_path = os.path.join(inchikey_dir, args['method_dir'], args['inchikey'] + '.log')
        self.add_messages_to_log(base_log_path, args['method_name'], ['status: ' + self.status])
        self.add_messages_to_log(method_log_path, args['method_name'], ['status: ' + self.status])
    
    def setup_properties(self):
        self.insert_property(
            'HEAT OF FORMATION', 'MOPAC property', 'float')
        self.insert_property(
            'TOTAL ENERGY', 'MOPAC property', 'float')
        self.insert_property(
            'ELECTRONIC ENERGY', 'MOPAC property', 'float')
        self.insert_property(
            'POINT GROUP', 'MOPAC property', 'str')
        self.insert_property(
            'CORE-CORE REPULSION', 'MOPAC property', 'float')
        self.insert_property(
            'COSMO AREA', 'MOPAC property', 'float')
        self.insert_property(
            'COSMO VOLUME', 'MOPAC property', 'float')
        self.insert_property(
            'GRADIENT NORM', 'MOPAC property', 'float')
        self.insert_property(
            'IONIZATION POTENTIAL', 'MOPAC property', 'float')
        self.insert_property(
            'HOMO', 'MOPAC property', 'float')
        self.insert_property(
            'LUMO', 'MOPAC property', 'float')
        self.insert_property(
            'FILLED LEVELS', 'MOPAC property', 'int')
        self.db.commit()
        
    def import_properties(self, inchikey, method_path_id, moo_out):
        try: # moo check
            with open(moo_out) as f:
                # parse out properties
                for line in f:
                    if ('FINAL HEAT OF FORMATION' in line):
                        s = line.split()
                        heat_of_formation_kcal = s[5]
                        heat_of_formation_kJ = s[8]
                    if ('TOTAL ENERGY' in line):
                        total_energy = line.split()[3]
                    if ('ELECTRONIC ENERGY' in line):
                        s = line.split()
                        electronic_energy = s[3]
                        point_group = s[7]
                    if ('CORE-CORE REPULSION' in line):
                        core_repulsion = line.split()[3]
                    if ('COSMO AREA' in line):
                        cosmo_area = line.split()[3]
                    if ('COSMO VOLUME' in line):
                        cosmo_volume = line.split()[3]
                    if ('GRADIENT NORM' in line):
                        gradient_norm = line.split()[3]
                    if ('IONIZATION POTENTIAL' in line):
                        ionization_potential = line.split()[3]
                    if ('HOMO LUMO ENERGIES' in line):
                        s = line.split()
                        homo = s[5]
                        lumo = s[6]
                    if ('NO. OF FILLED LEVELS' in line):
                        filled_levels = line.split()[5]
                    if ('COMPUTATION TIME' in line):
                        # no properties after, might as well skip ahead
                        break
        except IOError:
            sys.exit('Expected ' + moo_out + ' to be a valid path. Check method code.')
        # insert properties into db
        self.insert_property_value(
            inchikey, None, method_path_id,
            'HEAT OF FORMATION', 'MOPAC property', 'float',
            heat_of_formation_kcal, 'kcal/mol')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'HEAT OF FORMATION', 'MOPAC property', 'float',
            heat_of_formation_kJ, 'kJ/mol')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'TOTAL ENERGY', 'MOPAC property', 'float',
            total_energy, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'ELECTRONIC ENERGY', 'MOPAC property', 'float',
            electronic_energy, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'POINT GROUP', 'MOPAC property', 'str',
            point_group, '')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'CORE-CORE REPULSION', 'MOPAC property', 'float',
            core_repulsion, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'COSMO AREA', 'MOPAC property', 'float',
            cosmo_area, 'A^2')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'COSMO VOLUME', 'MOPAC property', 'float',
            cosmo_volume, 'A^3')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'GRADIENT NORM', 'MOPAC property', 'float',
            gradient_norm, '')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'IONIZATION POTENTIAL', 'MOPAC property', 'float',
            ionization_potential, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'HOMO', 'MOPAC property', 'float',
            homo, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'LUMO', 'MOPAC property', 'float',
            lumo, 'eV')
        self.insert_property_value(
            inchikey, None, method_path_id,
            'FILLED LEVELS', 'MOPAC property', 'int',
            filled_levels, '')
        self.db.commit()
