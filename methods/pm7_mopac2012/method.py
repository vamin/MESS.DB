# calculate mopac pm7 geometry and energies
# Victor Amin 2013

import math
import os
import sqlite3
import subprocess
import sys

from _methods import AbstractMethod

class Method(AbstractMethod):
    # method info
    method_name = 'pm7_mopac2012'
    method_description = ('precise pm7 geometry optimization '
                          'and energies using mopac')
    method_level = 'semiempirical'
    geop = 1
    # program info
    prog_name = 'MOPAC'
    prog_version = '2012'
    prog_url = 'http://openmopac.net/'
    
    def check_dependencies(self):
        try:
            mopac = subprocess.Popen(['MOPAC2012.exe'], stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
        except OSError:
            sys.exit('The ' + self.method_name + 
                     ' method requires MOPAC2012.exe (' + self.prog_url + 
                     ') to be installed and in PATH.')
        return True

    def execute(self, args):
        if args['parent_method_dir'] is None:
            sys.exit(('This method requires a parent path with a valid '
                      'xyz file (i.e., it cannot accept an InChI).'))
        inchikey = args['inchikey']
        method_dir = args['method_dir']
        parent_method_dir = args['parent_method_dir']
        p = args['path']
        inchikey_dir = self.get_inchikey_dir(inchikey)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, method_dir))
        self.setup_dir(out_dir)
        mop_file = os.path.join(out_dir, inchikey + '.mop')
        out_file = os.path.join(out_dir, inchikey + '.out')
        xyz_in = os.path.relpath(os.path.join(inchikey_dir, parent_method_dir, 
                                 inchikey + '.xyz').encode("ascii"))
        xyz_out = os.path.relpath(os.path.join(out_dir, 
                                  inchikey + '.xyz').encode("ascii"))
        if not self.check(out_file, xyz_out):
            keywords = ('PM7 PRECISE LARGE ALLVEC BONDS LOCALIZE AUX LET MMOK '
                        'DDMIN=0.0 T=1W')
            subprocess.Popen(['obabel', '-ixyz', xyz_in, '-omop', 
                              '-xk' + keywords], stdout=open(mop_file, 'w'), 
                             stderr=open(os.devnull, 'w'))
            pwd = os.getcwd()
            os.chdir(out_dir) # mopac unhappy if not run in same dir as input
            subprocess.Popen(['MOPAC2012.exe', inchikey + '.mop'])
            os.chdir(pwd)
            #subprocess.call(['obabel', '-imoo', os.path.relpath(os.path.join(out_dir, inchikey + '.out').encode("ascii")), '-oxyz'], stdout=open(xyz_out, 'w'))
            # Open Babel prints input geometry instead of 
            # the optimized for some reason. Ugh.
            self.moo_to_xyz(os.path.relpath(out_file).encode("ascii")), xyz_out)
            if self.check(out_file, xyz_out):
                self.import_properties(inchikey, p.path_id, out_file)
        else:
            self.status = 'calculation skipped'
            self.import_properties(inchikey, p.path_id, out_file)
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
            try:
                if (int(atoms) == i - 1):
                    xyz_check = True
            except (UnboundLocalError, ValueError):
                # 0 lines in file, i not set
                # or something more horrible
                pass
        except IOError:
            pass
        # moo check
        try:
            tail = subprocess.Popen(['tail', '-n 2', moo_out], 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE)
            if (tail.stdout.read().strip() == '== MOPAC DONE =='):
                moo_check = True
        except:
            pass
        if (moo_check and xyz_check):
            self.status = 'PM7 calculation completed successfully'
            return True
        else:
            self.status = 'PM7 calculation or xyz conversion failed'
            return False
    
    def log(self, args, inchikey_dir):
        base_log_path = os.path.join(inchikey_dir, args['inchikey'] + '.log')
        method_log_path = os.path.join(inchikey_dir, args['method_dir'], 
                                       args['inchikey'] + '.log')
        self.add_messages_to_log(base_log_path, self.method_name, 
                                 ['status: ' + self.status])
        self.add_messages_to_log(method_log_path, self.method_name, 
                                 ['status: ' + self.status])

    def setup_parameters(self):
        self.insert_method_parameter('PM7', '')
        self.insert_method_parameter('PRECISE', '')
        self.insert_method_parameter('LARGE', '')
        self.insert_method_parameter('ALLVEC', '')
        self.insert_method_parameter('BONDS', '')
        self.insert_method_parameter('LOCALIZE', '')
        self.insert_method_parameter('AUX', '')
        self.insert_method_parameter('LET', '')
        self.insert_method_parameter('MMOK', '')
        self.insert_method_parameter('DDMIN', 0.0)
        self.insert_method_parameter('T', '1W')
    
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
            sys.exit('Expected ' + moo_out + 
                     ' to be a valid path. Check method code.')
        # insert properties into db
        try:
            self.insert_property_value(
                inchikey, '', method_path_id,
                'HEAT OF FORMATION', 'MOPAC property', 'float',
                heat_of_formation_kcal, 'kcal/mol')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'HEAT OF FORMATION', 'MOPAC property', 'float',
                heat_of_formation_kJ, 'kJ/mol')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'TOTAL ENERGY', 'MOPAC property', 'float',
                total_energy, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'ELECTRONIC ENERGY', 'MOPAC property', 'float',
                electronic_energy, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'POINT GROUP', 'MOPAC property', 'str',
                point_group, '')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'CORE-CORE REPULSION', 'MOPAC property', 'float',
                core_repulsion, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'COSMO AREA', 'MOPAC property', 'float',
                cosmo_area, 'A^2')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'COSMO VOLUME', 'MOPAC property', 'float',
                cosmo_volume, 'A^3')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'GRADIENT NORM', 'MOPAC property', 'float',
                gradient_norm, '')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'IONIZATION POTENTIAL', 'MOPAC property', 'float',
                ionization_potential, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'HOMO', 'MOPAC property', 'float',
                homo, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'LUMO', 'MOPAC property', 'float',
                lumo, 'eV')
            self.insert_property_value(
                inchikey, '', method_path_id,
                'FILLED LEVELS', 'MOPAC property', 'int',
                filled_levels, '')
        except UnboundLocalError:
            self.status = 'could not parse properties'
            return 0
        self.db.commit()

    def moo_to_xyz(self, moo, xyz):
        with open(moo) as f:
            skip = True
            xyz_coords = []
            for line in f:
                if ('COMPUTATION TIME' in line):
                    skip = False
                if (skip):
                    continue
                # parse xyz
                s = line.split()
                try:
                    if (int(s[0]) + 1 > 0): # check if first col is a number
                        xyz_coords.append((s[1], s[2], s[4], s[6]))
                except (IndexError, ValueError, TypeError):
                    pass # not a coordinate line
                if ('Empirical Formula:' in line):
                    # check integrity
                    if not (len(xyz_coords) == int(line.split()[-2])):
                        print moo + ':'
                        print 'unrecoverable error in xyz conversion.' + "\n"
                        break
                    # print xyz
                    with open(xyz, 'w') as x:
                        x.write(str(len(xyz_coords)) + "\n")
                        x.write(moo + "\n")
                        for a in xyz_coords:
                            x.write("\t".join(a) + "\n")
                    break