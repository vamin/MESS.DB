from __future__ import print_function
from __future__ import unicode_literals

import codecs
import math
import os
import subprocess
import sys

from _method import AbstractMethod
from utils import get_inchikey_dir, setup_dir

class Mopac(AbstractMethod):
    # method info
    description = ('Semiempirical geometry optimization, energies, and bond '
                   'order analysis')
    geop = 1
    # program info
    prog_name = 'MOPAC'
    prog_version = '2012'
    prog_url = 'http://openmopac.net/'
    # parameters
    parameters = {'PM7': '',
                  'PRECISE': '',
                  'LARGE': '',
                  'ALLVEC': '',
                  'BONDS': '',
                  'LOCALIZE': '',
                  'AUX': '',
                  'LET': '',
                  'MMOK': '',
                  'DDMIN': 0.0,
                  'T': '1W',
                  'CYCLES': 2048}
    tags = ['PM7']
    threads = 1
    
    def check_dependencies(self):
        """Check that MOPAC2012.exe is installed and runnable.
        
        Returns:
            True if it's possible to run MOPAC2012.exe
        
        """
        try:
            mopac = subprocess.Popen(['MOPAC2012.exe'], stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
        except OSError:
            sys.exit(('The %s method requires MOPAC2012.exe (%s) to be '
                      'installed and in PATH.') % (self.method_name, 
                                                   self.prog_url))
        return True

    def execute(self, args):
        """Run mopac and import calculated properties for molecule."""
        p = args['path']
        parent_method_dir = p.parent_method_dir
        if parent_method_dir is '':
            sys.exit(('This method requires a parent path with a valid '
                      'xyz file (i.e., it cannot accept an InChI).'))
        method_dir = p.method_dir
        inchikey = args['inchikey']
        inchikey_dir = get_inchikey_dir(inchikey)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, method_dir))
        setup_dir(out_dir)
        mop_file = os.path.join(out_dir, '%s.mop' % inchikey)
        out_file = os.path.join(out_dir, '%s.out' % inchikey)
        xyz_in = os.path.abspath(os.path.join(inchikey_dir, parent_method_dir, 
                                 '%s.xyz' % inchikey))
        xyz_out = os.path.abspath(os.path.join(out_dir, 
                                  '%s.xyz' % inchikey))
        if not self.check(out_file, xyz_out):
            keywords = ''
            for k, v in self.parameters.items():
                if v:
                    keywords += '%s=%s ' % (k, v)
                else:
                    keywords += '%s ' % k
            keywords += 'THREADS=%s ' % self.threads
            q = ('SELECT result AS charge '
                 'FROM molecule_method_property mpp '
                 'JOIN property p ON mpp.property_id = p.property_id '
                 "WHERE p.name='charge' AND mpp.inchikey=?")
            r = self.c.execute(q, (inchikey,)).fetchone()
            keywords += 'CHARGE=%i' % r.charge
            babel = subprocess.Popen(['obabel', '-ixyz', xyz_in, '-omop', 
                                      '-xk' + keywords], 
                                     stdout=codecs.open(mop_file, 'w', 
                                                        'utf-8'), 
                                     stderr=subprocess.PIPE)
            babel_stderr = babel.stderr.read()
            pwd = os.getcwd()
            os.chdir(out_dir) # mopac unhappy if not run in same dir as input
            subprocess.Popen(['MOPAC2012.exe', '%s.mop' % inchikey]).wait()
            os.chdir(pwd)
            self.moo_to_xyz(os.path.abspath(out_file), xyz_out)
            if self.check(out_file, xyz_out):
                self.import_properties(inchikey, p.path_id, out_file)
            else:
                print(babel_stderr, file=sys.stderr)
        else:
            self.status = 'calculation skipped'
            self.import_properties(inchikey, p.path_id, out_file)
        self.log(args, inchikey_dir)
        return self.status

    def check(self, moo_out, xyz_out):
        """Check that a valid Mopac output and xyz file was generated.
        
        Args:
            moo_out: Path to Mopac output file.
            moo_out: Path to xyz output file.
        
        Sets:
            self.status: Message describing whether Mopac succeeded or failed.
        
        Returns:
            True if everything is fine, False otherwise.
        
        """
        xyz_check = False
        moo_check = False
        try: # xyz check
            with codecs.open(xyz_out, encoding='utf-8') as f:
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
        """Log messages to base log and method log.
        
        Args:
            args: The args parameter passed to the method (dict expected to
                  contain 'inchikey')
            inchikey_dir: Directory of molecule.
        
        """
        base_log_path = os.path.join(inchikey_dir, '%s.log' % args['inchikey'])
        method_log_path = os.path.join(inchikey_dir, args['path'].method_dir, 
                                       '%s.log' % args['inchikey'])
        self.add_messages_to_log(base_log_path, self.method_name, 
                                 ['status: %s' % self.status])
        self.add_messages_to_log(method_log_path, self.method_name, 
                                 ['status: %s' % self.status])
    
    def import_properties(self, inchikey, method_path_id, moo_out):
        """Load properties available in Mopac output into mess.db.
        
        Args:
            inchikey: The molecule InChIKey.
            method_path_id: Path id for this Mopac calculation.
            moo_out: Path to a mopac output file.
        
        """
        try: # moo check
            with codecs.open(moo_out, encoding='utf-8') as f:
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
                        try:
                            lumo = s[6]
                        except IndexError:
                            lumo = None
                            if homo.count('.') > 1:
                                s = homo.split('.')
                                homo = s[0] + '.' + s[1]
                    if ('NO. OF FILLED LEVELS' in line):
                        filled_levels = line.split()[5]
                    if ('COMPUTATION TIME' in line):
                        # no properties after, might as well skip ahead
                        break
        except IOError:
            sys.exit('Expected %s to be a valid path. Check method code.' % 
                     moo_out)
        # insert properties into db
        try:
            self.insert_property(
                inchikey, method_path_id,
                'HEAT OF FORMATION', 'MOPAC property', 'float',
                heat_of_formation_kcal, 'kcal/mol')
            self.insert_property(
                inchikey, method_path_id,
                'HEAT OF FORMATION', 'MOPAC property', 'float',
                heat_of_formation_kJ, 'kJ/mol')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'TOTAL ENERGY', 'MOPAC property', 'float',
                total_energy, 'eV')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'ELECTRONIC ENERGY', 'MOPAC property', 'float',
                electronic_energy, 'eV')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'POINT GROUP', 'MOPAC property', 'str',
                point_group, '')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'CORE-CORE REPULSION', 'MOPAC property', 'float',
                core_repulsion, 'eV')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'COSMO AREA', 'MOPAC property', 'float',
                cosmo_area, 'A^2')
            self.insert_property(
                inchikey, method_path_id,
                'COSMO VOLUME', 'MOPAC property', 'float',
                cosmo_volume, 'A^3')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'GRADIENT NORM', 'MOPAC property', 'float',
                gradient_norm, '')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'IONIZATION POTENTIAL', 'MOPAC property', 'float',
                ionization_potential, 'eV')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'HOMO', 'MOPAC property', 'float',
                homo, 'eV')
            self.insert_property(
                inchikey, method_path_id,
                'LUMO', 'MOPAC property', 'float',
                lumo, 'eV')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        try:
            self.insert_property(
                inchikey, method_path_id,
                'FILLED LEVELS', 'MOPAC property', 'int',
                filled_levels, '')
        except UnboundLocalError as e:
            print(e, file=sys.stderr)
            self.status = 'some properties not parsed'
        self.db.commit()

    def moo_to_xyz(self, moo, xyz):
        """Convert Mopac output to xyz.
        
        Args:
            moo: Path to Mopac output file.
            xyz: Path to xyz file to print to.
        
        """
        with codecs.open(moo, encoding='utf-8') as f:
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
                        print('%s:' % moo, file=sys.stderr)
                        print('unrecoverable error in xyz conversion.\n',
                              file=sys.stderr)
                        break
                    # print xyz
                    with codecs.open(xyz, 'w', 'utf-8') as x:
                        x.write('%i\n' % len(xyz_coords))
                        x.write('%s\n' % moo)
                        for a in xyz_coords:
                            x.write('%s\n' % '\t'.join(a))
                    break


def load(db):
    """Load Mopac(db)."""
    return Mopac(db)