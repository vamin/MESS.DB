# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB balloon method

This module contains the balloon method class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import binascii
import codecs
import os
import subprocess
import sys
import time

import pybel

from mess.method import AbstractMethod
from mess.utils import setup_dir


class Balloon(AbstractMethod):
    """This method uses balloon to generate 3D structures from 0D strings."""
    # method info
    description = 'Generate 3d structures from InChI with balloon'
    geop = 1
    # program info
    prog_name = 'Balloon'
    prog_version = ''  # set dynamically by property method
    prog_url = 'http://users.abo.fi/mivainio/balloon/'
    prog_citation = ('Mikko J. Vainio and Mark S. Johnson (2007) Generating '
                     'Conformer Ensembles Using a Multiobjective Genetic '
                     'Algorithm. Journal of Chemical Information and '
                     'Modeling, 47, 2462 - 2474.')
    # parameters
    parameters = {'--maxtime': '1024',
                  '--nGenerations': '1024',
                  '--maxFlipDistance': '32',
                  '--distanceDependent': '',
                  '--fullforce': '',
                  '--singleconf': '',
                  '--randomSeed': '#crc32(inchikey)'}
    
    @property
    def prog_version(self):
        """Get prog_version from call to balloon."""
        try:
            balloon = subprocess.Popen(['balloon'], stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            return balloon.stdout.read().split()[2]
        except OSError:
            sys.exit('The %s method requires Balloon (%s).' %
                     (self.method_name, self.prog_url))
    
    def check_dependencies(self):
        """Check that $BALLOON_FORCEFIELD is set."""
        # setting prog_version checks for balloon, so no need here
        if not os.environ.get('BALLOON_FORCEFIELD'):
            sys.exit(('You must set the $BALLOON_FORCEFIELD environment '
                      'variable to the path to MMFF94.mff.'))
        return True
    
    def map(self, inchikey, inchikey_dir):
        """Generate 3D structures with Balloon."""
        self.inchikey = inchikey
        start = time.time()
        out_dir = os.path.realpath(os.path.join(inchikey_dir, self.method_dir))
        setup_dir(out_dir)
        sdf_out = os.path.realpath(os.path.join(out_dir,
                                                '%s.sdf' % self.inchikey))
        xyz_out = os.path.join(out_dir, '%s.xyz' % self.inchikey)
        messages = []
        if not self.check(xyz_out):
            query = 'SELECT smiles FROM molecule WHERE inchikey=?'
            r = self.db.execute(query, (self.inchikey,)).next()
            # get positive 32-bit integer
            seed = binascii.crc32(inchikey) & 0xffffffff
            try:
                os.remove(sdf_out)
            except OSError:
                pass
            balloon_cmd = ['balloon']
            for k, v in self.parameters.items():
                if k.startswith('#') or v.startswith('#'):
                    continue
                balloon_cmd.append(k)
                if v:
                    balloon_cmd.append(v)
            balloon_cmd.extend(['--randomSeed', str(seed), r.smiles, sdf_out])
            balloon = subprocess.Popen(balloon_cmd, cwd=out_dir,
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            balloon.stdin.write('Y')  # in case balloon asks about overwrite
            messages.append(balloon.stdout.read())
            messages.append(balloon.stderr.read())
            forcefield = b'mmff94s'
            steps = 512
            try:
                mol = pybel.readfile('sdf', sdf_out).next()
                mol.write(b'xyz', str(xyz_out))
            except IOError:
                sdf_bad = os.path.join(out_dir, '%s_bad.sdf' % inchikey)
                try:
                    mol = pybel.readfile('sdf', sdf_bad).next()
                    mol.localopt(forcefield=forcefield, steps=steps)
                    self.log_all.info(('"bad" %s sdf cleaned up '
                                       'with %s forcefield '
                                       'and %i steps'),
                                      self.inchikey,
                                      forcefield,
                                      steps)
                    mol.write(b'xyz', str(xyz_out))
                except IOError:
                    pass
            if self.check(xyz_out):
                if abs(mol.molwt - pybel.readstring('smi',
                                                    r.smiles).molwt) > 0.001:
                    mol = pybel.readstring(b'smi', str(r.smiles))
                    mol.make3D(forcefield, steps)
                    mol.write(b'xyz', str(xyz_out), overwrite=True)
                    self.log_all.info(('%s 3D coordinates generation '
                                       'attempted by '
                                       'Open Babel rule-based algorithm '
                                       '(forcefields=%s steps=%i) instead of '
                                       'balloon due to hydrogen atom '
                                       'mismatch'),
                                      self.inchikey, forcefield, steps)
            else:
                mol = pybel.readstring(b'smi', str(r.smiles))
                mol.make3D(forcefield, steps)
                mol.write(b'xyz', str(xyz_out), overwrite=True)
                self.log_all.info(('%s 3D coordinates generation attempted by '
                                   'Open Babel rule-based algorithm '
                                   '(forcefields=%s steps=%i) instead of '
                                   'balloon due to unexpected failure'),
                                  self.inchikey, forcefield, steps)
            if self.check(xyz_out):
                self.log_all.info('%s 3D coordinates generated successfully',
                                  self.inchikey)
            else:
                self.log_all.warning('%s coordinate generation failed',
                                     self.inchikey)
                
            yield self.get_timing_query(self.inchikey, start)
        else:
            self.log_console.info('%s skipped', self.inchikey)
    
    def check(self, xyz_out):
        """Check that a valid xyz file was created.
        
        Args:
            xyz_out: Path to the xyz file generated by the balloon method.
        
        Returns:
            True if valid xyz, False otherwise.
        
        """
        try:
            with codecs.open(xyz_out, encoding='utf-8') as f:
                for i, l in enumerate(f):
                    if i == 0:
                        atoms = l.strip()
            lines = i + 1
            if int(atoms) == lines - 2:
                return True
            else:
                return False
        except IOError:
            return False


def load():
    """Load Balloon()."""
    return Balloon()
