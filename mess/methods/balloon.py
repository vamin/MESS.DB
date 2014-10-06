# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB balloon method

This module contains the balloon method class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import binascii
import codecs
import math
import os
import subprocess
import sys
import time
from distutils.version import LooseVersion

import pybel

from decorators import decorate, UnicodeDecorator
from method import AbstractMethod
from utils import setup_dir


class Balloon(AbstractMethod):
    # method info
    description = 'Generate 3d structures from InChI with balloon'
    geop = 1
    # program info
    prog_name = 'Balloon'
    prog_version = ''  # set dynamically by property method
    prog_url = 'http://users.abo.fi/mivainio/balloon/'
    # parameters
    parameters = {'-v': '1',
                  '--maxtime': '1024',
                  '--nGenerations': '1024',
                  '--singleconf': '',
                  '--randomSeed': '>>>crc32(inchikey)',
                  ">>>pybel mol.localopt(forcefield='uff', steps=128)": ''}
    tags = []
    
    @property
    def prog_version(self):
        """Get prog_version from call to balloon."""
        try:
            balloon = subprocess.Popen(['balloon'], stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            return balloon.stdout.read().split()[2]
        except OSError:
            sys.exit('The %s method requires Balloon (%s).' %\
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
                if k.startswith('>>>') or v.startswith('>>>'):
                    continue
                balloon_cmd.append(k)
                if (v):
                    balloon_cmd.append(v)
            balloon_cmd.extend(['--randomSeed', str(seed), r.smiles, sdf_out])
            balloon = subprocess.Popen(balloon_cmd, cwd=out_dir,
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            balloon.stdin.write('Y')  # in case balloon asks about overwrite
            messages.append(balloon.stdout.read())
            messages.append(balloon.stderr.read())
            try:
                mol = pybel.readfile('sdf', sdf_out).next()
            except IOError:
                sdf_bad = os.path.join(out_dir, '%s_bad.sdf' % inchikey)
                mol = pybel.readfile('sdf', sdf_bad).next()
            decorate(mol, UnicodeDecorator)
            mol.localopt(forcefield='mmff94s', steps=128)
            mol.write('xyz', xyz_out)
            if self.check(xyz_out):
                self.log_all.info('%s 3D coordinates generated',
                                      self.inchikey)
            else:
                self.log_all.warning('%s coordinate generation failed',
                                         self.inchikey)
            yield self.get_timing_query(inchikey, self.path_id, start)
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
                    if (i == 0):
                        atoms = l.strip()
            lines = i + 1
            if (int(atoms) == lines - 2):
                return True
            else:
                return False
        except IOError:
            return False


def load():
    """Load Balloon()."""
    return Balloon()
