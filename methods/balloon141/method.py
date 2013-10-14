#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

# generate 3d structure with balloon

import binascii
import math
import os
import sqlite3
import subprocess
import sys
from datetime import datetime

import pybel

from _methods import AbstractMethod

class Method(AbstractMethod):
    # method info
    method_name = 'balloon141'
    method_description = 'generate 3d structures from InChI with balloon'
    method_level = 'mm'
    geop = 1
    # program info
    prog_name = 'Balloon'
    prog_version = '1.4.1'
    prog_url = 'http://users.abo.fi/mivainio/balloon/'
    
    def execute(self, args):
        inchikey = args['inchikey']
        method_dir = args['method_dir']
        inchikey_dir = self.get_inchikey_dir(inchikey)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, method_dir))
        self.setup_dir(out_dir)
        xyz_out = os.path.join(out_dir, inchikey + '.xyz')
        sdf_out = os.path.realpath(os.path.join(inchikey_dir, method_dir, inchikey + '.sdf'))
        balloon_stdout = ''
        balloon_stderr = ''
        if not self.check(xyz_out):
            smiles = self.c.execute('SELECT smiles FROM molecule WHERE inchikey=?', (inchikey,)).next()
            seed = binascii.crc32(inchikey) & 0xffffffff # returns positive 32-bit integer
            try:
                os.remove(sdf_out)
            except OSError:
                pass
            balloon = subprocess.Popen(['/home/victor/software/balloon_x86_64/balloon', '-v', '1', '--maxtime', '256', '--randomSeed', str(seed), '--nGenerations', '1024', '--singleconf', '--fullforce', '--stripSalts', '-n', smiles['smiles'], sdf_out], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            balloon_stdout = balloon.stdout.read()
            balloon_stderr = balloon.stderr.read()
            try:
                mol = pybel.readfile('sdf', str(sdf_out)).next()
            except IOError:
                sdf_bad = os.path.join(os.path.dirname(__file__), '../../', inchikey + '_bad.sdf')
                sdf_out = os.path.join(out_dir, inchikey + '_bad.sdf')
                os.rename(sdf_bad, sdf_out)
                mol = pybel.readfile('sdf', str(sdf_out)).next()
            mol.localopt(forcefield='uff', steps=128)
            mol.write('xyz', str(xyz_out))
            self.check(xyz_out)
        else:
            self.status = 'skipped'
        self.log(args, inchikey_dir, [balloon_stdout, balloon_stderr])
        self.db.commit()
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
    
    def log(self, args, inchikey_dir, messages):
        base_log_path = os.path.join(inchikey_dir, args['inchikey'] + '.log')
        method_log_path = os.path.join(inchikey_dir, args['method_dir'], args['inchikey'] + '.log')
        self.add_messages_to_log(base_log_path, self.method_name, ['status: ' + self.status])
        if (messages):
            self.add_messages_to_log(method_log_path, self.method_name, messages)
    
    def setup_parameters(self):
        self.insert_method_parameter('-v', '1')
        self.insert_method_parameter('--maxtime', '256')
        self.insert_method_parameter('--randomSeed', 'crc32(inchikey)')
        self.insert_method_parameter('--nGenerations', '1024')
        self.insert_method_parameter('--singleconf', '')
        self.insert_method_parameter('--stripSalts', '')
        self.insert_method_parameter('-n', '')
        self.insert_method_parameter("pybel mol.localopt(forcefield='uff', steps=128)", '')
        
    def import_properties(self):
        pass