# generate 3d structure with balloon
# Victor Amin 2013

from __future__ import print_function
from __future__ import unicode_literals

import binascii
import codecs
import math
import os
import sqlite3
import subprocess
import sys
from datetime import datetime
from distutils.version import LooseVersion

import pybel

from _decorators import decorate, UnicodeDecorator
from _methods import AbstractMethod
from _utils import get_inchikey_dir, setup_dir

decorate(pybel, UnicodeDecorator)

class Balloon(AbstractMethod):
    # method info
    description = 'generate 3d structures from InChI with balloon'
    geop = 1
    # program info
    prog_name = 'Balloon'
    prog_version = '' # set dynamically by property method
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
        try:
            balloon = subprocess.Popen(['balloon'], stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
            return balloon.stdout.read().split()[2]
        except OSError:
            sys.exit('The ' + self.method_name + ' method requires balloon (' + 
                     self.prog_url + ').')

    def check_dependencies(self):
        # setting prog_version checks for Balloon
        return True
    
    def execute(self, args):
        inchikey = args['inchikey']
        method_dir = args['path'].method_dir
        inchikey_dir = get_inchikey_dir(inchikey)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, method_dir))
        setup_dir(out_dir)
        xyz_out = os.path.join(out_dir, inchikey + '.xyz')
        sdf_out = os.path.realpath(os.path.join(inchikey_dir, method_dir, 
                                                inchikey + '.sdf'))
        balloon_stdout = ''
        balloon_stderr = ''
        pwd = os.getcwd()
        os.chdir(os.path.join(os.path.dirname(__file__), '../../molecules'))
        if not self.check(xyz_out):
            q = 'SELECT smiles FROM molecule WHERE inchikey=?'
            row = self.c.execute(q, (inchikey,)).next()
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
            balloon_cmd += ['--randomSeed', str(seed), row.smiles, sdf_out]
            balloon = subprocess.Popen(balloon_cmd, 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
            balloon_stdout = balloon.stdout.read()
            balloon_stderr = balloon.stderr.read()
            try:
                mol = pybel.readfile('sdf', sdf_out).next()
            except IOError:
                sdf_bad = os.path.join(os.path.dirname(__file__), 
                                       '../../molecules', inchikey + '_bad.sdf')
                sdf_out = os.path.join(out_dir, inchikey + '_bad.sdf')
                os.rename(sdf_bad, sdf_out)
                mol = pybel.readfile('sdf', sdf_out).next()
            decorate(mol, UnicodeDecorator)
            mol.localopt(forcefield='mmff94s', steps=128)
            mol.write('xyz', xyz_out)
            self.check(xyz_out)
        else:
            self.status = 'skipped'
        self.log(args, inchikey_dir, [balloon_stdout, balloon_stderr])
        self.db.commit()
        os.chdir(pwd)
        return self.status

    def check(self, xyz_out):
        try:
            with codecs.open(xyz_out, encoding='utf-8') as f:
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
        method_log_path = os.path.join(inchikey_dir, args['path'].method_dir, 
                                       args['inchikey'] + '.log')
        self.add_messages_to_log(base_log_path, self.method_name, 
                                 ['status: ' + self.status])
        if (messages):
            self.add_messages_to_log(method_log_path, self.method_name, 
                                     messages)
    
    def import_properties(self):
        pass


def load(db):
    # loads the current plugin
    return Balloon(db)