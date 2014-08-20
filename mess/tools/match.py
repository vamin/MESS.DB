from __future__ import print_function
from __future__ import unicode_literals

import argparse
import csv
import codecs
import json
import os
import sys

import numpy as np
import pybel

from mess._db import MessDB
from mess._path import MethodPath
from mess._tool import AbstractTool
from mess.utils import get_inchikey_dir


class Match(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = ('Find matches to target molecules or SMARTS '
                            'patterns')
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help='a list of inchikeys (default: STDIN)')
        subparser.add_argument('-d', '--delimiter', type=str, default='\t',
                               help=('choose a delimiter for output files '
                                     '(default: tab)'))
        subparser.add_argument('-hd', '--headers', action='store_true',
                               help=('include headers in output, not '
                                     'recommended if piping to '
                                     "'mess calculate'"))
        smarts_group = subparser.add_argument_group('smarts matching')
        smarts_group.add_argument('-m', '--smarts', type=str,
                                  help=('matches to SMARTS pattern or file '
                                        'that lists a series of SMARTS '
                                        'patterns'))
        fp_group = subparser.add_argument_group('fingerprint matching')
        fp_group.add_argument('-t', '--target', type=str,
                              help=('report fingerprint distance relative to '
                                    'target'))
        fp_group.add_argument('-c', '--cutoff', type=float, default=0.0,
                              help=('minimum similarity distance to report '
                                    'match (default: 0.0)'))
        fp_group.add_argument('-f', '--fingerprint', type=str,
                              choices=[b'FP2', b'FP3', b'FP4',
                                       b'MACCS', b'MNA', b'MPD'],
                              help=('compare fingerprint'))
        fp_group.add_argument('-s', '--spectrophore', action='store_true',
                              help=('compare Spectrophore '
                                    'descriptor/fingerprint; '
                                    'requires 3D geometry (i.e., you must '
                                    'set a path to a method that has '
                                    'generated xyz coordinates)'))
        fp_group.add_argument('-p', '--path', type=int, default=None,
                              help=('path id to 3D geometry calc, '
                                    'only used for Spectrophore'))
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
        # parse args
        if not (args.smarts or args.fingerprint or args.spectrophore):
            sys.exit('No operations were selected, nothing to match.')
        if sum(bool(arg) for arg in (args.smarts,
                                     args.fingerprint,
                                     args.spectrophore)) > 1:
            sys.exit(('One thing at a time, please. The arguments --smarts, '
                      '--fingerprint, and --spectrophore are mutually '
                      'exclusive.'))
        if args.smarts and args.target:
            print('Warning: --target ignored, proceeding with SMARTS matching',
                  file=sys.stderr)
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
        # load target and target fingerprints
        target_mol = None
        target_fp = None
        target_sp = None
        if args.target:
            if os.path.exists(args.target):
                target_mol = pybel.readfile(args.target.split('.')[-1],
                                            args.target).next()
            else:
                target_mol = pybel.readstring('smi', args.target)
        if target_mol is not None:
            if args.fingerprint:
                target_fp = self.calculate_fingerprint(target_mol,
                                                       args.fingerprint)
            if args.spectrophore:
                target_sp = self.calculate_spectrophore(target_mol, sp_args)
        # match every input
        db = MessDB()
        inchi_query = 'SELECT inchi FROM molecule WHERE inchikey = ?'
        fp_query = ('SELECT fingerprint FROM molecule_fingerprint '
                    'WHERE inchikey = ? AND name = ? '
                    'AND settings = ? AND method_path_id = ?')
        writer = csv.writer(sys.stdout, delimiter=args.delimiter)
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            if args.smarts or args.fingerprint:
                inchi = db.execute(inchi_query, (inchikey,)).fetchone()[0]
                mol = pybel.readstring('inchi', 'InChI=%s' % inchi)
            if args.smarts:
                canonical = pybel.ob.OBOp.FindType(b"canonical")
                canonical.Do(mol.OBMol)
                for (smarts_obj,
                     smarts_str) in self.smarts_generator(args.smarts):
                    matches = [match for match in smarts_obj.findall(mol)]
                    if len(matches) > 0:
                        writer.writerow([inchikey, smarts_str] + matches)
            if args.fingerprint:
                try:
                    fp = db.execute(fp_query, (inchikey, args.fingerprint,
                                               '', '')).fetchone()[0]
                except TypeError:
                    fp = self.calculate_fingerprint(mol, args.fingerprint)
                if target_fp is not None:
                    similarity = self.calculate_similarity(target_fp, fp,
                                                           'tanimoto')
                    writer.writerow([inchikey, args.fingerprint,
                                     args.target, similarity])
                else:
                    writer.writerow([inchikey, args.fingerprint] + fp)
            if args.spectrophore:
                try:
                    sp = db.execute(fp_query, (inchikey, 'Spectrophore',
                                               json.dumps(sp_args,
                                                          sort_keys=True),
                                               args.path)).fetchone()[0]
                except TypeError:
                    xyz_file = os.path.join(get_inchikey_dir(inchikey),
                                            method_dir,
                                            '%s.xyz' % inchikey)
                    mol = pybel.readfile('xyz', xyz_file).next()
                    sp = Match.calculate_spectrophore(mol, sp_args)
                if target_sp is not None:
                    similarity = self.calculate_similarity(target_sp, sp,
                                                           'cos')
                    writer.writerow([inchikey, 'Spectrophore',
                                     args.target, similarity])
                else:
                    writer.writerow([inchikey, 'Spectrophore'] + sp)
    
    @classmethod
    def calculate_fingerprint(cls, mol, fp='FP2'):
        if fp.upper() == 'FP2':
            return mol.calcfp('FP2').bits
        elif fp.upper() == 'FP3':
            return mol.calcfp('FP3').bits
        elif fp.upper() == 'FP4':
            return mol.calcfp('FP4').bits
        elif fp.upper().startswith('MAC'):  # MACCS
            return mol.calcfp('MACCS').bits
        elif fp.upper().startswith('MN'):  # MNA
            mol.title = ''
            mna = mol.write('MNA')
            return mna.split()[3:]
        elif fp.upper().startswith('MP'):  # MPD
            mpd = mol.write('MPD')
            return mpd.split()[1:]
        else:  # not a valid fingerprint
            sys.exit('Unrecognized fingerprint option: %s.' % fp)
    
    @classmethod
    def calculate_spectrophore(cls, mol, args=None):
        if not mol.dim == 3:
            print(('Warning: molecule representation not 3D. Spectrophore '
                   'will return vector of nan and 0.0.'), file=sys.stderr)
        spec = pybel.ob.OBSpectrophore()
        if args is not None:
            for key in args:
                if key not in ('normalization', 'accuracy',
                               'stereo', 'resolution'):
                    print(('Warning: ignored unrecognized Spectrophore option '
                           '(%s).') % key, file=sys.stderr)
            if 'normalization' in args and args['normalization']:
                normalization = args['normalization']
                if normalization == 'No':
                    spec.SetNormalization(
                        pybel.ob.OBSpectrophore().NoNormalization)
                elif normalization == 'ZeroMean':
                    spec.SetNormalization(
                        pybel.ob.OBSpectrophore().NormalizationTowardsZeroMean)
                elif normalization == 'UnitStd':
                    spec.SetNormalization(
                        pybel.ob.OBSpectrophore().NormalizationTowardsUnitStd)
                elif normalization == 'ZeroMeanAndUnitStd':
                    spec.SetNormalization(
                        pybel.ob.OBSpectrophore()
                        .NormalizationTowardsZeroMeanAndUnitStd)
                else:
                    print(('Warning: unrecognized Spectrophore normalization '
                           'argument (%s), normalization '
                           'set to default (No).') % normalization,
                          file=sys.stderr)
                    spec.SetNormalization(
                        pybel.ob.OBSpectrophore().NoNormalization)
            if 'accuracy' in args and args['accuracy']:
                accuracy = args['accuracy']
                accuracy_map = {
                    1: pybel.ob.OBSpectrophore().AngStepSize1,
                    2: pybel.ob.OBSpectrophore().AngStepSize2,
                    5: pybel.ob.OBSpectrophore().AngStepSize5,
                    10: pybel.ob.OBSpectrophore().AngStepSize10,
                    15: pybel.ob.OBSpectrophore().AngStepSize15,
                    20: pybel.ob.OBSpectrophore().AngStepSize20,
                    30: pybel.ob.OBSpectrophore().AngStepSize30,
                    36: pybel.ob.OBSpectrophore().AngStepSize36,
                    45: pybel.ob.OBSpectrophore().AngStepSize45,
                    60: pybel.ob.OBSpectrophore().AngStepSize60
                }
                if accuracy not in accuracy_map:
                    print(('Warning: unrecognized Spectrophore accuracy '
                           'argument (%s), accuracy '
                           'set to default (20).') % accuracy,
                          file=sys.stderr)
                    accuracy = 20
                spec.SetAccuracy(accuracy_map[accuracy])
            if 'stereo' in args and args['stereo']:
                if args['stereo'] == 'No':
                    spec.SetStereo(0)
                elif args['stereo'] == 'Unique':
                    spec.SetStereo(pybel.ob.OBSpectrophore()
                                   .UniqueStereoSpecificProbes)
                elif args['stereo'] == 'Mirror':
                    spec.SetStereo(pybel.ob.OBSpectrophore()
                                   .MirrorStereoSpecificProbes)
                elif args['stereo'] == 'All':
                    spec.SetStereo(pybel.ob.OBSpectrophore()
                                   .AllStereoSpecificProbes)
                else:
                    print(('Warning: unrecognized Spectrophore stereospecificity '
                           'argument (%s), cage type '
                           'set to default (No).') % args['stereo'],
                          file=sys.stderr)
                    spec.SetStereo(0)
            if 'resolution' in args and args['resolution']:
                try:
                    if float(args['resolution']) > 0:
                        spec.SetResolution(args['resolution'])
                    else:
                        raise ValueError
                except ValueError:
                    print(('Warning: invalid Spectrophore resolution '
                           'argument (%s), resolution '
                           'set to default (3.0).') % args['resolution'],
                          file=sys.stderr)
                    spec.SetResolution(3.0)
        return [element for element in spec.GetSpectrophore(mol.OBMol)]
    
    @classmethod
    def calculate_similarity(cls, target, query, method='tanimoto'):
        if (method.lower().startswith('tan')
                or method.lower().startswith('jac')):
            # Tanimoto / Jaccard similarity
            tset = set(target)
            qset = set(query)
            return (float(len(tset.intersection(qset)))
                    / float(len(tset.union(qset))))
        elif method.lower().startswith('cos'):  # cosine similarity
            u = np.array(target)
            v = np.array(query)
            return np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        else:
            sys.exit('Unrecognized similarity method: %s.' % method)
    
    @classmethod
    def smarts_generator(cls, smarts):
        smarts_strings = []
        # check if file or string
        if os.path.exists(smarts):
            # if file, import to string
            with codecs.open(smarts, 'r', 'utf-8') as smarts_file:
                for line in smarts_file:
                    line = line.strip()
                    if line == '' or line.startswith('#'):
                        continue
                    cols = line.split()
                    if cols[0] == 'SEEDCHARGE':
                        continue
                    if cols[0] == 'TRANSFORM':
                        cols.pop(0)
                    smarts_strings.append(cols[0])
        else:
            smarts_strings.append(smarts)
        try:
            smarts_string = ''
            for smarts_string in smarts_strings:
                smarts_obj = pybel.Smarts(smarts_string)
                yield smarts_obj, smarts_string
        except:
            sys.exit('The SMARTS string \'%s\' is not valid.' % smarts_string)


def load():
    """Load Match()."""
    return Match()
