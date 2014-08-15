from __future__ import print_function
from __future__ import unicode_literals

import argparse
import csv
import codecs
import os
import sys

import pybel

from _db import MessDB
from _tool import AbstractTool


class Match(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Find matches to SMARTS patterns.'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help=('A list of inchikeys, file or passed in '
                                     'through STDIN'))
        subparser.add_argument('-s', '--smarts', type=str,
                               help=('A SMARTS pattern or file that '
                                     'lists a series of SMARTS patterns'))
        subparser.add_argument('-d', '--delimiter', type=str, default='\t',
                               help=('Choose a delimiter for output files, '
                                     'tab is default'))
    
    def execute(self, args):
        """Match molecules to SMARTS patterns."""
        if args.inchikeys.name == '<stdin>' and args.inchikeys.isatty():
            sys.exit('No input specified.')
        db = MessDB()
        query = 'SELECT inchi FROM molecule WHERE inchikey = ?'
        writer = csv.writer(sys.stdout, delimiter=args.delimiter)
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            inchi = db.execute(query, (inchikey,)).fetchone()[0]
            mol = pybel.readstring('inchi', 'InChI=%s' % inchi)
            canonical = pybel.ob.OBOp.FindType(b"canonical")
            canonical.Do(mol.OBMol)
            for (smarts_obj,
                 smarts_str) in self.smarts_generator(args.smarts):
                matches = [match for match in smarts_obj.findall(mol)]
                if len(matches) > 0:
                    writer.writerow([inchikey, smarts_str] + matches)
    
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
