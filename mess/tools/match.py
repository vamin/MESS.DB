from __future__ import print_function
from __future__ import unicode_literals

import codecs
import os
import sys
from distutils.version import LooseVersion

import pybel

from _decorators import decorate, UnicodeDecorator
from _methods import AbstractMethod
from _utils import get_inchikey_dir, setup_dir

decorate(pybel, UnicodeDecorator)

class Match(AbstractMethod):
    # method info
    description = 'Match molecules to SMARTS patterns'
    geop = 0
    # program info
    prog_name = 'Open Babel'
    prog_version = ''
    prog_url = 'http://openbabel.org/wiki/Main_Page'
    # parameters
    parameters = {}
    tags = []
        
    def check_dependencies(self):
        """Check for dependencies (Open Babel >=2.3).
        
        Returns:
            True if all dependencies are met.
        
        """
        try:
            if (LooseVersion(pybel.ob.OBReleaseVersion()) <
                LooseVersion('2.3.0')):
                sys.exit(('This tool requires Open Babel (and its python '
                          'module, pybel) version >=2.3.0.'))
        except AttributeError, OSError:
            sys.exit(('This tool requires Open Babel (and its python module, '
                      'pybel) version >=2.3.0.'))
        return True
    
    def execute(self, args):
        """Match molecules to SMARTS patterns."""
        inchikey = args['inchikey']
        method_dir = args['path'].method_dir
        inchikey_dir = get_inchikey_dir(inchikey)
        out_dir = os.path.realpath(os.path.join(inchikey_dir, method_dir))
        setup_dir(out_dir)
        import_smarts(args['smarts'])
        #messages = []
        mol = pybel.readfile('inchi', os.path.join(inchikey_dir, '%s.inch' % inchikey))

        #self.log(args, inchikey_dir, messages)
        return self.status
    
    def check(self):
        return True
    
    def log(self, args, inchikey_dir, messages):
        """Log messages to base log and method log.
        
        Args:
            args: The args parameter passed to the method (dict expected to
                  contain 'inchikey' and 'path')
            inchikey_dir: Directory of molecule.
            messages: List of messages to be written to logs.
        
        """
        base_log_path = os.path.join(inchikey_dir, '%s.log' % args['inchikey'])
        method_log_path = os.path.join(inchikey_dir, args['path'].method_dir,
                                       '%s.log' % args['inchikey'])
        self.add_messages_to_log(base_log_path, self.method_name,
                                 ['status: %s' % self.status])
        if (messages):
            self.add_messages_to_log(method_log_path, self.method_name,
                                     messages)
    
    def import_properties(self):
        pass

    def import_smarts(self, smarts):
        # check if file or string

        # if file, import to string
        with codecs.open(smarts, 'r', 'utf-8') as f:
            for l in f:
                l = l.strip()
                if l == '' or l.startswith('#'): continue
                cols = l.split()
                if cols[0] == 'SEEDCHARGE': continue
                if cols[0] == 'TRANSFORM': cols.pop(0)


        # load into parameters


def load(db):
    """Load Match(db)."""
    return Match(db)