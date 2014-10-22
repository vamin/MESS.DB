# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB reset module

This module contains the reset tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import glob
import os
import shutil
import sys

from mess.tool import AbstractTool


class Reset(AbstractTool):
    """Reset (clear) MESS.DB."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Reset (clear) MESS.DB instance'
        self.epilog = ('This tool removes the mess database, everything in '
                       'the molecules directory, and all logs. It does not '
                       'touch backups or sources.')
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('-f', '--force', action='store_true',
                               help='do not prompt before reset')
    
    def execute(self, args):
        """Remove db/mess.db, molecules/*/, and logs/*.log."""
        if not args.force:
            if not self.query_user(('Are you sure? It might be a good idea to '
                                    'backup first.')):
                sys.exit('Canceled.')
        db_files = glob.glob(os.path.join(os.path.dirname(__file__),
                                          '../../db/mess.db*'))
        mol_dirs = glob.glob(os.path.join(os.path.dirname(__file__),
                                          '../../molecules/*/'))
        log_files = glob.glob(os.path.join(os.path.dirname(__file__),
                                           '../../logs/*.log'))
        self.log_console.info('Removing db/mess.db')
        for f in db_files:
            os.remove(f)
        self.log_console.info('Removing molecules/*/')
        for d in mol_dirs:
            shutil.rmtree(d)
        self.log_console.info('Removing logs/*.log')
        for f in log_files:
            os.remove(f)
        self.log_console.info('Done')
    
    def query_user(self, question, default='no'):
        """Ask a yes/no question via raw_input() and return user answer.
        
        question: a string that is presented to the user.
        default: the presumed answer if the user just hits <Enter>.
            It must be "yes", "no" or None (meaning
            an answer is required of the user).
        
        Returns answer as a boolean, True (yes) or False (no).
        """
        valid = {"yes": True, "y": True, "ye": True,
                 "no": False, "n": False}
        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)
        while True:
            sys.stdout.write(question + prompt)
            choice = raw_input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")


def load():
    """Load Reset()."""
    return Reset()
