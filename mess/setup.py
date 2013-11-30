from __future__ import print_function
from __future__ import unicode_literals

import argparse
import codecs
import os
import sqlite3
import sys

from _tools import AbstractTool

class Setup(AbstractTool):
    """Setup MESS.DB."""
    def __init__(self):
        """Set description of tool."""
        self.description = 'Initialize MESS.DB for the first time.',
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('-rm', '--remove', action='store_true', 
                               help=('remove mess.db if it exists, '
                                     'destroying all data'))
    
    def check_dependencies(self):
        """Returns True, no dependencies to check."""
        return True
    
    def execute(self, args):
        """Create mess.db."""
        mess_db_path = os.path.realpath(os.path.join(os.path.dirname(__file__), 
                                    '../../db/mess.db'))
        if (args.remove):
            try:
                os.remove(mess_db_path)
            except OSError:
                print('%s not deleted because it already does not exist' % 
                      mess_db_path, file=sys.stderr)
        mess_db_conn = sqlite3.connect(mess_db_path)
        c = mess_db_conn.cursor()
        c.executescript(codecs.open(os.path.join(os.path.dirname(__file__), 
                                             '../../db/schema.sql'), 
                                encoding='utf-8').read())  


def load():
    """Load Setup()."""
    return Setup()