from __future__ import print_function
from __future__ import unicode_literals

import os
import sqlite3
import sys
from collections import namedtuple

class MessDB(object):
    def __init__(self):
        # connect to db
        try:
            self.conn = sqlite3.connect(os.path.join(os.path.dirname(__file__), 
                                        '../db/mess.db'))
        except IOError:
            sys.exit('could not find mess.db')
        self.conn.row_factory = self.namedtuple_factory
    
    def cursor(self):
        return self.conn.cursor()
    
    def commit(self):
        return self.conn.commit()

    def namedtuple_factory(self, cursor, row):
        """
        Usage:
        con.row_factory = namedtuple_factory
        """
        fields = [col[0] for col in cursor.description]
        Row = namedtuple('Row', fields)
        return Row(*row)
    
    def __del__(self):
        self.conn.commit()
        self.conn.close()