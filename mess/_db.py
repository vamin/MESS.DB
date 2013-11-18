from __future__ import print_function
from __future__ import unicode_literals

import os
import sqlite3
import sys
import traceback
from collections import namedtuple

class MessDB(object):
    """Set up a connection to mess.db."""
    def __init__(self):
        """Initialize db connection and row factory."""
        self.tries = 0
        self.max_tries = 3
        try:
            self.conn = sqlite3.connect(os.path.join(os.path.dirname(__file__), 
                                        '../db/mess.db'))
        except IOError:
            sys.exit('could not find mess.db')
        self.conn.row_factory = self.namedtuple_factory

    def close(self):
        """Close the db connection."""
        return self.conn.close()

    def commit(self):
        """Commit transactions to db."""
        while self.tries < self.max_tries:
            try:
                return self.conn.commit()
            except sqlite3.OperationalError:
                print('Database is locked, trying again.', file=sys.stderr)
                self.tries += 1
        traceback.print_stack()
        sys.exit('Database is still locked after %d tries.' % self.tries)

    def cursor(self):
        """Return cursor."""
        return self.conn.cursor()

    def namedtuple_factory(self, cursor, row):
        """
        Usage:
            con.row_factory = namedtuple_factory

        """
        fields = [col[0] for col in cursor.description]
        Row = namedtuple('Row', fields)
        return Row(*row)
    
    def __del__(self):
        """On deletion, commit transations and close the connection."""
        self.conn.commit()
        self.conn.close()