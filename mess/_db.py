from __future__ import print_function
from __future__ import unicode_literals

import codecs
import os
import sqlite3
import sys
import traceback
from collections import namedtuple
from distutils.version import LooseVersion

class MessDB(object):
    """Manage a connection to mess.db."""
    def __init__(self, *args, **kwargs):
        """Initialize attributes and db connection."""
        self.args = args
        self.kwargs = kwargs
        self.tries = 0
        self.max_tries = 3
        self.check_version()
        self.open()

    def check_version(self):
        if (LooseVersion(sqlite3.sqlite_version) < LooseVersion('3.7.0')):
            sys.exit(('The database requres SQLite version 3.7 or greater '
                      'due to the use of WAL mode.'))

    def open(self):
        """Open a connection to db/mess.db and set up row factory."""
        try:
            self.conn = sqlite3.connect(os.path.join(os.path.dirname(__file__),
                                                     '../db/mess.db'),
                                        timeout=120, *self.args, **self.kwargs)
        except IOError:
            sys.exit('could not find/create mess.db')
        self.conn.row_factory = self.namedtuple_factory
        if not self.check():
            self.initialize()
    
    def reopen(self):
        """Commit, close, and open db connection."""
        self.commit()
        self.close()
        self.open()
    
    def close(self):
        """Close the db connection."""
        return self.conn.close()
    
    def cursor(self):
        """Return cursor."""
        return self.conn.cursor()
    
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
    
    def namedtuple_factory(self, cursor, row):
        """
        Usage:
            conn.row_factory = namedtuple_factory
        
        """
        fields = [col[0] for col in cursor.description]
        Row = namedtuple('Row', fields)
        return Row(*row)
    
    def check(self):
        """Check that mess.db has the right number of tables in it."""
        q = 'SELECT count(*) count FROM sqlite_master WHERE type=?'
        r = self.cursor().execute(q, ('table',)).next()
        if (r.count != 15):
            return False
        return True
    
    def initialize(self):
        """Load the mess.db schema."""
        c = self.cursor()
        schema = os.path.join(os.path.dirname(__file__), '../db/schema.sql')
        c.executescript(codecs.open(schema, encoding='utf-8').read())
        r = c.execute('PRAGMA journal_mode=wal').next()
        if not r.journal_mode == 'wal':
            sys.exit(('Setting jounral mode to WAL failed. Run PRAGMA '
                      'journal_mode=wal in SQLite before continuing.'))
        print('New mess.db initialized.', file=sys.stderr)
    
    def __del__(self):
        """On deletion, commit transations and close the connection."""
        try:
            self.commit()
            self.close()
        except AttributeError:
            pass # there was no db connection