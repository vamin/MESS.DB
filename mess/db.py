# -*- coding: utf-8 -*-
"""MESS.DB database module

This module contains the MessDB class, which MESS.DB uses to connect to and
interact with the 'mess.db' sqlite3 database in the 'db' directory.
"""
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
    """Manage a connection to mess.db.
    
    Attributes:
        args (list of str): Positional arguments, passed to sqlite3.connect()
        kwargs (dict): Named arguments, passed to sqlite3.connect()
        tries (int): Number of times a db commit has failed
        max_tries (int): Maximum number of times to try db commit
        conn (obj): sqlite3 db connection object
    """
    def __init__(self, **kwargs):
        """Initialize attributes and db connection.
        
        Args:
            All named args will be passed to sqlite3.connect()
        """
        self.tries = 0
        self.max_tries = 3
        self.check_version()
        self.open(**kwargs)
    
    def check_version(self):
        """Ensure that all dependencies are satisfied."""
        if LooseVersion(sqlite3.sqlite_version) < LooseVersion('3.7.0'):
            sys.exit(('The database requres SQLite version 3.7 or greater '
                      'due to the use of WAL mode.'))
    
    def open(self, **kwargs):
        """Open a connection to db/mess.db and set up row factory."""
        try:
            if 'database' not in kwargs:
                kwargs['database'] = os.path.join(os.path.dirname(__file__),
                                                  '../db/mess.db')
            if 'timeout' not in kwargs:
                kwargs['timeout'] = 120
            self.conn = sqlite3.connect(**kwargs)
        except IOError:
            sys.exit('could not find/create mess.db')
        self.conn.row_factory = self.namedtuple_factory
        if not self.check():
            self.initialize()
    
    def reopen(self):
        """Commit, close, and open db connection."""
        self.close()
        self.open()
    
    def close(self):
        """Close the db connection."""
        self.commit()
        self.conn.close()
    
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
    
    def execute(self, query, values=()):
        """Execute a single query."""
        try:
            with self.conn:
                return self.conn.execute(query, values)
        except sqlite3.Error as err:
            print(('Error executing query:\n%s\n'
                   'With values:\n%s\n'
                   '%s') % (query, values, err), file=sys.stderr)
            
    def executemany(self, query, values):
        """Execute a single query with many values."""
        try:
            with self.conn:
                return self.conn.executemany(query, values)
        except sqlite3.Error as err:
            print(('Error executing query:\n%s\n'
                   'With many values.'
                   '%s') % (query, err), file=sys.stderr)
    
    def executescript(self, script):
        """Read and execute queries from file."""
        try:
            with self.conn:
                return self.conn.executescript(script)
        except sqlite3.Error as err:
            print(('Error executing script:\n%s\n'
                   '%s') % (script, err), file=sys.stderr)
    
    def namedtuple_factory(self, cursor, row):
        """
        Usage:
            conn.row_factory = namedtuple_factory
        """
        fields = [col[0] for col in cursor.description]
        row_tuple = namedtuple('Row', fields)
        return row_tuple(*row)
    
    def check(self):
        """Check that mess.db has the right number of tables in it."""
        query = 'SELECT count(*) count FROM sqlite_master WHERE type=?'
        result = self.cursor().execute(query, ('table',)).next()
        if result.count != 15:
            return False
        return True
    
    def initialize(self):
        """Load the mess.db schema."""
        tables = os.path.join(os.path.dirname(__file__),
                              '../db/schema/tables.sql')
        views = os.path.join(os.path.dirname(__file__),
                             '../db/schema/views.sql')
        triggers = os.path.join(os.path.dirname(__file__),
                                '../db/schema/triggers.sql')
        self.executescript(codecs.open(tables, encoding='utf-8').read())
        self.executescript(codecs.open(views, encoding='utf-8').read())
        self.executescript(codecs.open(triggers, encoding='utf-8').read())
        result = self.execute('PRAGMA journal_mode=wal').next()
        if not result.journal_mode == 'wal':
            sys.exit(('Setting journal mode to WAL failed. Run PRAGMA '
                      'journal_mode=wal in SQLite before continuing.'))
        print('New mess.db initialized.', file=sys.stderr)
    
    @property
    def total_changes(self):
        """Returns the total number of database rows that have been modified,
        inserted, or deleted since the db connection was opened."""
        return self.conn.total_changes
    
    def __del__(self):
        """On deletion, commit transations and close the connection."""
        try:
            self.close()
        except AttributeError:
            pass  # there was no db connection
