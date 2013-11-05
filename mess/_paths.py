from __future__ import print_function
from __future__ import unicode_literals

import sqlite3
import sys

class Path(object):
    # After calling setup, the following attributes will be set:
    # method_id
    # path_id
    # path_length
    # parent_path_id
    # parent_method_id
    # parent_method_name
    # parent_level
    # superparent_method_id
    # superparent_method_name
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()
    
    def setup(self, method_id, parent_path_id = '', parent_method_name = None):
        self.method_id = method_id
        self.parent_path_id = parent_path_id
        if (parent_method_name and parent_method_name == 'import'):
            q = 'SELECT method_id, name FROM method WHERE name = ?'
            import_method_row = self.c.execute(q, 
                                               (parent_method_name,)
                                               ).fetchone()
            try:
                self.path_length = 1
                self.parent_method_name = import_method_row.name
                self.parent_method_id = import_method_row.method_id
                self.superparent_method_id = ''
            except AttributeError:
                sys.exit('It seems like import has not been run yet.')
        elif (parent_path_id):
            q = ('SELECT '
                 'mp.length, '
                 'me.parent_method_id, me.child_method_id, '
                 'm.name method_name, '
                 'l.name level_name '
                 'FROM method_path mp '
                 'JOIN method_path_edge mpe '
                 'ON mpe.distance = mp.length AND '
                 'mpe.method_path_id = mp.method_path_id '
                 'JOIN method_edge me '
                 'ON me.method_edge_id = mpe.method_edge_id '
                 'JOIN method m ON me.child_method_id = m.method_id '
                 'JOIN level l ON l.level_id = m.level_id '
                 'WHERE mp.method_path_id= ?;')
            parent_path_row = self.c.execute(q, 
                                             (self.parent_path_id,)
                                             ).fetchone()
            try:
                self.path_length = parent_path_row.length + 1
                self.parent_method_name = parent_path_row.method_name
                self.parent_method_id = parent_path_row.child_method_id
                self.superparent_method_id = parent_path_row.parent_method_id
            except AttributeError:
                sys.exit(parent_path_id + (' is an invalid path id (i.e., '
                                           'it does not have a valid record '
                                           'in the database).'))
        else: # path root (i.e. import)
            self.path_length = 0
            self.parent_method_name = ''
            self.parent_method_id = method_id
            self.superparent_method_id = ''
        # check if path exists
        q = ('SELECT method_path_id FROM method_path_parent '
             'WHERE method_id=? AND parent_method_path_id=?')
        method_path_parent_row = self.c.execute(q, 
                                                (method_id, parent_path_id)
                                                ).fetchone()
        try:
            self.path_id = method_path_parent_row.method_path_id
        except AttributeError:
            # insert new path
            q = 'INSERT INTO method_path (length) VALUES (?);'
            self.c.execute(q, (self.path_length,))
            self.path_id = self.c.lastrowid
            q = ('INSERT INTO method_path_parent '
                 '(method_path_id, parent_method_path_id, method_id) '
                 'VALUES (?, ?, ?);')
            self.c.execute(q, 
                           (self.path_id, self.parent_path_id, self.method_id))
            # insert edges
            self.insert_edges(self.method_id, self.parent_method_id)
            # populate path edges
            q = ('INSERT OR IGNORE INTO method_path_edge '
                 '(method_path_id, method_edge_id, distance) '
                 'SELECT ?, method_edge_id, distance '
                 'FROM method_path_edge WHERE method_path_id = ? '
                 'UNION ALL SELECT ?, method_edge_id, ? '
                 'FROM method_edge '
                 'WHERE parent_method_id = ? AND child_method_id = ?')
            self.c.execute(q, (self.path_id, self.path_id, self.path_id, 
                               self.path_length, self.parent_method_id, 
                               self.method_id))
        # get parent method name and level of theory
        q = ('SELECT m.name method_name, l.name level_name '
             'FROM method m JOIN level l ON l.level_id = m.level_id '
             'WHERE m.method_id = ?')
        parent_method_row = self.c.execute(q, 
                                           (self.parent_method_id,)
                                           ).fetchone()
        q = 'SELECT name FROM method WHERE method_id = ?'
        superparent_method_row = self.c.execute(q, 
                                                (self.superparent_method_id,)
                                                ).fetchone()
        try:
            self.parent_method_name = parent_method_row.method_name
            self.parent_level = parent_method_row.level_name
        except AttributeError:
            self.parent_method_name = ''
            self.parent_level = ''
        try:
            self.superparent_method_name = superparent_method_row.name
        except AttributeError:
            self.superparent_method_name = ''
        self.db.commit()
    
    def insert_edges(self, child_method_id, parent_method_id):
        q = ('INSERT OR IGNORE INTO method_edge '
             '(parent_method_id, child_method_id) '
             'SELECT parent_method_id, ? FROM method_edge '
             'WHERE child_method_id = ? UNION ALL SELECT ?, ?;')
        return self.c.execute(q, (child_method_id, parent_method_id, 
                                  child_method_id, child_method_id))