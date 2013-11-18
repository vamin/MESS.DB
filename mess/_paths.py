from __future__ import print_function
from __future__ import unicode_literals

import sqlite3
import sys

class Path(object):
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()

    def setup(self, method_id, parent_path_id = None):
        method = self.setup_method(method_id)
        if (parent_path_id):
            (parent_method, 
             superparent_method,
             path_length) = self.setup_parent_path(parent_path_id)
        elif (method['name'] == 'import'):
            q = ('SELECT method_path_id FROM method_path_parent '
                 'WHERE method_id=? AND parent_method_path_id=method_path_id')
            r = self.c.execute(q, (method['id'], )).fetchone()
            try:
                parent_path_id = r.method_path_id
            except AttributeError:
                pass
            parent_method = method
            superparent_method = {}
            path_length = 0
        else: # assume parent is import
            q = 'SELECT method_path_id FROM method_path WHERE length=?'
            r = self.c.execute(q, (0, )).fetchone()
            parent_path_id = r.method_path_id
            q = 'SELECT method_id, name FROM method WHERE name = ?'
            r = self.c.execute(q, ('import',)).fetchone()
            parent_method = self.setup_method(r.method_id)
            superparent_method = {}
            path_length = 1
        # check if path exists, add if not
        self.path_id = self.insert_path(method, parent_method, parent_path_id, 
                                   path_length)
        self.db.commit()
        # set dir
        self.method_dir = self.get_dir(method, parent_method, self.path_id)
        # set parent dir
        self.parent_method_dir = self.get_dir(parent_method, 
                                              superparent_method,
                                              parent_path_id)

    def setup_method(self, method_id):
        q = 'SELECT name, hash FROM method WHERE method_id = ?'
        r = self.c.execute(q, (method_id, )).fetchone()
        method = {}
        method['id'] = method_id
        method['name'] = r.name
        method['hash'] = r.hash
        q = ('SELECT p.name, mp.setting FROM parameter p '
             'JOIN method_parameter mp ON p.parameter_id = mp.parameter_id '
             'JOIN method_tag mt ON p.parameter_id = mt.parameter_id '
             'WHERE mt.method_id = ?')
        method['tags'] = []
        for r in self.c.execute(q, (method_id, )).fetchall():
            if (r.setting):
                method['tags'].append(r.name + '=' + r.setting)
            else:
                method['tags'].append(r.name)
        method['tags'].sort()
        return method

    def setup_parent_path(self, parent_path_id):
        q = ('SELECT mp.length, me.parent_method_id, me.child_method_id, '
             'm.hash hash '
             'FROM method_path mp '
             'JOIN method_path_edge mpe ON mpe.distance = mp.length AND '
             'mpe.method_path_id = mp.method_path_id '
             'JOIN method_edge me ON me.method_edge_id = mpe.method_edge_id '
             'JOIN method m ON me.child_method_id = m.method_id '
             'WHERE mp.method_path_id= ?')
        r = self.c.execute(q, (parent_path_id, )).fetchone()
        try:
            path_length = r.length + 1
            parent_method = self.setup_method(r.child_method_id)
            superparent_method = self.setup_method(r.parent_method_id)
            return (parent_method, superparent_method, path_length)
        except AttributeError:
            sys.exit(parent_path_id + (' is an invalid path id (i.e., '
                                           'it does not have a valid record '
                                           'in the database).'))
    
    def insert_edges(self, child_method_id, parent_method_id):
        q = ('INSERT OR IGNORE INTO method_edge '
             '(parent_method_id, child_method_id) '
             'SELECT parent_method_id, ? FROM method_edge '
             'WHERE child_method_id = ? UNION ALL SELECT ?, ?;')
        return self.c.execute(q, (child_method_id, parent_method_id, 
                                  child_method_id, child_method_id))

    def populate_edges(self, method_id, parent_method_id, path_id, 
                       parent_path_id, path_length):
        q = ('INSERT OR IGNORE INTO method_path_edge '
             '(method_path_id, method_edge_id, distance) '
             'SELECT ?, method_edge_id, distance '
             'FROM method_path_edge WHERE method_path_id = ? '
             'UNION ALL SELECT ?, method_edge_id, ? '
             'FROM method_edge '
             'WHERE parent_method_id = ? AND child_method_id = ? '
             'UNION ALL SELECT ?, method_edge_id, distance '
             'FROM method_path_edge WHERE method_path_id = ?')
        return self.c.execute(q, (path_id, path_id, path_id, path_length, 
                                  parent_method_id, method_id, path_id, 
                                  parent_path_id))

    def insert_path(self, method, parent_method, parent_path_id, path_length):
        # check if path exists
        q = ('SELECT method_path_id FROM method_path_parent '
             'WHERE method_id=? AND parent_method_path_id=?')
        r = self.c.execute(q, (method['id'], parent_path_id)).fetchone()
        try:
            path_id = r.method_path_id
        except AttributeError:
            # insert new path
            q = 'INSERT INTO method_path (length) VALUES (?);'
            self.c.execute(q, (path_length,))
            path_id = self.c.lastrowid
            if not parent_path_id:
                parent_path_id = path_id
            q = ('INSERT INTO method_path_parent '
                 '(method_id, parent_method_path_id, method_path_id) '
                 'VALUES (?, ?, ?);')
            self.c.execute(q, (method['id'], parent_path_id, path_id))
            # insert edges
            self.insert_edges(method['id'], parent_method['id'])
            # populate path edges
            self.populate_edges(method['id'], parent_method['id'], path_id,
                                parent_path_id, path_length)
        return path_id

    def get_dir(self, method, parent_method, path_id):
        try:
            d = method['name']
            if (method['tags']):
                d += '_' + '_'.join(method['tags'])
            d += '_' + method['hash'][:7] + '_FROM_' + parent_method['name']
            if (parent_method['tags']):
                d += '_'.join(parent_method['tags'])
            d += '_' + parent_method['hash'][:7] + '_PATH_' + str(path_id)
            return d
        except KeyError:
            return ''