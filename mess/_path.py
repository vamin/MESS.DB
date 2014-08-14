# -*- coding: utf-8 -*-
"""MESS.DB path module

This module contains the Path class, which describes the relationships between
methods. Compute methods often take the outputs of other methods as input. The
parent-child relationships between methods for a directed acyclic graph, which
is described in mess.db and maintained by a Path object.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys


class Path(object):
    """This class describes the relationships between compute methods, whose
    parent-child relationships form a directed acyclic graph which is described
    in mess.db and maintained by the methods in this class.
    
    Attributes:
        db (obj): A MessDB object
        cur (obj): db.cursor()
        path_id (int): The method_path_id from the method_path table in mess.db
        method_dir (str): The name of the molecule subdirectory that contains
                          the results of the method calculation
        parent_method_dir (str): The name of the molecule subdirectory that
                                 contains the inputs to the method calculation
    """
    def __init__(self, db):
        """Initialize db cursor.
        
        Args:
            db (obj): A MessDB object
        """
        self.db = db
    
    def setup(self, method_id, parent_path_id=None):
        """Setup path in mess.db."""
        method = self.setup_method(method_id)
        if parent_path_id:
            (parent_method,
             superparent_method,
             path_length) = self.setup_parent_path(parent_path_id)
        elif method['name'] == 'import':
            query = ('SELECT method_path_id FROM method_path_parent '
                     'WHERE method_id=? '
                     'AND parent_method_path_id=method_path_id')
            result = self.db.execute(query, (method['id'], )).fetchone()
            try:
                parent_path_id = result.method_path_id
            except AttributeError:
                pass
            parent_method = method
            superparent_method = {}
            path_length = 0
        else:  # assume parent is import
            query = 'SELECT method_path_id FROM method_path WHERE length=?'
            result = self.db.execute(query, (0, )).fetchone()
            parent_path_id = result.method_path_id
            query = 'SELECT method_id, name FROM method WHERE name = ?'
            result = self.db.execute(query, ('import',)).fetchone()
            parent_method = self.setup_method(result.method_id)
            superparent_method = {}
            path_length = 1
        # check if path exists, add if not
        self.path_id = self.insert_path(method, parent_method, parent_path_id,
                                        path_length)
        # set dir
        self.method_dir = self.get_dir(method, parent_method, self.path_id)
        # set parent dir
        self.parent_method_dir = self.get_dir(parent_method,
                                              superparent_method,
                                              parent_path_id)
        self.db.close()
    
    def setup_method(self, method_id):
        """Set up method dict.
        Args:
            method_id: The method id.
        
        Returns:
            A dict containing the id, name, hash, and tags for the method.
        """
        query = 'SELECT name, hash FROM method WHERE method_id = ?'
        result = self.db.execute(query, (method_id, )).fetchone()
        method = {}
        method['id'] = method_id
        method['name'] = result.name
        method['hash'] = result.hash
        query = ('SELECT p.name, mp.setting FROM parameter p '
                 'JOIN method_parameter mp '
                 'ON p.parameter_id = mp.parameter_id '
                 'JOIN method_tag mt ON p.parameter_id = mt.parameter_id '
                 'WHERE mt.method_id = ?')
        method['tags'] = []
        for result in self.db.execute(query, (method_id, )).fetchall():
            if result.setting:
                method['tags'].append(result.name + '=' + result.setting)
            else:
                method['tags'].append(result.name)
        method['tags'].sort()
        return method
    
    def setup_parent_path(self, parent_path_id):
        """Get parent method, superparent method, and path length.
        
        Args:
            parent_path_id: The parent path id.
        
        Returns:
            A tuple containing a parent method dict, a superparent method
            dict, and a path length.
        """
        query = ('SELECT mp.length, me.parent_method_id, me.child_method_id, '
                 'm.hash hash '
                 'FROM method_path mp '
                 'JOIN method_path_edge mpe ON mpe.distance = mp.length AND '
                 'mpe.method_path_id = mp.method_path_id '
                 'JOIN method_edge me '
                 'ON me.method_edge_id = mpe.method_edge_id '
                 'JOIN method m ON me.child_method_id = m.method_id '
                 'WHERE mp.method_path_id= ?')
        result = self.db.execute(query, (parent_path_id, )).fetchone()
        try:
            path_length = result.length + 1
            parent_method = self.setup_method(result.child_method_id)
            superparent_method = self.setup_method(result.parent_method_id)
            return (parent_method, superparent_method, path_length)
        except AttributeError:
            sys.exit(('%s is an invalid path id (i.e., it does not have a '
                      'valid record in the database).') % parent_path_id)
    
    def insert_edges(self, child_method_id, parent_method_id):
        """Insert edge between vertecies (methods) in mess.db."""
        query = ('INSERT OR IGNORE INTO method_edge '
                 '(parent_method_id, child_method_id) '
                 'SELECT parent_method_id, ? FROM method_edge '
                 'WHERE child_method_id = ? UNION ALL SELECT ?, ?;')
        self.db.execute(query, (child_method_id, parent_method_id,
                                child_method_id, child_method_id))
    
    def populate_edges(self, method_id, parent_method_id, path_id,
                       parent_path_id, path_length):
        """Specify a path using edges from DAG described in method_edge.
        
        Args:
            method_id: Method id of the endpoint method.
            parent_method_id: Method id of the parent to the endpoint.
            path_id: Path id for the entire path.
            parent_path_id: Path id for the path up to the parent method.
            path_length: The length of the entire path.
        """
        query = ('INSERT OR IGNORE INTO method_path_edge '
                 '(method_path_id, method_edge_id, distance) '
                 'SELECT ?, method_edge_id, distance '
                 'FROM method_path_edge WHERE method_path_id = ? '
                 'UNION ALL SELECT ?, method_edge_id, ? '
                 'FROM method_edge '
                 'WHERE parent_method_id = ? AND child_method_id = ? '
                 'UNION ALL SELECT ?, method_edge_id, distance '
                 'FROM method_path_edge WHERE method_path_id = ?')
        self.db.execute(query, (path_id, path_id, path_id, path_length,
                                parent_method_id, method_id, path_id,
                                parent_path_id))
    
    def insert_path(self, method, parent_method, parent_path_id, path_length):
        """Insert a new path into mess.db.
        
        Args:
            method: An endpoint method tuple (e.g. from setup_method).
            parent_method: Method tuple for the parent of the endpoint method.
            parent_path_id: Path id for the path up to the parent method.
            path_length: The length of the entire path.
        """
        # check if path exists
        query = ('SELECT method_path_id FROM method_path_parent '
                 'WHERE method_id=? AND parent_method_path_id=?')
        result = self.db.execute(query,
                                 (method['id'], parent_path_id)).fetchone()
        try:
            path_id = result.method_path_id
        except AttributeError:
            # insert new path
            query = 'INSERT INTO method_path (length) VALUES (?);'
            cur = self.db.execute(query, (path_length,))
            path_id = cur.lastrowid
            if not parent_path_id:
                parent_path_id = path_id
            query = ('INSERT INTO method_path_parent '
                     '(method_id, parent_method_path_id, method_path_id) '
                     'VALUES (?, ?, ?);')
            self.db.execute(query, (method['id'], parent_path_id, path_id))
            # insert edges
            self.insert_edges(method['id'], parent_method['id'])
            # populate path edges
            self.populate_edges(method['id'], parent_method['id'], path_id,
                                parent_path_id, path_length)
        return path_id
    
    def get_dir(self, method, parent_method, path_id):
        """Generate a directory name for a path.
        
        Args:
            method: An endpoint method tuple (e.g. from setup_method)
            parent_method: Method tuple for the parent of the endpoint method.
            path_id: Path id for the entire path.
        
        Returns:
            Directory name string of the form
            'mname_mtags_FROM_parentmname_parentmtags_PATH_pid'.
        """
        try:
            directory = method['name']
            if method['tags']:
                directory += '_' + '_'.join(method['tags'])
            directory += ('_' + method['hash'][:7]
                          + '_FROM_' + parent_method['name'])
            if parent_method['tags']:
                directory += '_'.join(parent_method['tags'])
            directory += ('_' + parent_method['hash'][:7]
                          + '_PATH_' + str(path_id))
            return directory
        except KeyError:
            return ''
