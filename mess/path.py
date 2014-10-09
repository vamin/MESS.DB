# -*- coding: utf-8 -*-
"""MESS.DB path module

This module contains the MethodPath class, which describes the relationships
between methods. Compute methods often take the outputs of other methods as
input. The parent-child relationships between methods comprise a directed
graph, which is described in mess.db and maintained by a MethodPath object.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys

from mess.db import MessDB
from mess.log import Log


class Node(object):
    """Node class defines a node and its parent/child nodes.

    Attributes:
        _id: A unique identifier for the node
        _parents (dict): parent nodes => edge identifiers
        _children (dict): child nodes => edge identifiers
    """
    def __init__(self, node_id):
        """Initialize the node

        Args:
            node_id: A unique identifier for the node
        """
        self._id = node_id
        self._parents = {}
        self._children = {}
    
    def add_parent(self, parent, edge_id):
        """Add a node to the parents dict of this node.

        Args:
            parent: A node object
            edge_id: An identifier for the parent->node edge relationship
        """
        self._parents[parent] = edge_id
    
    def add_child(self, child, edge_id):
        """Add a node to the children dict of this node.

        Args:
            child: A node object
            edge_id: An identifier for the node->child edge relationship
        """
        self._children[child] = edge_id
    
    def __str__(self):
        """Return an ascii representation of the nodes parents and children.
        """
        parents_str = '\n'.join(['%s ->' % node.get_id()
                                 for node in self._parents])
        children_str = '\n'.join(['    -> %s' % node.get_id()
                                  for node in self._children])
        return '%s\n  %s\n%s' % (parents_str, self._id, children_str)
    
    def get_parents(self):
        """Return the parent node representations as a list."""
        return self._parents.keys()
    
    def get_children(self):
        """Return the child node representations as a list."""
        return self._children.keys()
    
    def get_id(self):
        """Return this node's id."""
        return self._id
    
    def get_edge_to(self, child):
        """Return the edge identifier for a node->child relationship, or None
        if no relationship exists.

        Args:
            child: A child object
        """
        try:
            return self._children[child]
        except KeyError:
            return None
    
    def get_edge_from(self, parent):
        """Return the edge identifier for a parent->node relationship, or None
        if no relationship exists.

        Args:
            parent: A parent object
        """
        try:
            return self._parents[parent]
        except KeyError:
            return None


class DirectedGraph(object):
    """A collection of Node objects forming a directed graph.

    Attributes:
        _nodes (dict): node id => node instance
        _node_count (int): The number of nodes in the graph
    """
    def __init__(self):
        """Initialize the _nodes and _node_count attributes."""
        self._nodes = {}
        self._node_count = 0
    
    def add_node(self, node_id):
        """Add a node to the graph.

        Args:
            node_id: A unique identifier for the node
        Returns:
            node instance, or None if the node_id was already taken
        """
        try:
            node = self._nodes[node_id]
            return None
        except KeyError:
            node = Node(node_id)
            self._nodes[node_id] = node
            self._node_count += 1
            return node
    
    def add_edge(self, edge_id, parent_id, child_id):
        """Add an edge to the graph.

        Args:
            edge_id: A unique identifier for the edge
            parent_id: A unique identifier for the parent node
            child_id: A unique identifier for the child node
        """
        try:
            parent = self._nodes[parent_id]
        except KeyError:
            parent = self.add_node(parent_id)
        try:
            child = self._nodes[child_id]
        except KeyError:
            child = self.add_node(child_id)
        parent.add_child(child, edge_id)
        child.add_parent(parent, edge_id)
    
    def get_edge(self, parent_id, child_id):
        """Return the edge id identifying the relationship between a parent and
        child node, or None if no relationship exists.

        Args:
            parent_id: A unique identifier for a parent node in the graph
            child_id: A unique identifier for a child node in the graph
        """
        try:
            for parent in self._nodes[child_id].get_parents():
                if parent_id == parent.get_id():
                    return self._nodes[child_id].get_edge_from(parent)
            return None
        except KeyError:
            return None
    
    def get_edge_node_ids(self, edge_id):
        """Return the parent and child node ids for an edge, or None if the
        edge is not in the graph.

        Args:
            edge_id: An edge identifier for a node in the graph
        """
        for child_id in self.get_node_ids():
            child_node = self.get_node(child_id)
            for parent_node in child_node.get_parents():
                parent_id = parent_node.get_id()
                if self.get_edge(parent_id, child_id) == edge_id:
                    return (parent_id, child_id)
        print('Edge does not exist: %s.' % edge_id, file=sys.stderr)
        return None
    
    def get_node(self, node_id):
        """Return a node instance, or None if the node instance is not in the
        graph.

        Args:
            node_id: A unique identifier for a node
        """
        if node_id in self._nodes:
            return self._nodes[node_id]
        else:
            return None
    
    def get_node_ids(self):
        """Return a list of node ids corresponding to node instances in the
        graph.
        """
        return self._nodes.keys()


class MethodPath(object):
    """This class describes the relationships between compute methods, whose
    parent-child relationships form a directed graph which is described in
    mess.db and maintained by the methods in this class.

    Attributes:
        _path_id: An identifier for the path
        _path (list): Edges comprising the path, mirrors method_path_edge
        _db: A MessDB instance
        _graph: A DirectedGraph instance that mirrors the method_edge DB table
    """
    def __init__(self):
        """Initialize all attributes and execute _load_graph."""
        self._path_id = None
        self._path = []
        self._db = MessDB()
        self._graph = DirectedGraph()
        self._log = Log('all')
        self._load_graph()
    
    def _load_graph(self):
        """Load graph from mess.db method_edge table into _graph
        representation.
        """
        for row in self._db.execute('SELECT * FROM method_edge').fetchall():
            self._graph.add_edge(row.method_edge_id,
                                 row.parent_method_id,
                                 row.child_method_id)

    def _insert_edge(self, parent_id, child_id):
        """Insert edge between vertecies (methods) in mess.db. Returns edge id.
        """
        if parent_id is None:
            parent_id = child_id
        total_changes = self._db.total_changes
        insert_query = ('INSERT OR IGNORE INTO method_edge '
                        '(parent_method_id, child_method_id) '
                        'SELECT parent_method_id, ? FROM method_edge '
                        'WHERE child_method_id = ? UNION ALL SELECT ?, ?')
        self._db.execute(insert_query,
                         (child_id, parent_id, child_id, child_id))
        select_query = ('SELECT method_edge_id FROM method_edge '
                        'WHERE parent_method_id = ? AND child_method_id = ?')
        edge_id = self._db.execute(select_query,
                                   (parent_id, child_id)).fetchone()[0]
        self._load_graph()
        if self._db.total_changes - total_changes > 0:
            self._log.info('edge added to method calculation graph')
        return edge_id
    
    def _get_path_id_query(self, path):
        """Return a query that will return a path id.

        Args:
            path (dict): A path dictionary, like self._path
        """
        distance_query = ('SELECT method_path_id FROM method_path_edge '
                          'GROUP BY method_path_id HAVING sum(distance) = ?')
        edge_query = ('SELECT method_path_id FROM method_path_edge '
                      'WHERE method_edge_id = ?')
        return (' INTERSECT '.join([distance_query]
                                   + [edge_query]
                                   * len(path)),
                (sum(range(len(path))), ) + tuple(path))
    
    def _get_directory(self, parent_id, method_id, path_id):
        """Return a descriptive directory name, or None if none exists.

        Args:
            parent_id: Method id for parent method
            method_id: Method id for endpoint of path
            path_id: Path id
        """
        method_query = 'SELECT name, hash FROM method WHERE method_id = ?'
        method_result = self._db.execute(method_query,
                                         (method_id, )).fetchone()
        parent_method_result = self._db.execute(method_query,
                                                (parent_id, )).fetchone()
        if method_result == parent_method_result:
            return None
        try:
            return '%s-%s_FROM_%s-%s_PATH_%s' % (method_result.name,
                                                 method_result.hash[:7],
                                                 parent_method_result.name,
                                                 parent_method_result.hash[:7],
                                                 path_id)
        except AttributeError:
            return None
    
    def get_length(self):
        """Returns the length of the path."""
        return len(self._path) - 1
    
    def get_method_id(self):
        """Returns the proximal method id, or None if the path is empty."""
        try:
            child_id = self._graph.get_edge_node_ids(self._path[-1])[1]
            return child_id
        except (IndexError, TypeError):
            return None
    
    def get_parent_method_id(self):
        """Returns the method id of the penultimate method id, or None if the
        path is less than two methods long."""
        try:
            parent_id = self._graph.get_edge_node_ids(self._path[-1])[0]
            return parent_id
        except (IndexError, TypeError):
            return None
    
    def get_superparent_method_id(self):
        """Returns the method id of the grandparent method to the end method of
        the current path, or None if the path is less than three methods long.
        """
        try:
            superparent_id = self._graph.get_edge_node_ids(self._path[-2])[0]
            return superparent_id
        except (IndexError, TypeError):
            return None
    
    def get_path_id(self):
        """Return the current path_id, or None if the path is empty."""
        query, values = self._get_path_id_query(self._path)
        try:
            return self._db.execute(query, values).fetchone()[0]
        except TypeError:
            return None
    
    def get_parent_path_id(self):
        """Return the path id for the path that does not contain the last
        method of the current path."""
        query, values = self._get_path_id_query(self._path[:-1])
        try:
            return self._db.execute(query, values).fetchone()[0]
        except TypeError:
            return None
    
    def get_path_directory(self):
        """Return a directory name for the current path."""
        return self._get_directory(self.get_parent_method_id(),
                                   self.get_method_id(),
                                   self._path_id)
    
    def get_parent_path_directory(self):
        """Return a directory name for the parent path."""
        return self._get_directory(self.get_superparent_method_id(),
                                   self.get_parent_method_id(),
                                   self.get_parent_path_id())
    
    def set_path(self, path_id):
        """Load path from method_path/method_path_edge DB tables.

        Args:
            path_id: A unique path identifer
        """
        if path_id is None:
            # assume parent is import
            select_query = ('SELECT method_path_id FROM method_path '
                            'WHERE length = ? ORDER BY method_path_id ASC')
            try:
                path_id = self._db.execute(select_query, (0, )).fetchone()[0]
            except TypeError:
                pass
        else:
            select_query = ('SELECT method_path_id FROM method_path '
                            'WHERE method_path_id = ?')
            try:
                path_id = self._db.execute(select_query,
                                           (path_id, )).fetchone()[0]
            except TypeError:
                sys.exit('Path %s does not exist.' % path_id)
        self._path_id = path_id
        self._path = []
        select_query = ('SELECT me.child_method_id FROM method_edge AS me '
                        'JOIN method_path_edge AS mpe '
                        'ON me.method_edge_id = mpe.method_edge_id '
                        'WHERE mpe.method_path_id = ? '
                        'ORDER BY mpe.distance')
        for row in self._db.execute(select_query, (path_id, )).fetchall():
            self.extend_path(row.child_method_id)
        assert self._path_id == self.get_path_id()
    
    def setup_path(self, method_id, parent_path_id=None):
        """Load parent path and then extend it to include a new method.

        Args:
            method_id: A unique method identifer
            parent_path_id: An existing path id
        """
        self.set_path(parent_path_id)
        self.extend_path(method_id)
        assert self._path_id == self.get_path_id()
    
    def extend_path(self, method_id):
        """Create a new path that is identital to the current path but has an
        additional step.

        Args:
            method_id: A unique method identifier
        """
        edge_id = self._insert_edge(self.get_method_id(), method_id)
        if len(self._path) == 0 or edge_id != self._path[-1]:
            self._path.append(edge_id)
        elif (len(self._path) > 1
              and self.get_method_id() == method_id
              and edge_id == self._path[-1]):
            self._path.append(edge_id)
        parent_path_id = self._path_id
        self._path_id = self.get_path_id()
        if self._path_id is None:
            cur = self._db.cursor()
            cur.execute('INSERT INTO method_path (length) VALUES (?)',
                        (self.get_length(), ))
            self._path_id = cur.lastrowid
            if parent_path_id is None:
                parent_path_id = self._path_id
            cur.execute(('INSERT OR IGNORE INTO method_path_edge '
                         '(method_path_id, method_edge_id, distance) '
                         'SELECT ?, method_edge_id, distance '
                         'FROM method_path_edge WHERE method_path_id = ? '
                         'UNION ALL SELECT ?, method_edge_id, ? '
                         'FROM method_edge '
                         'WHERE parent_method_id = ? '
                         'AND child_method_id = ? '
                         'UNION ALL SELECT ?, method_edge_id, distance '
                         'FROM method_path_edge '
                         'WHERE method_path_id = ?'),
                        (self._path_id, self._path_id, self._path_id,
                         self.get_length(),
                         self.get_parent_method_id(),
                         self.get_method_id(),
                         self._path_id, parent_path_id))
            self._db.commit()
            self._log.info('method path extended by one in calculation graph')
