"""
Unit tests for Node, DirectedGraph, and MethodPath mess/_path.py.
"""

import os
import shutil
import sys
import unittest

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../..'))
from mess.path import Node, DirectedGraph, MethodPath
from mess.db import MessDB
from mess.tests.helpers import suppress_stderr


class TestNode(unittest.TestCase):
    def setUp(self):
        self.node = Node(42)
    
    def test_node_init(self):
        self.assertEquals(self.node.get_id(), 42)
    
    def test_node_parents(self):
        parent_node = Node('parent')
        self.node.add_parent(parent_node, 'edge')
        count = 0
        for parent in self.node.get_parents():
            self.assertEquals(parent, parent_node)
            count += 1
        self.assertEquals(count, 1)
        another_node = Node('another')
        self.node.add_parent(parent_node, 'edge')
        self.node.add_parent(another_node, 'edge')
        for parent in self.node.get_parents():
            self.assertIn(parent, (parent_node, another_node))
            self.assertIn(parent.get_id(), ('parent', 'another'))
            count += 1
        self.assertEquals(count, 3)
    
    def test_node_children(self):
        child_node = Node('child')
        self.node.add_child(child_node, 'edge')
        count = 0
        for child in self.node.get_children():
            self.assertEquals(child, child_node)
            count += 1
        self.assertEquals(count, 1)
        another_node = Node('another')
        self.node.add_child(child_node, 'edge')
        self.node.add_child(another_node, 'edge')
        for child in self.node.get_children():
            self.assertIn(child, (child_node, another_node))
            self.assertIn(child.get_id(), ('child', 'another'))
            count += 1
        self.assertEquals(count, 3)
    
    def test_node_edges(self):
        parent_node = Node('parent')
        child_node = Node('child')
        self.node.add_parent(parent_node, 'parent_edge')
        self.node.add_child(child_node, 'child_edge')
        self.assertEquals(self.node.get_edge_to(child_node), 'child_edge')
        self.assertEquals(self.node.get_edge_from(parent_node), 'parent_edge')


class TestDirectedGraph(unittest.TestCase):
    def setUp(self):
        self.dg = DirectedGraph()
    
    def test_add_node(self):
        self.assertIsNone(self.dg.get_node(42))
        self.dg.add_node(42)
        self.assertEquals(self.dg.get_node(42).get_id(), 42)
        self.dg.add_node(24)
        count = 0
        for node_id in self.dg.get_node_ids():
            self.assertIn(node_id, (42, 24))
            count += 1
        self.assertEquals(count, 2)
        self.assertEquals(self.dg._node_count, 2)
    
    def test_add_edge(self):
        self.assertIsNone(self.dg.get_edge(42, 24))
        self.dg.add_edge('edge', 42, 24)
        self.assertEquals(self.dg.get_edge(42, 24), 'edge')
        parent_id, child_id = self.dg.get_edge_node_ids('edge')
        self.assertEquals(parent_id, 42)
        self.assertEquals(child_id, 24)


class TestMethodPath(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = './tmp'
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        self.path = MethodPath()
        with suppress_stderr():  # silence 'MESS.DB created' message
            self.path._db = MessDB(database='%s/test.db' % self.tmp_dir)
        self.path._graph = DirectedGraph()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def assert_path_consistency(self):
        self.assertEquals(self.path._path_id, self.path.get_path_id())

    def test_init(self):
        self.assertEquals(self.path._db.tries, 0)
        self.assertEquals(self.path._graph._node_count, 0)
        self.assertEquals(self.path._path, [])
        self.assertIsNone(self.path._path_id)
        self.assert_path_consistency()
    
    def test_load_graph(self):
        insert_query = 'INSERT INTO method_edge VALUES (?, ?, ?)'
        self.path._db.executemany(insert_query, ((1, 1, 1),
                                                 (2, 1, 2),
                                                 (3, 2, 3)))
        self.path._load_graph()
        self.assertEquals(sorted(self.path._graph.get_node_ids()), [1, 2, 3])
        self.assert_path_consistency()

    def test_setup_path(self):
        self.path.setup_path(1)
        self.assertEquals(self.path._graph._node_count, 1)
        self.assertEquals(self.path._path, [1])
        self.assertEquals(self.path.get_length(), 0)
        self.assert_path_consistency()
        # check that new path isn't added for same method
        new_path = MethodPath()
        new_path._db = self.path._db
        new_path._graph = DirectedGraph()
        new_path._load_graph()
        new_path.setup_path(1)
        self.assertEquals(new_path._graph._node_count, 1)
        self.assertEquals(new_path._path, [1])
        self.assertEquals(new_path.get_length(), 0)

        
    def test_extend_path(self):
        self.path.setup_path(42)
        self.path.extend_path(43)
        self.assertEquals(self.path._graph._node_count, 2)
        self.assertEquals(self.path._path, [1, 2])
        self.assertEquals(self.path.get_length(), 1)
        self.path.extend_path(44)
        self.assertEquals(self.path._graph._node_count, 3)
        self.assertEquals(self.path._path, [1, 2, 5])
        self.assertEquals(self.path.get_length(), 2)
        self.path.extend_path(42)
        self.assertEquals(self.path._graph._node_count, 3)
        self.assertEquals(self.path._path, [1, 2, 5, 8])
        self.assertEquals(self.path.get_length(), 3)
        self.assert_path_consistency()
        
    def test_set_path(self):
        self.path.setup_path(42)
        self.path.extend_path(43)
        self.path.extend_path(44)
        self.path.set_path(2)
        self.assertEquals(self.path._graph._node_count, 3)
        self.assertEquals(self.path._path, [1, 2])
        self.assertEquals(self.path.get_length(), 1)
        self.assert_path_consistency()
    
    def test_get_path_id(self):
        self.assertIsNone(self.path.get_path_id())
        self.path.setup_path(42)
        self.path.extend_path(43)
        self.path.extend_path(44)
        self.assertEquals(self.path.get_path_id(), 3)
    
    def test_get_parent_path_id(self):
        self.assertIsNone(self.path.get_parent_path_id())
        self.path.setup_path(42)
        self.path.extend_path(43)
        self.path.extend_path(44)
        self.assertEquals(self.path.get_parent_path_id(), 2)
    
    def test_get_method_id(self):
        self.assertIsNone(self.path.get_method_id())
        self.path.setup_path(42)
        self.assertEquals(self.path.get_method_id(), 42)
        self.path.extend_path(43)
        self.assertEquals(self.path.get_method_id(), 43)
    
    def test_get_parent_method_id(self):
        self.assertIsNone(self.path.get_parent_method_id())
        self.path.setup_path(42)
        self.assertEquals(self.path.get_parent_method_id(), 42)
        self.path.extend_path(43)
        self.assertEquals(self.path.get_parent_method_id(), 42)
        self.path.extend_path(44)
        self.assertEquals(self.path.get_parent_method_id(), 43)
    
    def test_get_superparent_method_id(self):
        self.assertIsNone(self.path.get_superparent_method_id())
        self.path.setup_path(42)
        self.assertEquals(self.path.get_parent_method_id(), 42)
        self.path.extend_path(43)
        self.assertEquals(self.path.get_parent_method_id(), 42)
        self.path.extend_path(44)
        self.assertEquals(self.path.get_superparent_method_id(), 42)
        self.path.extend_path(45)
        self.assertEquals(self.path.get_superparent_method_id(), 43)
    
    def test_get_directory(self):
        self.assertIsNone(self.path._get_directory(1, 2, 3))

    def test_get_path_directory(self):
        self.assertIsNone(self.path.get_path_directory())

    def test_get_parent_path_directory(self):
        self.assertIsNone(self.path.get_parent_path_directory())


if __name__ == '__main__':
    unittest.main()
