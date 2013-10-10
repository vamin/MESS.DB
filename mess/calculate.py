#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import glob
import imp
import os
import sqlite3
import sys
import time
from cStringIO import StringIO

import pybel

from _db import MessDB
from _tools import AbstractTool

class Calculate(AbstractTool):
    def __init__(self):
        self.description = 'applys a method to a set of molecules'
    
    def subparse(self, subparser):
        subparser.add_argument('inchikeys', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='A list of inchikeys, file or passed in through STDIN')
        subparser.add_argument("-m", "--method", required=True, help='A method name or method_id.')
        subparser.add_argument("-pp", "--parent-path", help="Specify a parent path id. If not set, start from InChI.")
        #subparser.add_argument("-s", "--state", help="Apply method to special state/condition (cation, anion, protonated, etc.)")
        #subparser.add_argument("-ss", "--parent-state", help="Apply method to special state/condition (cation, anion, protonated, etc.)")
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        # import method
        (method_id, method_name, level_of_theory) = self.parse_method(c, args.method)
        method = imp.load_source('method', os.path.join(os.path.dirname( __file__ ), '..', 'methods', level_of_theory, method_name, 'method.py'))
        m = method.Method(db)
        # parse path
        (path_id, parent_method_name, parent_level_of_theory, parent_path_id, superparent_method_name) = self.parse_path(c, method_id, args.parent_path)
        db.commit()
        # apply method to molecule
        molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
        # method dir
        method_dir = os.path.join(level_of_theory, method_name + '_FROM_' + parent_method_name + '_PATH_' + str(path_id))
        # parent method dir
        if (parent_path_id):
            parent_method_dir = os.path.join(parent_level_of_theory, parent_method_name + '_FROM_' + superparent_method_name + '_PATH_' + parent_path_id)
        else:
            parent_method_dir = None
        # to be implemented: conditions (state) checking
        #if (args.state):
        #    method_dir = os.path.join('_', args.state, method_dir)
            #parse_state(
        method_args = {}
        method_args['method_name'] = method_name
        method_args['method_dir'] = method_dir
        method_args['parent_method_dir'] = parent_method_dir
        method_args['path_id'] = path_id
        for row in args.inchikeys:
            method_args['inchikey'] = row.split()[0].strip()
            status = m.execute(method_args)
            print method_args['inchikey'] + ': ' + status

    def parse_path(self, c, method_id, path_id):
        if (path_id):
            # parent path specified
            path_check_row = c.execute('SELECT mp.method_path_id, mp.length, me.parent_method_id, me.child_method_id, m.name method_name, l.name level_name FROM method_path mp JOIN method_path_edge mpe ON mpe.distance = mp.length AND mpe.method_path_id = mp.method_path_id JOIN method_edge me ON me.method_edge_id = mpe.method_edge_id JOIN method m ON me.child_method_id = m.method_id JOIN level l ON l.level_id = m.level_id WHERE mp.method_path_id=? AND distance=1;', (path_id,)).fetchone()
            # check if valid
            if (path_check_row is None):
                sys.exit(path_id + " is an invalid path id (i.e., it does not have a valid record in the database).")
            path_length = path_check_row['length'] + 1
            parent_path_id = path_check_row['method_path_id']
            parent_method_id = path_check_row['child_method_id']
            parent_method_name = path_check_row['method_name']
            superparent_method_id = path_check_row['parent_method_id']
            # insert edges
            c.execute('INSERT OR IGNORE INTO method_edge (parent_method_id, child_method_id) SELECT parent_method_id, ? FROM method_edge WHERE child_method_id = ? UNION ALL SELECT ?, ?;', (method_id, parent_method_id, method_id, method_id))
            # check for path
            method_path_parent_row = c.execute('SELECT method_path_id FROM method_path_parent WHERE method_id=? AND parent_method_path_id=?', (method_id, parent_path_id)).fetchone()
            try:
                path_id = method_path_parent_row['method_path_id']
            except TypeError:
                # insert new path
                c.execute('INSERT INTO method_path (length) VALUES (?);', (path_length,))
                path_id = c.lastrowid
                c.execute('INSERT INTO method_path_parent (method_path_id, parent_method_path_id, method_id) VALUES (?, ?, ?);', (path_id, parent_path_id, method_id)) 
            # populate path edges
            c.execute('INSERT OR IGNORE INTO method_path_edge (method_path_id, method_edge_id, distance) SELECT ?, method_edge_id, distance FROM method_path_edge WHERE method_path_id = ? UNION ALL SELECT ?, method_edge_id, ? FROM method_edge WHERE parent_method_id = ? AND child_method_id = ?', (parent_path_id, path_id, parent_path_id, path_length, parent_method_id, method_id))
            # get parent method name and level of theory
            c.execute('SELECT m.name method_name, l.name level_name FROM method m JOIN level l ON l.level_id = m.level_id WHERE m.method_id = ?', (parent_method_id,))
            parent_method_row = c.fetchone()
            c.execute('SELECT name FROM method WHERE method_id = ?', (superparent_method_id,))
            superparent_method_row = c.fetchone()
            # return path row
            return (path_id, parent_method_row['method_name'], parent_method_row['level_name'], str(parent_path_id), superparent_method_row['name'])
            #(path_id, parent_method_name, parent_level_of_theory, parent_path_id, superparent_method_name)
        else:
            # parent is inchi
            parent_method_row = c.execute("SELECT method_id FROM method WHERE name='inchi';").fetchone()
            try:
                parent_method_row['method_id']
            except TypeError:
                sys.exit("INCHI is an invalid path id (i.e., it does not have a valid record in the database).")
            path_length = 1
            parent_method_id = parent_method_row['method_id']
            # insert edges
            c.execute('INSERT OR IGNORE INTO method_edge (parent_method_id, child_method_id) SELECT parent_method_id, ? FROM method_edge WHERE child_method_id = ? UNION ALL SELECT ?, ?;', (method_id, parent_method_id, method_id, method_id))
            # check for path
            method_path_parent_row = c.execute('SELECT method_path_id FROM method_path_parent WHERE method_id=? AND parent_method_path_id="";', (method_id,)).fetchone()
            try:
                path_id = method_path_parent_row['method_path_id']
            except TypeError:
                # insert new path
                c.execute('INSERT INTO method_path (length) VALUES (?);', (path_length,))
                path_id = c.lastrowid
                c.execute('INSERT INTO method_path_parent (method_path_id, parent_method_path_id, method_id) VALUES (?, ?, ?);', (path_id, '', method_id))     
            # populate path edges
            c.execute('INSERT OR IGNORE INTO method_path_edge (method_path_id, method_edge_id, distance) SELECT ?, method_edge_id, 1 FROM method_edge WHERE parent_method_id = ? AND child_method_id = ?;', (path_id, parent_method_id, method_id))
            # get parent method name and level of theory
            c.execute('SELECT name FROM method WHERE method_id = ?;', (parent_method_id,))
            parent_method_row = c.fetchone()
            # return path row
            return (path_id, parent_method_row['name'], None, None, None)
            #(path_id, parent_method_name, parent_level_of_theory, parent_path_id, superparent_method_name)

    def parse_method(self, c, method):
        try:
            [method_sql] = glob.glob(os.path.join(os.path.dirname(__file__), '../db/methods/*/', method + '.sql'))
        except ValueError:
            sys.exit(method +' is not a valid method with a valid sql file.')
        c.executescript(open(method_sql).read())
        c.execute("".join(['SELECT level.name l_name, method.name m_name, method.method_id FROM method JOIN level ON level.level_id = method.level_id WHERE method.name = \'', method, '\' OR method.method_id = \'', method, '\';']))
        method_row = c.fetchone()
        method_path = os.path.join(os.path.dirname(__file__), '../methods', method_row['l_name'])
        sys.path.append(method_path)
        return (method_row['method_id'], method_row['m_name'], method_row['l_name'])


def load():
    # loads the current plugin
    return Calculate()