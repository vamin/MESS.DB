#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import imp
import os
import sqlite3
import sys

import pybel

from _db import MessDB
from _paths import Path
from _tools import AbstractTool

class Calculate(AbstractTool):
    def __init__(self):
        self.description = 'applys a method to a set of molecules'
    
    def subparse(self, subparser):
        subparser.add_argument('inchikeys', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='A list of inchikeys, file or passed in through STDIN')
        subparser.add_argument("-m", "--method", required=True, help='A method name.')
        subparser.add_argument("-pp", "--parent-path", default='', help="Specify a parent path id. If not set, start from InChI.")
        #subparser.add_argument("-s", "--state", help="Apply method to special state/condition (cation, anion, protonated, etc.)")
        #subparser.add_argument("-ss", "--parent-state", help="Apply method to special state/condition (cation, anion, protonated, etc.)")
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        # setup import method
        try:
            method = imp.load_source('method', os.path.join(os.path.dirname( __file__ ), '..', 'methods', args.method, 'method.py'))
        except IOError:
            sys.exit(args.method + " is not a valid method.")
        m = method.Method(db)
        m.setup()
        method_id = m.get_method_id()
        # setup path
        p = Path(db)
        parent_method_name = None
        if not args.parent_path:
            parent_method_name = 'import'
        p.setup(m.get_method_id(), args.parent_path, parent_method_name)
        # apply method to molecule
        molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
        # method dir
        method_dir = os.path.join(m.method_name + '_FROM_' + p.parent_method_name + '_PATH_' + str(p.path_id))
        # parent method dir
        if (p.parent_path_id):
            parent_method_dir = os.path.join(p.parent_method_name + '_FROM_' + p.superparent_method_name + '_PATH_' + str(p.parent_path_id))
        else:
            parent_method_dir = None
        # to be implemented: conditions (state) checking
        #if (args.state):
        #    method_dir = os.path.join('_', args.state, method_dir)
        #    parse_state()
        pybel.ob.obErrorLog.StopLogging() 
        method_args = {}
        method_args['method_dir'] = method_dir
        method_args['parent_method_dir'] = parent_method_dir
        method_args['path'] = p
        for row in args.inchikeys:
            method_args['inchikey'] = row.split()[0].strip()
            status = m.execute(method_args)
            print method_args['inchikey'] + ': ' + status
        

def load():
    # loads the current plugin
    return Calculate()
