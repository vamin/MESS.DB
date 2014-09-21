# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB remove module

This module contains the remove tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import shutil
import os
import sys

from mess._db import MessDB
from mess._tool import AbstractTool
from mess.utils import get_inchikey_dir


class Remove(AbstractTool):
    """This tool removes molecules from the molecules directory and the
    database."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = ('Remove molecules from MESS.DB')
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('inchikeys', nargs='?',
                               type=argparse.FileType('r'), default=sys.stdin,
                               help='a list of inchikeys (default: STDIN)')
    
    def execute(self, args):
        """Remove specified elements."""
        db = MessDB()
        cur = db.cursor()
        for row in args.inchikeys:
            inchikey = row.split()[0].strip()
            try:
                inchikey_dir = get_inchikey_dir(inchikey)
                shutil.rmtree(inchikey_dir)
                self.log_all.info('%s dir removed', inchikey)
            except OSError:
                self.log_console.info('%s did not have a directory', inchikey)
            try:
                parent = os.path.relpath(os.path.join(inchikey_dir, '../'))
                os.removedirs(parent)
            except OSError:
                pass
            records = 0
            query = 'DELETE from molecule WHERE inchikey=?'
            cur.execute(query, (inchikey,))
            records += cur.rowcount
            query = 'DELETE from molecule_synonym WHERE inchikey=?'
            cur.execute(query, (inchikey,))
            records += cur.rowcount
            query = 'DELETE from molecule_source WHERE inchikey=?'
            cur.execute(query, (inchikey,))
            records += cur.rowcount
            query = ('DELETE from molecule_state_method_property '
                     'WHERE inchikey=?')
            cur.execute(query, (inchikey,))
            records += cur.rowcount
            db.commit()
            self.log_all.info('%i %s records removed from db',
                              records, inchikey)


def load():
    """Load Remove()."""
    return Remove()
