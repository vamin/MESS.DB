# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB inspect module

This module contains the inspect tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse

from mess.tool import AbstractTool


class Inspect(AbstractTool):
    """Get information about an element of MESS.DB."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Inspect the database, a molecule, or a calculation'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        pass
    
    def execute(self, args):
        """Inspect stuff."""
        sys.exit('not yet implemented')


def load():
    """Load Inspect()."""
    return Inspect()
