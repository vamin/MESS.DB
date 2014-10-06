# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB transform module

This module contains the transform tool class and load function.
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse

from mess.tool import AbstractTool


class Transform(AbstractTool):
    """Use reaction SMARTS to transform molecules."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Transform molecules by reaction SMARTS'
        self.epilog = ''
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        pass
    
    def execute(self, args):
        """Transform stuff."""
        sys.exit('not yet implemented')


def load():
    """Load Transform()."""
    return Transform()
