# -*- coding: utf-8 -*-
"""MESS.DB init module

This module initializes the __all__ array with the names of all of the tools
available in the 'tools' directory.
"""

from __future__ import unicode_literals

import os

__all__ = []
for filename in os.listdir(os.path.join(os.path.dirname(__file__), 'tools')):
    base, extension = os.path.splitext(filename)
    if extension == '.py' and not (filename.startswith('_')
                                   or filename.startswith('.')):
        __all__.append(base)
