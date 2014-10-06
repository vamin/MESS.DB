# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB tools module

This module contains the ToolManager class, which is used to load a tool for
execution, and the AbstractTool class, which all tool objects must inherit
from.
"""

from __future__ import print_function
from __future__ import unicode_literals

import os

from mess.log import Log
from mess.utils import CustomArgparseFormatter, is_inchikey


class ToolManager(object):
    """Handles the reading and loading of tool modules.
    
    Attributes:
    _tools (dict): Look up tool object by tool name
    _tool_names: Sorted list of tool names
    """
    def __init__(self):
        """Initialize the tools list."""
        self._tools = {}
        self._tool_names = []
        for filename in os.listdir(os.path.join(os.path.dirname(__file__),
                                                'tools')):
            base, extension = os.path.splitext(filename)
            if extension == '.py' and not (filename.startswith('_')
                                           or filename.startswith('.')):
                self._tool_names.append(base)
        self._tool_names.sort()  # so they show up alphabetized in help
    
    def populate_parser(self, parser):
        """Get argparse information from all of the tools."""
        subparsers = parser.add_subparsers(help='tools', dest='subparser_name')
        for tool_name in self._tool_names:
            tool = self.load_tool(tool_name)
            subparser = subparsers.add_parser(
                tool_name, help=tool.description,
                description=tool.description,
                epilog=tool.epilog,
                formatter_class=CustomArgparseFormatter)
            tool.subparse(subparser)
    
    def load_tool(self, tool_name):
        """Load a tool.
        
        Args:
            tool_name: The name of a valid tool.
        
        Returns:
            The tool object.
        """
        if tool_name not in self._tool_names:
            raise KeyError("tool '%s' not found" % tool_name)
        try:
            tool = self._tools[tool_name]
        except KeyError:
            # load the plugin only if not loaded yet
            module = __import__('mess.tools.%s' % tool_name,
                                fromlist=['tools'], level=0)
            tool = module.load()
            self._tools[tool_name] = tool
        return tool


# each tool must provide a load method at module level that will be
# used to instantiate the plugin (e.g. def load():\ return Tool())
class AbstractTool(object):
    """All tools must inherit from this class."""
    
    log_console = Log('console')
    log_all = Log('all')
    _inchikey = None
    
    def __init__(self):
        """Raise error if __init__ not implemented."""
        self.description = ''
        self.epilog = ''
        raise NotImplementedError(('every tool needs an __init__ that sets '
                                   'description and epilog attributes'))
    
    def subparse(self, subparser):
        """Set tool-specific arguments for argparser. Raise error if not
        implemented.
        """
        raise NotImplementedError("every tool needs a 'subparse' method")
    
    def execute(self, args):
        """Raise error if execute not implemented."""
        # all tools should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError("every tool needs an 'execute' method")
    
    @property
    def inchikey(self):
        """Get inchikey."""
        return self._inchikey
    
    @inchikey.setter
    def inchikey(self, inchikey):
        """Set inchikey, and update inchikey of logger."""
        if inchikey is not None and not is_inchikey(inchikey):
            raise RuntimeError('invalid inchikey: %s' % inchikey)
        self._inchikey = inchikey
        self.log_all.inchikey = inchikey
