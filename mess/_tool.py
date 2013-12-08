from __future__ import print_function
from __future__ import unicode_literals

from mess import __all__

class ToolManager(object):
    """Handles the reading and loading of tool modules."""
    def __init__(self):
        """Initialize the tools list."""
        self.__tools = {}
    
    def populate_parser(self, parser):
        """Get argparse information from all of the tools."""
        subparsers = parser.add_subparsers(help='tools', dest='subparser_name')
        __all__.sort() # so they show up alphabetized in help
        for tool_name in __all__:
            tool = self.load_tool(tool_name)
            subparser = subparsers.add_parser(tool_name, help=tool.description,
                                              description=tool.description,
                                              epilog=tool.epilog)
            tool.subparse(subparser)
    
    def load_tool(self, tool_name):
        """Load a tool.
        
        Args:
            tool_name: The name of a valid tool.
        
        Returns:
            The tool module.
        
        """
        if not tool_name in __all__:
            raise KeyError("tool '%s' not found" % tool_name)
        try:
            tool = self.__tools[tool_name]
        except KeyError:
            # load the plugin only if not loaded yet
            module = __import__('mess.tools.%s' % tool_name, 
                                fromlist=['tools'])
            tool = module.load()
            self.__tools[tool_name] = tool
        return tool
    
    def list(self):
        """Returns a list of the available tools."""
        return __all__


class AbstractTool(object):
    """All tools should inherit from this class."""
    def __init__(self):
        """Raise error if __init__ not implemented."""
        self.description = ''
        self.epilog = ''
        raise NotImplementedError(('every tool needs an __init__ with a '
                                   'description and epilog attribute'))
    
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


# each tool must provide a load method at module level that will be
# used to instantiate the plugin (e.g. def load():\ return Tool())