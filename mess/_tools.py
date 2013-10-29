from mess import __all__

class ToolsManager(object):
    def __init__(self):
        # initialize the tools list
        self.__tools = {}
    
    def populate_parser(self, parser):
        subparsers = parser.add_subparsers(help='tools', dest='subparser_name')
        for tool_name in __all__:
            tool = self.load_tool(tool_name)
            tool_subparser = subparsers.add_parser(tool_name, 
                                                   help=tool.description)
            tool.subparse(tool_subparser)
    
    def load_tool(self, tool_name):
        # loads a single tool given its name
        if not tool_name in __all__:
            raise KeyError("tool '" + tool_name + "' not found")
        try:
            tool = self.__tools[tool_name]
        except KeyError:
            # load the plugin only if not loaded yet
            module = __import__('mess.' + tool_name, fromlist=['tools'])
            tool = module.load()
            self.__tools[tool_name] = tool
        return tool
    
    def list(self):
        # returns a list of the available tools
        return __all__


class AbstractTool(object):
    # all tools should inherit from this class
    
    def subparse(self, subparser):
        # this method sets tool-specific arguments for argparser
        raise NotImplementedError("every tool needs a 'subparse' method")

    def check_dependencies(self):
        raise NotImplementedError(("every tool needs a 'check_dependencies'"
                                   'method'))
    
    def execute(self, args):
        # all tools should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError("every tool needs an 'execute' method")

# each tool must provide a load method at module level that will be
# used to instantiate the plugin (e.g. def load():\ return Tool())