from __future__ import print_function
from __future__ import unicode_literals

import collections
import functools
import inspect
import types

from _utils import unicode_replace

def decorate(object_, decorator, *args, **kwargs):
    """Apply a decorator to all callable functions of an object.
    
    Args:
        object_: An object containing callable functions to be decorated.
        decorator: A decorator.
    
    """
    for attr in dir(object_):
        if not '__' in attr:
            try:
                if inspect.isbuiltin(getattr(object_, attr)):
                    continue
                if callable(getattr(object_, attr)):
                    try:
                        setattr(object_, attr,
                                decorator(getattr(object_, attr)), *args,
                                          **kwargs)
                    except TypeError:
                        pass
            except AttributeError:
                pass
    return object_

class ExceptionDecorator(object):
    def __init__(self, func, exceptions, replacement=None):
        """Inits ExceptionDecorator with function to be decorated and updates
        wrapper.
        
        Args:
            func: A function to be decorated.
            exception: The exception(s) to be caught.
            replacement: A function that accepts the decorated function,
                         or a return value.
        
        """
        self.func = func
        self.exceptions = exceptions
        self.replacement = replacement
        functools.update_wrapper(self, func)
    
    def __call__(self, *args, **kwargs):
        if self.func is None:
            self.func = args[0]
            return self
        try:
            return self.func(*args, **kwargs)
        except self.exceptions:
            if (isinstance(self.replacement, types.FunctionType)):
                return replacement(self.func)
            elif (self.replacement):
                return self.replacement
            else:
                pass


class UnicodeDecorator(object):
    """Make sure that unicode inputs are encoded and that all return values
    are converted back to unicode.
    
    """
    def __init__(self, func):
        """Inits UnicodeDecorator with function to be decorated and updates
        wrapper.
        
        Args:
            func: A function to be decorated.
        
        """
        self.func = func
        functools.update_wrapper(self, func)
    
    def __call__(self, *args, **kwargs):
        """Call decorated (i.e. converted to unicode) function.
        
        Args:
            *args: Arguments of func which will be converted to bytes.
            **kwargs: Keyword arguments of func, which will be converted to
                      bytes.
        
        Returns:
            The output of func, converted to unicode if possible.
        
        """
        encoded_args = []
        encoded_kwargs = {}
        for a in args:
            if isinstance(a, unicode):
                encoded_args.append(a.encode())
            else:
                encoded_args.append(a)
        for k, v in kwargs.items():
            if isinstance(k, unicode):
                encoded_k = k.encode()
            else:
                encoded_k = k
            if isinstance(v, unicode):
                encoded_kwargs[encoded_k] = v.encode()
            else:
                encoded_kwargs[encoded_k] = v
        r = self.func(*encoded_args, **encoded_kwargs)
        if isinstance(r, (list, set)):
            r = type(r)(map(unicode_replace, r))
        elif isinstance(r, collections.Mapping):
            r = dict(map(unicode_replace, r.iteritems()))
        elif isinstance(r, str):
            r = unicode_replace(r)
        return r        