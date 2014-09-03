# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB decorators module

This module contains various decorator classes and functions that are used
throughout MESS.DB.
"""

from __future__ import print_function
from __future__ import unicode_literals

import collections
import functools
import inspect

from mess.utils import unicode_replace


def decorate(object_, decorator, *args, **kwargs):
    """Apply a decorator to all callable functions of an object.
    
    Args:
        object_: An object containing callable functions to be decorated.
        decorator: A decorator.
    """
    
    for attr in dir(object_):
        if '__' not in attr:
            try:
                if inspect.isbuiltin(getattr(object_, attr)):
                    continue
                if callable(getattr(object_, attr)):
                    try:
                        setattr(object_, attr,
                                decorator(getattr(object_, attr)),
                                *args, **kwargs)
                    except TypeError:
                        pass
            except AttributeError:
                pass
    return object_


class HandlerSelectiveFilter(object):
    """Filter messages depending on handler type.
    
    Attributes:
        _handler (str): The name of the handler being filtered.
    """
    
    def __init__(self, handler):
        """Initialize attributes.
        
        Args:
            handler (str): The name of the handler being filtered.
        """
        self._handler = handler
    
    def filter(self, record):
        """Returns true if record should be processed."""
        record.modifiedname = record.name
        if self._handler == 'StreamHandler':
            return self._stream_filter(record)
        elif self._handler == 'FileHandler':
            return self._file_filter(record)
        else:
            return True
    
    def _stream_filter(self, record):
        """Returns true if record name indicates console-specific."""
        if 'console' in record.name.split('.')[-1]:
            record.modifiedname = '.'.join(record.name.split('.')[:-1])
        return True
    
    def _file_filter(self, record):
        """Returns false if record name indicates console-specific."""
        if 'console' in record.name.split('.')[-1]:
            return False
        else:
            return True


class AddHandlerDecorator(object):
    """Decorates a logger addHandler method with handler selective filtering.
    """
    
    def __init__(self, func):
        """Inits AddHandlerDecorator.
        
        Args:
            func: A logger addHandler method.
        """
        self.func = func
        functools.update_wrapper(self, func)
    
    def __call__(self, handler):
        """Returns addHandler after adding HandlerSelectiveFilter."""
        handler.addFilter(HandlerSelectiveFilter(handler.__class__.__name__))
        return self.func(handler)


class GetLoggerDecorator(object):
    """Decorates all addHandler method in logging module with handler selective
    filtering.
    """
    
    def __init__(self, func):
        """Inits GetLoggerDecorator.
        
        Args:
            func: A logging module getLogger method.
        """
        self.func = func
        functools.update_wrapper(self, func)
    
    def __call__(self, name):
        """Decorates addHandler with AddHandlerDecorator and returns logger."""
        logger = self.func(name)
        logger.addHandler = AddHandlerDecorator(logger.addHandler)
        return logger


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
        for arg in args:
            if isinstance(arg, unicode):
                encoded_args.append(arg.encode())
            else:
                encoded_args.append(arg)
        for k, v in kwargs.items():
            if isinstance(k, unicode):
                encoded_k = k.encode()
            else:
                encoded_k = k
            if isinstance(v, unicode):
                encoded_kwargs[encoded_k] = v.encode()
            else:
                encoded_kwargs[encoded_k] = v
        result = self.func(*encoded_args, **encoded_kwargs)
        if isinstance(result, (list, set)):
            result = type(result)(map(unicode_replace, result))
        elif isinstance(result, collections.Mapping):
            result = dict(map(unicode_replace, result.iteritems()))
        elif isinstance(result, str):
            result = unicode_replace(result)
        return result
