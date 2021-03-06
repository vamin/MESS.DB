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
