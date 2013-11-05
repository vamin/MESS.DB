from __future__ import print_function
from __future__ import unicode_literals

import collections
import functools
import inspect

def decorate(object_, decorator):
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
                                decorator(getattr(object_, attr)))
                    except TypeError:
                        pass
            except AttributeError:
                pass


class UnicodeDecorator(object):
    """Make sure that unicode inputs are encoded and that all return values are
    converted back to unicode.

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
        """Decorated (i.e. converted to unicode) function output.

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
            r = type(r)(map(self.unicode_replace, r))
        elif isinstance(r, collections.Mapping):
            r = dict(map(self.unicode_replace, r.iteritems()))
        elif isinstance(r, str):
            r = self.unicode_replace(r)
        return r

    def unicode_replace(self, x, enc='utf-8', err='replace'):
        """Convert str to unicode.

        Args:
            x: A string.
            enc: Encoding of input string, defaults to 'utf-8'.
            err: What to do on unicode conversion error, defaults to 'replace'.

        Returns:
            Unicode string if x is str, x otherwise.

        """
        if isinstance(x, str):
            return unicode(x, enc, err)
        else:
            return x