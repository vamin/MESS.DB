# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB log module

This module contains classes relating to logging, including custom formatters
and filters for the 'logging' module. It also contains the Log class, which
wraps the 'logging' module to incorporate customizations.
"""

import functools
import inspect
import logging
import os
import sys
import textwrap
from datetime import date

from mess.utils import is_inchikey, get_inchikey_dir


class CustomLoggingFormatter(logging.Formatter):
    """This logging function overrides the default format method in
    logging.Formatter. The only change is the addition of a hanging indent
    using textwrap.
    """
    
    def format(self, record):
        """Format record with hanging indent."""
        record.message = record.getMessage()
        if self.usesTime():
            record.asctime = self.formatTime(record, self.datefmt)
        s = self._fmt % record.__dict__
        if record.exc_info:
            # Cache the traceback text to avoid converting it multiple times
            # (it's constant anyway)
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            try:
                s = s + record.exc_text
            except UnicodeError:
                # Sometimes filenames have non-ASCII chars, which can lead
                # to errors when s is Unicode and record.exc_text is str
                # See issue 8924.
                # We also use replace for when there are multiple
                # encodings, e.g. UTF-8 for the filesystem and latin-1
                # for a script. See issue 13232.
                s = s + record.exc_text.decode(sys.getfilesystemencoding(),
                                               'replace')
        hanging_indent = '                       '
        return textwrap.fill(s, width=79, initial_indent='',
                             subsequent_indent=hanging_indent)


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
    
    def filter(self, record):
        """Returns true if record should be processed."""
        record.modifiedname = record.name
        if self._handler == 'StreamHandler':
            return self._stream_filter(record)
        elif self._handler == 'FileHandler':
            return self._file_filter(record)
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


class Log(object):
    """Interface to logging module that incorporates customizations and makes
    it easy to broadcast to central and molecule logs."""
    
    inchikey = None
    
    def __init__(self, scope=None):
        """Sets the scope, which must be 'console' or 'all'."""
        if scope is None:
            raise RuntimeError(('Log object must be initialized with valid '
                                'scope: \'console\' or \'all\'.'))
        self.scope = scope
        self.context = self._context()
    
    @classmethod
    def setup(cls, mess_dir, verbose=False):
        """Setup logging environment, including stream and file handlers."""
        logging_formatter = CustomLoggingFormatter(('[%(asctime)s] '
                                                    '%(modifiedname)s '
                                                    '%(levelname)s: '
                                                    '%(message)s'),
                                                   '%Y-%m-%d %H:%M:%S')
        logging._defaultFormatter = logging_formatter
        logging.getLogger = GetLoggerDecorator(logging.getLogger)
        # setup main logger and stream handler
        logger = logging.getLogger('mess')
        stream_handler = logging.StreamHandler()
        if verbose:
            logger.setLevel(logging.DEBUG)
            stream_handler.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        logger.addHandler(stream_handler)
        # setup file handlers
        central_log_name = '%s.log' % date.today().strftime('%Y-%m-%d')
        central_log_path = os.path.join(mess_dir,
                                        '../logs/%s' % central_log_name)
        central_log_handler = logging.FileHandler(central_log_path)
        central_log_handler.setLevel(logging.INFO)
        logger.addHandler(central_log_handler)
        molecule_log_handler = logging.FileHandler('/dev/null')
        molecule_log_handler.setLevel(logging.INFO)
        logger.addHandler(molecule_log_handler)
    
    @classmethod
    def _context(cls):
        """Determine the name of the module being logged from."""
        frame = inspect.stack()[2][0]
        args, _, _, value_dict = inspect.getargvalues(frame)
        try:
            if args[0] == 'self':
                instance = value_dict.get('self', None)
                if instance:
                    return getattr(instance, '__class__', None).__name__
            else:
                return None
        except IndexError:
            return None
    
    def _to_all(self, inchikey=None):
        """Log to console, central log, and molecule log. Use to report things
        that change the database."""
        if inchikey is not None:
            if not is_inchikey(inchikey):
                sys.exit('invalid inchikey passed to logger')
            molecule_log_path = '%s/%s.log' % (get_inchikey_dir(inchikey),
                                               inchikey)
            molecule_log_handler = logging.FileHandler(molecule_log_path)
        logger = logging.getLogger('mess')
        
        for handler in logger.handlers:
            try:
                if ('molecules/' in handler.baseFilename
                        or '/dev/null' in handler.baseFilename):
                    if inchikey is not None:
                        logger.removeHandler(handler)
                        logger.addHandler(molecule_log_handler)
                        break
                    else:
                        logger.removeHandler(handler)
                        logger.addHandler('/dev/null')
                        break
            except AttributeError:
                continue
        if self.context is not None:
            return logging.getLogger('mess.%s' % self.context.lower())
        else:
            return logging.getLogger('mess')
    
    def _to_console(self):
        """Log to console only. Use to report things that don't result in a
        change in the database."""
        if self.context is not None:
            return logging.getLogger('mess.%s.consoleonly' %
                                     self.context.lower())
        else:
            return logging.getLogger('mess.consoleonly')
    
    def critical(self, *args):
        """Log a critical message."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().critical(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).critical(*args)
    
    def error(self, *args):
        """Log an error."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().error(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).error(*args)
    
    def warning(self, *args):
        """Log a warning."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().warning(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).warning(*args)
    
    def info(self, *args):
        """Log information."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().info(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).info(*args)
    
    def debug(self, *args):
        """Log a debug message."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().debug(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).debug(*args)
    
    def exception(self, *args):
        """Log a debug message."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().exception(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).exception(*args)
    
    def log(self, *args):
        """Log a debug message."""
        self.context = self._context()
        if self.scope == 'console':
            self._to_console().log(*args)
        elif self.scope == 'all':
            self._to_all(self.inchikey).log(*args)
