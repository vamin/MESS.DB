# -*- coding: utf-8 -*-
"""MESS.DB method module

This module contains the AbstractMethod class, which every MESS.DB calculation
method object must inherit from.
"""

from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import json
import sys
import time

from mess.db import MessDB
from mess.log import Log
from mess.path import MethodPath
from mess.utils import is_inchikey


class AbstractMethod(object):
    """All methods must inherit from this class.
    
    Attributes:
        db (obj): A MessDB object
        method_name (str): The name of the method
        description (str): Description of method
        geop (bool): Whether the method generates a new geometry
        prog_name (str): Program name
        prog_version (str): Program version
        prog_url (str): Program url
        parameters (dict): Parameters that affect program execution
    """
    parameters = dict()
    shortdesc = None
    method_citation = None
    prog_citation = None
    _inchikey = None
    _path_id = None
    _parent_path_id = None
    _method_dir = None
    _parent_method_dir = None
    
    def __init__(self):
        """Set up db, check for attributes, dependencies, and setup."""
        self.db = MessDB()
        self.path = MethodPath()
        self.log_console = Log('console')
        self.log_all = Log('all')
        self.method_name = self.get_method_name()
        try:
            self.parameters
            self.description
            self.geop  # flag indicates method results in new xyz coordinates
            self.prog_name
            self.prog_version
            self.prog_url
        except AttributeError as err:
            print(''.join([str(err), '\n']), file=sys.stderr)
            sys.exit(('Each method class needs to define description, geop, '
                      'prog_name, prog_version, prog_url, '
                      'parameters as attributes.'))
        self.check_dependencies()
    
    def __hash__(self):
        """Hash based on method name and parameters.
        
        Returns:
            A hex string of the sha1 hash of self.method_name plus
            JSON-serialized self.parameters. Keys are sorted.
        """
        return hashlib.sha1(self.method_name +
                            json.dumps(dict((str(k).lower(),
                                             str(v).lower())
                                            for k, v
                                            in self.parameters.iteritems()),
                                       sort_keys=True)).hexdigest()
    
    @property
    def hash(self):
        """Get hash."""
        return self.__hash__()
    
    @property
    def method_id(self):
        """Get the object's method_id attribute."""
        query = ('SELECT method_id FROM method '
                 'WHERE hash = ?;')
        row = self.db.execute(query, (self.hash,)).fetchone()
        return row.method_id
    
    @property
    def path_id(self):
        """Get the path id of the method."""
        if not self.path.get_method_id() == self.method_id:
            self._setup_path()
        return self._path_id
    
    @property
    def method_dir(self):
        """Get the directory name of the method."""
        if not self.path.get_method_id() == self.method_id:
            self._setup_path()
        return self._method_dir
    
    @property
    def parent_method_dir(self):
        """Get the parent directory name of the method."""
        if not self.path.get_method_id() == self.method_id:
            self._setup_path()
        return self._parent_method_dir
    
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
    
    @classmethod
    def get_method_name(cls):
        """Return the name of the method, derived from the subclass name."""
        return cls.__name__.replace('_', '').lower()
    
    def _setup_path(self):
        """Setup path given current method id and parent path."""
        self.path.setup_path(self.method_id, self._parent_path_id)
        self._path_id = self.path.get_path_id()
        self._method_dir = self.path.get_path_directory()
        self._parent_method_dir = self.path.get_parent_path_directory()
    
    def _insert_method(self):
        """Set insert program to db, set up hash, and insert method to db."""
        total_changes = self.db.total_changes
        query = ('INSERT OR IGNORE INTO method '
                 '(program_id, geop, name, shortdesc, citation, hash) '
                 'SELECT program.program_id, ?, ?, ?, ?, ? '
                 'FROM program '
                 'WHERE program.name=? AND program.version=?')
        self.db.execute(query, (self.geop, self.method_name, self.shortdesc,
                                self.method_citation, self.hash,
                                self.prog_name, self.prog_version))
        if self.db.total_changes - total_changes > 0:
            self.log_all.info('new %s method added to MESS.DB',
                              self.method_name)
    
    def _insert_program(self):
        """Adds row to program table in mess.db."""
        total_changes = self.db.total_changes
        query = ('INSERT OR IGNORE INTO program '
                 '(name, version, url, citation) '
                 'VALUES (?, ?, ?, ?)')
        self.db.execute(query,
                        (self.prog_name, self.prog_version, self.prog_url,
                         self.prog_citation))
        if self.db.total_changes - total_changes > 0:
            self.log_all.info('program %s %s added to MESS.DB',
                              self.prog_name, self.prog_version)
    
    def _insert_parameters(self):
        """Import paramaters dict to mess.db.
        
        Args:
            name: Name of parameter.
            setting: The value the parameter is set to.
        """
        added_parameters = 0
        for name, setting in self.parameters.items():
            query = ('INSERT OR IGNORE INTO parameter (name) VALUES (?)')
            self.db.execute(query, (name, ))
            total_changes = self.db.total_changes
            query = ('INSERT OR IGNORE INTO method_parameter '
                     '(method_id, parameter_id, setting) '
                     'SELECT ?, parameter.parameter_id, ? '
                     'FROM program, parameter '
                     'WHERE parameter.name=?')
            self.db.execute(query, (self.method_id, setting, name))
            added_parameters += (self.db.total_changes - total_changes)
        if added_parameters > 0:
            self.log_all.info('%i method parameters added to MESS.DB',
                              added_parameters)
    
    def get_insert_property_query(self, inchikey, name, description,
                                  format_, value, units=''):
        """Returns query to insert property value to mess.db.
        
        Args:
            inchikey: The inchikey of a molecule in MESS.DB.
            method_path_id: Path id for the calculations that generated the
                            property.
            name: The property name.
            description: A description of the property.
            format_: A description of the format the property is in.
            value: The calculated property.
            units: Units for the property value.
        """
        query = ('INSERT OR IGNORE INTO molecule_method_property_denorm '
                 'VALUES (?, ?, ?, ?, ?, ?, ?);')
        return (query, (inchikey, self.path_id, name, description,
                        format_, units, value))
    
    def get_insert_moldata_queries(self, inchikey, mol,
                                   description='', units=''):
        """Returns queries to insert molecule data values to mess.db."""
        for name, value in mol.data.iteritems():
            yield self.get_insert_property_query(inchikey,
                                                 name,
                                                 description,
                                                 type(value).__name__,
                                                 value,
                                                 units)
    
    def get_timing_query(self, inchikey, start):
        """Get a query to insert execution time property into db."""
        return self.get_insert_property_query(inchikey, 'runtime',
                                              'execution time',
                                              type(start).__name__,
                                              time.time() - start, 's')
    
    def set_parent_path(self, parent_path):
        """Set the parent path (e.g., path to method containing input
        geometry.)"""
        if parent_path > 0:
            self._parent_path_id = parent_path
    
    def has_parent_path(self, inchikey):
        """Returns True if molecule has had entire parent path calculated,
        False otherwise."""
        query = ('SELECT inchikey FROM molecule_method_property WHERE '
                 'inchikey = ? AND method_path_id = ?')
        try:
            self.db.execute(query,
                            (inchikey, self._parent_path_id)).fetchone()[0]
            return True
        except TypeError:
            return False
    
    def check_dependencies(self):
        """If check_dependencies is not implemented, raise error."""
        raise NotImplementedError(("every method needs a 'check_dependencies' "
                                   'method'))
    
    def check(self):
        """If check is not implemented, raise error."""
        # the check method should be called before a calculation (so
        # calculations are not repeated) and after (to verify success)
        raise NotImplementedError("every method needs a 'check' method")
    
    def map(self, inchikey, inchikey_dir):
        """Generally, maps molecule to calculation via method, emits
        query/value pairs.
        """
        raise NotImplementedError(("every method needs a 'map' method"))
    
    def reduce(self, query, values):
        """Run queries/values on the db."""
        total_changes = self.db.total_changes
        if query or values[0]:
            self.db.executemany(query, values)
            self.log_all.info('%i properties added to MESS.DB',
                              self.db.total_changes - total_changes)
            total_changes = self.db.total_changes
    
    def setup(self):
        """Set up method."""
        self._insert_program()
        self._insert_method()
        self._insert_parameters()
