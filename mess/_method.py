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


class AbstractMethod(object):
    """All methods must inherit from this class.
    
    Attributes:
        db (obj): A MessDB object
        cur (obj): db.cursor()
        method_name (str): The name of the method
        status (str): The current status of the method
        is_setup (bool): Indicates whether setup() has been run
        description (str): Description of method
        geop (bool): Whether the method generates a new geometry
        prog_name (str): Program name
        prog_version (str): Program version
        prog_url (str): Program url
        parameters (dict): Parameters that affect program execution
        tags (list of str): List of parameters that identify the method
    """
    parameters = dict()
    tags = list()
    status = None
    is_setup = False
    
    def __init__(self, db):
        """Set up db, check for attributes, dependencies, and setup.
        
        Args:
            db (obj): A MessDB object
        """
        self.db = db
        self.method_name = self.get_method_name()
        try:
            self.description
            self.geop  # flag indicates method results in new xyz coordinates
            self.prog_name
            self.prog_version
            self.prog_url
        except AttributeError as err:
            print(''.join([str(err), '\n']), file=sys.stderr)
            sys.exit(('Each method class needs to define description, geop, '
                      'prog_name, prog_version, prog_url, '
                      'parameters, tags as attributes.'))
        self.check_dependencies()
    
    def __hash__(self):
        """Hash based on method name and parameters.
        
        Returns:
            A hex string of the sha1 hash of self.method_name plus
            JSON-serialized self.parameters. Keys are sorted.
        """
        return hashlib.sha1(self.method_name +
                            json.dumps(self.parameters,
                                       sort_keys=True)).hexdigest()
    
    def check_dependencies(self):
        """If check_dependencies is not implemented, raise error."""
        raise NotImplementedError(("every method needs a 'check_dependencies' "
                                   'method'))
    
    def check(self, args):
        """If check is not implemented, raise error."""
        # the check method should be called before a calculation (so
        # calculations are not repeated) and after (to verify success)
        raise NotImplementedError("every method needs a 'check' method")
    
    def log(self, args):
        """If log is not implemented, raise error."""
        # the log method should be called at the end of every execute method
        # to record that a calculation has been attempted (in the main log)
        # and to record method-specific messages into the method log
        raise NotImplementedError("every method needs a 'log' method")
    
    def import_properties(self):
        """If import_properties is not implemented, raise error."""
        # this method reads molecule-proprty values from calc
        # into db
        raise NotImplementedError(("every method needs an 'import_properties' "
                                   'method'))
    
    def map():
        """Generally, maps molecule to calculation via method, emits
        query/value pairs.
        """
        raise NotImplementedError(("every method needs a 'map' method"))
    
    def reduce(self, query, values):
        """Run queries/values on the db."""
        if query or values[0]:
            self.db.executemany(query, values)
    
    def setup(self):
        """Set up method."""
        if not self.is_setup:
            self.insert_program()
            self.insert_method()
            self.insert_parameters()
            self.insert_tags()
            self.is_setup = True
    
    def insert_method(self): # WITH:
        """Set insert program to db, set up hash, and insert method to db."""
        query = ('INSERT OR IGNORE INTO method '
                 '(program_id, geop, name, hash) '
                 'SELECT program.program_id, ?, ?, ? '
                 'FROM program '
                 'WHERE program.name=? AND program.version=?')
        self.db.execute(query, (self.geop, self.method_name, self.hash,
                                self.prog_name, self.prog_version))
    
    def insert_program(self):
        """Adds row to program table in mess.db."""
        query = ('INSERT OR IGNORE INTO program (name, version, url) '
                 'VALUES (?, ?, ?)')
        self.db.execute(query,
                        (self.prog_name, self.prog_version, self.prog_url))
    
    def insert_parameters(self):
        """Import paramaters dict to mess.db.
        
        Args:
            name: Name of parameter.
            setting: The value the parameter is set to.
        """
        for name, setting in self.parameters.items():
            query = ('INSERT OR IGNORE INTO parameter (name) VALUES (?)')
            self.db.execute(query, (name, ))
            query = ('INSERT OR IGNORE INTO method_parameter '
                     '(method_id, parameter_id, setting) '
                     'SELECT ?, parameter.parameter_id, ? '
                     'FROM program, parameter '
                     'WHERE parameter.name=?')
            self.db.execute(query, (self.method_id, setting, name))
    
    def insert_tags(self):
        """Add tags to method_tag table in mess.db."""
        query = ('INSERT OR IGNORE INTO method_tag (method_id, parameter_id) '
                 'SELECT ?, parameter.parameter_id FROM parameter '
                 'WHERE parameter.name= ?')
        for tag in self.tags:
            self.db.execute(query, (self.method_id, tag))
    
    def get_insert_property_query(self, inchikey, method_path_id, name,
                                  description, format_, value, units):
        """Adds property value to mess.db.
        
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
        query = ('INSERT INTO molecule_method_property_denorm '
                 'VALUES (?, ?, ?, ?, ?, ?, ?);')
        return (query, (inchikey, method_path_id, name, description,
                        format_, units, value))
    
    def get_timing_query(self, inchikey, path_id, start):
        return self.get_insert_property_query(inchikey, path_id, 'runtime',
                                              'execution time',
                                              type(start).__name__,
                                              time.time()-start, 's')
    
    @classmethod
    def get_method_name(cls):
        """Return the name of the method, derived from the subclass name."""
        return cls.__name__.replace('_', '').lower()
    
    @property
    def method_id(self):
        """Get the object's method_id attribute."""
        query = ('SELECT method_id FROM method '
                 'WHERE hash = ?;')
        row = self.db.execute(query, (self.hash,)).fetchone()
        return row.method_id
    
    @property
    def hash(self):
        return self.__hash__()
