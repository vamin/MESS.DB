# -*- coding: utf-8 -*-
"""MESS.DB method module

This module contains the AbstractMethod class, which every MESS.DB calculation
method object must inherit from.
"""

from __future__ import print_function
from __future__ import unicode_literals

import codecs
import hashlib
import json
import sys
from datetime import datetime


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
    def __init__(self, db):
        """Set up db, check for attributes, dependencies, and setup.
        
        Args:
            db (obj): A MessDB object
        """
        self.db = db
        self.cur = db.cursor()
        self.method_name = self.get_method_name()
        self.status = ''
        self.is_setup = False
        try:
            self.description
            self.geop  # flag indicates method results in new xyz coordinates
            self.prog_name
            self.prog_version
            self.prog_url
            self.parameters
            self.tags
        except AttributeError as err:
            print(''.join([str(err), '\n']), file=sys.stderr)
            sys.exit(('Each method class needs to define description, geop, '
                      'prog_name, prog_version, prog_url, '
                      'parameters, tags as attributes.'))
        self.check_dependencies()
    
    def setup(self):
        """Set up method."""
        if not self.is_setup:
            self.setup_method()
            self.set_method_id()
            self.setup_parameters()
            self.is_setup = True
    
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
    
    def reduce(self, query, values):
        """Run queries/values on the db."""
        if query or values[0]:
            self.cur.executemany(query, values)
            self.db.commit()
        return len(values)
    
    def add_messages_to_log(self, log_path, method_name, messages):
        """Write messages to a log.
        
        Args:
            log_path: Path to the log to be written to.
            method_name: Name of the method doing the logging.
            messages: List of messages to write to the log.
        """
        with codecs.open(log_path, 'a', 'utf-8') as log:
            log.write(': '.join([datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                 method_name]))
            log.write('\n')
            log.write(' '.join(sys.argv))
            log.write('\n')
            for message in messages:
                log.write(message)
                log.write('\n')
            log.write('-' * 79)
            log.write('\n')
    
    def insert_program(self):
        """Adds row to program table in mess.db."""
        query = ('INSERT OR IGNORE INTO program (name, version, url) '
                 'VALUES (?, ?, ?)')
        self.cur.execute(query,
                         (self.prog_name, self.prog_version, self.prog_url))
        self.db.commit()
    
    def insert_parameter(self, name, setting):
        """Adds parameter to mess.db.
        
        Args:
            name: Name of parameter.
            setting: The value the parameter is set to.
        """
        query = ('INSERT OR IGNORE INTO parameter (name) VALUES (?)')
        self.cur.execute(query, (name, ))
        query = ('INSERT OR IGNORE INTO method_parameter '
                 '(method_id, parameter_id, setting) '
                 'SELECT ?, parameter.parameter_id, ? '
                 'FROM program, parameter '
                 'WHERE parameter.name=?')
        self.cur.execute(query, (self.method_id, setting, name))
        self.db.commit()
    
    def insert_property(self, inchikey, method_path_id, name, description,
                        format_, value, units):
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
        query = ('INSERT OR IGNORE INTO property (name, description, format) '
                 'VALUES (?, ?, ?);')
        self.cur.execute(query, (name, description, format_))
        query = ('INSERT OR REPLACE INTO molecule_method_property '
                 '(inchikey, method_path_id, property_id, units, result) '
                 'SELECT ?, ?, property.property_id, ?, ? '
                 'FROM property '
                 'WHERE '
                 'property.name=? AND property.description=? AND '
                 'property.format=?')
        self.cur.execute(query, (inchikey, method_path_id, units,
                                 value, name, description, format))
        self.db.commit()
    
    def insert_tags(self):
        """Add tags to method_tag table in mess.db."""
        query = ('INSERT OR IGNORE INTO method_tag (method_id, parameter_id) '
                 'SELECT ?, parameter.parameter_id FROM parameter '
                 'WHERE parameter.name= ?')
        for tag in self.tags:
            self.cur.execute(query, (self.method_id, tag))
        self.db.commit()
    
    def insert_property_query(self, name, description, format_):
        """
        Args:
            inchikey: The inchikey of a molecule in MESS.DB.
            method_path_id: Path id for the calculations that generated the
                            property.
            name: The property name.
            description: A description of the property.
            format: A description of the format the property is in.
            value: The calculated property.
            units: Units for the property value.
        """
        query = ('INSERT OR IGNORE INTO property (name, description, format) '
                 'VALUES (?, ?, ?);')
        return (query, (name, description, format_))
    
    def insert_property_value_query(self, inchikey, method_path_id, name,
                                    description, format_, value, units):
        """Insert property values into db."""
        query = ('INSERT OR REPLACE INTO molecule_method_property '
                 '(inchikey, method_path_id, property_id, units, result) '
                 'SELECT ?, ?, property.property_id, ?, ? '
                 'FROM property '
                 'WHERE '
                 'property.name=? AND property.description=? AND '
                 'property.format=?')
        return (query, (inchikey, method_path_id, units,
                        value, name, description, format_))
    
    def setup_parameters(self):
        """Import paramaters dict to mess.db."""
        for k, v in self.parameters.items():
            self.insert_parameter(k, v)
        self.insert_tags()
    
    def setup_method(self):
        """Set insert program to db, set up hash, and insert method to db."""
        self.insert_program()
        self.hash = self.__hash__()
        name = self.get_method_name()
        query = ('INSERT OR IGNORE INTO method '
                 '(program_id, geop, name, hash) '
                 'SELECT program.program_id, ?, ?, ? '
                 'FROM program '
                 'WHERE program.name=? AND program.version=?')
        self.cur.execute(query, (self.geop, name, self.hash,
                                 self.prog_name, self.prog_version))
        self.db.commit()
    
    def set_method_id(self):
        """Set the object's method_id attribute."""
        query = ('SELECT method_id FROM method '
                 'WHERE hash = ?;')
        row = self.cur.execute(query, (self.hash,)).fetchone()
        self.method_id = row.method_id
    
    @classmethod
    def get_method_name(cls):
        """Return the name of the method, derived from the subclass name."""
        return cls.__name__.replace('_', '').lower()
    
    def __hash__(self):
        """Hash based on method name and parameters.
        
        Returns:
            A hex string of the sha1 hash of self.method_name plus
            JSON-serialized self.parameters. Keys are sorted.
        """
        return hashlib.sha1(self.method_name +
                            json.dumps(self.parameters,
                                       sort_keys=True)).hexdigest()
