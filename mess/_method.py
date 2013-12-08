from __future__ import print_function
from __future__ import unicode_literals

import codecs
import hashlib
import json
import os
import sys
from datetime import datetime

class AbstractMethod(object):
    """All methods should inherit from this class."""
    def __init__(self, db):
        """Set up db, check for attributes, dependencies, and setup."""
        self.db = db
        self.c = db.cursor()
        self.method_name = self.get_method_name()
        self.status = ''
        self.is_setup = False
        try:
            self.description
            self.geop # flag indicates method results in new xyz coordinates
            self.prog_name
            self.prog_version
            self.prog_url
            self.parameters
            self.tags
        except AttributeError as e:
            print(''.join([str(e), '\n']), file=sys.stderr)
            sys.exit(('Each method class needs to define description, geop, '
                      'prog_name, prog_version, prog_url, '
                      'parameters, tags as attributes.'))
        self.check_dependencies()
        self.setup()
    
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
    
    def execute(self, args):
        """If execute is not implemented, raise error."""
        # all methods should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError("every method needs an 'execute' method")
        self.log(args) # all methods should record logs
        return self.status # should be set by self.check()
    
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
            for m in messages:
                log.write(m)
                log.write('\n')
            log.write('-' * 79)
            log.write('\n')
    
    def insert_program(self):
        """Adds row to program table in mess.db."""
        q = ('INSERT OR IGNORE INTO program (name, version, url) '
             'VALUES (?, ?, ?)')
        self.c.execute(q, (self.prog_name, self.prog_version, self.prog_url))
        self.db.commit()
    
    def insert_parameter(self, name, setting):
        """Adds parameter to mess.db.
        
        Args:
            name: Name of parameter.
            setting: The value the parameter is set to.
        
        """
        q = ('INSERT OR IGNORE INTO parameter (name) VALUES (?)')
        self.c.execute(q, (name, ))
        q = ('INSERT OR IGNORE INTO method_parameter '
             '(method_id, parameter_id, setting) '
             'SELECT ?, parameter.parameter_id, ? '
             'FROM program, parameter '
             'WHERE parameter.name=?')
        self.c.execute(q, (self.method_id, setting, name))
        self.db.commit()
    
    def insert_property(self, inchikey, method_path_id,
                              name, description,
                              format, value, units):
        """Adds property value to mess.db.
        
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
        q = ('INSERT OR IGNORE INTO property (name, description, format) '
             'VALUES (?, ?, ?);')
        self.c.execute(q, (name, description, format))
        q = ('INSERT OR REPLACE INTO molecule_method_property '
             '(inchikey, method_path_id, property_id, units, result) '
             'SELECT ?, ?, property.property_id, ?, ? '
             'FROM property '
             'WHERE '
             'property.name=? AND property.description=? AND '
             'property.format=?')
        self.c.execute(q, (inchikey, method_path_id, units, 
                           value, name, description, format))
        self.db.commit()
    
    def insert_tags(self):
        """Add tags to method_tag table in mess.db."""
        q = ('INSERT OR IGNORE INTO method_tag (method_id, parameter_id) '
             'SELECT ?, parameter.parameter_id FROM parameter '
             'WHERE parameter.name= ?')
        for t in self.tags:
            self.c.execute(q, (self.method_id, t))
        self.db.commit()
    
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
        q = ('INSERT OR IGNORE INTO method '
             '(program_id, geop, name, hash) '
             'SELECT program.program_id, ?, ?, ? '
             'FROM program '
             'WHERE program.name=? AND program.version=?')
        self.c.execute(q, (self.geop, name, self.hash,
                           self.prog_name, self.prog_version))
        self.db.commit()
    
    def set_method_id(self):
        """Set the object's method_id attribute."""
        q = ('SELECT method_id FROM method '
             'WHERE hash = ?;')
        row = self.c.execute(q, (self.hash,)).fetchone()
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