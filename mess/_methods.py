from __future__ import print_function
from __future__ import unicode_literals

import codecs
import os
import sys
from datetime import datetime

from _utils import hash_dict

class AbstractMethod(object):
    """All methods should inherit from this class."""
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()
        self.method_name = self.get_method_name()
        self.status = 'not set'
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
        if not self.is_setup:
            self.setup_method()
            self.db.commit()
            self.set_method_id()
            self.setup_parameters()
            self.db.commit()
            self.is_setup = True

    def check_dependencies(self):
        raise NotImplementedError(("every method needs a 'check_dependencies' "
                                   'method'))
    
    def execute(self, args):
        # all methods should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError("every method needs an 'execute' method")
        self.log(args) # all methods should record logs
        return self.status # should be set by self.check()
    
    def check(self, args):
        # the check method should be called before a calculation (so 
        # calculations are not repeated) and after (to verify success)
        raise NotImplementedError("every method needs a 'check' method")
    
    def log(self, args):
        # the log method should be called at the end of every execute method
        # to record that a calculation has been attempted (in the main log)
        # and to record method-specific messages into the method log
        raise NotImplementedError("every method needs a 'log' method")
    
    def import_properties(self):
        # this method reads molecule-proprty values from calc
        # into db
        raise NotImplementedError(("every method needs an 'import_properties' "
                                   'method'))
    
    def add_messages_to_log(self, log_path, method_name, messages):
        log = codecs.open(log_path, 'a', 'utf-8')
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
        log.close()

    def insert_program(self):
        q = ('INSERT OR IGNORE INTO program (name, version, url) '
             'VALUES (?, ?, ?)')
        return self.c.execute(q, (self.prog_name, self.prog_version, 
                                  self.prog_url))

    def insert_parameter(self, name, setting):
        q = ('INSERT OR IGNORE INTO parameter (name) VALUES (?)')
        r = self.c.execute(q, (name, ))
        q = ('INSERT OR IGNORE INTO method_parameter '
             '(method_id, parameter_id, setting) '
             'SELECT ?, parameter.parameter_id, ? '
             'FROM program, parameter '
             'WHERE parameter.name=?')
        return self.c.execute(q, (self.method_id, setting, name))

    def insert_property(self, inchikey, method_path_id, 
                              name, description, 
                              format, value, units):
        q = ('INSERT OR IGNORE INTO property (name, description, format) '
             'VALUES (?, ?, ?);')
        r = self.c.execute(q, (name, description, format))
        q = ('INSERT OR REPLACE INTO molecule_method_property '
             '(inchikey, method_path_id, property_id, units, result) '
             'SELECT ?, ?, property.property_id, ?, ? '
             'FROM property '
             'WHERE '
             'property.name=? AND property.description=? AND '
             'property.format=?')
        return self.c.execute(q, (inchikey, method_path_id, units, 
                                  value, name, description, format))

    def insert_tags(self):
        q = ('INSERT OR IGNORE INTO method_tag (method_id, parameter_id) '
             'SELECT ?, parameter.parameter_id FROM parameter '
             'WHERE parameter.name= ?')
        for t in self.tags:
            self.c.execute(q, (self.method_id, t))

    def setup_parameters(self):
        for k, v in self.parameters.items():
            self.insert_parameter(k, v)
        self.insert_tags()
    
    def setup_method(self):
        self.insert_program()
        self.param_hash = hash_dict(self.parameters)
        name = self.get_method_name()
        q = ('INSERT OR IGNORE INTO method '
             '(program_id, geop, name, hash) '
             'SELECT program.program_id, ?, ?, ? '
             'FROM program '
             'WHERE program.name=? AND program.version=?')
        return self.c.execute(q, (self.geop, name, self.param_hash, 
                                  self.prog_name, self.prog_version))
    
    def set_method_id(self):
        q = ('SELECT method_id FROM method '
             'WHERE hash = ?;')
        row = self.c.execute(q, (self.param_hash,)).fetchone()
        self.method_id = row.method_id

    @classmethod
    def get_method_name(cls):
        """Return the name of the method, derived from the subclass name."""
        return cls.__name__.replace('_', '').lower()