import os
import sys
from datetime import datetime

class AbstractMethod(object):
    # all tools should inherit from this class
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()
        self.status = 'not set'
        self.is_setup = False
        try:
            self.method_name
            self.method_description
            self.method_level # level of theory (e.g. dft, semiempirical, mm, etc.)
            self.geop # flag indicating whether the method results in new xyz coordinates
            self.prog_name
            self.prog_version
            self.prog_url
        except AttributeError as e:
            sys.stderr.write(''.join([str(e), "\n"]))
            sys.exit('Each method class needs to define method_name, method_description, method_level, geop, prog_name, prog_version, and prog_url as attributes.')
    
    def setup(self):
        if not self.is_setup:
            self.setup_method()
            self.setup_parameters()
            self.db.commit()
            self.is_setup = True
    
    def execute(self, args):
        # all methods should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError( "every method needs an execute method" )
        self.log(args) # all methods should record logs
        self.db.commit()
        return self.status # should be set by self.check()
    
    def check(self, args):
        # the check method should be called before a calculation (so calculations are
        # not repeated) and after (to verify success)
        raise NotImplementedError( "every method needs a check method" )
    
    def log(self, args):
        # the log method should be called at the end of every execute method
        # to record that a calculation has been attempted (in the main log)
        # and to record method-specific messages into the method log
        raise NotImplementedError( "every method needs a log method" )
    
    def setup_parameters(self):
        # this method sets up parameters
        raise NotImplementedError( "every method needs a setup_parameters method" )
        self.db.commit() # make sure you implement commit
    
    def import_properties(self):
        # this method reads molecule-proprty values from calc
        # into db
        raise NotImplementedError( "every method needs an import_properties method" )
        self.db.commit() # make sure you implement commit
    
    def get_inchikey_dir(self, inchikey):
        molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
        return os.path.relpath(os.path.join(molecules_dir, inchikey[:1], inchikey[1:3], inchikey[3:]))
    
    def setup_dir(self, directory):
        if not os.path.isdir(directory):
            os.makedirs(directory)
    
    def add_messages_to_log(self, log_path, method_name, messages):
        log = open(log_path, 'a')
        log.write(": ".join([datetime.now().strftime('%Y-%m-%d %H:%M:%S'), method_name]))
        log.write("\n")
        log.write(" ".join(sys.argv))
        log.write("\n")
        for m in messages:
            log.write(m)
            log.write("\n")
        log.write("-" * 80)
        log.write("\n")
        log.close()
    
    def setup_method(self):
        self.insert_program()
        self.insert_theory_level()
        return self.c.execute('INSERT OR IGNORE INTO method (level_id, program_id, geop, name, description) SELECT level.level_id, program.program_id, ?, ?, ? FROM level, program WHERE level.name=? AND program.name=? AND program.version=?;', (self.geop, self.method_name, self.method_description, self.method_level, self.prog_name, self.prog_version))
            
    def insert_program(self):
        return self.c.execute('INSERT OR IGNORE INTO program (name, version, url) \
            VALUES (?, ?, ?);',
            (self.prog_name, self.prog_version, self.prog_url))
    
    def insert_theory_level(self):
        return self.c.execute('INSERT OR IGNORE INTO level (name) VALUES (?);', (self.method_level,))
    
    def insert_parameter(self, parameter_name, parameter_description):
        return self.c.execute('INSERT OR IGNORE INTO parameter (name, description) VALUES (?, ?);', (parameter_name, parameter_description))
    
    def insert_method_parameter(self, parameter_name, parameter_setting, parameter_description=''):
        self.insert_parameter(parameter_name, parameter_description)
        return self.c.execute('INSERT OR IGNORE INTO method_parameter (method_id, parameter_id, setting) SELECT method.method_id, parameter.parameter_id, ? FROM method, parameter WHERE method.name=? AND parameter.name=?;', (parameter_setting, self.method_name, parameter_name))
    
    def insert_property(self, name, description, format):
        return self.c.execute('INSERT OR IGNORE INTO property (name, description, format) VALUES (?, ?, ?);', (name, description, format))
    
    def insert_property_value(self, inchikey, state_id, method_path_id, property_name, property_description, property_format, value, units):
        self.insert_property(property_name, property_description, property_format)
        return self.c.execute('INSERT OR REPLACE INTO molecule_state_method_property (inchikey, state_id, method_path_id, property_id, units, result) SELECT ?, ?, ?, property.property_id, ?, ? FROM property WHERE property.name=? AND property.description=? AND property.format=?', (inchikey, state_id, method_path_id, units, value, property_name, property_description, property_format))
    
    def get_method_id(self):
        row = self.c.execute('SELECT method_id FROM method WHERE method.name = ?;', (self.method_name,)).fetchone()
        return row['method_id']