import os
import sys
from datetime import datetime

class AbstractMethod(object):
    # all tools should inherit from this class
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()
        self.setup_properties()
        self.status = 'not set'
    
    def execute(self, args):
        # all methods should have a method that executes its tasks
        # based on the given commands
        raise NotImplementedError( "every method needs an execute method" )
        self.log(args) # all methods should record logs
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
    
    def setup_properties(self):
        # this method sets up properties (populates property table in db)
        # in preparation for molecule-property values
        raise NotImplementedError( "every method needs a setup_properties method" )
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
    
    def insert_property(self, name, description, format):
        return self.c.execute('INSERT OR IGNORE INTO property (name, description, format) VALUES (?, ?, ?);', (name, description, format))
    
    def insert_property_value(self, inchikey, state_id, method_path_id, property_name, property_description, property_format, value, units):
        return self.c.execute('INSERT OR REPLACE INTO molecule_state_method_property (inchikey, state_id, method_path_id, property_id, units, result) SELECT ?, ?, ?, property.property_id, ?, ? FROM property WHERE property.name=? AND property.description=? AND property.format=?', (inchikey, state_id, method_path_id, units, value, property_name, property_description, property_format))