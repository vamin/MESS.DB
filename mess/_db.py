import os
import sqlite3
import sys

class MessDB(object):
    def __init__(self):
        # connect to db
        try:
            self.conn = sqlite3.connect(os.path.join(os.path.dirname(__file__), 
                                        '../db/mess.db'))
        except IOError:
            sys.exit('could not find mess.db')
        self.conn.row_factory = sqlite3.Row
    
    def cursor(self):
        return self.conn.cursor()
    
    def commit(self):
        return self.conn.commit()
    
    def __del__(self):
        self.conn.close()