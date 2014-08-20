"""
Unit tests for methods in mess/_db.py.
"""

import os
import shutil
import sys
import unittest

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../..'))
from mess._db import MessDB
from mess.tests.helpers import suppress_stderr


class TestDB(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = './tmp'
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        with suppress_stderr():  # silence 'MESS.DB created' message
            self.db = MessDB(database='%s/test.db' % self.tmp_dir)
    
    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
    
    def test_reopen(self):
        self.db.reopen()
    
    def test_close(self):
        self.db.close()
    
    def test_cursor(self):
        cur = self.db.cursor()
        self.assertRegexpMatches(str(cur), 'sqlite3.Cursor')
    
    def test_insert_select(self):
        insert_query = 'INSERT INTO molecule (inchikey) VALUES (?)'
        self.db.execute(insert_query, ('TEST', ))
        count = 0
        for row in self.db.execute('SELECT * FROM molecule').fetchall():
            self.assertEquals(row.inchikey, 'TEST')
            count += 1
        self.assertEquals(count, 1)
    
    def test_executemany(self):
        insert_query = 'INSERT INTO method_path VALUES (?, ?)'
        values = ((1, 17), (2, 17), (3, 17))
        self.db.executemany(insert_query, values)
        select_query = 'SELECT length FROM method_path'
        count = 0
        for row in self.db.execute(select_query).fetchall():
            self.assertEquals(row.length, 17)
            count += 1
        self.assertEquals(count, 3)
    
    def test_executescript(self):
        pass
    
    def test_check(self):
        self.assertTrue(self.db.check())
        self.db.execute('DROP TABLE molecule')
        self.assertFalse(self.db.check())


if __name__ == '__main__':
    unittest.main()
