"""
Unit tests for methods in mess/_utils.py.
"""

import os
import shutil
import sys
import unittest

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../..'))
import mess.utils as utils
from mess.tests.helpers import suppress_stderr


class TestUtils(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = './tmp'
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
    
    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
    
    def test_get_inchikey_dir(self):
        molecules_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '../../molecules/'))
        test_dir = os.path.join(molecules_dir, 'B/QJ/CRHHNABKAKU-KBQPJGBKSA-N')
        self.assertEqual(utils.get_inchikey_dir('BQJCRHHNABKAKU-KBQPJGBKSA-N'),
                         test_dir)
    
    def test_get_mem_usage(self):
        self.assertRegexpMatches(utils.get_mem_usage(),
                                 '.+\d+.+\d+.+\d+.+mb\n$')
    
    def test_hash_dict(self):
        self.assertEqual(utils.hash_dict({'test': 'dictionary'}),
                         'fb2d5e04507b904ea5d45040d9110466397f7b66')
    
    def test_is_inchikey(self):
        self.assertTrue(utils.is_inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-N'))
        self.assertTrue(utils.is_inchikey('ADVPTQAUNPRNPO-UHFFFAOYSA-N', True))
        self.assertFalse(utils.is_inchikey('not an inchikey'))
    
    def test_load_method(self):
        
        class DB(object):
            def dummy():
                pass
        
        db = DB()
        with self.assertRaises(SystemExit) as context:
            with suppress_stderr():
                utils.load_method('not_a_valid_method', DB())
        self.assertRegexpMatches(context.exception.message,
                                 'not a valid method')
        sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../..'))
        method = utils.load_method('_import', DB())
        self.assertRegexpMatches(str(method), 'import')
    
    def test_setup_dir(self):
        test_dir = os.path.join(self.tmp_dir, 'test_dir/test_subdir')
        utils.setup_dir(test_dir)
        self.assertTrue(os.path.exists(test_dir))
    
    def test_touch(self):
        test_file = os.path.join(self.tmp_dir, 'temp.txt')
        utils.touch(test_file)
        self.assertTrue(os.path.exists(test_file))
    
    def test_unicode_replace(self):
        self.assertTrue(isinstance(utils.unicode_replace(b'test'), unicode))
        self.assertTrue(isinstance(utils.unicode_replace(u'test'), unicode))
    
    def test_write_to_log(self):
        test_log = os.path.join(self.tmp_dir, 'test.log')
        utils.write_to_log(test_log, ['test message'])
        self.assertRegexpMatches(open(test_log).read(), 'test message')
    
    def test_xstr(self):
        self.assertEqual(utils.xstr('string'), 'string')
        self.assertEqual(utils.xstr(None), '')


if __name__ == '__main__':
    unittest.main()
