"""
Unit tests for methods in mess/_utils.py.
"""

import os
import shutil
import sys
import unittest


class TestUtils(unittest.TestCase):
    def setUp(self):
        sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
        import utils
        self.utils = utils
        self.tmp_dir = './tmp'
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
    
    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
    
    def test_get_inchikey_dir(self):
        get_inchikey_dir = self.utils.get_inchikey_dir
        molecules_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '../../molecules/'))
        test_dir = os.path.join(molecules_dir, 'B/QJ/CRHHNABKAKU-KBQPJGBKSA-N')
        self.assertEqual(get_inchikey_dir('BQJCRHHNABKAKU-KBQPJGBKSA-N'),
                         test_dir)
    
    def test_get_mem_usage(self):
        get_mem_usage = self.utils.get_mem_usage
        self.assertRegexpMatches(get_mem_usage('test'),
                                 '^test.+\d+.+\d+.+\d+.+mb\n$')
    
    def test_hash_dict(self):
        hash_dict = self.utils.hash_dict
        self.assertEqual(hash_dict({'test': 'dictionary'}),
                         'fb2d5e04507b904ea5d45040d9110466397f7b66')
    
    def test_is_inchikey(self):
        is_inchikey = self.utils.is_inchikey
        self.assertTrue(is_inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-N'))
        self.assertTrue(is_inchikey('ADVPTQAUNPRNPO-UHFFFAOYSA-N', True))
        self.assertFalse(is_inchikey('not an inchikey'))
    
    def test_load_method(self):
        load_method = self.utils.load_method
        
        class DB(object):
            def dummy():
                pass
        
        db = DB()
        with self.assertRaises(SystemExit) as context:
            load_method('not_a_valid_method', DB())
        self.assertRegexpMatches(context.exception.message,
                                 'not a valid method')
        sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../..'))
        method = load_method('_import', DB())
        self.assertRegexpMatches(str(method), 'import')
    
    def test_setup_dir(self):
        setup_dir = self.utils.setup_dir
        test_dir = os.path.join(self.tmp_dir, 'test_dir/test_subdir')
        setup_dir(test_dir)
        self.assertTrue(os.path.exists(test_dir))
    
    def test_touch(self):
        touch = self.utils.touch
        test_file = os.path.join(self.tmp_dir, 'temp.txt')
        touch(test_file)
        self.assertTrue(os.path.exists(test_file))
    
    def test_unicode_replace(self):
        unicode_replace = self.utils.unicode_replace
        self.assertTrue(isinstance(unicode_replace('test'), unicode))
    
    def test_write_to_log(self):
        write_to_log = self.utils.write_to_log
        test_log = os.path.join(self.tmp_dir, 'test.log')
        write_to_log(test_log, ['test message'])
        self.assertRegexpMatches(open(test_log).read(), 'test message')
    
    def test_xstr(self):
        xstr = self.utils.xstr
        self.assertEqual(xstr('string'), 'string')
        self.assertEqual(xstr(None), '')


if __name__ == '__main__':
    unittest.main()
