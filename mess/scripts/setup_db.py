#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import os
import sqlite3
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Set up a fresh mess.db.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-rm', '--remove', action='store_true', 
                        help='remove mess.db if it exists, destroying all data')
    args = parser.parse_args()
    mess_db_path = os.path.realpath(os.path.join(os.path.dirname(__file__), 
                                    '../../db/mess.db'))
    if (args.remove):
        try:
            os.remove(mess_db_path)
        except OSError:
            print >> sys.stderr, mess_db_path + (' not deleted because '
                                                 'it already does not exist')
    mess_db_conn = sqlite3.connect(mess_db_path)
    c = mess_db_conn.cursor()
    c.executescript(open(os.path.join(os.path.dirname(__file__), 
                    '../../db/schema.sql')).read())  

if __name__ == '__main__':
    main()