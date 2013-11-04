import argparse
import os
import subprocess
import sys
import time
import StringIO
from distutils.version import LooseVersion

from _db import MessDB
from _tools import AbstractTool

class Backup(AbstractTool):
    def __init__(self):
        self.description = 'Backup or restore mess.db and the molecules dir '
        self.epilog = ('Backups consist of the db/mess.db file and the entire '
                       'molecules directory. Everything is tarred and gzipped '
                       'into messdb/backups with the filename '
                       'MESS.DB.timestamp.tgz. Backups can be restored '
                       'directly from thse files, overwriting whatever is '
                       'currently in MESS.DB. Backups are recommended before '
                       'any operation that writes to the database (i.e., '
                       'import, calculate, or remove).')

    
    def subparse(self, subparser):
        subparser.add_argument('-r', 
                               '--restore', 
                               help=('restore from the specified backup'))
        subparser.add_argument('-b', 
                               '--backups_path', 
                               help=('path to deposit backup, default is '
                                     'messdb/backups dir'))

    def check_dependencies(self):
        python_version = ".".join(str(x) for x in sys.version_info)
        if (LooseVersion(python_version) < LooseVersion('2.7')):
                sys.exit("The backup tool requres Python >=2.7.") 
        return True
    
    def execute(self, args):
        mess_db_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                       '../db/mess.db'))
        molecules_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                         '../molecules'))
        backups_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                       '../backups'))

        if (args.restore):
            mess_db_check = 0
            molecules_dir_check = 0
            try:
                output = subprocess.check_output(['tar', '-tzf', 
                                       os.path.relpath(args.restore)])
            except subprocess.CalledProcessError:
                sys.exit(args.restore + ' is not a valid backup file')
            for line in StringIO.StringIO(output):
                split = line.split('/')
                if (split[0].strip() == 'molecules'):
                    molecules_dir_check = 1
                try:
                    if (split[1].strip() == 'mess.db'):
                        mess_db_check = 1
                except IndexError:
                    pass
                if (mess_db_check and molecules_dir_check):
                    subprocess.call(['rm', '-rf', molecules_path])
                    subprocess.call(['tar', '-xzvf', os.path.relpath(args.restore)])
                    sys.exit('***backup restored***')
            sys.exit(args.restore + ' is not a valid backup file')
        else:
            timestamp = str(time.time())
            subprocess.call(['tar', '-czvf', 
                             backups_path + '/' + 'MESS.DB.' + timestamp + 
                             '.tgz', mess_db_path, molecules_path])


def load():
    # loads the current plugin
    return Backup()