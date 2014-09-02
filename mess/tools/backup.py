from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import glob
import inspect
import os
import shutil
import subprocess
import sys
import time
from cStringIO import StringIO

from mess._tool import AbstractTool

class Backup(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Backup or restore mess.db and the molecules dir '
        self.epilog = ('Backups consist of the db/mess.db file and the entire '
                       'molecules  and logs directories. Everything is tarred '
                       'and bzipped into messdb/backups with the filename '
                       'MESS.DB.timestamp.tbz2. Backups can be restored '
                       'directly from these files, overwriting whatever is '
                       'currently in MESS.DB. Backups are recommended before '
                       'any operation that writes to the database (i.e., '
                       'import, calculate, or remove).')
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        tools_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
        backups_dir = os.path.realpath(os.path.join(tools_dir,
                                                    '../../backups'))
        subparser.add_argument('-p', '--max-procs', type=int, default=1,
                               help=('number of CPUs to use for compression'))
        subparser.add_argument('-r', '--restore-path',
                               help=('restore from backup'))
        subparser.add_argument('-b', '--backups-dir', type=str,
                               default=backups_dir,
                               help=('path to deposit backup'))
    
    def execute(self, args):
        """Run backup/restore."""
        mess_db_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                    '../../db/mess.db'))
        mol_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                '../../molecules'))
        logs_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                 '../../logs'))
        backups_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                    '../../backups'))
        
        if args.restore:
            print('***validating integrity of %s***' % args.restore,
                  file=sys.stderr)
            mess_db_check = 0
            molecules_dir_check = 0
            try:
                output = subprocess.check_output(['tar', '-tjf', args.restore])
            except subprocess.CalledProcessError:
                sys.exit('%s is not a valid backup file' % args.restore)
            for line in StringIO(output):
                split = line.split('/')
                if split[0].strip() == 'molecules':
                    molecules_dir_check = 1
                try:
                    if split[1].strip() == 'mess.db':
                        mess_db_check = 1
                except IndexError:
                    pass
                if mess_db_check and molecules_dir_check:
                    print('***restore from %s initiated***' % args.restore,
                          file=sys.stderr)
                    try:
                        shutil.rmtree(mol_path)
                        shutil.rmtree(logs_path)
                    except OSError:
                        # already deleted
                        pass
                    mess_db_glob = glob.glob('%s*' % mess_db_path)
                    for file in mess_db_glob:
                        os.remove(file)
                    subprocess.call(['tar', '-jxvf',
                                     os.path.relpath(args.restore)])
                    sys.exit('***backup %s restored***' % args.restore)
            sys.exit('%s is not a valid backup file' % args.restore)
        else:
            backup_name = 'MESS.DB.%s.tbz2' % str(time.time())
            print('***backup to %s initiated***' % backup_name,
                  file=sys.stderr)
            tar = subprocess.Popen(['tar', '-cv', mess_db_path, mol_path,
                                    logs_path], stdout=subprocess.PIPE)
            subprocess.Popen(['parallel', '--gnu', '--pipe', '--recend', '',
                              '--max-procs', str(args.max_procs), '-k',
                              'bzip2', '--best', '-c'], stdin=tar.stdout,
                             stdout=open(os.path.join(backups_path,
                                                      backup_name), 'w'))
            tar.wait()
            print('***backup to %s completed***' % backup_name,
                  file=sys.stderr)


def load():
    """Load Backup()."""
    return Backup()
