from __future__ import print_function
from __future__ import unicode_literals

import os
import subprocess
import sys
import time
from distutils.version import LooseVersion
from StringIO import StringIO

from _db import MessDB
from _tools import AbstractTool

class Backup(AbstractTool):
    def __init__(self):
        """Set description of tool."""
        self.description = 'Backup or restore mess.db and the molecules dir '
        self.epilog = ('Backups consist of the db/mess.db file and the entire '
                       'molecules  and logs directories. Everything is tarred '
                       'and bzipped into messdb/backups with the filename '
                       'MESS.DB.timestamp.tbz2. Backups can be restored '
                       'directly from thse files, overwriting whatever is '
                       'currently in MESS.DB. Backups are recommended before '
                       'any operation that writes to the database (i.e., '
                       'import, calculate, or remove).')
    
    def subparse(self, subparser):
        """Set tool-specific argparse arguments."""
        subparser.add_argument('-p', '--max-procs', type=int, default=1,
                               help=('Number of CPUs to use for compression, '
                                     'default is 1'))
        subparser.add_argument('-r', '--restore',
                               help=('Restore from the specified backup'))
        subparser.add_argument('-b', '--backups_path',
                               help=('Path to deposit backup, default is '
                                     'messdb/backups dir'))
    
    def check_dependencies(self):
        """Check for dependencies (Python >=2.7).
        
        Returns:
            True if all dependencies are met.
        
        """
        python_version = '.'.join(str(x) for x in sys.version_info)
        if (LooseVersion(python_version) < LooseVersion('2.7')):
            sys.exit('The backup tool requres Python >=2.7.')
        try:
            if (LooseVersion(pybel.ob.OBReleaseVersion()) <
                LooseVersion('2.3.0')):
                sys.exit(('This tool requires Open Babel (and its python '
                          'module, pybel) version >=2.3.0.'))
        except AttributeError, OSError:
            sys.exit(('This tool requires Open Babel (and its python module, '
                      'pybel) version >=2.3.0.'))
        return True
    
    def execute(self, args):
        """Run backup/restore."""
        mess_db_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                       '../db/mess.db'))
        mol_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                         '../molecules'))
        logs_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                         '../logs'))
        backups_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                       '../backups'))
        
        if (args.restore):
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
                if (split[0].strip() == 'molecules'):
                    molecules_dir_check = 1
                try:
                    if (split[1].strip() == 'mess.db'):
                        mess_db_check = 1
                except IndexError:
                    pass
                if (mess_db_check and molecules_dir_check):
                    print('***restore from %s initiated***' % args.restore, 
                          file=sys.stderr)
                    subprocess.call(['rm', '-rf', mol_path])
                    subprocess.call(['rm', '-rf', logs_path])
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