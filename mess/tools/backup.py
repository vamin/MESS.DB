# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB backup module

This module contains the backup tool class and load function.
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import inspect
import os
import shutil
import subprocess
import sys
import time
from cStringIO import StringIO

from mess._tool import AbstractTool


class Backup(AbstractTool):
    """This tool backs up all data unique to this MESS.DB instance."""
    
    def __init__(self):
        """Set description of tool."""
        self.description = 'Backup or restore mess.db and the molecules dir '
        self.epilog = ('Backups consist of the db/mess.db file and the entire '
                       'molecules and logs directories. Everything is tarred '
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
        subparser.add_argument('-r', '--restore', metavar='PATH_TO_BACKUP',
                               help=('restore from backup'))
        subparser.add_argument('-b', '--backups-dir', type=str,
                               default=backups_dir,
                               help=('path to deposit backup'))
    
    def execute(self, args):
        """Run backup/restore."""
        db_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                               '../../db'))
        mol_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                '../../molecules'))
        logs_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                 '../../logs'))
        sources_path = os.path.relpath(os.path.join(os.path.dirname(__file__),
                                                    '../../sources'))
        if args.restore:
            self.log_consoleonly.info('validating integrity of %s',
                                      args.restore)
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
                    self.log_consoleonly.info('restore from %s initiated',
                                              args.restore)
                    try:
                        shutil.rmtree(db_path)
                        shutil.rmtree(mol_path)
                        shutil.rmtree(logs_path)
                        shutil.rmtree(sources_path)
                    except OSError:
                        pass  # already deleted
                    subprocess.call(['tar', '-jxvf',
                                     os.path.relpath(args.restore)])
                    self.log_consoleonly.info('backup %s restored',
                                              args.restore)
                    sys.exit()
            sys.exit('%s is not a valid backup file' % args.restore)
        else:
            backup_name = 'MESS.DB.%s.tbz2' % str(time.time())
            self.log_consoleonly.info('backup to %s initiated', backup_name)
            backup_out = os.path.join(args.backups_dir, backup_name)
            if args.max_procs > 1:
                self.check_for_parallel()
                tar = subprocess.Popen(['tar', '-cv',
                                        db_path, mol_path,
                                        logs_path, sources_path],
                                       stdout=subprocess.PIPE)
                subprocess.Popen(['parallel', '--gnu', '--pipe',
                                  '--recend', '',
                                  '--max-procs', str(args.max_procs),
                                  '-k', 'bzip2',
                                  '--best', '-c'], stdin=tar.stdout,
                                 stdout=open(backup_out, 'w'))
                tar.wait()
            else:
                tar = subprocess.Popen(['tar', '-cv',
                                        db_path, mol_path,
                                        logs_path, sources_path],
                                       stdout=subprocess.PIPE)
                subprocess.Popen(['bzip2', '--best', '-c'], stdin=tar.stdout,
                                 stdout=open(backup_out, 'w'))
                tar.wait()
            self.log_consoleonly.info('backup to %s completed', backup_name)
    
    def check_for_parallel(self):
        """Make sure GNU parallel is available on the system, otherwise exit.
        """
        try:
            check = subprocess.check_output(['which', 'parallel'],
                                            stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            self.log_consoleonly.info('backup cancelled')
            sys.exit(("'parallel' needs to be installed "
                      'to use more than one processor '
                      'for compression'))


def load():
    """Load Backup()."""
    return Backup()
