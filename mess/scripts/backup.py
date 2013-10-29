import argparse
import os
import subprocess
import sys
import time
import StringIO

def main():
    parser = argparse.ArgumentParser(
        description='Backup (or restore) a mess.db/molecules dir pair.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--restore', help=('restore from the specified '
                                                 'mess.db/molecules dir '
                                                 'tgz file'))
    args = parser.parse_args()
    mess_db_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                   '../../db/mess.db'))
    molecules_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                     '../../molecules'))
    backups_path = os.path.relpath(os.path.join(os.path.dirname(__file__), 
                                   '../../backups'))

    if (args.restore):
        mess_db_check = 0
        molecules_dir_check = 0
        try:
            output = check_output(['tar', '-tzf', 
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
                         backups_path + '/' + 'messdb.' + timestamp + '.tgz', 
                         mess_db_path, molecules_path])

def check_output(*popenargs, **kwargs):
    """Run command with arguments and return its output as a byte string.

    Backported from Python 2.7 as it's implemented as pure python on stdlib.

    >>> check_output(['/usr/bin/python', '--version'])
    Python 2.6.2
    """
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

if __name__ == '__main__':
    main()