from __future__ import print_function
from __future__ import unicode_literals

import codecs
import hashlib
import json
import os
import resource
import sys
from datetime import datetime


def get_inchikey_dir(inchikey):
    """Convert InChIKey into a path.
    
    Args:
        inchikey: An InChIKey string.
    
    Returns:
        Absolute path of the form e.g.:
        path/to/molecules/C/PE/LXLSAUQHCOX-UHFFFAOYSA-M/
    
    """
    molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
    return os.path.abspath(os.path.join(molecules_dir, inchikey[:1],
                                        inchikey[1:3], inchikey[3:]))


def get_mem_usage(point):
    """Get the current memory usage.
    
    Args:
        point: A name for the point in time memory usage is being requested.
    
    Returns:
        A human-readable string.
    
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return ('%s: usertime=%s systime=%s '
            'mem=%s mb\n') % (point, usage[0],usage[1],
                              (usage[2] * resource.getpagesize()) / 1000000.0)


def hash_dict(d):
    """Serialize and hash a dict.
    
    Args:
        d: A dict.
    
    Returns:
        A hex string of the sha1 hash of the JSON-serialized dict. Keys are
        sorted.
    
    """
    return hashlib.sha1(json.dumps(d, sort_keys=True)).hexdigest()


def is_inchikey(inchikey, enforce_standard=False):
    """Check if a string is a valid InChIKey.
    
    Args:
        inchikey: A supposed InChIkey string.
        enfore_standard: Make sure InChIKey is "standard". Default: False.
    
    Returns:
        boolean
    
    """
    if '=' in inchikey:
        inchikey = inchikey.split('=')[1]
    if len(inchikey) == 27:
        s = inchikey.split('-')
        try:
            if len(s[0]) == 14 and len(s[1]) == 10 and len(s[2]) == 1:
                if s[0].isalpha() and s[1].isalpha() and s[2].isalpha():
                    if not enforce_standard or s[1][-2] == 'S':
                        return True
        except IndexError:
            pass
    return False


def load_method(method_name, db):
    """Locate a method in mess/methods and return an instance of it."""
    try:
        module = __import__('mess.methods.%s' % method_name,
                            fromlist=['methods'])
        method = module.load(db)
    except ImportError as err:
        print('Error: %s;' % err)
        sys.exit('\'%s\' is not a valid method.' % method_name)
    return method


def setup_dir(directory):
    """If directory does not exist, create it and its parents."""
    if not os.path.isdir(directory):
        os.makedirs(directory)


def touch(fname, times=None):
    """Update the timestamp on a file."""
    fhandle = file(fname, 'a')
    try:
        os.utime(fname, times)
    finally:
        fhandle.close()


def unicode_replace(x, enc='utf-8', err='replace'):
    """Convert str to unicode.
    
    Args:
        x: A string.
        enc: Encoding of input string, defaults to 'utf-8'.
        err: What to do on unicode conversion error, defaults to 'replace'.
    
    Returns:
        Unicode string if x is str, x otherwise.
    
    """
    if isinstance(x, str):
        return unicode(x, enc, err)
    else:
        return x


def write_to_log(log_path, messages):
    """Write messages to a log.
    
    Args:
        log_path: Path to the log to be written to.
        messages: List of messages to write to the log.
    """
    with codecs.open(log_path, 'a', 'utf-8') as log:
        log.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        log.write('\n')
        log.write(' '.join(sys.argv))
        log.write('\n')
        try:
            for message in messages:
                log.write(message)
                log.write('\n')
        except TypeError:
            log.write('')
            log.write('\n')
        log.write('-' * 79)
        log.write('\n')


def xstr(s):
    """Return str(), except that None returns empty string."""
    if s is None:
        return ''
    return str(s)
