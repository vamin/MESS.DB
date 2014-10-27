# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB source module

This module contains the Source class, which is used to interact with and
import molecule source directories.
"""

from __future__ import print_function
from __future__ import unicode_literals

import codecs
import ConfigParser as cp
import csv
import os
import re
import sys

import pybel

from mess.db import MessDB
from mess.log import Log
from mess.utils import get_inchikey_dir, unicode_replace


class Source(object):
    """This class provides methods for importing sources and maintaining source
    information in mess.db and in the molecule directories.
    
    Attributes:
        db (obj): A MessDB object
        log (obj): A Log('all') object
        source_dir (str): A path to the source directory
        id (int): The source_id in the mess.db source table
        name (str): A name for the source
        dirname (str): The name of the source subdirectory in the 'sources'
                       directory
        url (str): Url where the source can be downloaded
        url_template (str): A url template that can be used to go to the web
                            page for a particular molecule in a source catalog
        last_update (str): Date when source was last downloaded
    """
    def __init__(self):
        """Initialize db cursor.
        
        Args:
            db (obj): A MessDB object
        """
        self.db = MessDB()
        self.log = Log('all')
        self.source_dir = None
        self.id = None
        self.name = None
        self.dirname = None
        self.url = None
        self.url_template = None
        self.last_update = None
    
    @classmethod
    def get_source_dir(cls, source):
        """Returns source directory from user-supplied source, which may be a
        directory or may be the name of a source in the source directory."""
        if not os.path.isdir(source):
            source_dir = os.path.join(os.path.dirname(__file__), '../',
                                      'sources/', source)
            if not os.path.isdir(source_dir):
                sys.exit('%s is not a valid source or source directory.' %
                         source)
        else:
            source_dir = source
        source_path_parent = os.path.abspath(os.path.join(source_dir, '..'))
        if not source_path_parent.split('/')[-1] == 'sources':
            sys.exit("All sources must reside in the 'sources' directory.")
        return os.path.abspath(source_dir)
    
    def files(self):
        """Returns a list of importable files in the source directory."""
        return [f for f in os.listdir(self.source_dir)
                if not (f.startswith('.') or f.startswith('~'))
                and f.split('.')[-1] in pybel.informats]
    
    def setup(self, source):
        """Setup source in mess.db.
        
        Args:
            source: A path to a source directory or a source basename.
        """
        self.source_dir = self.get_source_dir(source)
        # read attributes
        source_basename = os.path.basename(self.source_dir.split('/')[-1])
        source_ini = os.path.join(os.path.splitext(self.source_dir)[0],
                                  '%s.ini' % source_basename)
        if not os.path.isfile(source_ini):
            sys.exit(('All sources must have a corresponding INI file '
                      'in their directory.'))
        source_attributes = self.parse_ini(source_ini)
        # insert/update source in the database
        total_changes = self.db.total_changes
        insert_query = ('INSERT OR IGNORE INTO source '
                        '(name, dirname, url, url_template, '
                        'citation, last_update) '
                        'VALUES (?, ?, null, null, null, null)')
        update_query = ('UPDATE source '
                        'SET url=?, url_template=?, citation=?, last_update=? '
                        'WHERE dirname=?;')
        self.db.execute(insert_query, (source_attributes['name'],
                                       source_basename))
        self.db.execute(update_query, (source_attributes['url'],
                                       source_attributes['url_template'],
                                       source_attributes['citation'],
                                       source_attributes['last_update'],
                                       source_basename))
        if self.db.total_changes - total_changes > 0:
            self.log.info('%s added to sources in database', source_basename)
        select_query = ('SELECT source_id, name, dirname, '
                        'url, url_template, last_update '
                        'FROM source WHERE dirname=?')
        source_row = self.db.execute(select_query,
                                     (source_basename,)).fetchone()
        # set attributes
        self.id = source_row.source_id
        self.name = source_row.name
        self.dirname = source_row.dirname
        self.url = source_row.url
        self.url_template = source_row.url_template
        self.last_update = source_row.last_update
    
    def parse_ini(self, ini):
        """Parse source ini and return attributes."""
        source_attributes = {
            'url': None,
            'url_template': None,
            'citation': None
        }
        config = cp.ConfigParser(dict_type=CaseInsensitiveDict)
        config.read(ini)
        for section in config.sections():
            for option in config.options(section):
                source_attributes[option] = unicode_replace(config.get(section,
                                                                       option))
        required_attributes = ('name', 'last_update')
        if not all(att in source_attributes for att in required_attributes):
            sys.exit('Source INI missing required attributes: %s.'
                     % ' and/or '.join(required_attributes))
        return source_attributes
    
    def update_molecule_source_query(self, inchikey, identifier):
        """Update the source in mess.db.
        
        Args:
            inchikey: A molecule InChIKey.
            identifier: A source identifier (usually a catalog number).
        """
        query = ('INSERT OR IGNORE INTO molecule_source '
                 '(inchikey, source_id, identifier) '
                 'VALUES (?, ?, ?)')
        return (query, (inchikey, self.id, identifier))
    
    def update_source_tsv(self, inchikey, identifier):
        """Update the sources.tsv file.
        
        Args:
            inchikey_dir: Dir to a molecule in the molecules dir.
            identifier: A source identifier (usually a catalog number).
        
        """
        inchikey_dir = get_inchikey_dir(inchikey)
        name = self.name.encode('ascii', 'replace')
        dirname = self.dirname.encode('ascii', 'replace')
        identifier = identifier.encode('ascii', 'replace')
        sources_tsv = os.path.join(inchikey_dir, '%s.sources.tsv' % inchikey)
        with codecs.open(sources_tsv, 'r', 'ascii') as sources_in:
            with codecs.open(sources_tsv, 'a', 'ascii') as sources_out:
                sources_in = csv.reader(sources_in, delimiter=b'\t')
                sources_out = csv.writer(sources_out, delimiter=b'\t')
                # check if source has been recorded
                source_present = False
                for row in sources_in:
                    try:
                        if row[1] == dirname and row[2] == identifier:
                            source_present = True
                    except IndexError:
                        pass
                if not source_present:
                    if self.url_template:
                        url_split = re.split(r"\[|\]", self.url_template)
                        (match, replace) = re.split(r",\s?", url_split[1])
                        url_identifier = re.sub(match, replace, identifier)
                        source_url = url_split[0] + url_identifier
                        if 2 < len(url_split):
                            source_url += url_split[2]
                    else:
                        source_url = ''
                    sources_out.writerow([name, dirname, identifier,
                                          source_url.encode('ascii',
                                                            'replace')])
                    self.log.info('%s added to %s sources', name, inchikey)


class CaseInsensitiveDict(dict):
    """Case insensitive dict."""
    def __setitem__(self, key, value):
        """Make all keys lowercase when setting."""
        super(CaseInsensitiveDict, self).__setitem__(key.lower(), value)

    def __getitem__(self, key):
        """Request key as lowercase when getting."""
        return super(CaseInsensitiveDict, self).__getitem__(key.lower())
