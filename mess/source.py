# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB source module

This module contains the Source class, which is used to interact with and
import molecule source directories.
"""

from __future__ import print_function
from __future__ import unicode_literals

import codecs
import csv
import os
import re
import sys

from mess.log import Log
from mess.utils import get_inchikey_dir


class Source(object):
    """This class provides methods for importing sources and maintaining source
    information in mess.db and in the molecule directories.
    
    Attributes:
        db (obj): A MessDB object
        cur (obj): db.cursor()
        source_dir (str): The user-supplied source directory
        id (int): The source_id in the mess.db source table
        name (str): A name for the source
        dirname (str): The name of the source subdirectory in the 'sources'
                       directory
        url (str): Url where the source can be downloaded
        url_template (str): A url template that can be used to go to the web
                            page for a particular molecule in a source catalog
        last_update (str): Date when source was last downloaded
    """
    def __init__(self, db):
        """Initialize db cursor.
        
        Args:
            db (obj): A MessDB object
        """
        self.db = db
        self.log = Log('all')
    
    def setup(self, source):
        """Setup source in mess.db.
        
        Args:
            source_path: A path to a directory containing source molecules.
        """
        if not os.path.isdir(source):
            self.source_dir = os.path.join(os.path.dirname(__file__), '../',
                                           'sources/', source)
            if not os.path.isdir(self.source_dir):
                sys.exit('%s is not a valid source or source directory.' %
                         source)
        else:
            self.source_dir = source
        source_path_parent = os.path.abspath(os.path.join(self.source_dir,
                                                          '..'))
        if not source_path_parent.split('/')[-1] == 'sources':
            sys.exit("All sources must reside in the 'sources' directory.")
        source_basename = self.source_dir.strip('/').split('/')[-1]
        source_basename = os.path.basename(source_basename)
        source_sql = os.path.join(os.path.splitext(self.source_dir)[0],
                                  '%s.sql' % source_basename)
        if not os.path.isfile(source_sql):
            sys.exit(('All sources must have a corresponding sql file '
                      'in their directory.'))
        # insert/update source in the database
        self.db.executescript(codecs.open(source_sql,
                                          encoding='utf-8').read())
        # TODO if source is new: self.log.info('%s added to sources')
        # get source id
        query = ('SELECT source_id, name, dirname, '
                 'url, url_template, last_update '
                 'FROM source WHERE dirname=?')
        source_row = self.db.execute(query, (source_basename,)).fetchone()
        # set attributes
        self.id = source_row.source_id
        self.name = source_row.name
        self.dirname = source_row.dirname
        self.url = source_row.url
        self.url_template = source_row.url_template
        self.last_update = source_row.last_update
    
    def files(self):
        """Returns a list of files in the source directory."""
        return filter(self.is_source_file, os.listdir(self.source_dir))
    
    def is_source_file(self, file_):
        """Check if file is a valid .sources.tsv file. Returns boolean."""
        return not (file_.split('.')[-1] == 'sql'
                    or file_.split('.')[-1] == 'txt'
                    or file_.split('.')[-1] == 'bak'
                    or file_[-1] == '~')
    
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
                    if self.url_template and 'from:' not in identifier:
                        url_split = re.split(r"\[|\]", self.url_template)
                        (match, replace) = re.split(r",\s?", url_split[1])
                        url_identifier = re.sub(match, replace, identifier)
                        source_identifier_url = url_split[0] + url_identifier
                        if 2 < len(url_split):
                            source_identifier_url += url_split[2]
                    else:
                        source_identifier_url = ''
                    sid_url = source_identifier_url.encode('ascii', 'replace')
                    sources_out.writerow([name, dirname, identifier,
                                          sid_url])
                    self.log.info('%s added to %s sources', identifier,
                                                            inchikey)
