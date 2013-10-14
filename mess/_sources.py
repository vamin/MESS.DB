import csv
import os
import re
import sqlite3
import sys

class Source(object):
    def __init__(self, db):
        self.db = db
        self.c = db.cursor()
            
    def setup(self, source_path):
        self.source = source_path
        # check that the provided source is a real file in the right place
        if not (os.path.isdir(self.source)):
            sys.exit('This source is not a valid directory.')
        if not (os.path.abspath(os.path.join(self.source + '/..')).split('/')[-1] == 'sources'):
            sys.exit("All sources must reside in the 'sources' directory.")
        source_basename = os.path.basename(self.source.strip('/').split('/')[-1])
        source_sql = os.path.join(os.path.splitext(self.source)[0], source_basename + '.sql')
        if not (os.path.isfile(source_sql)):
            sys.exit('All sources must have a corresponding sql file in their directory.')
        # insert/update source in the database
        self.c.executescript(open(source_sql).read())
        # get source id
        source_row = self.c.execute("SELECT source_id, name, dirname, url, url_template, last_update FROM source WHERE dirname=?", (source_basename,)).fetchone()
        # set attributes
        self.id = source_row['source_id']
        self.name = source_row['name']
        self.dirname = source_row['dirname']
        self.url = source_row['url']
        self.url_template = source_row['url_template']
        self.last_update = source_row['last_update']
    
    def files(self):
        return os.listdir(self.source)
    
    def update_molecule_source(self, inchikey, identifier):
        return self.c.execute("INSERT OR REPLACE INTO molecule_source \
                    (inchikey, source_id, identifier) VALUES \
                    (?, ?, ?)", \
                    (inchikey, self.id, identifier))
    
    def update_source_tsv(self, inchikey_dir, identifier):
        with open(os.path.join(inchikey_dir, 'sources.tsv'), 'r') as sources_in:
            with open(os.path.join(inchikey_dir, 'sources.tsv'), 'a') as sources_out:
                sources_in = csv.reader(sources_in, delimiter='\t')
                sources_out = csv.writer(sources_out, delimiter='\t')
                # check if source has been recorded
                source_present = False
                for row in sources_in:
                    try:
                        if (row[1] == self.dirname):
                            source_present = True
                    except IndexError:
                        pass
                if not source_present:
                    if (self.url_template):
                        url_split = re.split(r"\[|\]", self.url_template)
                        (match, replace) = re.split(r",\s?", url_split[1])
                        url_identifier = re.sub(match, replace, identifier)
                        source_identifier_url = url_split[0] + url_identifier
                        if (2 < len(url_split)):
                            source_identifier_url += url_split[2]
                    else:
                        source_identifier_url = None
                    sources_out.writerow([self.name, self.dirname, identifier, source_identifier_url])