#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

import argparse
import csv
import os
import re
import requests
import sqlite3
import sys
import time
from cStringIO import StringIO
from datetime import datetime

import pybel

from _db import MessDB
from _tools import AbstractTool

class Import(AbstractTool):
    def __init__(self):
        self.description = 'imports molecule file, multi-molecule file, or dir of molecules into mess.db'
        self.cir = True
    
    def subparse(self, subparser):
        subparser.add_argument("source", help='A molecule source file or directory.')
    
    def execute(self, args):
        db = MessDB()
        c = db.cursor()
        # make sure inchi method is loaded
        if not os.path.isfile(os.path.join(os.path.dirname(__file__), '../db/methods/inchi.sql')):
            sys.exit('Could not find inchi.sql.')
        c.executescript(open(os.path.join(os.path.dirname(__file__), '../db/methods/inchi.sql')).read())
        # set up source and retrieve source id
        source_row = self.parse_source(c, args.source)
        # iterate over source file or dir and update folder and sql databases
        if os.path.isdir(args.source):
            files = os.listdir(args.source)
        else:
            files = [args.source]
        molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
        pybel.ob.obErrorLog.ClearLog()
        pybel.ob.obErrorLog.StopLogging() 
        for f in files:
            for mol in pybel.readfile(f.split('.')[-1], f):
                ob_logs = []
                # generate inchikey
                inchikey = mol.write('inchikey').rstrip()
                # check if molecule dir exists, create dir structure if not
                inchikey_dir = os.path.join(molecules_dir, inchikey[:1], inchikey[1:3], inchikey[3:])
                inchikey_basename = os.path.join(inchikey_dir, inchikey)
                # set up dir structure
                try:
                    os.makedirs(inchikey_dir)
                except OSError:
                    pass #print >> sys.stderr, inchikey + ' dir exists'
                if not os.path.exists(inchikey_basename + '.inchi'):
                    # catch ob log messages so they can be redirected to inchikey.log later
                    pybel.ob.obErrorLog.StartLogging()
                    mol.write('inchi', inchikey_basename + '.inchi')
                    for i in range(3):
                        ob_logs.append(pybel.ob.obErrorLog.GetMessagesOfLevel(i)) # (0-error, 1-warning, 2-info, 3-audit)
                    pybel.ob.obErrorLog.ClearLog()
                    pybel.ob.obErrorLog.StopLogging() # prevents warnings from being repeated on every conversion
                # capture and clear title (so it doesn't appear in png later)
                identifier = mol.title
                mol.title = ''
                if not os.path.exists(inchikey_basename + '.png'):
                    mol.write('_png2', inchikey_basename + '.png')
                mol.title = identifier
                self.touch(inchikey_basename + '.log')
                self.touch(inchikey_basename + '.notes')
                self.touch(os.path.join(inchikey_dir, 'sources.tsv'))
                # append any messages from inchi conversion to log
                log = open(inchikey_basename + '.log', 'a')
                log.write(": ".join([datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Import']))
                log.write("\n")
                log.write(" ".join(sys.argv))
                log.write("\n")
                for ll in ob_logs:
                    for l in ll:
                        log.write(l)
                log.write("-" * 80)
                log.write("\n")
                log.close()
                # insert molecule to db
                self.update_molecule(c, inchikey, mol)
                # update source catalog in db
                self.update_molecule_source(c, source_row['source_id'], inchikey, identifier)
                # update source list tsv
                with open(os.path.join(inchikey_dir, 'sources.tsv'), 'r') as sources_in:
                    with open(os.path.join(inchikey_dir, 'sources.tsv'), 'a') as sources_out:
                        sources_in = csv.reader(sources_in, delimiter='\t')
                        sources_out = csv.writer(sources_out, delimiter='\t')
                        # check if source has been recorded
                        source_present = False
                        for row in sources_in:
                            try:
                                if (row[1] == source_row['filename']):
                                    source_present = True
                            except IndexError:
                                pass
                        if not source_present:
                            if (source_row['url_template']):
                                url_split = re.split(r"\[|\]", source_row['url_template'])
                                (match, replace) = re.split(r",\s?", url_split[1])
                                url_identifier = re.sub(match, replace, identifier)
                                source_identifier_url = url_split[0] + url_identifier
                                if (2 < len(url_split)):
                                    source_identifier_url += url_split[2]
                            else:
                                source_identifier_url = None
                            sources_out.writerow([source_row['name'], source_row['filename'], identifier, source_identifier_url])   
                # commit changes
                db.commit()

    def touch(self, fname, times=None):
        fhandle = file(fname, 'a')
        try:
            os.utime(fname, times)
        finally:
            fhandle.close()

    def cir_request(self, inchikey, representation):
        if not self.cir:
            return None
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + inchikey + '/' + representation
        time.sleep(0.1) # protect cactus from hammering
        try:
            r = requests.get(url)
            if (r.status_code == 200):
                return r.text
        except requests.ConnectionError:
            print >> sys.stderr, 'CIR is down. Proceeding without importing IUPAC names or other synonyms.'
            self.cir = False
        return None

    def parse_source(self, c, source):
        # check that the provided source is a real file in the right place
        if not (os.path.isdir(source) or os.path.isfile(source)):
            sys.exit('This source is not a valid file or directory.')
        if not (os.path.abspath(os.path.join(source + '/..')).split('/')[-1] == 'sources'):
            sys.exit('All sources must reside in the "sources" directory.')
        # check for a valid source sql file
        source_basename = os.path.basename(source.strip('/').split('/')[-1])
        source_root = os.path.splitext(source_basename)[0]
        source_sql = os.path.join(os.path.dirname(__file__), '../db/sources/', source_root + '.sql')
        if not (os.path.isfile(source_sql)):
            sys.exit('All sources must have a corresponding sql file in "db/sources".')
        # insert/update source in the database
        c.executescript(open(source_sql).read())
        # get source id
        c.execute("SELECT * FROM source WHERE filename=?", (source_basename,))
        return c.fetchone()

    def update_molecule_source(self, c, source_id, inchikey, title):
        c.execute("INSERT OR REPLACE INTO molecule_source VALUES \
                    (?, ?, ?)", \
                    (inchikey, source_id, title))

    def update_molecule(self, c, inchikey, mol):
        # insert/update molecule synonyms
        synonyms = self.cir_request(inchikey, 'names')
        if (synonyms):
            for synonym in (synonyms.split("\n")):
                c.execute("INSERT OR IGNORE INTO molecule_synonym \
                        (inchikey, name) VALUES (?, ?)", \
                        (inchikey, synonym))
        # calculate identifiers with ob/cir, unless entry exists
        c.execute("SELECT inchi FROM molecule WHERE inchikey=?", (inchikey,))
        inchikey_check_row = c.fetchone()
        inchi = mol.write('inchi').rstrip().split('=')[1]
        if inchikey_check_row is not None and inchikey_check_row['inchi'] == inchi:
            return 0 # this molecule is already correct in the db
        smiles = mol.write('can').rstrip() # canonical smiles
        formula = mol.formula
        # get identifiers from CIR
        iupacs = self.cir_request(inchikey, 'iupac_name').splitlines(True)
        iupac = max(iupacs, key=len).rstrip() # if multiple iupacs, take the longest (most specific) one
        if (len(iupacs) > 1): # if multiple iupacs entered, add others as synonym
            for i in iupacs:
                if i != max(mylist, key=len):
                    c.execute("INSERT OR IGNORE INTO molecule_synonym \
                            (inchikey, name) VALUES (?, ?)", \
                            (inchikey, i.rstrip()))
        # insert molecule identifiers
        c.execute("INSERT OR IGNORE INTO molecule \
                (inchikey, inchi, smiles, formula, iupac) \
                VALUES \
                (?, ?, ?, ?, ?)", \
                (inchikey, inchi, smiles, formula, iupac))


def load():
    # loads the current plugin
    return Import()