#!/usr/bin/env python
# encoding: utf-8
# Victor Amin 2013

# imports multi-line inchi to sreeeningdb/molecules dir structure
# expects inchi file from ZINC database

import os, sys
import time
import argparse
import sqlite3
import csv
import requests
import pybel
import re
from datetime import datetime
from cStringIO import StringIO

def main():
    parser = argparse.ArgumentParser(
        description="Import molecules to screening db.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)#,
        #add_help=False)
    parser.add_argument("source", help='A molecule source file or directory.')
    args = parser.parse_args()
    # connect to db
    properties_db_file = os.path.join(os.path.dirname(__file__), '../properties.db')
    if not os.path.isfile(properties_db_file):
        sys.exit('Could not find properties.db.')
    properties_db_conn = sqlite3.connect(properties_db_file)
    properties_db_conn.row_factory = sqlite3.Row
    # set up source and retrieve source id
    source_row = parse_source(properties_db_conn, args.source)
    # iterate over source file or dir and update folder and sql databases
    if os.path.isdir(args.source):
        files = os.listdir(args.source)
    else:
        files = [args.source]
    molecules_dir = os.path.join(os.path.dirname(__file__), '../molecules/')
    pybel.ob.obErrorLog.ClearLog()
    for f in files:
        for mol in pybel.readfile(f.split('.')[-1], f):
            # generate inchikey
            inchikey = mol.write('inchikey').rstrip()
            # catch ob log messages so they can be redirected to inchikey.log later
            ob_logs = []
            for i in range(3):
                ob_logs.append(pybel.ob.obErrorLog.GetMessagesOfLevel(i)) # (0-error, 1-warning, 2-info, 3-audit)
            pybel.ob.obErrorLog.ClearLog()
            pybel.ob.obErrorLog.StopLogging() # prevents warnings from being repeated on every conversion
            # check if molecule dir exists, create dir structure if not
            inchikey_dir = os.path.join(molecules_dir, inchikey[:1], inchikey[1:3], inchikey[3:])
            # capture and clear title (so it doesn't appear in png later)
            identifier = mol.title
            mol.title = ''
            if not os.path.isdir(inchikey_dir):
                # set up dir structure
                os.makedirs(inchikey_dir)
                os.mkdir(os.path.join(inchikey_dir, '_cation'))
                os.mkdir(os.path.join(inchikey_dir, '_anion'))
                os.mkdir(os.path.join(inchikey_dir, '_triplet'))
                os.mkdir(os.path.join(inchikey_dir, '_stripped'))
                os.mkdir(os.path.join(inchikey_dir, '_protonated'))
            if not os.path.isfile(os.path.join(inchikey_dir, 'sources.tsv')):
                touch(os.path.join(inchikey_dir, 'sources.tsv'))
            if not os.path.isfile(os.path.join(inchikey_dir, inchikey + '.inchi')):
                mol.write('inchi', os.path.join(inchikey_dir, inchikey + '.inchi'))
            if not os.path.isfile(os.path.join(inchikey_dir, inchikey + '.png')):
                mol.write('_png2', os.path.join(inchikey_dir, inchikey + '.png'))
            if not os.path.isfile(os.path.join(inchikey_dir, inchikey + '.log')):
                touch(os.path.join(inchikey_dir, inchikey + '.log'))
                # append any messages from inchi conversion to log
                log = open(os.path.join(inchikey_dir, inchikey + '.log'), 'a')
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
            if not os.path.isfile(os.path.join(inchikey_dir, inchikey + '.notes')):
                touch(os.path.join(inchikey_dir, inchikey + '.notes'))
            # insert molecule to db
            update_molecule(properties_db_conn, inchikey, mol)
            # update source catalog in db
            update_molecule_source(properties_db_conn, source_row['source_id'], inchikey, identifier)
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
            properties_db_conn.commit()
            # restart logging for next iteration
            pybel.ob.obErrorLog.StartLogging()
        # close db
        properties_db_conn.close()

def touch(fname, times=None):
    fhandle = file(fname, 'a')
    try:
        os.utime(fname, times)
    finally:
        fhandle.close()

def cir_request(inchikey, representation):
    url = 'http://cactus.nci.nih.gov/chemical/structure/' + inchikey + '/' + representation
    time.sleep(0.1) # protect cactus from hammering
    r = requests.get(url)
    #if (r.status_code != 200):
    #    time.sleep(1)
    #    r = requests.get(url)
    if (r.status_code == 200):
        return r.text
    else:
        return None

def parse_source(conn, source):
    # check that the provided source is a real file in the right place
    if not (os.path.isdir(source) or os.path.isfile(source)):
        sys.exit('This source is not a valid file or directory.')
    if not (os.path.abspath(os.path.join(source + '/..')).split('/')[-1] == 'sources'):
        sys.exit('All sources must reside in the "sources" directory.')
    # check for a valid source sql file
    source_basename = os.path.basename(source.strip('/').split('/')[-1])
    source_root = os.path.splitext(source_basename)[0]
    source_sql = os.path.join(os.path.dirname(__file__), '../sql/sources/', source_root + '.sql')
    if not (os.path.isfile(source_sql)):
        sys.exit('All sources must have a corresponding sql file in "sql/sources".')
    # insert/update source in the database
    c = conn.cursor()
    c.executescript(open(source_sql).read())
    # get source id
    c.execute("SELECT * FROM source WHERE filename=?", (source_basename,))
    return c.fetchone()


def update_molecule_source(conn, source_id, inchikey, title):
    c = conn.cursor()
    c.execute("INSERT OR REPLACE INTO molecule_source VALUES \
                (?, ?, ?)", \
                (inchikey, source_id, title))

def update_molecule(conn, inchikey, mol):
    c = conn.cursor()
    # calculate identifiers with ob
    inchi = mol.write('inchi').rstrip().split('=')[1]
    smiles = mol.write('smi').rstrip()
    formula = mol.formula
    # get identifiers from CIR
    iupac = cir_request(inchikey, 'iupac_name')
    try:
        hashisy = ",".join(cir_request(inchikey, 'hashisy').split()) # can return multiple keys
    except AttributeError:
        hashisy = None
    # insert molecule synonyms
    synonyms = cir_request(inchikey, 'names')
    if (synonyms):
        for synonym in (synonyms.split("\n")):
            c.execute("INSERT OR IGNORE INTO molecule_synonym \
                    (inchikey, name) VALUES (?, ?)", \
                    (inchikey, synonym))
    
    # insert molecule identifiers
    c.execute("INSERT OR IGNORE INTO molecule \
            (inchikey, inchi, smiles, cactvs_hashisy, formula, iupac) \
            VALUES \
            (?, ?, ?, ?, ?, ?)", \
            (inchikey, inchi, smiles, hashisy, formula, iupac))
    

if __name__ == '__main__':
    main()
