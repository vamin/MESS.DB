# MESS.DB - Molecular Electronic Structure (and other Stuff) DB

This project provides a framework for organizing electronic structure (and other) calculations on organic molecules into a rational file structure and relational database.  

## Motivation and Design Philosophy

MESS.DB is designed to facilitate large scale screening. Molecules are represented in the database and file structure by their InChiKeys. MESS.DB is intended to be simple, portable, and human-friendly.

## How it Works
Suppose you import the molecule Morphine. Morphine's InChiKey will be calculated (BQJCRHHNABKAKU-KBQPJGBKSA-N) and a directory will be created for it:

molecules/B/QJ/CRHHNABKAKU-KBQPJGBKSA-N/
    BQJCRHHNABKAKU-KBQPJGBKSA-N.inchi <- the molecule in InChi format
    BQJCRHHNABKAKU-KBQPJGBKSA-N.log <- a log tracking what has been done to the molecule
    BQJCRHHNABKAKU-KBQPJGBKSA-N.notes <- a blank space for notes
    BQJCRHHNABKAKU-KBQPJGBKSA-N.png <- a 2D representation of the molecule
    sources.tsv <- a table of sources for the molecule, including where to buy if the source is commercial

In addition, morphine, along with its smiles, inchi, IUPAC name, synonyms, and basic properties (like MW, charge, etc.) will be imported to mess.db, an SQLite relational database. For the curious, the schema is in db/schema.sql.

Methods (which, as far as MESS is concerned, are plugins that describe how to run a particular calculation) can be run against the database (or a subset on it). If I apply the balloon141 method, which genereates 3D structures from smiles strings, a new folder appears in the molecules folder:

molecules/B/QJ/CRHHNABKAKU-KBQPJGBKSA-N/
    balloon141_FROM_import_PATH_2/ <- contains logs and output from running balloon
    BQJCRHHNABKAKU-KBQPJGBKSA-N.inchi
    BQJCRHHNABKAKU-KBQPJGBKSA-N.log
    BQJCRHHNABKAKU-KBQPJGBKSA-N.notes
    BQJCRHHNABKAKU-KBQPJGBKSA-N.png
    sources.tsv

If balloon generates any new properties that are not in the database, they are added. Now we can use the balloon 3D coordinates to run another calculation, and get:

molecules/B/QJ/CRHHNABKAKU-KBQPJGBKSA-N/
    balloon141_FROM_import_PATH_2/ <- contains logs and output from running balloon
    pm7_mopac2012_FROM_balloon141_PATH_3/ <- contains logs and output from running mopac
    BQJCRHHNABKAKU-KBQPJGBKSA-N.inchi
    BQJCRHHNABKAKU-KBQPJGBKSA-N.log
    BQJCRHHNABKAKU-KBQPJGBKSA-N.notes
    BQJCRHHNABKAKU-KBQPJGBKSA-N.png
    sources.tsv

Even though most relevant properties are imported into the database after a run, all output files are retained for your reading and copying pleasure.

MESS.DB scales happily to thousands, if not millions, of molecules.

## How to Install
```bash
git clone path/to/messdb.git
cd messdb
python setup.py
```

## Usage Examples
### import a set of molecules
```bash
mess import sources/fda
```
This imports the "FDA-approved drugs" data set into mess.db and the molecules dir.

### apply a method to all molecules in the database
```bash
mess select 'select * from molecule' | mess calculate -m balloon141
```
Balloon generates 3D structures from smiles.

### apply a method that uses previous series of methods (parent path) output
```bash
mess select 'select * from molecule' | mess calculate -m pm7_mopac2012 -pp 2
```
Run a semiempirical calculation using the output from path 2 (the balloon 3D structures in this case, if you've been following along)

## Current Features
-import from most common molecule formats (smi, inchi, xyz, sdf, etc.)
-rationial file structure with graceful duplicate handling
-relational database of all molecules, sources, methods, and properties
-source tracking
-select molecules based on sql queries of their properties
-apply calculations (methods) to any selection

## Planned Features
-report generation
-self-integrity checking
-handling of multiple molecular states (e.g cation, anion, triplet, conformers, etc.)
-database backup/restore
-database pruning

## Contributors
Victor Amin, 2013-
