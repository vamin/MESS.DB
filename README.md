# MESS.DB - Molecular Electronic Structure (and other Stuff) DB

This project provides a framework for organizing electronic structure (and other) calculations on organic molecules into a relational database.  

## Motivation and Design Philosophy

MESS.DB is designed to facilitate large scale screening. Molecules are represented in the database and file structure by their InChiKeys. When multiple computations (methods) are conducted sequentially, the path used to arrive at a result is tracked in the database. Methods can be applied to a subset of the database, facilitating 

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
### apply a method to all molecules in the database
```bash
mess select 'select * from molecule' | mess calculate -m balloon141
```
### apply a method that uses previous series of methods (parent path) output
```bash
mess select 'select * from molecule' | mess calculate -m dft_b3lyp_6-31g_nwchem -pp 3
```

## Current Features
-import
-select
-apply method
-duplicate handling
-source tracking
-method tracking
-queueing engine friendly

## Planned Features
-report generation
-improved self-validation and integrity checking
-handling of multiple molecular states (e.g cation, anion, triplet, etc.)

## Contributors
Victor Amin, 2013-