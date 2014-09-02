# Molecular Electronic Structure Screening DB (MESS.DB) #

This project provides a framework for organizing electronic structure (and
other) calculations on organic molecules into a rational file structure and
relational database. MESS.DB scales happily to millions of molecules.

## Quick Start ##

- Clone the repository:
  ```bash
  git clone git@github.com:vamin/MESS.DB.git
  cd messdb  
  ```

- Alias mess executable:
  ```bash
  alias mess=${PWD}/bin/mess
  ```

- Get help:
  ```bash
  mess -h
  mess [tool] -h
  ```

- Import the sample source molecules:
  ```bash
  mess import fda
  ```

- Prepare your own source for import:
  ```bash
  sources/setup_source.sh
  ```

- Generate 3D structures for all molecules in the database:
  ```bash
  mess select | mess calculate balloon
  ```

- Select only molecules in a particular molecular weight range:
  ```bash
  mess select -hd -n MW -o '>' -v 250
  ```

- Find matches to a [SMARTS][] string:
  ```bash
  mess select | mess match -hd -m [CX3]=[OX1]
  ```

[SMARTS]: http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

## Requirements ##

- Python 2.7+
- SQLite 3.7+
- [Open Babel 2.3+](http://openbabel.org/wiki/Main_Page) 
  with Python bindings

In addition, `mess calculate` provides access to various calculation methods,
each of which may have their own dependencies. You can learn about module
dependencies by running e.g. `mess calculate balloon -h`.

## Storage Structures ##

Every molecule is represented with its own directory and a collection of
records in a relational database.

### Molecule Directories ###

Each molecule is identified by its [inchikey][]. For example, morphine's
InChIKey is BQJCRHHNABKAKU-KBQPJGBKSA-N. Molecule directories are based on the
InChIKey, so all files related to morphine will be stored in
molecules/B/QJ/CRHHNABKAKU-KBQPJGBKSA-N/. In addition to any in/out files
generated during a calculation, each molecule directory contains:

- INCHIKEY.inchi -- the molecule in InChI format  
- INCHIKEY.log -- a log tracking the molecule's calculation history  
- INCHIKEY.notes -- a space for manual annotations of the molecule  
- INCHIKEY.png -- a 2D representation of the molecule  
- sources.tsv -- a table of sources for the molecule, including where to buy if
  the source is commercial

[inchikey]: http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#2.7

### mess.db ###

The relational database is stored in db/mess.db, an SQLite database which is
initialized on the first run of `mess`. The schema is described in db/schema.
It is possible to query the database directly via `mess select`, which also
provides command line interface for common queries.

## Workflow Examples ##

For each example, we will assume that a corpus of molecules to investigate has
already been imported with `mess import [source]`.

### Find low molecular weight molecules with small ionization potentials. ###

- Generate 3D structures with [Balloon][]:
  ```bash
  mess select | mess calculate balloon
  ```
  Running `mess select` with no options outputs a list of every molecule in the
  database.

- Calculate electronic structure with [Mopac][]:
  ```bash
  mess select | mess calculate mopac -p 2
  ```
  The `-p 2` specifies that the calculation should use the geometry in path id
  2 (the Balloon-generated geometry) as input.

- Select candidate molecules:
  ```bash
  mess select -n MW -o '<' -v 250 | mess select -n 'IONIZATION POTENTIAL' -o '<' -v '7'
  ```
  The resulting table can be sorted by piping to common unix tools like `sort`.

### Find molecules with conformations similar to aspirin that also contain a halogen. ###

- Generate 3D structures with [Balloon][]:
  ```bash
  mess select | mess calculate balloon
  ```

- Compare halogen-containing molecules to aspirin target geometry by
  [Spectrophore][]:
  ```bash
  mess match -m [F,Cl,Br,I] | mess match -t aspirin.xyz -s -p 2
  ```

[Balloon]: http://users.abo.fi/mivainio/balloon/
[Mopac]: http://openmopac.net/MOPAC2012.html
[Spectrophore]: http://openbabel.org/docs/dev/Fingerprints/spectrophore.html

## List of Tools ##

For more information about these tools, run `mess [tool] -h`.

- annotate
- backup
- calculate
- check
- import
- inspect
- match
- remove
- select
- transform

There is also a helper script for setting up sources for import,
sources/setup_source.sh.

## Known Bugs ##

Many. This software is currently in alpha, meaning not every feature has been
implemented and those that have been may behave unexpectedly. When the project
has progressed to the point that it is safe to use, this section will be
updatedd with a specific bug list.

## Copyright and License ##
Code and documentation copyright [Victor Amin][] 2013-2014 and made available
under the [AGPL][] license.

[Victor Amin]: http://vamin.net
[AGPL]: https://www.gnu.org/licenses/agpl-3.0.html
