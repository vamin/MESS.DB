# Molecule Directories #

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
