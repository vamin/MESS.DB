#!/bin/bash

cd ~
echo "module load python; module load openbabel; module load mopac; module load mess; mess select -s $1 'select inchikey from molecule;' | mess calculate -m pm7_mopac2012" | qsub -V -l nodes=1:ppn=1 -N 'pm7_mopac2012'
