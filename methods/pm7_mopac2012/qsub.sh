#!/bin/bash

cd ~
echo "module load openbabel; python /projects/victor/messdb/mess select -s $1 'select inchikey from molecule;' | python /projects/victor/messdb/mess calculate -m pm7_mopac2012" | qsub -V -l nodes=1:ppn=1 -N 'pm7_mopac2012'