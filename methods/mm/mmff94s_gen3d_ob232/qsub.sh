#!/bin/bash

cd ~
echo "module load openbabel; python /projects/victor/messdb/mess select -s $1 'select inchikey from molecule;' | python /projects/victor/messdb/mess calculate -m mmff94s_gen3d_ob232" | qsub -V -l nodes=1:ppn=1 -N 'mmff94s_gen3d_ob232'