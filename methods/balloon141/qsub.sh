#!/bin/bash

cd ~
echo "module load python; module load openbabel; python2.7 /projects/victor/messdb/mess select -s $1 'select inchikey from molecule;' | python2.7 /projects/victor/messdb/mess calculate -m mmff94s_gen3d_ob232" | qsub -V -l walltime=4:00:00 -l nodes=1:ppn=1 -N 'mmff94s_gen3d_ob232'