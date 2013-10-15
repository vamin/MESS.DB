#!/bin/bash

cd ~
echo "module load python; module load openbabel; module load balloon; module load mess; mess select -s $1 'select inchikey from molecule;' | mess calculate -m balloon141" | qsub -V -l walltime=2:00:00:00 -l nodes=1:ppn=1 -N 'balloon141'
