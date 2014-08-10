#!/bin/bash

module load python
module load sqlite
module load openbabel

THREADS=12
HOST=tim.selfip.org
DIR="$( cd "$( dirname "$0" )" && pwd )"
PASS=$(python ${DIR}/../../mess select | python ${DIR}/../../mess calculate -m balloon --mapreduce)

for i in $(seq 1 $THREADS); do
    python ${DIR}/../mincemeatpy/mincemeat.py -p $PASS $HOST
done