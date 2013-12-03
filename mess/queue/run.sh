#!/bin/bash

function usage {
    echo "
    This script will submit a MESS.DB job to a PBS queue.
    
    TYPICAL USAGE:
    $0 -t <threads> -m <method> -p <path_id>
    
    OPTIONS:
        FLAGS:
        -h Show this message
        
        INPUTS:
        -m <method> Name of method to calculate [required]
        -p <path_id> Path id of parent path [default=0]
        -t <threads> Number of processors to split job over [default=1]
        "
        exit 0
}

# Set the defaults
PATH_ID=0
THREADS=1

# Override defaults with options from command line if they are set
while getopts "m:p:t:h" OPT
do
    case $OPT in
        h)
            usage
            exit 1
            ;;
        m)
            METHOD=$OPTARG
            ;;
        p)
            PATH_ID=$OPTARG
            ;;
        t)
            THREADS=$OPTARG
            ;;
        \?) # if illegal options
            usage
            exit 1
            ;;
        :) # if option supplied without a value
            usage
            exit 1
            ;;
    esac
done

# Find the directory that this script resides in
PBSDIR=`dirname $0`
MESSDIR=$PBSDIR'/../..'
LOGS=$MESSDIR'/logs'
QSUB=$MESSDIR/calculate.qsub

for i in $(seq 1 $THREADS)
    do
        # Replace template stand-ins with specified variables and pipe to qsub
        sed "s,__METHOD__,$METHOD,g" $PBSDIR/calculate.pbs \
        | sed "s,__PATH_ID__,$PATH_ID,g" \
        | sed "s,__PART__,$i,g" \
        | sed "s,__TOTAL__,$THREADS,g"\
        | sed "s,__PATH_TO_LOGS__,$LOGS,g" >$QSUB
        qsub $QSUB
        sleep 0.5
        rm $QSUB
    done