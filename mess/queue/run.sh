#!/bin/bash

function usage {
    echo "
    This script will submit a MESS.DB job to a PBS queue.
    
    TYPICAL USAGE:
    $0 -t <tool> [tool-specific options]
    
    GLOBAL OPTIONS:
        -h Show this message
        -t <tool> import, backup, or calculate [required]

    IMPORT OPTIONS:
        -s <source> Base name of source (i.e., not a path) [required]

    BACKUP OPTIONS:
        -r <backup> Path of a backup to restore from
        -n <threads> # of threads (limited to single node) [default=1]

    CALCULATE OPTIONS:
        -m <method> Name of method to calculate [required]
        -p <path_id> Path id of parent path [default=0]
        -q <query> Query file
        -n <threads> # of threads [default=1]
        "
}

# Set the defaults
PATH_ID=0
THREADS=1

# Override defaults with options from command line if they are set
while getopts "t:s:r:m:p:q:n:h" OPT; do
    case $OPT in
        h)
            usage
            exit 0
            ;;
        t)
            TOOL=$OPTARG
            ;;
        s)
            SOURCE=$OPTARG
            ;;
        r)
            RESTORE='-r '$OPTARG
            ;;
        m)
            METHOD=$OPTARG
            ;;
        p)
            PATH_ID=$OPTARG
            ;;
        q)
            QUERY='-q '`readlink -f $OPTARG`
            ;;
        n)
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

# Check arguments
if [[ $TOOL != 'import' ]] && [[ $TOOL != 'backup' ]] && \
   [[ $TOOL != 'calculate' ]]; then
   echo "Tool must be set to import, backup, or calculate."
   usage
   exit 1
fi

if [[ $TOOL == 'import' ]]; then
    if [[ -z $SOURCE ]]; then
        echo "You must specify a source to run import."
        usage
        exit 1
    fi
fi

if [[ $TOOL == 'calculate' ]]; then
    if [[ -z $METHOD ]]; then
        echo "You must specify a method to run calculate."
        usage
        exit 1
    fi
fi

# Detect msub/qsub
if command -v msub >/dev/null 2>&1; then
    SUB_CMD=msub
elif command -v qsub >/dev/null 2>&1; then
    SUB_CMD=qsub
else
    echo "Neither msub or qsub are installed."
    exit 1
fi

# Find the directory that this script resides in
QUEUEDIR=`dirname $0`
MESSDIR=$QUEUEDIR'/../..'
LOGS=$MESSDIR'/logs'
SUB=$MESSDIR/$TOOL.sub

for i in $(seq 1 $THREADS); do
    # Replace template stand-ins with specified variables and pipe to sub
    sed "s,__PATH_TO_LOGS__,$LOGS,g" $QUEUEDIR/$TOOL.pbs \
    | sed "s,__SOURCE__,$SOURCE,g" \
    | sed "s,__RESTORE__,$RESTORE,g" \
    | sed "s,__METHOD__,$METHOD,g" \
    | sed "s,__PATH_ID__,$PATH_ID,g" \
    | sed "s,__QUERY__,$QUERY,g" \
    | sed "s,__PART__,$i,g" \
    | sed "s,__THREADS__,$THREADS,g" >$SUB
    $SUB_CMD $SUB
    sleep 0.2
    rm $SUB
    if [[ $TOOL == 'backup' ]]; then
        break
    fi 
done