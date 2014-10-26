#!/bin/bash
# Copyright 2013-2014 Victor Amin, http://vamin.net/

function usage {
    echo "
    This script will submit a 'mess calculate' job to a PBS or SGE queue.
    
    TYPICAL USAGE:
    mess select | $0 <method> [options]
    
    OPTIONS:
        <method> Name of method to calculate [required]
        -p <path> Parent path
        
        -N <threads> # of threads [default=3; 1 server, 1 mapper, 1 reducer]
        -H Show this message
        "
}

# Set and check the defaults
if [ -t 0 ]; then
    echo "You must provide a list of INCHIKEYs to apply the method to."
    usage
    exit 1
fi

METHOD=$1
PATH_ID=0
THREADS=3

if [[ $METHOD == '' ]] || [[ $METHOD == -* ]]; then
    echo "You must specify a method first."
    usage
    exit 1
else
    shift
fi

# Override defaults with options from command line if they are set
while getopts "p:N:H" OPT; do
    case $OPT in
        H)
            usage
            exit 0
            ;;
        N)
            THREADS=$OPTARG
            ;;
        p)
            PATH_ID=$OPTARG
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

# Check overridden arguments
if [[ $THREADS < 3 ]]; then
    echo "You must ask for at least 3 threads."
    echo "The first two threads are reserved for the server and the reducer."
    echo "Additional threads are assigned to mappers."
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
QUEUEDIR=$(cd $(dirname $0); pwd -P)
MESSDIR=$QUEUEDIR'/../..'
TEMPDIR=$MESSDIR'/temp/'
INPUTS=$TEMPDIR'inputs-'$METHOD'.txt'

# Create temp directory if it does not exist
if [ ! -d $TEMPDIR ]
    then
        mkdir $TEMPDIR
fi

echo "Caching inputs."
cat >$INPUTS

echo "Deploying server for '$METHOD' calculation."
SUB=$TEMPDIR'mess-server-'$METHOD'.qsub'
sed "s,__NAME__,MESS.DB-server-$METHOD,g" $QUEUEDIR/qsub/active_comments.qsub >$SUB
cat $QUEUEDIR/qsub/dependencies.qsub >>$SUB
sed "s,__INPUTS__,$INPUTS,g" $QUEUEDIR/qsub/mess-server.qsub \
	 | sed "s,__MESSDIR__,$MESSDIR,g" \
	 | sed "s,__METHOD__,$METHOD,g" \
	 | sed "s,__PARENT__,$PATH_ID,g" >>$SUB
$SUB_CMD $SUB
sleep 2

echo "Deploying reducer for '$METHOD' calculation."
SUB=$TEMPDIR'mess-reducer-'$METHOD'.qsub'
sed "s,__NAME__,MESS.DB-reducer-$METHOD,g" $QUEUEDIR/qsub/active_comments.qsub >$SUB
cat $QUEUEDIR/qsub/dependencies.qsub >>$SUB
sed "s,__MESSDIR__,$MESSDIR,g" $QUEUEDIR/qsub/mess-reducer.qsub \
	 | sed "s,__METHOD__,$METHOD,g" \
	 | sed "s,__PARENT__,$PATH_ID,g" >>$SUB
$SUB_CMD $SUB
sleep 0.2

echo "Deploying mappers ($(($THREADS-2))) for '$METHOD' calculation."
SUB=$TEMPDIR'mess-mapper-'$METHOD'.qsub'
sed "s,__NAME__,MESS.DB-mapper-$METHOD,g" $QUEUEDIR/qsub/active_comments.qsub >$SUB
cat $QUEUEDIR/qsub/dependencies.qsub >>$SUB
sed "s,__MESSDIR__,$MESSDIR,g" $QUEUEDIR/qsub/mess-mapper.qsub \
	 | sed "s,__METHOD__,$METHOD,g" \
	 | sed "s,__PARENT__,$PATH_ID,g" >>$SUB
for i in $(seq 1 $(($THREADS-2))); do
    $SUB_CMD $SUB
    sleep 0.1
done
