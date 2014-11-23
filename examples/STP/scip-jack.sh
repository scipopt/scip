#!/usr/bin/env bash

BINARY=bin/stp
SETTINGS=settings/stp3.set

FILENAME=$1
TIME=$2
THREADS=$3
OUTPUTFILE=$4

# check all variables defined
if [ -z ${OUTPUTFILE} ]
then
    echo Skipping SCIP-Jack run since not command line parameters were given
    echo "usage: scip-jack.sh FILENAME TIME THREADS OUTPUTFILE "
    exit 1;
fi


bash -c  "$BINARY -c 'set load $SETTINGS set stp logfile $OUTPUTFILE set limit time $TIME read $FILENAME o write stpsol disp stat q'"
