#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

APPLICATIONS=(Coloring Scheduler STP)

for APPLICATION in ${APPLICATIONS[@]}
do
    echo
    echo ===== $APPLICATION =====
    echo
    if (! make -C $APPLICATION doc )
    then
	exit $STATUS
    fi
    echo
done
