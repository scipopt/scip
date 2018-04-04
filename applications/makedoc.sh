#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
APPLICATIONS=(Coloring CycleClustering MinIISC PolySCIP Ringpacking Scheduler STP)

for APPLICATION in ${APPLICATIONS[@]}
do
    echo
    echo ===== $APPLICATION =====
    echo
    if (! make -C $DIR/$APPLICATION doc )
    then
	exit $STATUS
    fi
    echo
done
