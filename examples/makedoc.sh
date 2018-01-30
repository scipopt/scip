#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

EXAMPLES=(Binpacking CallableLibrary CAP Eventhdlr GMI LOP MIPSolver Queens TSP VRP)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo ===== $EXAMPLE =====
    echo
    if (! make -C $EXAMPLE doc )
    then
	exit $STATUS
    fi
    echo
done
