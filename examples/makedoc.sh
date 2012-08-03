#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

EXAMPLES=(Coloring Binpacking Eventhdlr LOP MIPSolver Queens Scheduler Testslack TSP VRP CallableLibrary)

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
