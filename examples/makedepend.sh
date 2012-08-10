#!/usr/bin/env bash
#
# This scripts generates the dependences for SCIP
#

EXAMPLES=(Coloring Binpacking Eventhdlr LOP MIPSolver Queens Scheduler Testslack TSP VRP CallableLibrary)
OPTS=(opt dbg opt-gccold)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo ===== $EXAMPLE =====
    echo
    cd $EXAMPLE
    for OPT in ${OPTS[@]}
    do
	make OPT=$OPT ZIMPL=false LPS=none depend
    done
    cd ..
    echo
done
