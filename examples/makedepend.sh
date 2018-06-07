#!/usr/bin/env bash
#
# This scripts generates the dependences for SCIP
#

# stop on error
set -e

EXAMPLES=(Binpacking CallableLibrary Eventhdlr GMI LOP MIPSolver Queens Relaxator SCFLP TSP VRP)
OPTS=(opt dbg opt-gccold)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo ===== $EXAMPLE =====
    echo
    pushd $EXAMPLE > /dev/null
    for OPT in ${OPTS[@]}
    do
        make OPT=$OPT ZIMPL=false LPS=none depend
    done
    popd > /dev/null
    echo
done
