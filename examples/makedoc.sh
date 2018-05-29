#!/usr/bin/env bash

# stop on error
set -e

EXAMPLES=(Binpacking CallableLibrary Eventhdlr GMI LOP MIPSolver Queens SCFLP TSP VRP)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo ===== $EXAMPLE =====
    echo
    make -C $EXAMPLE doc
    echo
done
