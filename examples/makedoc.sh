#!/usr/bin/env bash

# stop on error
set -e

EXAMPLES=(Binpacking CallableLibrary CAP Eventhdlr GMI LOP MIPSolver Queens TSP VRP)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo ===== $EXAMPLE =====
    echo
    make -C $EXAMPLE doc
    echo
done
