#!/usr/bin/env bash
#
# This scripts generates the dependencies for the SCIP applications
#

APPLICATIONS=(Coloring CycleClustering MinIISC PolySCIP Ringpacking Scheduler STP)
OPTS=(opt dbg opt-gccold)

for APPLICATION in ${APPLICATIONS[@]}
do
    echo ===== $APPLICATION =====
    echo
    cd $APPLICATION
    for OPT in ${OPTS[@]}
    do
	make OPT=$OPT ZIMPL=false LPS=none depend
    done
    cd ..
    echo
done
