#!/usr/bin/env bash
#
# This scripts generates the dependences for SCIP
#

APPLICATIONS=(Coloring MinIISC Scheduler STP)
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
