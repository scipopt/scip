#!/usr/bin/env bash
#
# This scripts generates the dependencies for the SCIP applications
#

APPLICATIONS=(Scheduler)

for APPLICATION in ${APPLICATIONS[@]}
do
    echo ===== $APPLICATION =====
    echo
    cd $APPLICATION
    make depend
    cd ..
    echo
done
