#!/usr/bin/env bash
#
# This scripts generates the dependencies for the MinIISC application
#

OPTS=(opt dbg opt-gccold)

for OPT in ${OPTS[@]}
do
    make OPT=$OPT depend
done
