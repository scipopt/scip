#!/usr/bin/env bash
#
# This scripts generates the dependencies for the CallableLibrary example
#

OPTS=(opt dbg opt-gccold)

for OPT in ${OPTS[@]}
do
    make OPT=$OPT depend
done
