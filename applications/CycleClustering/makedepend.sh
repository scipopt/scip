#!/usr/bin/env bash
#
# This scripts generates the dependences for the CycleClustering example
#

OPTS=(opt dbg opt-gccold)

for OPT in ${OPTS[@]}
do
    make OPT=$OPT depend
done
