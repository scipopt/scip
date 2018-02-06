#!/usr/bin/env bash
# 
# This scripts compiles SCIP for all LP solver for which the links exit
#

MAKEFLAG=$1

LPSS=(cpx spx spx2 xprs msk clp grb qso)
OPTS=(opt dbg prf)

# check if zimpl is available
if test -e lib/include/zimplinc/zimpl
then
    ZIMPL=true
else
    ZIMPL=false
fi

for OPT in ${OPTS[@]}
do
    make OPT=$OPT ZIMPL=$ZIMPL LPS=none  $@

    for LPS in ${LPSS[@]}
    do
        # check if the header for the LP solver are available
        if test -e lib/$LPS"inc"
        then
            make LPS=$LPS OPT=$OPT ZIMPL=$ZIMPL $MAKEFLAG $@
        fi
    done
done
