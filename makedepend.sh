#!/bin/sh
# 
# This scripts generates the dependences for SCIP 
#

LPSS=(cpx spx spx132 xprs msk clp grb qso none)
OPTS=(opt dbg prf)


for OPT in ${OPTS[@]}
do
    # dependencies of main SCIP source and objscip library
    # with ZIMPL disabled
    make OPT=$OPT ZIMPL=false LPS=none scipdepend

    # dependencies of cmain and cppmain
    make OPT=$OPT ZIMPL=false LPS=none LINKER=C   maindepend
    make OPT=$OPT ZIMPL=false LPS=none LINKER=CPP maindepend

    # dependencies of nlpi library (so far only for default config)
    make OPT=$OPT LPS=none nlpidepend
    
    for LPS in ${LPSS[@]}
    do
	# check if the header for the LP solver are available 
	if test -e lib/$LPS"inc"
	then
	    make LPS=$LPS OPT=$OPT lpidepend
	fi

	if test "$LPS" == "none"
	then
	    make LPS=$LPS OPT=$OPT lpidepend
	fi
    done
done
