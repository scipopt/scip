#!/bin/sh
# 
# This scripts generates the dependences for SCIP 
#

LPSS=(cpx spx spx132 xprs msk clp grb qso none)
OPTS=(opt dbg prf)

# disable ZIMPL for dependency generation

make ZIMPL=false LPS=none depend

for OPT in ${OPTS[@]}
do 
    make OPT=$OPT ZIMPL=false maindepend
    
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