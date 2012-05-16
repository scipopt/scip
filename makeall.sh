#!/usr/bin/env bash
#
# This scripts compiles SCIP for all LP solver for which the links exit
#

MAKEFLAG=$1

LPSS=(cpx spx spx132 xprs msk clp grb qso)
LPSEXS=(qsoex none)
OPTS=(opt dbg prf)
REDUCEDSOLVES=(true false)
BRANCHPLGSS=(true false)
EXACTSOLVES=(false true false)

# check if zimpl is available
if test -e lib/zimplinc/zimpl
then
    ZIMPL=true
else
    ZIMPL=false
fi

make clean

for OPT in ${OPTS[@]}
do
    for REDUCEDSOLVE in ${REDUCEDSOLVES[@]}
    do
	#make OPT=$OPT ZIMPL=$ZIMPL LPS=none REDUCEDSOLVE=$REDUCEDSOLVE $MAKEFLAG $@

	for LPS in ${LPSS[@]}
	do
	    for LPSEX in ${LPSEXS[@]}
	    do
                # check if the header for the LP solver are available,
                # or we are in the special case "none"
                # in the case "qso", the include directory is called qsinc
                # and
                # check if the header for the exact LP solver are available,
                # or we are in the special case "none"
                # in the case "qsoex", the include directory is called qsexinc
		if [ -e lib/$LPS"inc" ] || [ "$LPS" == "none" ] || [ "$LPS" == "qso" -a -e lib/qsinc ] || [ "$LPS" == "clp" -a -e lib/clp.*.opt ]
		then
		    if [ -e lib/$LPSEX"inc" ] || [ "$LPSEX" == "none" ] || [ "$LPSEX" == "qsoex" -a -e lib/qsexinc ]
		    then
			if [ "$REDUCEDSOLVE" == "true" ]
			then
			    for EXACTSOLVE in ${EXACTSOLVES[@]}
			    do
				for BRANCHPLGS in ${BRANCHPLGSS[@]}
				do
				    echo
				    echo $LPS $LPSEX $OPT REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=$EXACTSOLVE BRANCHPLGS=$BRANCHPLGS
				    echo --------------------------
				    make LPS=$LPS LPSEX=$LPSEX OPT=$OPT ZIMPL=$ZIMPL REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=$EXACTSOLVE BRANCHPLGS=$BRANCHPLGS $MAKEFLAG $@
				    #echo $LPS $LPSEX $OPT REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=$EXACTSOLVE BRANCHPLGS=$BRANCHPLGS
				    #bin/scip
				done
			    done
			else
			    echo
			    echo $LPS $LPSEX $OPT REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=false BRANCHPLGS=false
			    echo --------------------------
			    make LPS=$LPS LPSEX=$LPSEX OPT=$OPT ZIMPL=$ZIMPL REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=false BRANCHPLGS=false $MAKEFLAG $@
			    #echo $LPS $LPSEX $OPT REDUCEDSOLVE=$REDUCEDSOLVE EXACTSOLVE=false BRANCHPLGS=false
			    #bin/scip
			fi
		    fi
		fi
            done
	done
    done
done
