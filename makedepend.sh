#!/bin/sh
# 
# This scripts generates the dependences for SCIP 
#

LPSS=(cpx spx spx132 xprs msk clp grb qso none)
OPTS=(opt dbg prf)
EXPRINTS=(none cppad)

for OPT in ${OPTS[@]}
do
    # dependencies of main SCIP source and objscip library
    # with ZIMPL disabled
    make OPT=$OPT ZIMPL=false LPS=none scipdepend

    # dependencies of cmain and cppmain
    make OPT=$OPT ZIMPL=false LPS=none LINKER=C   maindepend
    make OPT=$OPT ZIMPL=false LPS=none LINKER=CPP maindepend

    for LPS in ${LPSS[@]}
    do
        # check if the header for the LP solver are available,
        # or we are in the special case "none"
        if test -e lib/$LPS"inc" -o "$LPS" == "none"
        then
             make LPS=$LPS OPT=$OPT lpidepend
        fi

    done

    # dependencies of nlpi libraries
    for EXPRINT in ${EXPRINTS[@]}
    do
        if test "$EXPRINT" == "none" -o -e lib/$EXPRINT -o -e lib/$EXPRINT"inc"
        then
            make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=false nlpidepend

            if test -e lib/ipoptinc
            then
                 make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=true nlpidepend
            fi
        fi
    done

done
