#!/usr/bin/env bash
#
# This scripts generates the dependences for SCIP
#

LPSS=(cpx spx spx132 xprs msk clp grb qso none)
LPSEXS=(qsoex none)
OPTS=(opt dbg prf opt-gccold)
EXPRINTS=(none cppad)

for OPT in ${OPTS[@]}
do
    # dependencies of main SCIP source and objscip library
    # with ZIMPL disabled
    # with EXACTSOLVE, REDUCEDSOLVE, BRANCHPLGS disabled
    make OPT=$OPT ZIMPL=false LPS=none LPSEX=none REDUCEDSOLVE=false EXACTSOLVE=false BRANCHPLGS=false scipdepend

    # dependencies of cmain and cppmain
    make OPT=$OPT ZIMPL=false LPS=none LPSEX=none REDUCEDSOLVE=false EXACTSOLVE=false BRANCHPLGS=false LINKER=C   maindepend
    make OPT=$OPT ZIMPL=false LPS=none LPSEX=none REDUCEDSOLVE=false EXACTSOLVE=false BRANCHPLGS=false LINKER=CPP maindepend

    for LPS in ${LPSS[@]}
    do
        # check if the header for the LP solver are available,
        # or we are in the special case "none"
        # in the case "qso", the include directory is called qsinc
        if [ -e lib/$LPS"inc" ] || [ "$LPS" == "none" ] || [ "$LPS" == "qso" -a -e lib/qsinc ] || [ "$LPS" == "clp" -a -e lib/clp.*.opt ]
        then
             make LPS=$LPS OPT=$OPT lpidepend
        fi

    done

    for LPSEX in ${LPSEXS[@]}
    do
        # check if the header for the exact LP solver are available,
        # or we are in the special case "none"
        # in the case "qsoex", the include directory is called qsexinc
        if [ -e lib/$LPSEX"inc" ] || [ "$LPSEX" == "none" ] || [ "$LPSEX" == "qsoex" -a -e lib/qsexinc ]
        then
             make LPSEX=$LPSEX OPT=$OPT lpiexdepend
        fi
    done

    # dependencies of nlpi libraries
    for EXPRINT in ${EXPRINTS[@]}
    do
        if test "$EXPRINT" == "none" -o "$EXPRINT" == "cppad" -o -e lib/$EXPRINT -o -e lib/$EXPRINT"inc"
        then
            make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=false nlpidepend

            if ls lib/ipopt.*.opt > /dev/null 2>&1;
            then
                make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=true IPOPTOPT=opt nlpidepend
            elif ls lib/ipopt.*.dbg > /dev/null 2>&1;
            then
                make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=true IPOPTOPT=dbg nlpidepend
            fi
        fi
    done
done
