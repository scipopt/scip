#!/usr/bin/env bash
#
# This scripts generates the dependencies for SCIP
#

LPSS=(cpx spx spx1 xprs msk clp grb qso none)
OPTS=(opt dbg prf opt-gccold)
TPIS=(omp tny none)
EXPRINTS=(none cppad)

for OPT in ${OPTS[@]}
do
   # dependencies of main SCIP source and objscip library
   # with ZIMPL disabled
   make OPT=$OPT ZIMPL=false LPS=none SYM=none scipdepend

   # dependencies of cmain and cppmain
   make OPT=$OPT ZIMPL=false LPS=none LINKER=C   maindepend
   make OPT=$OPT ZIMPL=false LPS=none LINKER=CPP maindepend

   for LPS in ${LPSS[@]}
   do
      # check if the header for the LP solver are available,
      # or we are in the special case "none"
      # in the case "qso", the include directory is called qsinc
      if [ -e lib/include/$LPS"inc" ] || [ "$LPS" == "none" ] || [ "${LPS:0:3}" == "spx" -a -e lib/include/spxinc ] || [ "$LPS" == "qso" -a -e lib/include/qsinc ] || [ "$LPS" == "clp" -a -e lib/static/clp.*.opt ]
      then
         make LPS=$LPS OPT=$OPT lpidepend
      fi
   done

   # dependencies of nlpi libraries
   for EXPRINT in ${EXPRINTS[@]}
   do
      if test "$EXPRINT" == "none" -o "$EXPRINT" == "cppad" -o -e lib/include/$EXPRINT -o -e lib/include/$EXPRINT"inc"
      then
         make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=false nlpidepend

         for libtype in shared static
         do
            for ipoptopt in opt dbg
            do
               if ls lib/$libtype/ipopt.*.$ipoptopt > /dev/null 2>&1;
               then
                  [ $libtype = "shared" ] && shared=true || shared=false
                  make OPT=$OPT LPS=none EXPRINT=$EXPRINT IPOPT=true SHARED=$shared IPOPTOPT=$ipoptopt nlpidepend
                  break 2
               fi
            done
         done
      fi
   done

   for TPI in ${TPIS[@]}
   do
      make OPT=$OPT ZIMPL=false LPS=none TPI=$TPI tpidepend
   done
done
