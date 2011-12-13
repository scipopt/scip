#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="1.2.0.8"
NAME="scip-$VERSION"
rm -f $NAME
ln -s . $NAME
if test ! -e release
then
    mkdir release
fi
rm -f release/$NAME.tgz
tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile \
$NAME/doc/* $NAME/doc/inc/faq.inc $NAME/doc/inc/faqcss.inc $NAME/doc/inc/authors.inc \
$NAME/make/make.* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/check_blis.sh $NAME/check/evalcheck_blis.sh $NAME/check/check_blis.awk \
$NAME/check/check_cbc.sh $NAME/check/evalcheck_cbc.sh $NAME/check/check_cbc.awk \
$NAME/check/check_cplex.sh $NAME/check/evalcheck_cplex.sh $NAME/check/check_cplex.awk \
$NAME/check/check_glpk.sh $NAME/check/evalcheck_glpk.sh $NAME/check/check_glpk.awk \
$NAME/check/check_gurobi.sh $NAME/check/evalcheck_gurobi.sh $NAME/check/check_gurobi.awk \
$NAME/check/check_mosek.sh $NAME/check/evalcheck_mosek.sh $NAME/check/check_mosek.awk \
$NAME/check/check_symphony.sh $NAME/check/evalcheck_symphony.sh $NAME/check/check_symphony.awk \
$NAME/check/checkcount.sh $NAME/check/evalcheckcount.sh $NAME/check/checkcount.awk \
$NAME/check/shortmiplib.test $NAME/check/shortmiplib.solu \
$NAME/check/cmpres.awk $NAME/check/allcmpres.sh \
$NAME/check/getlastprob.awk \
$NAME/mps2zpl.sh $NAME/mps2zpl.awk \
$NAME/settings/cuts/*.set $NAME/settings/emphasis/*.set $NAME/settings/heuristics/*.set $NAME/settings/presolving/*.set \
$NAME/release-notes/SCIP-* \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h \
$NAME/src/scip/*.c $NAME/src/scip/*.cpp $NAME/src/scip/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.cpp $NAME/src/blockmemshell/*.h \
$NAME/src/rectlu/*.c $NAME/src/rectlu/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.cpp $NAME/src/tclique/*.h \
$NAME/src/objscip/*.c $NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
$NAME/examples/Coloring/* $NAME/examples/Coloring/doc/* $NAME/examples/Coloring/data/* \
$NAME/examples/Coloring/src/depend.* \
$NAME/examples/Coloring/src/*.c $NAME/examples/Coloring/src/*.h \
$NAME/examples/LOP/* $NAME/examples/LOP/doc/* $NAME/examples/LOP/data/* \
$NAME/examples/LOP/src/depend.* \
$NAME/examples/LOP/src/*.c $NAME/examples/LOP/src/*.h \
$NAME/examples/MIPSolver/Makefile  $NAME/examples/MIPSolver/INSTALL $NAME/examples/MIPSolver/scipmip.set \
$NAME/examples/MIPSolver/doc/scipmip.dxy $NAME/examples/MIPSolver/doc/xternal.c \
$NAME/examples/MIPSolver/src/depend.* \
$NAME/examples/MIPSolver/src/*.c $NAME/examples/MIPSolver/src/*.cpp $NAME/examples/MIPSolver/src/*.h \
$NAME/examples/Queens/* \
$NAME/examples/Queens/src/depend.* \
$NAME/examples/Queens/src/*.c $NAME/examples/Queens/src/*.cpp \
$NAME/examples/Queens/src/*.h $NAME/examples/Queens/src/*.hpp \
$NAME/examples/SamplePricer/Makefile $NAME/examples/SamplePricer/INSTALL \
$NAME/examples/SamplePricer/doc/* \
$NAME/examples/SamplePricer/src/depend.* \
$NAME/examples/SamplePricer/src/*.c $NAME/examples/SamplePricer/src/*.cpp $NAME/examples/SamplePricer/src/*.h \
$NAME/examples/TSP/Makefile $NAME/examples/TSP/INSTALL \
$NAME/examples/TSP/runme.sh $NAME/examples/TSP/runviewer.sh \
$NAME/examples/TSP/sciptsp.set \
$NAME/examples/TSP/doc/* \
$NAME/examples/TSP/src/depend.* \
$NAME/examples/TSP/src/*.c $NAME/examples/TSP/src/*.cpp $NAME/examples/TSP/src/*.h \
$NAME/examples/TSP/tspviewer/*.java $NAME/examples/TSP/tspdata/*.tsp \
$NAME/examples/VRP/Makefile  $NAME/examples/VRP/INSTALL  \
$NAME/examples/VRP/doc/* $NAME/examples/VRP/data/* \
$NAME/examples/VRP/src/depend.* \
$NAME/examples/VRP/src/*.c $NAME/examples/VRP/src/*.cpp $NAME/examples/VRP/src/*.h \
$NAME/check/IP/miplib/bell3a.mps.gz \
$NAME/check/IP/miplib/bell5.mps.gz \
$NAME/check/IP/miplib/blend2.mps.gz \
$NAME/check/IP/miplib/dcmulti.mps.gz \
$NAME/check/IP/miplib/egout.mps.gz \
$NAME/check/IP/miplib/enigma.mps.gz \
$NAME/check/IP/miplib/flugpl.mps.gz \
$NAME/check/IP/miplib/gt2.mps.gz \
$NAME/check/IP/miplib/lseu.mps.gz \
$NAME/check/IP/miplib/misc03.mps.gz \
$NAME/check/IP/miplib/mod008.mps.gz \
$NAME/check/IP/miplib/modglob.mps.gz \
$NAME/check/IP/miplib/p0033.mps.gz \
$NAME/check/IP/miplib/p0201.mps.gz \
$NAME/check/IP/miplib/p0282.mps.gz \
$NAME/check/IP/miplib/p0548.mps.gz \
$NAME/check/IP/miplib/pp08a.mps.gz \
$NAME/check/IP/miplib/pp08aCUTS.mps.gz \
$NAME/check/IP/miplib/rgn.mps.gz \
$NAME/check/IP/miplib/stein27.mps.gz \
$NAME/check/IP/miplib/stein45.mps.gz \
$NAME/check/IP/miplib/vpm1.mps.gz \
$NAME/check/IP/miplib/vpm2.mps.gz
rm -f $NAME
echo ""
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
