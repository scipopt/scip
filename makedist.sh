#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.xx" with subversion number.
VERSION="1.1.0"
NAME="scip-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f release/$NAME.tgz
tar --no-recursion -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile \
$NAME/doc/* $NAME/doc/inc/faq.inc $NAME/doc/inc/faqcss.inc \
$NAME/lib \
$NAME/make/make.* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/check_cplex.sh $NAME/check/evalcheck_cplex.sh $NAME/check/check_cplex.awk \
$NAME/check/check_cbc.sh $NAME/check/evalcheck_cbc.sh $NAME/check/check_cbc.awk \
$NAME/check/checkcount.sh $NAME/check/evalcheckcount.sh $NAME/check/checkcount.awk \
$NAME/check/shortmiplib.test $NAME/check/shortmiplib.solu \
$NAME/check/cmpres.awk $NAME/check/allcmpres.sh \
$NAME/settings/cuts/*.set $NAME/settings/emphasis/*.set $NAME/settings/heuristics/*.set $NAME/settings/presolving/*.set \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h \
$NAME/src/scip/*.c $NAME/src/scip/*.cpp $NAME/src/scip/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.cpp $NAME/src/blockmemshell/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.cpp $NAME/src/tclique/*.h \
$NAME/src/objscip/*.c $NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
$NAME/examples/Coloring/* $NAME/examples/Coloring/doc/* $NAME/examples/Coloring/data/* \
$NAME/examples/Coloring/src/depend.* \
$NAME/examples/Coloring/src/*.c $NAME/examples/Coloring/src/*.h \
$NAME/examples/LOP/* $NAME/examples/LOP/doc/* $NAME/examples/LOP/data/* \
$NAME/examples/LOP/src/depend.* \
$NAME/examples/LOP/src/*.c $NAME/examples/LOP/src/*.h \
$NAME/examples/MIPSolver/* $NAME/examples/MIPSolver/doc/* \
$NAME/examples/MIPSolver/src/depend.* \
$NAME/examples/MIPSolver/src/*.c $NAME/examples/MIPSolver/src/*.cpp $NAME/examples/MIPSolver/src/*.h \
$NAME/examples/Queens/* \
$NAME/examples/Quenns/src/depend.* \
$NAME/examples/Queens/src/*.c $NAME/examples/Queens/src/*.cpp \
$NAME/examples/Queens/src/*.h $NAME/examples/Queens/src/*.hpp \
$NAME/examples/SamplePricer/* $NAME/examples/SamplePricer/doc/* \
$NAME/examples/SamplePricer/src/depend.* \
$NAME/examples/SamplePricer/src/*.c $NAME/examples/SamplePricer/src/*.cpp $NAME/examples/SamplePricer/src/*.h \
$NAME/examples/SamplePricer_C/* $NAME/examples/SamplePricer_C/doc/* \
$NAME/examples/SamplePricer_C/src/depend.* \
$NAME/examples/SamplePricer_C/src/*.c $NAME/examples/SamplePricer_C/src/*.h \
$NAME/examples/TSP/* $NAME/examples/TSP/doc/* \
$NAME/examples/TSP/src/depend.* \
$NAME/examples/TSP/src/*.c $NAME/examples/TSP/src/*.cpp $NAME/examples/TSP/src/*.h \
$NAME/examples/TSP/tspviewer/*.java $NAME/examples/TSP/tspdata/*.tsp \
$NAME/examples/VRP/* $NAME/examples/VRP/doc/* $NAME/examples/VRP/data/* \
$NAME/examples/VRP/src/depend.* \
$NAME/examples/VRP/src/*.c $NAME/examples/VRP/src/*.cpp $NAME/examples/VRP/src/*.h \
$NAME/check/IP/miplib/bell3a.mps \
$NAME/check/IP/miplib/bell5.mps \
$NAME/check/IP/miplib/blend2.mps \
$NAME/check/IP/miplib/dcmulti.mps \
$NAME/check/IP/miplib/egout.mps \
$NAME/check/IP/miplib/enigma.mps \
$NAME/check/IP/miplib/flugpl.mps \
$NAME/check/IP/miplib/gt2.mps \
$NAME/check/IP/miplib/lseu.mps \
$NAME/check/IP/miplib/misc03.mps \
$NAME/check/IP/miplib/mod008.mps \
$NAME/check/IP/miplib/modglob.mps \
$NAME/check/IP/miplib/p0033.mps \
$NAME/check/IP/miplib/p0201.mps \
$NAME/check/IP/miplib/p0282.mps \
$NAME/check/IP/miplib/p0548.mps \
$NAME/check/IP/miplib/pp08a.mps \
$NAME/check/IP/miplib/pp08aCUTS.mps \
$NAME/check/IP/miplib/rgn.mps \
$NAME/check/IP/miplib/stein27.mps \
$NAME/check/IP/miplib/stein45.mps \
$NAME/check/IP/miplib/vpm1.mps \
$NAME/check/IP/miplib/vpm2.mps
rm -f $NAME
echo ""
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
