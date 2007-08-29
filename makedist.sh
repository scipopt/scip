#!/bin/sh

# For release versions, only use VERSION="x.xx".
# For development versions, use VERSION="x.xx.xx" with subversion number.
VERSION="1.00.1"
NAME="scip-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f release/$NAME.tgz
tar --no-recursion -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* $NAME/lib \
$NAME/make/make.* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/check_cplex.sh $NAME/check/evalcheck_cplex.sh $NAME/check/check_cplex.awk \
$NAME/check/check_cbc.sh $NAME/check/evalcheck_cbc.sh $NAME/check/check_cbc.awk \
$NAME/check/shortmiplib.test $NAME/check/shortmiplib.solu \
$NAME/check/cmpres.awk $NAME/check/allcmpres.sh \
$NAME/settings/cuts/*.set $NAME/settings/emphasis/*.set $NAME/settings/heuristics/*.set $NAME/settings/presolving/*.set \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h \
$NAME/src/scip/*.c $NAME/src/scip/*.cpp $NAME/src/scip/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.cpp $NAME/src/blockmemshell/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.cpp $NAME/src/tclique/*.h \
$NAME/src/objscip/*.c $NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
$NAME/examples/MIPSolver/* $NAME/examples/MIPSolver/doc/* $NAME/examples/MIPSolver/lib \
$NAME/examples/MIPSolver/make/make.* \
$NAME/examples/MIPSolver/src/depend.* \
$NAME/examples/MIPSolver/src/*.c $NAME/examples/MIPSolver/src/*.cpp $NAME/examples/MIPSolver/src/*.h \
$NAME/examples/TSP/* $NAME/examples/TSP/doc/* $NAME/examples/TSP/lib \
$NAME/examples/TSP/make/make.* \
$NAME/examples/TSP/src/depend.* \
$NAME/examples/TSP/src/*.c $NAME/examples/TSP/src/*.cpp $NAME/examples/TSP/src/*.h \
$NAME/examples/TSP/tspviewer/*.java $NAME/examples/TSP/tspdata/*.tsp \
$NAME/examples/SamplePricer/* $NAME/examples/SamplePricer/doc/* $NAME/examples/SamplePricer/lib \
$NAME/examples/SamplePricer/make/make.* \
$NAME/examples/SamplePricer/src/depend.* \
$NAME/examples/SamplePricer/src/*.c $NAME/examples/SamplePricer/src/*.cpp $NAME/examples/SamplePricer/src/*.h \
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
