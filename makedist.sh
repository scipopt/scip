#!/bin/sh
VERSION="0.90.11"
NAME="scip-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f release/$NAME.zip
zip release/$NAME.zip $NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* $NAME/lib \
$NAME/make/make.* $NAME/scip.set \
$NAME/check/check.sh $NAME/check/check.awk $NAME/check/check_cplex.sh $NAME/check/check_cplex.awk \
$NAME/check/miplib3.test $NAME/check/miplib3.solu \
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
$NAME/examples/SamplePricer_C/* $NAME/examples/SamplePricer_C/doc/* $NAME/examples/SamplePricer_C/lib \
$NAME/examples/SamplePricer_C/make/make.* \
$NAME/examples/SamplePricer_C/src/depend.* \
$NAME/examples/SamplePricer_C/src/*.c $NAME/examples/SamplePricer_C/src/*.cpp $NAME/examples/SamplePricer_C/src/*.h \
-x ".*" -x "*~" -x "*/CVS/*"
rm -f $NAME
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
