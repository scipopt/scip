#!/bin/sh
VERSION="0.80a"
NAME="scip-$VERSION"
ln -s . $NAME
zip release/$NAME.zip $NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* $NAME/lib \
$NAME/lint/* $NAME/lint/cpp/* $NAME/lint/posix/* $NAME/lint/posix/sys/* $NAME/make/* $NAME/scip.set \
$NAME/src/* $NAME/src/scip/* $NAME/src/blockmemshell/* $NAME/src/tclique/* $NAME/src/objscip/* \
$NAME/examples/* \
$NAME/examples/MIPSolver/* $NAME/examples/MIPSolver/doc/* $NAME/examples/MIPSolver/lib \
$NAME/examples/MIPSolver/lint/* $NAME/examples/MIPSolver/lint/cpp/* $NAME/examples/MIPSolver/lint/posix/* \
$NAME/examples/MIPSolver/lint/posix/sys/* $NAME/examples/MIPSolver/make/* $NAME/examples/MIPSolver/src/* \
$NAME/examples/TSP/* $NAME/examples/TSP/doc/* $NAME/examples/TSP/lib \
$NAME/examples/TSP/lint/* $NAME/examples/TSP/lint/cpp/* $NAME/examples/TSP/lint/posix/* \
$NAME/examples/TSP/lint/posix/sys/* $NAME/examples/TSP/make/* $NAME/examples/TSP/src/* $NAME/examples/TSP/tspviewer/* \
$NAME/examples/SamplePricer/* $NAME/examples/SamplePricer/doc/* $NAME/examples/SamplePricer/lib \
$NAME/examples/SamplePricer/lint/* $NAME/examples/SamplePricer/lint/cpp/* $NAME/examples/SamplePricer/lint/posix/* \
$NAME/examples/SamplePricer/lint/posix/sys/* $NAME/examples/SamplePricer/make/* $NAME/examples/SamplePricer/src/* \
-x ".*" -x "*~" -x "*/CVS/*"
rm $NAME
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and zippack.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "VERSION" Makefile
