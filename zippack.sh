#!/bin/sh
VERSION="0.78j"
NAME="scip-$VERSION"
cd ..
ln -s scip $NAME
zip scip/release/$NAME.zip $NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* $NAME/lib $NAME/lint/* $NAME/lint/cpp/* $NAME/lint/posix/* $NAME/lint/posix/sys/* $NAME/make/* $NAME/scip.set $NAME/src/* $NAME/src/scip/* $NAME/src/objscip/* $NAME/examples/* $NAME/examples/MIPSolver/* $NAME/examples/MIPSolver/doc/* $NAME/examples/MIPSolver/lib $NAME/examples/MIPSolver/lint/* $NAME/examples/MIPSolver/lint/cpp/* $NAME/examples/MIPSolver/lint/posix/* $NAME/examples/MIPSolver/lint/posix/sys/* $NAME/examples/MIPSolver/make/* $NAME/examples/MIPSolver/src/* -x ".*" -x "*~" -x "*/CVS/*"
rm $NAME
cd scip
echo "check version numbers in src/def.h, doc/xternal.c and zippack.sh ($VERSION):"
grep "SCIP_VERSION" src/def.h
grep "@version" doc/xternal.c
