#!/bin/sh
NAME=scip076a
cd ..
ln -s scip $NAME
zip scip/release/$NAME $NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* $NAME/lib $NAME/lint/* $NAME/lint/cpp/* $NAME/lint/posix/* $NAME/lint/posix/sys/* $NAME/make/* $NAME/scip.set $NAME/src/* -x ".*" -x "*~" -x "*/CVS/*"
rm $NAME
cd scip
echo "check version numbers in src/def.h, doc/xternal.c and zippack.sh ($NAME):"
grep "SCIP_VERSION" src/def.h
grep "@version" doc/xternal.c
