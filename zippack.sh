zip scip060 INSTALL Makefile doc/* lib lint/* lint/cpp/* lint/posix/* lint/posix/sys/* make/* scip.set src/* -x ".*" -x "*~" -x "*/CVS/*"
echo "check version numbers in src/def.h, doc/xternal.c and zippack.sh:"
grep "SCIP_VERSION" src/def.h
grep "@version" doc/xternal.c
