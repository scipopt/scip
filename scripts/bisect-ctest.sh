#! /bin/bash
#
# Script to bisect the git history for the first fail of a ctest.
#
# Usage:
#
# 1. start git bisect (from the SCIP root directory), e.g. via
#
#      git bisect master v600
#
# 2. execute (from the SCIP root directory)
#
#      git bisect scripts/bisect-ctest.sh <BUILDDIR> <TEST>


BUILDDIR=$1
TEST=$2

make -C ${BUILDDIR} scip || exit 125

cd ${BUILDDIR} ; ctest -R ${TEST} -V --timeout 30

exit $?
