#!/bin/bash -e

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="7.0.2.4"
NAME="scip-$VERSION"
if test ! -e release
then
    mkdir release
fi
rm -f release/$NAME.tgz
rm -f $NAME.tar

echo "store git hash"
GITHASH=`git describe --always --dirty  | sed 's/^.*-g//'`
echo "#define SCIP_GITHASH \"$GITHASH\"" > src/scip/githash.c

# Before we create a tarball change the director and file rights in a command way
echo "adjust file modes"
git ls-files | xargs dirname | sort -u | xargs chmod 750
git ls-files | xargs chmod 640
git ls-files "*.sh" "scripts/split_scip/*.py" | xargs chmod 750

chmod 750 scripts/* interfaces/ampl/get.ASL check/cmpres.awk cmake/Modules/asan-wrapper applications/PolySCIP/mult_zimpl/mult_zimpl_to_mop.py

# pack files tracked by git and append $NAME to the front
echo "pack files"
git ls-files -c | xargs tar --transform "s|^|${NAME}/|" -chf $NAME.tar \
--exclude="*~" \
--exclude=".*" \
--exclude="scripts/*" \
--exclude="license/*" \
--exclude="lint/*" \
--exclude="make/local/*" \
--exclude="makedist.sh" \
--exclude="makepictures.sh" \
--exclude="suppressions.*" \
--exclude="tex/*" \
--exclude="interfaces/matlab" \
--exclude="check/check_*" \
--exclude="check/*cluster*" \
--exclude="check/checkcuts.py" \
--exclude="check/evalcheck_*" \
--exclude="check/schulz.sh" \
--exclude="check/testset/*"

# append additional files that were excluded before
tar --transform "s|^|${NAME}/|" -rf $NAME.tar \
check/check_count.awk \
check/check_count.sh \
check/check_coverage.sh \
check/evalcheck_count.sh \
check/evalcheck_cluster.sh \
check/testset/short.test \
check/testset/short.solu \
check/testset/coverage.test \
check/testset/stochastic.test \
applications/*/check/testset/* \
examples/*/check/testset/* \
interfaces/*/check/testset/* \
src/scip/githash.c \
scripts/trainEstimation/

# compress the archive
gzip -c $NAME.tar > release/$NAME.tgz

# remove temporary archive
rm -f $NAME.tar

echo ""
echo "check version numbers ($VERSION):"
grep -H "SCIP_VERSION" src/scip/def.h
grep -H "@version" doc/xternal.c
grep -H "^SCIP_VERSION" make/make.project
grep -H "SCIP_VERSION" CMakeLists.txt
echo ""
tail src/scip/githash.c
