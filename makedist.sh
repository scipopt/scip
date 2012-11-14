#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="3.0.1"
NAME="scip-$VERSION"
rm -f $NAME
ln -s . $NAME
if test ! -e release
then
    mkdir release
fi
rm -f release/$NAME.tgz

# run git status to clean the dirty git hash
git status

#echo generating default setting files
#make LPS=none OPT=opt-gccold READLINE=false ZLIB=false ZIMPL=false -j4
#bin/scip -c "set default set save doc/inc/parameters.set quit"

# Before we create a tarball change the director and file rights in a command way
echo adjust file modes
find ./ -type d -exec chmod 750 {} \;
find ./ -type f -exec chmod 640 {} \;
find ./ -name "*.sh" -exec chmod 750 {} \;
chmod 750 bin/* scripts/* interfaces/ampl/get.ASL

tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
--exclude="*xyz*" \
--exclude="nlpioracle.c" \
--exclude="nlpioracle.h" \
--exclude="exprinterpret_cppad.cpp " \
--exclude="nlpi_ipopt.cpp" \
--exclude="lpi_c*" \
--exclude="lpi_grb.c" \
--exclude="lpi_msk.c" \
--exclude="lpi_none.c" \
--exclude="lpi_qso.c" \
--exclude="lpi_spx121.cpp" \
--exclude="lpi_spx132.cpp" \
--exclude="lpi_xprs.c" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile \
$NAME/doc/scip* $NAME/doc/xternal.c $NAME/doc/inc/faq.inc \
$NAME/doc/inc/faqcss.inc $NAME/doc/inc/authors.inc $NAME/doc/inc/parameters.set \
$NAME/doc/pictures/miniscippy.png $NAME/doc/pictures/scippy.png \
$NAME/make/make.* \
$NAME/release-notes/SCIP-* \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp \
$NAME/src/scip/*.c $NAME/src/scip/*.cpp $NAME/src/scip/*.h \
$NAME/src/nlpi/*.c $NAME/src/nlpi/*.cpp $NAME/src/nlpi/*.h \
$NAME/src/xml/*.c $NAME/src/xml/*.h \
$NAME/src/dijkstra/*.c $NAME/src/dijkstra/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.h \
$NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
rm -f $NAME
echo ""
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
tail src/scip/githash.c
