#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="3.0.0"
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

echo generating default setting files
make LPS=none OPT=opt-gccold READLINE=false ZLIB=false ZIMPL=false -j4
bin/scip -c "set default set save doc/inc/parameters.set quit"

# Before we create a tarball change the director and file rights in a command way
echo adjust file modes
find ./ -type d -exec chmod 750 {} \;
find ./ -type f -exec chmod 640 {} \;
find ./ -name "*.sh" -exec chmod 750 {} \;
chmod 750 bin/* interfaces/ampl/get.ASL

tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile \
$NAME/doc/scip* $NAME/doc/xternal.c $NAME/doc/inc/faq.inc \
$NAME/doc/inc/faqcss.inc $NAME/doc/inc/authors.inc $NAME/doc/inc/parameters.set \
$NAME/doc/pictures/miniscippy.png $NAME/doc/pictures/scippy.png \
$NAME/make/make.* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/check_blis.sh $NAME/check/evalcheck_blis.sh $NAME/check/check_blis.awk \
$NAME/check/check_cbc.sh $NAME/check/evalcheck_cbc.sh $NAME/check/check_cbc.awk \
$NAME/check/check_cplex.sh $NAME/check/evalcheck_cplex.sh $NAME/check/check_cplex.awk \
$NAME/check/check_gamscluster.sh $NAME/check/rungamscluster.sh $NAME/check/finishgamscluster.sh \
$NAME/check/evalcheck_gamscluster.sh $NAME/check/check_gams.awk $NAME/check/schulz.sh \
$NAME/check/check_glpk.sh $NAME/check/evalcheck_glpk.sh $NAME/check/check_glpk.awk \
$NAME/check/check_gurobi.sh $NAME/check/evalcheck_gurobi.sh $NAME/check/check_gurobi.awk \
$NAME/check/check_mosek.sh $NAME/check/evalcheck_mosek.sh $NAME/check/check_mosek.awk \
$NAME/check/check_symphony.sh $NAME/check/evalcheck_symphony.sh $NAME/check/check_symphony.awk \
$NAME/check/check_count.sh $NAME/check/evalcheck_count.sh $NAME/check/check_count.awk \
$NAME/check/testset/short.test $NAME/check/testset/short.solu \
$NAME/check/cmpres.awk $NAME/check/allcmpres.sh \
$NAME/check/getlastprob.awk \
$NAME/release-notes/SCIP-* \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h \
$NAME/src/scip/*.c $NAME/src/scip/*.cpp $NAME/src/scip/*.h \
$NAME/src/nlpi/*.c $NAME/src/nlpi/*.cpp $NAME/src/nlpi/*.h \
$NAME/src/xml/*.c $NAME/src/xml/*.h \
$NAME/src/dijkstra/*.c $NAME/src/dijkstra/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.cpp $NAME/src/blockmemshell/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.cpp $NAME/src/tclique/*.h \
$NAME/src/objscip/*.c $NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
$NAME/src/cppad/* $NAME/src/cppad/local/* \
$NAME/examples/Binpacking/Makefile $NAME/examples/Binpacking/INSTALL \
$NAME/examples/Binpacking/doc/* $NAME/examples/Binpacking/doc/pics/binpacking.png \
$NAME/examples/Binpacking/check/testset/short.test $NAME/examples/Binpacking/check/testset/short.solu \
$NAME/examples/Binpacking/src/depend.* \
$NAME/examples/Binpacking/src/*.c $NAME/examples/Binpacking/src/*.h \
$NAME/examples/Binpacking/data/*.bpa \
$NAME/examples/Coloring/* $NAME/examples/Coloring/doc/* $NAME/examples/Coloring/data/* \
$NAME/examples/Coloring/check/testset/short.test $NAME/examples/Coloring/check/testset/short.solu \
$NAME/examples/Coloring/src/depend.* \
$NAME/examples/Coloring/src/*.c $NAME/examples/Coloring/src/*.h \
$NAME/examples/Eventhdlr/* $NAME/examples/Eventhdlr/doc/* \
$NAME/examples/Eventhdlr/src/depend.* \
$NAME/examples/Eventhdlr/src/*.c $NAME/examples/Eventhdlr/src/*.h \
$NAME/examples/LOP/* $NAME/examples/LOP/doc/* $NAME/examples/LOP/data/* \
$NAME/examples/LOP/check/check.sh $NAME/examples/LOP/check/testset/short.test $NAME/examples/LOP/check/testset/short.solu \
$NAME/examples/LOP/src/depend.* \
$NAME/examples/LOP/src/*.c $NAME/examples/LOP/src/*.h \
$NAME/examples/MIPSolver/Makefile  $NAME/examples/MIPSolver/INSTALL $NAME/examples/MIPSolver/scipmip.set \
$NAME/examples/MIPSolver/doc/scipmip.dxy $NAME/examples/MIPSolver/doc/xternal.c \
$NAME/examples/MIPSolver/src/depend.* \
$NAME/examples/MIPSolver/src/*.c $NAME/examples/MIPSolver/src/*.cpp $NAME/examples/MIPSolver/src/*.h \
$NAME/examples/Queens/* $NAME/examples/Queens/doc/scip_intro.tex \
$NAME/examples/Queens/src/depend.* \
$NAME/examples/Queens/src/*.c $NAME/examples/Queens/src/*.cpp \
$NAME/examples/Queens/src/*.h $NAME/examples/Queens/src/*.hpp \
$NAME/examples/TSP/Makefile $NAME/examples/TSP/INSTALL \
$NAME/examples/TSP/runme.sh $NAME/examples/TSP/runviewer.sh \
$NAME/examples/TSP/sciptsp.set \
$NAME/examples/TSP/doc/* \
$NAME/examples/TSP/check/testset/short.test $NAME/examples/TSP/check/testset/short.solu \
$NAME/examples/TSP/src/depend.* \
$NAME/examples/TSP/src/*.c $NAME/examples/TSP/src/*.cpp $NAME/examples/TSP/src/*.h \
$NAME/examples/TSP/tspviewer/*.java $NAME/examples/TSP/tspdata/*.tsp \
$NAME/examples/VRP/Makefile  $NAME/examples/VRP/INSTALL  \
$NAME/examples/VRP/doc/* $NAME/examples/VRP/data/* \
$NAME/examples/VRP/src/depend.* \
$NAME/examples/VRP/src/*.c $NAME/examples/VRP/src/*.cpp $NAME/examples/VRP/src/*.h \
$NAME/interfaces/matlab/* \
$NAME/interfaces/ampl/Makefile $NAME/interfaces/ampl/INSTALL $NAME/interfaces/ampl/get.ASL \
$NAME/interfaces/ampl/src/*.{c,h} $NAME/interfaces/ampl/src/depend.* \
$NAME/interfaces/ampl/check/check.sh $NAME/interfaces/ampl/check/testset/short.test \
$NAME/interfaces/ampl/check/instances/MINLP/*.{col,row,nl} \
$NAME/check/instances/CP/*.cip \
$NAME/check/instances/Indicator/*.lp \
$NAME/check/instances/MIP/*.fzn \
$NAME/check/instances/MIP/*.lp \
$NAME/check/instances/MIP/*.mps \
$NAME/check/instances/MIP/*.osil \
$NAME/check/instances/Orbitope/*.cip \
$NAME/check/instances/MINLP/*.cip \
$NAME/check/instances/MINLP/*.mps \
$NAME/check/instances/MINLP/*.osil \
$NAME/check/instances/MINLP/*.pip \
$NAME/check/instances/PseudoBoolean/*.opb \
$NAME/check/instances/PseudoBoolean/*.wbo \
$NAME/check/instances/SAT/*.cnf \
$NAME/check/instances/SOS/*.lp \
$NAME/check/instances/Semicontinuous/*.lp \
$NAME/check/instances/Semicontinuous/*.mps
rm -f $NAME
echo ""
echo "check version numbers in src/scip/def.h, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/scip/def.h
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
tail src/scip/githash.c
