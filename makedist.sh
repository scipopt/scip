#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="4.0.1"
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

# Before we create a tarball change the director and file rights in a command way
echo adjust file modes
find ./ -type d -exec chmod 750 {} \;
find ./ -type f -exec chmod 640 {} \;
find ./ -name "*.sh" -exec chmod 750 {} \;
chmod 750 bin/* scripts/* interfaces/ampl/get.ASL check/cmpres.awk

tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile \
$NAME/doc/scip* $NAME/doc/xternal.c $NAME/doc/builddoc.sh \
$NAME/doc/inc/codestyle/* $NAME/doc/inc/shelltutorial $NAME/doc/inc/debugexamples \
$NAME/doc/inc/faq/faqtext.txt $NAME/doc/inc/faq/*.py $NAME/doc/inc/faq/localfaq.php \
$NAME/doc/inc/makeexamples/* $NAME/doc/inc/shelltutorial/commands $NAME/doc/inc/shelltutorial/insertsnippetstutorial.py \
$NAME/doc/inc/debugexamples/example*.txt \
$NAME/doc/pictures/miniscippy.png $NAME/doc/pictures/scippy.png \
$NAME/make/make.* \
scip/applications/PolySCIP/doc/CMakeLists.txt  \
$NAME/applications/PolySCIP/src/CMakeLists.txt \
$NAME/applications/PolySCIP/CMakeLists.txt     \
$NAME/check/CMakeLists.txt                     \
$NAME/tests/CMakeLists.txt                     \
$NAME/src/CMakeLists.txt                       \
$NAME/CMakeLists.txt \
$NAME/scip-config.cmake.in \
$NAME/examples/TSP/check/CMakeLists.txt        \
$NAME/examples/TSP/CMakeLists.txt              \
$NAME/examples/Binpacking/check/CMakeLists.txt \
$NAME/examples/Binpacking/CMakeLists.txt       \
$NAME/examples/CallableLibrary/CMakeLists.txt  \
$NAME/examples/Relaxator/check/CMakeLists.txt  \
$NAME/examples/Relaxator/CMakeLists.txt        \
$NAME/examples/LOP/check/CMakeLists.txt        \
$NAME/examples/LOP/CMakeLists.txt              \
$NAME/examples/MIPSolver/CMakeLists.txt        \
$NAME/examples/GMI/CMakeLists.txt              \
$NAME/examples/VRP/check/CMakeLists.txt        \
$NAME/examples/VRP/CMakeLists.txt              \
$NAME/examples/Queens/CMakeLists.txt           \
$NAME/examples/CMakeLists.txt                  \
$NAME/examples/Eventhdlr/check/CMakeLists.txt  \
$NAME/examples/Eventhdlr/CMakeLists.txt        \

$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/check_gamscluster.sh $NAME/check/rungamscluster.sh $NAME/check/finishgamscluster.sh \
$NAME/check/evalcheck_gamscluster.sh $NAME/check/check_gams.awk $NAME/check/schulz.sh \
$NAME/check/check_count.sh $NAME/check/evalcheck_count.sh $NAME/check/check_count.awk \
$NAME/check/testset/short.test $NAME/check/testset/short.solu \
$NAME/check/cmpres.awk $NAME/check/allcmpres.sh \
$NAME/check/getlastprob.awk \
$NAME/check/configuration_set.sh $NAME/check/configuration_logfiles.sh \
$NAME/check/configuration_tmpfile_setup_scip.sh \
$NAME/check/run.sh $NAME/check/evalcheck_cluster.sh \
$NAME/release-notes/SCIP-* \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp \
$NAME/src/scip/*.c $NAME/src/scip/*.h \
$NAME/src/nlpi/*.c $NAME/src/nlpi/*.cpp $NAME/src/nlpi/*.h \
$NAME/src/lpi/*.c $NAME/src/lpi/*.cpp $NAME/src/lpi/*.h \
$NAME/src/tpi/*.c $NAME/src/tpi/*.h \
$NAME/src/xml/*.c $NAME/src/xml/*.h \
$NAME/src/tinycthread/*.c $NAME/src/tinycthread/*.h \
$NAME/src/dijkstra/*.c $NAME/src/dijkstra/*.h \
$NAME/src/blockmemshell/*.c $NAME/src/blockmemshell/*.h \
$NAME/src/tclique/*.c $NAME/src/tclique/*.h \
$NAME/src/objscip/*.cpp $NAME/src/objscip/*.h \
$NAME/src/cppad/* $NAME/src/cppad/local/* $NAME/src/cppad/utility/* \
$NAME/applications/Coloring/* $NAME/applications/Coloring/doc/* $NAME/applications/Coloring/data/* \
$NAME/applications/Coloring/check/testset/short.test $NAME/applications/Coloring/check/testset/short.solu \
$NAME/applications/Coloring/src/depend.* \
$NAME/applications/Coloring/src/*.c $NAME/applications/Coloring/src/*.h \
$NAME/applications/Scheduler/doc/* \
$NAME/applications/Scheduler/check/testset/short.test $NAME/applications/Scheduler/check/testset/short.solu \
$NAME/applications/Scheduler/src/depend.* \
$NAME/applications/Scheduler/src/*.c $NAME/applications/Scheduler/src/*.cpp $NAME/applications/Scheduler/src/*.h \
$NAME/applications/Scheduler/data/*.sm \
$NAME/applications/Scheduler/data/*.cmin \
$NAME/applications/Scheduler/Makefile \
$NAME/applications/MinIISC/Makefile $NAME/applications/MinIISC/INSTALL \
$NAME/applications/MinIISC/doc/* \
$NAME/applications/MinIISC/src/* \
$NAME/applications/MinIISC/data/* \
$NAME/applications/MinIISC/check/configuration_tmpfile_setup_miniisc.sh $NAME/applications/MinIISC/check/run.sh \
$NAME/applications/MinIISC/check/testset/short.* \
$NAME/applications/PolySCIP/doc/* \
$NAME/applications/PolySCIP/src/*.cpp $NAME/applications/PolySCIP/src/*.h \
$NAME/applications/PolySCIP/src/tclap/* \
$NAME/applications/PolySCIP/data/*.mop \
$NAME/applications/PolySCIP/data/AP_p-3_n-5.dat \
$NAME/applications/PolySCIP/Makefile \
$NAME/applications/PolySCIP/mult_zimpl/AP_p-3_n-5.zpl \
$NAME/applications/PolySCIP/mult_zimpl/README \
$NAME/applications/PolySCIP/mult_zimpl/tenfelde_podehl.zpl \
$NAME/applications/PolySCIP/mult_zimpl/*.py \
$NAME/applications/PolySCIP/INSTALL $NAME/applications/PolySCIP/LICENCE \
$NAME/applications/PolySCIP/README $NAME/applications/PolySCIP/scipmip.set \
$NAME/applications/STP/doc/* \
$NAME/applications/STP/src/depend.* \
$NAME/applications/STP/src/*.c $NAME/applications/STP/src/*.h \
$NAME/applications/STP/check/testset/*.test $NAME/applications/STP/check/testset/*.solu \
$NAME/applications/STP/data/short/* \
$NAME/applications/STP/Makefile $NAME/applications/STP/INSTALL \
$NAME/examples/xternal_examples.c \
$NAME/examples/Binpacking/Makefile $NAME/examples/Binpacking/INSTALL \
$NAME/examples/Binpacking/doc/* $NAME/examples/Binpacking/doc/pics/binpacking.png \
$NAME/examples/Binpacking/check/testset/short.test $NAME/examples/Binpacking/check/testset/short.solu \
$NAME/examples/Binpacking/src/depend.* \
$NAME/examples/Binpacking/src/*.c $NAME/examples/Binpacking/src/*.h \
$NAME/examples/Binpacking/data/*.bpa \
$NAME/examples/CallableLibrary/Makefile \
$NAME/examples/CallableLibrary/INSTALL \
$NAME/examples/CallableLibrary/doc/* \
$NAME/examples/CallableLibrary/src/depend.* $NAME/examples/CallableLibrary/src/*.c \
$NAME/examples/Eventhdlr/* $NAME/examples/Eventhdlr/doc/* \
$NAME/examples/Eventhdlr/src/depend.* \
$NAME/examples/Eventhdlr/src/*.c $NAME/examples/Eventhdlr/src/*.h \
$NAME/examples/GMI/Makefile \
$NAME/examples/GMI/src/Makefile \
$NAME/examples/GMI/INSTALL \
$NAME/examples/GMI/settings/scipdefault.set \
$NAME/examples/GMI/doc/* \
$NAME/examples/GMI/check/testset/short.* \
$NAME/examples/GMI/settings/gmi* $NAME/examples/GMI/src/depend.* \
$NAME/examples/GMI/src/*.c $NAME/examples/GMI/src/*.h \
$NAME/examples/LOP/* $NAME/examples/LOP/doc/* $NAME/examples/LOP/data/* \
$NAME/examples/LOP/check/testset/short.* \
$NAME/examples/LOP/src/depend.* $NAME/examples/LOP/src/Makefile \
$NAME/examples/LOP/settings/default.set \
$NAME/examples/LOP/src/*.c $NAME/examples/LOP/src/*.h \
$NAME/examples/MIPSolver/Makefile $NAME/examples/MIPSolver/INSTALL $NAME/examples/MIPSolver/scipmip.set \
$NAME/examples/MIPSolver/doc/* \
$NAME/examples/MIPSolver/src/depend.* \
$NAME/examples/MIPSolver/src/*.cpp \
$NAME/examples/Queens/* $NAME/examples/Queens/doc/scip_intro.tex \
$NAME/examples/Queens/src/depend.* \
$NAME/examples/Queens/src/*.cpp $NAME/examples/Queens/src/*.hpp \
$NAME/examples/Relaxator/INSTALL \
$NAME/examples/Relaxator/Makefile \
$NAME/examples/Relaxator/check/testset/short.test \
$NAME/examples/Relaxator/doc/xternal_relaxator.c \
$NAME/examples/Relaxator/makedepend.sh \
$NAME/examples/Relaxator/src/* \
$NAME/examples/TSP/Makefile $NAME/examples/TSP/INSTALL \
$NAME/examples/TSP/runme.sh $NAME/examples/TSP/runviewer.sh \
$NAME/examples/TSP/sciptsp.set \
$NAME/examples/TSP/doc/* \
$NAME/examples/TSP/check/testset/short.* \
$NAME/examples/TSP/src/depend.* \
$NAME/examples/TSP/src/*.cpp $NAME/examples/TSP/src/*.h \
$NAME/examples/TSP/tspviewer/*.java $NAME/examples/TSP/tspdata/*.tsp \
$NAME/examples/VRP/Makefile $NAME/examples/VRP/INSTALL \
$NAME/examples/VRP/doc/* $NAME/examples/VRP/data/* \
$NAME/examples/VRP/src/depend.* \
$NAME/examples/VRP/src/*.cpp $NAME/examples/VRP/src/*.h \
$NAME/examples/VRP/check/check.sh $NAME/examples/VRP/check/testset/* \
$NAME/interfaces/matlab/* \
$NAME/interfaces/ampl/Makefile $NAME/interfaces/ampl/INSTALL $NAME/interfaces/ampl/get.ASL \
$NAME/interfaces/ampl/src/* $NAME/interfaces/ampl/check/check.sh \
$NAME/interfaces/ampl/check/testset/short.test $NAME/interfaces/ampl/check/instances/MINLP/*.col \
$NAME/interfaces/ampl/check/instances/MINLP/*.row $NAME/interfaces/ampl/check/instances/MINLP/*.nl \
$NAME/interfaces/ampl/check/instances/SOS/*.col $NAME/interfaces/ampl/check/instances/SOS/*.row \
$NAME/interfaces/ampl/check/instances/SOS/*.nl $NAME/interfaces/ampl/check/testset/short.solu \
$NAME/interfaces/gams/Makefile $NAME/interfaces/gams/INSTALL $NAME/interfaces/gams/gamsinst.sh \
$NAME/interfaces/gams/test.sh $NAME/interfaces/gams/src/* \
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
$NAME/check/instances/PseudoBoolean/*.cip \
$NAME/check/instances/SAT/*.cnf \
$NAME/check/instances/SOS/*.lp \
$NAME/check/instances/Semicontinuous/*.lp \
$NAME/check/instances/Semicontinuous/*.mps \
$NAME/tests/Makefile \
$NAME/tests/README \
$NAME/tests/include/scip_test.* \
$NAME/tests/src/bugs/depthlevel.c \
$NAME/tests/src/cons/cons.c \
$NAME/tests/src/cons/initlp.c \
$NAME/tests/src/cons/knapsack/solveknapsackapprox.c \
$NAME/tests/src/cons/linear/parsing.c \
$NAME/tests/src/cons/nonlinear/getCoeffsAndConstantFromLinearExpr.c \
$NAME/tests/src/cons/nonlinear/issue1326.c \
$NAME/tests/src/cons/nonlinear/issue1326.cip \
$NAME/tests/src/cons/quadratic/gauge.c \
$NAME/tests/src/cons/quadratic/projection.c \
$NAME/tests/src/heur/multistart.c \
$NAME/tests/src/lpi/bases.c \
$NAME/tests/src/lpi/boundchg.c \
$NAME/tests/src/lpi/matrix.c \
$NAME/tests/src/lpi/solve.c \
$NAME/tests/src/misc/normaldistribution.c \
$NAME/tests/src/misc/regression.c \
$NAME/tests/src/misc/select.c \
$NAME/tests/src/presol/qpkktref.c \
$NAME/tests/src/prop/nlobbt.c \
$NAME/tests/src/sepa/convexproj.c \
$NAME/tests/src/sepa/gauge.c \
$NAME/tests/src/test/stages.c

rm -f $NAME
echo ""
echo "check version numbers in src/scip/def.h, doc/xternal.c, make.project, Makefile.nmake, and makedist.sh ($VERSION):"
grep -H "SCIP_VERSION" src/scip/def.h
grep -H "@version" doc/xternal.c
grep -H "^SCIP_VERSION" make/make.project
grep -H "^VERSION" Makefile.nmake
echo ""
tail src/scip/githash.c
