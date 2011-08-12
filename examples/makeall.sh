#
# run with bash -e makeall.sh to stop on errors
#

EXAMPLES=(Coloring Binpacking Eventhdlr LOP MIPSolver Queens SamplePricer SamplePricer_C TSP VRP)
LPSOLVERS=(clp cpx none spx)
OPTS=(opt dbg)

for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo
    echo ===== $EXAMPLE =====
    echo
    cd $EXAMPLE
    echo
    for OPT in ${OPTS[@]}
    do
	echo make OPT=$OPT LPS=$LPS ZIMPL=false depend
	make OPT=$OPT LPS=$LPS ZIMPL=false depend
	echo
	for LPS in ${LPSOLVERS[@]}
	do
	    echo make OPT=$OPT LPS=$LPS clean
	    make OPT=$OPT LPS=$LPS clean
	    echo
	    echo make OPT=$OPT LPS=$LPS
	    make OPT=$OPT LPS=$LPS
	    echo
	done
    done
    cd -
done
