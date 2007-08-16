for i in MIPSolver SamplePricer SamplePricer_C TSP
do
    echo
    echo ===== $i =====
    echo
    cd $i
    make OPT=opt depend
    make OPT=dbg depend
    make OPT=opt
    make OPT=dbg
    cd -
done
