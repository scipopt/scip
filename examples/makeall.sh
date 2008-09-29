for i in Coloring LOP MIPSolver Queens SamplePricer SamplePricer_C TSP VRP
do
    echo
    echo ===== $i =====
    echo
    cd $i
    make OPT=dbg depend
    make OPT=opt depend
    make OPT=dbg clean
    make OPT=dbg
    make OPT=opt clean
    make OPT=opt
    make OPT=dbg LPS=clp clean
    make OPT=dbg LPS=clp
    make OPT=opt LPS=clp clean
    make OPT=opt LPS=clp
    make OPT=dbg LPS=cpx clean
    make OPT=dbg LPS=cpx
    make OPT=opt LPS=cpx clean
    make OPT=opt LPS=cpx
    make OPT=dbg LPS=none clean
    make OPT=dbg LPS=none
    make OPT=opt LPS=none clean
    make OPT=opt LPS=none
    cd -
done
