for i in Coloring LOP MIPSolver Queens SamplePricer SamplePricer_C TSP VRP
do
    echo
    echo ===== $i =====
    echo
    cd $i
    echo make OPT=dbg depend
    make OPT=dbg depend
    echo
    echo make OPT=opt depend
    make OPT=opt depend
    echo
    echo make OPT=dbg clean
    make OPT=dbg clean
    echo
    echo make OPT=dbg
    make OPT=dbg
    echo
    echo make OPT=opt clean
    make OPT=opt clean
    echo
    echo make OPT=opt
    make OPT=opt
    echo
    echo make OPT=dbg LPS=clp clean
    make OPT=dbg LPS=clp clean
    echo
    echo make OPT=dbg LPS=clp
    make OPT=dbg LPS=clp
    echo
    echo make OPT=opt LPS=clp clean
    make OPT=opt LPS=clp clean
    echo
    echo make OPT=opt LPS=clp
    make OPT=opt LPS=clp
    echo
    echo make OPT=dbg LPS=cpx clean
    make OPT=dbg LPS=cpx clean
    echo
    echo make OPT=dbg LPS=cpx
    make OPT=dbg LPS=cpx
    echo
    echo make OPT=opt LPS=cpx clean
    make OPT=opt LPS=cpx clean
    echo
    echo make OPT=opt LPS=cpx
    make OPT=opt LPS=cpx
    echo
    echo make OPT=dbg LPS=none clean
    make OPT=dbg LPS=none clean
    echo
    echo make OPT=dbg LPS=none
    make OPT=dbg LPS=none
    echo
    echo make OPT=opt LPS=none clean
    make OPT=opt LPS=none clean
    echo
    echo make OPT=opt LPS=none
    make OPT=opt LPS=none
    cd -
done
