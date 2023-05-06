#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# set all solver specific data:
#  solver ["?"]
#  solverversion ["?"]
#  solverremark [""]
#  bbnodes [0]
#  db [-infty]
#  pb [+infty]
#  aborted [1]
#  timeout [0]

# solver name
BEGIN {
   solver = "MATLAB";
   solverremark = "";
   solver_objsense = +1;
   solver_type = "MIP";
   db = -infty;
   pb = +infty;
   timeout = 0;
}
# solver version
/^Release/ { solverversion = $2; }

# branch and bound nodes
/^BBnodes/ {
   bbnodes = $2;
   aborted = 0;
}

/exceeded the iteration limit/ {
   timeout = 1;
}

/exceeded the time limit/ {
   timeout = 1;
}

/exceeded its allocated memory/ {
   timeout = 1;
}

/reached the maximum number of nodes/ {
   timeout = 1;
}

# dual bound before BB starts
/LP: *Optimal objective/ {
   db = $6 + 0.0;
}

# dual and primal bound
/^PrimalBound/ {pb = $2};
/^DualBound/ {db = $2};

# solver detected infeasibility
/No feasible solution found./ {
   pb = +infty;
   db = +infty;
}
