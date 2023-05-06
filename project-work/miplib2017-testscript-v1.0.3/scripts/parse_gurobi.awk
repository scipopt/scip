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

BEGIN {
  solver = "Gurobi";
}
# solver version
/^Gurobi Optimizer version/ { solverversion = $4; }
# branch and bound nodes
/^Explored/ { bbnodes = $2; }
# infeasible model
/^Model is infeasible/ {
  db = pb;
}
# dual and primal bound
/^Best objective/ {
 if( $3 != "-," )
  pb = $3 + 0.0;
 if( $6 != "-," )
  db = $6 + 0.0;
}
/^Optimal objective/ {
 pb = $3;
 db = $3;
}
# solving status
/^Explored/ { aborted = 0; }
/^Solved in/ { aborted = 0; }
/^Time limit reached/ { timeout = 1; }
/^Solve interrupted/ { timeout = 1; }
/^Optimal solution found/ { timeout = 0; }
