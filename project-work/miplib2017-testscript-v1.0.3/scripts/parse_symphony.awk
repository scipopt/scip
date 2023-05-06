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
   solver = "SYMPHONY";
   optimal = 0;
}
# solver version
/^==  Version:/  {
   solverversion = $3;
}
#
# solution
#
/^\* Optimal Solution Found/ {
   aborted = 0;
   optimal = 1;
}
/^Number of analyzed nodes:/ {
   bbnodes = $5;
}
/^Current Upper Bound:/ {
   pb = $4;
}
/^Current Lower Bound:/ {
   db = $4;
}
/^Solution Cost:/ {
   pb = $3;
   if( optimal )
      db = pb;
}
/^Preprocessing detected infeasibility/ {
   db = +infty;
   pb = +infty;
   aborted = 0;
}
/^The problem is infeasible!/ {
   db = +infty;
   pb = +infty;
   aborted = 0;
}
/^\* Time Limit Reached/ {
   aborted = 0;
   timeout = 1;
}
