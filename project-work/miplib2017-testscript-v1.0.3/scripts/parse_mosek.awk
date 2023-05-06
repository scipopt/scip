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
   solver = "MOSEK";
}
# solver version
/^MIPLIBsolverversion/ { solverversion = $3; }
# branch and bound nodes
/^MIPLIBbbnodes/ { bbnodes = $3; }
# dual and primal bound
/^MIPLIBpb/ {
   if( $3 == "infeasible" )
      pb = +infty;
   else if( $3 == "+infty" )
      pb = +infty;
   else
      pb = $3;
}
/^MIPLIBdb/ {
   if( $3 == "-" )
      db = +infty;
   else if( $3 == "-infty" )
      db = -infty;
   else
      db = $3;
}
# solving status
/^MIPLIBaborted/ {
   if( $3 == "TRUE" )
      aborted = 1;
   else
      aborted = 0;
}
/^MIPLIBtimeout/ {
   if( $3 == "TRUE" )
      timeout = 1;
   else
      timeout = 0;
}
