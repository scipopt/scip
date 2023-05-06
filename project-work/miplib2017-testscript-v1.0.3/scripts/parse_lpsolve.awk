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
   solver = "lpsolve";
}
# solver version
/^ MEMO:/  { solverversion = $4; }

#
# solution
#
/^Optimal solution/ {
   pb = $3;
   bbnodes = $7;
   aborted = 0;
}
/^Relaxed solution/ {
   db = $3;
}
/^lp_solve optimization was stopped due to time-out./ {
   timeout = 1;
}
/^Suboptimal solution/ {
   timeout = 1;
}
END {
   if( !timeout )
      db = pb;
}
