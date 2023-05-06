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
  solver = "MIPCL";
}
# solver version
/^MIPCL version/ { solverversion = $3; }
# branch and bound nodes
/^Branch-and-Cut nodes:/ {
	bbnodes = $3;
	aborted = 0;
 }
# infeasible model
/^This problem is infeasible/ {
	db = pb;
}
# dual and primal bound
/^Objective value:/ {
	pb = $3 + 0.0;
	if ($6 == "proven") {
		db=pb;
	}
}

# solving status
/^Time limit reached/ {
 	timeout = 1;
}
/lower-bound:/ {
	db = $2 + 0.0;
}
