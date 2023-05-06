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
   solver = "NUOPT";
}
# solver version 
/^MSI Numerical Optimizer/ { solverversion = $4; }

# branch and bound nodes
/^PARTIAL_PROBLEM_COUNT/ { bbnodes = $2;}
# dual and primal bound 
/^VALUE_OF_OBJECTIVE/ { 
    pb = $2;
}
/^GAP/ { 
    gap = $2; 
}

/^SOLUTION_FILE/ { 
    aborted = 0;
    db = pb - gap;
}

# solving status
/^STATUS/ { 
    status = $2;
}

/^ERROR_TYPE/ { 
    type=$2;
    regex = "timeout"
    where = match( $0, regex)
    if (where != 0){
	timeout = 1
    }
}

#END {
#     print "ver:", solverversion;
#     print "pb:", pb;
#     print "db:", db;
#     print "bbnodes:", bbnodes;
#     print "aborted:", aborted;
#     print "timeout:", timeout;
#}
