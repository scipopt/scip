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
   solver = "XPRESS";
   solverremark = "";
   solver_objsense = +1;
   solver_type = "MIP";
}
# solver version for version <= 27.01.01
/^FICO Xpress Optimizer.* v[0-9][0-9]/ { solverversion = $(NF-2); }

# solver version for version > 27.01.01
/^FICO Xpress-Optimizer.* v[0-9][0-9]/ { solverversion = $(NF-2); }

# solver version for version >= 32.01 (Suite version >= 8.3)
/^FICO Xpress Solver.* v[0-9][0-9]/ { solverversion = $5; }

# objective sense
/^Minimizing MI[A-Z]*P/ {
   solver_objsense = +1;
}
/^Maximizing MI[A-Z]*P/ {
   solver_objsense = -1;
}

# Search completion
/Search completed/ {
  aborted = 0;
  timeout = 0;
}
/Search unfinished/ {
  aborted = 0;
  timeout = 1;
}


# branch and bound nodes
/^Nodes explored = / {
   bbnodes = $NF;
}

# Objective value
/^Objective value = / {
  pb = $NF
}

# Bound value
/^Best Bound = / {
  db = $NF
}
