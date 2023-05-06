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
   solver = "CPLEX";
   solverremark = "";
   solver_objsense = +1;
   solver_type = "MIP";
   timeout = 0;
}
# solver version
/^Welcome to.* Interactive Optimizer/ { solverversion = $NF; }

# objective sense
/^Problem is a minimization problem./ {
   solver_objsense = +1;
}
/^Problem is a maximization problem./ {
   solver_objsense = -1;
}

# branch and bound nodes
/^Solution time/ {
   bbnodes = $11;
   aborted = 0;
}

/^MIP - Error termination/ {
   timeout = 1;
}

# dual and primal bound, and solving status
/^MIP - / { solver_type = $1; $0 = substr($0, 7, length($0)-6); }
/^Barrier - / { solver_type = $1; $0 = substr($0, 11, length($0)-10); }
/^Primal simplex - / { solver_type = $1; $0 = substr($0, 18, length($0)-17); }
/^Dual simplex - / { solver_type = $1; $0 = substr($0, 16, length($0)-15); }
/^Populate - / { solver_type = "MIP"; $0 = substr($0, 12, length($0)-11); }
/^Presolve - / { solver_type = $1; $0 = substr($0, 12, length($0)-11); }
/^Integer /  {
   if ($2 ~ /infeasible/ )
   {
      db = solver_objsense * infty;
      pb = solver_objsense * infty;
   }
   else
   {
      db = $NF;
      pb = $NF;
   }
}
/^Optimal:  Objective = / {
   pb = $NF;
   db = $NF;
}
/^Non-optimal:  Objective = / {
   pb = $NF;
   db = $NF;
}
/^Unbounded or infeasible./ {
   db = solver_objsense * infty;
   pb = solver_objsense * infty;
}
/^Infeasible/ {
   db = solver_objsense * infty;
   pb = solver_objsense * infty;
}
/^Unbounded/ {
   db = solver_objsense * infty;
   pb = solver_objsense * infty;
}
/^Dual objective limit exceeded/ {
   db = solver_objsense * infty;
   pb = solver_objsense * infty;
}
/^Primal objective limit exceeded/ {
   db = solver_objsense * infty;
   pb = solver_objsense * infty;
}
/^Solution limit exceeded/ {
   pb = $NF;
}
/^Time limit exceeded/  {
   pb = solver_objsense * infty;
   if ( solver_type == "MIP" && $4 != "no" )
      pb = $8;
   timeout = 1;
}
/^Memory limit exceeded/ {
   pb = solver_objsense * infty;
   if ( solver_type == "MIP" && $4 != "no" )
      pb = $8;
   timeout = 1;
}
/^Node limit exceeded/ {
   pb = ($4 == "no") ? solver_objsense * infty : $8;
   timeout = 1;
}
/^Tree /  {
   pb = ($4 == "no") ? solver_objsense * infty : $8;
   timeout = 1;
}
/^Aborted, / {
   pb = solver_objsense * infty;
   if ( solver_type == "MIP" && $2 != "no" )
      pb = $6;
   aborted = 1;
}
/^Error /  {
   pb = solver_objsense * infty;
   if ( solver_type == "MIP" && $3 != "no" )
      pb = $7;
   aborted = 1;
}
/^Unknown status / {
   pb = solver_objsense * infty;
   if ( solver_type == "MIP" && $4 == "Objective" )
      pb = $6;
}
/^All reachable solutions enumerated, integer / {
   if ($6 ~ /infeasible/ )
   {
      db = solver_objsense * infty;
      pb = solver_objsense * infty;
   }
   else
   {
      db = $NF;
      pb = $NF;
   }
}
/^Populate solution limit exceeded/ {
   pb = $NF;
}
/^Current MIP best bound =/ {
   db = $6;
}
/^CPLEX Error  1001: Out of memory./ {
   timeout = 1;
}
