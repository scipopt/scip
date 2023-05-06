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
   solver = "FiberSCIP";
   solverversion = "";
   solverremark = "";
}
/timeout/ { timeout_timelimit = $2; }
# solver version
/^SCIP version/ { solverversion = $3; }
# LP solver name and version
/^SCIP version/ {
   if( $13 == "SoPlex" )
      lps = "spx";
   else if( $13 == "CPLEX" )
      lps = "cpx";
   else if( $13 == "NONE" )
      lps = "none";
   else if( $13 == "Clp" )
      lps = "clp";
   else if( $13 == "MOSEK" )
      lps = "msk";
   else if( $13 == "Gurobi" )
      lps = "grb";
   else if( $13 == "QSopt" )
      lps = "qso";

   if( NF >= 14 )
   {
      split($14, v, "]");
      lpsversion = v[1];
   }

   solverremark = lps "(" lpsversion ")";
}
# branch and bound nodes
/^  nodes \(total\)    :/ { bbnodes = $4; }
# dual and primal bound
/^  Primal Bound     :/ {
   if( $4 == "infeasible" )
   {
      pb = +infty;
      db = +infty;
   }
   else if( $4 == "-" )
      pb = +infty;
   else
      pb = $4;
}
/^  Dual Bound       :/ {
   if( $4 != "-" )
      db = $4;
}
# solving status
/^SCIP Status        :/ { aborted = 0; }
/solving was interrupted/ { timeout = 1; }
/problem is solved/ { timeout = 0; }
$1 ~ /^[0-9]+$/  && $1 > timeout_timelimit {
   timeout = 1;
   aborted = 0;
   if( $5 == "-" )
      pb = +infty;
   else
      pb = $5;
   if( $6 != "-" )
      db = $6;
}
$1 == "*" && $2 ~ /^[0-9]+$/  && $2 > timeout_timelimit {
   timeout = 1;
   aborted = 0;
   if( $6 == "-" )
      pb = +infty;
   else
      pb = $6;
   if( $7 != "-" )
      db = $7;
}
$1 ~ /^[0-9]+$/  {
   if( $5 == "-" )
      pb = +infty;
   else
      pb = $5;
   if( $6 != "-" )
      db = $6;
}
$1 == "*" && $2 ~ /^[0-9]+$/ {
   if( $6 == "-" )
      pb = +infty;
   else
      pb = $6;
   if( $7 != "-" )
      db = $7;
}
