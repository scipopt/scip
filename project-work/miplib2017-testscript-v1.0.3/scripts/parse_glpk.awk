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
   solver = "GLPK";
   feasible = 1;
}
# solver version
/^GLPSOL: GLPK LP\/MIP Solver,/ { solverversion = $5; }

#
# solution
#
/^PROBLEM HAS NO INTEGER FEASIBLE SOLUTION/ {
   db = +infty;
   pb = +infty;
   feasible = 0;
   aborted = 0;
}
/^PROBLEM HAS NO FEASIBLE SOLUTION/ {
   db = +infty;
   pb = +infty;
   feasible = 0;
   aborted = 0;
}
/^PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION/ {
   db = +infty;
   pb = +infty;
   feasible = 0;
   aborted = 0;
}
/^TIME LIMIT EXCEEDED/ {
   timeout = 1;
   aborted = 0;
}
/^INTEGER OPTIMAL SOLUTION FOUND/ {
   aborted = 0;
   db = pb; # fix missing display of optimal dual bound due to the message "tree is empty" (version 4.39)
}
/mip = / {
   if ( feasible == 1 )
   {
      if ( $1 == "+" )
      {
         if ( $5 == "not" )
	 {
            pb = +infty;
            if ( $9 != "tree" )
               db = ($9 == "-inf") ? -infty : $9;
            bbnodes = $11;
         }
         else {
            pb = $5;
            if ( $7 == "tree" )
               bbnodes = $12;
            else
	    {
               db = $7;
               bbnodes = $10;
            }
         }
      }
      else {
         if ( $4 == "not" )
	 {
            pb = +infty;
            if ( $8 != "tree" )
               db = ($8 == "-inf") ? -infty : $8;
            bbnodes = $10;
         }
         else
	 {
            pb = $4;
            if ( $6 == "tree" )
               bbnodes = $11;
            else
	    {
               db = $6;
               bbnodes = $9;
            }
         }
      }
   }
}
/>>>>>/ {
   if ( $1 == "+" )
   {
      if ( feasible == 1 )
      {
         if ( $4 == "not" )
	 {
            pb = +infty;
            db = ($8 == "-inf") ? -infty : $8;
            bbnodes = $10;
         }
         else
	 {
            pb = $4;
            if ( $6 == "tree" )
               bbnodes = $11;
            else
	    {
               db = $6;
               bbnodes = $9;
            }
         }
      }
      else
      {
         if ( $3 == "not" )
	 {
            pb = +infty;
            db = ($7 == "-inf") ? -infty : $7;
            bbnodes = $9;
         }
         else
	 {
            pb = $3;
            if ( $5 == "tree" )
               bbnodes = $10;
            else
	    {
               db = $5;
               bbnodes = $8;
            }
         }
      }
   }
}
