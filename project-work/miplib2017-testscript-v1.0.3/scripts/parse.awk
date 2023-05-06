#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

function abs(x)
{
   return x < 0 ? -x : x;
}
function min(x,y)
{
   return (x) < (y) ? (x) : (y);
}
function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}
BEGIN {
   printf("----------------------------+----------------+----------------+------+---------+-------+--------+---------\n");
   printf("Name                        |   Dual Bound   |  Primal Bound  | Gap%% |  Nodes  |  Time | Status | Solution \n");
   printf("----------------------------+----------------+----------------+------+---------+-------+--------+---------\n");

   infty = +1e+20;
   eps = 1e-04;
   largegap = 1e+04;
   reltol = 1e-06

   # initialize summary data
   nsolved = 0;
   nstopped = 0;
   nfailed = 0;

   # initialize data to be set in parse_<solver>.awk
   solver = "?";
   solverversion = "?";
   solverremark = "";

   # initialize MIPLIB script version
   miplibversion = "?";
}

# parse solution status from solution file (optional)
/=opt=/  { knownsolstatus[$2] = "opt"; knownsolvalue[$2] = $3; }   # get optimum
/=inf=/  { knownsolstatus[$2] = "inf"; }                      # problem infeasible (no feasible solution exists)
/=best=/ { knownsolstatus[$2] = "best"; knownsolvalue[$2] = $3; }  # get best known solution value
/=feas=/ { knownsolstatus[$2] = "feas"; }                     # no feasible solution known
/=unkn=/ { knownsolstatus[$2] = "unkn"; }                     # no feasible solution known

# MIPLIB script version
/^MIPLIB script version/ {
   miplibversion = $4;
}

# instance name 
/^@01/ { 
   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   n = split(prob, b, "#");
   prob = b[1];

   # initialize data to be set in parse.awk
   timelimit = 0;
   starttime = 0.0;
   endtime = 0.0;
   time = 0.0;

   # initialize data to be set parse_<solver>.awk
   bbnodes = 0;
   pb = +infty;
   db = -infty;
   aborted = 1;
   timeout = 0;
   solstatus = "error";
}
# time
/@03/ { starttime = $2; }
/@04/ { endtime = $2; }
/@05/ { timelimit = $2; }
# solution status
/Read SOL:/ { 
   solstatus = "--";
}
/Check SOL:/ { 
   intcheck = $4;
   conscheck = $6;
   objcheck = $8;
   if( intcheck && conscheck && objcheck ) 
      solstatus = "ok";
   else
      solstatus = "fail";
}
/^=ready=/ {
   # measure wallclock time externaly rounded up to the next second 
   time = max(1, endtime - starttime);

   # determine solving status
   status = "";
   if( aborted ) 
      status = "abort";
   else if( timeout )
      status = "stopped";
   else
      status = "ok";

   # fix solstatus w.r.t. consistency to solution file (optional) 
   if( knownsolstatus[prob] != "" && solstatus != "fail" && solstatus != "error" )
   {
      knownsolstat = knownsolstatus[prob];
      knownsol =  knownsolvalue[prob];
      abstol = LINTOL;

      # check if solver claims infeasible (NOTE currently only working with MINIMIZATION)
      infeasible = (pb >= infty && db >= infty);
      
      # We check different know solution status. 
      #
      # 1) in case of knowing an optimal solution, we check if the solution
      #    reported by the solver in consistent with tolerances
      # 2) in case of knowing infeasbility, the inconsistenty would be
      #    detected by the solution checker
      # 3) in case of knowing nothing (=unkn=) no consistency check
      #     possible
      # 4) in case of knowing a feasible solution, we check the feasibility
      #    status and optimal solution status for consistency
      # 5) in case of knowing feasibility, we check the feasibility
      #    status 
      if( knownsolstat == "opt" || knownsolstat == "best" )
      {
	 if( infeasible )
	 {
	    # solver claimed infeasible, but we know a solution
	    solstatus = "mismatch";
	 }
	 else if( !timeout && !aborted )
	 {
	    # solver solved instance to optimality but the solution value
	    # is inconsistent to the know (optimal) solution
	    if( abs(knownsol) <= 100 )
	    {
	       if( pb - knownsol > abstol )
		  solstatus = "mismatch";
	    }
	    else if( (pb - knownsol)/abs(knownsol) > reltol )
	       solstatus = "mismatch";
	 }
      }	    
      else if( knownsolstat == "feas" )
      {
	 if( infeasible )
	 {
	    # solver claimed infeasible, but we know a solution
	    solstatus = "mismatch";
	 }
      }
   }

   # determine overall status from solving status and solution status:

   # instance solved correctly (including case that no solution was found) 
   if( status == "ok" && (solstatus == "ok" || solstatus == "--") )
      nsolved++;
   # incorrect solving process or infeasible solution (including errors with solution checker)
   else if( status == "abort" || (solstatus == "fail" || solstatus == "error" || solstatus == "mismatch") )
      nfailed++;
   # stopped due to imposed limits
   else
      nstopped++;
 
  # compute gap
   temp = pb;
   pb = 1.0*temp;
   temp = db;
   db = 1.0*temp;

   if( abs(pb - db) < eps && pb < +infty ) 
      gap = 0.0;
   else if( abs(db) < eps || abs(pb) < eps)
      gap = -1.0;
   else if( pb*db < 0.0 )
      gap = -1.0;
   else if( abs(db) >= +infty )
      gap = -1.0;
   else if( abs(pb) >= +infty )
      gap = -1.0;
   else 
      gap = 100.0*abs((pb-db)/min(abs(db), abs(pb)));

   if( gap < 0.0 )
      gapstr = "    --";
   else if( gap < largegap )
      gapstr = sprintf("%6.1f", gap);
   else
      gapstr = " Large";

   printf("%-28s %16.9g %16.9g %6s %9d %7d %8s %9s\n", prob, db, pb, gapstr, bbnodes, time, status, solstatus);
}
END {
   printf("----------------------------+----------------+----------------+------+---------+-------+--------+---------\n");
   printf("\n");
   printf("solved/stopped/failed: %d/%d/%d\n", nsolved, nstopped, nfailed);
   printf("\n");
   printf("@03 MIPLIB script version %s\n", miplibversion);
   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 %s(%s)%s\n", solver, solverversion, solverremark);

}
