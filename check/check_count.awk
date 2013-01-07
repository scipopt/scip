#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    checkcount.awk
#@brief   SCIP Check Report Generator for counting tests
#@author  Stefan Heinz
#
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
   useshortnames = 0;  # should problem name be truncated to fit into column?
   printf("SCIP %s\n", LPS);
   printf("SETTING %s\n", SETTING);
   printf("------------------------+------+---------------+--------+--------+----------+-------+------------------+--------\n");
   printf("Name                    | Type | Conss |  Vars |Conftime| FeasST |   Nodes  |  Time |    Solutions     |\n");
   printf("------------------------+------+-------+-------+--------+--------+----------+-------+------------------+--------\n");
  
   nprobs = 0;
   sbab = 0;
   stottime = 0.0;
   timeouts = 0;
   memouts = 0;
   readerrors = 0; 
   failtime = 0.0;
   fail = 0;
   pass = 0;
   settings = "default";
   timelimit = 0.0;
   memlimit = 0;
}

# read number of solutions from the solution file for comparison
/^=sol=/  {  knownsols[$2] = $3; }

/^@01/ { 
   filename = $2;

   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > 24 )
      shortprob = substr(prob, length(prob)-23, 24);
   else
      shortprob = prob;

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];
 
   vars = 0;
   binvars = 0;
   intvars = 0;
   implvars = 0;
   contvars = 0;
   cons = 0;
   origvars = 0;
   origcons = 0;
   interrupted = 0;
   simpiters = 0;
   bbnodes = 0;
   tottime = 0.0;
   tottimestr = "";
   readerror = 0;
   counted = 0;
   timelimitreached = 0;
   memlimitreached = 0;
   timelimit = 0.0;
   starttime = 0.0;
   endtime = 0.0;
   inoriginalprob = 1;
   noncuts = 0;
   feasST = 0;
   sols = 0;
   discards = 0;
   maxdepth = 0;
   conftime = 0.0;
   confclauses = 0;
}
/@03/ { starttime = $2; }
/@04/ { endtime = $2; }
/^loaded parameter file/ { settings = $4; sub(/<..\/settings\//, "", settings); sub(/\.set>/, "", settings); }
/^SCIP> loaded parameter file/ { settings = $5; sub(/<..\/settings\//, "", settings); sub(/\.set>/, "", settings); }
/^parameter <limits\/time> set to/ { timelimit = $5; }
/^SCIP> parameter <limits\/time> set to/ { timelimit = $6; }
/^SCIP> set limits memory/ { memlimit = $5; }

#
# problem size
#
/^Presolved Problem  :/ { inoriginalprob = 0; }
/^  Variables        :/ {
   if( inoriginalprob )
      origvars = $3;
   else {
      vars = $3;
      intvars = $6;
      implvars = $8;
      contvars = $11;
      binvars = vars - intvars - implvars - contvars;
   }
}
/^  Constraints      :/ {
   if( inoriginalprob )
      origcons = $3;
   else
      cons = $3;
}
#
# solution
#
/^Original Problem   : no problem exists./ { readerror = 1; }
/solving was interrupted/  { interrupted = 1; counted = 0; }
/memory limit reached/ { memlimitreached = 1; }
/time limit reached/  { timelimitreached = 1; }
/problem is solved/   { interrupted = 0; }
/infeasible/    { counted = 1; }

#
# iterations
#
/^  primal LP        :/ { simpiters += $6; }
/^  dual LP          :/ { simpiters += $6; }
/^  barrier LP       :/ { simpiters += $6; }
/^  nodes \(total\)    :/ { bbnodes = $4 }
/^  max depth \(total\):/ { maxdepth = $4 }

#
# conflict analysis
#
/^Conflict Analysis  :/ { inconflict = 1; }
/^  propagation      :/ {
   if( inconflict == 1 ) {
      conftime += $3;
   }
}
/^  infeasible LP    :/ {
   if( inconflict == 1 ) {
      conftime += $4;
   }
}
/^  strong branching :/ {
   if( inconflict == 1 ) {
      conftime += $4; 
   }
}
/^  pseudo solution  :/ {
   if( inconflict == 1 ) {
      conftime += $4; 
   }
}
/^  applied globally :/ {
   if( inconflict == 1 ) {
      confclauses += $7; 
   }
}
/^  applied locally  :/ {
   if( inconflict == 1 ) {
      confclauses += $7; 
   }
}

/^Separators         :/ { inconflict = 0; }

#
# time
#
/^Solving Time       :/ { tottime = $4 } # for older scip version ( < 2.0.1.3 )
/^  solving          :/ { tottime = $3 }

#
# Solutions Colectd
#
/^Feasible Solutions :/ { 
   sols = $4;  
   feasST = substr($5, 2, length($5)-1); 
}

#
# Output
#
/^=ready=/ {
   nprobs++;
  
   if( !readerror ) {
      # figure problem type out 
      if( vars == 0 )
         probtype = "--";
      else if( binvars == 0 && intvars == 0 )
         probtype = "LP";
      else if( contvars == 0 ) {
         if( intvars == 0 && implvars == 0 )
            probtype = "BP";
         else
            probtype = "IP";
      }
      else {
         if( intvars == 0 )
            probtype = "01MIP";
         else
            probtype = "MIP";
      }

      if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

      stottime += tottime;
      sbab += bbnodes;
      tottimestr = sprintf("%.1f", tottime);

      if( interrupted ) {
         fail++;
         if(memlimitreached) {
            memouts++;
         }
         if( timelimitreached ) {
            tottimestr = "timeout";
            timeouts++;
         }
      }      
      else if( counted ) {
         pass++;
      }
      else {
         tottimestr = "error";
      }
  
      if( knownsols[prob] == "" ) {
         status = "unknown";
      }
      else if( (interrupted && knownsols[prob] >= sols) ) {
         status = "ok";
      }
      else if( !interrupted &&  knownsols[prob] == sols ) {
         status = "ok";
      }
      else {
         status = "fail";
         swrong++;
      }

      if( sols > 1e+16 )
         solsstr = sprintf("%e", sols);
      else
         solsstr = sprintf("%d", sols);

      if( interrupted ) {
         texsolsstr = sprintf("$\\geq$%s", solsstr);
         solsstr = sprintf(">%s", solsstr);
      }
      else
         texsolsstr = sprintf("%s", solsstr);
     
      printf("%-25s %-5s %7d %7d %8.1f %8d %10s %7.1f %18s %7s\n",
         shortprob, probtype, origcons, origvars,
         conftime, feasST, bbnodes, tottimestr, solsstr, status);
   }
   else { 
      # read error
      readerrors++;
     
      printf("%-25s %-5s %7d %7d %8d %9d %7s %8d %8d %8d %18s %6s\n",
         shortprob, "-", "-", "-", 
         "-", "-", "-", "-", "-", "-", "readerror", "-");
   }
}
END {
   printf("------------------------+------+-------+-------+--------+--------+----------+-------+------------------+--------\n");

   printf("\n");
   printf("-----------+------[Out]------+-----------\n");
   printf("  Cnt  Pass  Time   Mem  Read  Total Time\n");
   printf("-----------+-----+-----+-----+-----------\n");
   printf("%5d %5d %5d %5d %5d %11.1f \n",
      nprobs, pass, timeouts, memouts, readerrors, stottime);
   printf("-----------+-----+-----+-----+-----------\n");

   printf("time limit: %g\n", timelimit);
   printf("memory limit: %g\n", memlimit);
   printf("SCIP:%s.set\n", settings);
}
