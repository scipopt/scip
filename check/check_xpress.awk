#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    check_xpress.awk
#@brief   XPRESS Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Marc Pfetsch
#@author  Robert Waniek
#@author  Michael Winkler
#@author  Marc Pfetsch
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
   timegeomshift = 60.0;
   nodegeomshift = 1000.0;
   onlyinsolufile = 0;  # should only instances be reported that are included in the .solu file?
   useshortnames = 1;   # should problem name be truncated to fit into column?
   writesolufile = 0;   # should a solution file be created from the results
   NEWSOLUFILE = "new_solufile.solu";
   infty = +1e+20;
   headerprinted = 0;
   namelength = 18;             # maximal length of instance names (can be increased)
   timeline = 0;        # only used as flag to read next line once

   nprobs   = 0;
   sbab     = 0;
   scut     = 0;
   stottime = 0.0;
   sgap     = 0.0;
   nodegeom = 0.0;
   timegeom = 0.0;
   shiftednodegeom = nodegeomshift;
   shiftedtimegeom = timegeomshift;
   failtime = 0.0;
   timeouttime = 0.0;
   fail     = 0;
   pass     = 0;
   timeouts = 0;
   settings = "default";
   version = "?";
   threads = 1;
   starttime = 0.0;
   endtime = 0.0;
   timelimit  = 0.0;
}
/=opt=/  { solstatus[$2] = "opt"; sol[$2] = $3; }   # get optimum
/=inf=/  { solstatus[$2] = "inf"; }                 # problem infeasible (no feasible solution exists)
/=best=/ { solstatus[$2] = "best"; sol[$2] = $3; }  # get best known solution value
/=unkn=/ { solstatus[$2] = "unkn"; }                # no feasible solution known
#
# problem name
#
/^@01/ {
   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > namelength )
      shortprob = substr(prob, length(prob)-namelength-1, namelength);
   else
      shortprob = prob;

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];

   # for marking section (is repeated with bounds in newer CPLEX versions)
   origdata = 0;
   presdata = 0;

   origvars   = 0;
   origcons   = 0;
   vars       = 0;
   cons       = 0;
   timeout    = 0;
   memout     = 0;
   nodeout    = 0;
   opti       = 0;
   feasible   = 1;
   cuts       = 0;
   pb         = +infty;
   db         = -infty;
   mipgap     = 0;
   absmipgap  = 1e-6;
   bbnodes    = 0;
   iters      = 0;
   primlps    = 0;
   primiter   = 0;
   duallps    = 0;
   dualiter   = 0;
   sblps      = 0;
   sbiter     = 0;
   tottime    = 0.0;
   aborted    = 1;
}

/@03/ { starttime = $2; }
/@04/ { endtime = $2; }

/^FICO Xpress Solver 64bit/ { version = $5; }

/^timelimit/ {
   timelimit = abs($2);
}
/^mipgap/ {
   mipgap = abs($2);
}

#
# problem size
#
/^Original problem has:/ {
   origdata = 1;
}
/^Presolved problem has:/ {
   presdata = 1;
}

# read lines that start with up to 10⁹-1 rows
/^  / {
   if ( origdata == 1 ) {
      origcons = $1;
      origvars = $3;
      origdata = 0;
   }
   else if ( presdata == 1 ) {
      cons = $1;
      vars = $3;
      # no binary variables
      intvars = $7;
      contvars = vars - intvars;
      presdata = 0;
   }
}
/^Best integer solution found is/ {
   pb = $6;
   feasible = 1;
}
/^Best bound is/ {
   if( feasible == 1 )
      db = $4;
}
/^Problem is integer infeasible/ {
   db = infty;
   pb = infty;
   feasible = 0;
}
/^Optimal solution found/ {
   aborted = 0;
   feasible = 1;
}
/^Barrier solved problem/ {
   aborted = 0;
   feasible = 1;
   timeline = 1;
   next;
}
// {if( timeline )
{
   tottime = substr($5,0,length($5)-1) * 1.0;
   timeline = 0;
}
}
# this is not "per simplex iteration" but "per Xpress iteration"
/ microseconds per iteration/ {
   if( tottime == 0.0 )
      tottime = $1 / 1000000;
}
/^ \*\*\* Search completed \*\*\*     Time:/ {
   bbnodes = $8;
   tottime = $6;
   aborted = 0;
}
/^ \*\*\* Search unfinished \*\*\*    Time:/ {
   bbnodes = $8;
   tottime = $6;
   aborted = 0;
   timeout = 1;
}
/^simplexiter/{
   iters = 1;
   next;
}

#
# evaluation
#
# solver status overview (in order of priority):
# 1) solver broke before returning solution => abort
# 2) solver cut off the optimal solution (solu-file-value is not between primal and dual bound) => fail
#    (especially of problem is claimed to be solved but solution is not the optimal solution)
# 3) solver solved problem with the value in solu-file (if existing) => ok
# 4) solver solved problem which has no (optimal) value in solu-file => solved
# 5) solver found solution better than known best solution (or no solution was known so far) => better
# 6) solver reached any other limit (like time or nodes) => timeout
# 7) otherwise => unknown
#
/^=ready=/ {

   #since the header depends on the parameter settings it is no longer possible to print it in the BEGIN section
   if( !headerprinted ) {
      printf("\\documentclass[leqno]{article}\n")                      >TEXFILE;
      printf("\\usepackage{a4wide}\n")                                 >TEXFILE;
      printf("\\usepackage{amsmath,amsfonts,amssymb,booktabs}\n")      >TEXFILE;
      printf("\\usepackage{supertabular}\n")                           >TEXFILE;
      printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
      printf("\\begin{document}\n\n")                                  >TEXFILE;
      printf("\\begin{center}\n")                                      >TEXFILE;
      printf("\\setlength{\\tabcolsep}{2pt}\n")                        >TEXFILE;
      printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n")    >TEXFILE;
      printf("\\tablehead{\n\\toprule\n")                              >TEXFILE;
      printf("Name                &  Conss &   Vars &     Dual Bound &   Primal Bound &  Gap\\%% &     Nodes &     Time \\\\\n") >TEXFILE;
      printf("\\midrule\n}\n")                                         >TEXFILE;
      printf("\\tabletail{\n\\midrule\n")                              >TEXFILE;
      printf("\\multicolumn{%d}{r} \\; continue next page \\\\\n", ntexcolumns) >TEXFILE;
      printf("\\bottomrule\n}\n")                                      >TEXFILE;
      printf("\\tablelasttail{\\bottomrule}\n")                        >TEXFILE;
      printf("\\tablecaption{CPLEX with %s settings}\n",settings)      >TEXFILE;
      printf("\\begin{supertabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrr@{}}\n") >TEXFILE;

      # prepare header
      hyphenstr = "";
      for (i = 0; i < namelength; ++i)
         hyphenstr = sprintf("%s-", hyphenstr);

      # first part: name of given length
      tablehead1 = hyphenstr;
      tablehead2 = sprintf("Name%*s", namelength-4, " ");
      tablehead3 = hyphenstr;

      # append rest of header
      tablehead1 = tablehead1"+------+--- Original --+-- Presolved --+----------------+----------------+------+---------+--------+-------+";
      tablehead2 = tablehead2"| Type | Conss |  Vars | Conss |  Vars |   Dual Bound   |  Primal Bound  | Gap%% |  Iters  |  Nodes |  Time |";
      tablehead3 = tablehead3"+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+";

      tablehead1 = tablehead1"--------\n";
      tablehead2 = tablehead2"       \n";
      tablehead3 = tablehead3"--------\n";

      printf(tablehead1);
      printf(tablehead2);
      printf(tablehead3);

      headerprinted = 1;
   }

   if( !onlyinsolufile || solstatus[prob] != "" ) {

      #avoid problems when comparing floats and integer (make everything float)
      temp = pb;
      pb = 1.0*temp;
      temp = db;
      db = 1.0*temp;

      bbnodes = max(bbnodes, 1);  # in case solver reports 0 nodes if the primal heuristics find the optimal solution in the root node

      nprobs++;

      optimal = 0;
      markersym = "\\g";
      if( abs(pb - db) < 1e-06 && pb < infty ) {
         gap = 0.0;
         optimal = 1;
         markersym = "  ";
      }
      else if( abs(db) < 1e-06 )
         gap = -1.0;
      else if( abs(pb) < 1e-06 )
         gap = -1.0;
      else if( pb*db < 0.0 )
         gap = -1.0;
      else if( abs(db) >= +infty )
         gap = -1.0;
      else if( abs(pb) >= +infty )
         gap = -1.0;
      else
         gap = 100.0*abs((pb-db)/min(abs(db),abs(pb)));
      if( gap < 0.0 )
         gapstr = "  --  ";
      else if( gap < 1e+04 )
         gapstr = sprintf("%6.1f", gap);
      else
         gapstr = " Large";

      if( vars == 0 )
         probtype = "--";
      else if( binvars == 0 && intvars == 0 )
         probtype = "   LP";
      else if( contvars == 0 ) {
         if( intvars == 0 && implvars == 0 )
            probtype = "   BP";
         else
            probtype = "   IP";
      }
      else {
         if( intvars == 0 )
            probtype = "  MBP";
         else
            probtype = "  MIP";
      }

      if( aborted && endtime - starttime > timelimit && timelimit > 0.0 ) {
         timeout = 1;
         aborted = 0;
         tottime = endtime - starttime;
      }
      if( aborted && tottime == 0.0 )
         tottime = timelimit;
      if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

      printf("%-*s & %6d & %6d & %14.9g & %14.9g & %6s &%s%8d &%s%7.1f \\\\\n",
	     namelength, pprob, cons, vars, db, pb, gapstr, markersym, bbnodes, markersym, tottime) >TEXFILE;

      printf("%-*s  %-5s %7d %7d %7d %7d %16.9g %16.9g %6s %9d %8d %7.1f ",
	     namelength, shortprob, probtype, origcons, origvars, cons, vars, db, pb, gapstr, iters, bbnodes, tottime);

      if( aborted ) {
         printf("abort\n");
         failtime += tottime;
         fail++;
      }
      else if( solstatus[prob] == "opt" ) {
         reltol = max(mipgap, 1e-5) * max(abs(pb),1.0);
         abstol = max(absmipgap, 1e-4);

         if( ( pb-db > max(abstol,reltol) && (db-sol[prob] > reltol || sol[prob]-pb > reltol)) || ( db-pb > max(reltol,abstol) && (sol[prob]-db > reltol || pb-sol[prob] > reltol)) ) {
            printf("fail\n");
            failtime += tottime;
            fail++;
         }
         else {
            if (timeout) {
               printf("timeout\n");
               timeouttime += tottime;
               timeouts++;
            }
            else if (memout) {
               printf("memlimit\n");
               timeouttime += tottime;
               timeouts++;
            }
            else if (nodeout) {
               printf("nodelimit\n");
               timeouttime += tottime;
               timeouts++;
            }
            else {
               if( (abs(pb - db) <= max(abstol, reltol)) && abs(pb - sol[prob]) <= reltol ) {
                  printf("ok\n");
                  pass++;
               }
               else {
                  printf("fail\n");
                  failtime += tottime;
                  fail++;
               }
            }
         }
      }
      else if( solstatus[prob] == "best" ) {
         reltol = max(mipgap, 1e-5) * max(abs(pb),1.0);
         abstol = max(absmipgap, 1e-4);

         if( ( pb-db > max(abstol,reltol) && db-sol[prob] > reltol) || ( db-pb > max(reltol,abstol) && sol[prob]-db > reltol) ) {
            printf("fail\n");
            failtime += tottime;
            fail++;
         }
         else {
            if (timeout || memout || nodeout) {
               if( (pb-db > max(abstol,reltol) && sol[prob]-pb > reltol) || (db-pb > max(abstol,reltol) && pb-sol[prob] > reltol) ) {
                  printf("better\n");
                  timeouttime += tottime;
                  timeouts++;
               }
	       else if (memout) {
		  printf("memlimit\n");
		  timeouttime += tottime;
		  timeouts++;
	       }
	       else if (nodeout) {
		  printf("nodelimit\n");
		  timeouttime += tottime;
		  timeouts++;
	       }
               else {
                  printf("timeout\n");
                  timeouttime += tottime;
                  timeouts++;
               }
            }
            else {
               if( abs(pb - db) <= max(abstol, reltol) ) {
                  printf("solved not verified\n");
                  pass++;
               }
               else {
                  printf("fail\n");
                  failtime += tottime;
                  fail++;
               }
            }
         }
      }
      else if( solstatus[prob] == "unkn" ) {
         reltol = max(mipgap, 1e-5) * max(abs(pb),1.0);
         abstol = max(absmipgap, 1e-4);

         if( timeout || memout || nodeout ) {
            if( abs(pb) < infty )
               printf("better\n");
            else if( timeout )
	       printf("timeout\n");
	    else if (memout)
	       printf("memlimit\n");
	    else if (nodeout)
	       printf("nodelimit\n");
	    timeouttime += tottime;
	    timeouts++;
         }
         else if( abs(pb - db) <= max(abstol, reltol) ) {
            printf("solved not verified\n");
            pass++;
         }
	 else
	    printf("unknown\n");
      }
      else if( solstatus[prob] == "inf" ) {
	 if( timeout || memout || nodeout ) {
            if( timeout )
               printf("timeout\n");
	    else if (memout)
               printf("memlimit\n");
            else if (nodeout)
               printf("nodelimit\n");
	    timeouttime += tottime;
	    timeouts++;
	 }
	 else if( !feasible ) {
	    printf("ok\n");
	    pass++;
	 }
	 else {
	    printf("fail\n");
	    failtime += tottime;
	    fail++;
         }
      }
      else {
         reltol = max(mipgap, 1e-5) * max(abs(pb),1.0);
         abstol = max(absmipgap, 1e-4);

         if( timeout || memout || nodeout ) {
            if( timeout )
	       printf("timeout\n");
	    else if (memout)
	       printf("memlimit\n");
	    else if (nodeout)
	       printf("nodelimit\n");
	    timeouttime += tottime;
	    timeouts++;
         }
         else if( abs(pb - db) <= max(abstol, reltol) ) {
            printf("solved not verified\n");
            pass++;
         }
	 else
	    printf("unknown\n");
      }

      if( writesolufile ) {
         if( pb == +infty && db == +infty )
            printf("=inf= %-18s\n",prob)>NEWSOLUFILE;
         else if( pb == db )
            printf("=opt= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else if( pb < +infty )
            printf("=best= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else
            printf("=unkn= %-18s\n",prob)>NEWSOLUFILE;
      }

      sbab     += bbnodes;
      scut     += cuts;
      stottime += tottime;
      timegeom = timegeom^((nprobs-1)/nprobs) * max(tottime, 1.0)^(1.0/nprobs);
      nodegeom = nodegeom^((nprobs-1)/nprobs) * max(bbnodes, 1.0)^(1.0/nprobs);
      shiftedtimegeom = shiftedtimegeom^((nprobs-1)/nprobs) * max(tottime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftednodegeom = shiftednodegeom^((nprobs-1)/nprobs) * max(bbnodes+nodegeomshift, 1.0)^(1.0/nprobs);
   }
}
END {
   shiftednodegeom -= nodegeomshift;
   shiftedtimegeom -= timegeomshift;

   printf("\\midrule\n")                                                 >TEXFILE;
   printf("%-14s (%2d) &        &        &                &                &        & %9d & %8.1f \\\\\n",
	  "Total", nprobs, sbab, stottime) >TEXFILE;
   printf("%-14s      &        &        &                &                &        & %9d & %8.1f \\\\\n",
	  "Geom. Mean", nodegeom, timegeom) >TEXFILE;
   printf("%-14s      &        &        &                &                &        & %9d & %8.1f \\\\\n",
	  "Shifted Geom.", shiftednodegeom, shiftedtimegeom) >TEXFILE;
   printf("\\bottomrule\n")                                              >TEXFILE;
   printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
   printf("\\end{supertabular*}\n")                                      >TEXFILE;
   printf("\\end{center}\n")                                             >TEXFILE;
   printf("\\end{document}\n")                                           >TEXFILE;

   printf("------------------+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+-------\n");

   printf("\n");
   printf("------------------------------[Nodes]---------------[Time]------\n");
   printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom. \n");
   printf("----------------------------------------------------------------\n");
   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f\n",
	  nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom);
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f\n",
	  nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   printf("----------------------------------------------------------------\n");

   printf("@02 threads: %g\n", threads);
   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 XPRESS(%s):%s\n", version, settings);
}
