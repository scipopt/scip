#!/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check_cbc.awk,v 1.7 2009/04/06 13:06:47 bzfberth Exp $
#
#@file    check_cbc.awk
#@brief   CBC Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#
function abs(x)
{
    return x < 0 ? -x : x;
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

   printf("\\documentclass[leqno]{article}\n")                      >TEXFILE;
   printf("\\usepackage{a4wide}\n")                                 >TEXFILE;
   printf("\\usepackage{amsmath,amsfonts,amssymb,booktabs}\n")      >TEXFILE;
   printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
   printf("\\begin{document}\n\n")                                  >TEXFILE;
   printf("\\begin{table}[p]\n")                                    >TEXFILE;
   printf("\\begin{center}\n")                                      >TEXFILE;
   printf("\\setlength{\\tabcolsep}{2pt}\n")                        >TEXFILE;
   printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n")    >TEXFILE;
   printf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrrrrrrr@{}}\n") >TEXFILE;
   printf("\\toprule\n")                                            >TEXFILE;
   printf("Name                &  Conss &   Vars &     Dual Bound &   Primal Bound &  Gap\\% &     Nodes &     Time \\\\\n") > TEXFILE;
   printf("\\midrule\n")                                            >TEXFILE;

   printf("------------------+-------+------+----------------+----------------+------+---------+--------+-------+-------\n");
   printf("Name              | Conss | Vars |   Dual Bound   |  Primal Bound  | Gap% |   Iters |  Nodes |  Time |       \n");
   printf("------------------+-------+------+----------------+----------------+------+---------+--------+-------+-------\n");

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
}
/=opt=/  { solstatus[$2] = "opt"; sol[$2] = $3; }  # get optimum
/=inf=/  { solstatus[$2] = "inf"; sol[$2] = 0.0; } # problem infeasible
/=best=/ { solstatus[$2] = "best"; sol[$2] = $3; } # get best known solution value
/=unkn=/ { solstatus[$2] = "unkn"; }               # no feasible solution known

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

   if( useshortnames && length(prob) > 18 )
      shortprob = substr(prob, length(prob)-17, 18);
   else
      shortprob = prob;

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];
   vars       = 0;
   cons       = 0;
   timeout    = 0;
   opti       = 0;
   feasible   = 1;
   cuts       = 0;
   pb         = +1e+75;
   db         = -1e+75;
   mipgap     = 1e-4;
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
   timelimit  = 0.0;
   cbcversion = "";
}

#
# remove Coin:
#
/^Coin:/ {
   $0 = substr($0, 6, length($0)-5);
}

#
# Cbc version
#
/^Coin Cbc and Clp/ {
   split($7, v, ",");
   cbcversion = v[1];
}

#
# settings
#
/^@03 SETTINGS:/ {
   settings = $3;
}
/^ratioGap has value/ {
   mipgap = $4;
}
/^allowableGap has value/ {
   absmipgap = $4;
}
/^seconds has value/ {
   timelimit = $4;
}

#
# problem size
#
/^Problem has/ {
   cons = $3;
   vars = $5;
}

#
# node log (dual bound)
#
#/^Cbc0010I After [0-9]* nodes, [0-9]* on tree, [-+.0-9eE]* best solution, best possible [-+.0-9eE]* ([-+.0-9eE]* seconds)/ {
/^Cbc0010I After .* nodes, .* on tree, .* best solution, best possible .* (.* seconds)/ {
   bbnodes = $3;
   pb = $8;
   db = $13;
   tottime = substr($14, 2, length($14)-1);
}

#
# solution
#
/^Result - Finished objective/ {
   pb       = $5;
   db       = pb;
   bbnodes  = $7;
   dualiter = $10;
   iters    = dualiter;
   tottime  = $14;
   feasible = (pb < 1e+20);
   opti     = 1;
   aborted  = 0;
}
/^Result - Finishedgap objective/ {
   pb       = $5;
   db       = pb;
   bbnodes  = $7;
   dualiter = $10;
   iters    = dualiter;
   tottime  = $14;
   feasible = (pb < 1e+20);
   opti     = 1;
   aborted  = 0;
}
/^Result - Stopped on time objective/ {
   pb       = $7;
   bbnodes  = $9;
   dualiter = $12;
   iters    = dualiter;
   tottime  = $16;
   feasible = (pb < 1e+20);
   opti     = 0;
   timeout  = 1;
   aborted  = 0;
}
/^Result - Stopped on nodes objective/ {
   pb       = $7;
   bbnodes  = $9;
   dualiter = $12;
   iters    = dualiter;
   tottime  = $16;
   feasible = (pb < 1e+20);
   opti     = 0;
   timeout  = 1;
   aborted  = 0;
}
/^Cputime limit exceeded/ {
   feasible = (pb < 1e+20);
   opti     = 0;
   timeout  = 1;
   aborted  = 0;
   tottime  = max(tottime, timelimit);
}

#
# evaluation
#
/^=ready=/ {
   if( !onlyinsolufile || solstatus[prob] != "" )
   {
      bbnodes = max(bbnodes, 1);

      nprobs++;
    
      optimal = 0;
      markersym = "\\g";
      if( abs(pb - db) < 1e-06 )
      {
         gap = 0.0;
         optimal = 1;
         markersym = "  ";
      }
      else if( abs(db) < 1e-06 )
         gap = -1.0;
      else if( pb*db < 0.0 )
         gap = -1.0;
      else if( abs(db) >= 1e+20 )
         gap = -1.0;
      else if( abs(pb) >= 1e+20 )
         gap = -1.0;
      else
         gap = 100.0*abs((pb-db)/db);
      if( gap < 0.0 )
         gapstr = "  --  ";
      else if( gap < 1e+04 )
         gapstr = sprintf("%6.1f", gap);
      else
         gapstr = " Large";

      printf("%-19s & %6d & %6d & %14.9g & %14.9g & %6s &%s%8d &%s%7.1f \\\\\n",
         pprob, cons, vars, db, pb, gapstr, markersym, bbnodes, markersym, tottime) >TEXFILE;

      printf("%-19s %6d %6d %16.9g %16.9g %6s %9d %8d %7.1f ",
         shortprob, cons, vars, db, pb, gapstr, iters, bbnodes, tottime);

      if( aborted )
      {
         printf("abort\n");
         failtime += tottime;
         fail++;
      }
      else if( solstatus[prob] == "opt" )
      {
 	 reltol = max(mipgap, 1e-5) * max(abs(pb),1.0);
	 abstol = max(absmipgap, 1e-4);
         if( (abs(pb - db) > max(abstol, reltol)) || (abs(pb - sol[prob]) > reltol) )
         {
            if (timeout)
            {
               printf("timeout\n");
               timeouttime += tottime;
               timeouts++;
            }
            else
            {
               printf("fail\n");
               failtime += tottime;
               fail++;
            }
         }
         else
         {
            printf("ok\n");
            pass++;
         }
      }
      else if( solstatus[prob] == "inf" )
      {
         if (feasible)
         {
            if (timeout)
            {
               printf("timeout\n");
               timeouttime += tottime;
               timeouts++;
            }
            else
            {
               printf("fail\n");
               failtime += tottime;
               fail++;
            }
         }
         else
         {
            printf("ok\n");
            pass++;
         }
      }
      else
      {
         if (timeout)
         {
            printf("timeout\n");
            timeouttime += tottime;
            timeouts++;
         }
         else
            printf("unknown\n");
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
   printf("\\end{tabular*}\n")                                           >TEXFILE;
   printf("\\caption{CPLEX with default settings}\n")                    >TEXFILE;
   printf("\\end{center}\n")                                             >TEXFILE;
   printf("\\end{table}\n")                                              >TEXFILE;
   printf("\\end{document}\n")                                           >TEXFILE;

   printf("------------------+-------+------+----------------+----------------+------+---------+--------+-------+-------\n");

   printf("\n");
   printf("------------------------------[Nodes]---------------[Time]------\n");
   printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom. \n");
   printf("----------------------------------------------------------------\n");
   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f\n",
      nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom);
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f\n",
      nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   printf("----------------------------------------------------------------\n");

   printf("@01 CBC(%s):%s\n", cbcversion, settings);
}
