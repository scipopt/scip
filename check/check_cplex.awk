#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2006 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2006 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check_cplex.awk,v 1.17 2006/02/23 12:40:31 bzfpfend Exp $
#
#@file    check_cplex.awk
#@brief   CPLEX Check Report Generator
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

    printf("-------------------------+-------+------+--------------+--------------+------+---------------+-------+------+-------\n");
    printf("Name                     | Conss | Vars |   Dual Bound | Primal Bound | Gap% |               | Nodes | Time |       \n");
    printf("-------------------------+-------+------+--------------+--------------+------+---------------+-------+------+-------\n");

    nprobs   = 0;
    sbab     = 0;
    scut     = 0;
    stottime = 0.0;
    sgap     = 0.0;
    nodegeom = 1.0;
    timegeom = 1.0;
    failtime = 0.0;
    timeouttime = 0.0;
    fail     = 0;
    pass     = 0;
    timeouts = 0;
}
/=opt=/ { solfeasible[$2] = 1; sol[$2] = $3; }  # get optimum
/=inf=/ { solfeasible[$2] = 0; sol[$2] = 0.0; } # problem infeasible
#
# problem name
#
/^@01/ { 
    n  = split ($2, a, "/");
    split(a[n], b, ".");
    prob = b[1];
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
    pb         = 0.0;
    db         = 0.0;
    bbnodes    = 0;
    primlps    = 0;
    primiter   = 0;
    duallps    = 0;
    dualiter   = 0;
    sblps      = 0;
    sbiter     = 0;
    tottime    = 0.0;
}
#
# problem size
#
/^Reduced MIP has/ { cons = $4; vars = $6; }
#
# solution
#
/^Integer/  {
    if ($2 == "infeasible." || $2 == "infeasible")
    {
	db = 1e+20;
	pb = 1e+20;
	absgap = 0.0;
        feasible = 0;
    }
    else
    {
	db = $NF;
	pb = $NF;
	absgap = 0.0;
        feasible = 1;
    }
    opti = ($2 == "optimal") ? 1 : 0;
}
/^MIP - Integer/  { # since CPLEX 10.0
    if ($4 == "infeasible." || $4 == "infeasible")
    {
	db = 1e+20;
	pb = 1e+20;
	absgap = 0.0;
        feasible = 0;
    }
    else
    {
	db = $NF;
	pb = $NF;
	absgap = 0.0;
        feasible = 1;
    }
    opti = ($4 == "optimal") ? 1 : 0;
}
/^Solution limit exceeded, integer feasible:/ {
   pb = $8;
   timeout = 1;
}
/^MIP - Solution limit exceeded, integer feasible:/ { # since CPLEX 10.0
   pb = $10;
   timeout = 1;
}
/^Time/  {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^MIP - Time/  { # since CPLEX 10.0
    pb = ($6 == "no") ? 1e+75 : $10;
    timeout = 1;
}
/^Node/ {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^MIP - Node/ { # since CPLEX 10.0
    pb = ($6 == "no") ? 1e+75 : $10;
    timeout = 1;
}
/^Tree/  {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^MIP - Tree/  { # since CPLEX 10.0
    pb = ($6 == "no") ? 1e+75 : $10;
    timeout = 1;
}
/^Error/  {
    pb = ($3 == "no") ? 1e+75 : $7;
}
/^MIP - Error/  { # since CPLEX 10.0
    pb = ($5 == "no") ? 1e+75 : $9;
}
/cuts applied/ { 
    cuts += $NF;
}
/^Current MIP best bound =/ {
    db = $6;
    split ($9, a, ",");
    absgap = a[1];
    if (opti == 1) 
	absgap = 0.0;
}
/^Solution/ {
   tottime   = $4;
   bbnodes   = $11;
}
/^=ready=/ {
   sbab     += bbnodes;
   scut     += cuts;
   stottime += tottime;
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
   
   printf("%-26s %6d %6d %14.9g %14.9g %6s                 %7d %6.1f ",
      prob, cons, vars, db, pb, gapstr, bbnodes, tottime);
   
   if (sol[prob] == "")
      printf("unknown\n");
   else
   {
      if (solfeasible[prob])
      {
         if ((abs(pb - db) > 1e-4) || (abs(pb - sol[prob]) > 1e-6*max(abs(pb),1.0)))
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
   }
   
   if( tottime < 1.0 )
      tottime = 1.0;
   timegeom = timegeom^((nprobs-1)/nprobs) * tottime^(1.0/nprobs);
   if( bbnodes < 1 )
      bbnodes = 1;
   nodegeom = nodegeom^((nprobs-1)/nprobs) * bbnodes^(1.0/nprobs);
}
END {   
    printf("\\midrule\n")                                                 >TEXFILE;
    printf("%-14s (%2d) &        &        &                &                &        & %9d & %8.1f \\\\\n",
       "Total", nprobs, sbab, stottime) >TEXFILE;
    printf("%-14s      &        &        &                &                &        & %9d & %8.1f \\\\\n",
       "Geom. Mean", nodegeom, timegeom) >TEXFILE;
    printf("\\bottomrule\n")                                              >TEXFILE;
    printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
    printf("\\end{tabular*}\n")                                           >TEXFILE;
    printf("\\caption{CPLEX with default settings}\n")                    >TEXFILE;
    printf("\\end{center}\n")                                             >TEXFILE;
    printf("\\end{table}\n")                                              >TEXFILE;
    printf("\\end{document}\n")                                           >TEXFILE;
    
    printf("------------------+-------+------+--------------+--------------+------+---------------+-------+------+-------\n");
    
    printf("\n");
    printf("------------------------------[Nodes]---------------[Time]------\n");
    printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom. \n");
    printf("----------------------------------------------------------------\n");
    printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f\n",
       nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom);
    printf("----------------------------------------------------------------\n");
}
