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
# $Id: check.awk,v 1.29 2006/01/03 12:22:37 bzfpfend Exp $
#
#@file    check.awk
#@brief   SCIP Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Alexander Martin
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
    printf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrrrrrrr@{}}\n")  >TEXFILE;
#    printf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrr@{}}\n")  >TEXFILE;
    printf("\\toprule\n")                                            >TEXFILE;
    printf("Name                &  Conss &   Vars &     Dual Bound &   Primal Bound &  Gap\\% &    Confs &    Lits &     Nodes &     Time &   BTime &   OTime &   CTime \\\\\n") > TEXFILE;
#    printf("Name                &  Conss &   Vars &     Dual Bound &   Primal Bound &  Gap\\% &     Nodes &     Time \\\\\n") > TEXFILE;
    printf("\\midrule\n")                                            >TEXFILE;

    printf("------------------+-------+------+--------------+--------------+------+-------+-------+-------+------+------+------+------+-------\n");
    printf("Name              | Conss | Vars |   Dual Bound | Primal Bound | Gap% | Confs |  Lits | Nodes | Time | BTim | OTim | CTim |       \n");
    printf("------------------+-------+------+--------------+--------------+------+-------+-------+-------+------+------+------+------+-------\n");
#    printf("------------------+-------+------+--------------+--------------+------+-------+------+-------\n");
#    printf("Name              | Conss | Vars |   Dual Bound | Primal Bound | Gap% | Nodes | Time |       \n");
#    printf("------------------+-------+------+--------------+--------------+------+-------+------+-------\n");

    nprobs = 0;
    sbab = 0;
    slp = 0;
    ssim = 0;
    ssblp = 0;
    stottime = 0.0;
    nodegeom = 1.0;
    timegeom = 1.0;
    sblpgeom = 1.0;
    timeouttime = 0.0;
    timeouts = 0;
    failtime = 0.0;
    fail = 0;
    pass = 0;
    settings = "default";
    conftottime = 0.0;
    conftimegeom = 1.0;
    overheadtottime = 0.0;
    overheadtimegeom = 1.0;
    basictimegeom = 1.0;
}
/=opt=/ { solfeasible[$2] = 1; sol[$2] = $3; }  # get optimum
/=inf=/ { solfeasible[$2] = 0; sol[$2] = 0.0; } # problem infeasible
/^@01/ { 
    n  = split ($2, a, "/");
    split(a[n], b, ".");
    prob = b[1];
    if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
       m--;
    for( i = 2; i < m; ++i )
       prob = prob "." b[i];

    if( length(prob) > 18 )
       shortprob = substr(prob, length(prob)-17, 18);
    else
       shortprob = prob;

    # Escape _ for TeX
    n = split(prob, a, "_");
    pprob = a[1];
    for( i = 2; i <= n; i++ )
       pprob = pprob "\\_" a[i];
    vars = 0;
    cons = 0;
    timeout = 0;
    feasible = 1;
    pb = +1e20;
    db = -1e20;
    bbnodes = 0;
    primlps = 0;
    primiter = 0;
    duallps = 0;
    dualiter = 0;
    sblps = 0;
    sbiter = 0;
    tottime = 0.0;
    inconflict = 0;
    inconstime = 0;
    confclauses = 0;
    confliterals = 0.0;
    conftime = 0.0;
    overheadtime = 0.0;
    aborted = 1;
}
/^SCIP> loaded parameter file/ { settings = $5; sub(/<settings\//, "", settings); sub(/.set>/, "", settings); }
#
# conflict analysis
#
/^Conflict Analysis  :/ { inconflict = 1; }
/^  propagation      :/ {
   if( inconflict == 1 )
   {
      conftime += $3; #confclauses += $5 + $7; confliterals += $5 * $6 + $7 * $8;
   }
}
/^  infeasible LP    :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  strong branching :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  pseudo solution  :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  applied globally :/ {
   if( inconflict == 1 )
   {
      confclauses += $6; confliterals += $6 * $7;
   }
}
/^  applied locally  :/ {
   if( inconflict == 1 )
   {
      confclauses += $6; confliterals += $6 * $7;
   }
}
/^Separators         :/ { inconflict = 0; }
/^Constraint Timings :/ { inconstime = 1; }
#/^  logicor          :/ { if( inconstime == 1 ) { overheadtime += $3; } }
/^Propagators        :/ { inconstime = 0; }
/^  switching time   :/ { overheadtime += $4; }
#
# problem size
#
/^  Variables        :/ { vars = $3; }
/^  Constraints      :/ { cons = $3; }
#
# solution
#
/^SCIP Status        :/ { aborted = 0; }
/solving was interrupted/  { timeout = 1; }
/problem is solved/    { timeout = 0; }
/^  Primal Bound     :/ {
   if( $4 == "infeasible" )
      feasible = 0;
   else if( $4 != "-" )
      pb = $4;
}
/^  Dual Bound       :/ { if( $4 != "-" ) db = $4; }
#
# iterations
#
/^  nodes            :/ { bbnodes = $3 }
/^  primal LP        :/ { primlps = $5; primiter = $6; }
/^  dual LP          :/ { duallps = $5; dualiter = $6; }
/^  strong branching :/ { sblps = $5; sbiter = $6; }
#
# time
#
/^Solving Time       :/ { tottime = $4 }
#
# Output
#
/^=ready=/ {
   lps = primlps + duallps;
   simplex = primiter + dualiter;
   stottime += tottime;
   sbab += bbnodes;
   slp += lps;
   ssim += simplex;
   ssblp += sblps;
   conftottime += conftime;
   overheadtottime += overheadtime;
   nprobs++;

   optimal = 0;
   markersym = "\\g";
   if( abs(pb - db) < 1e-06 || !feasible )
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

   printf("%-19s & %6d & %6d & %14.9g & %14.9g & %6s & %8d & %7.1f &%s%8d &%s%7.1f & %7.1f & %7.1f & %7.1f \\\\\n",
      pprob, cons, vars, db, pb, gapstr, confclauses, (confclauses > 0 ? confliterals / confclauses : 0.0), 
      markersym, bbnodes, markersym, tottime, tottime - conftime - overheadtime, overheadtime, conftime) >TEXFILE;
#   printf("%-19s & %6d & %6d & %14.9g & %14.9g & %6s &%s%8d &%s%7.1f \\\\\n",
#      pprob, cons, vars, db, pb, gapstr, markersym, bbnodes, markersym, tottime) >TEXFILE;
   
   printf("%-19s %6d %6d %14.9g %14.9g %6s %7d %7.1f %7d %6.1f %6.1f %6.1f %6.1f ",
      shortprob, cons, vars, db, pb, gapstr, confclauses, (confclauses > 0 ? confliterals / confclauses : 0.0), 
      bbnodes, tottime, tottime - conftime - overheadtime, overheadtime, conftime);
#   printf("%-19s %6d %6d %14.9g %14.9g %6s %7d %6.1f ",
#      shortprob, cons, vars, db, pb, gapstr, bbnodes, tottime);
   
   if( aborted )
   {
      printf("abort\n");
      failtime += tottime;
      fail++;
   }
   else if (sol[prob] == "")
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
   else
   {
      if( solfeasible[prob] )
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
   
   basictime = tottime - conftime - overheadtime;
   if( tottime < 1.0 )
      tottime = 1.0;
   timegeom = timegeom^((nprobs-1)/nprobs) * tottime^(1.0/nprobs);
   if( bbnodes < 1 )
      bbnodes = 1;
   nodegeom = nodegeom^((nprobs-1)/nprobs) * bbnodes^(1.0/nprobs);
   if( sblps < 1 )
      sblps = 1;
   sblpgeom = sblpgeom^((nprobs-1)/nprobs) * sblps^(1.0/nprobs);
   if( conftime < 1.0 )
      conftime = 1.0;
   conftimegeom = conftimegeom^((nprobs-1)/nprobs) * conftime^(1.0/nprobs);
   if( overheadtime < 1.0 )
      overheadtime = 1.0;
   overheadtimegeom = overheadtimegeom^((nprobs-1)/nprobs) * overheadtime^(1.0/nprobs);
   if( basictime < 1.0 )
      basictime = 1.0;
   basictimegeom = basictimegeom^((nprobs-1)/nprobs) * basictime^(1.0/nprobs);
}
END {   
    printf("\\midrule\n")                                                 >TEXFILE;
    printf("%-14s (%2d) &        &        &                &                &        &          &         & %9d & %8.1f & %7.1f & %7.1f & %7.1f \\\\\n",
       "Total", nprobs, sbab, stottime, stottime - conftottime - overheadtottime, overheadtottime, conftottime) >TEXFILE;
    printf("%-14s      &        &        &                &                &        &          &         & %9d & %8.1f & %7.1f & %7.1f & %7.1f \\\\\n",
       "Geom. Mean", nodegeom, timegeom, basictimegeom, overheadtimegeom, conftimegeom) >TEXFILE;
#    printf("%-14s (%2d) &        &        &                &                &        & %9d & %8.1f \\\\\n",
#       "Total", nprobs, sbab, stottime) >TEXFILE;
#    printf("%-14s      &        &        &                &                &        & %9d & %8.1f \\\\\n",
#       "Geom. Mean", nodegeom, timegeom) >TEXFILE;
    printf("\\bottomrule\n")                                              >TEXFILE;
    printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
    printf("\\end{tabular*}\n")                                           >TEXFILE;
    printf("\\caption{%s}\n", settings)                                   >TEXFILE;
    printf("\\end{center}\n")                                             >TEXFILE;
    printf("\\end{table}\n")                                              >TEXFILE;
    printf("\\end{document}\n")                                           >TEXFILE;
    
    printf("------------------+-------+------+--------------+--------------+------+-------+-------+-------+------+------+------+------+-------\n");
#    printf("------------------+-------+------+--------------+--------------+------+-------+------+-------\n");
    
    printf("\n");
    printf("------------------------------[Nodes]---------------[Time]-----------[Basic Time]-------[Overhead Time]-----[Conflict Time]-\n");
    printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom.     total     geom.     total     geom.     total     geom. \n");
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
    printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f\n",
       nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom, 
       stottime - conftottime - overheadtottime, basictimegeom,
       overheadtottime, overheadtimegeom, conftottime, conftimegeom);
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
}
