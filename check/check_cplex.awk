#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2004 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2004 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the SCIP Academic Licence.        *
#*                                                                           *
#*  You should have received a copy of the SCIP Academic License             *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check_cplex.awk,v 1.3 2004/09/28 10:50:29 bzfpfend Exp $
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
    printf("\\documentclass[leqno]{article}\n") >TEXFILE;
    printf("\\usepackage{a4wide}\n") >TEXFILE;
    printf("\\usepackage{amsmath,amsfonts,amssymb}\n") >TEXFILE;
    printf("\\font\\alex=cmr10 at 9truept\n") >TEXFILE;
    printf("\\begin{document}\n\n") >TEXFILE;
    printf("\\begin{table}[p]\n") >TEXFILE;
    printf("\\renewcommand{\\arraystretch}{0.818}\n") >TEXFILE;
    printf("\\begin{alex}\n") >TEXFILE;
    printf("\\begin{center}\n") >TEXFILE;
    printf("\\begin{tabular}{l@{\\quad\\enspace}rrrrrr}\n") >TEXFILE;
    printf("\\hline\n") >TEXFILE;
    printf("\\noalign{\\smallskip}\n") >TEXFILE;
    printf("Example & B \\& B & Cuts & Dual Bound & Primal Bound & Time & Gap \\% \\\\\n") >TEXFILE;
    printf("\\hline\n") >TEXFILE;
    printf("\\noalign{\\smallskip}\n") >TEXFILE;
    nprobs   = 0;
    sbab     = 0;
    scut     = 0;
    stottime = 0.0;
    sgap     = 0.0;
    nodegeom = 1.0;
    timegeom = 1.0;
    failtime = 0.0;
    fail     = 0;
    pass     = 0;

    printf("------------------+-------+------+-------+--------------+--------------+------+------+-------\n");
    printf("Name              | Conss | Vars | Nodes |   Dual Bound | Primal Bound |  Gap | Time |\n");
    printf("------------------+-------+------+-------+--------------+--------------+------+------+-------\n");
}
/=opt=/ { sol[$2] = $3; }  # get optimum
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
    if ($2 == "infeasible.") {
	db = 1e+20;
	pb = 1e+20;
	absgap = 0.0;
    }
    else {
	db = $NF;
	pb = $NF;
	absgap = 0.0;
    }
    opti = ($2 == "optimal") ? 1 : 0;
}
/^Time/  {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^Node/ {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^Tree/  {
    pb = ($4 == "no") ? 1e+75 : $8;
    timeout = 1;
}
/^Error/  {
    pb = ($3 == "no") ? 1e+75 : $7;
}
/cuts applied/ { 
    cuts += $NF;
}
/^Current MIP best bound/ {
    db = $6;
    split ($9, a, ",");
    absgap = a[1];
    if (opti == 1) 
	absgap = 0.0;
}
/^Solution/ {
    tottime   = $4;
    bbnodes   = $11;
    sbab     += bbnodes;
    scut     += cuts;
    stottime += tottime;
    nprobs++;
    
    if (pb > 1e+70  ||  pb < -1e+70) {
	printf ("%-12s & %7d & %5d & %14.10g & %14s & %6.1f & %7s \\\\\n", pprob, bbnodes, cuts, db, "       -      ", tottime, "   -   ") >TEXFILE;
    }
    else {
	if (absgap < 0.0000001 && absgap > -0.0000001) 
	    gap = 0.0;
	else if ( pb > db ) {
	    if  (db >  0.001)  
		gap = 100.0 * (pb - db) / (1.0 * db);
	    else if (db < -0.001)  
		gap = 100.0 * (pb - db) / (-1.0 * db);
	    else
		gap = 0.0;
	}
	else {
	    if (db >  0.001)  
		gap = 100.0 * (db - pb) / (1.0 * db);
	    else if (db < -0.001)  
		gap = 100.0 * (db - pb) / (-1.0 * db);
	    else
		gap = 0.0;
	}
        sgap     += gap;

	printf ("%-12s & %7d & %5d & %14.10g & %14.10g & %6.1f & %7.3f \\\\\n", pprob, bbnodes, cuts, db, pb, tottime, gap) >TEXFILE;
    }
    printf("%-19s %6d %6d %7d %14.9g %14.9g %6.1f %6.1f ",
       prob, cons, vars, bbnodes, db, pb, gap, tottime);

   if (sol[prob] == "")
      printf("unknown\n");
   else {
      if ((abs(pb - db) > 1e-4) || (abs(pb - sol[prob]) > 1e-6*max(abs(pb),1.0))) {
	 if (timeout)
	    printf("timeout\n");
	 else
	    printf("fail\n");
	 failtime += tottime;
	 fail++;
      }
      else {
	 printf("ok\n");
	 pass++;
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
    printf("\\hline\n") >TEXFILE;
    printf("\\noalign{\\vspace{1.5pt}}\n") >TEXFILE;
    printf ("%-12s (%d) & %8d & %5d & & & %8.1f & %7.3f \\\\\n", "Total", nprobs, sbab, scut, stottime, sgap) >TEXFILE;
    printf("\\hline\n") >TEXFILE;
    printf("\\end{tabular}\n") >TEXFILE;
    printf("\\caption{{\\tt CPLEX} with default parameter settings}\n") >TEXFILE;
    printf("\\end{center}\n") >TEXFILE;
    printf("\\end{alex}\n") >TEXFILE;
    printf("\\end{table}\n") >TEXFILE;
    printf("\\end{document}") >TEXFILE;

    printf("------------------+-------+------+-------+--------------+--------------+------+------+-------\n");
    
    printf("\n----------------------------------------------------------------\n");
    printf("  Cnt  Pass  Fail  kNodes FailTime  TotTime  NodeGeom  TimeGeom\n");
    printf("----------------------------------------------------------------\n");
    printf("%5d %5d %5d %7d %8.0f %8.0f %9.1f %9.1f\n",
	   nprobs, pass, fail, sbab / 1000, failtime, stottime, nodegeom, timegeom);
    printf("----------------------------------------------------------------\n");
}
