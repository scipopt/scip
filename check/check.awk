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
# $Id: check.awk,v 1.9 2004/09/02 09:19:41 bzfpfend Exp $
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
    printf("\\usepackage{amsmath,amsfonts,amssymb}\n")               >TEXFILE;
    printf("\\font\\alex=cmr10 at 9truept\n")                        >TEXFILE;
    printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
    printf("\\begin{document}\n\n")                                  >TEXFILE;
    printf("\\begin{table}[p]\n")                                    >TEXFILE;
    printf("\\renewcommand{\\arraystretch}{0.818}\n")                >TEXFILE;
    printf("\\begin{alex}\n")                                        >TEXFILE;
    printf("\\begin{center}\n")                                      >TEXFILE;
    printf("\\begin{tabular}{l@{\\quad\\enspace}rrrrrr}\n")          >TEXFILE;
    printf("\\hline\n")                                              >TEXFILE;
    printf("\\noalign{\\smallskip}\n")                               >TEXFILE;
    printf("Example & B \\& B & StrBr & Dual Bound & ")              >TEXFILE;
    printf("Primal Bound & Time & Gap \\% \\\\\n")                   >TEXFILE;
    printf("\\hline\n")                                              >TEXFILE;
    printf("\\noalign{\\smallskip}\n")                               >TEXFILE;

    nprobs   = 0;
    sbab     = 0;
    slp      = 0;
    ssim     = 0;
    ssblp    = 0;
    stottime = 0.0;
    sgap     = 0.0;
    nodegeom = 1.0;
    timegeom = 1.0;
    sblpgeom = 1.0;
    failtime = 0.0;
    fail     = 0;
    pass     = 0;
    settings = "default";

    printf("------------------+-------+------+-------+--------------+------+------+-------\n");
    printf("Name              | Conss | Vars | Nodes |   Upperbound |  Gap | Time |\n");
    printf("------------------+-------+------+-------+--------------+------+------+-------\n");
}
/=opt=/ { sol[$2] = $3; }  # get optimum
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
/^SCIP> loaded parameter file:/ { settings = $5; }
#
# problem size
#
/^  Variables        :/ { vars = $3; }
/^  Constraints      :/ { cons = $3; }
#
# solution
#
/solving was interrupted/  { timeout = 1; }
/problem is solved/    { timeout = 0; }
/^  Primal Bound     :/ { pb = $4; }
/^  Dual Bound       :/ { db = $4; }
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
   lps       = primlps + duallps;
   simplex   = primiter + dualiter;
   stottime += tottime;
   sbab     += bbnodes;
   slp      += lps;
   ssim     += simplex;
   ssblp    += sblps;
   nprobs++;

   if (pb > 1e+19  ||  pb < -1e+19) {
      printf ("%-19s & %7d & %5d & %14.9g & %14s & %6.1f & %s \\\\\n", 
	      pprob, bbnodes, sblps, db, 
	      "\\multicolumn{1}{c}{   -   }", tottime, 
	      "\\multicolumn{1}{c}{   -   }")                      >TEXFILE;
   }
   else {
      if (abs(db) > 1e-4)
      {
	 gap = abs((pb-db)/db);
      }
      else
      {
	 gap = 0.0;
      }
      sgap += gap;
      
      printf ("%-19s & %7d & %5d & %14.9g & %14.9g & %6.1f & %7.3f \\\\\n",
	      pprob, bbnodes, sblps, db, pb, tottime, gap)           >TEXFILE;
   }
   printf("%-19s %6d %6d %7d %14.9g %6.1f %6.1f ",
	  prob, cons, vars, bbnodes, pb, gap, tottime);
   
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
   if( sblps < 1 )
      sblps = 1;
   sblpgeom = sblpgeom^((nprobs-1)/nprobs) * sblps^(1.0/nprobs);
}
END {   
    printf("\\hline\n")                                                   >TEXFILE;
    printf("\\noalign{\\vspace{1.5pt}}\n")                                >TEXFILE;
    printf ("%-14s (%d) & %8d & %5d &                &                & %8.1f & %7.3f \\\\\n", 
	    "Total", nprobs, sbab, ssblp, stottime, sgap)                 >TEXFILE;
    printf ("%-14s      & %8.1f & %5.1f &                &                & %8.1f &         \\\\\n",
            "Geom.\\ Mean", nodegeom, sblpgeom, timegeom)                 >TEXFILE;
    printf("\\hline\n")                                                   >TEXFILE;
    printf("\\end{tabular}\n")                                            >TEXFILE;
    printf("\\caption{%s}\n", settings)                                   >TEXFILE;
    printf("\\end{center}\n")                                             >TEXFILE;
    printf("\\end{alex}\n")                                               >TEXFILE;
    printf("\\end{table}\n")                                              >TEXFILE;
    printf("\\end{document}")                                             >TEXFILE;
    
    printf("------------------+-------+------+-------+--------------+------+------+-------\n");
    
    printf("\n----------------------------------------------------------------\n");
    printf("  Cnt  Pass  Fail  kNodes FailTime  TotTime  NodeGeom  TimeGeom\n");
    printf("----------------------------------------------------------------\n");
    printf("%5d %5d %5d %7d %8.0f %8.0f %9.1f %9.1f\n",
	   nprobs, pass, fail, sbab / 1000, failtime, stottime, nodegeom, timegeom);
    printf("----------------------------------------------------------------\n");
}
