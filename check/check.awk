# $Id: check.awk,v 1.1 2002/10/23 14:31:35 bzfpfend Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: check.awk                                                     *
#*   Name....: SIP Check Report Generator                                    *
#*   Author..: Alexander Martin, Thorsten Koch, Tobias Pfender               *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# head line
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
    printf("Example & B \\& B & Cuts & Dual Bound & ")               >TEXFILE;
    printf("Primal Bound & Time & Gap \\% \\\\\n")                   >TEXFILE;
    printf("\\hline\n")                                              >TEXFILE;
    printf("\\noalign{\\smallskip}\n")                               >TEXFILE;

    total    = 0;
    sbab     = 0;
    slp      = 0;
    ssim     = 0;
    scut     = 0;
    stottime = 0.0;
    sgap     = 0.0;
    brrule   = "default";

    printf("------------------+-------+------+-------+--------------+------+------+-------\n");
    printf("Name              |  Rows | Cols | Nodes |   Upperbound |  Gap | Time |\n");
    printf("------------------+-------+------+-------+--------------+------+------+-------\n");
}
/=opt=/ { sol[$2] = $3; }  # get optimum
/New variable selection:/ { brrule = $4 " " $5 " " $6; } 
#/^Problem/ { 
#    if ($2 != "infeasible") {
#	n  = split($2, a, "/");
#	prob = a[n];
#	n = split(prob, a, "_");
#	if (n == 1) 
#	    pprob = prob;
#	else 
#	    pprob = a[1] "\\_" a[2];
#	
#	separation = 0;
#	cut        = 0;
#    }
#}
/^\(SIP\) Problem/ { 
    n  = split ($3, a, "/");
    split(a[n], b, ".");
    prob = b[1];
    # Escape _ for TeX
    n = split(prob, a, "_");
    pprob = a[1];
    for( i = 2; i <= n; i++ )
       pprob = pprob "\\_" a[i];
#    if (n == 1) 
#	pprob = prob;
#    else 
#	pprob = a[1] "\\_" a[2];

    cols       = $5;
    rows       = $13;
    nzos       = $15;    
    separation = 0;
    cut        = 0;
}
/^Problem/ { 
    n  = split ($2, a, "/");
    split(a[n], b, ".");
    prob = b[1];
    # Escape _ for TeX
    n = split(prob, a, "_");
    pprob = a[1];
    for( i = 2; i <= n; i++ )
       pprob = pprob "\\_" a[i];
#    if (n == 1) 
#	pprob = prob;
#    else 
#	pprob = a[1] "\\_" a[2];

    cols       = $4;
    rows       = $12;
    nzos       = $14;    
    separation = 0;
    cut        = 0;
}
/time limit reached/ { timeout = 1; }
/node limit reached/ { timeout = 1; }
#
# solution
#
/^   Optimum Value  :/ { ub = db = $4; }
/^   Primal Bound   :/ { ub = $4; }
/^   Dual Bound     :/ { db = $4; }
#
# iterations
#
/^   B & B  nodes   :/ { bbnodes = $6 }
/^   LPs            :/ { lps = $3     }
/^   Simplex        :/ { simplex = $3 }
#
# separation
#
/^Separation/          { separation = 1 }
/^   Gomory's MIR in:/ { if (separation == 1 && $4 != "off")  cut += $4 }
/^   c-MIR inequalit:/ { if (separation == 1 && $3 != "off")  cut += $3 }
/^   generalized kna:/ { if (separation == 1 && $3 != "off")  cut += $3 }
/^   knapsack       :/ { if (separation == 1 && $3 != "off")  cut += $3 }
/^   SOS knapsack   :/ { if (separation == 1 && $4 != "off")  cut += $4 }
/^   aggregate knaps:/ { if (separation == 1 && $3 != "off")  cut += $3 }
/^   pool separation:/ { if (separation == 1 && $3 != "off")  cut += $3 }
#
# time
#
/^Time/                { separation = 0 }
/^   Total          :/ { 
    tottime = $3;
    sbab += bbnodes;
    scut += cut;
    slp  += lps;
    ssim += simplex;
    total++;
    stottime += tottime;
    if (ub > 1e+19  ||  ub < -1e+19) {
	printf ("%-19s & %7d & %5d & %14.10g & %14s & %6.1f & %s \\\\\n", 
		pprob, bbnodes, cut, db, 
		"\\multicolumn{1}{c}{   -   }", tottime, 
		"\\multicolumn{1}{c}{   -   }")                      >TEXFILE;
    }
    else {
	if ( ub > db ) {
	    if (db >  0.001)  
		gap = 100.0 * (ub - db) / (1.0 * db);
	    else if (db < -0.001)  
		gap = 100.0 * (ub - db) / (-1.0 * db);
	    else      
		gap = 0.0;
	}
	else {
	    if (db >  0.001)  
		gap = 100.0 * (db - ub) / (1.0 * db);
	    else if (db < -0.001)  
		gap = 100.0 * (db - ub) / (-1.0 * db);
	    else                   
		gap = 0.0;
	}
	sgap += gap;

	printf ("%-19s & %7d & %5d & %14.10g & %14.10g & %6.1f & %7.3f \\\\\n",
		pprob, bbnodes, cut, db, ub, tottime, gap)           >TEXFILE;
    }
    printf("%-19s %6d %6d %7d %14.10g %6.1f %6.1f ",
	   prob, rows, cols, bbnodes, ub, gap, tottime);
    
    probs++;
    sumtime  += tottime;
    sumnodes += bbnodes;
    
    if (sol[prob] == "")
      printf("unknown\n");
    else {
      if ((abs(ub - db) > 1e-4) || (abs(ub - sol[prob])/max(abs(ub),1.0) > 1e-6)) {
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
}
END {   
    printf("\\hline\n")                                              >TEXFILE;
    printf("\\noalign{\\vspace{1.5pt}}\n")                           >TEXFILE;
    printf ("%-19s (%d) & %8d & %5d & & & %8.1f & %7.3f \\\\\n", 
	    "Total", total, sbab, scut, stottime, sgap)              >TEXFILE;
    printf("\\hline\n")                                              >TEXFILE;
    printf("\\end{tabular}\n")                                       >TEXFILE;
    printf("\\caption{%s}\n", brrule)                                >TEXFILE;
    printf("\\end{center}\n")                                        >TEXFILE;
    printf("\\end{alex}\n")                                          >TEXFILE;
    printf("\\end{table}\n")                                         >TEXFILE;
    printf("\\end{document}")                                        >TEXFILE;

    printf("------------------+-------+------+-------+--------------+------+------+-------\n");

    printf("\n-------------------------------------------\n");
    printf("  Cnt  Pass  Fail  kNodes FailTime  TotTime\n");
    printf("-------------------------------------------\n");
    printf("%5d %5d %5d %7d %8.0f %8.0f\n",
       probs, pass, fail, sumnodes / 1000, failtime, sumtime);
    printf("-------------------------------------------\n");
}






