function extract_name(fullname,  a,b,c,n,pprob) {
    n = split (fullname, a, "/");
    n = split (a[n], b, "'");
    split(b[n-1], c, ".");
    n = split (c[1], a, "_");
    if (n == 1) {
	pprob = sprintf ("%s", a[1]);
    }
    else {
	pprob = sprintf ("%s\\_%s", a[1], a[2]);
    }
    return pprob;
}
#
# head line
#
BEGIN {
    printf("\\documentclass[leqno]{article}\n");
    printf("\\usepackage{a4wide}\n");
    printf("\\usepackage{amsmath,amsfonts,amssymb}\n");
    printf("\\font\\alex=cmr10 at 9truept\n");
    printf("\\begin{document}\n\n");
    printf("\\begin{table}[p]\n");
    printf("\\renewcommand{\\arraystretch}{0.818}\n");
    printf("\\begin{alex}\n");
    printf("\\begin{center}\n");
    printf("\\begin{tabular}{l@{\\quad\\enspace}rrrrrr}\n");
    printf("\\hline\n");;
    printf("\\noalign{\\smallskip}\n");
    printf("Example & B \\& B & Cuts & Dual Bound & Primal Bound & Time & Gap \\% \\\\\n");;
    printf("\\hline\n");
    printf("\\noalign{\\smallskip}\n");
    total    = 0;
    sbab     = 0;
    scut     = 0;
    spb      = 0;
    sdb      = 0;
    stottime = 0.0;
    sgap     = 0.0;
}
#
# problem name
#
/^Problem/ { 
    pprob = extract_name($2);
    opti = 0;
    cuts = 0;
}
/^CPLEX> Problem/ { 
    pprob = extract_name($3);
    opti = 0;
    cuts = 0;
}
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
}
/^Node/ {
    pb = ($4 == "no") ? 1e+75 : $8;
}
/^Tree/  {
    pb = ($4 == "no") ? 1e+75 : $8;
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
    babnodes  = $11;
    sbab     += babnodes;
    scut     += cuts;
    stottime += tottime;
    sgap     += gap;
    total++;
    
    if (pb > 1e+70  ||  pb < -1e+70) {
	printf ("%-12s & %7d & %5d & %14.10g & %14s & %6.1f & %7s \\\\\n", pprob, babnodes, cuts, db, "       -      ", tottime, "   -   ");
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
	printf ("%-12s & %7d & %5d & %14.10g & %14.10g & %6.1f & %7.3f \\\\\n", pprob, babnodes, cuts, db, pb, tottime, gap);
    }
}
END {   
    printf("\\hline\n");
    printf("\\noalign{\\vspace{1.5pt}}\n");
    printf ("%-12s (%d) & %8d & %5d & & & %8.1f & %7.3f \\\\\n", "Total", total, sbab, scut, stottime, sgap);
    printf("\\hline\n");
    printf("\\end{tabular}\n");
    printf("\\caption{{\\tt CPLEX} with default parameter settings}\n");
    printf("\\end{center}\n");
    printf("\\end{alex}\n");
    printf("\\end{table}\n");
    printf("\\end{document}");
}



