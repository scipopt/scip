#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
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
#@file    cmpres.awk
#@brief   SCIP Check Comparison Report Generator
#@author  Tobias Achterberg
#@author  Robert Waniek
#@author  Marc Pfetsch
#@author  Timo Berthold
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

function ceil(x)
{
   return (x == int(x) ? x : (x < 0 ? int(x) : int(x+1)));
}

function floor(x)
{
   return (x == int(x) ? x : (x < 0 ? int(x-1) : int(x)));
}

function fracceil(x,f)
{
   return ceil(x/f)*f;
}

function fracfloor(x,f)
{
   return floor(x/f)*f;
}

function mod(x,m)
{
   return (x - m*floor(x/m));
}

function printhline(nsolver,short, printsoltimes)
{
   for( s = 0; s < nsolver; ++s )
   {
      
      if( s == 0 )
         printf("--------------------+-+---------+--------+");
      else
      {
         if( !short )
            printf("-+---------+--------+------+------+");
         else
            printf("-+---------+--------+");
      }
      if( printsoltimes )
      {
         if( s == 0 )
            printf("---------+--------+");
         else
            printf("------+------+");
      }
   }
   printf("-------------\n");
}

function isfaster(t,reft,tol)
{
   return (t < 1.0/tol*reft && t <= reft - 0.2);
}

function isslower(t,reft,tol)
{
   return isfaster(reft, t, tol);
}

function texcompstr(val,refval, x,s,t)
{
   x = floor(100*(val/refval-1.0)+0.5);
   s = "";
   t = "";
   if( x < 0 )
   {
      if( x <= -texcolorlimit )
      {
         s = "\\textcolor{red}{\\raisebox{0.25ex}{\\tiny $-$}";
         t = "}";
      }
      else
         s = "\\raisebox{0.25ex}{\\tiny $-$}";
   }
   else if( x > 0 )
   {
      if( x >= +texcolorlimit )
      {
         s = "\\textcolor{blue}{\\raisebox{0.25ex}{\\tiny $+$}";
         t = "}";
      }
      else
         s = "\\raisebox{0.25ex}{\\tiny $+$}";
   }

   return sprintf("%s%d%s", s, abs(x), t);
}

function texstring(s, ts)
{
   ts = s;
   gsub(/_/, "\\_", ts);

   return ts;
}

function texint(x, ts,r)
{
   ts = "";
   x = floor(x);
   while( x != 0 )
   {
      r = mod(x, 1000);
      x = floor(x/1000);
      if( ts != "" )
         ts = "\\," ts;

      if ( x > 0 )
         ts = sprintf("%03d", r) "" ts;
      else
         ts = r "" ts;
   }

   return ts;
}

function texsolvername(s, sname)
{
   sname = solvername[s];
   if( setname[sname] != "" )
      sname = setname[sname];
   else
   {
      sub(/.*:/, "", sname);
      sub(/.*_/, "", sname);
      if( length(sname) > 12 )
         sname = substr(sname, length(sname)-11, 12);
   }

   return sname;
}

BEGIN {

   short = 0;  #for each non reference solver, only absolute time and number of nodes are printed 
   printsoltimes = 0; # for reference solver, absolute time to first and best solution are printed, for other solvers the corresponding ratios
                      #! please NOTE that this additional output is currently only available for SCIP .res-files created with the evalcheck.sh script and
                      #  the flag printsoltimes = 1 set in check.awk. If other solvers are involved, leave this flag set to 0.
   printgap = 0; # if timeout, then print absolute gap at termination in time column, if gap is finite
   printsoltimes = !short && printsoltimes; # short deactivates the detailed solution times output
   infinity = 1e+20;
   timegeomshift = 10.0;
   nodegeomshift = 100.0;
   mintime = 0.5;
   wintolerance = 1.1;
   markbettertime = 1.1;
   markworsetime  = 1.1;
   markbetternodes = 5.0;
   markworsenodes  = 5.0;
   onlymarked = 0;
   onlyprocessed = 0;
   maxscore = 10.0;
   consistency = 1;
   onlyfeasible = 0;
   onlyinfeasible = 0;
   onlyfail = 0;
   exclude = "";
   texfile = "";
   texincfile = "";
   texsummaryfile = "";
   texsummaryheader = 0;
   texsummaryweight = 0;
   texsummaryshifted = 0;
   texcolorlimit = 5;
   textestset = "";
   texcmpfile = "";
   texcmpfiledir = "";
   texcmpfilename = "";
   diagramfile = "";
   diagramnsteps = 5;        # number of steps at the time line
   diagramyellowbg = 0;      # Should the background be colored in SCIP-HP-yellow ?
   diagramgerman = 0;        # Soll die Beschriftung deutsch sein?
   thesisnames = 0;
   nsetnames = 0;
   onlygroup = 0;
   group = "default";

   problistlen = 0;
   nsolver = 0;
   nprobs[nsolver] = 0;
   fulltotaltime = 0.0;
}
/^=group=/ {
   group = $2;
}
/^=setname= / {
   if( setorder[$2] == 0 )
   {
      nsetnames++;
      setorder[$2] = nsetnames;
      setname[$2] = $3;
      for( i = 4; i <= NF; i++ )
         setname[$2] = setname[$2]" "$i;
   }
   setingroup[$2,group] = 1;
}
/^@02 timelimit: / {
   timelimit[nsolver] = $3;
}
/^@01 / {
   if( onlygroup == 0 || setingroup[$2,onlygroup] )
   {
      solvername[nsolver] = $2;
      nsolver++;
   }
   nprobs[nsolver] = 0;
}
// {
   statuses["ok"];
   statuses["timeout"];
   statuses["unknown"];
   statuses["abort"];
   statuses["fail"];
   statuses["readerror"];
   statuses["better"];
   statuses["solved"];
   statuses["sollimit"];
   statuses["gaplimit"];
   statuses["memlimit"];
   statuses["nodelimit"];
   
   name[nsolver,nprobs[nsolver]] = $1;
   validline = 0;
   # check if this is a useable line
   if( $10 in statuses ) # BLIS, SYMPHONY
   {
      # collect data (line with problem size and simplex iterations)
      type[nsolver,nprobs[nsolver]] = "?";
      conss[nsolver,nprobs[nsolver]] = $2;
      vars[nsolver,nprobs[nsolver]] = $3;
      dualbound[nsolver,nprobs[nsolver]] = max(min($4, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($5, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $6;
      iters[nsolver,nprobs[nsolver]] = $7;
      nodes[nsolver,nprobs[nsolver]] = max($8,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($9,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $10;
      printsoltimes = 0; # additional output is only available for SCIP-.res files
      validline = 1;
   }
   if( $11 in statuses ) # from NLP-trace-files
   {
      # collect data (line with problem type, problem size and simplex iterations)
       type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $3;
      vars[nsolver,nprobs[nsolver]] = $4;
      dualbound[nsolver,nprobs[nsolver]] = max(min($5, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($6, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $7;
      iters[nsolver,nprobs[nsolver]] = $8;
      nodes[nsolver,nprobs[nsolver]] = max($9,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($10,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $11;
      printsoltimes = 0; # additional output is only available for SCIP-.res files
      validline = 1;
   }
   if( $12 in statuses ) # GUROBI, CBC
   {
      # collect data (line with original and presolved problem size and simplex iterations)
        type[nsolver,nprobs[nsolver]] = "?";
      conss[nsolver,nprobs[nsolver]] = $4;
      vars[nsolver,nprobs[nsolver]] = $5;
      dualbound[nsolver,nprobs[nsolver]] = max(min($6, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $8;
      iters[nsolver,nprobs[nsolver]] = $9;
      nodes[nsolver,nprobs[nsolver]] = max($10,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($11,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $12;
      printsoltimes = 0; # additional output is only available for SCIP-.res files
      validline = 1;
   }
   if( $13 in statuses ) # GLPK, CPLEX, SCIP without columns displaying times to first and best solution
   {
      # collect data (line with problem type, original and presolved problem size and simplex iterations)
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($8, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $9;
      iters[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($12,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $13;
      printsoltimes = 0; # additional output is only available for SCIP-.res files
      validline = 1;
   }

   if( $15 in statuses ) # SCIP with solution times to first/last
   {
      # collect data (line with problem type, original and presolved problem size and simplex iterations)
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($8, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $9;
      iters[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($12,mintime),0.1);
      timetofirst[nsolver,nprobs[nsolver]] = fracceil(max($13,mintime),0.1);
      timetobest[nsolver, nprobs[nsolver]] = fracceil(max($14, mintime), 0.1);
      status[nsolver,nprobs[nsolver]] = $15;
      validline = 1;
   }

   if( validline )
   {
      # postprocessing of information
      if( status[nsolver,nprobs[nsolver]] == "better" )
	 status[nsolver,nprobs[nsolver]] = "timeout";
      if( status[nsolver,nprobs[nsolver]] == "sollimit" || status[nsolver,nprobs[nsolver]] == "gaplimit" || status[nsolver,nprobs[nsolver]] == "solved" )
	 status[nsolver,nprobs[nsolver]] = "ok";
   
      if( status[nsolver,nprobs[nsolver]] == "timeout" || status[nsolver,nprobs[nsolver]] == "nodelimit" ||  status[nsolver,nprobs[nsolver]] == "memlimit") 
	 hitlimit[nsolver,nprobs[nsolver]] = 1;
      else
	 hitlimit[nsolver,nprobs[nsolver]] = 0;
      probidx[$1,nsolver] = nprobs[nsolver];
      probcnt[$1]++;
      nprobs[nsolver]++;
      if( probcnt[$1] == 1 )
      {
	 problist[problistlen] = $1;
	 problistlen++;
      }
   }
}
END {
   if( onlygroup > 0 && nsolver == 1 && solvername[1] == "SCIP:default" )
   {
      printf("only SCIP:default setting found\n");
      exit 1;
   }
   if( nsolver == 0 )
   {
      printf("no instances found in log file\n");
      exit 1;
   }

   # tex comparison file: either directly as 'texcmpfile' or as pair 'texcmpfiledir/texcmpfilename'
   if( texcmpfile == "" && texcmpfiledir != "" && texcmpfilename != "" )
      texcmpfile = texcmpfiledir "/" texcmpfilename;

   # process exclude string
   n = split(exclude, a, ",");
   for( i = 1; i <= n; i++ )
      excluded[a[i]] = 1;

   # initialize means
   for( s = 0; s < nsolver; ++s )
   {
      # cat: -1 - all within time limit, 0 - all, 1 - different path, 2 - equal path, 3 - all timeout
      for( cat = -1; cat <= 3; cat++ )
      {
         nevalprobs[s,cat] = 0;
         nprocessedprobs[s,cat] = 0;
         timetotal[s,cat] = 0.0;
         nodetotal[s,cat] = 0.0;
         timegeom[s,cat] = 1.0;
         nodegeom[s,cat] = 1.0;
         timeshiftedgeom[s,cat] = timegeomshift;
         timetofirstgeom[s,cat] = 1.0;
         timetofirstshiftedgeom[s,cat] = timegeomshift;
         timetobestgeom[s,cat] = 1.0;
         timetobestshiftedgeom[s,cat] = timegeomshift;
         nodeshiftedgeom[s,cat] = nodegeomshift;
         reftimetotal[s,cat] = 0.0;
         refnodetotal[s,cat] = 0.0;
         reftimegeom[s,cat] = 1.0;
         reftimetofirstgeom[s,cat] = 1.0;
         reftimetobestgeom[s,cat] = 1.0;
         refnodegeom[s,cat] = 1.0;
         reftimeshiftedgeom[s,cat] = timegeomshift;
         refnodeshiftedgeom[s,cat] = nodegeomshift;
         reftimetofirstshiftedgeom[s,cat] = timegeomshift;
         reftimetobestshiftedgeom[s,cat] = timegeomshift;
         
         wins[s,cat] = 0;
         nsolved[s,cat] = 0;
         ntimeouts[s,cat] = 0;
         nfails[s,cat] = 0;
         better[s,cat] = 0;
         worse[s,cat] = 0;
         betterobj[s,cat] = 0;
         worseobj[s,cat] = 0;
         feasibles[s,cat] = 0;
         score[s,cat] = 1.0;
      }
   }
   besttimegeom = 1.0;
   besttimetofirstgeom = 1.0;
   besttimetobestgeom = 1.0;
   bestnodegeom = 1.0;
   besttimeshiftedgeom = timegeomshift;
   besttimetofirstshiftedgeom = timegeomshift;
   besttimetobestshiftedgeom = timegeomshift;
   bestnodeshiftedgeom = nodegeomshift;
   bestnsolved = 0;
   bestntimeouts = 0;
   bestnfails = 0;
   bestbetter = 0;
   bestbetterobj = 0;
   bestfeasibles = 0;

   # calculate the order in which the columns should be printed: CPLEX < SCIP, default < non-default
   for( s = 0; s < nsolver; ++s )
   {
      sname = solvername[s];
      for( o = 0; o < s; ++o )
      {
         i = printorder[o];
         iname = solvername[i];
         if( nsetnames > 0 )
         {
            # use order given by =setname= entries
            if( setorder[sname] < setorder[iname] )
               break;
         }
         else
         {
            # use alphabetical order, but put CPLEX before SCIP and "default" before all others
            if( substr(sname, 1, 5) == "CPLEX" && substr(iname, 1, 5) != "CPLEX" )
               break;
            if( substr(sname, 1, 5) == substr(iname, 1, 5) &&
               match(sname, "default") != 0 && match(iname, "default") == 0 )
               break;
            if( substr(sname, 1, 5) == substr(iname, 1, 5) &&
               (match(sname, "default") == 0) == (match(iname, "default") == 0) &&
               sname < iname )
               break;
         }
      }
      for( j = s-1; j >= o; --j )
         printorder[j+1] = printorder[j];
      printorder[o] = s;
   }

   # print headers
   for( o = 0; o < nsolver; ++o )
   {
       s = printorder[o];
       sname = solvername[s];
       if( o == 0 )
       {
	   if( printsoltimes )
	   {
	       if ( length(sname) <= 58 )
		   printf(" %58s |", sname);
	       else
		   printf(" *%57s |", substr(sname, length(sname)-58));
	   }
	   else
	   {
	       if ( length(sname) <= 39 )
		   printf(" %39s |", sname)
	       else
		   printf(" *%38s |", substr(sname, length(sname)-39));
	   }
       }
       else
       {
	   if( short )
	   {
	       if( length(sname) <= 19 )
		   printf("%19s |", sname);
	       else
		   printf("*%18s |", substr(sname, length(sname)-19));
	   }
	   else if( printsoltimes )
	   {
	       if( length(sname) <= 47 )
		   printf("%47s |", sname);
	       else
		   printf("*%46s |", substr(sname, length(sname)-47));
	   }
	   else
	   {            
	       if( length(sname) <= 33 )
		   printf("%31s |", sname);
	       else
		   printf("*%30s |", substr(sname, length(sname)-31));
	   }
       }
   }
   printf("\n");
   printhline(nsolver,short, printsoltimes);
   printf("  Name              |");
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 || short )
         printf("F|   Nodes |   Time |");
      else
         printf("F|   Nodes |   Time | NodQ | TimQ |");
      if( printsoltimes )
      {
        if( s == 0 )
          printf(" ToFirst | ToLast |");
        else 
          printf(" FirQ | LasQ |");
      }
   }
   printf(" bounds check\n");
   printhline(nsolver,short, printsoltimes);

   # tex comparison headers
   if( texcmpfile != "" )
   {
      printf("{\\sffamily\n") > texcmpfile;
      printf("\\scriptsize\n") > texcmpfile;
      printf("\\setlength{\\extrarowheight}{1pt}\n") > texcmpfile;
      printf("\\setlength{\\tabcolsep}{2pt}\n") > texcmpfile;
      printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n") > texcmpfile;
      printf("\\newcommand{\\spc}{\\hspace{2em}}\n") > texcmpfile;

      # add names of solvers
      for ( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         solverextension = solverextension "i";
         printf("\\newcommand{\\solvername%s}{%s}\n", solverextension, solvername[s]) > texcmpfile;
      }

      printf("\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}l") > texcmpfile;
      for( s = 0; s < nsolver; ++s )
         printf("@{\\spc}rr") > texcmpfile;
      printf("@{}}\n") > texcmpfile;
      printf("\\toprule\n") > texcmpfile;
      solverextension = "";
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         solverextension = solverextension "i";
         printf("& \\multicolumn{2}{@{\\spc}c%s}{\\solvername%s} ", o < nsolver-1 ? "@{\\spc}" : "", solverextension) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
         printf("& Nodes & Time ") > texcmpfile;
      printf("\\\\\n") > texcmpfile;
      printf("\\midrule\n") > texcmpfile;
   }

   # display the problem results and calculate mean values
   for( i = 0; i < problistlen; ++i )
   {
      p = problist[i];
      if( length(p) > 18 )
         shortp = substr(p, length(p)-17, 18);
      else
         shortp = p;

      line = sprintf("%-18s", shortp);
      fail = 0;
      readerror = 0;
      unprocessed = 0;
      mindb = +infinity;
      maxdb = -infinity;
      minpb = +infinity;
      maxpb = -infinity;
      itercomp = -1;
      nodecomp = -1;
      timecomp = -1;
      timetofirstcomp = -1;
      timetobestcomp = -1;
      besttime = +infinity;
      besttimetofirst = +infinity;
      besttimetobest = +infinity;
      bestnodes = +infinity;
      worsttime = -infinity;
      worstnodes = -infinity;
      worstiters = -infinity;
      worsttimetofirst = -infinity;
      worsttimetobest = -infinity;
      nthisunprocessed = 0;
      nthissolved = 0;
      nthistimeouts = 0;
      nthisfails = 0;
      ismini = 0;
      ismaxi = 0;
      mark = " ";
      countprob = 1;
      notimeout = 1;

      # check for exclusion
      if( excluded[p] )
      {
         unprocessed = 1;
         countprob = 0;
      }

      # find best and worst run and check whether this instance should be counted in overall statistics
      for( s = 0; s < nsolver; ++s )
      {
         pidx = probidx[p,s];
         processed = (pidx != "");

         # make sure, nodes and time are non-zero for geometric means
         nodes[s,pidx] = max(nodes[s,pidx], 1);
         time[s,pidx] = max(time[s,pidx], mintime);
         fulltotaltime += time[s,pidx];

         # If we got a timeout although the time limit has not been reached (e.g., due to a memory limit),
         # we assume that the run would have been continued with the same nodes/sec.
         # Set the time to the time limit and increase the nodes accordingly.
         # if( status[s,pidx] == "timeout" && time[s,pidx] < timelimit[s] )
         # {
         #    nodes[s,pidx] *= timelimit[s]/time[s,pidx];
         #    time[s,pidx] = timelimit[s];
         # }

         # if the solver exceeded the timelimit, set status accordingly
         if( (status[s,pidx] == "ok" || status[s,pidx] == "unknown") && timelimit[s] > 0.0 && time[s,pidx] > timelimit[s] )
         {
            status[s,pidx] = "timeout";
            time[s,pidx] = timelimit[s];
         }

         # check if all solvers processed the problem
         if( !processed )
         {
            marker = "?";
            unprocessed = 1;
         }

         # check if solver ran successfully (i.e., no abort nor fail)
         if( processed && (status[s,pidx] == "ok" || status[s,pidx] == "unknown" || status[s,pidx] == "timeout" || status[s,pidx] == "nodelimit" || status[s,pidx] == "memlimit") )
         {
            besttime = min(besttime, time[s,pidx]);
            bestnodes = min(bestnodes, nodes[s,pidx]);
            besttimetofirst = min(besttimetofirst, timetofirst[s,pidx]);
            besttimetobest = min(besttimetobest, timetobest[s,pidx]);
            worsttime = max(worsttime, time[s,pidx]);
            worstnodes = max(worstnodes, nodes[s,pidx]);
            worstiters = max(worstiters, iters[s,pidx]);
            worsttimetofirst = max(worsttimetofirst, timetofirst[s, pidx]);
            worsttimetobest = max(worsttimetobest, timetobest[s, pidx]);
         }
         else
            countprob = 0;
      }
      worsttime = max(worsttime, mintime);
      worsttimetofirst = max(worsttimetofirst, mintime);
      worsttimetobest = max(worsttimetobest, mintime);
      worstnodes = max(worstnodes, 1);
      worstiters = max(worstiters, 0);

      # check for each solver if it has same path as reference solver -> category
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");

         if( !processed )
         {
            notimeout = 0;
            continue;
         }
         else
         {
            if ( status[s,pidx] == "timeout" )
            {
               # If memory limit was exceeded or we hit a hard time/memory limit,
               # replace time and nodes by worst time and worst nodes of all runs.
               # Note this also takes action if the time limits of the runs are
               # different: in this case we set the values to the worst case.
               # if ( time[s,pidx] < 0.99*worsttime || nodes[s,pidx] <= 1 )
               # {
               #    iters[s,pidx] = worstiters+s; # make sure this is not treated as equal path
               #    nodes[s,pidx] = worstnodes;
               #    time[s,pidx] = worsttime;
               # }
            }
         }

         if( nodecomp == -1 )
         {
            itercomp = iters[s,pidx];
            nodecomp = nodes[s,pidx];
            timecomp = time[s,pidx];
            timeoutcomp = hitlimit[s,pidx];
            timetofirstcomp = max(mintime, timetofirst[s,pidx]);
            timetobestcomp = max(mintime, timetobest[s,pidx]);
         }
         iseqpath = (iters[s,pidx] == itercomp && nodes[s,pidx] == nodecomp);
         hastimeout = timeoutcomp || hitlimit[s,pidx];
         notimeout = notimeout && !timeoutcomp &&  !hitlimit[s,pidx];

         # which category?
         if( hastimeout )
            category[s] = 3;
         else if( iseqpath )
            category[s] = 2;
         else
            category[s] = 1;
      }

      # evaluate instance for all solvers
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");

         if( processed && name[s,pidx] != p )
            printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

         # check if solver ran successfully (i.e., no abort nor fail)
         if( processed )
         {
            if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
            {
               marker = " ";
               if( !unprocessed )
               {
                  if ( notimeout )
                     nsolved[s,-1]++;
                  nsolved[s,0]++;
                  nsolved[s,category[s]]++;
                  nthissolved++;
               }
            }
            else if( hitlimit[s,pidx] )
            {
               marker = ">";
               notimeout = 0;
               if( !unprocessed )
               {
                  if( countprob )
                  {
                     ntimeouts[s,0]++;
                     ntimeouts[s,category[s]]++;
                     nthistimeouts++;
                  }
               }
            }
            else
            {
               marker = "!";
               fail = 1;
               notimeout = 0;
               if( status[s,pidx] == "readerror" )
                  readerror = 1;
               if( !unprocessed )
               {
                  nfails[s,0]++;
                  nfails[s,category[s]]++;
                  nthisfails++;
               }
            }
         }

         if( primalbound[s,pidx] < infinity )
            feasmark = " ";
         else
            feasmark = "#";

         if( processed && !fail )
         {
            mindb = min(mindb, dualbound[s,pidx]);
            maxdb = max(maxdb, dualbound[s,pidx]);
            minpb = min(minpb, primalbound[s,pidx]);
            maxpb = max(maxpb, primalbound[s,pidx]);
            ismini = ismini || (primalbound[s,pidx] > dualbound[s,pidx] + 1e-06);
            ismaxi = ismaxi || (primalbound[s,pidx] < dualbound[s,pidx] - 1e-06);
         }

         # print statistics
         if( !processed )
            line = sprintf("%s           -        -", line);
         else
         {
            if( printgap && hitlimit[s,pidx] && gap[s,pidx] != "--" && gap[s,pidx] != "Large" )
              line = sprintf("%s %s%10d %7.2f%%", line, feasmark, nodes[s,pidx], gap[s,pidx]);
            else
              line = sprintf("%s %s%10d %s%7.1f", line, feasmark, nodes[s,pidx], marker, time[s,pidx]);
            if( printsoltimes && o == 0 )
               line = sprintf("%s  %8.1f %8.1f", line, timetofirst[s,pidx], timetobest[s, pidx] );
         }
         if( o > 0 && !short )
         {
            if( !processed )
               line = sprintf("%s      -", line);
            else if( nodes[s,pidx]/nodecomp > 999.99 )
               line = sprintf("%s  Large", line);
            else
               line = sprintf("%s %6.2f", line, nodes[s,pidx]/nodecomp);
            if( !processed )
               line = sprintf("%s      -", line);
            else if( time[s,pidx]/timecomp > 999.99 )
               line = sprintf("%s  Large", line);
            else
               line = sprintf("%s %6.2f", line, time[s,pidx]/timecomp);
            if( processed &&
		(timeoutcomp != hitlimit[s,pidx] ||
                 nodes[s,pidx] > markworsenodes * nodecomp ||
                 nodes[s,pidx] < 1.0/markbetternodes * nodecomp ||
                 isfaster(time[s,pidx], timecomp, markbettertime) ||
                 isslower(time[s,pidx], timecomp, markworsetime)) )
               mark = "*";
         }
         if( o > 0 && printsoltimes )
         {
            if( !processed )
               line = sprintf("%s      -", line);
            else if( timetofirst[s,pidx]/timetofirstcomp > 999.99 )
               line = sprintf("%s  Large", line);
            else
               line = sprintf("%s %6.2f", line, timetofirst[s,pidx]/timetofirstcomp);
            if( !processed )
               line = sprintf("%s      -", line);
            else if( timetobest[s,pidx]/timetobestcomp > 999.99 )
               line = sprintf("%s  Large", line);
            else
               line = sprintf("%s %6.2f", line, timetobest[s,pidx]/timetobestcomp);
         }
      }

      # update the best status information
      if( nthissolved > 0 )
        bestnsolved++;
      else if( nthistimeouts > 0 )
        bestntimeouts++;
      else if( nthisfails == nsolver - nthisunprocessed )
        bestnfails++;

      # check for inconsistency in the primal and dual bounds
      if( readerror )
      {
         line = sprintf("%s  readerror", line);
         mark = " ";
      }
      else if( fail )
      {
         line = sprintf("%s  fail", line);
         mark = " ";
      }
      else if( consistency &&
         ((ismini && ismaxi) ||
            (ismini && maxdb - minpb > 1e-5 * max(max(abs(maxdb), abs(minpb)), 1.0)) ||
            (ismaxi && maxpb - mindb > 1e-5 * max(max(abs(maxpb), abs(mindb)), 1.0)) ||
            (!ismini && !ismaxi && abs(maxpb - minpb) > 1e-5 * max(abs(maxpb), 1.0))) )
      {
         line = sprintf("%s  inconsistent", line);
         fail = 1;
         mark = " ";
      }
      else if( excluded[p] )
      {
         line = sprintf("%s  excluded", line);
         mark = " ";
      }
      else if( unprocessed )
      {
         line = sprintf("%s  unprocessed", line);
         mark = " ";
      }
      else
         line = sprintf("%s  ok", line);

      # calculate number of instances for which feasible solution has been found
      hasfeasible = 0;
      if( !unprocessed )
      {
         for( s = 0; s < nsolver; ++s )
         {
            if ( notimeout )
               nprocessedprobs[s,-1]++;
            nprocessedprobs[s,0]++;
            nprocessedprobs[s,category[s]]++;
            pidx = probidx[p,s];
            if( primalbound[s,pidx] < infinity ) 
            {
               if ( notimeout )
                  feasibles[s,-1]++;
               feasibles[s,0]++;
               feasibles[s,category[s]]++;
               hasfeasible = 1;
            }
         }
         if( hasfeasible )
           bestfeasibles++;
      }

      if( (!onlymarked || mark == "*") && (!onlyprocessed || !unprocessed) &&
          (!onlyfeasible || hasfeasible) && (!onlyinfeasible || !hasfeasible) &&
          (!onlyfail || fail) )
      {
         printf("%s %s\n", mark, line);

         # tex comparison file
         if( texcmpfile != "" )
         {
            printf("%s ", texstring(p)) > texcmpfile;
            ref = printorder[0];
            refnodes = nodes[ref,probidx[p,ref]];
            reftime = time[ref,probidx[p,ref]];
            refstatus = status[ref,probidx[p,ref]];
            refhitlimit = hitlimit[ref,probidx[p,ref]];
            for( o = 0; o < nsolver; ++o )
            {
               s = printorder[o];
               pidx = probidx[p,s];

               if( hitlimit[s,pidx] )
                  timeoutmarker = "\\g";
               else
                  timeoutmarker = "  ";

               if( nodes[s,pidx] <= 0.5*refnodes && !hitlimit[s,pidx] )
                  nodecolor = "red";
               else if( nodes[s,pidx] >= 2.0*refnodes && !refhitlimit )
                  nodecolor = "blue";
               else
                  nodecolor = "black";

               if( (time[s,pidx] <= 0.5*reftime && !hitlimit[s,pidx]) ||
                  (!hitlimit[s,pidx] && refhitlimit) )
                  timecolor = "red";
               else if( (time[s,pidx] >= 2.0*reftime && !refhitlimit) ||
                  (hitlimit[s,pidx] && !refhitlimit) )
                  timecolor = "blue";
               else
                  timecolor = "black";

               if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" || status[s,pidx] == "timeout" || status[s,pidx] == "memlimit" || status[s,pidx] == "nodelimit" )
                  printf("&\\textcolor{%s}{%s %8s} &\\textcolor{%s}{%s %8.1f} ",
                     nodecolor, timeoutmarker, texint(nodes[s,pidx]), timecolor, timeoutmarker, time[s,pidx]) > texcmpfile;
               else
                  printf("&        --- &        --- ") > texcmpfile;
            }
            printf("\\\\\n") > texcmpfile;


         }
      }

      # calculate totals and means for instances where no solver failed
      if( !fail && !unprocessed &&
          (!onlyfeasible || hasfeasible) && (!onlyinfeasible || !hasfeasible) )
      {
         reftime = time[printorder[0],probidx[p,printorder[0]]];
         refhitlimit = hitlimit[printorder[0],probidx[p,printorder[0]]];
         refnodes = nodes[printorder[0],probidx[p,printorder[0]]];
         refobj = primalbound[printorder[0],probidx[p,printorder[0]]];
         reftimetofirst = timetofirst[printorder[0],probidx[p,printorder[0]]];
         reftimetobest = timetobest[printorder[0],probidx[p,printorder[0]]];
         hasbetter = 0;
         hasbetterobj = 0;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            for( cat = 0; cat <= 3; cat = 3*cat + category[s] )
            {
               nevalprobs[s,cat]++;
               nep = nevalprobs[s,cat];
               timetotal[s,cat] += time[s,pidx];
               timetofirsttotal[s,cat] += timetofirst[s,pidx];
               timetobesttotal[s, cat] += timetobest[s, pidx];
               nodetotal[s,cat] += nodes[s,pidx];
               timegeom[s,cat] = timegeom[s,cat]^((nep-1)/nep) * time[s,pidx]^(1.0/nep);
               timetofirstgeom[s,cat] = timetofirstgeom[s,cat]^((nep-1)/nep) * max(timetofirst[s,pidx], mintime)^(1.0/nep);
               timetobestgeom[s,cat] = timetobestgeom[s,cat]^((nep-1)/nep) * max(timetobest[s,pidx])^(1.0/nep);
               nodegeom[s,cat] = nodegeom[s,cat]^((nep-1)/nep) * nodes[s,pidx]^(1.0/nep);
               timeshiftedgeom[s,cat] = timeshiftedgeom[s,cat]^((nep-1)/nep) * (time[s,pidx]+timegeomshift)^(1.0/nep);
               timetofirstshiftedgeom[s,cat] = timetofirstshiftedgeom[s,cat]^((nep-1)/nep) * max(timetofirst[s,pidx]+timegeomshift, 1.0)^(1.0/nep);
               timetobestshiftedgeom[s,cat] = timetobestshiftedgeom[s,cat]^((nep-1)/nep) * max(timetobest[s,pidx]+timegeomshift, 1.0)^(1.0/nep);
               nodeshiftedgeom[s,cat] = nodeshiftedgeom[s,cat]^((nep-1)/nep) * (nodes[s,pidx]+nodegeomshift)^(1.0/nep);
               reftimetotal[s,cat] += reftime;
               reftimetofirsttotal[s,cat] += reftimetofirst;
               reftimetobesttotal[s,cat] += reftimetobest; 
               refnodetotal[s,cat] += refnodes;
               reftimegeom[s,cat] = reftimegeom[s,cat]^((nep-1)/nep) * reftime^(1.0/nep);
               reftimetofirstgeom[s,cat] = reftimetofirstgeom[s,cat]^((nep-1)/nep) * reftimetofirst^(1.0/nep);
               reftimetobestgeom[s,cat] = reftimetobestgeom[s,cat]^((nep-1)/nep) * reftimetobest^(1.0/nep);
               refnodegeom[s,cat] = refnodegeom[s,cat]^((nep-1)/nep) * refnodes^(1.0/nep);
               reftimeshiftedgeom[s,cat] = reftimeshiftedgeom[s,cat]^((nep-1)/nep) * (reftime+timegeomshift)^(1.0/nep);
               reftimetofirstshiftedgeom[s,cat] = reftimetofirstshiftedgeom[s,cat]^((nep-1)/nep) * (reftimetofirst+timegeomshift)^(1.0/nep);
               reftimetobestshiftedgeom[s,cat] = reftimetobestshiftedgeom[s,cat]^((nep-1)/nep) * (reftimetobest+timegeomshift)^(1.0/nep);
               refnodeshiftedgeom[s,cat] = refnodeshiftedgeom[s,cat]^((nep-1)/nep) * (refnodes+nodegeomshift)^(1.0/nep);
               if( time[s,pidx] <= wintolerance*besttime )
                  wins[s,cat]++;
               if( !hitlimit[s,pidx] && (isfaster(time[s,pidx], reftime, wintolerance) || refhitlimit))
               {
                  better[s,cat]++;
                  hasbetter = 1;
               }
               else if( !refhitlimit && (isslower(time[s,pidx], reftime, wintolerance) || (hitlimit[s,pidx])))
                  worse[s,cat]++;
               pb = primalbound[s,pidx];
               if( (ismini && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
                   (ismaxi && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ) {
                  betterobj[s,cat]++;
                  hasbetterobj = 1;
               }
               else if( (ismini && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
                        (ismaxi && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
                  worseobj[s,cat]++;
               thisscore = reftime/time[s,pidx];
               thisscore = max(thisscore, 1/maxscore);
               thisscore = min(thisscore, maxscore);
               score[s,cat] = score[s,cat]^((nep-1)/nep) * thisscore^(1.0/nep);
            }
         }
         s = printorder[0];
         besttimegeom = besttimegeom^((nevalprobs[s,0]-1)/nevalprobs[s,0]) * besttime^(1.0/nevalprobs[s,0]);
         bestnodegeom = bestnodegeom^((nevalprobs[s,0]-1)/nevalprobs[s,0]) * bestnodes^(1.0/nevalprobs[s,0]);
         besttimeshiftedgeom = besttimeshiftedgeom^((nevalprobs[s,0]-1)/nevalprobs[s,0]) * (besttime+timegeomshift)^(1.0/nevalprobs[s,0]);
         bestnodeshiftedgeom = bestnodeshiftedgeom^((nevalprobs[s,0]-1)/nevalprobs[s,0]) * (bestnodes+nodegeomshift)^(1.0/nevalprobs[s,0]);
         if( hasbetter )
            bestbetter++;
         if( hasbetterobj )
           bestbetterobj++;

         # once again for the case in which all instances have been solved to optimality
         if ( notimeout )
         {
            for( s = 0; s < nsolver; ++s )
            {
               pidx = probidx[p,s];
               cat = -1;
               nevalprobs[s,cat]++;
               nep = nevalprobs[s,cat];
               timetotal[s,cat] += time[s,pidx];
               nodetotal[s,cat] += nodes[s,pidx];
               timegeom[s,cat] = timegeom[s,cat]^((nep-1)/nep) * time[s,pidx]^(1.0/nep);
               nodegeom[s,cat] = nodegeom[s,cat]^((nep-1)/nep) * nodes[s,pidx]^(1.0/nep);
               timeshiftedgeom[s,cat] = timeshiftedgeom[s,cat]^((nep-1)/nep) * (time[s,pidx]+timegeomshift)^(1.0/nep);
               nodeshiftedgeom[s,cat] = nodeshiftedgeom[s,cat]^((nep-1)/nep) * (nodes[s,pidx]+nodegeomshift)^(1.0/nep);
               reftimetotal[s,cat] += reftime;
               refnodetotal[s,cat] += refnodes;
               reftimegeom[s,cat] = reftimegeom[s,cat]^((nep-1)/nep) * reftime^(1.0/nep);
               refnodegeom[s,cat] = refnodegeom[s,cat]^((nep-1)/nep) * refnodes^(1.0/nep);
               reftimeshiftedgeom[s,cat] = reftimeshiftedgeom[s,cat]^((nep-1)/nep) * (reftime+timegeomshift)^(1.0/nep);
               refnodeshiftedgeom[s,cat] = refnodeshiftedgeom[s,cat]^((nep-1)/nep) * (refnodes+nodegeomshift)^(1.0/nep);
               if( time[s,pidx] <= wintolerance*besttime )
                  wins[s,cat]++;
               if( isfaster(time[s,pidx], reftime, wintolerance) )
                  better[s,cat]++;
               else if( isslower(time[s,pidx], reftime, wintolerance) )
                  worse[s,cat]++;
               pb = primalbound[s,pidx];
               if( (ismini && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
                   (ismaxi && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ) {
                  betterobj[s,cat]++;
               }
               else if( (ismini && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
                        (ismaxi && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
                  worseobj[s,cat]++;
               thisscore = reftime/time[s,pidx];
               thisscore = max(thisscore, 1/maxscore);
               thisscore = min(thisscore, maxscore);
               score[s,cat] = score[s,cat]^((nep-1)/nep) * thisscore^(1.0/nep);
            }
         }
      }
   }
   printhline(nsolver,short, printsoltimes);

   # make sure total time and nodes is not zero
   for( s = 0; s < nsolver; ++s )
   {
      for( cat = -1; cat <= 3; cat++ )
      {
         nodetotal[s,cat] = max(nodetotal[s,cat], 1);
         refnodetotal[s,cat] = max(refnodetotal[s,cat], 1);
         timetotal[s,cat] = max(timetotal[s,cat], mintime);
         reftimetotal[s,cat] = max(reftimetotal[s,cat], mintime);
      }
   }

   # print solvers' overall statistics
   probnumstr = "("nevalprobs[printorder[0],0]")";
   printf("%-14s %5s", "total", probnumstr);
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 || short )
      {
          printf(" %11d %8d", nodetotal[s,0], timetotal[s,0]);
          if( o == 0 && printsoltimes )
             printf(" %9d %8d" , timetofirsttotal[s,0], timetobesttotal[s,0]);
      }
      else
      {
          printf(" %11d %8d              ", nodetotal[s,0], timetotal[s,0]);
          if( printsoltimes )
          {
             referencesolvername = printorder[0];
             comptimetofirst = timetofirsttotal[s,0]/(max(timetofirsttotal[referencesolvername,0], 1));
             comptimetobest = timetobesttotal[s,0]/(max(timetobesttotal[referencesolvername,0], 1));
             printf("%7.2f %6.2f", comptimetofirst, comptimetobest);
          }
      }
   }
   printf("\n");
   printf("%-20s", "geom. mean");
   nodegeomcomp = -1;
   timegeomcomp = -1;
   nodetotalcomp = -1;
   timetotalcomp = -1;
   timetofirstcomp = -1;
   timetobestcomp = -1;
   timetofirstgeomcomp = -1;
   timetobestgeomcomp = -1;

   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 || short )
      {
          printf(" %11d %8.1f", nodegeom[s,0], timegeom[s,0]);
          if( o == 0 && printsoltimes )
             printf(" %9.1f %8.1f", timetofirstgeom[s,0], timetobestgeom[s,0]);

         if( nodegeomcomp < 0 )
            nodegeomcomp = nodegeom[s,0];
         if( timegeomcomp < 0 )
            timegeomcomp = timegeom[s,0];
         if( nodetotalcomp < 0 )
            nodetotalcomp = nodetotal[s,0];
         if( timetotalcomp < 0 )
            timetotalcomp = timetotal[s,0];
         if( timetofirstcomp < 0 )
             timetofirstcomp = timetofirsttotal[s,0];
         if( timetobestcomp < 0 )
             timetobestcomp = timetobesttotal[s,0];
         if( timetofirstgeomcomp < 0 )
             timetofirstgeomcomp = timetofirstgeom[s,0];
         if( timetobestgeomcomp < 0 )
             timetobestgeomcomp = timetobestgeom[s,0];
      }
      else
      {
          printf(" %11d %8.1f %6.2f %6.2f", nodegeom[s,0], timegeom[s,0], nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp);

          if( printsoltimes )
             printf(" %6.2f %6.2f", timetofirstgeom[s,0]/timetofirstgeomcomp, timetobestgeom[s,0]/timetobestgeomcomp);
      }
   }
   printf("\n");
   printf("%-20s", "shifted geom.");
   nodeshiftedgeomcomp = -1;
   timeshiftedgeomcomp = -1;
   timetofirstshiftedgeomcomp = -1;
   timetobestshiftedgeomcomp = -1;
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      for( cat = -1; cat <= 3; cat++ )
      {
         nodeshiftedgeom[s,cat] -= nodegeomshift;
         timeshiftedgeom[s,cat] -= timegeomshift;
         timetofirstshiftedgeom[s,cat] -= timegeomshift;
         timetobestshiftedgeom[s,cat] -= timegeomshift;
         nodeshiftedgeom[s,cat] = max(nodeshiftedgeom[s,cat], 1);
         timeshiftedgeom[s,cat] = max(timeshiftedgeom[s,cat], mintime);
         timetofirstshiftedgeom[s,cat] = max(timetofirstshiftedgeom[s,cat], mintime);
         timetobestshiftedgeom[s,cat] = max(timetobestshiftedgeom[s,cat], mintime);
         refnodeshiftedgeom[s,cat] -= nodegeomshift;
         reftimeshiftedgeom[s,cat] -= timegeomshift;
         refnodeshiftedgeom[s,cat] = max(refnodeshiftedgeom[s,cat], mintime);
         reftimeshiftedgeom[s,cat] = max(reftimeshiftedgeom[s,cat], mintime);
      }
      if( o == 0 || short )
      {
          printf(" %11d %8.1f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0]);

          if( o == 0 && printsoltimes )
             printf(" %9.1f %8.1f", timetofirstshiftedgeom[s,0], timetobestshiftedgeom[s,0]);

         if( nodeshiftedgeomcomp < 0 )
            nodeshiftedgeomcomp = nodeshiftedgeom[s,0];
         if( timeshiftedgeomcomp < 0 )
            timeshiftedgeomcomp = timeshiftedgeom[s,0];
         if( timetofirstshiftedgeomcomp < 0 )
             timetofirstshiftedgeomcomp = timetofirstshiftedgeom[s,0];
         if( timetobestshiftedgeomcomp < 0 )
             timetobestshiftedgeomcomp = timetobestshiftedgeom[s,0];
      }
      else
      {
         printf(" %11d %8.1f %6.2f %6.2f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0],
                nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp);

         if( printsoltimes )
            printf(" %6.2f %6.2f", timetofirstshiftedgeom[s,0]/timetofirstshiftedgeomcomp, timetobestshiftedgeom[s,0]/timetobestshiftedgeomcomp);
      }
   }
   bestnodeshiftedgeom -= nodegeomshift;
   besttimeshiftedgeom -= timegeomshift;
   bestnodeshiftedgeom = max(bestnodeshiftedgeom, 1.0);
   besttimeshiftedgeom = max(besttimeshiftedgeom, 1.0);
   
   printf("\n");

   #since the rows of the quotients are not printed, print the quotients of the geometric means
   if( short )
   {
      printf("quot. geom. mean                         ");
      for( o = 0; o < nsolver; ++o )
      {
         if( o > 0 )
         {
            s = printorder[o];
            printf("      %6.2f   %6.2f",nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp);
         }
      }
      printf("\n");
      printf("quot. sh. geom. mean                     ");
      for( o = 0; o < nsolver; ++o )
      {
         if( o > 0 )
         {
            s = printorder[o];
            printf("      %6.2f   %6.2f",nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp);
         }
      }
      printf("\n");
      printf("percent not solved                 ");
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("%6.2f",100-100*nsolved[s,0]/nprocessedprobs[s,0]);
         printf("               ");
      }
      printf("\n");
   }
   printhline(nsolver,short, printsoltimes);
   
   # tex comparison footer
   if( texcmpfile != "" )
   {
      # all problems
      printf("\\midrule\n") > texcmpfile;
      printf("geom. mean     ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodegeom[s,0]), timegeom[s,0]) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      printf("sh. geom. mean ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodeshiftedgeom[s,0]), timeshiftedgeom[s,0]) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      printf("arithm. mean   ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodetotal[s,0]/max(1,nevalprobs[s,0])), timetotal[s,0]/max(1,nevalprobs[s,0])) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      
      # add statistics for problems solved to optimality
      printf("\\midrule\n") > texcmpfile;
      printf("\\multicolumn{%d}{@{}l}{all optimal}\\\\\n", 1 + 2 * nsolver) > texcmpfile;
      printf("geom. mean     ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodegeom[s,-1]), timegeom[s,-1]) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      printf("sh. geom. mean ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodeshiftedgeom[s,-1]), timeshiftedgeom[s,-1]) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      printf("arithm. mean   ") > texcmpfile;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("& %8s & %8.1f", texint(nodetotal[s,-1]/max(1,nevalprobs[s,-1])), timetotal[s,-1]/max(1,nevalprobs[s,-1])) > texcmpfile;
      }
      printf("\\\\\n") > texcmpfile;
      printf("\\bottomrule\n") > texcmpfile;
      printf("\\end{tabular*}\n") > texcmpfile;
      printf("}\n") > texcmpfile;
   }

   if( !short )
   {
      for( cat = 0; cat <= 3; cat++ )
      {   
#         if( nprocessedprobs[cat] == 0 )
#            continue;
      
         header = (cat == -1 ? "optimal" : (cat == 0 ? "all" : (cat == 1 ? "diff" : (cat == 2 ? "equal" : "timeout"))));
         printf("\n");
         printf("%-7s                                            proc eval fail time solv wins bett wors bobj wobj feas    gnodes   shnodes   gnodesQ  shnodesQ   gtime  shtime  gtimeQ shtimeQ   score\n",
            header);
   
         for( o = 0; o < nsolver; ++o )
         {
            s = printorder[o];
	    sname = solvername[s];
            if( o == 0 )
            {
               nodegeomcomp = nodegeom[s,cat];
               timegeomcomp = timegeom[s,cat];
               nodeshiftedgeomcomp = nodeshiftedgeom[s,cat];
               timeshiftedgeomcomp = timeshiftedgeom[s,cat];
            }
            if( (o > 0 || cat == 0 || cat == -1) && nevalprobs[s,cat] > 0 )
            {
	       if ( length(sname) <= 50 )
		  printf("%-50s %4d %4d %4d %4d %4d %4d", sname, nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat],
			 ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
	       else
		  printf("*%-49s %4d %4d %4d %4d %4d %4d", substr(sname, length(sname)-48), nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat],
			 ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
 
               printf(" %4d %4d", better[s,cat], worse[s,cat]);
               printf(" %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", 
                  betterobj[s,cat], worseobj[s,cat], feasibles[s,cat],
                  nodegeom[s,cat], nodeshiftedgeom[s,cat], nodegeom[s,cat]/refnodegeom[s,cat],
                  nodeshiftedgeom[s,cat]/refnodeshiftedgeom[s,cat],
                  timegeom[s,cat], timeshiftedgeom[s,cat], timegeom[s,cat]/reftimegeom[s,cat],
                  timeshiftedgeom[s,cat]/reftimeshiftedgeom[s,cat], score[s,cat]);
            }
         }
         if( cat == 0 )
         {
            printf("%-50s           %4d %4d %4d %4s", "optimal auto settings", bestnfails, bestntimeouts, bestnsolved, "");
            printf(" %4d %4s", bestbetter, "");
            printf(" %4d %4s %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7s\n",
                   bestbetterobj, "", bestfeasibles,
                   bestnodegeom, bestnodeshiftedgeom, bestnodegeom/nodegeomcomp, bestnodeshiftedgeom/nodeshiftedgeomcomp,
                   besttimegeom, besttimeshiftedgeom, besttimegeom/timegeomcomp, besttimeshiftedgeom/timeshiftedgeomcomp,
                   "");
         }
      }

      # output the all optimal case
      header = "all optimal";
      printf("\n");
      printf("%-11s                                        proc eval fail time solv wins bett wors bobj wobj feas    gnodes   shnodes   gnodesQ  shnodesQ   gtime  shtime  gtimeQ shtimeQ   score\n", header);
      cat = -1;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
	 sname = solvername[s];
         if( o == 0 )
         {
            nodegeomcomp = nodegeom[s,cat];
            timegeomcomp = timegeom[s,cat];
            nodeshiftedgeomcomp = nodeshiftedgeom[s,cat];
            timeshiftedgeomcomp = timeshiftedgeom[s,cat];
         }
         if( (o > 0 || cat == 0 || cat == -1) && nevalprobs[s,cat] > 0 )
         {
	    if ( length(sname) <= 50 )
	       printf("%-50s %4d %4d %4d %4d %4d %4d", sname, nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat],
		      ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
	    else
	       printf("*%-49s %4d %4d %4d %4d %4d %4d", substr(sname, length(sname)-48), nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat],
		      ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
            printf(" %4d %4d", better[s,cat], worse[s,cat]);
            printf(" %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", 
                   betterobj[s,cat], worseobj[s,cat], feasibles[s,cat],
                   nodegeom[s,cat], nodeshiftedgeom[s,cat], nodegeom[s,cat]/refnodegeom[s,cat],
                   nodeshiftedgeom[s,cat]/refnodeshiftedgeom[s,cat],
                   timegeom[s,cat], timeshiftedgeom[s,cat], timegeom[s,cat]/reftimegeom[s,cat],
                   timeshiftedgeom[s,cat]/reftimeshiftedgeom[s,cat], score[s,cat]);
         }
      }
      if( cat == 0 )
      {
         printf("%-50s           %4d %4d %4d %4s", "optimal auto settings", bestnfails, bestntimeouts, bestnsolved, "");
         printf(" %4d %4s", bestbetter, "");
         printf(" %4d %4s %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7s\n",
                bestbetterobj, "", bestfeasibles,
                bestnodegeom, bestnodeshiftedgeom, bestnodegeom/nodegeomcomp, bestnodeshiftedgeom/nodeshiftedgeomcomp,
                besttimegeom, besttimeshiftedgeom, besttimegeom/timegeomcomp, besttimeshiftedgeom/timeshiftedgeomcomp,
                "");
      }
   }

   printf("\n");
   printf("total time over all settings: %.1f sec = %.1f hours = %.1f days = %.1f weeks = %.1f months\n",
      fulltotaltime, fulltotaltime/3600.0, fulltotaltime/(3600.0*24), fulltotaltime/(3600.0*24*7),
      fulltotaltime/(3600.0*24*30));

   # generate tex file
   if( texfile != "" )
   {
      hasequalpath = 0;
      for( o = 0; o < nsolver; ++o )
      {
         if( nevalprobs[s,2] > 0 )
         {
            hasequalpath = 1;
            break;
         }
      }

      printf("generating tex file <%s>\n", texfile);
      printf("{\\sffamily\n") > texfile;
      printf("\\scriptsize\n") > texfile;
      printf("\\setlength{\\extrarowheight}{1pt}\n") > texfile;
      printf("\\setlength{\\tabcolsep}{2pt}\n") > texfile;
      printf("\\newcommand{\\spc}{\\hspace{%dem}}\n", 2-hasequalpath) > texfile;

      printf("\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lrrr@{\\spc}rrrrrr@{\\spc}rrrr") > texfile;
      if( hasequalpath )
         printf("@{\\spc}rrrr") > texfile;
      printf("@{}}\n") > texfile;

      printf("\\toprule\n") > texfile;

      printf("& & & & \\multicolumn{6}{c@{\\spc}}{all instances (%d)} & \\multicolumn{4}{c@{\\spc}}{different path}",
         nevalprobs[printorder[0],0]) > texfile;
      if( hasequalpath )
         printf("& \\multicolumn{4}{c}{equal path}") > texfile;
      printf("\\\\\n") > texfile;

      printf("setting & T & fst & slw & $\\textit{n}_\\textit{gm}$ & $\\textit{n}_\\textit{sgm}$ & $\\textit{n}_\\textit{tot}$ & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$ & \\# & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$") > texfile;
      if( hasequalpath )
         printf("& \\# & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$") > texfile;
      printf("\\\\\n") > texfile;

      printf("\\midrule\n") > texfile;

      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         printf("%-45s & %4d & %3d & %3d", texsolvername(s), ntimeouts[s,0],  better[s,0], worse[s,0]) > texfile;
         printf(" & %5s & %5s & %5s & %5s & %5s & %5s",
            texcompstr(nodegeom[s,0], refnodegeom[s,0]),
            texcompstr(nodeshiftedgeom[s,0], refnodeshiftedgeom[s,0]),
            texcompstr(nodetotal[s,0], refnodetotal[s,0]),
            texcompstr(timegeom[s,0], reftimegeom[s,0]),
            texcompstr(timeshiftedgeom[s,0], reftimeshiftedgeom[s,0]),
            texcompstr(timetotal[s,0], reftimetotal[s,0])) > texfile;
         if( nevalprobs[s,1] > 0 )
         {
            printf(" & %2d & %5s & %5s & %5s",
               nevalprobs[s,1],
               texcompstr(timegeom[s,1], reftimegeom[s,1]),
               texcompstr(timeshiftedgeom[s,1], reftimeshiftedgeom[s,1]),
               texcompstr(timetotal[s,1], reftimetotal[s,1])) > texfile;
         }
         else
            printf(" &  0 &     --- &     --- &     ---") > texfile;
         if( hasequalpath )
         {
            if( nevalprobs[s,2] > 0 )
            {
               printf(" & %2d & %5s & %5s & %5s",
                  nevalprobs[s,2],
                  texcompstr(timegeom[s,2], reftimegeom[s,2]),
                  texcompstr(timeshiftedgeom[s,2], reftimeshiftedgeom[s,2]),
                  texcompstr(timetotal[s,2], reftimetotal[s,2])) > texfile;
            }
            else
               printf(" &  0 &     --- &     --- &     ---") > texfile;
         }
         printf(" \\\\\n") > texfile;
      }

      printf("\\bottomrule\n") > texfile;
      printf("\\end{tabular*}\n") > texfile;
      printf("}\n") > texfile;

      # extend tex include file
      if( texincfile != "" )
      {
         n = split(texfile, a, "/");
         texpath = "";
         for( i = 1; i < n; i++ )
            texpath = texpath a[i] "/";
         texbase = a[n];
         sub(/\.tex/, "", texbase);
         n = split(texbase, a, "_");
         texsetname = a[2];
         texgroupname = a[3];
         textestname = a[4];

         printf("\n") >> texincfile;
         printf("\\begin{table}[hp]\n") > texincfile;
         printf("\\input{Tables/mip/auto/%s}\n", texbase) > texincfile;
         printf("\\smalltabcaption{\\label{table_%s_%s_%s}\n", texsetname, texgroupname, textestname) > texincfile;
         printf("Evaluation of \\setting%s%s on test set \\testset{%s}.}\n", texsetname, texgroupname, textestname) > texincfile;
         printf("\\end{table}\n") > texincfile;
         printf("\n") > texincfile;
      }
   }

   # generate (or extend) tex summary file
   if( texsummaryfile != "" )
   {
      n = split(texsummaryfile, a, "/");
      texsummarypath = "";
      for( i = 1; i < n; i++ )
         texsummarypath = texsummarypath a[i] "/";
      texsummarybase = a[n];
      sub(/\.tex/, "", texsummarybase);
      texsummaryfiletime = texsummarypath texsummarybase "_time.tex";
      texsummaryfilenodes = texsummarypath texsummarybase "_nodes.tex";
      if( texsummaryheader > 0 )
      {
         printf("{\\sffamily\n") > texsummaryfile;
         printf("\\scriptsize\n") > texsummaryfile;
         printf("\\setlength{\\extrarowheight}{1pt}\n") > texsummaryfile;
         printf("\\setlength{\\tabcolsep}{2pt}\n") > texsummaryfile;
         printf("\\newcommand{\\spc}{\\hspace{1em}}\n") > texsummaryfile;
         for( si = 0; si <= 2; si++ )
         {
            printf("\\ifthenelse{\\summaryinfo = %d}{\n", si) > texsummaryfile;
            printf("\\begin{tabular*}{\\columnwidth}{@{}ll@{\\extracolsep{\\fill}}") > texsummaryfile;
            for( o = 1; o < nsolver; o++ )
               printf("r") > texsummaryfile;
            printf("@{}}\n") > texsummaryfile;
            printf("\\toprule\n") > texsummaryfile;
            printf("& test set") > texsummaryfile;
            if( nsolver >= 9 )
            {
               for( o = 1; o < nsolver; o++ )
                  printf(" & %s", texsolvername(printorder[o])) > texsummaryfile;
            }
            else
            {
               for( o = 1; o < nsolver; o++ )
                  printf(" & \\makebox[0em][r]{%s}", texsolvername(printorder[o])) > texsummaryfile;
            }
            printf(" \\\\\n") > texsummaryfile;
            if( si == 0 || si == 1 )
            {
               printf("\\midrule\n") > texsummaryfile;
               printf("\\input{Tables/mip/auto/%s_time}\n", texsummarybase) > texsummaryfile;
            }
            if( si == 0 || si == 2 )
            {
               printf("\\midrule\n") > texsummaryfile;
               printf("\\input{Tables/mip/auto/%s_nodes}\n", texsummarybase) > texsummaryfile;
            }
            printf("\\bottomrule\n") >> texsummaryfile;
            printf("\\end{tabular*}\n") >> texsummaryfile;
            printf("}{}\n") >> texsummaryfile;
         }
         printf("}\n") > texsummaryfile;
         printf("\\raisebox{-%.1fex}[0em][0em]{\\rotatebox{90}{\\makebox[3em]{time}}}",
            1.5*(texsummaryheader+1)) > texsummaryfiletime;
         printf("\\raisebox{-%.1fex}[0em][0em]{\\rotatebox{90}{\\makebox[3em]{nodes}}}",
            1.5*(texsummaryheader+1)) > texsummaryfilenodes;
      }
      printf("& \\testset{%s}", textestset) >> texsummaryfiletime;
      for( o = 1; o < nsolver; o++ )
      {
         s = printorder[o];
         if( texsummaryshifted )
            printf(" & %s", texcompstr(timeshiftedgeom[s,0], reftimeshiftedgeom[s,0])) > texsummaryfiletime;
         else
            printf(" & %s", texcompstr(timegeom[s,0], reftimegeom[s,0])) > texsummaryfiletime;
      }
      printf("\\\\\n") > texsummaryfiletime;
      printf("& \\testset{%s}", textestset) >> texsummaryfilenodes;
      for( o = 1; o < nsolver; o++ )
      {
         s = printorder[o];
         if( texsummaryshifted )
            printf(" & %s", texcompstr(nodeshiftedgeom[s,0], refnodeshiftedgeom[s,0])) > texsummaryfilenodes;
         else
            printf(" & %s", texcompstr(nodegeom[s,0], refnodegeom[s,0])) > texsummaryfilenodes;
      }
      printf("\\\\\n") > texsummaryfilenodes;

      # add tex comment to summary file which is later be used to generate overall statistics
      for( o = 1; o < nsolver; o++ )
      {
         s = printorder[o];
         weight = (texsummaryweight == 0 ? nevalprobs[s,0] : texsummaryweight);
         if( texsummaryshifted )
            printf("%% =mean=  %s %.4f %.4f %g\n", solvername[s], 
               timeshiftedgeom[s,0]/reftimeshiftedgeom[s,0], nodeshiftedgeom[s,0]/refnodeshiftedgeom[s,0],
               weight) >> texsummaryfile;
         else
            printf("%% =mean=  %s %.4f %.4f %g\n", solvername[s], 
               timegeom[s,0]/reftimegeom[s,0], nodegeom[s,0]/refnodegeom[s,0], weight) >> texsummaryfile;
      }
   }

   if( diagramfile != "" )
   {
      refsolver = 0;
      mintime = timeshiftedgeom[0,0];
      maxtime = timeshiftedgeom[0,0];
      for( s = 0; s < nsolver; ++s )
      {
         if( substr(solvername[s], 1, 4) == "SCIP" )
            refsolver = s;

         if( timeshiftedgeom[s,0] < mintime )
            mintime = timeshiftedgeom[s,0];
         else if( timeshiftedgeom[s,0] > mintime )
            maxtime = timeshiftedgeom[s,0];
      }
      # bound, over which solvers get cut off
      upperbound = diagramnsteps*timeshiftedgeom[refsolver,0];

      yscale = 3000/upperbound;

      printf("\\documentclass{article}\n\\usepackage{tikz}\n\\pagestyle{empty}\n\\begin{document}\n\n") > diagramfile;

      printf("\\definecolor{c1}{HTML}{000060}\n") > diagramfile;
      printf("\\definecolor{c2}{HTML}{0000FF}\n") > diagramfile;
      printf("\\definecolor{c3}{HTML}{36648B}\n") > diagramfile;
      printf("\\definecolor{c4}{HTML}{4682B4}\n") > diagramfile;
      printf("\\definecolor{c5}{HTML}{5CACEE}\n") > diagramfile;
      printf("\\definecolor{c6}{HTML}{00FFFF}\n") > diagramfile;
      #printf("\\definecolor{c7}{HTML}{A0FFFF}\n") > diagramfile;
      printf("\\definecolor{c7}{HTML}{008888}\n") > diagramfile;
      printf("\\definecolor{c8}{HTML}{00DD99}\n") > diagramfile;
      printf("\\definecolor{c9}{HTML}{527B10}\n") > diagramfile;
      printf("\\definecolor{c10}{HTML}{7BC618}\n") > diagramfile;
      printf("\\definecolor{c11}{HTML}{A00060}\n") > diagramfile;
      printf("\\definecolor{c12}{HTML}{A000FF}\n") > diagramfile;
      printf("\\definecolor{c13}{HTML}{E6648B}\n") > diagramfile;
      printf("\\definecolor{c14}{HTML}{F682B4}\n") > diagramfile;
      printf("\\definecolor{c15}{HTML}{DCACEE}\n") > diagramfile;
      printf("\\definecolor{c16}{HTML}{A0FFFF}\n") > diagramfile;
      #printf("\\definecolor{c7}{HTML}{A0FFFF}\n") > diagramfile;
      printf("\\definecolor{c17}{HTML}{A08888}\n") > diagramfile;
      printf("\\definecolor{c18}{HTML}{A0DD99}\n") > diagramfile;
      printf("\\definecolor{c19}{HTML}{D27B10}\n") > diagramfile;
      printf("\\definecolor{c20}{HTML}{FBC618}\n") > diagramfile;

      printf("\\definecolor{darkgreen}{HTML}{006600}\n") > diagramfile;
      if( diagramyellowbg )
        printf("\\definecolor{background}{HTML}{FFFFE6}\n\n") > diagramfile;
      else
        printf("\\definecolor{background}{HTML}{FFFFFF}\n\n") > diagramfile;

      printf("\\begin{tikzpicture}[auto,scale=0.8,yscale=%1.2f]\n",yscale) > diagramfile;
      printf("\n%% tikz styles\n") > diagramfile;
      printf("\\tikzstyle{box} +=[draw=black,rectangle,inner sep=0mm, minimum size = 2.5mm,left];\n") > diagramfile;
      printf("\\tikzstyle{legend} +=[inner sep=1mm, right];\n") > diagramfile;
      printf("\\tikzstyle{tic} +=[left];\n") > diagramfile;
      printf("\\tikzstyle{bel} +=[below,inner sep=0.8mm];\n") > diagramfile;
      printf("\\tikzstyle{abo} +=[above,inner sep=0.3mm];\n\n") > diagramfile;

      #extendedub = (1+1/(2*diagramnsteps-2))*upperbound;
      extendedub = 1.2*upperbound;
      printf("\\draw[background,fill=background] (-1.5,-0.3) rectangle (%1.2f,%1.2f);\n",nsolver+6.5,extendedub/1000.0+0.1) > diagramfile;

      printf("%% Beschriftungen und Hoehenlinien\n") > diagramfile;
      for( i = 0; i < diagramnsteps; ++i )
      {
         perc = i/(diagramnsteps-1.0)*upperbound/1000.0;
         printf("\\node () at (0.2,%1.2f) [tic] {\\small %d};\n",perc,perc*1000) > diagramfile;
         printf("\\draw[dotted] (0.2,%1.2f) -- (%1.2f,%1.2f);\n",perc,nsolver+0.8,perc) > diagramfile;
      }
      if( diagramgerman )
         printf("\\node () at (-1,%1.2f) [rotate=90]{\\footnotesize Zeit in Sekunden};\n",extendedub/2000.0) > diagramfile;
      else
         printf("\\node () at (-1,%1.2f) [rotate=90]{\\footnotesize time in seconds};\n",extendedub/2000.0) > diagramfile;
      printf("\\draw[] (0.2,0) rectangle (%1.2f,%1.2f);\n\n",nsolver+0.8,extendedub/1000.0) > diagramfile;

      printf("%% BALKEN\n") > diagramfile;
      for( i = 1; i <= nsolver; ++i )
      {
         pos[i] = i-1;
         if( timeshiftedgeom[pos[i],0] >= upperbound)
            printf("\\draw[fill=c%d] (%d.5,0) rectangle (%d.5,%1.3f);\n",i,i-1,i,extendedub/1000.0) > diagramfile;
         else
            printf("\\draw[fill=c%d] (%d.5,0) rectangle (%d.5,%1.3f);\n",i,i-1,i,timeshiftedgeom[pos[i],0]/1000.0) > diagramfile;
      }
      printf("\n") > diagramfile;

      printf("%% TRENNLINIEN\n") > diagramfile;
      upperbreak = 0.92 * extendedub/1000.0;
      lowerbreak = 0.90 * extendedub/1000.0;
      printf("\\fill[background] (0.1,%1.4f) rectangle (%d.5,%1.4f);\n",lowerbreak,nsolver,upperbreak) > diagramfile;
      printf("\\draw[black] (0.1,%1.4f) --  (0.3,%1.4f);\n",lowerbreak,lowerbreak) > diagramfile;
      printf("\\draw[black] (0.1,%1.4f) --  (0.3,%1.4f);\n",upperbreak,upperbreak) > diagramfile;

      printf("%% BALKENBESCHRIFTUNG \n") > diagramfile;
      printf("\\draw[c%d] (0.2,%1.3f) -- (%1.2f,%1.3f);\n",refsolver+1,timeshiftedgeom[refsolver,0]/1000.0,
         nsolver+0.8,timeshiftedgeom[refsolver,0]/1000.0) > diagramfile;

      for( i = 1; i <= nsolver; ++i )
      {
        if( timeshiftedgeom[pos[i],0] >= upperbound )
           printf("\\node () at (%d,%1.3f) [bel,inner sep=0.3mm] {\\footnotesize\\textcolor{white}{%2.1fx}};\n",i,
              extendedub/1000.0,timeshiftedgeom[pos[i],0]/timeshiftedgeom[refsolver,0]) > diagramfile;
        else if( pos[i] == refsolver )
           printf("\\node () at (%d,%1.3f) [abo] {\\footnotesize\\textbf{%1.2fx}};\n",i,
              timeshiftedgeom[pos[i],0]/1000.0,timeshiftedgeom[pos[i],0]/timeshiftedgeom[refsolver,0]) > diagramfile;
        else
           printf("\\node () at (%d,%1.3f) [abo] {\\footnotesize %1.2fx};\n",i,
              timeshiftedgeom[pos[i],0]/1000.0,timeshiftedgeom[pos[i],0]/timeshiftedgeom[refsolver,0]) > diagramfile;
      }
      printf("\n") > diagramfile;

      printf("%% NOT SOLVED\n") > diagramfile;
      if( diagramgerman )
         printf("\\node () at (-0.4, %1.3f) [bel]{\\footnotesize\\textcolor{blue}{Nicht gel\\\"ost}};\n", -upperbound/6000.0) > diagramfile;
      else
         printf("\\node () at (-0.4, %1.3f) [bel]{\\footnotesize\\textcolor{blue}{not solved}};\n", -upperbound/6000.0) > diagramfile;
      for( i = 1; i <= nsolver; ++i )
         printf("\\node () at (%d, %1.3f) [bel]{\\footnotesize\\textcolor{blue}{%2.0f\\%}};\n",i,-upperbound/6000.0,100.0-100.0*((nsolved[pos[i],0]+0.0)/(nprocessedprobs[pos[i],0]+0.0))) > diagramfile;
      printf("\n") > diagramfile;

      printf("%% LEGEND\n") > diagramfile;
      for( i = 1; i <= nsolver; ++i )
      {
         perc = (nsolver-i)/(nsolver-1.0)*extendedub/1000.0;
         if( pos[i] == refsolver )
            printf("\\node () at (%1.4f,%1.4f) [legend]{\\small \\textbf{%s}};\n",nsolver+1.4,perc,solvername[pos[i]]) > diagramfile;
         else
            printf("\\node () at (%1.4f,%1.4f) [legend]{\\small %s};\n",nsolver+1.4,perc,solvername[pos[i]]) > diagramfile;
         printf("\\node () at (%1.4f,%1.4f) [box,color=c%d,fill=c%d]{};\n",nsolver+1.4,perc,i,i) > diagramfile;
      }
      printf("\\node () at (%1.4f,%1.3f) [legend]{\\footnotesize timings w.r.t. %d instances};\n",nsolver+1.4,-upperbound/6000.0,nevalprobs[printorder[0],0]) > diagramfile;
      printf("\\end{tikzpicture}\n") > diagramfile;
      printf("\n\\end{document}\n") > diagramfile;
   }
}
