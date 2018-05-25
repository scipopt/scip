#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    permcmpresall.awk
#@brief   compare different versions of runs with permuations
#@author  Jan Kuske
#@author  Marc Pfetsch

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

function printhline(nsolver,short)
{
   for (s = 0; s < nsolver; ++s)
   {
      if ( s == 0 )
      {
	 printf("--------------------+-+---------+--------+");
	 if (printdifftime)
	    printf("--------+");
      }
      else
      {
	 if ( !short )
	 {
	    if (!printdifftime)
	       printf("-+---------+--------+------+------+");
	    else
	       printf("-+---------+--------+--------+------+------+------+");
	 }
	 else
	 {
	    printf("-+---------+--------+");
	 }
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


BEGIN {
   short = 0;           # for each non reference solver, only absolute time and number of nodes are printed
   printgap = 0;        # if timeout, then print absolute gap at termination in time column, if gap is finite
   printdifftime=0;     # additional column with the difference of maximum time and minimum time

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
   if ( setorder[$2] == 0 )
   {
      nsetnames++;
      setorder[$2] = nsetnames;
      setname[$2] = $3;
      for (i = 4; i <= NF; i++)
	 setname[$2] = setname[$2]" "$i;
   }
   setingroup[$2,group] = 1;
}

/^@03 / {
   permutations[nsolver] = $3;
}

/^@02 timelimit: / {
   timelimit[nsolver] = $3;
}

/^@01 / {
   if ( onlygroup == 0 || setingroup[$2,onlygroup] )
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

   validline = 0;
   # check if this is a useable line

   if ( $14 in statuses ) # SCIP with max-min time difference
   {
      # collect data (line with problem type, original and presolved problem size and simplex iterations)
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($8, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $9;
      iters[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($12,mintime),0.1);
      timediff[nsolver,nprobs[nsolver]] = fracceil(max($13,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $14;
      statcount[nsolver,nprobs[nsolver]] = $15;
      validline = 1;
   }

   if ( $15 in statuses ) # SCIP with median or stddtime (both will not be used)
   {
      # collect data (line with problem type, original and presolved problem size and simplex iterations)
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($8, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $9;
      iters[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($12,mintime),0.1);
      timediff[nsolver,nprobs[nsolver]] = fracceil(max($13,mintime),0.1);
      status[nsolver,nprobs[nsolver]] = $15;
      statcount[nsolver,nprobs[nsolver]] = $16;
      validline = 1;
   }

   if ( $16 in statuses ) # SCIP with median and stddtime
   {
      # collect data (line with problem type, original and presolved problem size and simplex iterations)
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = max(min($7, +infinity), -infinity);
      primalbound[nsolver,nprobs[nsolver]] = max(min($8, +infinity), -infinity);
      gap[nsolver,nprobs[nsolver]] = $9;
      iters[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = fracceil(max($12,mintime),0.1);
      timediff[nsolver,nprobs[nsolver]] = fracceil(max($13,mintime),0.1);
      stdd[nsolver,nprobs[nsolver]] = $14;
      median[nsolver,nprobs[nsolver]] = $15;
      status[nsolver,nprobs[nsolver]] = $16;
      statcount[nsolver,nprobs[nsolver]] = $17;
      validline = 1;
   }

   if ( validline )
   {
      # postprocessing of information
      if ( status[nsolver,nprobs[nsolver]] == "better" )
	 status[nsolver,nprobs[nsolver]] = "timeout";
      if ( status[nsolver,nprobs[nsolver]] == "sollimit" || status[nsolver,nprobs[nsolver]] == "gaplimit" || status[nsolver,nprobs[nsolver]] == "solved" )
	 status[nsolver,nprobs[nsolver]] = "ok";

      if ( status[nsolver,nprobs[nsolver]] == "timeout" || status[nsolver,nprobs[nsolver]] == "nodelimit" ||  status[nsolver,nprobs[nsolver]] == "memlimit")
	 hitlimit[nsolver,nprobs[nsolver]] = 1;
      else
	 hitlimit[nsolver,nprobs[nsolver]] = 0;

      probidx[$1,nsolver] = nprobs[nsolver];
      probcnt[$1]++;
      nprobs[nsolver]++;

      if ( probcnt[$1] == 1 )
      {
	 problist[problistlen] = $1;
	 problistlen++;
      }
   }
}

END {
   printdifftime = !short && printdifftime;   # if short=1 set printdifftime=0

   if ( onlygroup > 0 && nsolver == 1 && solvername[1] == "SCIP:default" )
   {
      printf("only SCIP:default setting found.\n");
      exit 1;
   }

   if ( nsolver == 0 )
   {
      printf("No instances found in log files.\n");
      exit 1;
   }

   # process exclude string
   n = split(exclude, a, ",");
   for (i = 1; i <= n; i++)
      excluded[a[i]] = 1;

   permutation = permutations[0];

   # initialize means
   for (s = 0; s < nsolver; ++s)
   {
      if ( permutation != permutations[s] )
      {
	 printf("Permutations are not equal between runs.\n");
	 exit 1;
      }

      # cat: -1 - all within time limit, 0 - all, 1 - different path, 2 - equal path, 3 - all timeout
      for (cat = -1; cat <= 3; cat++)
      {
	 nevalprobs[s,cat] = 0;
	 nprocessedprobs[s,cat] = 0;

	 timetotal[s,cat] = 0.0;
	 mmttotal[s,cat] = 0.0;
	 nodetotal[s,cat] = 0.0;

	 timegeom[s,cat] = 1.0;
	 mmtgeom[s,cat] = 1.0;
	 nodegeom[s,cat] = 1.0;

	 timeshiftedgeom[s,cat] = timegeomshift;
	 mmtshiftedgeom[s,cat] = timegeomshift;
	 nodeshiftedgeom[s,cat] = nodegeomshift;

	 reftimetotal[s,cat] = 0.0;
	 refnodetotal[s,cat] = 0.0;

	 reftimegeom[s,cat] = 1.0;
	 refnodegeom[s,cat] = 1.0;

	 reftimeshiftedgeom[s,cat] = timegeomshift;
	 refnodeshiftedgeom[s,cat] = nodegeomshift;

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
   for (s = 0; s < nsolver; ++s)
   {
      sname = solvername[s];
      for (o = 0; o < s; ++o)
      {
	 i = printorder[o];
	 iname = solvername[i];
	 if ( nsetnames > 0 )
	 {
	    # use order given by =setname= entries
	    if ( setorder[sname] < setorder[iname] )
	       break;
	 }
	 else
	 {
	    # use alphabetical order, but put CPLEX before SCIP and "default" before all others
	    if ( substr(sname, 1, 5) == "CPLEX" && substr(iname, 1, 5) != "CPLEX" )
	       break;
	    if ( substr(sname, 1, 5) == substr(iname, 1, 5) && match(sname, "default") != 0 && match(iname, "default") == 0 )
	       break;
	    if ( substr(sname, 1, 5) == substr(iname, 1, 5) && (match(sname, "default") == 0) == (match(iname, "default") == 0) && sname < iname )
	       break;
	 }
      }
      for (j = s-1; j >= o; --j)
	 printorder[j+1] = printorder[j];
      printorder[o] = s;
   }

   # print headers
   for (o = 0; o < nsolver; ++o)
   {
      s = printorder[o];
      sname = solvername[s];

      if ( o == 0 )
      {
	 if ( length(sname) <= 39 )
	    printf(" %39s |", sname)
	 else
	    printf(" *%38s |", substr(sname, length(sname)-39));
      }
      else
      {
	 if ( short )
	 {
	    if ( length(sname) <= 19 )
	       printf("%19s |", sname);
	    else
	       printf("*%18s |", substr(sname, length(sname)-19));
	 }
	 else
	 {
	    if ( length(sname) <= 33 )
	       printf("%31s |", sname);
	    else
	       printf("*%30s |", substr(sname, length(sname)-31));
	 }
      }
   }
   printf("\n");
   printhline(nsolver,short);#, printsoltimes);

   printf("  Name              |");
   for (s = 0; s < nsolver; ++s)
   {
      if ( s == 0 || short )
      {
	 if ( !printdifftime )
	    printf("F|   Nodes |   Time |");
	 else
	    printf("F|   Nodes |   Time |  DifT  |");
      }
      else
      {
	 if ( ! printdifftime )
	    printf("F|   Nodes |   Time | NodQ | TimQ |");
	 else
	    printf("F|   Nodes |   Time |  DifT  | NodQ | TimQ | DifQ |");
      }
   }
   printf(" bounds check\n");
   printhline(nsolver,short);#, printsoltimes);

   # display the problem results and calculate mean values
   for (i = 0; i < problistlen; ++i)
   {
      p = problist[i];
      if ( length(p) > 18 )
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
      if ( excluded[p] )
      {
	 unprocessed = 1;
	 countprob = 0;
      }

      # find best and worst run and check whether this instance should be counted in overall statistics
      for (s = 0; s < nsolver; ++s)
      {
	 pidx = probidx[p,s];
	 processed = (pidx != "");

	 # make sure, nodes and time are non-zero for geometric means
	 nodes[s,pidx] = max(nodes[s,pidx], 1);
	 time[s,pidx] = max(time[s,pidx], mintime);
	 fulltotaltime += time[s,pidx];

	 # if the solver exceeded the timelimit, set status accordingly
	 if ( (status[s,pidx] == "ok" || status[s,pidx] == "unknown") && timelimit[s] > 0.0 && time[s,pidx] > timelimit[s] )
	 {
	    status[s,pidx] = "timeout";
	    time[s,pidx] = timelimit[s];
	 }

	 # check if all solvers processed the problem
	 if ( !processed )
	 {
	    marker = "?";
	    unprocessed = 1;
	 }

	 # check if solver ran successfully (i.e., no abort nor fail)
	 if ( processed && (status[s,pidx] == "ok" || status[s,pidx] == "unknown" || status[s,pidx] == "timeout" || status[s,pidx] == "nodelimit" || status[s,pidx] == "memlimit") )
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
      for (o = 0; o < nsolver; ++o)
      {
	 s = printorder[o];
	 pidx = probidx[p,s];
	 processed = (pidx != "");

	 if ( ! processed )
	 {
	    notimeout = 0;
	    continue;
	 }

	 if ( nodecomp == -1 )
	 {
	    itercomp = iters[s,pidx];
	    nodecomp = nodes[s,pidx];
	    timecomp = time[s,pidx];
	    maxmindiffcomp = timediff[s,pidx];
	    timeoutcomp = hitlimit[s,pidx];
	    timetofirstcomp = max(mintime, timetofirst[s,pidx]);
	    timetobestcomp = max(mintime, timetobest[s,pidx]);
	 }
	 iseqpath = (iters[s,pidx] == itercomp && nodes[s,pidx] == nodecomp);
	 hastimeout = timeoutcomp || hitlimit[s,pidx];
	 notimeout = notimeout && !timeoutcomp &&  !hitlimit[s,pidx];

	 # which category?
	 if ( hastimeout )
	    category[s] = 3;
	 else if ( iseqpath )
	    category[s] = 2;
	 else
	    category[s] = 1;
      }

      # evaluate instance for all solvers
      for (o = 0; o < nsolver; ++o)
      {
	 s = printorder[o];
	 pidx = probidx[p,s];
	 processed = (pidx != "");

	 if ( processed && name[s,pidx] != p )
	    printf("Error: solver %s, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

	 # check if solver ran successfully (i.e., no abort nor fail)
	 if ( processed )
	 {
	    if ( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
	    {
	       marker = " ";
	       if ( ! unprocessed )
	       {
		  if ( notimeout )
		     nsolved[s,-1]++;
		  nsolved[s,0]++;
		  nsolved[s,category[s]]++;
		  nthissolved++;
	       }
	    }
	    else if ( hitlimit[s,pidx] )
	    {
	       marker = ">";
	       notimeout = 0;
	       if ( ! unprocessed )
	       {
		  if ( countprob )
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
	       if ( status[s,pidx] == "readerror" )
		  readerror = 1;
	       if ( !unprocessed )
	       {
		  nfails[s,0]++;
		  nfails[s,category[s]]++;
		  nthisfails++;
	       }
	    }
	 }

	 if ( primalbound[s,pidx] < infinity )
	    feasmark = " ";
	 else
	    feasmark = "#";

	 if ( processed && !fail )
	 {
	    mindb = min(mindb, dualbound[s,pidx]);
	    maxdb = max(maxdb, dualbound[s,pidx]);
	    minpb = min(minpb, primalbound[s,pidx]);
	    maxpb = max(maxpb, primalbound[s,pidx]);
	    ismini = ismini || (primalbound[s,pidx] > dualbound[s,pidx] + 1e-06);
	    ismaxi = ismaxi || (primalbound[s,pidx] < dualbound[s,pidx] - 1e-06);
	 }

	 # print statistics
	 if ( !processed )
	    line = sprintf("%s           -        -", line);
	 else
	 {
	    if (!printdifftime)
	       line = sprintf("%s %s%10d %s%7.1f", line, feasmark, nodes[s,pidx], marker, time[s,pidx]);
	    else
	       line = sprintf("%s %s%10d %s%7.1f %8.1f", line, feasmark, nodes[s,pidx], marker, time[s,pidx],timediff[s,pidx]);
	 }

	 if ( o > 0 && !short )
	 {
	    if ( !processed )
	       line = sprintf("%s      -", line);
	    else if ( nodes[s,pidx]/nodecomp > 999.99 )
	       line = sprintf("%s  Large", line);
	    else
	       line = sprintf("%s %6.2f", line, nodes[s,pidx]/nodecomp);
	    if ( ! processed )
	       line = sprintf("%s      -", line);
	    else if ( time[s,pidx]/timecomp > 999.99 )
	       line = sprintf("%s  Large", line);
	    else
	       line = sprintf("%s %6.2f", line, time[s,pidx]/timecomp);
	    if ( ! processed && printdifftime )
	       line = sprintf("%s      -", line);
	    else if ( timediff[s,pidx]/maxmindiffcomp > 999.99 && printdifftime )
	       line = sprintf("%s  Large", line);
	    else if ( printdifftime )
	       line = sprintf("%s %6.2f", line, timediff[s,pidx]/maxmindiffcomp);

	    if ( processed &&
	       (timeoutcomp != hitlimit[s,pidx] ||
		nodes[s,pidx] > markworsenodes * nodecomp ||
		nodes[s,pidx] < 1.0/markbetternodes * nodecomp ||
		isfaster(time[s,pidx], timecomp, markbettertime) ||
		isslower(time[s,pidx], timecomp, markworsetime)) )
	       mark = "*";
	 }
      }

      # update the best status information
      if ( nthissolved > 0 )
	 bestnsolved++;
      else if ( nthistimeouts > 0 )
	 bestntimeouts++;
      else if ( nthisfails == nsolver - nthisunprocessed )
	 bestnfails++;

      # check for inconsistency in the primal and dual bounds
      if ( readerror )
      {
	 line = sprintf("%s  readerror", line);
	 mark = " ";
      }
      else if ( fail )
      {
	 line = sprintf("%s  fail", line);
	 mark = " ";
      }
      else if ( consistency &&
       ((ismini && ismaxi) ||
	  (ismini && maxdb - minpb > 1e-5 * max(max(abs(maxdb), abs(minpb)), 1.0)) ||
	  (ismaxi && maxpb - mindb > 1e-5 * max(max(abs(maxpb), abs(mindb)), 1.0)) ||
	  (!ismini && !ismaxi && abs(maxpb - minpb) > 1e-5 * max(abs(maxpb), 1.0))) )
      {
	 line = sprintf("%s  inconsistent", line);
	 fail = 1;
	 mark = " ";
      }
      else if ( excluded[p] )
      {
	 line = sprintf("%s  excluded", line);
	 mark = " ";
      }
      else if ( unprocessed )
      {
	 line = sprintf("%s  unprocessed", line);
	 mark = " ";
      }
      else
	 line = sprintf("%s  ok", line);

      # calculate number of instances for which feasible solution has been found
      hasfeasible = 0;
      if ( ! unprocessed )
      {
	 for ( s = 0; s < nsolver; ++s )
	 {
	    if ( notimeout )
	       nprocessedprobs[s,-1]++;
	    nprocessedprobs[s,0]++;
	    nprocessedprobs[s,category[s]]++;
	    pidx = probidx[p,s];
	    if ( primalbound[s,pidx] < infinity )
	    {
	       if ( notimeout )
		  feasibles[s,-1]++;
	       feasibles[s,0]++;
	       feasibles[s,category[s]]++;
	       hasfeasible = 1;
	    }
	 }
	 if ( hasfeasible )
	    bestfeasibles++;
      }

      if ( (!onlymarked || mark == "*") && (!onlyprocessed || !unprocessed) &&
	(!onlyfeasible || hasfeasible) && (!onlyinfeasible || !hasfeasible) &&
	(!onlyfail || fail) )
      {
	 printf("%s %s\n", mark, line);
      }

      # calculate totals and means for instances where no solver failed
      if ( ! fail && ! unprocessed && (! onlyfeasible || hasfeasible) && (! onlyinfeasible || ! hasfeasible) )
      {
	 reftime = time[printorder[0],probidx[p,printorder[0]]];
	 refhitlimit = hitlimit[printorder[0],probidx[p,printorder[0]]];
	 refnodes = nodes[printorder[0],probidx[p,printorder[0]]];
	 refobj = primalbound[printorder[0],probidx[p,printorder[0]]];
	 reftimetofirst = timetofirst[printorder[0],probidx[p,printorder[0]]];
	 reftimetobest = timetobest[printorder[0],probidx[p,printorder[0]]];
	 hasbetter = 0;
	 hasbetterobj = 0;
	 for (s = 0; s < nsolver; ++s)
	 {
	    pidx = probidx[p,s];
	    for (cat = 0; cat <= 3; cat = 3*cat + category[s])
	    {
	       nevalprobs[s,cat]++;
	       nep = nevalprobs[s,cat];
	       timetotal[s,cat] += time[s,pidx];
	       mmttotal[s,cat] += timediff[s,pidx];
	       timetofirsttotal[s,cat] += timetofirst[s,pidx];
	       timetobesttotal[s, cat] += timetobest[s, pidx];
	       nodetotal[s,cat] += nodes[s,pidx];
	       timegeom[s,cat] = timegeom[s,cat]^((nep-1)/nep) * time[s,pidx]^(1.0/nep);
	       mmtgeom[s,cat] = mmtgeom[s,cat]^((nep-1)/nep) * timediff[s,pidx]^(1.0/nep);
	       timetofirstgeom[s,cat] = timetofirstgeom[s,cat]^((nep-1)/nep) * max(timetofirst[s,pidx], mintime)^(1.0/nep);
	       timetobestgeom[s,cat] = timetobestgeom[s,cat]^((nep-1)/nep) * max(timetobest[s,pidx])^(1.0/nep);
	       nodegeom[s,cat] = nodegeom[s,cat]^((nep-1)/nep) * nodes[s,pidx]^(1.0/nep);
	       timeshiftedgeom[s,cat] = timeshiftedgeom[s,cat]^((nep-1)/nep) * (time[s,pidx]+timegeomshift)^(1.0/nep);
	       mmtshiftedgeom[s,cat] = mmtshiftedgeom[s,cat]^((nep-1)/nep) * (timediff[s,pidx]+timegeomshift)^(1.0/nep);
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
	       if ( time[s,pidx] <= wintolerance*besttime )
		  wins[s,cat]++;
	       if ( !hitlimit[s,pidx] && (isfaster(time[s,pidx], reftime, wintolerance) || refhitlimit))
	       {
		  better[s,cat]++;
		  hasbetter = 1;
	       }
	       else if ( !refhitlimit && (isslower(time[s,pidx], reftime, wintolerance) || (hitlimit[s,pidx])))
		  worse[s,cat]++;
	       pb = primalbound[s,pidx];
	       if ( (ismini && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) || (ismaxi && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ) {
		  betterobj[s,cat]++;
		  hasbetterobj = 1;
	       }
	       else if ( (ismini && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) || (ismaxi && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
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
	 if ( hasbetter )
	    bestbetter++;
	 if ( hasbetterobj )
	    bestbetterobj++;

	 # once again for the case in which all instances have been solved to optimality
	 if ( notimeout )
	 {
	    for (s = 0; s < nsolver; ++s)
	    {
	       pidx = probidx[p,s];
	       cat = -1;
	       nevalprobs[s,cat]++;
	       nep = nevalprobs[s,cat];
	       timetotal[s,cat] += time[s,pidx];
	       mmttotal[s,cat] += timediff[s,pidx];
	       nodetotal[s,cat] += nodes[s,pidx];
	       timegeom[s,cat] = timegeom[s,cat]^((nep-1)/nep) * time[s,pidx]^(1.0/nep);
	       mmtgeom[s,cat] = mmtgeom[s,cat]^((nep-1)/nep) * timediff[s,pidx]^(1.0/nep);
	       nodegeom[s,cat] = nodegeom[s,cat]^((nep-1)/nep) * nodes[s,pidx]^(1.0/nep);
	       timeshiftedgeom[s,cat] = timeshiftedgeom[s,cat]^((nep-1)/nep) * (time[s,pidx]+timegeomshift)^(1.0/nep);
	       mmtshiftedgeom[s,cat] = mmtshiftedgeom[s,cat]^((nep-1)/nep) * (timediff[s,pidx]+timegeomshift)^(1.0/nep);
	       nodeshiftedgeom[s,cat] = nodeshiftedgeom[s,cat]^((nep-1)/nep) * (nodes[s,pidx]+nodegeomshift)^(1.0/nep);
	       reftimetotal[s,cat] += reftime;
	       refnodetotal[s,cat] += refnodes;
	       reftimegeom[s,cat] = reftimegeom[s,cat]^((nep-1)/nep) * reftime^(1.0/nep);
	       refnodegeom[s,cat] = refnodegeom[s,cat]^((nep-1)/nep) * refnodes^(1.0/nep);
	       reftimeshiftedgeom[s,cat] = reftimeshiftedgeom[s,cat]^((nep-1)/nep) * (reftime+timegeomshift)^(1.0/nep);
	       refnodeshiftedgeom[s,cat] = refnodeshiftedgeom[s,cat]^((nep-1)/nep) * (refnodes+nodegeomshift)^(1.0/nep);
	       if ( time[s,pidx] <= wintolerance*besttime )
		  wins[s,cat]++;
	       if ( isfaster(time[s,pidx], reftime, wintolerance) )
		  better[s,cat]++;
	       else if ( isslower(time[s,pidx], reftime, wintolerance) )
		  worse[s,cat]++;
	       pb = primalbound[s,pidx];
	       if ( (ismini && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) || (ismaxi && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
	       {
		  betterobj[s,cat]++;
	       }
	       else if ( (ismini && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) || (ismaxi && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
		  worseobj[s,cat]++;

	       thisscore = reftime/time[s,pidx];
	       thisscore = max(thisscore, 1/maxscore);
	       thisscore = min(thisscore, maxscore);
	       score[s,cat] = score[s,cat]^((nep-1)/nep) * thisscore^(1.0/nep);
	    }
	 }
      }
   }
   printhline(nsolver,short);#, printsoltimes);

   # make sure total time and nodes is not zero
   for (s = 0; s < nsolver; ++s)
   {
      for (cat = -1; cat <= 3; cat++)
      {
	 nodetotal[s,cat] = max(nodetotal[s,cat], 1);
	 refnodetotal[s,cat] = max(refnodetotal[s,cat], 1);
	 timetotal[s,cat] = max(timetotal[s,cat], mintime);
	 mmttotal[s,cat] = max(mmttotal[s,cat], mintime);
	 reftimetotal[s,cat] = max(reftimetotal[s,cat], mintime);
      }
   }

   # print solvers' overall statistics
   probnumstr = "("nevalprobs[printorder[0],0]")";
   printf("%-14s %5s", "total", probnumstr);
   for (o = 0; o < nsolver; ++o)
   {
      s = printorder[o];
      if ( o == 0 || short )
      {
	 printf(" %11d %8d", nodetotal[s,0], timetotal[s,0]);
	 if (printdifftime)
	    printf(" %8d",mmttotal[s,0]);
      }
      else
      {
	 if (! printdifftime)
	    printf(" %11d %8d              ", nodetotal[s,0], timetotal[s,0]);
	 else
	    printf(" %11d %8d %8d                     ", nodetotal[s,0], timetotal[s,0],mmttotal[s,0]);
      }
   }
   printf("\n");
   printf("%-20s", "geom. mean");

   nodegeomcomp = -1;
   timegeomcomp = -1;
   mmtgeomcomp = -1;
   nodetotalcomp = -1;
   timetotalcomp = -1;
   mmttotalcomp = -1;

   for (o = 0; o < nsolver; ++o)
   {
      s = printorder[o];
      if ( o == 0 || short )
      {
	 printf(" %11d %8.1f", nodegeom[s,0], timegeom[s,0]);
	 if ( printdifftime )
	    printf(" %8.1f",mmtgeom[s,0])

	 if ( nodegeomcomp < 0 )
	    nodegeomcomp = nodegeom[s,0];
	 if ( timegeomcomp < 0 )
	    timegeomcomp = timegeom[s,0];
	 if ( mmtgeomcomp < 0 )
	    mmtgeomcomp = mmtgeom[s,0];
	 if ( nodetotalcomp < 0 )
	    nodetotalcomp = nodetotal[s,0];
	 if ( timetotalcomp < 0 )
	    timetotalcomp = timetotal[s,0];
	 if ( mmttotalcomp < 0 )
	    mmttotalcomp = mmttotal[s,0];
      }
      else
      {
	 if (!printdifftime)
	    printf(" %11d %8.1f %6.2f %6.2f", nodegeom[s,0], timegeom[s,0], nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp);
	 else
	    printf(" %11d %8.1f %8.1f %6.2f %6.2f %6.2f", nodegeom[s,0], timegeom[s,0], mmtgeom[s,0], nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp, mmtgeom[s,0]/mmtgeomcomp);
      }
   }
   printf("\n");
   printf("%-20s", "shifted geom.");
   nodeshiftedgeomcomp = -1;
   timeshiftedgeomcomp = -1;
   mmtshiftedgeomcomp = -1;

   for (o = 0; o < nsolver; ++o)
   {
      s = printorder[o];
      for (cat = -1; cat <= 3; cat++)
      {
	 nodeshiftedgeom[s,cat] -= nodegeomshift;
	 timeshiftedgeom[s,cat] -= timegeomshift;
	 mmtshiftedgeom[s,cat] -= timegeomshift;
	 nodeshiftedgeom[s,cat] = max(nodeshiftedgeom[s,cat], 1);
	 timeshiftedgeom[s,cat] = max(timeshiftedgeom[s,cat], mintime);
	 mmtshiftedgeom[s,cat] = max(mmtshiftedgeom[s,cat], mintime);
	 refnodeshiftedgeom[s,cat] -= nodegeomshift;
	 reftimeshiftedgeom[s,cat] -= timegeomshift;
	 refnodeshiftedgeom[s,cat] = max(refnodeshiftedgeom[s,cat], mintime);
	 reftimeshiftedgeom[s,cat] = max(reftimeshiftedgeom[s,cat], mintime);
      }

      if ( o == 0 || short )
      {
	 printf(" %11d %8.1f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0]);

	 if (printdifftime)
	   printf(" %8.1f", mmtshiftedgeom[s,0]);

	 if ( nodeshiftedgeomcomp < 0 )
	    nodeshiftedgeomcomp = nodeshiftedgeom[s,0];
	 if ( timeshiftedgeomcomp < 0 )
	    timeshiftedgeomcomp = timeshiftedgeom[s,0];
	 if ( mmtshiftedgeomcomp < 0 )
	    mmtshiftedgeomcomp = mmtshiftedgeom[s,0];
      }
      else
      {
	 if (!printdifftime)
	    printf(" %11d %8.1f %6.2f %6.2f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0], nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp);
	 else
	    printf(" %11d %8.1f %8.1f %6.2f %6.2f %6.2f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0], mmtshiftedgeom[s,0],
	       nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp , mmtshiftedgeom[s,0]/mmtshiftedgeomcomp);
      }
   }
   bestnodeshiftedgeom -= nodegeomshift;
   besttimeshiftedgeom -= timegeomshift;
   bestnodeshiftedgeom = max(bestnodeshiftedgeom, 1.0);
   besttimeshiftedgeom = max(besttimeshiftedgeom, 1.0);

   printf("\n");

   #since the rows of the quotients are not printed, print the quotients of the geometric means
   if ( short )
   {
      printf("quot. geom. mean                         ");
      for (o = 0; o < nsolver; ++o)
      {
	 if ( o > 0 )
	 {
	    s = printorder[o];
	    printf("      %6.2f   %6.2f",nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp);
	 }
      }
      printf("\n");
      printf("quot. sh. geom. mean                     ");
      for (o = 0; o < nsolver; ++o)
      {
	 if ( o > 0 )
	 {
	    s = printorder[o];
	    printf("      %6.2f   %6.2f",nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp);
	 }
      }
      printf("\n");
      printf("percent not solved                 ");
      for (o = 0; o < nsolver; ++o)
      {
	 s = printorder[o];
	 printf("%6.2f",100-100*nsolved[s,0]/nprocessedprobs[s,0]);
	 printf("               ");
      }
      printf("\n");
   }

   # print number of permutations
   permstr = "("permutations[0]")";
   printf("%-14s %5s\n", "permutations", permstr);
   printhline(nsolver,short);#, printsoltimes);

   if ( ! short )
   {
      for (cat = 0; cat <= 3; cat++)
      {

	 header = (cat == -1 ? "optimal" : (cat == 0 ? "all" : (cat == 1 ? "diff" : (cat == 2 ? "equal" : "timeout"))));
	 printf("\n");
	 printf("%-7s                                            proc eval fail time solv wins bett wors bobj wobj feas    gnodes   shnodes   gnodesQ  shnodesQ   gtime  shtime  gtimeQ shtimeQ   score\n", header);

	 for (o = 0; o < nsolver; ++o)
	 {
	    s = printorder[o];
	    sname = solvername[s];
	    if ( o == 0 )
	    {
	       nodegeomcomp = nodegeom[s,cat];
	       timegeomcomp = timegeom[s,cat];
	       nodeshiftedgeomcomp = nodeshiftedgeom[s,cat];
	       timeshiftedgeomcomp = timeshiftedgeom[s,cat];
	    }

	    if ( (o > 0 || cat == 0 || cat == -1) && nevalprobs[s,cat] > 0 )
	    {
	       if ( length(sname) <= 50 )
		  printf("%-50s %4d %4d %4d %4d %4d %4d", sname, nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat], ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
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
	 if ( cat == 0 )
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
      for (o = 0; o < nsolver; ++o)
      {
	 s = printorder[o];
	 sname = solvername[s];
	 if ( o == 0 )
	 {
	    nodegeomcomp = nodegeom[s,cat];
	    timegeomcomp = timegeom[s,cat];
	    nodeshiftedgeomcomp = nodeshiftedgeom[s,cat];
	    timeshiftedgeomcomp = timeshiftedgeom[s,cat];
	 }

	 if ( (o > 0 || cat == 0 || cat == -1) && nevalprobs[s,cat] > 0 )
	 {
	    if ( length(sname) <= 50 )
	       printf("%-50s %4d %4d %4d %4d %4d %4d", sname, nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat], ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
	    else
	       printf("*%-49s %4d %4d %4d %4d %4d %4d", substr(sname, length(sname)-48), nprocessedprobs[s,cat], nevalprobs[s,cat], nfails[s,cat], ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
	    printf(" %4d %4d", better[s,cat], worse[s,cat]);
	    printf(" %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n",
		   betterobj[s,cat], worseobj[s,cat], feasibles[s,cat],
		   nodegeom[s,cat], nodeshiftedgeom[s,cat], nodegeom[s,cat]/refnodegeom[s,cat],
		   nodeshiftedgeom[s,cat]/refnodeshiftedgeom[s,cat],
		   timegeom[s,cat], timeshiftedgeom[s,cat], timegeom[s,cat]/reftimegeom[s,cat],
		   timeshiftedgeom[s,cat]/reftimeshiftedgeom[s,cat], score[s,cat]);
	 }
      }
      if ( cat == 0 )
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
}
