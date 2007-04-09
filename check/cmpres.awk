#!/bin/gawk -f
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
# $Id: cmpres.awk,v 1.32 2007/04/09 20:41:53 bzfpfend Exp $
#
#@file    cmpres.awk
#@brief   SCIP Check Comparison Report Generator
#@author  Tobias Achterberg
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
   return (10*x) == int(10*x) ? (x) : int(10*x+1.0)/10.0;
}
function printhline(nsolver)
{
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf("--------------------+-+---------+--------+");
      else
         printf("-+---------+--------+------+------+");
   }
   printf("-------------\n");
}
BEGIN {
   timegeomshift = 60.0;
   nodegeomshift = 1000.0;
   mintime = 0.1;
   nwintolerances = 2;
   wintolerances[1] = 1.1;
   wintolerances[2] = 2.0;
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

   problistlen = 0;
   nsolver = 0;
   nprobs[nsolver] = 0;
   fulltotaltime = 0.0;
}
/^@01 / {
   solvername[nsolver] = $2;
   nsolver++;
   nprobs[nsolver] = 0;
}
// {
   # check if this is a useable line
   if( $10 == "ok" || $10 == "timeout" || $10 == "unknown" || $10 == "abort" || $10 == "fail" || $10 == "readerror" )
   {
      # collect data
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $3;
      vars[nsolver,nprobs[nsolver]] = $4;
      dualbound[nsolver,nprobs[nsolver]] = $5;
      primalbound[nsolver,nprobs[nsolver]] = $6;
      gap[nsolver,nprobs[nsolver]] = $7;
      iterations[nsolver,nprobs[nsolver]] = -nsolver; # different values for each solver -> category "diff"
      nodes[nsolver,nprobs[nsolver]] = max($8,1);
      time[nsolver,nprobs[nsolver]] = ceil(max($9,mintime));
      status[nsolver,nprobs[nsolver]] = $10;
      probidx[$1,nsolver] = nprobs[nsolver];
      probcnt[$1]++;
      nprobs[nsolver]++;
      if( probcnt[$1] == 1 )
      {
         problist[problistlen] = $1;
         problistlen++;
      }
   }
   else if( $12 == "ok" || $12 == "timeout" || $12 == "unknown" || $12 == "abort" || $12 == "fail" || $12 == "readerror" )
   {
      # collect data (line with original and presolved problem size)
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = $7;
      primalbound[nsolver,nprobs[nsolver]] = $8;
      gap[nsolver,nprobs[nsolver]] = $9;
      iterations[nsolver,nprobs[nsolver]] = -nsolver; # different values for each solver -> category "diff"
      nodes[nsolver,nprobs[nsolver]] = max($10,1);
      time[nsolver,nprobs[nsolver]] = ceil(max($11,mintime));
      status[nsolver,nprobs[nsolver]] = $12;
      probidx[$1,nsolver] = nprobs[nsolver];
      probcnt[$1]++;
      nprobs[nsolver]++;
      if( probcnt[$1] == 1 )
      {
         problist[problistlen] = $1;
         problistlen++;
      }
   }
   else if( $13 == "ok" || $13 == "timeout" || $13 == "unknown" || $13 == "abort" || $13 == "fail" || $13 == "readerror" )
   {
      # collect data (line with original and presolved problem size and simplex iterations)
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = $7;
      primalbound[nsolver,nprobs[nsolver]] = $8;
      gap[nsolver,nprobs[nsolver]] = $9;
      iterations[nsolver,nprobs[nsolver]] = $10;
      nodes[nsolver,nprobs[nsolver]] = max($11,1);
      time[nsolver,nprobs[nsolver]] = ceil(max($12,mintime));
      status[nsolver,nprobs[nsolver]] = $13;
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
   if( nsolver == 0 )
   {
      printf("no instances found in log file\n");
      exit 1;
   }

   # process exclude string
   n = split(exclude, a, ",");
   for( i = 1; i <= n; i++ )
      excluded[a[i]] = 1;

   # initialize means
   for( s = 0; s < nsolver; ++s )
   {
      timetotal[s] = 0.0;
      nodetotal[s] = 0.0;
      timegeom[s] = 1.0;
      nodegeom[s] = 1.0;
      timeshiftedgeom[s] = timegeomshift;
      nodeshiftedgeom[s] = nodegeomshift;
      wins[s] = 0;
      for( wt = 1; wt <= nwintolerances; wt++ )
      {
         better[s,wt] = 0;
         worse[s,wt] = 0;
      }
      betterobj[s] = 0;
      worseobj[s] = 0;
      feasibles[s] = 0;
      score[s] = 1.0;
   }
   besttimegeom = 1.0;
   bestnodegeom = 1.0;
   besttimeshiftedgeom = timegeomshift;
   bestnodeshiftedgeom = nodegeomshift;
   bestnsolved = 0;
   bestntimeouts = 0;
   bestnfails = 0;
   for( wt = 1; wt <= nwintolerances; wt++ )
      bestbetter[wt] = 0;
   bestbetterobj = 0;
   bestfeasibles = 0;

   # calculate the order in which the columns should be printed: CPLEX < SCIP, default < non-default
   for( s = 0; s < nsolver; ++s )
   {
      for( o = 0; o < s; ++o )
      {
         i = printorder[o];
         if( substr(solvername[s], 1, 5) == "CPLEX" && substr(solvername[i], 1, 5) != "CPLEX" )
            break;
         if( substr(solvername[s], 1, 5) == substr(solvername[i], 1, 5) &&
            match(solvername[s], "default") != 0 && match(solvername[i], "default") == 0 )
            break;
         if( substr(solvername[s], 1, 5) == substr(solvername[i], 1, 5) &&
            (match(solvername[s], "default") == 0) == (match(solvername[i], "default") == 0) &&
            solvername[s] < solvername[i] )
            break;
      }
      for( j = s-1; j >= o; --j )
         printorder[j+1] = printorder[j];
      printorder[o] = s;
   }

   # print headers
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
         printf(" %39s |", solvername[s]);
      else
         printf(" %32s |", solvername[s]);
   }
   printf("\n");
   printhline(nsolver);
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf("  Name              |F|   Nodes |   Time |");
      else
         printf("F|   Nodes |   Time | NodQ | TimQ |");
   }
   printf(" bounds check\n");
   printhline(nsolver);
   
   # display the problem results and calculate mean values
   nevalprobs = 0;
   nprocessedprobs = 0;
   for( i = 0; i < problistlen; ++i )
   {
      p = problist[i];
      line = sprintf("%-18s", p);
      fail = 0;
      readerror = 0;
      unprocessed = 0;
      mindb = +1e+100;
      maxdb = -1e+100;
      minpb = +1e+100;
      maxpb = -1e+100;
      nodecomp = -1;
      timecomp = -1;
      besttime = +1e+100;
      bestnodes = +1e+100;
      worsttime = -1e+100;
      worstnodes = -1e+100;
      nthisunprocessed = 0;
      nthissolved = 0;
      nthistimeouts = 0;
      nthisfails = 0;
      ismini = 0;
      ismaxi = 0;
      mark = " ";
      countprob = 1;

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

         # check if solver ran successfully (i.e., no abort nor fail)
         if( processed && (status[s,pidx] == "ok" || status[s,pidx] == "unknown" || status[s,pidx] == "timeout") )
         {
            besttime = min(besttime, time[s,pidx]);
            bestnodes = min(bestnodes, nodes[s,pidx]);
            worsttime = max(worsttime, time[s,pidx]);
	    worstnodes = max(worstnodes, nodes[s,pidx]);
         }
         else
            countprob = 0;
      }
      worsttime = max(worsttime, mintime);
      worstnodes = max(worstnodes, 1);

      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");
         if( processed && name[s,pidx] != p )
            printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

         # check if solver ran successfully (i.e., no abort nor fail)
         if( !processed )
	 {
            marker = "?";
	    unprocessed = 1;
            nthisunprocessed++;
	 }
         else if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
         {
            marker = " ";
            nsolved[s]++;
            nthissolved++;
         }
         else if( status[s,pidx] == "timeout" )
         {
            # if memory limit was exceeded, replace time and nodes by worst time and worst nodes of all runs
            if( time[s,pidx] < 0.99*worsttime )
            {
               nodes[s,pidx] = worstnodes;
               time[s,pidx] = worsttime;
            }
            marker = ">";
            if( countprob )
            {
               ntimeouts[s]++;
               nthistimeouts++;
            }
         }
         else
         {
            marker = "!";
            if( status[s,pidx] == "readerror" )
               readerror = 1;
            fail = 1;
            nfails[s]++;
            nthisfails++;
         }

	 if( primalbound[s,pidx] < 1e+20 )
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
            line = sprintf("%s %s%10d %s%7.1f", line, feasmark, nodes[s,pidx], marker, time[s,pidx]);
            if( nodecomp == -1 )
            {
               nodecomp = nodes[s,pidx];
               timecomp = time[s,pidx];
            }
         }
         if( o > 0 )
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
		(nodes[s,pidx] > markworsenodes * nodecomp ||
		 nodes[s,pidx] < 1.0/markbetternodes * nodecomp ||
		 time[s,pidx] > markworsetime * timecomp ||
		 time[s,pidx] < 1.0/markbettertime * timecomp ) )
	       mark = "*";
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
         nprocessedprobs++;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
	    if( primalbound[s,pidx] < 1e+20 ) {
	       feasibles[s]++;
	       hasfeasible = 1;
	    }
	 }
	 if( hasfeasible )
	   bestfeasibles++;
      }

      if( (!onlymarked || mark == "*") && (!onlyprocessed || !unprocessed) &&
          (!onlyfeasible || hasfeasible) && (!onlyinfeasible || !hasfeasible) &&
          (!onlyfail || fail) )
         printf("%s %s\n", mark, line);

      # calculate totals and means for instances where no solver failed
      if( !fail && !unprocessed &&
          (!onlyfeasible || hasfeasible) && (!onlyinfeasible || !hasfeasible) )
      {
         nevalprobs++;
	 reftime = time[printorder[0],probidx[p,printorder[0]]];
	 refobj = primalbound[printorder[0],probidx[p,printorder[0]]];
         for( wt = 1; wt <= nwintolerances; wt++ )
            hasbetter[wt] = 0;
	 hasbetterobj = 0;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            timetotal[s] += time[s,pidx];
            nodetotal[s] += nodes[s,pidx];
            timegeom[s] = timegeom[s]^((nevalprobs-1)/nevalprobs) * time[s,pidx]^(1.0/nevalprobs);
            nodegeom[s] = nodegeom[s]^((nevalprobs-1)/nevalprobs) * nodes[s,pidx]^(1.0/nevalprobs);
            timeshiftedgeom[s] = timeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (time[s,pidx]+timegeomshift)^(1.0/nevalprobs);
            nodeshiftedgeom[s] = nodeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (nodes[s,pidx]+nodegeomshift)^(1.0/nevalprobs);
            if( time[s,pidx] <= wintolerances[1]*besttime )
               wins[s]++;
            for( wt = 1; wt <= nwintolerances; wt++ )
            {
               if( time[s,pidx] < 1.0/wintolerances[wt]*reftime ) {
                  better[s,wt]++;
                  hasbetter[wt] = 1;
               }
               else if( time[s,pidx] > wintolerances[wt]*reftime )
                  worse[s,wt]++;
            }
	    pb = primalbound[s,pidx];
	    if( (ismini && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
		(ismaxi && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ) {
 	       betterobj[s]++;
	       hasbetterobj = 1;
	    }
	    else if( (ismini && pb - refobj > +0.01 * max(max(abs(refobj), abs(pb)), 1.0)) ||
		     (ismaxi && pb - refobj < -0.01 * max(max(abs(refobj), abs(pb)), 1.0)) )
 	       worseobj[s]++;
	    thisscore = reftime/time[s,pidx];
	    thisscore = max(thisscore, 1/maxscore);
	    thisscore = min(thisscore, maxscore);
	    score[s] = score[s]^((nevalprobs-1)/nevalprobs) * thisscore^(1.0/nevalprobs);
         }
	 besttimegeom = besttimegeom^((nevalprobs-1)/nevalprobs) * besttime^(1.0/nevalprobs);
	 bestnodegeom = bestnodegeom^((nevalprobs-1)/nevalprobs) * bestnodes^(1.0/nevalprobs);
	 besttimeshiftedgeom = besttimeshiftedgeom^((nevalprobs-1)/nevalprobs) * (besttime+timegeomshift)^(1.0/nevalprobs);
	 bestnodeshiftedgeom = bestnodeshiftedgeom^((nevalprobs-1)/nevalprobs) * (bestnodes+nodegeomshift)^(1.0/nevalprobs);
         for( wt = 1; wt <= nwintolerances; wt++ )
         {
            if( hasbetter[wt] )
               bestbetter[wt]++;
         }
	 if( hasbetterobj )
	   bestbetterobj++;
      }
   }
   printhline(nsolver);

   # print solvers' overall statistics
   probnumstr = "("nevalprobs")";
   printf("%-14s %5s", "total", probnumstr);
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
         printf(" %11d %8d", nodetotal[s], timetotal[s]);
      else
         printf(" %11d %8d              ", nodetotal[s], timetotal[s]);
   }
   printf("\n");
   printf("%-20s", "geom. mean");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
      {
         printf(" %11d %8.1f", nodegeom[s], timegeom[s]);
         nodegeomcomp = nodegeom[s];
         timegeomcomp = timegeom[s];
         nodetotalcomp = nodetotal[s];
         timetotalcomp = timetotal[s];
      }
      else
         printf(" %11d %8.1f %6.2f %6.2f", nodegeom[s], timegeom[s], nodegeom[s]/nodegeomcomp, timegeom[s]/timegeomcomp);
   }
   printf("\n");
   printf("%-20s", "shifted geom.");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      nodeshiftedgeom[s] -= nodegeomshift;
      timeshiftedgeom[s] -= timegeomshift;
      nodeshiftedgeom[s] = max(nodeshiftedgeom[s], 1.0);
      timeshiftedgeom[s] = max(timeshiftedgeom[s], 1.0);
      if( o == 0 )
      {
         printf(" %11d %8.1f", nodeshiftedgeom[s], timeshiftedgeom[s]);
         nodeshiftedgeomcomp = nodeshiftedgeom[s];
         timeshiftedgeomcomp = timeshiftedgeom[s];
      }
      else
         printf(" %11d %8.1f %6.2f %6.2f", nodeshiftedgeom[s], timeshiftedgeom[s],
            nodeshiftedgeom[s]/nodeshiftedgeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp);
   }
   bestnodeshiftedgeom -= nodegeomshift;
   besttimeshiftedgeom -= timegeomshift;
   bestnodeshiftedgeom = max(bestnodeshiftedgeom, 1.0);
   besttimeshiftedgeom = max(besttimeshiftedgeom, 1.0);
   
   printf("\n");
   printhline(nsolver);

   printf("\n");
   printf("%55s", "");
   for( wt = 1; wt <= nwintolerances; wt++ )
      printf("   %3d%%   ", 100.0*(wintolerances[wt]));
   printf("    obj   \n");
   printf("solver (%4d proc, %4d eval)       fail time solv wins", nprocessedprobs, nevalprobs);
   for( wt = 1; wt <= nwintolerances; wt++ )
      printf(" bett wors");
   printf(" bett wors feas     nodes   shnodes    nodesQ  shnodesQ    time  shtime   timeQ shtimeQ   score\n");
          
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      printf("%-35s %4d %4d %4d %4d", solvername[s], nfails[s], ntimeouts[s], nsolved[s], wins[s]);
      for( wt = 1; wt <= nwintolerances; wt++ )
         printf(" %4d %4d", better[s,wt], worse[s,wt]);
      printf(" %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", 
         betterobj[s], worseobj[s], feasibles[s],
         nodegeom[s], nodeshiftedgeom[s], nodegeom[s]/nodegeomcomp, nodeshiftedgeom[s]/nodeshiftedgeomcomp,
         timegeom[s], timeshiftedgeom[s], timegeom[s]/timegeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp,
	 score[s]);
   }
   printf("%-35s %4d %4d %4d %4s", "optimal auto settings", bestnfails, bestntimeouts, bestnsolved, "");
   for( wt = 1; wt <= nwintolerances; wt++ )   
      printf(" %4d %4s", bestbetter[wt], "");
   printf(" %4d %4s %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7s\n",
      bestbetterobj, "", bestfeasibles,
      bestnodegeom, bestnodeshiftedgeom, bestnodegeom/nodegeomcomp, bestnodeshiftedgeom/nodeshiftedgeomcomp,
      besttimegeom, besttimeshiftedgeom, besttimegeom/timegeomcomp, besttimeshiftedgeom/timeshiftedgeomcomp,
      "");

   printf("total time over all settings: %.1f sec = %.1f hours = %.1f days = %.1f weeks = %.1f months\n",
      fulltotaltime, fulltotaltime/3600.0, fulltotaltime/(3600.0*24), fulltotaltime/(3600.0*24*7),
      fulltotaltime/(3600.0*24*30));

   # generate tex summary file
   if( texfile != "" )
   {
      printf("generating tex file <%s>\n", texfile);
      printf("{\\sffamily\n") > texfile;
      printf("\\scriptsize\n") > texfile;
      printf("\\setlength{\\extrarowheight}{1pt}\n") > texfile;
      printf("\\setlength{\\tabcolsep}{2pt}\n") > texfile;
      printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n") > texfile;
      printf("\\newcommand{\\spc}{\\;\\;\\;\\;\\;\\;}\n") > texfile;
      printf("\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}l@{\\spc}rrr@{\\spc}rrrrrr@{\\spc}rrr@{\\spc}rrr@{}}\n") > texfile;
      printf("\\toprule\n") > texfile;
      printf("& & \\multicolumn{2}{c@{\\spc}}{fast/slow} & \\multicolumn{6}{c@{\\spc}}{all instances} & \\multicolumn{3}{c@{\\spc}}{different path} & \\multicolumn{3}{c}{equal path} \\\\\n") > texfile;
      printf("setting & >T") > texfile;
      for( wt = 1; wt <= nwintolerances; wt++ )
         printf(" & %d\\%", 100.0*wintolerances[wt]) > texfile;
      printf(" & $\\textit{n}_\\textit{gm}$ & $\\textit{n}_\\textit{sgm}$ & $\\textit{n}_\\textit{tot}$ & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$ & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$ & $\\textit{t}_\\textit{gm}$ & $\\textit{t}_\\textit{sgm}$ & $\\textit{t}_\\textit{tot}$ \\\\\n") > texfile;
      printf("\\midrule\n") > texfile;

      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         sname = solvername[s];
         sub(/.*:/, "", sname);
         sub(/.*_/, "", sname);
         printf("%-35s & %4d", sname, ntimeouts[s]) > texfile;
         for( wt = 1; wt <= nwintolerances; wt++ )
         {
            bws = sprintf("%d/%d", better[s,wt], worse[s,wt]);
            printf(" & %5s", bws) > texfile;
         }
         printf(" & %9.2f & %9.2f & %9.2f & %7.2f & %7.2f & %7.2f \\\\\n", 
            nodegeom[s]/nodegeomcomp, nodeshiftedgeom[s]/nodeshiftedgeomcomp, nodetotal[s]/nodetotalcomp,
            timegeom[s]/timegeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp, timetotal[s]/timetotalcomp,
            score[s]) > texfile;
      }

      printf("\\bottomrule\n") > texfile;
      printf("\\end{tabular*}\n") > texfile;
      printf("}\n") > texfile;
   }
}
