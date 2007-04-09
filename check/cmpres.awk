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
# $Id: cmpres.awk,v 1.33 2007/04/09 21:12:05 bzfpfend Exp $
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

function isfaster(t,reft,tol)
{
   return (t < 1.0/tol*reft && t <= reft - 0.2);
}

function isslower(t,reft,tol)
{
   return isfaster(reft, t, tol);
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
      iters[nsolver,nprobs[nsolver]] = -nsolver; # different values for each solver -> category "diff"
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
      iters[nsolver,nprobs[nsolver]] = -nsolver; # different values for each solver -> category "diff"
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
      iters[nsolver,nprobs[nsolver]] = $10;
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
      # cat: 0 - all, 1 - different path, 2 - equal path, 3 - all timeout
      for( cat = 0; cat <= 3; cat++ )
      {
         timetotal[s,cat] = 0.0;
         nodetotal[s,cat] = 0.0;
         timegeom[s,cat] = 1.0;
         nodegeom[s,cat] = 1.0;
         timeshiftedgeom[s,cat] = timegeomshift;
         nodeshiftedgeom[s,cat] = nodegeomshift;
         wins[s,cat] = 0;
         nsolved[s,cat] = 0;
         ntimeouts[s,cat] = 0;
         nfails[s,cat] = 0;
         for( wt = 1; wt <= nwintolerances; wt++ )
         {
            better[s,wt,cat] = 0;
            worse[s,wt,cat] = 0;
         }
         betterobj[s,cat] = 0;
         worseobj[s,cat] = 0;
         feasibles[s,cat] = 0;
         score[s,cat] = 1.0;
      }
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
   for( cat = 0; cat <= 3; cat++ )
   {
      nevalprobs[cat] = 0;
      nprocessedprobs[cat] = 0;
   }
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
      itercomp = -1;
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
      iseqpath = 1;
      alltimeout = 1;
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

         # check if all solvers processed the problem
         if( !processed )
	 {
            marker = "?";
	    unprocessed = 1;
	 }

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

      # check if all solvers have same path
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");

         if( processed )
         {
            if( nodecomp == -1 )
            {
               itercomp = iters[s,pidx];
               nodecomp = nodes[s,pidx];
               timecomp = time[s,pidx];
            }
            iseqpath = (iseqpath && iters[s,pidx] == itercomp && nodes[s,pidx] == nodecomp);
            alltimeout = alltimeout && (status[s,pidx] == "timeout");
         }
      }

      # which category?
      if( alltimeout )
         category = 3;
      else if( iseqpath )
         category = 2;
      else
         category = 1;

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
                  nsolved[s,0]++;
                  nsolved[s,category]++;
                  nthissolved++;
               }
            }
            else if( status[s,pidx] == "timeout" )
            {
               marker = ">";
               if( !unprocessed )
               {
                  # if memory limit was exceeded, replace time and nodes by worst time and worst nodes of all runs
                  if( time[s,pidx] < 0.99*worsttime )
                  {
                     nodes[s,pidx] = worstnodes;
                     time[s,pidx] = worsttime;
                  }
                  if( countprob )
                  {
                     ntimeouts[s,0]++;
                     ntimeouts[s,category]++;
                     nthistimeouts++;
                  }
               }
            }
            else
            {
               marker = "!";
               fail = 1;
               if( status[s,pidx] == "readerror" )
                  readerror = 1;
               if( !unprocessed )
               {
                  nfails[s,0]++;
                  nfails[s,category]++;
                  nthisfails++;
               }
            }
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
            line = sprintf("%s %s%10d %s%7.1f", line, feasmark, nodes[s,pidx], marker, time[s,pidx]);
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
		 isfaster(time[s,pidx], timecomp, markbettertime) ||
		 isslower(time[s,pidx], timecomp, markworsetime)) )
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
	       feasibles[s,0]++;
	       feasibles[s,category]++;
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
         nevalprobs[0]++;
         nevalprobs[category]++;
	 reftime = time[printorder[0],probidx[p,printorder[0]]];
	 refobj = primalbound[printorder[0],probidx[p,printorder[0]]];
         for( wt = 1; wt <= nwintolerances; wt++ )
            hasbetter[wt] = 0;
	 hasbetterobj = 0;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            for( cat = 0; cat <= 3; cat = 3*cat + category )
            {
               nep = nevalprobs[cat];
               timetotal[s,cat] += time[s,pidx];
               nodetotal[s,cat] += nodes[s,pidx];
               timegeom[s,cat] = timegeom[s,cat]^((nep-1)/nep) * time[s,pidx]^(1.0/nep);
               nodegeom[s,cat] = nodegeom[s,cat]^((nep-1)/nep) * nodes[s,pidx]^(1.0/nep);
               timeshiftedgeom[s,cat] = timeshiftedgeom[s,cat]^((nep-1)/nep) * (time[s,pidx]+timegeomshift)^(1.0/nep);
               nodeshiftedgeom[s,cat] = nodeshiftedgeom[s,cat]^((nep-1)/nep) * (nodes[s,pidx]+nodegeomshift)^(1.0/nep);
               if( time[s,pidx] <= wintolerances[1]*besttime )
                  wins[s,cat]++;
               for( wt = 1; wt <= nwintolerances; wt++ )
               {
                  if( isfaster(time[s,pidx], reftime, wintolerances[wt]) )
                  {
                     better[s,wt,cat]++;
                     hasbetter[wt] = 1;
                  }
                  else if( isslower(time[s,pidx], reftime, wintolerances[wt]) )
                     worse[s,wt,cat]++;
               }
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
	 besttimegeom = besttimegeom^((nevalprobs[0]-1)/nevalprobs[0]) * besttime^(1.0/nevalprobs[0]);
	 bestnodegeom = bestnodegeom^((nevalprobs[0]-1)/nevalprobs[0]) * bestnodes^(1.0/nevalprobs[0]);
	 besttimeshiftedgeom = besttimeshiftedgeom^((nevalprobs[0]-1)/nevalprobs[0]) * (besttime+timegeomshift)^(1.0/nevalprobs[0]);
	 bestnodeshiftedgeom = bestnodeshiftedgeom^((nevalprobs[0]-1)/nevalprobs[0]) * (bestnodes+nodegeomshift)^(1.0/nevalprobs[0]);
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
   probnumstr = "("nevalprobs[0]")";
   printf("%-14s %5s", "total", probnumstr);
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
         printf(" %11d %8d", nodetotal[s,0], timetotal[s,0]);
      else
         printf(" %11d %8d              ", nodetotal[s,0], timetotal[s,0]);
   }
   printf("\n");
   printf("%-20s", "geom. mean");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
      {
         printf(" %11d %8.1f", nodegeom[s,0], timegeom[s,0]);
         nodegeomcomp = nodegeom[s,0];
         timegeomcomp = timegeom[s,0];
         nodetotalcomp = nodetotal[s,0];
         timetotalcomp = timetotal[s,0];
      }
      else
         printf(" %11d %8.1f %6.2f %6.2f", nodegeom[s,0], timegeom[s,0], nodegeom[s,0]/nodegeomcomp, timegeom[s,0]/timegeomcomp);
   }
   printf("\n");
   printf("%-20s", "shifted geom.");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      for( cat = 0; cat <= 3; cat++ )
      {
         nodeshiftedgeom[s,cat] -= nodegeomshift;
         timeshiftedgeom[s,cat] -= timegeomshift;
         nodeshiftedgeom[s,cat] = max(nodeshiftedgeom[s,cat], 1.0);
         timeshiftedgeom[s,cat] = max(timeshiftedgeom[s,cat], 1.0);
      }
      if( o == 0 )
      {
         printf(" %11d %8.1f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0]);
         nodeshiftedgeomcomp = nodeshiftedgeom[s,0];
         timeshiftedgeomcomp = timeshiftedgeom[s,0];
      }
      else
         printf(" %11d %8.1f %6.2f %6.2f", nodeshiftedgeom[s,0], timeshiftedgeom[s,0],
            nodeshiftedgeom[s,0]/nodeshiftedgeomcomp, timeshiftedgeom[s,0]/timeshiftedgeomcomp);
   }
   bestnodeshiftedgeom -= nodegeomshift;
   besttimeshiftedgeom -= timegeomshift;
   bestnodeshiftedgeom = max(bestnodeshiftedgeom, 1.0);
   besttimeshiftedgeom = max(besttimeshiftedgeom, 1.0);
   
   printf("\n");
   printhline(nsolver);

   for( cat = 0; cat <= 3; cat++ )
   {
      if( nprocessedprobs[cat] == 0 )
         continue;

      header = (cat == 0 ? "all" : (cat == 1 ? "diff" : (cat == 2 ? "equal" : "timeout")));
      printf("\n");
      printf("%55s", "");
      for( wt = 1; wt <= nwintolerances; wt++ )
         printf("   %3d%%   ", 100.0*(wintolerances[wt]));
      printf("    obj   \n");
      printf("%-7s (%4d proc, %4d eval)      fail time solv wins", header, nprocessedprobs, nevalprobs);
      for( wt = 1; wt <= nwintolerances; wt++ )
         printf(" bett wors");
      printf(" bett wors feas     nodes   shnodes    nodesQ  shnodesQ    time  shtime   timeQ shtimeQ   score\n");
          
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         if( o == 0 )
         {
            nodegeomcomp = nodegeom[s,cat];
            timegeomcomp = timegeom[s,cat];
            nodeshiftedgeomcomp = nodeshiftedgeom[s,cat];
            timeshiftedgeomcomp = timeshiftedgeom[s,cat];
         }
         printf("%-35s %4d %4d %4d %4d", solvername[s], nfails[s,cat], ntimeouts[s,cat], nsolved[s,cat], wins[s,cat]);
         for( wt = 1; wt <= nwintolerances; wt++ )
            printf(" %4d %4d", better[s,wt,cat], worse[s,wt,cat]);
         printf(" %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", 
                betterobj[s,cat], worseobj[s,cat], feasibles[s,cat],
                nodegeom[s,cat], nodeshiftedgeom[s,cat], nodegeom[s,cat]/nodegeomcomp, nodeshiftedgeom[s,cat]/nodeshiftedgeomcomp,
                timegeom[s,cat], timeshiftedgeom[s,cat], timegeom[s,cat]/timegeomcomp, timeshiftedgeom[s,cat]/timeshiftedgeomcomp,
                score[s,cat]);
      }
      if( cat == 0 )
      {
         printf("%-35s %4d %4d %4d %4s", "optimal auto settings", bestnfails, bestntimeouts, bestnsolved, "");
         for( wt = 1; wt <= nwintolerances; wt++ )   
            printf(" %4d %4s", bestbetter[wt], "");
         printf(" %4d %4s %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7s\n",
                bestbetterobj, "", bestfeasibles,
                bestnodegeom, bestnodeshiftedgeom, bestnodegeom/nodegeomcomp, bestnodeshiftedgeom/nodeshiftedgeomcomp,
                besttimegeom, besttimeshiftedgeom, besttimegeom/timegeomcomp, besttimeshiftedgeom/timeshiftedgeomcomp,
                "");
      }
   }

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
