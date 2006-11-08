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
# $Id: cmpres.awk,v 1.17 2006/11/08 23:22:43 bzfpfend Exp $
#
#@file    compare.awk
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
         printf("--------------------+----------+-------+");
      else
         printf("----------+-------+------+------+");
   }
   printf("-------------\n");
}
BEGIN {
   timegeomshift = 60.0;
   nodegeomshift = 1000.0;
   mintime = 0.1;
   wintolerance = 1.1;
   markbettertime = 1.1;
   markworsetime  = 1.1;
   markbetternode = 5.0;
   markworsenode  = 5.0;
   onlymarked = 0;
   maxscore = 10.0;

   problistlen = 0;
   nsolver = 0;
   nprobs[nsolver] = 0;
}
/^@01 / {
   solvername[nsolver] = $2;
   nsolver++;
   nprobs[nsolver] = 0;
}
// {
   # check if this is a useable line
   if( $12 == "ok" || $12 == "timeout" || $12 == "unknown" || $12 == "abort" || $12 == "fail" || $12 == "readerror" )
   {
      # collect data
      name[nsolver,nprobs[nsolver]] = $1;
      type[nsolver,nprobs[nsolver]] = $2;
      conss[nsolver,nprobs[nsolver]] = $5;
      vars[nsolver,nprobs[nsolver]] = $6;
      dualbound[nsolver,nprobs[nsolver]] = $7;
      primalbound[nsolver,nprobs[nsolver]] = $8;
      gap[nsolver,nprobs[nsolver]] = $9;
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
}
END {
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
      better[s] = 0;
      worse[s] = 0;
      score[s] = 1.0;
   }

   # calculate the order in which the columns should be printed: CPLEX < SCIP, default < non-default
   for( s = 0; s < nsolver; ++s )
   {
      for( i = 0; i < s; ++i )
      {
         if( substr(solvername[s], 1, 5) == "CPLEX" && substr(solvername[i], 1, 5) != "CPLEX" )
            break;
         if( substr(solvername[s], 1, 5) == substr(solvername[i], 1, 5) &&
            match(solvername[s], "default") != 0 && match(solvername[i], "default") == 0 )
            break;
      }
      for( j = s-1; j >= i; --j )
         printorder[j+1] = printorder[j];
      printorder[i] = s;
   }

   # print headers
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
         printf(" %37s |", solvername[s]);
      else
         printf(" %30s |", solvername[s]);
   }
   printf("\n");
   printhline(nsolver);
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf("  Name              |    Nodes |  Time |");
      else
         printf("    Nodes |  Time | NodQ | TimQ |");
   }
   printf(" bounds check\n");
   printhline(nsolver);
   
   # display the problem results and calculate mean values
   nevalprobs = 0;
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
      ismini = 0;
      ismaxi = 0;
      mark = " ";
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");
         if( processed && name[s,pidx] != p )
            printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

	 # make sure, nodes and time are non-zero for geometric means
	 nodes[s,pidx] = max(nodes[s,pidx], 1);
	 time[s,pidx] = max(time[s,pidx], mintime);

         # check if solver ran successfully (i.e., no abort nor fail)
         if( !processed )
	 {
            marker = "?";
	    unprocessed = 1;
	 }
         else if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
         {
            nsolved[s]++;
            marker = " ";
            besttime = min(besttime, time[s,pidx]);
         }
         else if( status[s,pidx] == "timeout" )
         {
            ntimeouts[s]++;
            marker = ">";
         }
         else
         {
            if( status[s,pidx] == "readerror" )
               readerror = 1;
            nfails[s]++;
            fail = 1;
            marker = "!";
         }

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
            line = sprintf("%s          -       -", line);
         else
         {
            line = sprintf("%s %10d %s%6.1f", line, nodes[s,pidx], marker, time[s,pidx]);
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
		(nodes[s,pidx] > markworsenode * nodecomp ||
		 nodes[s,pidx] < 1.0/markbetternode * nodecomp ||
		 time[s,pidx] > markworsetime * timecomp ||
		 time[s,pidx] < 1.0/markbettertime * timecomp ) )
	       mark = "*";
         }
      }

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
      else if( (ismini && ismaxi) ||
	       (ismini && maxdb - minpb > 1e-5 * max(max(abs(maxdb), abs(minpb)), 1.0)) ||
	       (ismaxi && maxpb - mindb > 1e-5 * max(max(abs(maxpb), abs(mindb)), 1.0)) ||
	       (!ismini && !ismaxi && abs(maxpb - minpb) > 1e-5 * max(abs(maxpb), 1.0)) )
      {
         line = sprintf("%s  inconsistent", line);
         fail = 1;
	 mark = " ";
      }
      else if( unprocessed )
      {
         line = sprintf("%s  unprocessed", line);
	 mark = " ";
      }
      else
         line = sprintf("%s  ok", line);
      if( !onlymarked || mark == "*" )
	printf("%s %s\n", mark, line);

      # calculate totals and means for instances where no solver failed
      if( !fail && !unprocessed )
      {
         nevalprobs++;
	 reftime = time[printorder[0],probidx[p,printorder[0]]];
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            timetotal[s] += time[s,pidx];
            nodetotal[s] += nodes[s,pidx];
            timegeom[s] = timegeom[s]^((nevalprobs-1)/nevalprobs) * time[s,pidx]^(1.0/nevalprobs);
            nodegeom[s] = nodegeom[s]^((nevalprobs-1)/nevalprobs) * nodes[s,pidx]^(1.0/nevalprobs);
            timeshiftedgeom[s] = timeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (time[s,pidx]+timegeomshift)^(1.0/nevalprobs);
            nodeshiftedgeom[s] = nodeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (nodes[s,pidx]+nodegeomshift)^(1.0/nevalprobs);
            if( time[s,pidx] <= wintolerance*besttime )
               wins[s]++;
            if( time[s,pidx] < 1.0/wintolerance*reftime )
               better[s]++;
            else if( time[s,pidx] > wintolerance*reftime )
               worse[s]++;
	    thisscore = reftime/time[s,pidx];
	    thisscore = max(thisscore, 1/maxscore);
	    thisscore = min(thisscore, maxscore);
	    score[s] = score[s]^((nevalprobs-1)/nevalprobs) * thisscore^(1.0/nevalprobs);
         }
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
         printf(" %10d %7d", nodetotal[s], timetotal[s]);
      else
         printf(" %10d %7d              ", nodetotal[s], timetotal[s]);
   }
   printf("\n");
   printf("%-20s", "geom. mean");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
      {
         printf(" %10d %7.1f", nodegeom[s], timegeom[s]);
         nodegeomcomp = nodegeom[s];
         timegeomcomp = timegeom[s];
      }
      else
         printf(" %10d %7.1f %6.2f %6.2f", nodegeom[s], timegeom[s], nodegeom[s]/nodegeomcomp, timegeom[s]/timegeomcomp);
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
         printf(" %10d %7.1f", nodeshiftedgeom[s], timeshiftedgeom[s]);
         nodeshiftedgeomcomp = nodeshiftedgeom[s];
         timeshiftedgeomcomp = timeshiftedgeom[s];
      }
      else
         printf(" %10d %7.1f %6.2f %6.2f", nodeshiftedgeom[s], timeshiftedgeom[s],
            nodeshiftedgeom[s]/nodeshiftedgeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp);
   }
   printf("\n");
   printhline(nsolver);

   printf("\n");
   printf("solver                            fail time solv wins bett wors     nodes   shnodes    nodesQ  shnodesQ    time  shtime   timeQ shtimeQ   score\n");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      printf("%-33s %4d %4d %4d %4d %4d %4d %9d %9d %9.2f %9.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", 
         solvername[s], nfails[s], ntimeouts[s], nsolved[s], wins[s], better[s], worse[s],
         nodegeom[s], nodeshiftedgeom[s], nodegeom[s]/nodegeomcomp, nodeshiftedgeom[s]/nodeshiftedgeomcomp,
         timegeom[s], timeshiftedgeom[s], timegeom[s]/timegeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp,
	 score[s]);
   }
}
