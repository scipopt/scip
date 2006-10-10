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
# $Id: cmpres.awk,v 1.8 2006/10/10 14:09:08 bzfpfend Exp $
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
         printf("------------------+----------+-------+");
      else
         printf("----------+-------+------+------+");
   }
   printf("-------------\n");
}
BEGIN {
   timegeomshift = 60.0;
   nodegeomshift = 1000.0;
   wintolerance = 0.05;

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
      nodes[nsolver,nprobs[nsolver]] = max($8,1);
      time[nsolver,nprobs[nsolver]] = ceil(max($9,1.0));
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
         printf(" %35s |", solvername[s]);
      else
         printf(" %30s |", solvername[s]);
   }
   printf("\n");
   printhline(nsolver);
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf("Name              |    Nodes |  Time |");
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
      printf("%-18s", p);
      fail = 0;
      readerror = 0;
      maxdb = -1e+100;
      minpb = +1e+100;
      nodecomp = -1;
      timecomp = -1;
      besttime = +1e+100;
      for( o = 0; o < nsolver; ++o )
      {
         s = printorder[o];
         pidx = probidx[p,s];
         processed = (pidx != "");
         if( processed && name[s,pidx] != p )
            printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

         # check if solver ran successfully (i.e., no abort nor fail)
         if( !processed )
            marker = "?";
         else if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
         {
            nsolved[s]++;
            marker = " ";
            maxdb = max(maxdb, dualbound[s,pidx]);
            minpb = min(minpb, primalbound[s,pidx]);
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

         # print statistics
         if( !processed )
            printf("          -       -");
         else
         {
            printf(" %10d %s%6.1f", nodes[s,pidx], marker, time[s,pidx]);
            if( nodecomp == -1 )
            {
               nodecomp = nodes[s,pidx];
               timecomp = time[s,pidx];
            }
         }
         if( o > 0 )
         {
            if( !processed )
               printf("      -");
            else if( nodes[s,pidx]/nodecomp > 999.99 )
               printf("  Large");
            else
               printf(" %6.2f", nodes[s,pidx]/nodecomp);
            if( !processed )
               printf("      -");
            else if( time[s,pidx]/timecomp > 999.99 )
               printf("  Large");
            else
               printf(" %6.2f", time[s,pidx]/timecomp);
         }
      }

      # check for inconsistency in the primal and dual bounds
      if( readerror )
         printf("  readerror");
      else if( fail )
         printf("  fail");
      else if( maxdb - minpb > 1e-5 * max(max(abs(maxdb), abs(minpb)), 1.0) )
      {
         printf("  inconsistent");
         fail = 1;
      }
      else if( probcnt[p] != nsolver )
         printf("  unprocessed");
      else
         printf("  ok");
      printf("\n");

      # calculate totals and means for instances where no solver failed
      if( !fail && probcnt[p] == nsolver )
      {
         nevalprobs++;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            timetotal[s] += time[s,pidx];
            nodetotal[s] += nodes[s,pidx];
            timegeom[s] = timegeom[s]^((nevalprobs-1)/nevalprobs) * time[s,pidx]^(1.0/nevalprobs);
            nodegeom[s] = nodegeom[s]^((nevalprobs-1)/nevalprobs) * nodes[s,pidx]^(1.0/nevalprobs);
            timeshiftedgeom[s] = timeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (time[s,pidx]+timegeomshift)^(1.0/nevalprobs);
            nodeshiftedgeom[s] = nodeshiftedgeom[s]^((nevalprobs-1)/nevalprobs) * (nodes[s,pidx]+nodegeomshift)^(1.0/nevalprobs);
            if( time[s,pidx] <= (1.0+wintolerance)*besttime )
               wins[s]++;
            if( time[s,pidx] < (1.0-wintolerance)*time[0,pidx] )
               better[s]++;
            else if( time[s,pidx] > (1.0+wintolerance)*time[0,pidx] )
               worse[s]++;
         }
      }
   }
   printhline(nsolver);

   # print solvers' overall statistics
   probnumstr = "("nevalprobs")";
   printf("%-12s %5s", "total", probnumstr);
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      if( o == 0 )
         printf(" %10d %7d", nodetotal[s], timetotal[s]);
      else
         printf(" %10d %7d              ", nodetotal[s], timetotal[s]);
   }
   printf("\n");
   printf("%-18s", "geom. mean");
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
   printf("%-18s", "shifted geom.");
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
   printf("solver                           fails timeout  solved    wins  better   worse    time  shtime   timeQ shtimeQ\n");
   for( o = 0; o < nsolver; ++o )
   {
      s = printorder[o];
      printf("%-30s %7d %7d %7d %7d %7d %7d %7.1f %7.1f %7.2f %7.2f\n", 
         solvername[s], nfails[s], ntimeouts[s], nsolved[s], wins[s], better[s], worse[s], timegeom[s], timeshiftedgeom[s],
         timegeom[s]/timegeomcomp, timeshiftedgeom[s]/timeshiftedgeomcomp);
   }
}
