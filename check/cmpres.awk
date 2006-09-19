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
# $Id: cmpres.awk,v 1.2 2006/09/19 00:48:52 bzfpfend Exp $
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
   return (x) == int(x) ? (x) : int(x+1.0);
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
   
   nsolver = 0;
   nprobs[nsolver] = 0;
}
/^@01 / {
   solvername[nsolver] = $2;
   nsolver++;
}
// {
   # check if this is a useable line
   if( $10 == "ok" || $10 == "timeout" || $10 == "unknown" || $10 == "abort" || $10 == "fail" )
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
      if( nsolver == 0 )
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
   }

   # print headers
   for( s = 0; s < nsolver; ++s )
   {
       if( s == 0 )
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
      if( probcnt[p] == nsolver )
      {
         printf("%-18s", p);
         fail = 0;
         maxdb = -1e+100;
         minpb = +1e+100;
         for( s = 0; s < nsolver; ++s )
         {
            pidx = probidx[p,s];
            if( name[s,pidx] != p )
               printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

            # check if solver ran successfully (i.e., no abort nor fail)
            if( status[s,pidx] == "ok" || status[s,pidx] == "unknown" )
            {
               nsolved[s]++;
               marker = " ";
               maxdb = max(maxdb, dualbound[s,pidx]);
               minpb = min(minpb, primalbound[s,pidx]);
            }
            else if( status[s,pidx] == "timeout" )
            {
               ntimeouts[s]++;
               marker = ">";
            }
            else
            {
               nfails[s]++;
               fail = 1;
               marker = "!";
            }

            # print statistics
            printf(" %10d %s%6d", nodes[s,pidx], marker, time[s,pidx]);
            if( s == 0 )
            {
               nodecomp = nodes[s,pidx];
               timecomp = time[s,pidx];
            }
            else
            {
               if( nodes[s,pidx]/nodecomp > 999.99 )
                  printf("  Large");
               else
                  printf(" %6.2f", nodes[s,pidx]/nodecomp);
               if( time[s,pidx]/timecomp > 999.99 )
                  printf("  Large");
               else
                  printf(" %6.2f", time[s,pidx]/timecomp);
            }
         }

         # check for inconsistency in the primal and dual bounds
         if( fail )
            printf("  fail");
         else if( maxdb - minpb > 1e-5 * max(max(abs(maxdb), abs(minpb)), 1.0) )
         {
            printf("  inconsistent");
            fail = 1;
         }
         else
            printf("  ok");
         printf("\n");

         # calculate totals and means for instances where no solver failes
         if( !fail )
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
            }
         }
      }
   }
   printhline(nsolver);

   # print solvers' overall statistics
   probnumstr = "("nevalprobs")";
   printf("%-12s %5s", "total", probnumstr);
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf(" %10d %7d", nodetotal[s], timetotal[s]);
      else
         printf(" %10d %7d              ", nodetotal[s], timetotal[s]);
   }
   printf("\n");
   printf("%-18s", "geom. mean");
   for( s = 0; s < nsolver; ++s )
   {
      if( s == 0 )
         printf(" %10d %7d", nodegeom[s], timegeom[s]);
      else
         printf(" %10d %7d %6.2f %6.2f", nodegeom[s], timegeom[s], nodegeom[s]/nodegeom[0], timegeom[s]/timegeom[0]);
   }
   printf("\n");
   printf("%-18s", "shifted geom.");
   for( s = 0; s < nsolver; ++s )
   {
      nodeshiftedgeom[s] -= nodegeomshift;
      timeshiftedgeom[s] -= timegeomshift;
      if( s == 0 )
         printf(" %10d %7d", nodeshiftedgeom[s], timeshiftedgeom[s]);
      else
         printf(" %10d %7d %6.2f %6.2f", nodeshiftedgeom[s], timeshiftedgeom[s],
            nodeshiftedgeom[s]/nodeshiftedgeom[0], timeshiftedgeom[s]/timeshiftedgeom[0]);
   }
   printf("\n");
   printhline(nsolver);

   printf("\n");
   for( s = 0; s < nsolver; ++s )
   {
      printf("%-30s fails: %4d  timeouts: %4d  solved: %4d  time: %6.1f  shtime: %6.1f  timeQ: %5.2f  shtimeQ: %5.2f\n", 
         solvername[s], nfails[s], ntimeouts[s], nsolved[s], timegeom[s], timeshiftedgeom[s],
         timegeom[s]/timegeom[0], timeshiftedgeom[s]/timeshiftedgeom[0]);
   }
}
