#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
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
#@file    average.awk
#@brief   compute averages of several SCIP result files - can be used for getting result for permuted instances
#@author  Marc Pfetsch
#
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

function fracceil(x,f)
{
   return ceil(x/f)*f;
}

BEGIN {
   printsoltime = 1;            # output time to first and optimal solution
   printrounded = 0;            # output rounded average values for integral values
   namelength = 18;             # maximal length of instance names (can be increased)

   infinity = 1e+20;
   mintime = 0.5;
   timegeomshift = 10.0;
   nodegeomshift = 100.0;

   nruns = 0;
   nprobs[nruns] = 0;
   problistlen = 0;
}
/^@02 timelimit: / {
   timelimit[nruns] = $3;
}
/^@01/ {
   # determine solver and settings
   githash[nruns] = $4;

   solver[nruns] = $2;
   sub(/:.*/, "", solver[nruns]);

   settings[nruns] = $2;
   sub(/.*:/, "", settings[nruns]);

   nruns++;
   nprobs[nruns] = 0;
}
// {
   if ( NF >= 13 )
   {
      statuses["ok"] = 1;
      statuses["timeout"] = 1;
      statuses["unknown"] = 1;
      statuses["abort"] = 1;
      statuses["fail"] = 1;
      statuses["readerror"] = 1;
      statuses["better"] = 1;
      statuses["solved"] = 1;
      statuses["sollimit"] = 1;
      statuses["gaplimit"] = 1;
      statuses["memlimit"] = 1;
      statuses["nodelimit"] = 1;

      name[nruns,nprobs[nruns]] = $1;
      validline = 0;

      if ( NF >= 13 && $13 in statuses ) # GLPK, CPLEX, SCIP without columns displaying times to first and best solution
      {
	 # collect data (line with problem type, original and presolved problem size and simplex iterations)
	 type[nruns,nprobs[nruns]] = $2;
	 origconss[nruns,nprobs[nruns]] = $3;
	 origvars[nruns,nprobs[nruns]] = $4;
	 conss[nruns,nprobs[nruns]] = $5;
	 vars[nruns,nprobs[nruns]] = $6;
	 dualbound[nruns,nprobs[nruns]] = max(min($7, +infinity), -infinity);
	 primalbound[nruns,nprobs[nruns]] = max(min($8, +infinity), -infinity);
	 gap[nruns,nprobs[nruns]] = $9;
	 iters[nruns,nprobs[nruns]] = $10;
	 nodes[nruns,nprobs[nruns]] = max($11,1);
	 time[nruns,nprobs[nruns]] = fracceil(max($12,mintime),0.1);
	 status[nruns,nprobs[nruns]] = $13;
	 printsoltimes = 0; # additional output is only available for SCIP-.res files
	 validline = 1;
      }

      if ( NF >= 15 && $15 in statuses ) # SCIP with solution times to first/last
      {
	 # collect data (line with problem type, original and presolved problem size and simplex iterations)
	 type[nruns,nprobs[nruns]] = $2;
	 origconss[nruns,nprobs[nruns]] = $3;
	 origvars[nruns,nprobs[nruns]] = $4;
	 conss[nruns,nprobs[nruns]] = $5;
	 vars[nruns,nprobs[nruns]] = $6;
	 dualbound[nruns,nprobs[nruns]] = max(min($7, +infinity), -infinity);
	 primalbound[nruns,nprobs[nruns]] = max(min($8, +infinity), -infinity);
	 gap[nruns,nprobs[nruns]] = $9;
	 iters[nruns,nprobs[nruns]] = $10;
	 nodes[nruns,nprobs[nruns]] = max($11,1);
	 time[nruns,nprobs[nruns]] = fracceil(max($12,mintime),0.1);
	 timetofirst[nruns,nprobs[nruns]] = fracceil(max($13,mintime),0.1);
	 timetobest[nruns, nprobs[nruns]] = fracceil(max($14, mintime), 0.1);
	 status[nruns,nprobs[nruns]] = $15;
	 validline = 1;
      }

      if ( validline )
      {
	 # postprocessing of information
	 if ( status[nruns,nprobs[nruns]] == "better" )
	    status[nruns,nprobs[nruns]] = "timeout";
	 if ( status[nruns,nprobs[nruns]] == "sollimit" || status[nruns,nprobs[nruns]] == "gaplimit" || status[nruns,nprobs[nruns]] == "solved" )
	    status[nruns,nprobs[nruns]] = "ok";

	 probidx[$1,nruns] = nprobs[nruns];
	 if ( $1 in probcnt )
	    probcnt[$1]++;
	 else
	    probcnt[$1] = 1;
	 nprobs[nruns]++;
	 if ( probcnt[$1] == 1 )
	 {
	    problist[problistlen] = $1;
	    problistlen++;
	 }
      }
   }
}
END {
   if ( nruns == 0 )
   {
      printf("No instances found in log files.\n");
      exit 1;
   }

   # check whether time limits, solvers, and git hashes are the same
   t = timelimit[0];
   s = solver[0];
   h = githash[0];
   for (i = 1; i < nruns; ++i)
   {
      if ( timelimit[i] != t )
      {
	 printf("Time limits of the runs are different.\n");
	 exit 1;
      }
      if ( solver[i] != s )
      {
	 printf("Solvers of the runs are different.\n");
	 exit 1;
      }
      if ( githash[i] != h )
      {
	 printf("Git hashes of the runs are different.\n");
	 exit 1;
      }
   }

   # prepare header
   hyphenstr = "";
   for (i = 0; i < namelength; ++i)
      hyphenstr = sprintf("%s-", hyphenstr);

   # first part: name of given length
   tablehead1 = hyphenstr;
   tablehead2 = sprintf("Name%*s", namelength-4, " ");
   tablehead3 = hyphenstr;

   # append rest of header
   if ( printrounded )
   {
      tablehead1 = tablehead1"+------+--- Original --+-- Presolved --+----------------+----------------+------+---------+--------+-------+";
      tablehead2 = tablehead2"| Type | Conss |  Vars | Conss |  Vars |   Dual Bound   |  Primal Bound  | Gap%% |  Iters  |  Nodes |  Time |";
      tablehead3 = tablehead3"+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+";
   }
   else
   {
      tablehead1 = tablehead1"+------+--- Original --+---- Presolved ----+----------------+----------------+------+-----------+----------+-------+";
      tablehead2 = tablehead2"| Type | Conss |  Vars |   Conss |    Vars |   Dual Bound   |  Primal Bound  | Gap%% |   Iters   |   Nodes  |  Time |";
      tablehead3 = tablehead3"+------+-------+-------+---------+---------+----------------+----------------+------+-----------+----------+-------+";
   }

   tablehead1 = tablehead1"--------\n";
   tablehead2 = tablehead2"       \n";
   tablehead3 = tablehead3"--------\n";

   printf(tablehead1);
   printf(tablehead2);
   printf(tablehead3);

   # init averages over all instances
   nodegeom = 0.0;
   timegeom = 0.0;
   shiftednodegeom = nodegeomshift;
   shiftedtimegeom = timegeomshift;
   stotnodes = 0.0;
   stottime = 0.0;

   passes = 0;
   timeouts = 0;
   fails = 0;

   # display the mean values
   for (i = 0; i < problistlen; ++i)
   {
      prob = problist[i];

      if ( length(prob) > namelength )
	 shortprob = substr(prob, length(prob)-namelength-1, namelength);
      else
	 shortprob = prob;

      line = sprintf("%-18s", shortprob);

      # check whether instance has been processed
      if ( ! ((prob,0) in probidx) )
      {
	 printf("Problem <%s> not processed by run 0.\n", prob);
	 exit 1;
      }

      # get data for first instance to compare with
      pidx = probidx[prob,0];

      # check whether the following values are the same for each run
      typefirst = type[0,pidx];
      origconssfirst = origconss[0,pidx];
      origvarsfirst = origvars[0,pidx];
      statusfirst = status[0,pidx];
      finalstatus = statusfirst;
      finalstatusnr = 1;

      # initialize average values (note that the presolve # of conss and vars might differ between runs)
      avgconss = conss[0,pidx];
      avgvars = vars[0,pidx];
      avgdb = dualbound[0,pidx];
      avgpb = primalbound[0,pidx];
      avggap = gap[0,pidx];
      avgiters = iters[0,pidx];
      avgnodes = nodes[0,pidx];
      avgtime = time[0,pidx];
      if ( printsoltimes )
      {
	 avgtimetobest = timetobest[0,pidx];
	 avgtimetofirst = timetofirst[0,pidx];
      }

      # loop through runs
      for (s = 1; s < nruns; ++s)
      {
	 pidx = probidx[prob,s];

	 # check whether instance has been processed
	 if ( pidx == "" )
	 {
	    printf("Problem <%s> not processed by instance %d.\n", prob, s);
	    exit 1;
	 }

	 # compare statistics to first instance
	 # problem types are a bit fuzzy and dervided from presolve instance, so they can change from one permutation to another
	 #if ( type[s,pidx] != typefirst )
	 #{
	 #   printf("Warning: Problem <%s> type not equal between runs.\n", prob);
	 #   exit 1;
	 #}
	 if ( origconss[s,pidx] != origconssfirst )
	 {
	    printf("Error: Problem <%s> number of constraints not equal between runs (%d != %d (%d)).\n", prob, conssfirst, conss[s,pidx], s);
	    exit 1;
	 }
	 if ( origvars[s,pidx] != origvarsfirst )
	 {
	    printf("Error: Problem <%s> number of variables not equal between runs.\n", prob);
	    exit 1;
	 }
	 if ( status[s,pidx] == finalstatus )
	    finalstatusnr += 1;
	 else
	 {
	    if ( finalstatus != "fail" && finalstatus != "abort" && finalstatus != "readererror" )
	    {
	       if ( status[s,pidx] == "fail" || status[s,pidx] == "abort" || status[s,pidx] == "readererror" )
	       {
		  finalstatus = status[s,pidx];
		  finalstatusnr = 1;
	       }
	       else
	       {
		  if ( finalstatus != "sollimit" && finalstatus != "gaplimit" && finalstatus != "memlimit" && finalstatus != "nodelimit" )
		  {
		     if ( status[s,pidx] == "sollimit" || status[s,pidx] == "gaplimit" || status[s,pidx] == "memlimit" || status[s,pidx] == "nodelimit" )
		     {
			finalstatus = status[s,pidx];
			finalstatusnr = 1;
		     }
		     else
		     {
			if ( finalstatus != "timeout" )
			{
			   if ( status[s,pidx] == "timeout" )
			   {
			      finalstatus = status[s,pidx];
			      finalstatusnr = 1;
			   }
			}
		     }
		  }
	       }
	    }
	 }

	 # ---------------------------------------------
	 # compute average values

	 # take care of infinite values
	 if ( avgdb > -infinity && avgdb < infinity )
	 {
	    if ( dualbound[s,pidx] <= -infinity )
	       avgdb = -infinity;
	    else if ( dualbound[s,pidx] >= infinity )
	       avgdb = -infinity;
	    else
	       avgdb += dualbound[s,pidx];
	 }

	 if ( avgpb > -infinity && avgpb < infinity )
	 {
	    if ( primalbound[s,pidx] <= -infinity )
	       avgpb = infinity;
	    else if ( primalbound[s,pidx] >= infinity )
	       avgpb = infinity;
	    else
	       avgpb += primalbound[s,pidx];
	 }

	 # take care of gap
	 if ( gap[s,pidx] == "" || gap[s,pidx] == "--" || gap[s,pidx] == "Large" )
	    avggap = infinity;
	 else if ( gap[s,pidx] < infinity )
	    avggap += gap[s,pidx];
	 else
	    avggap = infinity;

	 # the other values should all be finite
	 avgconss += conss[s,pidx];
	 avgvars += vars[s,pidx];
	 avgiters += iters[s,pidx];
	 avgnodes += nodes[s,pidx];
	 avgtime += time[s,pidx];
	 if ( printsoltimes )
	 {
	    avgtimetobest += timetobest[s,pidx];
	    avgtimetofirst += timetofirst[s,pidx];
	 }
      }

      # final computation of average values
      if ( avgdb > -infinity && avgdb < infinity )
	 avgdb /= nruns;

      if ( avgpb > -infinity && avgpb < infinity )
	 avgpb /= nruns;

      if ( avggap > -infinity && avggap < infinity )
	 avggap /= nruns;

      avgconss /= nruns;
      avgvars /= nruns;
      avgiters /= nruns;
      avgnodes /= nruns;
      avgtime /= nruns;
      if ( printsoltimes )
      {
	 avgtimetobest /= nruns;
	 avgtimetofirst /= nruns;
      }

      # output
      if ( avggap < 0.0 )
	 gapstr = "  --  ";
      else if( avggap < 1e+04 )
	 gapstr = sprintf("%6.1f", avggap);
      else
	 gapstr = " Large";

      if ( printrounded )
      {
	 printf("%-*s  %-5s %7d %7d %7d %7d %16.9g %16.9g %6s %9d %8d %7.1f ",
		namelength, shortprob, typefirst, origconssfirst, origvarsfirst, avgconss, avgvars, avgdb, avgpb, gapstr, avgiters, avgnodes, avgtime);
      }
      else
      {
	 printf("%-*s  %-5s %7d %7d %9.1f %9.1f %16.9g %16.9g %6s %11.1f %10.1f %7.1f ",
		namelength, shortprob, typefirst, origconssfirst, origvarsfirst, avgconss, avgvars, avgdb, avgpb, gapstr, avgiters, avgnodes, avgtime);
      }

      if ( printsoltimes )
	 printf(" %9.1f %9.1f ", abgtimetofirst, avgtimetobest);

      printf("%s (%d)\n", finalstatus, finalstatusnr);

      # determine count as fails/timeout/passed
      if ( finalstatus == "timeout" || finalstatus == "sollimit" || finalstatus == "gaplimit" || finalstatus == "memlimit" || finalstatus == "nodelimit" )
	 timeouts += 1;

      if ( finalstatus == "ok" || finalstatus == "better" || finalstatus == "unknown" || finalstatus == "solved" )
	 passes += 1;

      if ( finalstatus == "fail" || finalstatus == "abort" || finalstatus == "readererror" )
	 fails += 1;

      # compute averages over all instances
      nodegeom = nodegeom^(i/(i+1)) * max(avgnodes, 1.0)^(1.0/(i+1));
      timegeom = timegeom^(i/(i+1)) * max(avgtime, 1.0)^(1.0/(i+1));
      shiftednodegeom = shiftednodegeom^(i/(i+1)) * max(avgnodes + nodegeomshift, 1.0)^(1.0/(i+1));
      shiftedtimegeom = shiftedtimegeom^(i/(i+1)) * max(avgtime + timegeomshift, 1.0)^(1.0/(i+1));
      stotnodes += avgnodes;
      stottime += avgtime;
   }
   if ( printrounded )
      printf("------------------+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+--------\n\n");
   else
      printf("------------------+------+-------+-------+---------+---------+----------------+----------------+------+---------+----------+---------+--------\n\n");

   # output averages over all instances
   tablebottom1 = "------------------------------[Nodes]---------------[Time]------";
   tablebottom2 = "  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom.";
   tablebottom3 = "----------------------------------------------------------------";

   tablebottom1 = tablebottom1"\n";
   tablebottom2 = tablebottom2"\n";
   tablebottom3 = tablebottom3"\n";

   printf(tablebottom1);
   printf(tablebottom2);
   printf(tablebottom3);

   shiftednodegeom -= nodegeomshift;
   shiftedtimegeom -= timegeomshift;

   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f ", problistlen, passes, timeouts, fails, stotnodes / 1000.0, nodegeom, stottime, timegeom);

   printf("\n");
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f ",
	  nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   printf("\n");
   printf(tablebottom3);

   # try to clean up settings (if used for permuted runs)
   setting = settings[1];
   sub(/-p([0-9])*$/, "", setting);

   printf("@02 timelimit: %g\n", timelimit[0]);
   printf("@01 %s:%s [GitHash: %s\n", solver[0], setting, githash[0]);
}
