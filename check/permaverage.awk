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
#@file    permaverage.awk
#@brief   compute averages over instances for different permuations
#@author  Jan Kuske
#@author  Marc Pfetsch
#
# This script creates the average values between different
# permutations over the same setting, solver, and instances. There are
# some switches to give additional information. They are listed and
# documented in the begin-block. This script will not interrupt if the
# type of a problem changes between the permutations. Instead it will
# print "**" in the column "Type". The status is shown as "status (n)",
# the number n tells how many times the status appeared. Only
# the Status with the highest priority will be counted.
#
# The priorities are:
#  1. fail / abort /readerror
#  2. some limit is reached except timeout
#  3. timeout
#  4. ok
#
# If some permutations are missing, the script will proceed and print
# a list with all instances with missing permutations. The same
# happens with failed/abort/readerror status permutations

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

function median(c,v,  j)
{
   asort(v,j);
   if (c % 2)
      return j[(c+1)/2];
   else
      return (j[c/2+1]+j[c/2])/2.0;
}

BEGIN {
   printrounded = 0;             # output rounded average values for integral values
   printmediantime = 0;          # output median time
   printstdtime = 0;             # output standard devation time
   namelength = 35;              # maximal length of instance names (can be increased)

   skipfails = 0;                # if true: the script will proceed, even if some permutations failed
				 # if false: if one permutation failed, the whole instance will be ignored
   missingpem = 0;               # to verify if there is at least one permutation missing
   failedpem = 0;                # to verify if there is at least one permutation failed

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

      validline = 0;

      if ( NF >= 13 && $13 in statuses ) # GLPK, CPLEX, SCIP without columns displaying times to first and best solution
      {
	 # collect data (line with problem type, original and presolved problem size and simplex iterations)
	 name[nruns,nprobs[nruns]] = $1;
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
	    probfirstrun[$1] = nruns;
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

   tablehead1 = tablehead1"------------+";
   tablehead2 = tablehead2" Delta Time |";
   tablehead3 = tablehead3"------------+";

   if (printstdtime)
   {
      tablehead1 = tablehead1"----------+";
      tablehead2 = tablehead2" StddTime |";
      tablehead3 = tablehead3"----------+";
   }

   if (printmediantime)
   {
      tablehead1 = tablehead1"----------+";
      tablehead2 = tablehead2"  Median  |";
      tablehead3 = tablehead3"----------+";
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
   k=0; #this is a counter is used to calculate the geom. mean.
   for (i = 0; i < problistlen; ++i)
   {
      prob = problist[i];

      if ( length(prob) > namelength )
	 shortprob = substr(prob, length(prob)-namelength-1, namelength);
      else
	 shortprob = prob;

      line = sprintf("%-18s", shortprob);

      firstrun = probfirstrun[prob];
      pemmiscnt[prob]=firstrun;
      if (firstrun != 0)
	 missingpem = 1;

      # get data for first instance to compare with
      pidx = probidx[prob,firstrun];

      # check whether the following values are the same for each run
      typefirst = type[firstrun,pidx];
      origconssfirst = origconss[firstrun,pidx];
      origvarsfirst = origvars[firstrun,pidx];
      statusfirst = status[firstrun,pidx];
      finalstatus = statusfirst;
      finalstatusnr = 1;
      if ( statusfirst == "fail" || statusfirst == "abort" || statusfirst == "readererror" )
      {
	 failedpem = 1;
	 finalstatus="ok";
	 finalstatusnr = 0;
	 pemfailedcnt[prob]=1;
      }

      # initialize average values (note that the presolve # of conss and vars might differ between runs)
      avgconss = conss[firstrun,pidx];
      avgvars = vars[firstrun,pidx];
      avgdb = dualbound[firstrun,pidx];
      avgpb = primalbound[firstrun,pidx];
      avggap = gap[firstrun,pidx];
      avgiters = iters[firstrun,pidx];
      avgnodes = nodes[firstrun,pidx];
      avgtime = time[firstrun,pidx];
      mintime = time[firstrun,pidx];
      maxtime = time[firstrun,pidx];

      # loop through runs
      for (s = firstrun+1; s < nruns; ++s)
      {
	 pidx = probidx[prob,s];

	 # check whether instance has been processed
	 if ( pidx == "" )
	 {
	    pemmiscnt[prob]++;
	    missingpem = 1;
	    continue;
	 }

	 if ( status[s,pidx] == "fail" || status[s,pidx] == "abort" || status[s,pidx] == "readererror" )
	 {
	    failedpem = 1;
	    pemfailedcnt[prob]++;
	    if (skipfails)
	    {
	       continue;
	    }
	    else
	    {
	       break;
	    }
	 }

	 # compare statistics to first instance
	 if ( type[s,pidx] != typefirst )
	 {
	    typerfirst = "**";   #if the type changes between runs, we write "**" into the type column.
	 }

	 if ( origconss[s,pidx] != origconssfirst )
	 {
	    printf("Error: Problem <%s> number of constraints not equal between runs (%d != %d (%d)).\n", prob, origconssfirst, origconss[s,pidx], s);
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

	 if (maxtime < time[s,pidx])
	    maxtime = time[s,pidx];
	 if (mintime > time[s,pidx])
	    mintime = time[s,pidx];
      }

      if (pemfailedcnt[prob] != 0 && !skipfails)
      {
	 continue;
      }

      if (pemfailedcnt[prob] != probcnt[prob])
      {
	 probcnt[prob]=probcnt[prob]-pemfailedcnt[prob];
      }
      else
      {
	 continue;
      }

      # final computation of average values
      if ( avgdb > -infinity && avgdb < infinity )
	 avgdb /=probcnt[prob];

      if ( avgpb > -infinity && avgpb < infinity )
	 avgpb /= probcnt[prob];

      if ( avggap > -infinity && avggap < infinity )
	 avggap /= probcnt[prob];

      avgconss /= probcnt[prob];
      avgvars /= probcnt[prob];
      avgiters /= probcnt[prob];
      avgnodes /= probcnt[prob];
      avgtime /= probcnt[prob];
      difftime = maxtime - mintime;

      stddtime = 0;
      tempcnt=0;
      for (s = 0; s < nruns; ++s)
      {
	 pidx = probidx[prob,s];
	 if ( pidx != "" )
	 {
	    stddtime = stddtime + (avgtime - time[s,pidx])^2;
	    medtimeList[tempcnt] = time[s,pidx];
	    tempcnt++;
	 }
      }
      if ( probcnt[prob] > 1 )
	 stddtime = (1/(probcnt[prob]-1))*stddtime;

      stddtime = sqrt(stddtime);
      medtime =  median(tempcnt,medtimeList);

      # output
      if ( avggap < 0.0 )
	 gapstr = "  --  ";
      else if ( avggap < 1e+04 )
	 gapstr = sprintf("%6.1f", avggap);
      else
	 gapstr = " Large";

      if ( printrounded )
      {
	 printf("%-*s  %-5s %7d %7d %7d %7d %16.9g %16.9g %6s %9d %8d %7.1f %12.1f ",
	 namelength, shortprob, typefirst, origconssfirst, origvarsfirst, avgconss, avgvars, avgdb, avgpb, gapstr, avgiters, avgnodes, avgtime,difftime);
      }
      else
      {
	 printf("%-*s  %-5s %7d %7d %9.1f %9.1f %16.9g %16.9g %6s %11.1f %10.1f %7.1f %12.1f ",
	 namelength, shortprob, typefirst, origconssfirst, origvarsfirst, avgconss, avgvars, avgdb, avgpb, gapstr, avgiters, avgnodes, avgtime,difftime);
      }

      if (printstdtime)
      {
	 printf("%10.1f ",stddtime);
      }
      if (printmediantime)
      {
	 printf("%10.1f ",medtime);
      }

      printf("%s (%d)\n", finalstatus, finalstatusnr);

      # determine count as fails/timeout/passed
      if ( finalstatus == "timeout" || finalstatus == "sollimit" || finalstatus == "gaplimit" || finalstatus == "memlimit" || finalstatus == "nodelimit" )
	 timeouts += 1;

      if ( finalstatus == "ok" || finalstatus == "better" || finalstatus == "unknown" || finalstatus == "solved" )
	 passes += 1;

      if ( finalstatus == "fail" || finalstatus == "abort" || finalstatus == "readererror" )
	 fails += 1;

      # compute averages over all instances
      nodegeom = nodegeom^(k/(k+1)) * max(avgnodes, 1.0)^(1.0/(k+1));
      timegeom = timegeom^(k/(k+1)) * max(avgtime, 1.0)^(1.0/(k+1));
      shiftednodegeom = shiftednodegeom^(k/(k+1)) * max(avgnodes + nodegeomshift, 1.0)^(1.0/(k+1));
      shiftedtimegeom = shiftedtimegeom^(k/(k+1)) * max(avgtime + timegeomshift, 1.0)^(1.0/(k+1));
      k++;
      stotnodes += avgnodes;
      stottime += avgtime;
   }
   printf(tablehead3)

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

   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f ",passes+timeouts+fails , passes, timeouts, fails, stotnodes / 1000.0, nodegeom, stottime, timegeom);

   printf("\n");
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f ",nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   printf("\n");
   printf(tablebottom3);

   if (missingpem)
   {
      printf("\n");
      printf(tablebottom3);
      printf("Missing Permutations \n");
      printf(tablebottom3);
      printf("%-18s %s\n","Instance","Permutations");
      printf(tablebottom3);

      for (i = 0; i < problistlen; ++i)
      {
	 prob = problist[i];
	 if ( length(prob) > namelength )
	    shortprob = substr(prob, length(prob)-namelength-1, namelength);
	 else
	    shortprob = prob;
	 if (pemmiscnt[prob] != 0)
	    printf("%-18s  %2d \n",shortprob,pemmiscnt[prob]);
      }
   }
   if (!skipfails && failedpem)
   {
      printf("\n");
      printf(tablebottom3);
      printf("Ignored Instances with at least one failed permutation \n");
      printf(tablebottom3);
      printf("%-18s \n","Instance");
      printf(tablebottom3);
      for (i = 0; i < problistlen; ++i)
      {
	 prob = problist[i];
	 if ( length(prob) > namelength )
	    shortprob = substr(prob, length(prob)-namelength-1, namelength);
	 else
	    shortprob = prob;
	 if (pemfailedcnt[prob] != 0)
	    printf("%-18s\n",shortprob);
      }
   }
   else if (skipfails && failedpem)
   {
      printf("\n");
      printf(tablebottom3);
      printf("Skipped failed permutations \n");
      printf(tablebottom3);
      printf("%-18s %s\n","Instance","failed");
      printf(tablebottom3);
      for (i = 0; i < problistlen; ++i)
      {
	 prob = problist[i];
	 if ( length(prob) > namelength )
	    shortprob = substr(prob, length(prob)-namelength-1, namelength);
	 else
	    shortprob = prob;

	 if (pemfailedcnt[prob] != 0)
	 {
	    if (pemfailedcnt[prob] != probcnt[prob])
	       printf("%-18s %2d\n",shortprob,pemfailedcnt[prob]);
	    else
	       printf("%-18s %2d (all)\n",shortprob,pemfailedcnt[prob]);
	 }
      }
   }
   printf("\n");

   # try to clean up settings (if used for permuted runs)
   setting = settings[1];
   sub(/-p([0-9])*$/, "", setting);

   printf("@03 Permutations: %d\n",nruns);
   printf("@02 timelimit: %g\n", timelimit[0]);
   printf("@01 %s:%s [GitHash: %s\n", solver[0], setting, githash[0]);
}
