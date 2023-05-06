#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    average.awk
#@brief   compute averages of several MIPLIB result files - can be used for getting result for permuted instances
#@author  Gerald Gamrath
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
   infinity = 1e+20;

   # initialize summary data
   nsolved = 0;
   nstopped = 0;
   nfailed = 0;

   nruns = 0;
   nprobs[nruns] = 0;
   problistlen = 0;
}
/^@03 MIPLIB script version: / {
   version[nruns] = $5;
}
/^@02 timelimit: / {
   timelimit[nruns] = $3;
}
/^@01/ {
   solver[nruns] = $2;
   nruns++;
   nprobs[nruns] = 0;
}
// {
   if ( NF >= 8 )
   {
      statuses["ok"] = 1;
      statuses["abort"] = 1;
      statuses["stopped"] = 1;

      solstatuses["ok"] = 1;
      solstatuses["fail"] = 1;
      solstatuses["mismatch"] = 1;
      solstatuses["error"] = 1;
      solstatuses["--"] = 1;

      probname = $1;
      sub(/\#p([0-9])*$/, "", probname);
      name[nruns,nprobs[nruns]] = probname;

      if ( $8 in solstatuses )
      {
	 # collect data (line with problem type, original and presolved problem size and simplex iterations)
	 dualbound[nruns,nprobs[nruns]] = $2;
	 primalbound[nruns,nprobs[nruns]] = $3;
	 gap[nruns,nprobs[nruns]] = $4;
	 nodes[nruns,nprobs[nruns]] = $5;
	 time[nruns,nprobs[nruns]] = $6;
	 status[nruns,nprobs[nruns]] = $7;
	 solstatus[nruns,nprobs[nruns]] = $8;

	 probidx[probname,nruns] = nprobs[nruns];
	 if ( probname in probcnt )
	    probcnt[probname]++;
	 else
	 {
	    probcnt[probname] = 1;
	    problist[problistlen] = probname;
	    problistlen++;
	 }
	 nprobs[nruns]++;
      }
   }
}
END {
   if ( nruns == 0 )
   {
      printf("No instances found in log files.\n");
      exit 1;
   }

   # check whether time limits, solvers, and MIPLIB script versions are the same
   t = timelimit[0];
   s = solver[0];
   v = version[0];
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
      if ( version[i] != v )
      {
	 printf("MIPLIB script versions of the runs are different.\n");
	 exit 1;
      }
   }

   printf("----------------------------+----------------+----------------+------+---------+-------+-------------+--------------\n");
   printf("Name                        |   Dual Bound   |  Primal Bound  | Gap%% |  Nodes  |  Time |    Status   |    Solution \n");
   printf("----------------------------+----------------+----------------+------+---------+-------+-------------+--------------\n");

   # display the mean values
   for (i = 0; i < problistlen; ++i)
   {
      prob = problist[i];

      # check whether instance has been processed
      if ( ! ((prob,0) in probidx) )
      {
	 printf("Problem <%s> not processed by run 0.\n", prob);
	 exit 1;
      }

      # get data for first instance to compare with
      pidx = probidx[prob,0];

      # initialize status and solution status
      statusfirst = status[0,pidx];
      finalstatus = statusfirst;
      finalstatusnr = 1;

      solstatusfirst = solstatus[0,pidx];
      finalsolstatus = solstatusfirst;
      finalsolstatusnr = 1;

      # initialize average values (note that the presolve # of conss and vars might differ between runs)
      avgdb = dualbound[0,pidx];
      avgpb = primalbound[0,pidx];
      avggap = gap[0,pidx];
      avgnodes = nodes[0,pidx];
      avgtime = time[0,pidx];

      # loop through runs
      for (s = 1; s < nruns; ++s)
      {
	 pidx = probidx[prob,s];

	 # check whether instance has been processed
	 if ( pidx == "" )
	 {
	    printf("Problem <%s> not processed by run %d.\n", prob, s);
	    exit 1;
	 }

	 # update status
	 if ( status[s,pidx] == finalstatus )
	    finalstatusnr += 1;
	 else
	 {
	    if ( finalstatus != "abort" )
	    {
	       if ( status[s,pidx] == "abort" )
	       {
		  finalstatus = status[s,pidx];
		  finalstatusnr = 1;
	       }
	       else
	       {
		  if ( finalstatus != "stopped" )
		  {
		     if ( status[s,pidx] == "stopped" )
		     {
			finalstatus = status[s,pidx];
			finalstatusnr = 1;
		     }
		  }
	       }
	    }
	 }

	 # update solution status
	 if ( solstatus[s,pidx] == finalsolstatus )
	    finalsolstatusnr += 1;
	 else
	 {
	    if ( finalsolstatus != "error" )
	    {
	       if ( solstatus[s,pidx] == "error" )
	       {
		  finalsolstatus = solstatus[s,pidx];
		  finalsolstatusnr = 1;
	       }
	       else
	       {
		  if ( finalsolstatus != "mismatch" )
		  {
		     if ( solstatus[s,pidx] == "mismatch" )
		     {
			finalsolstatus = solstatus[s,pidx];
			finalsolstatusnr = 1;
		     }
		  }
		  else
		  {
		     if ( finalsolstatus != "fail" )
		     {
			if ( solstatus[s,pidx] == "fail" )
			{
			   finalsolstatus = solstatus[s,pidx];
			   finalsolstatusnr = 1;
			}
		     }
		     else
		     {
			if ( finalsolstatus != "ok" )
			{
			   if ( solstatus[s,pidx] == "ok" )
			   {
			      finalsolstatus = solstatus[s,pidx];
			      finalsolstatusnr = 1;
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
	 if ( gap[s,pidx] == "Large" )
	    avggap = infinity;
	 else if ( gap[s,pidx] < infinity )
	 {
	    if ( gap[s,pidx] == "" || gap[s,pidx] == "--" )
	       avggap = "--";
	    else if ( avggap != "--" )
	       avggap += gap[s,pidx];
	 }
	 else
	    avggap = infinity;

	 # the other values should all be finite
	 avgnodes += nodes[s,pidx];
	 avgtime += time[s,pidx];
      }

      # final computation of average values
      if ( avgdb > -infinity && avgdb < infinity )
	 avgdb /= nruns;

      if ( avgpb > -infinity && avgpb < infinity )
	 avgpb /= nruns;

      if ( avggap > -infinity && avggap < infinity )
	 avggap /= nruns;

      avgnodes /= nruns;
      avgtime /= nruns;

      # output
      if ( avggap < 0.0 )
	 gapstr = "  --  ";
      else if( avggap < 1e+04 )
	 gapstr = sprintf("%6.1f", avggap);
      else
	 gapstr = " Large";

      statusnrstr = sprintf("(%d)", finalstatusnr);
      solstatusnrstr = sprintf("(%d)", finalsolstatusnr);

      printf("%-28s %16.9g %16.9g %6s %9d %7d %8s %4s %9s %4s\n", prob, avgdb, avgpb, gapstr, avgnodes, avgtime, finalstatus, statusnrstr, finalsolstatus, solstatusnrstr);

      # determine overall status from solving status and solution status:

      # instance solved correctly (including case that no solution was found)
      if( finalstatus == "ok" && (finalsolstatus == "ok" || finalsolstatus == "--") )
	 nsolved++;
      # incorrect solving process or infeasible solution (including errors with solution checker)
      else if( finalstatus == "abort" || (finalsolstatus == "fail" || finalsolstatus == "error" || finalsolstatus == "mismatch") )
	 nfailed++;
      # stopped due to imposed limits
      else
	 nstopped++;
   }
   printf("----------------------------+----------------+----------------+------+---------+-------+-------------+--------------\n");
   printf("\n");
   printf("solved/stopped/failed: %d/%d/%d\n", nsolved, nstopped, nfailed);
   printf("\n");
   printf("@04 aggregated results of %d runs\n", nruns);
   printf("@03 MIPLIB script version %s\n", version[0]);
   printf("@02 timelimit: %g\n", timelimit[0]);
   printf("@01 %s\n", solver[0]);
}
