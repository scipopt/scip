#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}

function printhline(nsolver)
{
   printf("-----------------------------------");
   for( s = 0; s < nsolver; ++s )
   {
      printf("+");

      for( c = 0; c <= colwidth; ++c )
	 printf("-");
   }
   printf("+\n");
}

function printsolvers(nsolver)
{
   printf("  Name                             |");
   #printf("%31s", "|");
   for( s = 0; s < nsolver; ++s )
   {
      nsplits = split(solvername[s], splitarray, "(");
      shortname[s] = splitarray[1];
      if( nsplits > 2 )
      {
	 nsplits = split(splitarray[2], splitarray2, ")");
	 shortname[s] = shortname[s] "-"  splitarray2[2];
      }

      printf("%"colwidth"s |", substr(shortname[s], 1, colwidth));
   }
   printf("\n");
}

BEGIN {

   infinity = 1e+20;
   timegeomshift = 10.0;
   nodegeomshift = 100.0;
   mintime = 1.0;

   problistlen = 0;
   nsolver = 0;
   nprobs[nsolver] = 0;

   colwidth = 14;
}

# collect used time limit
/^@02 timelimit: / {
   timelimit[nsolver] = $3;
}

# collect solver information
/^@01 / {
   solvername[nsolver] = $2;
   nsolver++;
   nprobs[nsolver] = 0;
}

# collect number of permutations
/^@04 aggregated results of/ {
   nperms[nsolver] = $5;
}

# collect instance results
// {
   statuses["ok"];
   statuses["stopped"];
   statuses["abort"];

   solstatuses["ok"];
   solstatuses["fail"];
   solstatuses["error"];
   solstatuses["mismatch"];
   solstatuses["--"];

   # check if this is a useable line
   validline = 0;

   # result file for a single run
   if( $7 in statuses && $8 in solstatuses )
   {
      # collect data from res-file
      name[nsolver,nprobs[nsolver]] = $1;
      time[nsolver,nprobs[nsolver]] = $6;
      status[nsolver,nprobs[nsolver]] = $7;
      solstatus[nsolver,nprobs[nsolver]] = $8;
      validline = 1;
   }

   # result file for average over multiple permutations
   if( $7 in statuses && $9 in solstatuses )
   {
      # collect data from res-file
      name[nsolver,nprobs[nsolver]] = $1;
      time[nsolver,nprobs[nsolver]] = $6;
      status[nsolver,nprobs[nsolver]] = $7;
      solstatus[nsolver,nprobs[nsolver]] = $9;
      validline = 1;
   }

   if( validline )
   {
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

# collect overall result
/^solved\/stopped\/failed:/ {
   overall[nsolver] = $2
}

END {

   # print headers
   printhline(nsolver);
   printsolvers(nsolver);
   printhline(nsolver);

   # display the problem results
   for( i = 0; i < problistlen; ++i )
   {
      p = problist[i];
      if( length(p) > 28 )
         shortp = substr(p, length(p)-32, 33);
      else
         shortp = p;

      line = sprintf("%-33s  ", shortp);

      # write results of all solvers for the current instance
      for( s = 0; s < nsolver; ++s )
      {
         pidx = probidx[p,s];
         processed = (pidx != "");

	 # if the solver exceeded the timelimit, set status accordingly ???????
         if( (status[s,pidx] == "ok") && timelimit[s] > 0.0 && time[s,pidx] > timelimit[s] )
         {
            status[s,pidx] = "stopped";
            time[s,pidx] = timelimit[s];
         }

         if( processed && name[s,pidx] != p )
            printf("Error: solver %d, probidx %d, <%s> != <%s>\n", solvername[s], pidx, name[s,pidx], p);

         if( processed )
         {
	    if( solstatus[s,pidx] == "ok" || solstatus[s,pidx] == "--" )
	    {
	       if( status[s,pidx] == "ok" )
		  line = sprintf("%s %"colwidth"d ", line, time[s,pidx]);
	       else if( status[s,pidx] == "stopped" )
		  line = sprintf("%s %"colwidth"s ", line, "timeout");
	       else
		  line = sprintf("%s %"colwidth"s ", line, status[s,pidx]);
	    }
	    else
	       line = sprintf("%s %"colwidth"s ", line, solstatus[s,pidx]);
	 }
	 else
	    line = sprintf("%s %"colwidth"s ", line, "-");

      }

      # display instance information
      printf("%s\n", line);
   }

   printhline(nsolver);

   # write results of all solvers for the current instance
   line = sprintf("%-33s  ", "solved/stopped/failed");

   for( s = 0; s < nsolver; ++s )
   {
      line = sprintf("%s|%"colwidth"s ", line, overall[s]);
   }
   printf("%s|\n", line);

   line = sprintf("%-33s  ", "timelimit [sec]");

   for( s = 0; s < nsolver; ++s )
   {
      line = sprintf("%s|%"colwidth"s ", line, timelimit[s]);
   }
   printf("%s|\n", line);

   printhline(nsolver);
}
