#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


# set all solver specific data:
#  solver ["?"]
#  solverversion ["?"]
#  solverremark [""]
#  bbnodes [0]
#  db [-infty]
#  pb [+infty]
#  aborted [1]
#  timeout [0]

# solver name
BEGIN {
    solver = "SAS";
    infeas = 0;
}

/^NOTE: SAS \(r\) Proprietary Software/ {
    solverversion = $6;
}

/SOLUTION_STATUS=/ {
    status = getValue("SOLUTION_STATUS");
    if( status ~ "^TIME_LIM" )
    {
        timeout = 1;
        aborted = 0;
	infeas = 0;
    }
    else if( status ~ "^OPTIMAL" )
    {
        timeout = 0;
        aborted = 0;
	infeas = 0;
    }
    # optimal with slight violations of the tolerances; the solution checker will check the solution anyway
    else if( status ~ "^OPTIMAL_COND" )
    {
        timeout = 0;
        aborted = 0;
	infeas = 0;
    }
    else if( status ~ "^INFEASIBLE" )
    {
        timeout = 0;
        aborted = 0;
	infeas = 1;
	pb = +infty;
	db = +infty;
    }
    else
    {
        timeout = 0;
        aborted = 1;
	infeas = 0;
    }
}

/NODES=/ {
    bbnodes = getValue("NODES");
}

/OBJECTIVE=/ {
   if( infeas == 0 )
   {
      pb = getValue("OBJECTIVE");
   }
}

/BEST_BOUND=/ {
   if( infeas == 0 )
   {
      db = getValue("BEST_BOUND");
   }
}

function getValue(key) {
    for( i = 1; i <= NF; i++ )
    {
        split($i, pair, "=");
        if ( pair[1] == key )
        {
            return pair[2];
        }
    }
    printf("Can't get value for key '%s'", key) > "/dev/stderr";
    exit -1;
}
