#!/usr/bin/env bash
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

### resets and fills a settings file SETFILE to run SCIP with
### sets correct limits, reads in settings, and controls
### display of the solving process

# environment variables passed as arguments
INSTANCE=$1      #  instance name to solve
SCIPPATH=$2      # - path to working directory for test (usually, the check subdirectory)
SCIP_INSTANCEPATH=$3 # instance path
TMPFILE=$4       # - the batch file to control SCIP
SETNAME=$5       # - specified basename of settings-file, or 'default'
SETFILE=$6       # - instance/settings specific set-file
THREADS=$7       # - the number of LP solver threads to use
SETCUTOFF=$8     # - should optimal instance value be used as objective limit (0 or 1)?
FEASTOL=$9       # - feasibility tolerance, or 'default'
TIMELIMIT=${10}  # - time limit for the solver
MEMLIMIT=${11}   # - memory limit for the solver
NODELIMIT=${12}  # - node limit for the solver
LPS=${13}        # - LP solver to use
DISPFREQ=${14}   # - display frequency for chronological output table
REOPT=${15}      # - true if we use reoptimization, i.e., using a difflist file instead if an instance file
OPTCOMMAND=${16} # - command that should per executed after reading the instance, e.g. optimize, presolve or count
CLIENTTMPDIR=${17}
SOLBASENAME=${18}
SETCUTOFF=${19}
VISUALIZE=${20}  # - true if visualization data should be recorded
SOLUFILE=${21}   # - solu file, only necessary if $SETCUTOFF is 1
#args=("$@")
#for ((i=0; i < $#; i++)) {
#   echo "argument $((i+1)): ${args[$i]}"
#}

# new environment variables after running this script
# -None

#set solfile
SOLFILE=$CLIENTTMPDIR/${USER}-tmpdir/$SOLBASENAME.sol

# reset SETFILE
echo > $SETFILE

# read in settings (even when using default, see bugzilla 600)
SETTINGS=$SCIPPATH/../settings/$SETNAME.set
if test $SETNAME == "default"
then
   # create empty settings file
   test -e $SETTINGS || touch $SETTINGS
fi
cat $SETTINGS >> $SETFILE

# set non-default feasibility tolerance
if test $FEASTOL != "default"
then
    echo numerics/feastol = $FEASTOL >> $SETFILE
fi

# if permutation counter is positive add permutation seed (0 = default)
if test $p -gt 0
then
    echo misc/permutationseed = $p   >> $SETFILE
fi

# avoid solving LPs in case of LPS=none
if test "$LPS" = "none"
then
    echo lp/solvefreq = -1           >> $SETFILE
fi
echo limits/time = $TIMELIMIT        >> $SETFILE
echo limits/nodes = $NODELIMIT       >> $SETFILE
echo limits/memory = $MEMLIMIT       >> $SETFILE
echo lp/threads = $THREADS           >> $SETFILE
echo timing/clocktype = 1            >> $SETFILE
echo display/freq = $DISPFREQ        >> $SETFILE
# avoid switching to dfs - better abort with memory error
echo memory/savefac = 1.0            >> $SETFILE
echo display/statistics = TRUE       >> $SETFILE

# set objective limit: optimal solution value from solu file, if existent
if test $SETCUTOFF = 1
then
    if test $SOLUFILE == ""
    then
        echo Exiting test because no solu file can be found for this test
        exit
    fi
    CUTOFF=`grep "$SHORTPROBNAME " $SOLUFILE | grep -v =feas= | grep -v =inf= | tail -n 1 | awk '{print $3}'`
    if test ""$CUTOFF != ""
    then
        echo limits/objective = $CUTOFF     >> $SETFILE
        echo heuristics/emphasis = \"off\"  >> $SETFILE
    fi
fi

if test "$VISUALIZE" = true
then
    BAKFILENAME="`basename $SETFILE .set`.dat"
    echo visualization output set to "$BAKFILENAME"
    echo visual/bakfilename = "results/${BAKFILENAME}" >> $SETFILE
fi
