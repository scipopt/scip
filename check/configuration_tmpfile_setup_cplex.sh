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

### resets and fills a batch file TMPFILE to run CPLEX with
### sets correct limits, reads in settings, and controls
### display of the solving process

# environment variables passed as arguments
INSTANCE=$1      #  instance name to solve
SCIPPATH=$2      # - path to working directory for test (usually, the check subdirectory)
TMPFILE=$3       # - the batch file to control CPLEX
SETNAME=$4       # - specified basename of settings-file, or 'default'
SETFILE=$5       # - instance/settings specific set-file
THREADS=$6       # - the number of LP solver threads to use
SETCUTOFF=$7     # - should optimal instance value be used as objective limit (0 or 1)?
FEASTOL=$8       # - feasibility tolerance, or 'default'
TIMELIMIT=$9     # - time limit for the solver
MEMLIMIT=${10}   # - memory limit for the solver
NODELIMIT=${11}  # - node limit for the solver
LPS=${12}        # - LP solver to use
DISPFREQ=${13}   # - display frequency for chronological output table
REOPT=${14}      # - true if we use reoptimization, i.e., using a difflist file instead if an instance file
OPTCOMMAND=${15} # - command that should per executed after reading the instance, e.g. optimize, presolve or count
CLIENTTMPDIR=${16}
SOLBASENAME=${17}
VISUALIZE=${18}
SOLUFILE=${19}   # - solu file, only necessary if $SETCUTOFF is 1

#args=("$@")
#for ((i=0; i < $#; i++)) {
#   echo "argument $((i+1)): ${args[$i]}"
#}

# new environment variables after running this script
# -None

#set solfile
SOLFILE=$CLIENTTMPDIR/${USER}-tmpdir/$SOLBASENAME.sol

# reset TMPFILE
echo > $TMPFILE
echo ""                              > $TMPFILE

# read in settings (even when using default, see bugzilla 600)
SETTINGS=$SCIPPATH/../settings/$SETNAME.prm
if test $SETNAME != "default"
then
    echo read $SETTINGS                  >> $TMPFILE
    echo disp settings changed           >> $TMPFILE
fi

# set non-default feasibility tolerance
if test $FEASTOL != "default"
then
    echo set simplex tolerances feas $FEASTOL    >> $TMPFILE
    echo set mip tolerances integrality $FEASTOL >> $TMPFILE
fi

# if permutation counter is positive add permutation seed (0 = default)
if test $p -gt 0
then
    echo "Warning: CPlex configuration currently cannot handle instance permutation"
    exit 1
fi

if test "$REOPT" = true
then
    # exit because reoptimization feature is not supported here
    echo "Warning: CPlex configuration currently cannot handle reoptimization"
    exit 1
fi

if test "$VISUALIZE" = true
then
    # exit because visualization feature is not supported here
    echo "Warning: CPlex configuration currently cannot handle visualization"
    exit 1
fi

# set objective limit: optimal solution value from solu file, if existent
if test $SETCUTOFF = 1 || test $SETCUTOFF = true
then
    # TODO setting cutoff requires knowledge about whether the objective sense is minimization or maximization
    echo "Warning: Setting a cutoff is currently not supported for Cplex configuration"
    exit 1
fi

echo set timelimit $TIMELIMIT           >> $TMPFILE
echo set clocktype 0                    >> $TMPFILE
echo set mip display 3                  >> $TMPFILE
echo set mip interval $DISPFREQ         >> $TMPFILE
echo set mip tolerances mipgap 0.0      >> $TMPFILE
echo set mip limits nodes $NODELIMIT    >> $TMPFILE
echo set mip limits treememory $MEMLIMIT >> $TMPFILE
echo set threads $THREADS               >> $TMPFILE
echo set parallel 1                     >> $TMPFILE
echo set lpmethod 4                     >> $TMPFILE
echo set barrier crossover -1           >> $TMPFILE
#echo write $SETFILE                     >> $TMPFILE
echo read $INSTANCE                  >> $TMPFILE
echo display problem stats              >> $TMPFILE
echo $OPTCOMMAND                        >> $TMPFILE
echo display solution quality           >> $TMPFILE
echo quit                               >> $TMPFILE

# currently, the solution checker only supports .mps-files.
# compare instance name (without .gz) to instance name stripped by .mps.
#if they are unequal, we have an mps-file
# TMPINSTANCE=`basename $INSTANCE .gz`
# TMPINSTANCEB=`basename $TMPINSTANCE .mps`
# if test "$TMPINSTANCEB" != "$TMPINSTANCE"
# then
   # TODO Solution checker implementation for CPLEX
   #echo write $SOLFILE             >> $TMPFILE
# fi
echo quit                              >> $TMPFILE
