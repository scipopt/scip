#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# configures environment variables for test runs both on the cluster and locally
# to be invoked inside a check(_cluster.sh) script
# This script cancels the process if required variables are not correctly set

# new environment variables defined by this script:
#    SCIPPATH - absolute path to invocation working directory
#    SETTINGSLIST - array of setting file basenames. script will abort if any of them doesn't exist
#    SETCUTOFF - should optimal solution value (from solu file) be passed as objective limit?
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values
#    VALGRINDCMD - the valgrind command to use

# input environment - these environment variables should be set before invoking this script
BINNAME=$1       # name of the binary
TSTNAME=$2       # name of the test set
SETNAMES=$3      # a comma-separated string of setting file basenames (without .set extension)
TIMELIMIT=$4     # the time limit in seconds
TIMEFORMAT=$5    # the format for the time (sec or format)
MEMLIMIT=$6      # the memory limit in MB
MEMFORMAT=$7     # the format for hard memory limit (kB or MB)
VALGRIND=$8      # should valgrind be used?
SETCUTOFF=$9     # set this to 1 if you want the scripts to (try to) pass a best known primal bound (from .solu file) to the solver

# get current SCIP path
SCIPPATH=`pwd`

# check if binary exists
if test ! -e $SCIPPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# create results directory if it doesn't already exist
if test ! -e $SCIPPATH/results
then
    mkdir $SCIPPATH/results
fi

# create settings directory if non-existent
if test ! -d $SCIPPATH/../settings/
then
    echo Create directory settings
    mkdir $SCIPPATH/../settings
fi

# check if all settings files exist
SETTINGSLIST=(${SETNAMES//,/ })
for SETNAME in ${SETTINGSLIST[@]}
do
    SETTINGS=$SCIPPATH/../settings/$SETNAME.set
    if test $SETNAME != "default" && test ! -e $SETTINGS
    then
        echo Skipping test since the settings file $SETTINGS does not exist.
        exit
    fi
done

# if cutoff should be passed, check for solu file
if test $SETCUTOFF = 1
then
    SOLUFILE=""
    for SOLU in testset/$TSTNAME.solu testset/all.solu
    do
        if test -e $SOLU
        then
            SOLUFILE=$SOLU
            break
        fi
    done
    if test $SOLUFILE = ""
    then
        echo "Skipping test: SETCUTOFF=1 set, but no .solu file (testset/$TSTNAME.solu or testset/all.solu) available"
        exit
    fi
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

if test $TIMEFORMAT = "format"
then
    #format is (d-)HH:MM:SS
    TMP=`expr $HARDTIMELIMIT`
    HARDTIMELIMIT=""
    DIVISORS=(60 60 24)
    for((i=0; i<=2; i++))
    do
        printf -v HARDTIMELIMIT "%02d${HARDTIMELIMIT}" `expr ${TMP} % ${DIVISORS[i]}`
        # separate the numbers by colons except for the last (HH hours)
        if test $i -lt 2
        then
            HARDTIMELIMIT=":${HARDTIMELIMIT}"
        fi
        TMP=`expr ${TMP} / ${DIVISORS[i]}`
    done
    if test ${TMP} -gt 0
    then
        HARDTIMELIMIT=${TMP}-${HARDTIMELIMIT}
    fi
fi

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

# qsub requires the memory limit to be displayed in kB
if test "$MEMFORMAT" = "kB"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`
elif test "$MEMFORMAT" = "B"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
fi
# check if the test run should be processed in the valgrind environment
if test "$VALGRIND" = "true"
then
    VALGRINDCMD="valgrind --log-fd=1 --leak-check=full "
else
    VALGRINDCMD=""
fi

#check if additional instance paths are given
POSSIBLEPATHS=$SCIPPATH
if test -e paths.txt
then
    POSSIBLEPATHS="${POSSIBLEPATHS} `cat paths.txt`"
fi
POSSIBLEPATHS="${POSSIBLEPATHS} / DONE"
echo $POSSIBLEPATHS
