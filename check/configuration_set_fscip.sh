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

# configures environment variables for test runs both on the cluster and locally
# to be invoked inside a check(_cluster.sh) script
# This script cancels the process if required variables are not correctly set

# new environment variables defined by this script:
#    SCIPPATH - absolute path to invocation working directory
#    SETTINGSLIST - array of setting file basenames. script will abort if any of them doesn't exist
#    SETCUTOFF - should optimal solution value (from solu file) be passed as objective limit?
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values
#    DEBUGTOOLCMD - a debug tool command to use
#    INSTANCELIST - list of all instances with complete path

# function to capitalize a whole string
function capitalize {
    echo "$1" | tr '[:lower:]' '[:upper:]'
}

# function to strip version of, e.g., scip-3.2... to only scip and scipampl.* to scipampl
function stripversion {
    NAMENOPATH=`basename $1`
    # by '%%', Trim the longest match from the end
    NAMENOVERSION=${NAMENOPATH%%-*}
    NAMENOVERSION=${NAMENOVERSION%%\.*}
    echo $NAMENOVERSION
}

# input environment - these environment variables should be set before invoking this script
BINNAME=$1       # name of the binary
TSTNAME=$2       # name of the test set
SETNAMES=$3      # a comma-separated string of setting file basenames (without .set extension)
TIMELIMIT=$4     # the time limit in seconds
TIMEFORMAT=$5    # the format for the time (sec or format)
MEMLIMIT=$6      # the memory limit in MB
MEMFORMAT=$7     # the format for hard memory limit (kB or MB)
DEBUGTOOL=$8      # which debug tool should be used, if any?
SETCUTOFF=$9     # set this to 1 if you want the scripts to (try to) pass a best known primal bound (from .solu file) to the solver

# get current SCIP path
SCIPPATH=`pwd`

# check if binary exists. The second condition checks whether there is a binary of that name directly available
# independent of whether it is a symlink, file in the working directory, or application in the path
if test ! -e $SCIPPATH/../bin/$BINNAME && ! type $BINNAME >/dev/null 2>&1
then
   echo "ERROR: \"$SCIPPATH/../bin/$BINNAME\" not found."
   echo "       This is needed by ${0} to work. Check your"
   echo "       \$PATH variable or install the tool \"$BINNAME\"."
   echo Skipping test since the binary $BINNAME does not exist.
   exit 2
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

# figure out the correct settings file extension
if test $BINNAME = cplex
then
    SETEXTEXTENSION="prm"
else
    SETEXTEXTENSION="set"
fi


# check if all settings files exist
SETTINGSLIST=(${SETNAMES//,/ })
for SETNAME in ${SETTINGSLIST[@]}
do
    SETTINGS="${SCIPPATH}/../../ug/settings/${SETNAME}.${SETEXTEXTENSION}"
    if test $SETNAME != "default" && test ! -e $SETTINGS
    then
        echo Skipping test since the settings file $SETTINGS does not exist.
        exit
    fi
done

SOLUFILE=""
for SOLU in testset/$TSTNAME.solu testset/all.solu
do
    if test -e $SOLU
    then
        SOLUFILE=$SOLU
        break
    fi
done

# if cutoff should be passed, solu file must exist
if test $SETCUTOFF = 1 || test $SETCUTOFF = true
then
    if test $SOLUFILE = ""
    then
        echo "Skipping test: SETCUTOFF=1 set, but no .solu file (testset/$TSTNAME.solu or testset/all.solu) available"
        exit
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
# check if the test run should be processed in a debug tool environment
if test "$DEBUGTOOL" = "valgrind"
then
    DEBUGTOOLCMD="valgrind --log-fd=1 --leak-check=full --suppressions=${SCIPPATH}/../suppressions.valgrind "
elif test "$DEBUGTOOL" = "gdb"
then
    #  set a gdb command, but leave a place holder for the error file we want to log to, which gets replaced in 'run.sh'
    DEBUGTOOLCMD='gdb -batch-silent -return-child-result -ex "run" -ex "set logging file ERRFILE_PLACEHOLDER" -ex "set logging on" -ex "thread apply all bt full" --args '
else
    DEBUGTOOLCMD=""
fi

#check if additional instance paths are given
POSSIBLEPATHS=$SCIPPATH
if test -e paths.txt
then
    POSSIBLEPATHS="${POSSIBLEPATHS} `cat paths.txt`"
fi
POSSIBLEPATHS="${POSSIBLEPATHS} / DONE"
# echo $POSSIBLEPATHS

#check if we use a ttest or a test file
if [ -f testset/$TSTNAME.ttest ];
then
    FULLTSTNAME="testset/$TSTNAME.ttest"
    TIMEFACTOR=$TIMELIMIT
else
    FULLTSTNAME="testset/$TSTNAME.test"
    TIMEFACTOR=1
fi

#write instance names to an array
COUNT=0
for INSTANCE in `cat $FULLTSTNAME | awk '{print $1}'`
do
    # if the key word DONE appears in the test file, skip the remaining test file
    if test "$INSTANCE" = "DONE"
    then
        break
    fi
    # check if problem instance exists
    for IPATH in ${POSSIBLEPATHS[@]}
    do
        #echo $IPATH
        if test "$IPATH" = "DONE"
        then
            echo "input file $INSTANCE not found!"
        elif test -f $IPATH/$INSTANCE
        then
            INSTANCELIST[$COUNT]="${IPATH}/${INSTANCE}"
            break
        fi
    done
    COUNT=$(( $COUNT + 1 ))
done
INSTANCELIST[$COUNT]="DONE"
COUNT=$(( $COUNT + 1 ))

#write timelimits to an array
#if no second column with timelimits exists in the test file the normal timelimit will be returned by the awk command
COUNT=0
for TL in `cat $FULLTSTNAME | awk '{print match($2, /[^ ]/) ? $2 : "'$TIMELIMIT'"}'`
do
    TMPTL=$(( $TL * $TIMEFACTOR ))
    TIMELIMLIST[$COUNT]=$TMPTL
    # we add 100% to the hard time limit and additional 600 seconds in case of small time limits
    HARDTIMELIMIT=`expr \`expr $TMPTL + 600\` + $TMPTL`

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
    HARDTIMELIMLIST[$COUNT]=$HARDTIMELIMIT
    COUNT=$(( $COUNT + 1 ))
done
