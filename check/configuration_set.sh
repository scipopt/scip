#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# configures environment variables for test runs both on the cluster and locally
# to be invoked inside a check(_cluster.sh) script
# This script cancels the process if required variables are not correctly set

# new environment variables defined by this script:
#    SCIPPATH - absolute path to invocation working directory
#    SETTINGSLIST - array of setting file basenames. script will abort if any of them doesn't exist
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values
#    HARDMEMLIMIT - hard memory limit for the optimization call, that is given to slurm or the shell
#    DEBUGTOOLCMD - a debug tool command to use
#    INSTANCELIST - list of all instances with complete path
#    TIMELIMLIST - list of time limits for the individual instances
#    HARDTIMELIMLIST - list of hard time limits for the individual instances, that are given to slurm
#    FULLTSTNAME - path to .test file either in testset/ or in instancedata/testset

# function to capitalize a whole string
capitalize() {
    echo "${1}" | tr '[:lower:]' '[:upper:]'
}

# function to strip version of, e.g., scip-3.2... to only scip and scipampl.* to scipampl
stripversion() {
    NAMENOPATH=$(basename "${1}")
    # by '%%', Trim the longest match from the end
    NAMENOVERSION=${NAMENOPATH%%-*}
    NAMENOVERSION=${NAMENOVERSION%%\.*}
    echo "${NAMENOVERSION}"
}

# format time in seconds into format  dd-hh:mm:ss
formattime() {
    local T="${1}"
    local D=$((T / 60 / 60 / 24))
    local H=$((T / 60 / 60 % 24))
    local M=$((T / 60 % 60))
    local S=$((T % 60))
    local F="%02d"
    printf -v PRETTYTIME "${F}-${F}:${F}:${F}" "${D}" "${H}" "${M}" "${S}"
    echo "${PRETTYTIME}"
}

# input environment - these environment variables should be set before invoking this script
BINNAME="${1}"       # name of the binary
TSTNAME="${2}"       # name of the test set
SETNAMES="${3}"      # a comma-separated string of setting file basenames (without .set extension)
TIMELIMIT="${4}"     # the time limit in seconds
TIMEFORMAT="${5}"    # the format for the time (sec or format)
MEMLIMIT="${6}"      # the memory limit in MB
MEMFORMAT="${7}"     # the format for hard memory limit (kB or MB)
DEBUGTOOL="${8}"     # which debug tool should be used, if any?
SETCUTOFF="${9}"     # set this to 1 if you want the scripts to (try to) pass a best known primal bound (from .solu file) to the solver

# get current SCIP path
SCIPPATH=$(pwd)

# check if binary exists. The second condition checks whether there is a binary of that name directly available
# independent of whether it is a symlink, file in the working directory, or application in the path
if test ! -e "${SCIPPATH}/../${BINNAME}" && ! type "${BINNAME}" >/dev/null 2>&1
then
    echo "ERROR: \"${BINNAME}\" not found."
    echo "       This is needed by ${0} to work. Check your"
    echo "       \${PATH} variable or install the tool \"${BINNAME}\"."
    echo Skipping test since the binary ${BINNAME} does not exist.
    exit 1
fi

# create ${OUTPUTDIR} directory if it doesn't already exist
if test ! -e "${SCIPPATH}/${OUTPUTDIR}"
then
    mkdir "${SCIPPATH}/${OUTPUTDIR}"
fi

# create settings directory if non-existent
if test ! -d "${SCIPPATH}/../settings/"
then
    echo "Create directory settings"
    mkdir "${SCIPPATH}/../settings"
fi

# find out the solver that should be used
SOLVER=$(stripversion "${BINNAME}")

# figure out the correct settings file extension
if test "${SOLVER}" = cplex
then
    SETEXTEXTENSION="prm"
else
    SETEXTEXTENSION="set"
fi


# check if all settings files exist
SETTINGSLIST=(${SETNAMES//,/ })
for SETNAME in "${SETTINGSLIST[@]}"
do
    if test "${SOLVER}" = fscip
    then
        SETTINGS="${SCIPPATH}/../../ug/settings/${SETNAME}.${SETEXTEXTENSION}"
    else
        SETTINGS="${SCIPPATH}/../settings/${SETNAME}.${SETEXTEXTENSION}"
    fi
    if test "${SETNAME}" != "default" && test ! -e "${SETTINGS}"
    then
        echo "Skipping test since the settings file ${SETTINGS} does not exist."
        exit 1
    fi
done

# call method to obtain solution file
# defines the following environment variable: SOLUFILE
. ./configuration_solufile.sh "${TSTNAME}"

# if cutoff should be passed, solu file must exist
if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
then
    if test "${SOLUFILE}" = ""
    then
        echo "Skipping test: SETCUTOFF=1 set, but no .solu file (${TSTNAME}.solu or all.solu in testset/ or instancedata/testsets/) available"
        exit 1
    fi
fi

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))

# qsub requires the memory limit to be displayed in kB
if test "${MEMFORMAT}" = "kB"
then
    HARDMEMLIMIT=$((HARDMEMLIMIT * 1024))
elif test "${MEMFORMAT}" = "B"
then
    HARDMEMLIMIT=$((HARDMEMLIMIT * 1024000))
fi
# check if the test run should be processed in a debug tool environment
if test "${DEBUGTOOL}" = "valgrind"
then
    DEBUGTOOLCMD="valgrind --log-fd=1 --leak-check=full --suppressions=${SCIPPATH}/../suppressions.valgrind "
elif test "${DEBUGTOOL}" = "rr"
then
    DEBUGTOOLCMD="rr record --chaos -o RRTRACEFOLDER_PLACEHOLDER.rr "
elif test "${DEBUGTOOL}" = "gdb"
then
    #  set a gdb command, but leave a place holder for the error file we want to log to, which gets replaced in 'run.sh'
    DEBUGTOOLCMD='gdb -batch-silent -return-child-result -ex "run" -ex "set logging file ERRFILE_PLACEHOLDER" -ex "set logging on" -ex "thread apply all bt full" --args '
    if test "${SOLVER}" = "fscip"
    then
        DEBUGTOOLCMD='gdb -batch-silent -return-child-result -ex "run" -ex "set logging file ERRFILE_PLACEHOLDER" -ex "set logging on" -ex "set verbose on" -ex "thread apply all bt full" --args '
    fi
else
    DEBUGTOOLCMD=""
fi

#search for test file and check if we use a ttest or a test file
if [ -f "testset/${TSTNAME}.ttest" ]
then
    FULLTSTNAME="testset/${TSTNAME}.ttest"
    TIMEFACTOR=${TIMELIMIT}
elif [ -f "testset/${TSTNAME}.test" ]
then
    FULLTSTNAME="testset/${TSTNAME}.test"
    TIMEFACTOR=1
elif [ -f "instancedata/testsets/${TSTNAME}.ttest" ]
then
    FULLTSTNAME="instancedata/testsets/${TSTNAME}.ttest"
    TIMEFACTOR=${TIMELIMIT}
elif [ -f "instancedata/testsets/${TSTNAME}.test" ]
then
    FULLTSTNAME="instancedata/testsets/${TSTNAME}.test"
    TIMEFACTOR=1
else
    echo "Skipping test: no ${TSTNAME}.(t)test file found in testset/ or instancedata/testsets/"
    exit 1
fi

if [ "${TIMEFACTOR}" -gt 5 ]
then
    echo "ERROR: Time factor ${TIMEFACTOR} is too large--did you accidentally use a time limit in seconds in combination with a TTEST file?"
    echo "Exiting, try again, SCIP script novice"
    exit 1
fi

#write instance names to an array
COUNT=0
for INSTANCE in $(cat "${FULLTSTNAME}" | awk '{print $1}')
do
    # if the key word DONE appears in the test file, skip the remaining test file
    if test "${INSTANCE}" = "DONE"
    then
        break
    fi
    # check if problem instance exists
    if test -f "${SCIPPATH}/${INSTANCE}"
    then
        INSTANCELIST[${COUNT}]="${SCIPPATH}/${INSTANCE}"
        COUNT=$(( COUNT + 1 ))
    else
        echo "input file ${SCIPPATH}/${INSTANCE} not found!"
    fi
done
INSTANCELIST[${COUNT}]="DONE"

#write timelimits to an array
#if no second column with timelimits exists in the test file the normal timelimit will be returned by the awk command
COUNT=0
for TL in $(cat "${FULLTSTNAME}" | awk '{print match($2, /[^ ]/) ? $2 : "'${TIMELIMIT}'"}')
do
    TMPTL=$(( TL * TIMEFACTOR ))
    TIMELIMLIST[${COUNT}]=${TMPTL}
    # we add 100% to the hard time limit and additional 600 seconds in case of small time limits
    HARDTIMELIMIT=$(((TMPTL + 600) + TMPTL))

    if test "${TIMEFORMAT}" = "format"
    then
        #format is (d-)HH:MM:SS
        HARDTIMELIMIT=$(formattime "${HARDTIMELIMIT}")
        # echo "${HARDTIMELIMIT}"
    fi
    HARDTIMELIMLIST[${COUNT}]="${HARDTIMELIMIT}"
    COUNT=$(( COUNT + 1 ))
done
