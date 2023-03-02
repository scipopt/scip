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

# Start a local scip testrun.
#
# To be invoked by Makefile 'make test'.

TSTNAME="${1}"
BINNAME="${2}" #EXECUTABLE from makefile
SETNAMES="${3}"
BINID="${4}"
OUTPUTDIR="${5}"
TIMELIMIT="${6}"
NODELIMIT="${7}"
MEMLIMIT="${8}"
THREADS="${9}"
FEASTOL="${10}"
DISPFREQ="${11}"
CONTINUE="${12}"
LOCK="${13}"
VERSION="${14}"
LPS="${15}"
DEBUGTOOL="${16}"
CLIENTTMPDIR="${17}"
REOPT="${18}"
OPTCOMMAND="${19}"
SETCUTOFF="${20}"
MAXJOBS="${21}"
VISUALIZE="${22}"
PERMUTE="${23}"
SEEDS="${24}"
GLBSEEDSHIFT="${25}"
STARTPERM="${26}"
PYTHON="${27}"
EMPHBENCHMARK="${28}"

# check if all variables defined (by checking the last one)
if test -z "${PYTHON}"
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = ${TSTNAME}"
    echo "BINNAME       = ${BINNAME}"
    echo "SETNAMES      = ${SETNAMES}"
    echo "BINID         = ${BINID}"
    echo "OUTPUTDIR     = ${OUTPUTDIR}"
    echo "TIMELIMIT     = ${TIMELIMIT}"
    echo "NODELIMIT     = ${NODELIMIT}"
    echo "MEMLIMIT      = ${MEMLIMIT}"
    echo "THREADS       = ${THREADS}"
    echo "FEASTOL       = ${FEASTOL}"
    echo "DISPFREQ      = ${DISPFREQ}"
    echo "CONTINUE      = ${CONTINUE}"
    echo "LOCK          = ${LOCK}"
    echo "VERSION       = ${VERSION}"
    echo "LPS           = ${LPS}"
    echo "DEBUGTOOL     = ${DEBUGTOOL}"
    echo "CLIENTTMPDIR  = ${CLIENTTMPDIR}"
    echo "REOPT         = ${REOPT}"
    echo "OPTCOMMAND    = ${OPTCOMMAND}"
    echo "SETCUTOFF     = ${SETCUTOFF}"
    echo "MAXJOBS       = ${MAXJOBS}"
    echo "VISUALIZE     = ${VISUALIZE}"
    echo "PERMUTE       = ${PERMUTE}"
    echo "SEEDS         = ${SEEDS}"
    echo "GLBSEEDSHIFT  = ${GLBSEEDSHIFT}"
    echo "STARTPERM     = ${STARTPERM}"
    echo "PYTHON        = ${PYTHON}"
    echo "EMPHBENCHMARK = ${EMPHBENCHMARK}"
    exit 1;
fi

# call routines for creating the result directory, checking for existence
# of passed settings, etc
# defines the following environment variables: SCIPPATH, SETTINGSLIST, SOLUFILE, HARDMEMLIMIT, DEBUGTOOLCMD, INSTANCELIST,
#                                              TIMELIMLIST, HARDTIMELIMLIST
TIMEFORMAT="sec"
MEMFORMAT="kB"
. ./configuration_set.sh "${BINNAME}" "${TSTNAME}" "${SETNAMES}" "${TIMELIMIT}" "${TIMEFORMAT}" "${MEMLIMIT}" "${MEMFORMAT}" "${DEBUGTOOL}" "${SETCUTOFF}"

if test -e "${SCIPPATH}/../${BINNAME}"
then
    EXECNAME="${SCIPPATH}/../${BINNAME}"

    # check if we can set hard memory limit (address, leak, or thread sanitzer don't like ulimit -v)
    if [ $(uname) == Linux ] && (ldd "${EXECNAME}" | grep -q lib[alt]san); then
        # skip hard mem limit if using AddressSanitizer (libasan), LeakSanitizer (liblsan), or ThreadSanitizer (libtsan)
        HARDMEMLIMIT="none"
    elif [ $(uname) == Linux ] && (nm "${EXECNAME}" | grep -q __[alt]san); then
        # skip hard mem limit if using AddressSanitizer, LeakSanitizer, or ThreadSanitizer linked statitically (__[alt]san symbols)
        HARDMEMLIMIT="none"
    else
        ULIMITMEM="ulimit -v ${HARDMEMLIMIT} k;"
    fi
else
    EXECNAME="${BINNAME}"
    ULIMITMEM="ulimit -v ${HARDMEMLIMIT} k;"
fi

# use specified python version if the binary ends with ".py"
EXT="${BINNAME##*.}"
if test "${EXT}" = "py"
then
   EXECNAME="${PYTHON} ${BINNAME}"
fi

EXECNAME="${DEBUGTOOLCMD}${EXECNAME}"

INIT="true"
COUNT=0
for idx in "${!INSTANCELIST[@]}"
do
    # retrieve instance and timelimits from arrays set in the configuration_set.sh script
    INSTANCE=${INSTANCELIST[${idx}]}
    TIMELIMIT=${TIMELIMLIST[${idx}]}
    HARDTIMELIMIT=${HARDTIMELIMLIST[${idx}]}

    COUNT=$((COUNT + 1))

    # run different random seeds
    for ((s = 0; ${s} <= ${SEEDS}; s++))
    do

        # permute transformed problem
        for ((p = 0; ${p} <= ${PERMUTE}; p++))
        do

            # loop over settings
            for SETNAME in "${SETTINGSLIST[@]}"
            do
                # waiting while the number of jobs has reached the maximum
                if [ "${MAXJOBS}" -ne 1 ]
                then
                    while [ $(jobs -r | wc -l) -ge "${MAXJOBS}" ]
                    do
                        sleep 10
                        echo "Waiting for jobs to finish."
                    done
                fi

                # infer the names of all involved files from the arguments
                QUEUE=$(hostname)

                # infer the names of all involved files from the arguments
                # defines the following environment variables: OUTFILE, ERRFILE, EVALFILE, OBJECTIVEVAL, SHORTPROBNAME,
                #                                              FILENAME, SKIPINSTANCE, BASENAME, TMPFILE, SETFILE

                . ./configuration_logfiles.sh "${INIT}" "${COUNT}" "${INSTANCE}" "${BINID}" "${PERMUTE}" "${SEEDS}" "${SETNAME}" \
                    "${TSTNAME}" "${CONTINUE}" "${QUEUE}" "${p}" "${s}" "${THREADS}" "${GLBSEEDSHIFT}" "${STARTPERM}"

                if test "${INSTANCE}" = "DONE"
                then
                    wait
                    ./evalcheck_cluster.sh "${EVALFILE}"
                    continue
                fi

                if test "${SKIPINSTANCE}" = "true"
                then
                    continue
                fi

                # find out the solver that should be used
                SOLVER=$(stripversion "${BINNAME}")

                CONFFILE="configuration_tmpfile_setup_${SOLVER}.sh"

                # we don't have separate configuration files for most examples and applications, use SCIP configuration file instead
                if ! test -f "${CONFFILE}"
                then
                    CONFFILE="configuration_tmpfile_setup_scip.sh"
                fi

                # call tmp file configuration for SCIP
                # this may modify the EXECNAME environment variable
                . "./${CONFFILE}" "${INSTANCE}" "${SCIPPATH}" "${TMPFILE}" "${SETNAME}" "${SETFILE}" "${THREADS}" "${SETCUTOFF}" \
                    "${FEASTOL}" "${TIMELIMIT}" "${MEMLIMIT}" "${NODELIMIT}" "${LPS}" "${DISPFREQ}"  "${REOPT}" "${OPTCOMMAND}" \
                    "${CLIENTTMPDIR}" "${FILENAME}" "${VISUALIZE}" "${SOLUFILE}" "${EMPHBENCHMARK}"

                # additional environment variables needed by run.sh
                export SOLVERPATH="${SCIPPATH}"
                export BASENAME="${FILENAME}"
                export FILENAME="${INSTANCE}"
                export SOLNAME="${SOLCHECKFILE}"
                export TIMELIMIT
                export CLIENTTMPDIR
                export OUTPUTDIR
                export EXECNAME
                export CHECKERPATH="${SCIPPATH}/solchecker"

                WRITESETTINGS="false"
                if test "$INIT" = "true"
                then
                    if test "$SOLVER" = "scip"
                    then
                        WRITESETTINGS="true"
                    fi
                fi
                if test "$WRITESETTINGS" = "true"
                then
                    echo -e "#!/usr/bin/env bash \n $EXECNAME -s $SETTINGS -c 'set save $CHECKSETFILE quit'" > write-settings.sh
                fi

                if test "$WRITESETTINGS" = "true"
                then
                    if [ $MAXJOBS -eq 1 ]
                    then
                        bash write-settings.sh
                    else
                        bash write-settings.sh &
                    fi
                fi

                echo "Solving instance ${INSTANCE} with settings ${SETNAME}, hard time ${HARDTIMELIMIT}, hard mem ${HARDMEMLIMIT}"
                if [ "${MAXJOBS}" -eq 1 ]
                then
                    bash -c "ulimit -t ${HARDTIMELIMIT} s; ${ULIMITMEM} ulimit -f 200000; ./run.sh"
                else
                    bash -c "ulimit -t ${HARDTIMELIMIT} s; ${ULIMITMEM} ulimit -f 200000; ./run.sh" &
                fi

                if test "$WRITESETTINGS" = "true"
                then
                    rm write-settings.sh
                fi

            done # end for SETNAME
        done # end for PERMUTE
    done # end for SEEDS

    # after the first termination of the set loop, no file needs to be initialized anymore
    INIT="false"
done # end for TSTNAME
