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

### resets and fills a batch file TMPFILE to run SCIP with
### sets correct limits, reads in settings, and controls
### display of the solving process

# environment variables passed as arguments
INSTANCE="${1}"        #  instance name to solve
SCIPPATH="${2}"        # - path to working directory for test (usually, the check subdirectory)
TMPFILE="${3}"         # - the batch file to control XPRESS
SETNAME="${4}"         # - specified basename of settings-file, or 'default'
SETFILE="${5}"         # - instance/settings specific set-file
THREADS="${6}"         # - the number of LP solver threads to use
SETCUTOFF="${7}"       # - should optimal instance value be used as objective limit (0 or 1)?
FEASTOL="${8}"         # - feasibility tolerance, or 'default'
TIMELIMIT="${9}"       # - time limit for the solver
MEMLIMIT="${10}"       # - memory limit for the solver
NODELIMIT="${11}"      # - node limit for the solver
LPS="${12}"            # - LP solver to use
DISPFREQ="${13}"       # - display frequency for chronological output table
REOPT="${14}"          # - true if we use reoptimization, i.e., using a difflist file instead if an instance file
OPTCOMMAND="${15}"     # - command that should per executed after reading the instance, e.g. optimize, presolve or count
CLIENTTMPDIR="${16}"   # - directory for temporary files
SOLBASENAME="${17}"    # - base name for solution file
VISUALIZE="${18}"      # - true, if the branch-and-bound search should be visualized
SOLUFILE="${19}"       # - solu file, only necessary if ${SETCUTOFF} is 1
EMPHBENCHMARK="${20}"  # - use set emphasis benchmark

#args=("$@")
#for ((i=0; i < $#; i++)) {
#   echo "argument $((i+1)): ${args[${i}]}"
#}

# new environment variables after running this script
# -None

#set solfile
SOLFILE="${CLIENTTMPDIR}/${USER}-tmpdir/${SOLBASENAME}.sol"

# reset TMPFILE
echo > "${TMPFILE}"

# read in settings (even when using default, see bugzilla 600)
SETTINGS="${SCIPPATH}/../settings/${SETNAME}.set"
if test "${SETNAME}" == "default"
then
    # create empty settings file
    test -e "${SETTINGS}" || touch "${SETTINGS}"
fi
echo set load "${SETTINGS}"            >>  "${TMPFILE}"

# set non-default feasibility tolerance
if test "${FEASTOL}" != "default"
then
    echo "set numerics feastol ${FEASTOL}" >> "${TMPFILE}"
fi

# if permutation counter is positive add permutation seed (0 = default)
PERM=$((p + STARTPERM))
if test "${PERM}" -gt 0
then
    echo "set randomization permutationseed ${PERM}"   >> "${TMPFILE}"
fi

# set random seed shift
SEED=$((s + GLBSEEDSHIFT))
if test "${SEED}" -gt 0
then
    echo "set randomization randomseedshift ${SEED}" >> "${TMPFILE}"
fi

# avoid solving LPs in case of LPS=none
if test "${LPS}" = "none"
then
    echo "set lp solvefreq -1"           >> "${TMPFILE}"
fi

# set reference value
if test "${OBJECTIVEVAL}" != ""
then
    #echo "Reference value ${OBJECTIVEVAL}"
    echo "set misc referencevalue ${OBJECTIVEVAL}"      >> "${TMPFILE}"
fi

INSTANCENAME=${INSTANCE%%.gz}
for i in gz mps cip lp
do
    INSTANCENAME=${INSTANCENAME%%.${i}}
done
INSTANCESETTINGSFILE="${INSTANCENAME}.set"

if test -f "${INSTANCESETTINGSFILE}"
then
    echo set load "${INSTANCESETTINGSFILE}"                      >> "${TMPFILE}"
fi

if  [ "${EMPHBENCHMARK}" = true ] ; then
    echo "set emphasis benchmark"                                >> "${TMPFILE}" # avoid switching to dfs etc. - better abort with memory error; this has to be first
fi
echo "set limits time ${TIMELIMIT}"                              >> "${TMPFILE}"
echo "set limits nodes ${NODELIMIT}"                             >> "${TMPFILE}"
echo "set limits memory ${MEMLIMIT}"                             >> "${TMPFILE}"
echo "set lp advanced threads ${THREADS}"                        >> "${TMPFILE}"
echo "set timing clocktype 1"                                    >> "${TMPFILE}"
echo "set display freq ${DISPFREQ}"                              >> "${TMPFILE}"
echo "set save ${SETFILE}"                                       >> "${TMPFILE}"

if test "${VISUALIZE}" = true
then
    BAKFILENAME="$(basename ${TMPFILE} .tmp).dat"
    echo "visualization output set to ${BAKFILENAME}"
    echo "set visual bakfilename ${OUTPUTDIR}/${BAKFILENAME}"    >> ${TMPFILE}
fi

if test "${REOPT}" = false
then
    # read and solve the instance
    echo "read ${INSTANCE}"                                      >> "${TMPFILE}"
    INSTANCENAME=${INSTANCE%%.gz}
    # if a decomposition in gzipped format (.dec.gz) with the basename of the instance lies in the same directory,
    # read it into SCIP, as well
    DECOMP="${INSTANCENAME}.dec.gz"
    if test -f "${DECOMP}"
    then
        echo "read ${DECOMP}"                                    >> "${TMPFILE}"
    fi
    # set objective limit: optimal solution value from solu file, if existent
    if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
    then
        if test "${OBJECTIVEVAL}" != ""
        then
            echo "set limits objective ${OBJECTIVEVAL}"          >> "${TMPFILE}"
        fi
        echo "set heur emph off"                                 >> "${TMPFILE}"
    fi

    echo "display parameters"                                    >> "${TMPFILE}"
    echo "${OPTCOMMAND}"                                         >> "${TMPFILE}"
    echo "display statistics"                                    >> "${TMPFILE}"
    echo "checksol"                                              >> "${TMPFILE}"
else
    # read the difflist file
    cat "${INSTANCE}"                                            >> "${TMPFILE}"
fi

# currently, the solution checker only supports .mps-files.
# compare instance name (without .gz) to instance name stripped by .mps.
#if they are unequal, we have an mps-file
TMPINSTANCE=$(basename "${INSTANCE}" .gz)
TMPINSTANCEB=$(basename "${TMPINSTANCE}" .mps)
if test "${TMPINSTANCEB}" != "${TMPINSTANCE}"
then
    echo "write sol ${SOLFILE}"                                  >> "${TMPFILE}"
fi
echo "quit"                                                      >> "${TMPFILE}"
