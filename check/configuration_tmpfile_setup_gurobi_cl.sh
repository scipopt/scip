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

### resets and fills a batch file TMPFILE to run CBC with
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

# new environment variables after running this script
# -None

# updated environment variables after running this script
# EXECNAME

# init an empty list of command line settings to be used by gurobi_cl
CLSETTINGSLIST=""

#append feasibility tolerance in case of non-default values
if test "${FEASTOL}" != "default"
then
    CLSETTINGSLIST="${CLSETTINGSLIST} FeasibilityTol=${FEASTOL} IntFeasTol=${FEASTOL}"
fi
CLSETTINGSLIST="${CLSETTINGSLIST} TimeLimit=${TIMELIMIT} NodeLimit=${NODELIMIT} DisplayInterval=${DISPFREQ} MIPGap=0.0 Threads=${THREADS} Crossover=0 Method=2"

# parse settings from settings file via awk
if test "${SETNAME}" != "default"
then
    echo $(pwd)
    CLSETTINGSLIST="$(awk 'BEGIN { finalstr=""} {finalstr=finalstr " "$1"="$2} END {print finalstr}' ${SETTINGS}) ${CLSETTINGSLIST}"
fi

#have a look if Gurobi is invoked with the settings you are asking for
echo "Gurobi will be invoked with arguments: ${CLSETTINGSLIST}"

EXECNAME="${EXECNAME} ${CLSETTINGSLIST} ${INSTANCE}"
