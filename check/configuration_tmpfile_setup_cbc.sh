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
MEMLIMIT="${10}"     # - memory limit for the solver
NODELIMIT="${11}"    # - node limit for the solver
LPS="${12}"          # - LP solver to use
DISPFREQ="${13}"     # - display frequency for chronological output table
REOPT="${14}"        # - true if we use reoptimization, i.e., using a difflist file instead if an instance file
OPTCOMMAND="${15}"   # - command that should per executed after reading the instance, e.g. optimize, presolve or count
CLIENTTMPDIR="${16}" # - directory for temporary files
SOLBASENAME="${17}"  # - base name for solution file
VISUALIZE="${18}"    # - true, if the branch-and-bound search should be visualized
SOLUFILE="${19}"     # - solu file, only necessary if "${SETCUTOFF}" is 1
EMPHBENCHMARK="${20}"  # - use set emphasis benchmark

# new environment variables after running this script
# -None

# set solfile
SOLFILE="${CLIENTTMPDIR}/${USER}-tmpdir/${SOLBASENAME}.sol"

if test "${p}" -gt 0
then
    echo "Warning: CBC configuration currently cannot handle instance permutation"
    exit 1
fi

if test "${SETNAME}" != "default"
then
    echo "Warning: CBC configuration currently cannot handle non-default settings"
    exit 1
fi

if test "${REOPT}" = true
then
    echo "Warning: CBC configuration currently cannot handle reoptimization"
    exit 1
fi

if test "${VISUALIZE}" = true
then
    echo "Warning: CBC configuration currently cannot handle visualization"
    exit 1
fi

if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
then
    echo "Warning: Setting a cutoff is currently not supported for CBC configuration"
    exit 1
fi

# The following variables are ignored:
# "${DISPFREQ}", "${OPTCOMMAND}"

# workaround: since CBC only looks at cpu-time, we multiply the timelimit with the number of threads
TIMELIMIT=$((TIMELIMIT * THREADS))
echo "seconds ${TIMELIMIT}"                 >  "${TMPFILE}"
# set mip gap:
echo "ratioGap 0.0"                         >> "${TMPFILE}"
if test "${FEASTOL}" != "default"
then
    echo "primalTolerance ${FEASTOL}"       >> "${TMPFILE}"
    echo "integerTolerance ${FEASTOL}"      >> "${TMPFILE}"
fi
echo "maxNodes ${NODELIMIT}"                >> "${TMPFILE}"
echo "threads ${THREADS}"                   >> "${TMPFILE}"

echo "import ${INSTANCE}"                   >> "${TMPFILE}"
echo "ratioGap"                             >> "${TMPFILE}"
echo "allowableGap"                         >> "${TMPFILE}"
echo "seconds"                              >> "${TMPFILE}"
echo "stat"                                 >> "${TMPFILE}"
echo "solve"                                >> "${TMPFILE}"
echo "quit"                                 >> "${TMPFILE}"
