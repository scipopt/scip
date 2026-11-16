#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      *
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

### sets EXECNAME to run lj with appropriate arguments
### creates instance/settings-specific settings file

# environment variables passed as arguments
INSTANCE="${1}"        # - instance name to solve
SCIPPATH="${2}"        # - path to working directory for test (usually, the check subdirectory)
TMPFILE="${3}"         # - the batch file to control lj - won't be created
SETNAME="${4}"         # - specified basename of settings-file, or 'default'
SETFILE="${5}"         # - instance/settings specific settings-file
THREADS="${6}"         # - the number of LP solver threads to use
SETCUTOFF="${7}"       # - should optimal instance value be used as objective limit (0 or 1)?
FEASTOL="${8}"         # - feasibility tolerance, or 'default'
TIMELIMIT="${9}"       # - time limit for the solver
MEMLIMIT="${10}"       # - memory limit for the solver
NODELIMIT="${11}"      # - node limit for the solver
LPS="${12}"            # - LP solver to use
DISPFREQ="${13}"       # - display frequency for chronological output table
REOPT="${14}"          # - true if we use reoptimization, i.e., using a difflist file instead of an instance file
OPTCOMMAND="${15}"     # - command that should be executed after reading the instance, e.g. optimize, presolve or count
CLIENTTMPDIR="${16}"   # - directory for temporary files
SOLBASENAME="${17}"    # - base name for solution file
VISUALIZE="${18}"      # - true, if the branch-and-bound search should be visualized
SOLUFILE="${19}"       # - solu file, only necessary if ${SETCUTOFF} is 1
EMPHBENCHMARK="${20}"  # - use set emphasis benchmark
CLOCKTYPE="${21}"      # - clocktype (1 = CPU, 2 = wallclock)
WITHCERTIFICATE="${22}" # - true, if a certificate file should be created
KEEPSOL="${23}"        # - true, if file with solution values should be written and kept

#args=("$@")
#for ((i=0; i < $#; i++)) {
#   echo "argument $((i+1)): ${args[${i}]}"
#}

# new environment variables after running this script
# -None

# updated environment variables after running this script
# EXECNAME

if [ "${REOPT}" = true ] ; then
  echo "REOPT=true not supported for Lennard-Jones Cluster example" > /dev/stderr
  exit 1
fi

if [ "${OPTCOMMAND}" != "optimize" ] ; then
  echo "OPTCOMMAND=${OPTCOMMAND} not supported for Lennard-Jones Cluster example" > /dev/stderr
  exit 1
fi

if [ "${WITHCERTIFICATE}" = true ] ; then
  echo "WITHCERTIFICATE=true not supported for Lennard-Jones Cluster example" > /dev/stderr
  exit 1
fi

if [ "${KEEPSOL}" = true ] ; then
  echo "KEEPSOL=true not supported for Lennard-Jones Cluster example; ignored"
fi

# create setting file
echo -n > "${SETFILE}"

# emphasis needs to be before other option settings
[ "${EMPHBENCHMARK}" = true ] && echo "emphasis: benchmark" >> "${SETFILE}"

ORIGSETFILE="${SCIPPATH}/../settings/${SETNAME}.set"
if [ -e "$ORIGSETFILE" ] ; then
  cat "${ORIGSETFILE}" >> "${SETFILE}"
elif [ "$SETNAME" != default ] ; then
  echo "Error: Parameter file $ORIGSETFILE not found."
  exit 1
fi

# set non-default feasibility tolerance
[ "${FEASTOL}" != "default" ] && echo "limits/feastol = ${FEASTOL}" >> "${SETFILE}"

# if permutation counter is positive, add permutation seed (0 = default)
PERM=$((p + STARTPERM))
[ "${PERM}" -gt 0 ] && echo "randomization/permutationseed = ${PERM}" >> "${SETFILE}"

# set random seed shift if nonzero
SEED=$((s + GLBSEEDSHIFT))
[ "${SEED}" -gt 0 ] && echo "randomization/randomseedshift = ${SEED}" >> "${SETFILE}"

# avoid solving LPs in case of LPS=none
[ "${LPS}" = "none" ] && echo "lp/solvefreq = -1" >> "${SETFILE}"

# set reference value
[ "${OBJECTIVEVAL}" != "" ] && echo "misc/referencevalue = ${OBJECTIVEVAL}" >> "${SETFILE}"

echo "limits/time = ${TIMELIMIT}" >> "${SETFILE}"
echo "limits/nodes = ${NODELIMIT}" >> "${SETFILE}"
echo "limits/memory = ${MEMLIMIT}" >> "${SETFILE}"
echo "lp/threads = ${THREADS}" >> "${SETFILE}"
echo "timing/clocktype = ${CLOCKTYPE}" >> "${SETFILE}"
echo "display/freq = ${DISPFREQ}" >> "${SETFILE}"

if [ "${VISUALIZE}" = true ] ; then
  BAKFILENAME="$(basename ${TMPFILE} .tmp).dat"
  echo "visualization output set to ${BAKFILENAME}"
  echo "visual/bakfilename = ${OUTPUTDIR}/${BAKFILENAME}" >> ${SETFILE}
fi

# set objective limit: optimal solution value from solu file, if existent
[ "( ${SETCUTOFF} = 1 -o ${SETCUTOFF} = true ) -a -n \"${OBJECTIVEVAL}\"" ] && echo "limits/objective = ${OBJECTIVEVAL}" >> "${SETFILE}"

# parse instance: the first word in the file gives the number of particles
NPARTICLES=`awk '{ print $1 ; exit; }' ${INSTANCE}`

# define how to call instance
EXECNAME="${EXECNAME} ${NPARTICLES} ${SETFILE}"
