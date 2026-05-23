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

# launches multiple solver instances in parallel on one exclusive node.
# used by check_cluster.sh when EXCLUSIVE=auto to pack instances onto one node.
# core assignment by task rank via single srun -n N with --distribution=cyclic:cyclic.
#
# usage: run_group.sh SRUN_PREFIX NUM_INSTANCES [FIELDS...]
#
# FIELDS contains NUM_INSTANCES groups of 14 positional arguments each:
#   SOLVERPATH BASENAME FILENAME CLIENTTMPDIR OUTPUTDIR HARDTIMELIMIT HARDMEMLIMIT
#   CHECKERPATH SETFILE TIMELIMIT EXECNAME VIPRCHECKNAME VIPRCOMPNAME VIPRCOMPRESSNAME

SRUN_PREFIX="$1"
NUM_INSTANCES="$2"
shift 2
# encode per-instance parameters as indexed environment variables
for ((i = 0; i < NUM_INSTANCES; i++))
do
    printf -v "INST_${i}" 'export SOLVERPATH=%q BASENAME=%q FILENAME=%q CLIENTTMPDIR=%q OUTPUTDIR=%q HARDTIMELIMIT=%q HARDMEMLIMIT=%q CHECKERPATH=%q SETFILE=%q TIMELIMIT=%q EXECNAME=%q VIPRCHECKNAME=%q VIPRCOMPNAME=%q VIPRCOMPRESSNAME=%q' \
        "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}"
    export "INST_${i}"
    shift 14
done

# hold cores until all steps submitted
TICK=0.03125
HOLD=$(awk "BEGIN {printf \"%.5f\", ${NUM_INSTANCES} * ${TICK}}")

# delay launches so that srun steps bind stable cores
for ((i = 0; i < NUM_INSTANCES; i++))
do
    sleep ${TICK}
    (n=INST_${i}; eval "${!n}"; EXECNAME="bash -c 'sleep ${HOLD}; ${EXECNAME}'"; export SRUN="${SRUN_PREFIX}"; ./run.sh) &
done
wait
