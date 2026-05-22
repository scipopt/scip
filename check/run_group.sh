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
# core assignment is handled by Slurm via --distribution=cyclic:cyclic on sbatch.
#
# usage: run_group.sh SRUN_PREFIX NUM_INSTANCES [FIELDS...]
#
# FIELDS contains NUM_INSTANCES groups of 14 positional arguments each:
#   SOLVERPATH BASENAME FILENAME CLIENTTMPDIR OUTPUTDIR HARDTIMELIMIT HARDMEMLIMIT
#   CHECKERPATH SETFILE TIMELIMIT EXECNAME VIPRCHECKNAME VIPRCOMPNAME VIPRCOMPRESSNAME

SRUN_PREFIX="$1"
NUM_INSTANCES="$2"
shift 2

for ((i = 0; i < NUM_INSTANCES; i++))
do
    # delay launches so that srun steps bind stable cores
    sleep 0.0078125
    (
    export SOLVERPATH="$1" BASENAME="$2" FILENAME="$3" CLIENTTMPDIR="$4" OUTPUTDIR="$5" \
        HARDTIMELIMIT="$6" HARDMEMLIMIT="$7" CHECKERPATH="$8" SETFILE="$9" TIMELIMIT="${10}" \
        EXECNAME="${11}" VIPRCHECKNAME="${12}" VIPRCOMPNAME="${13}" VIPRCOMPRESSNAME="${14}"
    export SRUN="${SRUN_PREFIX} --cpu_bind=verbose,cores"
    ./run.sh
    ) &
    shift 14
done
wait
