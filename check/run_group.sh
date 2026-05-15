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
#
# usage: run_group.sh SRUN_PREFIX NUM_INSTANCES THREADS_SAFE [FIELDS...]
#
# THREADS_SAFE is the number of physical cores per instance.
# Node topology (sockets, cores, NUMA) is queried from the compute node at runtime.
#
# FIELDS contains NUM_INSTANCES groups of 14 positional arguments each:
#   SOLVERPATH BASENAME FILENAME CLIENTTMPDIR OUTPUTDIR HARDTIMELIMIT HARDMEMLIMIT
#   CHECKERPATH SETFILE TIMELIMIT EXECNAME VIPRCHECKNAME VIPRCOMPNAME VIPRCOMPRESSNAME

SRUN_PREFIX="$1"
NUM_INSTANCES="$2"
THREADS_SAFE="$3"
shift 3

# query node topology
NSOCKETS=$(lscpu | awk '/^Socket/ {print $NF}')
CORES_PER_SOCKET=$(lscpu | awk '/Core.s. per socket/ {print $NF}')

# collect physical core CPU indices per socket (first thread of each core)
for ((s = 0; s < NSOCKETS; s++))
do
    SOCKET_CORES[$s]=$(lscpu --parse=CPU,Core,Socket | grep -v '^#' | awk -F, -v sock="${s}" '$3 == sock && !seen[$2]++ {printf "%s ", $1}')
done
# instances per socket and stride
IPS=$(( (NUM_INSTANCES + NSOCKETS - 1) / NSOCKETS ))
STRIDE=$(( CORES_PER_SOCKET / IPS ))

# assign instances sequentially across sockets with common stride
instance=0
for ((s = 0; s < NSOCKETS && instance < NUM_INSTANCES; s++))
do
    for ((i = 0; i < IPS && instance < NUM_INSTANCES; i++))
    do
        cores=(${SOCKET_CORES[$s]})
        CPU_MAP=""
        for ((t = 0; t < THREADS_SAFE; t++))
        do
            CPU_MAP="${CPU_MAP:+${CPU_MAP},}${cores[$((i * STRIDE + t))]}"
        done
        (
        export SOLVERPATH="$1" BASENAME="$2" FILENAME="$3" CLIENTTMPDIR="$4" OUTPUTDIR="$5" \
            HARDTIMELIMIT="$6" HARDMEMLIMIT="$7" CHECKERPATH="$8" SETFILE="$9" TIMELIMIT="${10}" \
            EXECNAME="${11}" VIPRCHECKNAME="${12}" VIPRCOMPNAME="${13}" VIPRCOMPRESSNAME="${14}"
        export SRUN="${SRUN_PREFIX} --cpu_bind=verbose,map_cpu:${CPU_MAP}"
        ./run.sh
        ) &
        shift 14
        instance=$((instance + 1))
    done
done
wait
