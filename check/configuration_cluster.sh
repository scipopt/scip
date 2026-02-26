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

# configures environment variables for cluster runs.
# to be invoked inside a check_cluster.sh script
# This script cancels the process if required variables are not correctly set

# input variables - should be passed to this script
QUEUE="${1}"     # the name of the cluster queue (e.g., M640, Gold6338, moskito)
PPN="${2}"       # number of cluster nodes to use
EXCLUSIVE="${3}" # should cluster nodes be blocked for other users while the jobs are running?
QUEUETYPE="${4}" # either 'srun' or 'qsub'

# check if queue has been defined
if test "${QUEUE}" = ""
then
    echo "Skipping test since the queue ${QUEUE} name has not been defined."
    exit 1
fi

# check if number of nodes has been defined
if test "${PPN}" = ""
then
    echo "Skipping test since the number ${PPN} of nodes has not been defined."
    exit 1
fi

# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "${QUEUETYPE}" = "srun"
then
    HOSTMEM=$(ulimit -m)
    if test "${HOSTMEM}" != "unlimited"
    then
        if [ "$((HARDMEMLIMIT * 1024))" -gt "${HOSTMEM}" ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit 1
        fi
    fi
fi

CLUSTERQUEUE="${QUEUE}"
CONSTRAINT=""
NICE=""
ACCOUNT=""

if [ "${CLUSTERQUEUE}" = "moskito" ] || [ "${CLUSTERQUEUE}" = "prio" ]; then
    ACCOUNT="dopt"
elif [[ "$(uname -n)" =~ htc ]]; then
    # z1 cluster
    ACCOUNT="optimi_integer"
fi

if test "${CLUSTERQUEUE}" = "M640-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M640"
elif test "${CLUSTERQUEUE}" = "M640v2-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M640"
fi

if test "${CLUSTERQUEUE}" = "Gold6338"
then
    CONSTRAINT="Gold6338"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "Gold6342"
then
    CONSTRAINT="Gold6342"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "M640v2"
then
    CONSTRAINT="Gold5222"
    CLUSTERQUEUE="opt_int"
elif test "${CLUSTERQUEUE}" = "M640"
then
    CONSTRAINT="Gold5122"
    CLUSTERQUEUE="opt_int"
fi

# check if the slurm blades should be used exclusively
if test "${EXCLUSIVE}" = "true"
then
    EXCLUSIVE=" --exclusive"
else
    EXCLUSIVE=""
fi
