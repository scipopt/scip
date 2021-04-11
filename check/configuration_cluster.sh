#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# configures environment variables for cluster runs.
# to be invoked inside a check_cluster.sh script
# This script cancels the process if required variables are not correctly set

# input variables - should be passed to this script
QUEUE="${1}"     # the name of the cluster (M620, dbg, telecom-dbg, mip-dbg, opt-low, opt)
PPN="${2}"       # number of cluster nodes to use
EXCLUSIVE="${3}" # should cluster nodes be blocked for other users while the jobs are running?
QUEUETYPE="${4}" # either 'srun' or 'qsub'

# new environment variables defined by this script:
NICE=""
ACCOUNT="mip"
CLUSTERQUEUE="${QUEUE}"

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

#define clusterqueue, which might not be the QUEUE, because this might be an alias for a bunch of QUEUEs
if test "${CLUSTERQUEUE}" = "dbg"
then
    CLUSTERQUEUE="mip-dbg,telecom-dbg"
    ACCOUNT="mip-dbg"
elif test "${CLUSTERQUEUE}" = "telecom-dbg"
then
    ACCOUNT="mip-dbg"
elif test "${CLUSTERQUEUE}" = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test "${CLUSTERQUEUE}" = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_opt
elif test "${CLUSTERQUEUE}" = "M620-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M620"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620
elif test "${CLUSTERQUEUE}" = "M620v3-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M620v3"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620v3
elif test "${CLUSTERQUEUE}" = "M630-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M630"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M630
elif test "${CLUSTERQUEUE}" = "M620x"
then
    CLUSTERQUEUE="M620,M620v2,M620v3"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620
    make --makefile=wakeup-slurm wake_M620v2
    make --makefile=wakeup-slurm wake_M620v3
elif test "${CLUSTERQUEUE}" = "M640-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M640"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M640
elif test "${CLUSTERQUEUE}" = "moskito"
then
    ACCOUNT="dopt"
fi

# check if the slurm blades should be used exclusively
if test "${EXCLUSIVE}" = "true"
then
    EXCLUSIVE=" --exclusive"
    if test "${CLUSTERQUEUE}" = "opt"
    then
        CLUSTERQUEUE="M640"
    fi
else
    EXCLUSIVE=""
fi
