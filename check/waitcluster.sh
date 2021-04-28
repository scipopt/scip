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

# In order to not overload the cluster, no jobs are submitted if the queue is too full
# instead, this script waits until the queue load falls under a threshold and returns
# for the calling script to continue submitting jobs.

# The cluster queue QUEUE has an upper bound of LIMIT jobs; if this limit is
# reached the submitted jobs are dumped; to avoid that we check the total
# load of the cluster and wait until it is safe (total load not more than
# LIMIT jobs or we have submitted less than LIMITQUEUE jobs) to submit the next job.

# Usage example: ./waitcluster.sh 1600 "${QUEUE}" 200

LIMIT="${1}"
QUEUE="${2}"
LIMITQUEUE="${3}"

while true
do
    ALLQUEUED=$(qstat | grep -c "")
    QUEUED=$(qstat -u "${USER}" | grep -c " ${QUEUE} ")
    RUNNING=$(qstat -u "${USER}" | grep -c " R ")

    # display current user load and total load
    echo "jobs in progress: ${RUNNING} / ${QUEUED} (${ALLQUEUED})"

    if test "${ALLQUEUED}" -le 1990
    then

        if test "${QUEUED}" -le "${LIMITQUEUE}"
        then
            break
        fi

        if test "${ALLQUEUED}" -le "${LIMIT}"
        then
            break
        fi
    fi

    sleep 30
done
