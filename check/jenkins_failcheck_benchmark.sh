#!/bin/bash -x
#
# Usage:
# clusterbench.sh QUEUES EXECUTABLE

# This script constructs the testrun files for each queuenode.
# note: on call we are in the check directory

echo "This is jenkins_failcheck_benchmark.sh running."

sleep 5

# evaluate the run
echo "Evaluating the benchmark run."

cd check

# SCIP check files are check.clusterbench.SCIPVERSION.otherstuff.SETTING.{out,err,res,meta}
EVALFILES=`ls ${OUTPUTDIR}/*.eval`

# evaluate with evalcheck_cluster.sh
for EVALFILE in ${EVALFILES};
do
  ./evalcheck_cluster.sh -r ${EVALFILE}
done

cd ..
