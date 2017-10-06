#! /bin/bash

#
# Usage:
# make testcluster | VERSION=scipbinversion PERF=performance check/jenkins_check_results.sh TESTSET SETTING

# This script reads stdout from make testcluster, parses the slurm job ids, and starts a
# job after the previously queued slurm jobs finish. This job waits for 5 seconds, then
# runs ./evalcheck_cluster.sh and greps for fails.
# To know which results to process with evalcheck_cluster, TESTSET and SETTING must be provided explicitly.
# optional: VERSION specifies the scip binary which should be used, PERF=performance enables rubberband support in jenkins_failcheck.sh.
# The results are uploaded to rubberband with rbcli and if there are fails, an email is sent to the admin.
#

# read from stdin
i=0
while read line
do
  if [[ $line == Submitted[[:space:]]batch[[:space:]]job[[:space:]]* ]]
  then
    stringarray=($line)
    slurmjobids[$i]=${stringarray[-1]}
    ((i++))
  fi
done < /dev/stdin

# build job ids string for sbatch dependency
jobidsstr=$(printf ",%s" "${slurmjobids[@]}")
jobidsstr=${jobidsstr:1}

export TESTSET=$1
export SETTING=$2

# process optional variables
if [ "${PERF}" != "" ]; then
  export PERFORMANCE=${PERF}
fi
# get some relevant information
if [ "${VERSION}" == "" ]; then
  export SCIPVERSION="scip"
  SCIPVERSIONOUTPUT=`bin/scip -v | sed -e 's/$/@/'`
else
  export SCIPVERSION="scip-${VERSION}"
  SCIPVERSIONOUTPUT=`bin/${SCIPVERSION}* -v | sed -e 's/$/@/'`
fi

# execute checker after all jobs completed
export GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`
export GITBRANCH=`git ls-remote --heads origin | grep $(git rev-parse HEAD)| cut -d / -f 3`
export OPT=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* OPT=\([^@]*\).*/\1/'`
export LPS=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* LPS=\([^@]*\).*/\1/'`

sbatch --dependency=afterany:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=100 --time=10 --partition=mip-dbg --account=mip check/jenkins_failcheck.sh
