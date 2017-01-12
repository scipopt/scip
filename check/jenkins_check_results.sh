#! /bin/bash

#
# Usage:
# make testcluster | check/jenkins_check_results.sh TESTSET

# This script reads stdout from make testcluster, parses the slurm job ids, and starts a
# job after the previously queued slurm jobs finish. This job waits for 5 seconds, then 
# runs ./evalcheck_cluster.sh and greps for fails. 
# To know which results to process with evalcheck_cluster, TESTSET must be provided explicitly.
# If there is a fail, an email is sent to the admin. Otherwise, the results are uploaded
# to rubberband with rbcli.
# 
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

# execute checker after all jobs completed
export TESTSET=$1
sbatch --dependency=afterany:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=100 --time=10 --partition=mip-dbg --account=mip check/jenkins_failcheck.sh
