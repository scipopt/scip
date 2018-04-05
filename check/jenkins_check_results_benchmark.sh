#! /bin/bash -x

#
# Usage:
# clusterbench.sh QUEUES EXECUTABLE | jenkins_check_results_benchmark.sh

# This script reads stdout from clusterbench.sh, parses the slurm job ids, and starts a
# job after the previously queued slurm jobs finish. This job waits for 5 seconds, then
# runs ./evalcheck_cluster_benchmark.sh and greps for fails.
# To know which results to process with evalcheck_cluster_benchmark
# The results are uploaded to rubberband with rbcli and if there are fails, an email is sent to the admin.
#

echo "This is jenkins_check_results_benchmark.sh running."

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

export OUTPUTDIR="results/clusterbench$(date +%Y%m%d)"

# build job ids string for sbatch dependency
jobidsstr=$(printf ",%s" "${slurmjobids[@]}")
jobidsstr=${jobidsstr:1}

# execute checker after all jobs completed
sbatch --dependency=afterany:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=4000 --time=500 --partition=mip-dbg --account=mip jenkins_failcheck_benchmark.sh
