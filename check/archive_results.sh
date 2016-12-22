#! /bin/bash

#
# Usage:
# make testcluster | check/archive_results.sh
#
# This script reads stdout from make testcluster, parses the slurm job ids, and starts a
# job after the previously queued slurm jobs finish. This job runs ./evalcheck_cluster.sh 
# and rbcli to send the results to rubberband.
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
done < "${1:-/dev/stdin}"

# build job ids string for sbatch dependency
jobidsstr=$(printf ",%s" "${slurmjobids[@]}")
jobidsstr=${jobidsstr:1}

# execute rbcliwrapper if all jobs succeed
sbatch --dependency=afterok:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=100 --time=2 --partition=mip-dbg --account=mip check/rbcliwrapper.sh

sbatch --dependency=afternotok:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=100 --time=2 --partition=mip-dbg --account=mip --mail-type=BEGIN --mail-user=timo-admin@zib.de --output=/dev/null --job-name=failed-run <<< "#! /bin/bash"
