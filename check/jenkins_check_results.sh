#! /bin/bash -x

#
# Usage:
# make testcluster | PERMUTE=permutations VERSION=scipbinversion PERF=performance check/jenkins_check_results.sh TESTSET SETTING

# This script is supposed to be used when you compiled your scip with make.
# This script reads stdout from make testcluster, parses the slurm job ids, and queues jenkins_failcheck.sh
# to run after the make testcluster jobs finish. The jenkins_failcheck script waits for 5 seconds, then
# runs ./evalcheck_cluster.sh and greps for fails, among other things
# To know which results to process with evalcheck_cluster, TESTSET and SETTING must be provided explicitly.
# optional: VERSION specifies the scip binary which should be used, PERF=performance enables rubberband support in jenkins_failcheck.sh.
# The results are uploaded to rubberband with rbcli and if there are fails, an email is sent to the admin.

echo "This is jenkins_check_results.sh running."

# set up environment for jenkins_failcheck.sh
export TESTSET=$1
export SETTING=$2

# get some relevant information
# process optional variables
if [ "${PERF}" != "" ]; then
  export PERFORMANCE=${PERF}
fi
if [ "${VERSION}" == "" ]; then
  # when version not given explicitly, recover it from executable
  SCIPVERSIONOUTPUT=`bin/scip -v | sed -e 's/$/@/'`
  export SCIPVERSION=scip-`echo $SCIPVERSIONOUTPUT | sed -e 's/.* VERSION=\([^@]*\).*/\1/'`
else
  export SCIPVERSION="scip-${VERSION}"
  SCIPVERSIONOUTPUT=`bin/${SCIPVERSION}* -v | sed -e 's/$/@/'`
fi

# if PERMUTE is not a number, set it to 0
re='^[0-9]+$'
if ! [[ $PERMUTE =~ $re ]] ; then
  PERMUTE="0"
fi
export PERMUTE

export GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`

# GIT_BRANCH is a jenkins variable, if not present, try to get it from the git repository. The second thing is not robust because there may be more branches that this HEAD is present in.
export GITBRANCH=`echo ${GIT_BRANCH} | cut -d / -f 2`
if [ "${GITBRANCH}" = "" ];
then
  export GITBRANCH=`git show -s --pretty=%D | cut -d , -f 2 | cut -d / -f 2 `
fi

export OPT=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* OPT=\([^@]*\).*/\1/'`
export LPS=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* LPS=\([^@]*\).*/\1/'`

# read from stdin
CANCEL_FILE=${TESTSET}_${SETTING}_${LPS}_cancel.txt
i=0
while read line
do
  if [[ $line == Submitted[[:space:]]batch[[:space:]]job[[:space:]]* ]]
  then
    stringarray=($line)
    slurmjobids[$i]=${stringarray[-1]}
    ((i++))
    echo "${stringarray[-1]}" >> "${CANCEL_FILE}"
  fi
done < /dev/stdin

echo "To cancel the jobs run"
echo 'for jobid in `cat '$CANCEL_FILE'`; do scancel $jobid; done'
echo "This is an experimental feature, use with caution. In particular, make sure no two jobs have the same TESTSET, SETTING and LPS combination!"

# build job ids string for sbatch dependency
jobidsstr=$(printf ",%s" "${slurmjobids[@]}")
jobidsstr=${jobidsstr:1}

# execute checker after all jobs completed
sbatch --dependency=afterany:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=4000 --time=500 --partition=mip-dbg --account=mip check/jenkins_failcheck.sh
