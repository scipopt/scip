#! /bin/bash -x

# Usage:
# make testcluster | TESTSET=testset SETTING=setting PERMUTE=permutations EXECUTABLE=build/bin/scip PERF=performance check/jenkins_check_results.sh
# make testcluster | TESTSET=testset SETTING=setting PERMUTE=permutations VERSION=scipbinversion PERF=performance check/jenkins_check_results.sh
# or export the above mentioned variables and simply run
# make testcluster | check/jenkins_check_results.sh
# NOTE:
# Provide either a VERSION or an EXECUTABLE (EXECUTABLE will be used over VERSION)

# This script reads stdout from make testcluster, parses the slurm job ids, and queues jenkins_failcheck.sh
# to run after the make testcluster jobs finish. The jenkins_failcheck script waits for 5 seconds, then
# runs ./evalcheck_cluster.sh and greps for fails, among other things
# optional: VERSION specifies the scip binary which should be used, PERF=performance enables rubberband support in jenkins_failcheck.sh.
# The results are uploaded to rubberband with rbcli and if there are fails, an email is sent to the admin.

# set up environment for jenkins_failcheck.sh
if [ "${TESTSET}" == "" ]; then
  TESTSET=short
  echo "No testset provided, defaulting to '${TESTSET}'."
fi
if [ "${OUTPUTDIR}" == "" ]; then
  OUTPUTDIR=results
  echo "No no outputdir provided, defaulting to '${OUTPUTDIR}'."
fi
if [ "${SETTING}" == "" ]; then
  SETTING=default
  echo "No setting provided, defaulting to '${SETTING}'."
fi
if [ "${EXECUTABLE}" == "" ]; then
  if [ "${VERSION}" == "" ]; then
    EXECUTABLE=bin/scip
    echo "Neither version nor executable provided, defaulting to '${EXECUTABLE}'."
  else
    EXECUTABLE=`ls bin/scip-${VERSION}*|head -n 1`
  fi
fi
echo "Using executable '${EXECUTABLE}'."

# exporting the variables to the environment for check/jenkins_failcheck.sh to use
export TESTSET
export OUTPUTDIR
export SETTING
export EXECUTABLE

# get some relevant information
# process optional variables
if [ "${PERF}" != "" ]; then
  export PERFORMANCE=${PERF}
fi

export SCIPVERSIONOUTPUT=`${EXECUTABLE} -v | sed -e 's/$/@/'`
export SCIPVERSION=scip-`echo ${SCIPVERSIONOUTPUT} | sed -e 's/.* VERSION=\([^@]*\).*/\1/'`
if [ "${OPT}" == "" ]; then
  OPT=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* OPT=\([^@]*\).*/\1/'`
fi
if [ "${LPS}" == "" ]; then
  LPS=`echo $SCIPVERSIONOUTPUT | sed -e 's/.* LPS=\([^@]*\).*/\1/'`
fi
export OPT
export LPS

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
    export GITBRANCH=`git show -s --pretty=%D | cut -d , -f 2 | cut -d / -f 2 | `
fi

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

env

# build job ids string for sbatch dependency
jobidsstr=$(printf ",%s" "${slurmjobids[@]}")
jobidsstr=${jobidsstr:1}

# execute checker after all jobs completed
#sbatch --dependency=afterany:${jobidsstr} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=4000 --time=100 --partition=mip-dbg --account=mip check/jenkins_failcheck.sh
