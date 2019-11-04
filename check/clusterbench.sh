#!/bin/bash

# Usage:
# clusterbench.sh QUEUES EXECUTABLE

# benchmarks the performance of all nodes of a given slurm queue or more (comma separated)
# on each node run exclusively a predefined test set with a specified SCIP binary;

# Careful: WE ARE IN check/ !

# arguments and defaults
: ${CB_ID:=$(date +%Y%m%d%H%M)}
: ${CB_OUTPUTDIR:="results/clusterbench_${CB_ID}"}

if [ $# -ge 3 ] || [ $# -eq 0 ]; then
  echo 'usage: ./clusterbench.sh QUEUES EXECUTABLE'
  exit 1
fi;

if test -z ${USER}; then
  echo "The USER environment variable cannot be empty."
  exit 1
fi

if [ ! "$(readlink -f `pwd`)" = "$(dirname $(readlink -f ${BASH_SOURCE[0]}))" ]; then
  echo 'cd to scip/check in order to run script'
  exit 1
fi;

if ! test -d "IP"; then
  echo "> Cannot find directory IP in check/, exiting."
  exit 1
fi;

# check if the test set file exists
TEST="clusterbench"
if [ ! -e "./testset/${TEST}.test" ]; then
  echo "> Cannot find ${TEST}.test in $(pwd)/testset/, exiting."
  exit 1
fi;


# communicate some settings to configuration_logfiles.sh that this is a clusterbenchmark,
# so that important metainfo can be written
export CLUSTERBENCHMARK=yes
export CB_ID=${CB_ID}
export CB_OUTPUTDIR=${CB_OUTPUTDIR}

export EXECUTABLE=$(pwd)/$2
export QUEUE=$1
export SETTINGS="clusterbench"

# from scip/check switch to the SCIP base directory to be able to execute make targets
cd ..

mkdir -p settings
touch settings/${SETTINGS}.set

# we need a counter for the jobids array
j=0

for Q in $(echo ${QUEUE} | sed -e 's/,/ /g'); do
  echo "> Preparing for queue $Q."

  nodes=( $(sinfo --format="%n" -p $Q|tail -n +2) )
  echo "> Found ${#nodes[@]} nodes (${nodes[*]})."

  echo "> Start testing ..."

  echo "make testcluster"
  echo "     EXECUTABLE=${EXECUTABLE}"
  echo "     QUEUE=${Q}"
  echo "     EXCLUSIVE=true"
  echo "     TEST=${TEST}"
  echo "     TIME=600"
  echo "     SETTINGS=${SETTINGS}"

  # execute the testing script on all nodes
  for n in ${nodes[*]}; do
    # OUTPUTDIR is needed in jenkins_failcheck_benchmark.sh
    export OUTPUTDIR="${CB_OUTPUTDIR}/${Q}/${n}"
    export CB_QUEUENODE=${n}

    echo "     CLUSTERNODES=${n}"
    echo "     OUTPUTDIR=${OUTPUTDIR}"

    mkdir -p check/${OUTPUTDIR}

    # run full test set on each node
    TESTRUN_OUTPUT=$(make testcluster EXECUTABLE=${EXECUTABLE} QUEUE=${Q} CLUSTERNODES=${n} SLURMACCOUNT=scip EXCLUSIVE=true TEST=${TEST} SETTINGS=${SETTINGS} TIME=601 OUTPUTDIR="${OUTPUTDIR}" | check/jenkins_check_results_benchmark.sh)

    # use for debugging
    #TESTRUN_OUTPUT=$(make testcluster EXECUTABLE=${EXECUTABLE} QUEUE=${Q} CLUSTERNODES=${n} SLURMACCOUNT=mip EXCLUSIVE=true TEST=${TEST} SETTINGS=${SETTINGS} TIME=10 OUTPUTDIR="${OUTPUTDIR}" | check/jenkins_check_results_benchmark.sh)
    echo $TESTRUN_OUTPUT

    # collect evaljobids; array declared implicitly
    JOBIDS[$j]=$(echo ${TESTRUN_OUTPUT} | cut -d ' ' -f 4)
    ((j++))

  done

done

# go back to scip/check
cd check

# build job ids string for sbatch dependency
JOBIDSSTR=$(printf ",%s" "${JOBIDS[@]}")
JOBIDSSTR=${JOBIDSSTR:1}

# execute script to glue and upload after all eval jobs completed
sbatch --dependency=afterany:${JOBIDSSTR} --kill-on-invalid-dep=yes --cpus-per-task=1 --mem=4000 --time=500 --partition=opt --account=scip jenkins_collect_benchmark.sh
