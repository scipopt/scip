#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

# NOTES:
#  - We use associative arrays, this requires bash4.

###################
### Description ###
###################

# This script is used by cijenkins.zib.de.
# Depending on the day of the week this script will start different testruns on the cluster.

# Usage: from scip root execute
#        GITBRANCH=master ./performance_runs.sh

# Arguments | defaultvalue                             | possibilities
# ----------|------------------------------------------|--------------
# GITBRANCH | master                                   | master, bugfix

echo "This is performance_runs.sh running."

# arguments and defaults
# Find out what day of week it is: mon-1 .. sun-7
: ${DAY_OF_WEEK:=$(date +%u)}

######################################
### evaluate commandline arguments ###
######################################

# set default arguments
# If no branch is given, try to guess the branch based on the current directory
if [ "${GITBRANCH}" == "" ]; then
  # GIT_BRANCH is a jenkins variable, if not present, try to get it from the git repository. The second thing is not robust because there may be more branches that this HEAD is present in.
  GITBRANCH=`echo ${GIT_BRANCH} | cut -d / -f 2`
  if [ "${GITBRANCH}" == "" ]; then
      GITBRANCH=`git show -s --pretty=%D | cut -s -d , -f 2 | cut -d / -f 2`
  fi
fi

if [ "${GITBRANCH}" != "master" ]; then
  if [ "${GITBRANCH}" != "consexpr" ]; then
    if [[ ${GITBRANCH} =~ "bugfix" ]]; then
      GITBRANCH=bugfix
    else
      echo "Branch is neither 'master' nor 'bugfix' nor 'consexpr'. Something is wrong. Exiting."
      exit 1
    fi
  fi
fi

export FULLGITHASH=$(git show -s --pretty=%H)
export GITBRANCH
export MODE=performance

# This soplex there is installed on pushes to soplex by the jenkins job SOPLEX_install_${GITBRANCH}.
# We have to export these variables to make them available to cmake.
# Scripts will also use nonexported variables correctly.
if [ "${GITBRANCH}" == "consexpr" ]; then
  export SOPLEX_DIR=/nfs/OPTI/adm_timo/soplex_master_Release/
else
  export SOPLEX_DIR=/nfs/OPTI/adm_timo/soplex_${GITBRANCH}_Release/
fi

export CRITERION_DIR=""
export BLISS_DIR=/nfs/OPTI/bzfgleix/software/bliss-0.73p-Ubuntu18.04
export IPOPT_DIR=/nfs/optimi/usr/sw/Ipopt-3.12.11~ub18.04
export ZIMPL_DIR=/nfs/OPTI/jenkins/workspace/ZIMPL_monthly/build-gnu-Release/

# create required directory
mkdir -p settings

#######################
### Update Branches ###
#######################

BRANCHNAME=${GITBRANCH}
if [ "${GITBRANCH}" == "bugfix" ]; then
  BRANCHNAME="v60-bugfix"
fi
if [ "${DAY_OF_WEEK}" == "6" ]; then
  git checkout ${BRANCHNAME}
  git pull
  git checkout performance-${GITBRANCH}
  git merge ${BRANCHNAME} --ff-only
  git push
  git checkout ${BRANCHNAME}
fi

####################################
### jobs configuration variables ###
####################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - Jobs indices start at DAY_OF_WEEK,1 and not at zero.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.
#  - To add settings please visit the section 'setup testruns'. This can only happen after compilation.
#  - Only 10 runs per day will be executed. If you need more you should overthink you overall concept.
#  - The check/jenkins_*_cmake.sh evaluation scripts don't work yet if you use a global seed shift.
# FORMAT:
#    JOBS[x,y]="EXCLUSIVE=true EXECUTABLE=scipoptspx/bin/scip BINID=scipoptspx-${GITBRANCH} MEM=100 QUEUE=opt TEST=short TIME=10 PERMUTE=2 SETTINGS=default PERFORMANCE=performance"

RANDOMSEED=$(date +%Y%m%d)

# use associative arrays, this requires bash4
# declaration
declare -A JOBS
declare -A TRIGGER

# SAPSETTINGS
SAPSETTINGS=sap-next-release-pure-diff
if [ "${GITBRANCH}" != "master" ]; then
  SAPSETTINGS=sap-600-pure-diff
fi

# symlink to SAP settings for the next release settings
ln -fs ~/sap-next-release-pure-diff.set settings/.
ln -fs ~/sap-600-pure-diff.set settings/.

# for descriptions on the testsets see scip/check/testsets/README.md

if [ "${GITBRANCH}" == "master" ]; then
  # running on saturday
  JOBS[6,1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M620v3 TEST=mipdev12merged-solvable TIME=7200 SETTINGS=default PERFORMANCE=performance SEEDS=4"
  TRIGGER[6,1]="https://adm_timo:11d1846ee478c8ff7b7e116b4dd0ddbe86@cijenkins.zib.de/job/SCIP_SAP_perfrun_${GITBRANCH}_weekly/build?token=weeklysaptoken"
  TRIGGER[6,2]="https://adm_timo:11d1846ee478c8ff7b7e116b4dd0ddbe86@cijenkins.zib.de/job/SCIP_SAP_perfrun_shootout_weekly/build?token=weeklysaptoken"
  TRIGGER[6,3]="https://adm_timo:11d1846ee478c8ff7b7e116b4dd0ddbe86@cijenkins.zib.de/job/SCIP_SAP_perfrun_presolve_master_weekly/build?token=weeklysaptoken"

  # jobs running on sunday
  #JOBS[7,1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M630v2 TEST=sapdev-solvable TIME=3600 SETTINGS=${SAPSETTINGS} PERFORMANCE=performance SEEDS=2"

elif [ "${GITBRANCH}" == "consexpr" ]; then
  # running on saturday
  JOBS[6,1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M640 TEST=minlpdev-solvable TIME=3600 SETTINGS=minlp_default PERFORMANCE=performance PERMUTE=4"

else # on bugfix
  # running on saturday
  JOBS[6,1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M620v3 TEST=mipdev12merged-solvable TIME=7200 SETTINGS=default PERFORMANCE=performance SEEDS=4"
  JOBS[6,2]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M640 TEST=minlpdev-solvable TIME=3600 SETTINGS=minlp_default PERFORMANCE=performance PERMUTE=4"
  TRIGGER[6,1]="https://adm_timo:11d1846ee478c8ff7b7e116b4dd0ddbe86@cijenkins.zib.de/job/SCIP_SAP_perfrun_${GITBRANCH}_weekly/build?token=weeklysaptoken"

  # jobs running on sunday
  #JOBS[7,1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M630v2 TEST=sapdev-solvable TIME=3600 SETTINGS=${SAPSETTINGS} PERFORMANCE=performance SEEDS=2"

fi

#########################
### process variables ###
#########################

# To improve accessibility move todays jobs into separate array
TODAYS_N_JOBS=0
TODAYS_N_TRIGGERS=0

# NOTE: only check up to 10 runs. If there are more there is something wrong...
for i in `seq 1 10`; do
  if [ "${JOBS[${DAY_OF_WEEK},$i]}" == "" ]; then
    break
  fi
  TODAYS_N_JOBS=$i
done

declare -A TODAYS_JOBS

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  TODAYS_JOBS[$i]="${JOBS[${DAY_OF_WEEK},$i]} OUTPUTDIR=results${RANDOMSEED}_${i}"
done

# Print some information about what is happening
echo "Today is `date +%A`. Running the following ${TODAYS_N_JOBS} jobs (index ${DAY_OF_WEEK},*):"
for i in `seq 1 ${TODAYS_N_JOBS}`; do
  echo "- job configuration: '${TODAYS_JOBS[$i]}'"
done
echo "Triggering the following jobs:"
for i in `seq 1 10`; do
  if [ "${TRIGGER[${DAY_OF_WEEK},$i]}" == "" ]; then
    break
  fi
  echo "- ${TRIGGER[${DAY_OF_WEEK},$i]}"
  TODAYS_N_TRIGGERS=$i
done

# exit if nothing to do
if [ "${TODAYS_N_JOBS}" == "0" -a "${TODAYS_N_TRIGGERS}" == "0" ]; then
  echo "No schedules for today! Exiting."
  exit 0
fi

if [ "${TODAYS_N_JOBS}" != "0" ]; then
  ###################
  ### Compilation ###
  ###################

  make solchecker

  # build with soplex
  BUILD_DIR=scipoptspx_${GITBRANCH}_${RANDOMSEED}
  mkdir -p ${BUILD_DIR}
  cd ${BUILD_DIR}
  cmake .. -DCMAKE_BUILD_TYPE=Release -DLPS=spx -DSOPLEX_DIR=${SOPLEX_DIR}
  make -j4
  cd ..

  # TODO tag this


  ######################
  ### Setup testruns ###
  ######################

  SCIP_BINARY=${BUILD_DIR}/bin/scip

  # NOTES:
  #  - When building a default setting with random seed, use a capital D. No setting name should be a prefix of another!

  # MIP settings

  # MINLP settings
  ${SCIP_BINARY} -c "set numerics checkfeastolfac 1000.0 set limits gap 1e-4 set diffsave settings/minlp_default.set q"

  # create more required symlinks
  ln -fs /nfs/optimi/kombadon/IP check/
  ln -fs /nfs/optimi/kombadon/MINLP check/

  # get testset files to the correct place
  cp check/IP/instancedata/testsets/*.test check/testset/

  #######################
  ### Submit Testruns ###
  #######################

  for i in `seq 1 ${TODAYS_N_JOBS}`; do
    FLAGS=${TODAYS_JOBS[$i]}
    for j in "SEEDS EXECUTABLE BINID MEM QUEUE TEST TIME PERMUTE PERFORMANCE EXCLUSIVE SETTINGS OUTPUTDIR"; do
      unset $j
    done
    export ${FLAGS}
    export PERFORMANCE=performance

    echo "Submitting job with configuration:\n- compilation: ${SCIPFLAGS}'\n- make testcluster: ${FLAGS}"
    make testcluster ${FLAGS} | check/jenkins_check_results_cmake.sh
  done
fi

set +e
if [ "${TODAYS_N_TRIGGERS}" != "0" ]; then
  # NOTE: only check up to 10 triggers. If there are more there is something wrong...
  echo "Will trigger the following jobs:"
  for i in `seq 1 ${TODAYS_N_TRIGGERS}`; do
    curl -f -I "${TRIGGER[${DAY_OF_WEEK},$i]}"
  done
fi
set -e
