#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

# NOTES:
#  - We use associative arrays, this requires bash4.

###################
### Description ###
###################

# This script is used by cijenkins.zib.de.

# Usage: from scip root execute
#        GITBRANCH=master ./performance_mergerequest.sh

# Arguments | defaultvalue                             | possibilities
# ----------|------------------------------------------|--------------
# GITBRANCH | master                                   | master, bugfix

echo "This is performance_mergerequest.sh running."
: ${TESTMODE:="all"}

if [ "${TESTMODE}" == "all" ]; then
  #echo "Testing mip, minlp and sapdev-solvable"
  echo "Testing mip and minlp"
elif [ "${TESTMODE}" == "short" ]; then
  echo "Testing short"
elif [ "${TESTMODE}" == "mip" ]; then
  echo "Testing mip"
elif [ "${TESTMODE}" == "minlp" ]; then
  echo "Testing minlp"
# elif [ "${TESTMODE}" == "sapdev-solvable" ]; then
#   echo "Testing sapdev-solvable"
else
  echo "Nothing to do, exiting."
  exit 0
fi

######################################
### evaluate commandline arguments ###
######################################

export GITBRANCH=${gitlabTargetBranch}

if [ "${GITBRANCH}" != "master" ]; then
  if [ "${GITBRANCH}" != "consexpr" ]; then
    if [[ ${GITBRANCH} =~ "bugfix" ]]; then
      GITBRANCH=bugfix
    else
      echo "Branch is neither 'master' nor 'bugfix'. Something is wrong. Exiting."
      exit 1
    fi
  fi
fi

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

####################################
### jobs configuration variables ###
####################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - Jobs indices start at 1 and not at zero.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.
#  - To add settings please visit the section 'setup testruns'. This can only happen after compilation.
#  - Only 10 runs will be executed. If you need more you should overthink you overall concept.
#  - The check/jenkins_*_cmake.sh evaluation scripts don't work yet if you use a global seed shift.
# FORMAT:
#    JOBS[x,y]="EXCLUSIVE=true EXECUTABLE=scipoptspx/bin/scip BINID=scipoptspx-${GITBRANCH} MEM=100 QUEUE=opt TEST=short TIME=10 PERMUTE=2 SETTINGS=default PERFORMANCE=mergerequest"

RANDOMSEED=$(date +%Y%m%d%H%M)

# use associative arrays, this requires bash4
# declaration
declare -A JOBS

# for descriptions on the testsets see scip/check/testsets/README.md
# jobs running

if [ "${TESTMODE}" == "all" ]; then
  JOBS[1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M620v3 TEST=mipdev12merged-solvable TIME=7200 SETTINGS=default PERFORMANCE=mergerequest SEEDS=4"
  JOBS[2]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M640 TEST=minlpdev-solvable TIME=3600 SETTINGS=minlp_default PERFORMANCE=mergerequest PERMUTE=4"
  #JOBS[3]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M630v2 TEST=sapdev-solvable TIME=3600 SETTINGS=${SAPSETTINGS} PERFORMANCE=mergerequest SEEDS=2"
elif [ "${TESTMODE}" == "short" ]; then
  JOBS[1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} EXCLUSIVE=false MEM=5000 QUEUE=opt TEST=short TIME=60 SETTINGS=default PERFORMANCE=mergerequest SEEDS=0"
elif [ "${TESTMODE}" == "mip" ]; then
  JOBS[1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M620v3 TEST=mipdev12merged-solvable TIME=7200 SETTINGS=default PERFORMANCE=mergerequest SEEDS=4"
elif [ "${TESTMODE}" == "minlp" ]; then
  JOBS[1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M640 TEST=minlpdev-solvable TIME=3600 SETTINGS=minlp_default PERFORMANCE=mergerequest PERMUTE=4"
#elif [ "${TESTMODE}" == "sap" ]; then
#  JOBS[1]="EXECUTABLE=scipoptspx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptspx_${GITBRANCH}_${RANDOMSEED} SLURMACCOUNT=scip EXCLUSIVE=true MEM=50000 QUEUE=M630v2 TEST=sapdev-solvable TIME=3600 SETTINGS=${SAPSETTINGS} PERFORMANCE=mergerequest SEEDS=2"
fi

SAPSETTINGS=sap-next-release-pure-diff
if [ "${GITBRANCH}" != "master" ]; then
  SAPSETTINGS=sap-600-pure-diff
fi

# symlink to SAP settings for the next release settings
ln -fs ~/sap-next-release-pure-diff.set settings/.
ln -fs ~/sap-600-pure-diff.set settings/.

#########################
### process variables ###
#########################

# To improve accessibility move todays jobs into separate array
N_JOBS=0

# NOTE: only check up to 10 runs. If there are more there is something wrong...
for i in `seq 1 10`; do
  if [ "${JOBS[$i]}" == "" ]; then
    break
  fi
  N_JOBS=$i
done

declare -A TODAYS_JOBS

for i in `seq 1 ${N_JOBS}`; do
  TODAYS_JOBS[$i]="${JOBS[$i]} OUTPUTDIR=results${RANDOMSEED}_${i}"
done

# Print some information about what is happening
echo "Today is `date +%A`. Running the following ${N_JOBS} jobs:"
for i in `seq 1 ${N_JOBS}`; do
  echo "- job configuration: '${TODAYS_JOBS[$i]}'"
done

if [ "${N_JOBS}" != "0" ]; then
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

  for i in `seq 1 ${N_JOBS}`; do
    FLAGS=${TODAYS_JOBS[$i]}
    for j in "SEEDS EXECUTABLE BINID MEM QUEUE TEST TIME PERMUTE PERFORMANCE EXCLUSIVE SETTINGS OUTPUTDIR"; do
      unset $j
    done
    export ${FLAGS}
    export PERFORMANCE=mergerequest

    echo "Submitting job with configuration:\n- compilation: ${SCIPFLAGS}'\n- make testcluster: ${FLAGS}"
    make testcluster ${FLAGS} | check/jenkins_check_results_cmake.sh
  done
fi

