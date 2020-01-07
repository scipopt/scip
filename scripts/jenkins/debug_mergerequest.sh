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
#        TESTMODE=short GITBRANCH=master ./scripts/jenkins/performance_mergerequest.sh

# Arguments | defaultvalue                             | possibilities
# ----------|------------------------------------------|--------------
# GITBRANCH | master                                   | master, v60-bugfix, consexpr
# TESTMODE  | ""                                       | minlp, short

echo "This is performance_mergerequest.sh running."
: ${TESTMODE:="minlp"}
: ${GITBRANCH:=${gitlabTargetBranch}}

ORIGBRANCH=${GITBRANCH}

if [ "${TESTMODE}" == "short" ]; then
  echo "Testing short"
elif [ "${TESTMODE}" == "minlp" ]; then
  echo "Testing minlp"
else
  echo "Nothing to do, exiting."
  exit 0
fi

######################################
### evaluate commandline arguments ###
######################################

if [ "${GITBRANCH}" != "master" ]; then
  if [ "${GITBRANCH}" != "consexpr" ]; then
    if [[ ${GITBRANCH} =~ "bugfix" ]]; then
      GITBRANCH=bugfix
    else
      export FAILREASON="Branch is neither 'master', 'bugfix' nor 'consexpr'. Something is wrong. Exiting."
      echo ${FAILREASON}
      exit 1
    fi
  fi
fi

export FULLGITHASH=$(git show -s --pretty=%H)
export MODE=debug
export MRSETTINGS="MR-${gitlabMergeRequestIid}"

####################################
### jobs configuration variables ###
####################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.
#  - To add settings please visit the section 'setup testruns'. This can only happen after compilation.
#  - Don't add LPS=xxx and LPSOPT=xxx but instead use EXECUTABLE=[scipdbgspx|scipdbgcpx].
#  - Only 10 runs will be executed. If you need more you should overthink you overall concept.
#  - The check/jenkins_*_cmake.sh evaluation scripts don't work yet if you use a global seed shift.
# FORMAT:
#    JOBS[x,y]="EXCLUSIVE=true EXECUTABLE=scipoptspx/bin/scip BINID=scipoptspx-${GITBRANCH} MEM=100 QUEUE=opt TEST=short TIME=10 PERMUTE=2 SETTINGS=default PERFORMANCE=mergerequest"

RANDOMSEED=$(date +%Y%m%d%H%M)

# for descriptions on the testsets see scip/check/testsets/README.md
# jobs running

if [ "${TESTMODE}" == "short" ]; then
  JOB="EXECUTABLE=scipoptcpx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptcpx_${GITBRANCH}_${RANDOMSEED} EXCLUSIVE=false MEM=50 QUEUE=opt TEST=short TIME=10 SETTINGS=${MRSETTINGS} PERFORMANCE=mergerequest SEEDS=0"
elif [ "${TESTMODE}" == "minlp" ]; then
  JOB="EXECUTABLE=scipoptcpx_${GITBRANCH}_${RANDOMSEED}/bin/scip BINID=scipoptcpx_${GITBRANCH}_${RANDOMSEED} EXCLUSIVE=false MEM=6000 QUEUE=M620,M620v2,M630,M630v2 TEST=MINLP_debug TIME=60 SETTINGS=${MRSETTINGS} PERFORMANCE=mergerequest"
fi


# create required directory
mkdir -p settings
JOB="${JOB} OUTPUTDIR=results${RANDOMSEED}"
FLAGS=${JOB}
export ${FLAGS}
export PERFORMANCE=mergerequest

# from check/jenkins_failcheck.sh
SCIP_BUILDDIR=$(echo ${EXECUTABLE}| cut -d '/' -f 1|cut -d '_' -f 1)

# get the comparison db for githash
if [ "${TESTMODE}" == "short" ]; then
  COMPARERBDB="/nfs/OPTI/adm_timo/databases/rbdb/${GITBRANCH}_${MODE}_${TEST}_${SETTINGS}_${SCIP_BUILDDIR}_rbdb.txt"
elif [ "${TESTMODE}" == "minlp" ]; then
  COMPARERBDB="/nfs/OPTI/adm_timo/databases/rbdb/${GITBRANCH}_${MODE}_${TEST}_${SETTINGS}_${SCIP_BUILDDIR}_rbdb.txt"
fi

export CRITERION_DIR=""
export BLISS_DIR=/nfs/OPTI/bzfgleix/software/bliss-0.73p-Ubuntu18.04
export IPOPT_DIR=/nfs/optimi/usr/sw/ipopt-static
export ZIMPL_DIR=/nfs/OPTI/jenkins/workspace/ZIMPL_monthly/build-gnu-Release/

###################
### Compilation ###
###################

make solchecker

# build with cplex
BUILD_DIR=scipoptcpx_${GITBRANCH}_${RANDOMSEED}
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake .. -DCMAKE_BUILD_TYPE=Debug -DLPS=cpx
make -j4
cd ..

######################
### Setup testruns ###
######################

SCIP_BINARY=${BUILD_DIR}/bin/scip

# NOTES:
#  - When building a default setting with random seed, use a capital D. No setting name should be a prefix of another!

# MIP settings
touch "settings/${MRSETTINGS}.set"

# MINLP settings
${SCIP_BINARY} -c "set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_${MRSETTINGS}.set q"

# create more required symlinks
ln -fs /nfs/optimi/kombadon/IP check/
ln -fs /nfs/optimi/kombadon/MINLP check/

# get testset files to the correct place
cp check/IP/instancedata/testsets/*.test check/testset/

#######################
### Submit Testrun ###
#######################

echo "Submitting job with configuration:\n '\n- make testcluster: ${FLAGS}"
make testcluster ${FLAGS} | check/jenkins_check_results_cmake.sh

