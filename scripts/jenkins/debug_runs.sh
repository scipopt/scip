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
#        ./debug_runs.sh GITBRANCH=master

# Arguments | defaultvalue                             | possibilities
# ----------|------------------------------------------|--------------
# GITBRANCH | master                                   | master, bugfix
# SPX_DIR   | /OPTI/adm_timo/soplex_${GITBRANCH}_Debug | *

echo "This is debug_runs.sh running."

######################################
### evaluate commandline arguments ###
######################################

# set default arguments
# If no branch is given, try to guess the branch based on the current directory
if [ "${GITBRANCH}" == "" ]; then
  # GIT_BRANCH is a jenkins variable, if not present, try to get it from the git repository. The second thing is not robust because there may be more branches that this HEAD is present in.
  GITBRANCH=`echo ${GIT_BRANCH} | cut -d / -f 2`
  if [ "${GITBRANCH}" == "" ]; then
      GITBRANCH=`git show -s --pretty=%D | cut -d , -f 2 | cut -d / -f 2 | `
  fi
fi

if [ "${GITBRANCH}" != "master" ]; then
  if [[ ${GITBRANCH} =~ "bugfix" ]]; then
    echo "Branch is neither 'master' nor 'bugfix'. Something is wrong. exiting."
    exit 1
  fi
fi

# This soplex there is installed on pushes to soplex by the jenkins job SOPLEX_install_${GITBRANCH}.
SOPLEX_DIR=/OPTI/adm_timo/soplex_${GITBRANCH}_Debug/

CRITERION_DIR=/optimi/usr/sw/criterion
IPOPT_DIR=/optimi/usr/sw/ipopt
CPLEX_DIR=/optimi/usr/sw/cplex
BLISS_DIR=/optimi/usr/sw/bliss

# evaluate commandline arguments
for i in $@
do
  eval $i
done

# Find out what day of week it is: mon-1 .. sun-7
DAY_OF_WEEK=`date +%u`

####################################
### jobs configuration variables ###
####################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - Jobs indices start at DAY_OF_WEEK,1 and not at zero.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.
#  - To add settings please visit the section 'setup testruns'. This can only happen after compilation.
#  - Don't add LPS=xxx and LPSOPT=xxx but instead use EXECUTABLE=[scipspx|scipcpx].
# FORMAT:
#    JOBS[x,y]="EXECUTABLE=scipspx MEM=100 QUEUE=opt TESTSET=short TIME=10 PERMUTE=2 PERFORMANCE=performance"

RANDOMSEED=`date +%Y%m%d`

# use associative arrays, this requires bash4
# declaration
declare -A JOBS

# # jobs running on monday
# JOBS[1,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
# JOBS[1,2]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"
#
# # jobs running on tuesday
# JOBS[2,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
# JOBS[2,2]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
#
# # jobs running on wednesday
# JOBS[3,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=presolaggr_sepaaggr_heuroff_${RANDOMSEED}"
# JOBS[3,2]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_presolaggr_sebaaggr_heuroff_${RANDOMSEED}"
#
# # jobs running on thursday
# JOBS[4,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=heuraggr_${RANDOMSEED}"
# JOBS[4,2]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_heuraggr_${RANDOMSEED}"
#
# # jobs running on friday
# JOBS[5,1]="EXECUTABLE=scipcpx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
# JOBS[5,2]="EXECUTABLE=scipcpx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"
#
# # jobs running on saturday
# JOBS[6,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=mipdev-solvable TIME=7200 SETTING=default"
#
# # jobs running on sunday
# JOBS[7,1]="EXECUTABLE=scipspx MEM=6000 QUEUE=opt TESTSET=minlpdev-solvable TIME=7200 SETTING=default"

JOBS[4,1]="EXECUTABLE=scipspx MEM=100 QUEUE=opt TESTSET=short TIME=10 PERMUTE=2 PERFORMANCE=performance"
JOBS[4,2]="EXECUTABLE=scipcpx MEM=100 QUEUE=opt TESTSET=short TIME=10"

#########################
### process variables ###
#########################

# To improve accessibility move todays jobs into seperate array
TODAYS_N_JOBS=0

# NOTE: only check up to 10 runs. If there are more there is something wrong...
for i in `seq 1 10`; do
  if [ "${JOBS[${DAY_OF_WEEK},$i]}" == "" ]; then
    break
  fi
  TODAYS_N_JOBS=$i
done
if [ "${TODAYS_N_JOBS}" == "0" ]; then
  echo "No schedules for today! Exiting."
  exit 0
fi

declare -A TODAYS_JOBS

LPSVERSIONS=""
for i in `seq 1 ${TODAYS_N_JOBS}`; do
  TODAYS_JOBS[$i]="${JOBS[${DAY_OF_WEEK},$i]} OUTPUTDIR=results$i"

  # Collect versions to know which lpsolvers to compile
  LPSVERSION=`echo ${TODAYS_JOBS[$i]} |grep -o 'scip[a-z]*'`
  LPSVERSIONS="${LPSVERSIONS} ${LPSVERSION}"

  # append /bin/scip to executable
  TODAYS_JOBS[$i]=`echo ${TODAYS_JOBS[$i]}|sed "s@\(scip[cs]px\)@\1/bin/scip@"`
done

# Print some information about what is happening
echo "Today is `date +%A`. Running the following ${TODAYS_N_JOBS} jobs (index ${DAY_OF_WEEK},*):"
for i in `seq 1 ${TODAYS_N_JOBS}`; do
  echo "- job configuration: '${TODAYS_JOBS[$i]}'"
done

#########################
### Setup compilation ###
#########################

# create all required directories
mkdir -p settings

###################
### Compilation ###
###################

# compile, while doing that say 'yes' to all questions

BUILD_DIR=""
# build with soplex only if today we have some soplex runs scheduled
# that is the case if $LPSVERSIONS contains 'scipspx'
if [[ ${LPSVERSIONS} =~ "scipspx" ]] ; then
  BUILD_DIR=scipspx
  mkdir -p ${BUILD_DIR}
  cd ${BUILD_DIR}
  cmake .. -DCMAKE_BUILD_TYPE=Debug -DLPS=spx
  make -j4
  cd ..
fi
# build with cplex only if today we have some cplex runs scheduled
# that is the case if $LPSVERSIONS contains 'scipcpx'
if [[ ${LPSVERSIONS} =~ "scipcpx" ]] ; then
  BUILD_DIR=scipcpx
  mkdir -p ${BUILD_DIR}
  cd ${BUILD_DIR}
  cmake .. -DCMAKE_BUILD_TYPE=Debug -DLPS=cpx
  make -j4
  cd ..
fi

######################
### Setup testruns ###
######################

SCIP_BINARY=${BUILD_DIR}/bin/scip

# NOTES:
#  - When building a default setting with random seed, use a capital D. No setting name should be a prefix of another!

# MIP settings
${SCIP_BINARY} -c "set rand rand ${RANDOMSEED} set diffsave settings/Default_${RANDOMSEED}.set q"
${SCIP_BINARY} -c "set heur emph aggr set rand rand ${RANDOMSEED} set diffsave settings/heuraggr_${RANDOMSEED}.set q"
${SCIP_BINARY} -c "set sepa emph aggr set presol emph aggr set heur emph off set rand rand ${RANDOMSEED} set diffsave settings/presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

# MINLP settings
${SCIP_BINARY} -c "set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_default.set q"
${SCIP_BINARY} -c "set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_Default_${RANDOMSEED}.set q"
${SCIP_BINARY} -c "set heur emph aggr set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_heuraggr_${RANDOMSEED}.set q"
${SCIP_BINARY} -c "set sepa emph aggr set presol emph aggr set heur emph off set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

# create more required symlinks
ln -fs /optimi/kombadon/IP check/
ln -fs /optimi/kombadon/MINLP check/

#######################
### Submit Testruns ###
#######################

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  FLAGS=${TODAYS_JOBS[$i]}
  export ${FLAGS}
  echo "Submitting job with configuration:\n- compilation: ${SCIPFLAGS}'\n- make testcluster: ${FLAGS}"
  make testcluster ${FLAGS} | check/jenkins_check_results_cmake.sh
done

