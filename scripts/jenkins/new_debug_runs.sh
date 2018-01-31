#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

# NOTES:
#  - We use associative arrays, this requires bash4.

###################
### Description ###
###################

# This script is used by cijenkins.zib.de.
# Depending on the day of the week this script will start different testruns on the cluster
# TODO

# Usage: from scip root execute
#        ./script.sh GITBRANCH=master

# Arguments | defaultvalue                             | possibilities
# ----------|------------------------------------------|--------------
# GITBRANCH | master                                   | master, bugfix
# SPX_DIR   | /OPTI/adm_timo/soplex_${GITBRANCH}_Debug | *

######################################
### evaluate commandline arguments ###
######################################

# set default arguments
GITBRANCH=master

# This soplex there is installed on pushes to soplex by the jenkins job SOPLEX_install_${GITBRANCH}.
# SPX_DIR=/nfs/OPTI/jenkins/workspace/SOPLEX_COMP=gnu_OPT=dbg_nightly #TODO
SPX_DIR=/OPTI/adm_timo/soplex_${GITBRANCH}_Debug

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
#  - Don't add LPS=xxx and LPSOPT=xxx but instead use VERSION=[scipspxopt|scipspxdbg|scipcpx].

SCIP_FLAGS="COMP=gnu IPOPT=true OPT=dbg SYM=bliss USRFLAGS=-Werror ZIMPL=false"
RANDOMSEED=`date +%Y%m%d`

# use associative arrays, this requires bash4
# declaration
declare -A JOBS

# jobs running on monday
JOBS[1,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
JOBS[1,2]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"

# jobs running on tuesday
JOBS[2,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[2,2]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"

# jobs running on wednesday
JOBS[3,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=presolaggr_sepaaggr_heuroff_${RANDOMSEED}"
JOBS[3,2]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_presolaggr_sebaaggr_heuroff_${RANDOMSEED}"

# jobs running on thursday
JOBS[4,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=heuraggr_${RANDOMSEED}"
JOBS[4,2]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_heuraggr_${RANDOMSEED}"

# jobs running on friday
JOBS[5,1]="VERSION=scipcpx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
JOBS[5,2]="VERSION=scipcpx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"

# jobs running on saturday
JOBS[6,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=mipdev-solvable TIME=7200 SETTING=default"

# jobs running on sunday
JOBS[7,1]="VERSION=scipspxdbg MEM=6000 QUEUE=opt TESTSET=minlpdev-solvable TIME=7200 SETTING=default"

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

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  TODAYS_JOBS[$i]="${JOBS[${DAY_OF_WEEK},$i]} OUTPUTDIR=results$i"
done

# Print some information about what is happening
# Collect versions to know which combinations to compile
LPSVERSIONS=""
echo "Today is `date +%A`. Running the following ${TODAYS_N_JOBS} jobs (index ${DAY_OF_WEEK},*):"
for i in `seq 1 ${TODAYS_N_JOBS}`; do
  echo "- job configuration: '${TODAYS_JOBS[$i]}'"
  LPSVERSIONS="${LPSVERSIONS}`echo ${TODAYS_JOBS[$i]} |grep -o 'VERSION=scip[a-z]*'` "
done

#########################
### Setup compilation ###
#########################

# create all required directories
mkdir -p lib/include
mkdir -p lib/include/zimplinc
mkdir -p lib/static
mkdir -p settings

# create all required symlinks
# symlink spx
ln -s ${SPX_DIR}/src lib/include/spxinc
ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.opt.a
ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.dbg.a lib/static/libsoplex.linux.x86_64.gnu.dbg.a

# symlink cpx
ln -s /optimi/usr/sw/cplex/include/ilcplex lib/include/cpxinc
ln -s /optimi/usr/sw/cplex/lib/x86-64_linux/static_pic/libcplex.a lib/static/libcplex.linux.x86_64.gnu.a

# symlink ipopt
ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt

# symlink bliss
ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a

# symlink zimpl
ln -s /optimi/usr/sw/zimpl/src lib/include/zimplinc/zimpl
ln -s /optimi/usr/sw/zimpl/lib/libzimpl.linux.x86_64.gnu.opt.a lib/static/libzimpl.linux.x86_64.gnu.opt.a

###################
### Compilation ###
###################

# compile, while doing that say 'yes' to all questions

if [[ ${LPSVERSIONS} =~ "scipspxopt" ]] ; then
  yes "" | make ${SCIP_FLAGS} LPS=spx LPSOPT=opt VERSION=scipspxopt -j4
fi
if [[ ${LPSVERSIONS} =~ "scipspxdbg" ]] ; then
  yes "" | make ${SCIP_FLAGS} LPS=spx LPSOPT=dbg VERSION=scipspxdbg -j4
fi
if [[ ${LPSVERSIONS} =~ "scipcpx" ]] ; then
  yes "" | make ${SCIP_FLAGS} LPS=cpx VERSION=scipcpx -j4
fi

######################
### Setup testruns ###
######################

# NOTES:
#  - When building a default setting with random seed, use a capital D. No setting name should be a prefix of another!

# MIP settings
./bin/scip -c "set rand rand ${RANDOMSEED} set diffsave settings/Default_${RANDOMSEED}.set q"
./bin/scip -c "set heur emph aggr set rand rand ${RANDOMSEED} set diffsave settings/heuraggr_${RANDOMSEED}.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set rand rand ${RANDOMSEED} set diffsave settings/presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

# MINLP settings
./bin/scip -c "set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_default.set q"
./bin/scip -c "set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_Default_${RANDOMSEED}.set q"
./bin/scip -c "set heur emph aggr set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_heuraggr_${RANDOMSEED}.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

# create more required symlinks
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

#######################
### Submit Testruns ###
#######################

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  FLAGS="${TODAYS_JOBS[$i]}"
  echo "Submitting job with configuration:\n- compilation: '${SCIPFLAGS}'\n- make testcluster: ${FLAGS}'"
  echo "make testcluster ${FLAGS} | ${FLAGS} check/jenkins_check_results.sh"
done
