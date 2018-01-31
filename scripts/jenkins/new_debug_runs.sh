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

# Arguments | defaultvalue                          | possibilities
# ----------|---------------------------------------|--------------
# GITBRANCH    | master                                | master, bugfix
# SPX_DIR   | /OPTI/adm_timo/soplex_${GITBRANCH}_Debug | *

######################################
### evaluate commandline arguments ###
######################################

# set default arguments
GITBRANCH=master

# This soplex there is installed on pushes to soplex by the jenkins job SOPLEX_install_${GITBRANCH}.
SPX_DIR=/OPTI/adm_timo/soplex_${GITBRANCH}_Debug

# evaluate commandline arguments
for i in $@
do
  eval $i
done

# Find out what day of week it is: mon-1 .. sun-7
DAY_OF_WEEK=`date +%u`

##############################################
### jobs configuration variables ###
##############################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - Jobs indices start at DAY_OF_WEEK,1 and not at zero.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.

SCIP_FLAGS="IPOPT=true SYM=bliss ZIMPL=false COMP=gnu OPT=opt"
RANDOMSEED=`datei +%Y%m%d`

# use associative arrays, this requires bash4
# declaration
declare -A JOBS

# jobs running on monday
JOBS[1,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
JOBS[1,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"
N_JOBS[1]=2

# jobs running on tuesday
JOBS[2,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[2,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[2]=2

# jobs running on wednesday
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=presolaggr_sepaaggr_heuroff_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_presolaggr_sebaaggr_heuroff_${RANDOMSEED}"
N_JOBS[3]=2

# jobs running on thursday
JOBS[4,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=heuraggr_${RANDOMSEED}"
JOBS[4,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_heuraggr_${RANDOMSEED}"
N_JOBS[4]=2

# jobs running on friday
JOBS[5,1]="LPS=cpx MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default"
JOBS[5,2]="LPS=cpx MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default"
N_JOBS[5]=2

# jobs running on saturday
JOBS[6,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdev-solvable TIME=7200 SETTING=default"
N_JOBS[6]=2

# jobs running on sunday
JOBS[7,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=minlpdev-solvable TIME=7200 SETTING=default"
N_JOBS[7]=3

# MIP settings
./bin/scip -c "set rand rand ${RANDOMSEED} set diffsave settings/default_${RANDOMSEED}.set q"
./bin/scip -c "set heur emph aggr set rand rand ${RANDOMSEED} set diffsave settings/heuraggr_${RANDOMSEED}.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set rand rand ${RANDOMSEED} set diffsave settings/presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

# MINLP settings
./bin/scip -c "set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_default.set q"
./bin/scip -c "set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_Default_${RANDOMSEED}.set q"
./bin/scip -c "set heur emph aggr set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_heuraggr_${RANDOMSEED}.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set numerics checkfeastolfac 1000.0 set rand rand ${RANDOMSEED} set diffsave settings/minlp_presolaggr_sepaaggr_heuroff_${RANDOMSEED}.set q"

##############################################
### process variables ###
##############################################

# To improve accessibility move todays jobs into seperate array
TODAYS_N_JOBS=${N_JOBS[$DAY_OF_WEEK]}

declare -A TODAYS_JOBS

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  TODAYS_JOBS[$i]=${JOBS[${DAY_OF_WEEK},$i]}
done


# Print some information about what is happening
echo "Today is `date +%A`. Running the following jobs (index ${DAY_OF_WEEK},*):"
for i in `seq 1 ${TODAYS_N_JOBS}`; do
  echo "- job configuration: '${TODAYS_JOBS[$i]}'"
done

exit 0

##############################################
### Setup compilation ###
##############################################

# create all required directories
mkdir -p lib/include
mkdir -p lib/static
mkdir -p settings

# create all required symlinks
ln -s ${SPX_DIR}/src lib/include/spxinc
ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.opt.a

ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt

ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a

###################
### Compilation ###
###################

# compile, while doing that say 'yes' to all questions
yes "" | make $SCIP_FLAGS -j4 USRFLAGS=-Werror

#########################
### Setup testruns ###
#########################

# create more required symlinks
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

##########################
### Testruns ###
##########################

for i in `seq 1 ${TODAYS_N_JOBS}`; do
  FLAGS="${SCIP_FLAGS} ${TODAYS_JOBS[$i]}"
  echo "Submitting job with configuration: '${FLAGS}'"
  make testcluster ${FLAGS} | ${FLAGS} check/jenkins_check_results.sh
done
