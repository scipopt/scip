#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

# NOTE: We use associative arrays, this requires bash4.

###################
### Description ###
###################
# This script is used by cijenkins.zib.de.
# Depending on the day of the week this script will start different testruns on the cluster
# TODO

# Usage: from scip root execute
#        ./TODO

# Arguments | defaultvalue | possibilities
# ----------|--------------|--------------
# TODO

######################################
### evaluate commandline arguments ###
######################################
# TODO
# for i in $@
# do
#   eval $i
# done

# Find out what day of week it is: mon-1 .. sun-7
DAY_OF_WEEK=`date +%u`
BRANCH=master #TODO

##############################################
### jobs configuration variables ###
##############################################
# NOTES:
#  - If you change the configuration, you have to make sure that you update the number of jobs in the N_JOBS array.
#  - Jobs indices start at DAY_OF_WEEK,1 and not at zero.
#  - For all jobs the calls to 'make' and 'make testcluster' the flags are concatenated from
#      the given flags and the SCIP_FLAGS.

SCIP_FLAGS="IPOPT=true SYM=bliss ZIMPL=false COMP=gnu OPT=opt"
SPX_DIR=/optimi/home/adm_timo/spx_master_Release #TODO
RANDOMSEED=`datei +%Y%m%d` #TODO

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
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[3]=3

# jobs running on thursday
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[3]=3

# jobs running on friday
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[3]=3

# jobs running on saturday
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[3]=3

# jobs running on sunday
JOBS[3,1]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=mipdebug TIME=60 SETTING=default_${RANDOMSEED}"
JOBS[3,2]="LPS=spx LPSOPT=dbg MEM=6000 QUEUE=opt TESTSET=MINLP TIME=60 SETTING=minlp_default_${RANDOMSEED}"
N_JOBS[3]=3

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
  echo "- jobconfig: '${TODAYS_JOBS[$i]}'"
done

exit 0

##############################################
### Setup compilation ###
##############################################
mkdir -p lib/include
mkdir -p lib/static
mkdir -p settings

ln -s ${SPX_DIR}/src lib/include/spxinc
ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.opt.a

ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt

ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a

###################
### Compilation ###
###################
yes "" | make $SCIP_FLAGS -j4 USRFLAGS=-Werror

######################
### setup of testruns ###
######################
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

##########################
### execution of testruns ###
##########################
for testset in mipdev-solvable; do
    for setting in default; do
        make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting EXCLUSIVE=true TIME=7200 MEM=35000 QUEUE=M620v3-low | PERF=performance check/jenkins_check_results.sh $testset $setting
    done
done
