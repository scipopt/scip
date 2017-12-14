#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

###################
### Description ###
###################
# This script is used by cijenkins.zib.de.
# It builds SCIP with Make and an LPSolver and runs mipdebug and MINLP on different settings.

# Usage: from scip root execute
#        ./nightly_testruns.sh LPS=spx

# Arguments | defaultvalue | possibilities
# ----------|--------------|--------------
# LPS       | spx          | cpx, spx
# SPX_DIR   | --           | --

##########################
### evaluate arguments ###
##########################
for i in $@
do
    eval $i
done

#########################
### Compilation setup ###
#########################
mkdir -p lib/include
mkdir -p lib/include/zimplinc
mkdir -p lib/static
mkdir -p settings
ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a
ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt
ln -s /optimi/usr/sw/zimpl/src lib/include/zimplinc/zimpl
ln -s /optimi/usr/sw/zimpl/lib/libzimpl.linux.x86_64.gnu.opt.a lib/static/libzimpl.linux.x86_64.gnu.opt.a

SCIP_FLAGS="ZIMPL=true COMP=gnu OPT=dbg IPOPT=true SYM=bliss"

# deal with different LPSolvers
if [ "${LPS}" == "cpx" ];
then
    # cplex is globally installed
    ln -s /optimi/usr/sw/cplex/include/ilcplex lib/include/cpxinc
    ln -s /optimi/usr/sw/cplex/lib/x86-64_linux/static_pic/libcplex.a lib/static/libcplex.linux.x86_64.gnu.a
    SCIP_FLAGS="$SCIP_FLAGS LPS=cpx"
else
    # soplex is in adm_timos jenkins workspace
    ln -s ${SPX_DIR}/src lib/include/spxinc
    ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.dbg.a lib/static/libsoplex.linux.x86_64.gnu.dbg.a
    SCIP_FLAGS="$SCIP_FLAGS LPS=spx LPSOPT=dbg"
fi

###################
### Compilation ###
###################
yes "" | make ${SCIP_FLAGS} -j4 USRFLAGS=-Werror

######################
### Testrun: setup ###
######################
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

##########################
### Testrun: execution ###
##########################

# NOTE: when building a default setting with random seed, use a capital D.
# No setting name should be a prefix of another!

# MIP
./bin/scip -c "set rand rand 20171124 set diffsave settings/Default_20171124.set q"
./bin/scip -c "set heur emph aggr set rand rand 20171124 set diffsave settings/heuraggr_20171124.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set rand rand 20171124 set diffsave settings/presolaggr_sepaaggr_heuroff_20171124.set q"
for testset in mipdebug; do
  for setting in Default_20171124 heuraggr_20171124 presolaggr_sepaaggr_heuroff_20171124 default; do
    make testcluster $SCIP_FLAGS TEST=$testset MEM=6000 TIME=300 SETTINGS=$setting QUEUE=opt-low | check/jenkins_check_results.sh $testset $setting
  done
done

# MINLP
./bin/scip -c "set numerics checkfeastolfac 1000.0 set set diffsave settings/minlp_default.set q"
./bin/scip -c "set numerics checkfeastolfac 1000.0 set rand rand 20171124 set diffsave settings/minlp_Default_20171124.set q"
./bin/scip -c "set heur emph aggr set numerics checkfeastolfac 1000.0 set rand rand 20171124 set diffsave settings/minlp_heuraggr_20171124.set q"
./bin/scip -c "set sepa emph aggr set presol emph aggr set heur emph off set numerics checkfeastolfac 1000.0 set rand rand 20171124 set diffsave settings/minlp_presolaggr_sepaaggr_heuroff_20171124.set q"
for testset in MINLP; do
  for setting in minlp_Default_20171124 minlp_heuraggr_20171124 minlp_presolaggr_sepaaggr_heuroff_20171124 minlp_default; do
    make testcluster $SCIP_FLAGS TEST=$testset MEM=6000 TIME=300 SETTINGS=$setting QUEUE=opt-low | check/jenkins_check_results.sh $testset $setting
  done
done
