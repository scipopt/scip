#!/bin/bash -x
# the -x is for writing each command to standard error (preceded by a ‘+ ’) before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'

###################
### Description ###
###################
# This script is used by cijenkins.zib.de
# Builds SCIP with Make and an LPSolver and runs mipdebug and MINLP on different settings

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

SCIP_FLAGS="ZIMPL=true COMP=gnu OPT=dbg IPOPT=true SYM=bliss "

# deal with different LPSolvers
if [ "${LPS}" == "spx" ];
then
    # soplex is in adm_timos jenkins workspace
    ln -s /nfs/OPTI/jenkins/workspace/SOPLEX_COMP=gnu_OPT=dbg_nightly/src lib/include/spxinc
    ln -s /nfs/OPTI/jenkins/workspace/SOPLEX_COMP=gnu_OPT=dbg_nightly/lib/libsoplex.linux.x86_64.gnu.dbg.a lib/static/libsoplex.linux.x86_64.gnu.dbg.a
    SCIP_FLAGS+="LPS=spx LPSOPT=dbg"
elif [ "${LPS}" == "cpx" ];
then
    # cplex is globally installed
    ln -s /optimi/usr/sw/cplex/include/ilcplex lib/include/cpxinc
    ln -s /optimi/usr/sw/cplex/lib/x86-64_linux/static_pic/libcplex.a lib/static/libcplex.linux.x86_64.gnu.a
    SCIP_FLAG+="LPS=cpx"
fi

###################
### Compilation ###
###################
yes "" | make $SCIP_FLAGS -j4 USRFLAGS=-Werror

#####################
### Testrun setup ###
#####################
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

#########################
### Testrun execution ###
#########################
# MIP
./bin/scip -c "set heur emph aggr set diffsave settings/heuraggr.set q"
./bin/scip -c "set sepa emph aggr set diffsave settings/sepaaggr.set q"
./bin/scip -c "set presol emph aggr set diffsave settings/presolaggr.set q"
for testset in mipdebug; do
  for setting in default heuraggr sepaaggr presolaggr; do
    make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting MEM=6000 TIME=60 QUEUE=opt-low | check/jenkins_check_results.sh $testset $setting
  done
done

# MINLP
./bin/scip -c "set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_default.set q"
./bin/scip -c "set heur emph aggr set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_heuraggr.set q"
./bin/scip -c "set sepa emph aggr set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_sepaaggr.set q"
./bin/scip -c "set presol emph aggr set numerics checkfeastolfac 1000.0 set diffsave settings/minlp_presolaggr.set q"
for testset in MINLP; do
  for setting in minlp_default minlp_heuraggr minlp_sepaaggr minlp_presolaggr; do
    make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting MEM=6000 TIME=60 QUEUE=opt-low | check/jenkins_check_results.sh $testset $setting
  done
done
