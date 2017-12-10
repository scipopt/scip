#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

###################
### Description ###
###################
# This script is used by cijenkins.zib.de.
# It builds optimized SCIP/SoPlex with make and runs mipdev-solvable, sapdev-solvable and minlpdev-solvable

# Usage: from scip root execute
#        ./weeklyperf_testruns.sh SPX_DIR=path/to/soplex

# Arguments | defaultvalue | possibilities
# ----------|--------------|--------------
# SPX_DIR   | --           | --


##########################
### evaluate arguments ###
##########################
for i in $@
do
    eval $i
done

#############
### Setup ###
#############
mkdir -p lib/include
mkdir -p lib/static
mkdir -p settings

ln -s ${SPX_DIR}/src lib/include/spxinc
ln -s ${SPX_DIR}/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.opt.a

ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt

ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a

SCIP_FLAGS="IPOPT=true SYM=bliss ZIMPL=false COMP=gnu LPS=spx OPT=opt"

###################
### Compilation ###
###################
yes "" | make $SCIP_FLAGS -j4 USRFLAGS=-Werror

######################
### Testrun: setup ###
######################
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

##########################
### Testrun: execution ###
##########################
for testset in mipdev-solvable; do
    for setting in default; do
        make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting EXCLUSIVE=true TIME=7200 MEM=35000 QUEUE=M620v3-low | PERF=performance check/jenkins_check_results.sh $testset $setting
    done
done

cp ~/sap-400-pure.set settings/.
for testset in sapdev-solvable; do
    for setting in sap-400-pure; do
        make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting EXCLUSIVE=true TIME=7200 MEM=35000 QUEUE=M630v2 | PERF=performance check/jenkins_check_results.sh $testset $setting
    done
done

echo "limits/gap = 0.0001" > settings/minlp.set
echo "numerics/checkfeastolfac = 10.0" >> settings/minlp.set
for testset in minlpdev-solvable; do
    for setting in minlp; do
        make testcluster $SCIP_FLAGS TEST=$testset SETTINGS=$setting EXCLUSIVE=true TIME=3600 MEM=35000 QUEUE=M630-low | PERF=performance check/jenkins_check_results.sh $testset $setting
    done
done
