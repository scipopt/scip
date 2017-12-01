#!/bin/bash -ex
# the -x is for writing each command to standard error (preceded by a '+') before it is executed.
# other flags for debugging: -Canvu, please also consult 'man sh'.

###################
### Description ###
###################
# This script is used by cijenkins.zib.de.
# It builds SCIP with Make and an LPSolver and runs mipdebug and MINLP on different settings.

# Usage: from scip root execute
#        ./weeklyperf_testruns.sh

#############
### Setup ###
#############
mkdir -p lib/include
mkdir -p lib/static
mkdir -p settings
ln -s ../../../SOPLEX_COMP=gnu_OPT=opt_weekly/src lib/include/spxinc
ln -s ../../../SOPLEX_COMP=gnu_OPT=opt_weekly/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.opt.a

ln -s /OPTI/adm_cple/cplex/include/ilcplex lib/include/cpxinc
ln -s /OPTI/adm_cple/cplex/lib/x86-64_linux/static_pic/libcplex.a lib/static/libcplex.linux.x86_64.gnu.a

ln -s /OPTI/adm_cple/ipopt lib/static/ipopt.linux.x86_64.gnu.opt

ln -s /optimi/usr/sw/bliss lib/include/bliss
ln -s /optimi/usr/sw/bliss/libbliss.a lib/static/libbliss.linux.x86_64.gnu.a

###################
### Compilation ###
###################
yes "" | make LPS=spx LPSOPT=opt VERSION=perfspx IPOPT=true SYM=bliss ZIMPL=false COMP=gnu OPT=opt -j4 USRFLAGS=-Werror
yes "" | make LPS=cpx LPSOPT=opt VERSION=perfcpx IPOPT=true SYM=bliss ZIMPL=false COMP=gnu OPT=opt -j4 USRFLAGS=-Werror

######################
### Testrun: setup ###
######################
ln -s /optimi/kombadon/IP check/
ln -s /optimi/kombadon/MINLP check/

##########################
### Testrun: execution ###
##########################
# soplex
for testset in mipdev-solvable; do
    for setting in default; do
        make testcluster TEST=$testset SETTINGS=$setting LPS=spx OPT=opt VERSION=perfspx IPOPT=true SYM=bliss EXCLUSIVE=true TIME=7200 MEM=35000 QUEUE=M620v3-low | PERF=performance VERSION=perfspx check/jenkins_check_results.sh $testset $setting
    done
done

# copy sap settings
cp ~/sap-400-pure.set settings/.

for testset in sapdev-solvable; do
    for setting in sap-400-pure; do
        make testcluster TEST=$testset SETTINGS=$setting LPS=spx OPT=opt VERSION=perfspx IPOPT=true SYM=bliss EXCLUSIVE=true TIME=7200 MEM=35000 QUEUE=M630v2 | PERF=performance VERSION=perfspx check/jenkins_check_results.sh $testset $setting
    done
done

# cplex
echo "limits/gap = 0.0001" > settings/minlp.set
echo "numerics/checkfeastolfac = 10.0" >> settings/minlp.set
for testset in minlpdev-solvable; do
    for setting in minlp; do
        make testcluster TEST=$testset SETTINGS=$setting LPS=cpx OPT=opt VERSION=perfcpx IPOPT=true SYM=bliss EXCLUSIVE=true TIME=3600 MEM=35000 QUEUE=M630-low | PERF=performance VERSION=perfcpx check/jenkins_check_results.sh $testset $setting
    done
done
