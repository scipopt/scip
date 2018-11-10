#!/usr/bin/env bash

VERSION=local-rl
OUTPUTDIR=results-local-rapidlearning

QUEUE=M640
MEM=35000
TIME=3600

SETTINGS=rapidlearning-freq-5-exp4-nsolsF,rapidlearning-freq-10-exp4-nsolsF #rapidlearning-freq-15-full-exp4,rapidlearning-freq-20-full-exp4,rapidlearning-freq-5-full-exp4,rapidlearning-freq-10-full-exp4,default,rapidlearning-glb,rapidlearning-freq-5-full,rapidlearning-freq-10-full

# FEASIBILITY
make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster

# MMM-IP
make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# MMMC
#make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMMc OUTPUTDIR=$OUTPUTDIR testcluster
