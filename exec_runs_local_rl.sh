#!/usr/bin/env bash

VERSION=local-rl-64c5928    #git hash c0ca752b9 (SoPlex 82cab95)
VERSION2=local-rl-cb58065fd    #git hash cb58065fd (SoPlex 82cab95)

OUTPUTDIR=results-local-rapidlearning

QUEUE=M640
MEM=35000
TIME=3600

#further settings: rapidlearning-freq-15-full-exp4,rapidlearning-freq-20-full-exp4,rapidlearning-freq-5-full-exp4,rapidlearning-freq-10-full-exp4,default,rapidlearning-freq-5-full,rapidlearning-freq-10-full,rapidlearning-freq-5-exp4-nsolsF,rapidlearning-freq-10-exp4-nsolsF


#### FEASIBILITY
# SETTINGS=rapidlearning-freq-5-exp4-no-checks,rapidlearning-freq-10-exp4-no-checks,rapidlearning-freq-15-exp4-no-checks
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-degeneracy,rapidlearning-freq-10-exp4-degeneracy,rapidlearning-freq-15-exp4-degeneracy
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-dualbound,rapidlearning-freq-10-exp4-dualbound,rapidlearning-freq-15-exp4-dualbound
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-leaves,rapidlearning-freq-10-exp4-leaves,rapidlearning-freq-15-exp4-leaves
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-localobj,rapidlearning-freq-10-exp4-localobj,rapidlearning-freq-15-exp4-localobj
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-sblps,rapidlearning-freq-15-exp4-sblps,rapidlearning-freq-10-exp4-sblps
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-nsols,rapidlearning-freq-10-exp4-nsols,rapidlearning-freq-15-exp4-nsols
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-glb
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster

SETTINGS=rapidlearning-freq-5-exp4-full,rapidlearning-freq-10-exp4-full,rapidlearning-freq-15-exp4-full
make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster

SETTINGS=rapidlearning-glb-full
make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster

#
# #### MMM-IP
# SETTINGS=rapidlearning-freq-5-exp4-no-checks,rapidlearning-freq-10-exp4-no-checks,rapidlearning-freq-15-exp4-no-checks
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-degeneracy,rapidlearning-freq-10-exp4-degeneracy,rapidlearning-freq-15-exp4-degeneracy
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-dualbound,rapidlearning-freq-10-exp4-dualbound,rapidlearning-freq-15-exp4-dualbound
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-leaves,rapidlearning-freq-10-exp4-leaves,rapidlearning-freq-15-exp4-leaves
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-localobj,rapidlearning-freq-10-exp4-localobj,rapidlearning-freq-15-exp4-localobj
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-sblps,rapidlearning-freq-15-exp4-sblps,rapidlearning-freq-10-exp4-sblps
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-freq-5-exp4-nsols,rapidlearning-freq-10-exp4-nsols,rapidlearning-freq-15-exp4-nsols
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
#
# SETTINGS=rapidlearning-glb
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

SETTINGS=rapidlearning-freq-5-exp4-full,rapidlearning-freq-10-exp4-full,rapidlearning-freq-15-exp4-full
make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

SETTINGS=rapidlearning-glb-full
make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster


#### MMMC
# make VERSION=$VERSION CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMMc OUTPUTDIR=$OUTPUTDIR testcluster



