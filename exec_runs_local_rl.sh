#!/usr/bin/env bash

VERSION=local-rl-64c5928       #git hash c0ca752b9 (SoPlex 82cab95)
VERSION2=local-rl-cb58065fd    #git hash cb58065fd (SoPlex 82cab95)
VERSION3=local-rl-625d18a6f    #git hash 625d18a6f (SoPlex 82cab95)

OUTPUTDIR=results-local-rapidlearning

QUEUE=M640
MEM=35000
TIME=3600

#further settings: rapidlearning-freq-15-full-exp4,rapidlearning-freq-20-full-exp4,rapidlearning-freq-5-full-exp4,rapidlearning-freq-10-full-exp4,default,rapidlearning-freq-5-full,rapidlearning-freq-10-full,rapidlearning-freq-5-exp4-nsolsF,rapidlearning-freq-10-exp4-nsolsF

#### FEASIBILITY
# SETTINGS=rapidlearning-freq-5-exp4-degeneracy,rapidlearning-freq-5-exp4-dualbound,rapidlearning-freq-5-exp4-leaves,rapidlearning-freq-5-exp4-localobj,rapidlearning-freq-5-exp4-sblps,rapidlearning-freq-5-exp4-nsols
# make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-glb-degeneracy,rapidlearning-glb-localobj,rapidlearning-glb-nsols
# make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=feasibility-timo-diss OUTPUTDIR=$OUTPUTDIR testcluster


#### MMM-IP

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy-applybdchgsF,rapidlearning-freq-5-exp4-degeneracy-applyinfervalsF,rapidlearning-freq-5-exp4-degeneracy-applyconflictsF,rapidlearning-freq-5-exp4-degeneracy-applyprimalsolF
# make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy,rapidlearning-freq-5-exp4-dualbound,rapidlearning-freq-5-exp4-leaves,rapidlearning-freq-5-exp4-localobj,rapidlearning-freq-5-exp4-sblps,rapidlearning-freq-5-exp4-nsols
# make VERSION=$VERSION2 GLBSEEDSHIFT=1 SEEDS=3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy-leaves,rapidlearning-freq-5-exp4-degeneracy-leaves-localobj,rapidlearning-freq-5-exp4-leaves-obj
# make VERSION=$VERSION2 GLBSEEDSHIFT=1 SEEDS=3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-glb-degeneracy,rapidlearning-glb-localobj,rapidlearning-glb-nsols
# make VERSION=$VERSION2 GLBSEEDSHIFT=1 SEEDS=3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy-leaves-localobj-applybdchgsF,rapidlearning-freq-5-exp4-degeneracy-leaves-localobj-applyinfervalsF,rapidlearning-freq-5-exp4-degeneracy-leaves-localobj-applyconflictsF,rapidlearning-freq-5-exp4-degeneracy-leaves-localobj-applyprimalsolF
# make VERSION=$VERSION2 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster


### FINAL TEST

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy,rapidlearning-freq-5-exp4-dualbound,rapidlearning-freq-5-exp4-leaves,rapidlearning-freq-5-exp4-localobj,rapidlearning-freq-5-exp4-sblps,rapidlearning-freq-5-exp4-nsols
SETTINGS=rapidlearning-freq-5-exp4-no-checks
make VERSION=$VERSION3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-glb,rapidlearning-glb-full,rapidlearning-glb-degeneracy,rapidlearning-glb-localobj,rapidlearning-glb-nsols
# make VERSION=$VERSION3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy-applybdchgsF,rapidlearning-freq-5-exp4-degeneracy-applyinfervalsF,rapidlearning-freq-5-exp4-degeneracy-applyconflictsF,rapidlearning-freq-5-exp4-degeneracy-applyprimalsolF
# make VERSION=$VERSION3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster

# SETTINGS=rapidlearning-freq-5-exp4-degeneracy-leaves,rapidlearning-freq-5-exp4-degeneracy-leaves-localobj,rapidlearning-freq-5-exp4-leaves-obj
# make VERSION=$VERSION3 CONTINUE=true TIME=$TIME MEM=$MEM QUEUE=$QUEUE EXCLUSIVE=true SETTINGS=$SETTINGS TEST=MMM-IP OUTPUTDIR=$OUTPUTDIR testcluster
