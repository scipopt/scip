#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# Call with "make testcluster"
#
# The queue is passed via $QUEUE (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call "./evalcheck_cluster.sh
# results/check.$TSTNAME.$BINNAME.$SETNAME.eval in directory check/
# This leads to result files
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.out
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.res
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAMES=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
LPS=${10}
DISPFREQ=${11}
CONTINUE=${12}
QUEUETYPE=${13}
QUEUE=${14}
PPN=${15}
CLIENTTMPDIR=${16}
NOWAITCLUSTER=${17}
EXCLUSIVE=${18}
PERMUTE=${19}
SEEDS=${20}
DEBUGTOOL=${21}
REOPT=${22}
OPTCOMMAND=${23}
SETCUTOFF=${24}
VISUALIZE=${25}

SOLVER=fscip

# check if all variables defined (by checking the last one)
if test -z $VISUALIZE
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAMES      = $SETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "LPS           = $LPS"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "QUEUETYPE     = $QUEUETYPE"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "EXCLUSIVE     = $EXCLUSIVE"
    echo "PERMUTE       = $PERMUTE"
    echo "SEEDS         = $SEEDS"
    echo "DEBUGTOOL      = $DEBUGTOOL"
    echo "REOPT         = $REOPT"
    echo "OPTCOMMAND    = $OPTCOMMAND"
    echo "SETCUTOFF     = $SETCUTOFF"
    echo "VISUALIZE     = $VISUALIZE"
    exit 1;
fi

# configure cluster-related environment variables
. ./configuration_cluster.sh $QUEUE $PPN $EXCLUSIVE $QUEUETYPE

# the srun queue requires a format duration HH:MM:SS (and optionally days),
# whereas the qsub requires the memory limit in kB
if test "$QUEUETYPE" != "qsub"
then
    TIMEFORMAT="format"
    MEMFORMAT="MB"
else
    TIMEFORMAT="sec"
    MEMFORMAT="B"
fi
# call routines for creating the result directory, checking for existence
# of passed settings, etc
. ./configuration_set_fscip.sh $BINNAME $TSTNAME $SETNAMES $TIMELIMIT $TIMEFORMAT $MEMLIMIT $MEMFORMAT $DEBUGTOOL $SETCUTOFF


# at the first time, some files need to be initialized. set to "" after the innermost loop
# finished the first time
INIT="true"

# counter to define file names for a test set uniquely
COUNT=0
# loop over permutations
# loop over testset
for idx in ${!INSTANCELIST[@]}
do
    # retrieve instance and timelimits from arrays set in the configuration_set.sh script
    INSTANCE=${INSTANCELIST[$idx]}
    TIMELIMIT=${TIMELIMLIST[$idx]}
    HARDTIMELIMIT=${HARDTIMELIMLIST[$idx]}
    # increase the index for the instance tried to solve, even if the filename does not exist
    COUNT=`expr $COUNT + 1`

    # we need the DONE keyword for the check.sh script to automatically run evalcheck, here it is not needed
    if test "$INSTANCE" = "DONE"
    then
        continue
    fi

    # run different random seeds
    for ((s = 0; $s <= $SEEDS; s++))
    do
        # permute transformed problem
	for ((p = 0; $p <= $PERMUTE; p++))
	do
	    # the cluster queue has an upper bound of 2000 jobs; if this limit is
	    # reached the submitted jobs are dumped; to avoid that we check the total
	    # load of the cluster and wait until it is save (total load not more than
	    # 1600 jobs) to submit the next job.
	    if test "${NOWAITCLUSTER}" -eq "0" && test "$QUEUETYPE" = "qsub"
	    then
		./waitcluster.sh 1600 $QUEUE 200
	    elif test "${NOWAITCLUSTER}" -eq "0"
	    then
		echo "waitcluster does not work on slurm cluster"
	    fi
	    # loop over settings
	    for SETNAME in ${SETTINGSLIST[@]}
	    do
		# infer the names of all involved files from the arguments
		. ./configuration_logfiles_fscip.sh $INIT $COUNT $INSTANCE $BINID $PERMUTE $SEEDS $SETNAME $TSTNAME $CONTINUE $QUEUE $THREADS $p $s

		# skip instance if log file is present and we want to continue a previously launched test run
		if test "$SKIPINSTANCE" = "true"
		then
		    continue
		fi

		# find out the solver that should be used
		SOLVER=`stripversion $BINNAME`

		JOBNAME="`capitalize ${SOLVER}`${SHORTPROBNAME}"

		export EXECNAME=$SCIPPATH/../bin/$BINNAME

		# check queue type
		if test  "$QUEUETYPE" = "srun"
		then
		# additional environment variables needed by run.sh
		    export SOLVERPATH=$SCIPPATH
		    export BASENAME=$FILENAME
		    export FILENAME=$INSTANCE
		    export CLIENTTMPDIR
		    export HARDTIMELIMIT
		    export HARDMEMLIMIT
		    export CHECKERPATH=$SCIPPATH/solchecker
		    export SETFILE
                    export SETNAME
                    export THREADS
		    export TIMELIMIT
		    # the space at the end is necessary
		    export SRUN="srun --cpu_bind=verbose,cores -v -v "
		    sbatch --ntasks=1 --cpus-per-task=`expr $THREADS + 1` --job-name=${JOBNAME} --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT $NICE --time=${HARDTIMELIMIT} --cpu-freq=highm1 ${EXCLUSIVE} --output=/dev/null run_fscip.sh
		else
		    # -V to copy all environment variables
		    qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N ${JOBNAME} \
			-v SOLVERPATH=$SCIPPATH,EXECNAME=${EXECNAME},BASENAME=$FILENAME,FILENAME=$INSTANCE,CLIENTTMPDIR=$CLIENTTMPDIR \
			-V -q $CLUSTERQUEUE -o /dev/null -e /dev/null run_fscip.sh
		fi
	    done # end for SETNAME
	done # end for PERMUTE
    done #end for SEEDS

    # after the first termination of the set loop, no file needs to be initialized anymore
    INIT="false"
done # end for TSTNAME
