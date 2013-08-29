#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
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
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
DISPFREQ=${10}
CONTINUE=${11}
QUEUETYPE=${12}
QUEUE=${13}
PPN=${14}
CLIENTTMPDIR=${15}
NOWAITCLUSTER=${16}
EXCLUSIVE=${17}
PERMUTE=${18}

# check if all variables defined (by checking the last one)
if test -z $PERMUTE
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAME       = $SETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "QUEUETYPE     = $QUEUETYPE"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "EXCLUSIVE     = $EXCLUSIVE"
    echo "PERMUTE       = $PERMUTE"
    exit 1;
fi

# get current SCIP path
SCIPPATH=`pwd`

if test ! -e $SCIPPATH/results
then
    mkdir $SCIPPATH/results
fi

SETTINGS=$SCIPPATH/../settings/$SETNAME.set

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e $SETTINGS
    then
        echo Skipping test since the settings file $SETTINGS does not exist.
        exit
    fi
fi

# check if binary exists
if test ! -e $SCIPPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
    exit
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "$QUEUETYPE" = "srun"
then
    HOSTMEM=`ulimit -m`
    if test "$HOSTMEM" != "unlimited"
    then
        if [ `expr $HARDMEMLIMIT \* 1024` -gt $HOSTMEM ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit
        fi
    fi
fi

# in case of qsub queue the memory is measured in kB and in case of srun the time needs to be formatted
if test  "$QUEUETYPE" = "qsub"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
else
    MYMINUTES=0
    MYHOURS=0
    MYDAYS=0

    #calculate seconds, minutes, hours and days
    MYSECONDS=`expr $HARDTIMELIMIT % 60`
    TMP=`expr $HARDTIMELIMIT / 60`
    if test "$TMP" != "0"
    then
	MYMINUTES=`expr $TMP % 60`
	TMP=`expr $TMP / 60`
	if test "$TMP" != "0"
	then
	    MYHOURS=`expr $TMP % 24`
	    MYDAYS=`expr $TMP / 24`
	fi
    fi
    #format seconds to have two characters
    if test ${MYSECONDS} -lt 10
    then
	MYSECONDS=0${MYSECONDS}
    fi
    #format minutes to have two characters
    if test ${MYMINUTES} -lt 10
    then
	MYMINUTES=0${MYMINUTES}
    fi
    #format hours to have two characters
    if test ${MYHOURS} -lt 10
    then
	MYHOURS=0${MYHOURS}
    fi
    #format HARDTIMELIMT
    if test ${MYDAYS} = "0"
    then
	HARDTIMELIMIT=${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    else
	HARDTIMELIMIT=${MYDAYS}-${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    fi
fi

#define clusterqueue, which might not be the QUEUE, cause this might be an alias for a bunch of QUEUEs
CLUSTERQUEUE=$QUEUE

NICE=""
ACCOUNT="mip"

if test $CLUSTERQUEUE = "dbg"
then
    CLUSTERQUEUE="mip-dbg,telecom-dbg"
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "telecom-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"
fi

# check if the slurm blades should be used exclusively
if test "$EXCLUSIVE" = "true"
then
    EXCLUSIVE=" --exclusive"
    if test $CLUSTERQUEUE = "opt"
    then
        CLUSTERQUEUE="M610"
    fi
else
    EXCLUSIVE=""
fi


# counter to define file names for a test set uniquely
COUNT=0

# loop over permutations
for ((p = 0; $p <= $PERMUTE; p++))
do
    # if number of permutations is positive, add postfix
    if test $PERMUTE -gt 0
    then
	EVALFILE=$SCIPPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"#p"$p.eval
    else
	EVALFILE=$SCIPPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
    fi
    echo > $EVALFILE

    # loop over testset
    for i in `cat testset/$TSTNAME.test` DONE
    do
	if test "$i" = "DONE"
	then
	    break
	fi

        # increase the index for the instance tried to solve, even if the filename does not exist
	COUNT=`expr $COUNT + 1`

        # check if problem instance exists
	if test -f $SCIPPATH/$i
	then

            # the cluster queue has an upper bound of 2000 jobs; if this limit is
            # reached the submitted jobs are dumped; to avoid that we check the total
            # load of the cluster and wait until it is save (total load not more than
            # 1600 jobs) to submit the next job.
	    if test "$NOWAITCLUSTER" != "1"
	    then
		if test  "$QUEUETYPE" != "qsub"
		then
		    echo "waitcluster does not work on slurm cluster"
		fi
		./waitcluster.sh 1600 $QUEUE 200
	    fi

	    SHORTPROBNAME=`basename $i .gz`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .mps`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .lp`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .opb`

	    # if number of permutations is positive, add postfix
	    if test $PERMUTE -gt 0
	    then
		FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME#"p"$p
	    else
		FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME
	    fi

	    BASENAME=$SCIPPATH/results/$FILENAME

	    TMPFILE=$BASENAME.tmp
	    SETFILE=$BASENAME.set

	    echo $BASENAME >> $EVALFILE

            # in case we want to continue we check if the job was already performed
	    if test "$CONTINUE" != "false"
	    then
		if test -e results/$FILENAME.out
		then
		    echo skipping file $i due to existing output file $FILENAME.out
		    continue
		fi
	    fi

	    echo > $TMPFILE
	    if test $SETNAME != "default"
	    then
		echo set load $SETTINGS            >>  $TMPFILE
	    fi
	    if test $FEASTOL != "default"
	    then
		echo set numerics feastol $FEASTOL >> $TMPFILE
	    fi

	    # if permutation counter is positive add permutation seed (0 = default)
	    if test $p -gt 0
	    then
		echo set misc permutationseed $p   >> $TMPFILE
	    fi

	    echo set limits time $TIMELIMIT        >> $TMPFILE
	    echo set limits nodes $NODELIMIT       >> $TMPFILE
	    echo set limits memory $MEMLIMIT       >> $TMPFILE
	    echo set lp advanced threads $THREADS  >> $TMPFILE
	    echo set timing clocktype 1            >> $TMPFILE
	    echo set display freq $DISPFREQ        >> $TMPFILE
	    echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
	    if test "$LPS" = "none"
	    then
		echo set lp solvefreq -1           >> $TMPFILE # avoid solving LPs in case of LPS=none
	    fi
	    echo set save $SETFILE                 >> $TMPFILE
	    echo read $SCIPPATH/$i                 >> $TMPFILE
#           echo presolve                          >> $TMPFILE
	    echo optimize                          >> $TMPFILE
	    echo display statistics                >> $TMPFILE
#           echo display solution                  >> $TMPFILE
	    echo checksol                          >> $TMPFILE
	    echo quit                              >> $TMPFILE

            # check queue type
	    if test  "$QUEUETYPE" = "srun"
	    then
                # additional environment variables needed by runcluster.sh
		export SOLVERPATH=$SCIPPATH
		export EXECNAME=$SCIPPATH/../$BINNAME
		export BASENAME=$FILENAME
		export FILENAME=$i
		export CLIENTTMPDIR=$CLIENTTMPDIR
		sbatch --job-name=SCIP$SHORTPROBNAME --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT $NICE --time=${HARDTIMELIMIT} ${EXCLUSIVE} --output=/dev/null runcluster.sh
	    else
                # -V to copy all environment variables
		qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N SCIP$SHORTPROBNAME -v SOLVERPATH=$SCIPPATH,EXECNAME=$SCIPPATH/../$BINNAME,BASENAME=$FILENAME,FILENAME=$i,CLIENTTMPDIR=$CLIENTTMPDIR -V -q $CLUSTERQUEUE -o /dev/null -e /dev/null runcluster.sh
	    fi
	else
	    echo "input file "$SCIPPATH/$i" not found!"
	fi
    done
done
