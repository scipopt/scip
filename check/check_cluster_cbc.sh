#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de       *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# Call with "make testclustercbc"
#
# The queue is passed via $QUEUE (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way 
# that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
# of each computer. If course, the value depends on the specific computer/queue.
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
BINID=$BINNAME.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
CONTINUE=${10}
QUEUE=${11}
PPN=${12}

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
        echo skipping test due to not existes of the settings file $SETTINGS
        exit
    fi
fi

# we add 10% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "exrpess" queue; this queue has only 4 CPUs
#       available 
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# since bash counts cpu time we need the time limit for each thread
HARDTIMELIMIT=`expr $HARDTIMELIMIT \* $THREADS`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`

EVALFILE=$SCIPPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
echo > $EVALFILE

for i in `cat testset/$TSTNAME.test` DONE
do
  if test "$i" = "DONE"
      then
      break
  fi

  echo adding instance $COUNT to queue

  # the cluster queue has an upper bound of 2000 jobs; if this limit is
  # reached the submitted jobs are dumped; to avoid that we check the total
  # load of the cluster and wait until it is save (total load not more than
  # 1900 jobs) to submit the next job.
  ./waitcluster.sh 1500 $QUEUE 200

  SHORTFILENAME=`basename $i .gz`
  SHORTFILENAME=`basename $SHORTFILENAME .mps`
  SHORTFILENAME=`basename $SHORTFILENAME .lp`
  SHORTFILENAME=`basename $SHORTFILENAME .opb`

  FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTFILENAME.$QUEUE.$BINID.$SETNAME
  BASENAME=$SCIPPATH/results/$FILENAME

  TMPFILE=$BASENAME.tmp
  SETFILE=$BASENAME.set
  
  echo $BASENAME >> $EVALFILE

  COUNT=`expr $COUNT + 1`

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
  if test $FEASTOL != "default"
      then
      echo primalTolerance $FEASTOL       >> $TMPFILE
      echo integerTolerance $FEASTOL      >> $TMPFILE
  fi
#workaround: since CBC only looks at cpu-time, we multiply the timelimit with the number of threads
  TIMELIMIT=`expr $TIMELIMIT \* $THREADS`
  echo seconds $TIMELIMIT                 >> $TMPFILE
  echo threads $THREADS                   >> $TMPFILE
  echo ratioGap 0.0                       >> $TMPFILE
  echo maxNodes $NODELIMIT                >> $TMPFILE
  echo import $SCIPPATH/$i                >> $TMPFILE
  echo ratioGap                           >> $TMPFILE
  echo allowableGap                       >> $TMPFILE
  echo seconds                            >> $TMPFILE
  echo stat                               >> $TMPFILE
  echo solve                              >> $TMPFILE
  echo quit                               >> $TMPFILE
 
 
  qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N CBC$SHORTFILENAME -v SOLVERPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$i,BASENAME=$FILENAME -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
 
done
