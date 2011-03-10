#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            *
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
DISPFREQ=${10}
CONTINUE=${11}
QUEUE=${12}
PPN=${13}
CLIENTTMPDIR=${14}
NOWAITCLUSTER=${15}


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

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "express" queue; this queue has only 4 CPUs
#       available 
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`

EVALFILE=$SCIPPATH/results/check.$QUEUE.$TSTNAME.$BINID.$SETNAME.eval
echo > $EVALFILE

# counter to define file names for a test set uniquely 
COUNT=1

for i in `cat $TSTNAME.test` DONE
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
  if test "$NOWAITCLUSTER" != "1"
  then
      ./waitcluster.sh 1500 $QUEUE 200
  fi

  SHORTFILENAME=`basename $i .gz`
  SHORTFILENAME=`basename $SHORTFILENAME .mps`
  SHORTFILENAME=`basename $SHORTFILENAME .lp`
  SHORTFILENAME=`basename $SHORTFILENAME .opb`

  FILENAME=$USER.$QUEUE.$TSTNAME.$COUNT"_"$SHORTFILENAME.$BINID.$SETNAME
  BASENAME=$SCIPPATH/results/$FILENAME

  TMPFILE=$BASENAME.tmp
  SETFILE=$BASENAME.prm
  
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

  if test -e $SETFILE
  then
      rm -f $SETFILE
  fi
  
  echo > $TMPFILE
  echo ""                              > $TMPFILE
  if test $FEASTOL != "default"
  then
      echo set simplex tolerances feas $FEASTOL    >> $TMPFILE
      echo set mip tolerances integrality $FEASTOL >> $TMPFILE
  fi
  echo set timelimit $TIMELIMIT           >> $TMPFILE
  echo set clocktype 0                    >> $TMPFILE
  echo set mip display 3                  >> $TMPFILE
  echo set mip interval $DISPFREQ         >> $TMPFILE
  echo set mip tolerances mipgap 0.0      >> $TMPFILE
  echo set mip limits nodes $NODELIMIT    >> $TMPFILE
  echo set mip limits treememory $MEMLIMIT >> $TMPFILE
  echo set threads $THREADS               >> $TMPFILE
  echo set parallel 1                     >> $TMPFILE
  echo set mip strategy kappastats 2      >> $TMPFILE
  echo write $SETFILE                     >> $TMPFILE
  echo read $SCIPPATH/$i                  >> $TMPFILE
  echo display problem stats              >> $TMPFILE
  echo optimize                           >> $TMPFILE
  echo display solution quality           >> $TMPFILE
  echo quit                               >> $TMPFILE

  qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N CPLEX$SHORTFILENAME -v SOLVERPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$i,BASENAME=$FILENAME,CLIENTTMPDIR=$CLIENTTMPDIR -q $QUEUE -o /dev/null -e /dev/null runcluster.sh

done

