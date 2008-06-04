#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check_cluster.sh,v 1.11 2008/06/04 17:03:34 bzfheinz Exp $
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8
DISPFREQ=$9
CONTINUE=${10}
LOCK=${11}
VERSION=${12}
OPT=${13}

# get cuurent SCIP path
SCIPPATH=`pwd`

SETDIR=../settings

# choose a queue for the cluster run 
if test "$OPT" == "opt"
    then
    QUEUE="gbe"
else
    QUEUE="ib"
fi

if test ! -e results
then
    mkdir results
fi

SETTINGS=$SETDIR/$SETNAME.set

# the jobs should have a hard running time of more than 5 minutes; if not so, these
# jobs get automatically assigned in the "exrpess" queue; this queue has only 4 CPU
# available 
HARDTIMELIMIT=`echo "($TIMELIMIT*1.1)+600" | bc` 
HARDMEMLIMIT=`echo "($MEMLIMIT+100)*1024" | bc`

USRPATH=`pwd`

EVALFILE=$USRPATH/results/check.$TSTNAME.$BINID.$SETNAME.eval
echo > $EVALFILE

for i in `cat $TSTNAME.test` DONE
do
  if test "$i" == "DONE"
      then
      break
  fi

  SHORTFILENAME=`basename $i .gz`
  SHORTFILENAME=`basename $SHORTFILENAME .mps`
  SHORTFILENAME=`basename $SHORTFILENAME .lp`

  DIR=`dirname $i`
  DIR=$(echo $DIR|sed 's/\//_/g')

  BASENAME=$USRPATH/results/check.$TSTNAME.$DIR"_"$SHORTFILENAME.$BINID.$SETNAME

  TMPFILE=$BASENAME.tmp
  ERRFILE=$BASENAME.err
  SETFILE=$BASENAME.set
  
  echo $BASENAME.out >> $EVALFILE
  
  echo > $TMPFILE
  if test $SETTINGS != "default"
      then
      echo set load $SETTINGS                >>  $TMPFILE
  fi
  if test $FEASTOL != "default"
      then
      echo set numerics feastol $FEASTOL    >> $TMPFILE
  fi
  echo set limits time $TIMELIMIT        >> $TMPFILE
  echo set limits nodes $NODELIMIT       >> $TMPFILE
  echo set limits memory $MEMLIMIT       >> $TMPFILE
  echo set timing clocktype 1            >> $TMPFILE
  echo set display verblevel 4           >> $TMPFILE
  echo set display freq $DISPFREQ        >> $TMPFILE
  echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
  echo set save $SETFILE                 >> $TMPFILE
  echo read $i                           >> $TMPFILE
  echo optimize                          >> $TMPFILE
  echo display statistics                >> $TMPFILE
#	    echo display solution                  >> $TMPFILE
  echo checksol                          >> $TMPFILE
  echo quit                              >> $TMPFILE

  echo $i                                > $ERRFILE
  date                                   >> $ERRFILE

  qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -N SCIP$SHORTFILENAME -v SCIPPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$i,BASENAME=$BASENAME -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
done
