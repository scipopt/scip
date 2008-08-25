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
# $Id: check_cluster_cbc.sh,v 1.2 2008/08/25 12:25:47 bzfwanie Exp $

TSTNAME=$1
BINID=$2
SETNAME=$3
BINNAME="check/$BINID"
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8
MIPGAP=$9
CONTINUE=${10}

# get cuurent SCIP path
SCIPPATH=`pwd`

SETDIR=../settings

# choose a queue for the cluster run 
if test "$OPT" = "opt"
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
HARDMEMLIMIT=`echo "($MEMLIMIT+100)*1024*1000" | bc`


EVALFILE=$SCIPPATH/results/check.$TSTNAME.$BINID.$SETNAME.eval
echo > $EVALFILE

for i in `cat $TSTNAME.test` DONE
do
  if test "$i" = "DONE"
      then
      break
  fi

  SHORTFILENAME=`basename $i .gz`
  SHORTFILENAME=`basename $SHORTFILENAME .mps`
  SHORTFILENAME=`basename $SHORTFILENAME .lp`

  DIR=`dirname $i`
  DIR=$(echo $DIR|sed 's/\//_/g')

  BASENAME=$SCIPPATH/results/check.$TSTNAME.$DIR"_"$SHORTFILENAME.$BINID.$4.$SETNAME

  TMPFILE=$BASENAME.tmp
  ERRFILE=$BASENAME.err
  SETFILE=$BASENAME.set
  
  echo $BASENAME.out >> $EVALFILE

  echo > $TMPFILE
  if test $FEASTOL != "default"
      then
      echo primalTolerance $FEASTOL       >> $TMPFILE
      echo integerTolerance $FEASTOL      >> $TMPFILE
  fi
  echo seconds $TIMELIMIT                 >> $TMPFILE
  if test $MIPGAP != "default"
      then
      echo ratioGap $MIPGAP               >> $TMPFILE
  fi
  echo maxNodes $NODELIMIT                >> $TMPFILE
  echo import $i                          >> $TMPFILE
  echo ratioGap                           >> $TMPFILE
  echo allowableGap                       >> $TMPFILE
  echo seconds                            >> $TMPFILE
  echo stat                               >> $TMPFILE
  echo solve                              >> $TMPFILE
  echo quit                                >> $TMPFILE
  echo $i                                 > $ERRFILE
  date                                    >> $ERRFILE
  
  qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -N CBC$SHORTFILENAME -v SCIPPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$i,BASENAME=$BASENAME -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
done
