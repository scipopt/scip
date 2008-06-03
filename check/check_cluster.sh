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
# $Id: check_cluster.sh,v 1.4 2008/06/03 12:14:31 bzfheinz Exp $
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

echo $SCIPPATH


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
HARDMEMLIMIT=`echo "($MEMLIMIT*1.1+10)*1024" | bc`

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

  OUTFILE=$BASENAME.out
  ERRFILE=$BASENAME.err
  TMPFILE=$BASENAME.tmp
  SETFILE=$BASENAME.set
  
  echo $OUTFILE >> $EVALFILE
  
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

  echo $i                                >> $ERRFILE
  date                                   >> $ERRFILE

  export SCIPPATH=$SCIPPATH
  export BINNAME=$BINNAME
  export TMPFILE=$TMPFILE
  export FILENAME=$i

  qsub -l walltime=$HARDTIMELIMIT -N SCIP$SHORTFILENAME -V -o $OUTFILE -e $ERRFILE -q $QUEUE runcluster.sh
done
