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
# $Id: check_cluster.sh,v 1.20 2009/01/20 13:31:54 bzfheinz Exp $
#
# Call with make testcluster
# Cluster nodes have 8 cores (queue "ib") and 16 GB RAM
# If no time is measured, change to PPN=1 (see below) in order to allow parallel runs
# For more information, see "http://www.zib.de/cluster-user/view/Main/Hardware"
#
# To get the result files call  "./evalcheck_cluster.sh results/check.$TSTNAME.$BINNMAE.$SETNAME.eval 
# in directory check/
# This leads to result files 
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

# number of needed core at a certain cluster node
#  - PPN=8 means we need all core at a node, therefore time measuring is possible
#  - PPN=1 means we need one core at a node, therefore time measuring is not possible
PPN=8

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

# get current SCIP path
SCIPPATH=`pwd`

SETDIR=../settings

# choose a queue for the cluster run 
if test "$OPT" = "opt"
    then
#   QUEUE="gbe"
    QUEUE="ib"
else
    QUEUE="ib"
fi

if test ! -e results
then
    mkdir results
fi

SETTINGS=$SETDIR/$SETNAME.set

# we add 10% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "exrpess" queue; this queue has only 4 CPUs
#       available 
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + \`expr $TIMELIMIT / 10\``

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`

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

  FILENAME=check.$TSTNAME.$DIR"_"$SHORTFILENAME.$BINID.$SETNAME
  BASENAME=$SCIPPATH/results/$FILENAME

  TMPFILE=$BASENAME.tmp
  SETFILE=$BASENAME.set
  
  echo $BASENAME >> $EVALFILE
  
  echo > $TMPFILE
  if test $SETTINGS != "default"
      then
      echo set load $SETTINGS            >>  $TMPFILE
  fi
  if test $FEASTOL != "default"
      then
      echo set numerics feastol $FEASTOL >> $TMPFILE
  fi
  echo set limits time $TIMELIMIT        >> $TMPFILE
  echo set limits nodes $NODELIMIT       >> $TMPFILE
  echo set limits memory $MEMLIMIT       >> $TMPFILE
  echo set timing clocktype 1            >> $TMPFILE
  echo set display verblevel 4           >> $TMPFILE
  echo set display freq $DISPFREQ        >> $TMPFILE
  echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
  echo set save $SETFILE                 >> $TMPFILE
  echo read /workbig/$i                  >> $TMPFILE
  echo optimize                          >> $TMPFILE
  echo display statistics                >> $TMPFILE
#	    echo display solution                  >> $TMPFILE
  echo checksol                          >> $TMPFILE
  echo quit                              >> $TMPFILE

  qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N SCIP$SHORTFILENAME -v SCIPPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$i,BASENAME=$FILENAME -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
done
