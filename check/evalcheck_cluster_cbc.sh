#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
export LANG=C

AWKARGS=""
DIR=`dirname $1`
EVALFILE=`basename $1 .eval`

OUTFILE=$DIR/$EVALFILE.out 
ERRFILE=$DIR/$EVALFILE.err
RESFILE=$DIR/$EVALFILE.res
TEXFILE=$DIR/$EVALFILE.tex
PAVFILE=$DIR/$EVALFILE.pav

echo > $OUTFILE
echo > $ERRFILE
echo create overall output and error file
for i in `cat $DIR/$EVALFILE.eval` DONE
  do
  if test "$i" = "DONE"
      then
      break
  fi

  FILE=$i.out
  if test -e $FILE
      then
      cat $FILE >> $OUTFILE
  fi

  FILE=$i.err
  if test -e $FILE
      then
      cat $FILE >> $ERRFILE
  fi
done


TSTNAME=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`

# look for solufiles under the name of the test, the name of the test with everything after the first "_" stripped, and all
SOLUFILE=""
for f in $TSTNAME ${TSTNAME%%_*} ${TSTNAME%%-*} all
do
  if test -f testset/${f}.solu
  then
      SOLUFILE=testset/${f}.solu
      break
  fi
done

awk -f check_cbc.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $SOLUFILE $OUTFILE | tee $RESFILE


