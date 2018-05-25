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

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-04
# absolut tolerance for checking integrality constraints
INTTOL=1e-04

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR/${USER}-tmpdir
then
    mkdir $CLIENTTMPDIR/${USER}-tmpdir
    echo Creating directory $CLIENTTMPDIR/${USER}-tmpdir for temporary outfile
fi

OUTFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.err
SOLFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.sol
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
echo "hard time limit: $HARDTIMELIMIT">>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT" >>$OUTFILE
echo @01 $FILENAME ===========      >> $OUTFILE
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE


SETTINGSSTR=""
if test $SETTINGS != ""
then
    SETTINGSSTR="-s "$SETTINGS
fi

#if we use a debugger command, we need to replace the errfile place holder by the actual err-file for logging
EXECNAME=${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}
echo $EXECNAME $FILENAME >> $ERRFILE

echo $EXECNAME $FILENAME $SETTINGSSTR -t $TIMELIMIT -m $MEMLIMIT -d $DISPFREQ 2>>$ERRFILE | tee -a $OUTFILE
bash -c "$EXECNAME $FILENAME $SETTINGSSTR -t $TIMELIMIT -m $MEMLIMIT -d $DISPFREQ" 2>>$ERRFILE | tee -a $OUTFILE
retcode=${PIPESTATUS[0]}
if test $retcode != 0
then
  echo "$EXECNAME returned with error code $retcode." >>$ERRFILE
fi

if test -e $SOLFILE
then
    # translate SCIP solution format into format for solution checker. The
    # SOLFILE format is a very simple format where in each line we have a
    # <variable, value> pair, separated by spaces.  A variable name of
    # =obj= is used to store the objective value of the solution, as
    # computed by the solver. A variable name of =infeas= can be used to
    # indicate that an instance is infeasible.
    sed ' /solution status:/d;
	    s/objective value:/=obj=/g;
	    s/no solution available//g' $SOLFILE > $TMPFILE
    mv $TMPFILE $SOLFILE

    # check if the link to the solution checker exists
    if test -f "$CHECKERPATH/bin/solchecker"
    then
      echo
      $SHELL -c " $CHECKERPATH/bin/solchecker $FILENAME $SOLFILE $LINTOL $INTTOL" 2>>$ERRFILE | tee -a $OUTFILE
      echo
    fi
fi

date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE

mv $OUTFILE $SOLVERPATH/results/$BASENAME.out
mv $ERRFILE $SOLVERPATH/results/$BASENAME.err

rm -f $TMPFILE
rm -f $SOLFILE
#chmod g+r $ERRFILE
#chmod g+r $SCIPPATH/results/$BASENAME.out
#chmod g+r $SCIPPATH/results/$BASENAME.set
