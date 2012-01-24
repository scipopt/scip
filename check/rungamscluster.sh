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

# check if tmp-path exists 
if test ! -d $CLIENTTMPDIR
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
LSTFILE=$CLIENTTMPDIR/$BASENAME.lst
TRCFILE=$CLIENTTMPDIR/$BASENAME.trc
SCRDIR=$CLIENTTMPDIR/$BASENAME.scr

# setup scratch directory
mkdir -p $SCRDIR
GAMSOPTS="$GAMSOPTS SCRDIR=$SCRDIR"

# initialize trace file
echo "* Trace Record Definition" > $TRCFILE
echo "* GamsSolve" >> $TRCFILE
echo "* InputFileName,ModelType,SolverName,OptionFile,Direction,NumberOfEquations,NumberOfVariables,NumberOfDiscreteVariables,NumberOfNonZeros,NumberOfNonlinearNonZeros," >> $TRCFILE
echo "* ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime,ETSolver,NumberOfIterations,NumberOfNodes" >> $TRCFILE

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
echo @01 $FILENAME ===========      >> $OUTFILE 
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE

# run GAMS and check return code
$GAMSBIN $GMSFILE $GAMSOPTS output=$LSTFILE inputdir=$INPUTDIR $MODTYPE=$SOLVER $GDXFILE traceopt=3 trace=$TRCFILE >> $OUTFILE 2>>$ERRFILE
gamsrc=$?
if test $gamsrc != 0
then
  echo "GAMS returned with error code $gamsrc. There was some problem." >>$ERRFILE
  #TODO write 13/13 into trace file
fi

date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE

mv $OUTFILE results/$BASENAME.out
mv $LSTFILE results/$BASENAME.lst
mv $ERRFILE results/$BASENAME.err
mv $TRCFILE results/$BASENAME.trc

rm -r $SCRDIR
