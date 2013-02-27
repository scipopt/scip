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
WORKDIR=$CLIENTTMPDIR/$BASENAME.scr

# setup scratch directory
mkdir -p $WORKDIR
GAMSOPTS="$GAMSOPTS curdir=$WORKDIR"

# ensure scratch directory is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap "
  rm -r $WORKDIR;
  test -e $OUTFILE && mv $OUTFILE results/$BASENAME.out
  test -e $LSTFILE && mv $LSTFILE results/$BASENAME.lst
  test -e $ERRFILE && mv $ERRFILE results/$BASENAME.err
  test -e $TRCFILE && mv $TRCFILE results/$BASENAME.trc
" EXIT

# initialize trace file
echo "* Trace Record Definition" > $TRCFILE
echo "* GamsSolve" >> $TRCFILE
echo "* InputFileName,ModelType,SolverName,OptionFile,Direction,NumberOfEquations,NumberOfVariables,NumberOfDiscreteVariables,NumberOfNonZeros,NumberOfNonlinearNonZeros," >> $TRCFILE
echo "* ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime,ETSolver,NumberOfIterations,NumberOfNodes" >> $TRCFILE

# setup gams file that sets cutoff
# this only works for models that include %gams.u1% and where the model name is m (e.g., MINLPLib instances)
if test -n "$CUTOFF"
then
  echo "m.cutoff = $CUTOFF;" > $WORKDIR/include.u1
  GAMSOPTS="$GAMSOPTS u1=$WORKDIR/include.u1"
fi

# add commands to .u1 file to read start solutions from available gdx files, if any
if test "$PASSSTARTSOL" = 1 ; then
  for sol in $INPUTDIR/${GMSFILE/%.gms/}*.gdx ;
  do
    if test -e $sol ; then
      # create .u1 file if not existing yet and add u1 command to GAMS options
      if test ! -e $WORKDIR/include.u1 ;
      then
        touch $WORKDIR/include.u1
        GAMSOPTS="$GAMSOPTS u1=$WORKDIR/include.u1"
      fi
      echo "execute_loadpoint '$sol'" >> $WORKDIR/include.u1
    fi
  done
fi

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
