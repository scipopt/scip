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

export LANG=C

REMOVE=0
AWKARGS=""
FILES=""

for i in $@
do
  if test ! -e $i
  then
      if test "$i" = "-r"
      then
          REMOVE=1
      else
          AWKARGS="$AWKARGS $i"
      fi
  else	
      FILES="$FILES $i"
  fi
done

for FILE in $FILES
do
 
  DIR=`dirname $FILE`
  EVALFILE=`basename $FILE .eval`
  EVALFILE=`basename $EVALFILE .out`

  OUTFILE=$DIR/$EVALFILE.out 
  ERRFILE=$DIR/$EVALFILE.err
  TRCFILE=$DIR/$EVALFILE.trc
  LSTFILE=$DIR/$EVALFILE.lst
  RESFILE=$DIR/$EVALFILE.res
  TEXFILE=$DIR/$EVALFILE.tex
  PAVFILE=$DIR/$EVALFILE.pav
  
  # check if the eval file exists; if this is the case construct the overall solution files
  if test -e $DIR/$EVALFILE.eval
  then
    # in case an output file exists, copy it away to save the results
    DATEINT=`date +"%s"`
    if test -e $OUTFILE
    then
        cp $OUTFILE $OUTFILE.old-$DATEINT
    fi
    if test -e $ERRFILE
    then
        cp $ERRFILE $ERRFILE.old-$DATEINT
    fi
    if test -e $TRCFILE
    then
        cp $TRCFILE $TRCFILE.old-$DATEINT
    fi
    if test -e $LSTFILE
    then
        cp $LSTFILE $LSTFILE.old-$DATEINT
    fi

    echo > $OUTFILE
    echo > $ERRFILE
    # initialize gams trace file
    echo "* Trace Record Definition" > $TRCFILE
    echo "* GamsSolve" >> $TRCFILE
    echo "* InputFileName,ModelType,SolverName,OptionFile,Direction,NumberOfEquations,NumberOfVariables,NumberOfDiscreteVariables,NumberOfNonZeros,NumberOfNonlinearNonZeros," >> $TRCFILE
    echo "* ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime,ETSolver,NumberOfIterations,NumberOfNodes" >> $TRCFILE
    echo "*" >> $TRCFILE

    echo "create overall output, error, and trace file for $EVALFILE"

    for i in `cat $DIR/$EVALFILE.eval` DONE
    do
      if test "$i" = "DONE"
      then
        break
      fi
     
      # pass auxiliary lines about solver and limits form eval file to trace file
      case $i in *,* )
        echo "* $i" >> $TRCFILE
        continue
        ;;
      esac

      FILE=$i.out
      if test -e $FILE
      then
        cat $FILE >> $OUTFILE
        if test "$REMOVE" = "1"
        then
          rm -f $FILE
        fi
      else
        echo Missing $i
      fi
      
      FILE=$i.err
      if test -e $FILE
      then
        cat $FILE >> $ERRFILE
        if test "$REMOVE" = "1"
        then
          rm -f $FILE
        fi
      fi
      
      FILE=$i.trc
      if test -e $FILE
      then
        grep -v "^*" $FILE >> $TRCFILE
        if test "$REMOVE" = "1"
        then
          rm -f $FILE
        fi
      else
        echo Missing $i
      fi
      
      FILE=$i.lst
      if test -e $FILE
      then
        cat $FILE >> $LSTFILE
        if test "$REMOVE" = "1"
        then
          rm -f $FILE
        fi
      fi
    done
      
    if test "$REMOVE" = "1"
    then
      rm -f $DIR/$EVALFILE.eval
    fi
  fi

  # check if the out file exists
  if test -e $DIR/$EVALFILE.out
  then
    echo create results for $EVALFILE

    # detect test set
    TSTNAME=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`
    echo "Testset " $TSTNAME

    if test -f testset/$TSTNAME.test
    then
      TESTFILE=testset/$TSTNAME.test
    else
      TESTFILE=""
    fi

    if test -f testset/$TSTNAME.solu
    then
      SOLUFILE=testset/$TSTNAME.solu
    else 
      if test -f testset/all.solu
      then
        SOLUFILE=testset/all.solu
      else
        SOLUFILE=""
      fi
    fi

    awk -f check_gams.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $SOLUFILE $TRCFILE | tee $RESFILE
  fi
done
