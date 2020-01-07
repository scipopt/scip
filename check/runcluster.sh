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

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    echo Skipping test since the path $CLIENTTMPDIR for the tmp-dir does not exist.
    exit
fi

OUTFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.err
TMPFILE=$SOLVERPATH/$OUTPUTDIR/$BASENAME.tmp

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE

# ensure TMPFILE is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap "
  mv $OUTFILE $SOLVERPATH/$OUTPUTDIR/$BASENAME.out
  mv $ERRFILE $SOLVERPATH/$OUTPUTDIR/$BASENAME.err
  rm -f $TMPFILE
" EXIT

echo "hard time limit: $HARDTIMELIMIT">>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT" >>$OUTFILE
echo @01 $FILENAME ===========      >> $OUTFILE
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE
$EXECNAME                < $TMPFILE >> $OUTFILE 2>>$ERRFILE
date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE
