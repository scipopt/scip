#!/bin/sh
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2005 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2005 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check.sh,v 1.19 2005/05/31 17:20:06 bzfpfend Exp $
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINID.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

SETTINGS=settings/$SETNAME.set

uname -a >$OUTFILE
uname -a >$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

for i in `cat $TSTNAME.test`
do
    if [ -f $i ]
    then
        echo @01 $i ===========
	echo @01 $i ===========                >> $ERRFILE
	echo set load $SETTINGS                >  $TMPFILE
	if [ $FEASTOL != "default" ]
	then
	    echo set numerics feastol $FEASTOL >> $TMPFILE
	fi
	echo set limits time $TIMELIMIT        >> $TMPFILE
	echo set limits nodes $NODELIMIT       >> $TMPFILE
	echo set limits memory $MEMLIMIT       >> $TMPFILE
	echo set timing clocktype 1            >> $TMPFILE
	echo set display verblevel 3           >> $TMPFILE
	echo set display freq 10000            >> $TMPFILE
	echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
	echo set save $SETFILE                 >> $TMPFILE
	echo read $i                           >> $TMPFILE
	echo optimize                          >> $TMPFILE
	echo display statistics                >> $TMPFILE
	echo display solution                  >> $TMPFILE
	echo checksol                          >> $TMPFILE
	echo quit                              >> $TMPFILE
	echo -----------------------------
	date
	echo -----------------------------
	../$2 < $TMPFILE 2>>$ERRFILE
	echo -----------------------------
	date
	echo -----------------------------
	echo
	echo =ready=
    else
	echo @02 FILE NOT FOUND: $i ===========
	echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
    fi
done | tee -a $OUTFILE

rm $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

gawk -f check.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
