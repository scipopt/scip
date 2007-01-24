#!/bin/sh
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2006 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2006 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: check.sh,v 1.38 2007/01/24 10:17:23 bzfpfend Exp $
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

if [ ! -e results ]
then
    mkdir results
fi
if [ ! -e settings ]
then
    mkdir settings
fi
if [ ! -e locks ]
then
    mkdir locks
fi

LOCKFILE=locks/$TSTNAME.$SETNAME.$VERSION.lock
RUNFILE=locks/$TSTNAME.$SETNAME.$VERSION.run.$BINID
DONEFILE=locks/$TSTNAME.$SETNAME.$VERSION.done

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINID.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

SETTINGS=settings/$SETNAME.set

if [ "$LOCK" == "true" ]
then
    if [ -e $DONEFILE ]
    then
	echo skipping test due to existing done file $DONEFILE
	exit
    fi
    if [ -e $LOCKFILE ]
    then
	if [ -e $RUNFILE ]
        then
	    echo continuing aborted run with run file $RUNFILE
	else
	    echo skipping test due to existing lock file $LOCKFILE
	    exit
	fi
    fi
    date > $LOCKFILE
    date > $RUNFILE
fi

if [ ! -e $OUTFILE ]
then
    CONTINUE=false
fi

if [ "$CONTINUE" == "true" ]
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if [ -e $OUTFILE ]
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if [ -e $ERRFILE ]
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if [ "$CONTINUE" == "true" ]
then
    LASTPROB=`getlastprob.awk $OUTFILE`
    echo Continuing benchmark. Last solved instance: $LASTPROB
    echo "" >> $OUTFILE
    echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
    echo "" >> $OUTFILE
else
    LASTPROB=""
fi

uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

HARDTIMELIMIT=`echo $TIMELIMIT*1.1 | bc`
HARDMEMLIMIT=`echo $MEMLIMIT*1.2 | bc`
echo hard time limit: $HARDTIMELIMIT >>$OUTFILE
echo hard mem limit: $HARDMEMLIMIT >>$OUTFILE

for i in `cat $TSTNAME.test` DONE
do
    if [ "$i" == "DONE" ]
    then
	date > $DONEFILE
	break
    fi

    if [ "$LASTPROB" == "" ]
    then
	LASTPROB=""
	if [ -f $i ]
	then
	    echo @01 $i ===========
	    echo @01 $i ===========                >> $ERRFILE
	    echo set load $SETTINGS                >  $TMPFILE
	    if [ $FEASTOL != "default" ]
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

	    echo -----------------------------
	    date
	    date >>$ERRFILE
	    echo -----------------------------
	    date +"@03 %s"
	    tcsh -c "limit cputime $HARDTIMELIMIT s; limit memoryuse $HARDMEMLIMIT M; limit filesize 200 M; ../$2 < $TMPFILE" 2>>$ERRFILE
#	    ../$2 < $TMPFILE 2>>$ERRFILE
	    date +"@04 %s"
	    echo -----------------------------
	    date
	    date >>$ERRFILE
	    echo -----------------------------
	    echo
	    echo =ready=
	else
	    echo @02 FILE NOT FOUND: $i ===========
	    echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
	fi
    else
	echo skipping $i
	if [ "$LASTPROB" == "$i" ]
	then
	    LASTPROB=""
        fi
    fi
done | tee -a $OUTFILE

rm -f $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

if [ -e $DONEFILE ]
then
    if [ -f $TSTNAME.solu ]
    then
	gawk -f check.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
    else
	gawk -f check.awk -vTEXFILE=$TEXFILE $OUTFILE | tee $RESFILE
    fi

    if [ "$LOCK" == "true" ]
    then
	rm -f $RUNFILE
    fi
fi
