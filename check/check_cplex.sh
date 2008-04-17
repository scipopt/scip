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

TSTNAME=$1
CPLEXBIN=$2
SETNAME=$3
BINNAME=$CPLEXBIN.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8
MIPGAP=$9
CONTINUE=${10}

if test ! -e results
then
    mkdir results
fi
if test ! -e settings
then
    mkdir settings
fi

OUTFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.prm

SETTINGS=settings/$SETNAME.cpxset

if test "$CONTINUE" == "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if test "$CONTINUE" == "true"
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

for i in `cat $TSTNAME.test`
do
    if test "$LASTPROB" == ""
    then
	LASTPROB=""
	if test -f $i
	then
	    rm -f $SETFILE
	    echo @01 $i ===========
	    echo @01 $i ===========                 >> $ERRFILE
	    if test $SETNAME != "default"
	    then
#	        echo read $SETTINGS                  > $TMPFILE
		cp $SETTINGS $TMPFILE
	    else
		echo ""                              > $TMPFILE
	    fi
	    if test $FEASTOL != "default"
	    then
		echo set simplex tolerances feas $FEASTOL    >> $TMPFILE
		echo set mip tolerances integrality $FEASTOL >> $TMPFILE
	    fi
	    echo set timelimit $TIMELIMIT           >> $TMPFILE
	    echo set clocktype 1                    >> $TMPFILE
	    echo set mip display 3                  >> $TMPFILE
	    echo set mip interval 10000             >> $TMPFILE
	    if test $MIPGAP != "default"
	    then
		echo set mip tolerances mipgap $MIPGAP >> $TMPFILE
	    fi
	    echo set mip limits nodes $NODELIMIT    >> $TMPFILE
	    echo set mip limits treememory $MEMLIMIT >> $TMPFILE
	    echo write $SETFILE                     >> $TMPFILE
	    echo read $i                            >> $TMPFILE
	    echo display problem stats              >> $TMPFILE
	    echo optimize                           >> $TMPFILE
	    echo quit                               >> $TMPFILE
	    echo -----------------------------
	    date
	    echo -----------------------------
	    tcsh -c "limit cputime $HARDTIMELIMIT s; limit memoryuse $HARDMEMLIMIT M; limit filesize 1000 M; $CPLEXBIN < $TMPFILE" 2>>$ERRFILE
	    echo -----------------------------
	    date
	    echo -----------------------------
	    echo =ready=
	else
	    echo @02 FILE NOT FOUND: $i ===========
	    echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
	fi
    else
	echo skipping $i
	if test "$LASTPROB" == "$i"
	then
	    LASTPROB=""
        fi
    fi
done | tee -a $OUTFILE

rm -f $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

if test -f $TSTNAME.solu
then
    gawk -f check_cplex.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
else
    gawk -f check_cplex.awk -vTEXFILE=$TEXFILE $OUTFILE | tee $RESFILE
fi
