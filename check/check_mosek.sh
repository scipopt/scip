#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

TSTNAME=$1
MOSEKBIN=$2
SETNAME=$3
BINNAME=$MOSEKBIN.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8
MIPGAP=$9
DISPFREQ=${10}
CONTINUE=${11}

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

SETTINGS=settings/$SETNAME.mskset

if test "$CONTINUE" = "true"
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

if test "$CONTINUE" = "true"
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

# we add 10% to the hard time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 10\` + \`expr $TIMELIMIT / 10\``

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`

echo "hard time limit: $HARDTIMELIMIT s" >>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT k" >>$OUTFILE

for i in `cat $TSTNAME.test`
do
    if test "$LASTPROB" = ""
    then
	LASTPROB=""
	if test -f $i
	then
	    rm -f $SETFILE
	    echo @01 $i ===========
	    echo @01 $i ===========                     >> $ERRFILE

	    echo BEGIN MOSEK                            > $TMPFILE
	    echo MSK_IPAR_LOG_MIO_FREQ $DISPFREQ        >> $TMPFILE
            echo MSK_DPAR_MIO_MAX_TIME $TIMELIMIT       >> $TMPFILE
#	    echo set mip limits nodes $NODELIMIT    >> $TMPFILE
#	    echo set mip limits treememory $MEMLIMIT >> $TMPFILE
	    if test $FEASTOL != "default"
	    then
		echo MSK_DPAR_MIO_TOL_REL_GAP $FEASTOL  >> $TMPFILE
	    fi
	    echo MSK_DPAR_MIO_NEAR_TOL_REL_GAP $MIPGAP  >> $TMPFILE
	    echo END MOSEK                              >> $TMPFILE
	    cp $TMPFILE $SETFILE
	    echo -----------------------------
	    date
	    echo -----------------------------
	    bash -c "ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ulimit -f 1000000; $MOSEKBIN -p $TMPFILE $i" 2>>$ERRFILE
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
	if test "$LASTPROB" = "$i"
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
    awk -f check_mosek.awk -v "TEXFILE=$TEXFILE" $TSTNAME.solu $OUTFILE | tee $RESFILE
else
    awk -f check_mosek.awk -v "TEXFILE=$TEXFILE"$OUTFILE | tee $RESFILE
fi
