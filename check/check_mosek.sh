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

TSTNAME=$1
MOSEKBIN=$2
SETNAME=$3
BINNAME=$MOSEKBIN.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
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
    LASTPROB=`awk -f getlastprob.awk $OUTFILE`
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

for i in `cat testset/$TSTNAME.test`
do
    if test "$LASTPROB" = ""
    then
        LASTPROB=""
        if test -f $i
        then
            echo @01 $i ===========
            echo @01 $i ===========                     >> $ERRFILE
            echo @05 TIMELIMIT: $TIMELIMIT
            echo @05 THREADS: $THREADS

            ENVPARAM=""
            ENVPARAM="$ENVPARAM -d MSK_IPAR_LOG_MIO_FREQ $DISPFREQ" 
            ENVPARAM="$ENVPARAM -d MSK_DPAR_MIO_MAX_TIME $TIMELIMIT"
            ENVPARAM="$ENVPARAM -d MSK_DPAR_OPTIMIZER_MAX_TIME $TIMELIMIT"
            ENVPARAM="$ENVPARAM -d MSK_IPAR_MIO_MAX_NUM_BRANCHES $NODELIMIT"
#            echo set mip limits treememory $MEMLIMIT >> $TMPFILE
#threads are only used for interior point method, not in b&b tree (version 6)
            ENVPARAM="$ENVPARAM -d MSK_IPAR_INTPNT_NUM_THREADS $THREADS"
 
            if test $FEASTOL != "default"
            then
                ENVPARAM="$ENVPARAM -d  MSK_DPAR_MIO_TOL_REL_GAP $FEASTOL"
            fi
            ENVPARAM="$ENVPARAM -d MSK_DPAR_MIO_NEAR_TOL_REL_GAP 0.0"
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            date +"@03 %s"
            bash -c "ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ulimit -f 1000000; $MOSEKBIN $ENVPARAM -pari $SETTINGS -paro $SETFILE $i" 2>>$ERRFILE
            date +"@04 %s"
            echo -----------------------------
            date
            date >>$ERRFILE
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

./evalcheck_mosek.sh $OUTFILE