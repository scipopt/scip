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
# $Id: check_cplex.sh,v 1.18 2006/08/21 20:13:17 bzfpfend Exp $
TSTNAME=$1
CPLEXBIN=$2
SETNAME=$3
BINNAME=$CPLEXBIN.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
FEASTOL=$8

OUTFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.prm

SETTINGS=settings/$SETNAME.cpxset

uname -a >$OUTFILE
uname -a >$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

for i in `cat $TSTNAME.test`
do
    if [ -f $i ]
    then
        rm -f $SETFILE
	echo @01 $i ===========
	echo @01 $i ===========                 >> $ERRFILE
	if [ $SETNAME != "default" ]
	then
	    cp $SETTINGS $TMPFILE
	else
	    echo ""                              > $TMPFILE
	fi
	if [ $FEASTOL != "default" ]
	then
	    echo set simplex tolerances feas $FEASTOL    >> $TMPFILE
	    echo set mip tolerances integrality $FEASTOL >> $TMPFILE
	fi
	echo set timelimit $TIMELIMIT           >> $TMPFILE
	echo set clocktype 1                    >> $TMPFILE
	echo set mip display 3                  >> $TMPFILE
	echo set mip interval 10000             >> $TMPFILE
	echo set mip tolerances absmipgap 1e-9  >> $TMPFILE
	echo set mip tolerances mipgap 0.0      >> $TMPFILE
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
	$CPLEXBIN < $TMPFILE 2>>$ERRFILE
	echo -----------------------------
	date
	echo -----------------------------
	echo =ready=
    else
	echo @02 FILE NOT FOUND: $i ===========
	echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
    fi
done | tee -a $OUTFILE

rm $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

if [ -f $TSTNAME.solu ]
then
    gawk -f check_cplex.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
else
    gawk -f check_cplex.awk -vTEXFILE=$TEXFILE $OUTFILE | tee $RESFILE
fi
