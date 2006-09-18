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
# $Id: evalcheck.sh,v 1.1 2006/09/18 04:27:57 bzfpfend Exp $

for i in $@
do
    NAME=`basename $i .out`
    DIR=`dirname $i`
    OUTFILE=$DIR/$NAME.out
    RESFILE=$DIR/$NAME.res
    TEXFILE=$DIR/$NAME.tex
    TMPFILE=$NAME.tmp

    echo $NAME >$TMPFILE
    TSTNAME=`sed 's/check.\([a-zA-Z0-9_]*\).*/\1/g' $TMPFILE`
    rm $TMPFILE

    if [ -f $TSTNAME.solu ]
    then
	gawk -f check.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
    else
	gawk -f check.awk -vTEXFILE=$TEXFILE $OUTFILE | tee $RESFILE
    fi
done
