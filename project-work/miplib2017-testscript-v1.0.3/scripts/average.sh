#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e $i
    then
	AWKARGS="$AWKARGS $i"
    else
	f1=`basename $i .res`
	f2=`basename $i`
	if test "$f1" != "$f2"
	then
	    if test ! -n "$FILES"
	    then
		DIR=`dirname $i`
		AVGFILE="$DIR/$f1-average.res"
	    fi
	    FILES="$FILES $i"
	fi
    fi
done

export LC_NUMERIC=C

if test -n "$FILES"
then
    awk -f scripts/average.awk $AWKARGS $FILES  | tee $AVGFILE
fi

