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
# $Id: allcmpres.sh,v 1.6 2006/10/20 23:02:10 bzfpfend Exp $

AWKARGS=""
FILES=""
for i in $@
do
    if [ ! -e $i ]
    then
	AWKARGS="$AWKARGS $i"
    else
	FILES="$FILES $i"
    fi
done

for i in `ls -1 --color=none $FILES | sed 's!check\.\([^ .]*\)\.\([^ ]*\)\.res!check.\1!g' | sort -u`
do
    TESTSET=`echo $i | sed 's!results/check\.\([^ .]*\)\..*!\1!g'`
    echo
    echo ====vvvv==== $TESTSET ====vvvv====
    cmpres.awk $AWKARGS `ls -1 --color=none $FILES | grep "$i\..*\.res"`
    echo ====^^^^==== $TESTSET ====^^^^====
done
