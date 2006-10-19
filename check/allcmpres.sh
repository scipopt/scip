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
# $Id: allcmpres.sh,v 1.5 2006/10/19 19:43:29 bzfpfend Exp $

for i in `ls -1 --color=none $@ | sed 's!check\.\([^ .]*\)\.\([^ ]*\)\.res!check.\1!g' | sort -u`
do
    TESTSET=`echo $i | sed 's!results/check\.\([^ .]*\)\..*!\1!g'`
    echo
    echo ====vvvv==== $TESTSET ====vvvv====
    cmpres.awk `ls -1 --color=none $@ | grep "$i\..*\.res"`
    echo ====^^^^==== $TESTSET ====^^^^====
done
