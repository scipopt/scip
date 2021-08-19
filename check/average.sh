#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# compares averages of several SCIP result files
#
# Usage: average.sh <awkargs> <check* files>

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e "${i}"
    then
        AWKARGS="${AWKARGS} ${i}"
    else
        f1=$(basename "${i}" .res)
        f2=$(basename "${i}")
        if test "${f1}" != "${f2}"
        then
            FILES="${FILES} ${i}"
        fi
    fi
done

export LC_NUMERIC=C

if test -n "${FILES}"
then
    # the variables AWKARGS and FILES need to be without quotation marks here
    awk -f average.awk ${AWKARGS} ${FILES}
fi

