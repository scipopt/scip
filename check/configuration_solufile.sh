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

# configures SOLUFILE environment variable for test runs and evaluation runs

# input environment - these environment variables should be set before invoking this script
TSTNAME="${1}"       # name of the test set

# new environment variables defined by this script:
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values

# look for solufiles under the name of the test, the name of the test with everything after the first "_" or "-" stripped, and all;
# prefer more specific solufile names over general ones and the instance database solufiles over those in testset/
SOLUFILE=""
for f in "${TSTNAME}" ${TSTNAME%%_*} ${TSTNAME%%-*} all
do
    for d in instancedata/testsets testset
    do
        if test -f "${d}/${f}.solu"
        then
            SOLUFILE="${d}/${f}.solu"
            break
        fi
    done
    if ! test "${SOLUFILE}" = ""
    then
        break
    fi
done
