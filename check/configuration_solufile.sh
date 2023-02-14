#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
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
    for d in testset instancedata/testsets
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
