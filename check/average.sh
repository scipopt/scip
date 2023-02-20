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

