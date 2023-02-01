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

# Generate a comparison of two testruns, i.e. with different settings.
#
# Usage: 'allcmpres.sh check.run1.res check.run2.res'
#        'allcmpres.sh check.run1.*.res'

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e "${i}"
    then
        AWKARGS="${AWKARGS} ${i}"
    else
        FILES="${FILES} ${i}"
    fi
done

TESTSETS=""
for i in $(ls -1 ${FILES} | sed 's!\(.*\)check\.\([^ .]*\)\.\([^ ]*\)\.res!\2!g' | sort -u)
do
    TESTSETS="${TESTSETS} ${i}"
done

export LC_NUMERIC=C

for i in ${TESTSETS}
do
    echo
    echo "====vvvv==== ${i} ====vvvv===="
    # the variable AWKARGS needs to be without quotation marks here
    awk -f cmpres.awk ${AWKARGS} texcmpfile="cmpres.${i}.tex" $(ls -1f ${FILES} | grep "${i}\..*\.res")
    echo "====^^^^==== ${i} ====^^^^===="
done
