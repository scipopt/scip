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

# Calls check_count.awk on the testrun files and writes the output in a .res file.
#
# To be invoked by 'check_count.sh'.

export LANG=C

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e ${i}
    then
        AWKARGS="${AWKARGS} ${i}"
    else
        FILES="${FILES} ${i}"
    fi
done

for i in ${FILES}
do
    NAME=$(basename ${i} .out)
    DIR=$(dirname ${i})
    OUTFILE="${DIR}/${NAME}.out"
    RESFILE="${DIR}/${NAME}.res"
    TEXFILE="${DIR}/${NAME}.tex"
    PAVFILE="${DIR}/${NAME}.pav"

    TSTNAME=$(echo "${NAME}" | sed 's/checkcount.\([a-zA-Z0-9_]*\).*/\1/g')

    if test -f "testset/${TSTNAME}.test"
    then
        TESTFILE="testset/${TSTNAME}.test"
    else
        TESTFILE=""
    fi

    # call method to obtain solution file
    # defines the following environment variable: SOLUFILE
    . ./configuration_solufile.sh "${TSTNAME}"

    # the variable AWKARGS needs to be without quotation marks here
    awk -f check_count.awk ${AWKARGS} "${TESTFILE}" "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
done
